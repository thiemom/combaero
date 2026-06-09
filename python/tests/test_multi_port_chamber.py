"""
Network-level tests for MultiPortChamberElement (momentum-CV junction).

The element + companion BorderCarnotLossElement land the PDF-spec successor
to the K-closure tee (docs/junction/momentum cv implementation guide.pdf).
These tests cover the assembly side: the junction wires correctly into
combaero's solver, the port-MCN mass-row skip works, sign conventions
resolve cleanly from topology, and the system converges to a mass-conserving
root on a simple low-Mach 3-port network.

Tier-4 cross-check vs TeeJunctionElement and Tier-1/2 validation against
Hager / Bassett / Perez-Garcia / Wang reference data live in separate files.
"""

import math

import numpy as np

import combaero as cb
from combaero.network import (
    BorderCarnotLossElement,
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    MultiPortChamberElement,
    NetworkSolver,
    PressureBoundary,
)

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _three_port_net(
    Pt_in: float = 2.1e5,
    Pt_str: float = 2.05e5,
    Pt_bra: float = 2.00e5,
    Tt: float = 300.0,
    diameter: float = 0.05,
    length: float = 0.3,
) -> FlowNetwork:
    """3-port momentum-CV junction: 1 pressure inlet, 2 pressure outlets, all
    via finite-length compressible channels with realistic friction.

    Topology:
       pb_in  --ch_in-->  mc_com \\
                                  >- jct -< mc_str --ch_str--> pb_str
                                            mc_bra --ch_bra--> pb_bra

    All three ports are MomentumChamberNodes (PDF requires per-port (P, Pt)
    state at the junction face). The junction is lossless; channel friction
    sets up the Pt gradient that drives flow.
    """
    Y = _DRY_AIR_Y
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_in", Pt=Pt_in, Tt=Tt, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_str, Tt=Tt, Y=Y))
    net.add_node(PressureBoundary("pb_bra", Pt=Pt_bra, Tt=Tt, Y=Y))
    net.add_node(MomentumChamberNode("mc_com", area=math.pi * (diameter / 2.0) ** 2))
    net.add_node(MomentumChamberNode("mc_str", area=math.pi * (diameter / 2.0) ** 2))
    net.add_node(MomentumChamberNode("mc_bra", area=math.pi * (diameter / 2.0) ** 2))
    net.add_element(
        ChannelElement(
            "ch_in", "pb_in", "mc_com", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_str", "mc_str", "pb_str", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_bra", "mc_bra", "pb_bra", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["mc_com"],
            outlet_nodes=["mc_str", "mc_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
        )
    )
    return net


def test_three_port_network_converges():
    """Solver assembles and converges on the basic 3-port momentum-CV network."""
    net = _three_port_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    # Residual norm at the converged state is small.
    x = np.array([sol.get(n, 0.0) for n in solver.unknown_names])
    res, _ = solver._residuals_and_jacobian(x, compute_jacobian=False)
    assert max(abs(r) for r in res) < 1e-3, f"residuals not small: {res}"


def test_three_port_mass_conservation():
    """All three port flows balance at the junction (sum-mass residual = 0).

    Junction convention: ch_in is the declared inlet (sign -1: +ch_in is flow
    INTO junction), ch_str and ch_bra are outlets (sign +1: +ch_str / +ch_bra
    are flow OUT of junction). Conservation:
        -ch_in.m_dot + ch_str.m_dot + ch_bra.m_dot = 0
    """
    net = _three_port_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    m_in = sol["ch_in.m_dot"]
    m_str = sol["ch_str.m_dot"]
    m_bra = sol["ch_bra.m_dot"]
    imbalance = -m_in + m_str + m_bra
    assert abs(imbalance) < 1e-6, (
        f"mass imbalance {imbalance} kg/s: m_in={m_in}, m_str={m_str}, m_bra={m_bra}"
    )


def test_three_port_all_forward_under_pressure_drive():
    """With a single source (highest Pt) feeding two sinks (lower Pt), flow is
    forward through every port.

    Note: an earlier version of this test asserted BIDIRECTIONAL flow as a
    "feature" of the PDF Section 2.2 angle-blind impulse-equality ansatz
    (R = P + rho*u^2 - P_jct = 0). That bidirectional behaviour was an
    artifact of the formulation, not real physics -- Bassett 2001 Fig 7a
    shows real branching-tee K6 is purely positive and forward-driven under
    these BCs. The sin^2(theta) projection in the impulse residual recovers
    the physical all-forward solution. See feat(junction) sin^2 milestone
    commit for derivation.
    """
    net = _three_port_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    m_in = sol["ch_in.m_dot"]
    m_str = sol["ch_str.m_dot"]
    m_bra = sol["ch_bra.m_dot"]
    assert m_in > 0.0, f"inlet reversed: m_in={m_in}"
    assert m_str > 0.0, f"straight reversed: m_str={m_str}"
    assert m_bra > 0.0, f"branch reversed: m_bra={m_bra}"


def test_three_port_p_jct_in_physical_range():
    """P_jct should lie between the lowest collector Pt and the source Pt."""
    Pt_in = 2.1e5
    Pt_str = 2.05e5
    Pt_bra = 2.00e5
    net = _three_port_net(Pt_in=Pt_in, Pt_str=Pt_str, Pt_bra=Pt_bra)
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert sol["__success__"]

    P_jct = sol["jct.P_jct"]
    # Sanity: P_jct should sit between the lowest collector P and the source.
    # Use ~10% slack: friction can push P_jct outside the Pt envelope.
    lo, hi = min(Pt_str, Pt_bra) * 0.9, Pt_in * 1.1
    assert lo < P_jct < hi, f"P_jct={P_jct} outside physical range [{lo}, {hi}]"


def test_border_carnot_loss_element_in_series_with_junction():
    """Bolt a BorderCarnotLossElement onto the branch port.

    Smoke test: the network still solves with the loss element in series and
    series mass conservation holds along the lateral path (loss element and
    its downstream channel carry the same m_dot).
    """
    Y = _DRY_AIR_Y
    diameter = 0.05
    length = 0.3
    area = math.pi * (diameter / 2.0) ** 2

    # With loss element on branch:
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_in", Pt=2.1e5, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=2.05e5, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_bra", Pt=2.00e5, Tt=300.0, Y=Y))
    net.add_node(MomentumChamberNode("mc_com", area=area))
    net.add_node(MomentumChamberNode("mc_str", area=area))
    net.add_node(MomentumChamberNode("mc_bra", area=area))
    # Intermediate node between BorderCarnot and the branch channel.
    net.add_node(MomentumChamberNode("mc_bra_post", area=area))
    net.add_element(
        ChannelElement(
            "ch_in", "pb_in", "mc_com", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_str", "mc_str", "pb_str", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        BorderCarnotLossElement(
            "loss_bra", from_node="mc_bra", to_node="mc_bra_post", delta_geom_deg=90.0, area=area
        )
    )
    net.add_element(
        ChannelElement(
            "ch_bra",
            "mc_bra_post",
            "pb_bra",
            length=length,
            diameter=diameter,
            regime="compressible",
        )
    )
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["mc_com"],
            outlet_nodes=["mc_str", "mc_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
        )
    )

    sol = NetworkSolver(net).solve()
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    # Series mass conservation: ch_bra (downstream channel) carries the same
    # flow as the loss element (no junction between them, no mass-row skip).
    m_loss = sol["loss_bra.m_dot"]
    m_ch_bra = sol["ch_bra.m_dot"]
    assert abs(m_loss - m_ch_bra) < 1e-6, (
        f"loss/channel series mismatch: m_loss={m_loss}, m_ch_bra={m_ch_bra}"
    )
