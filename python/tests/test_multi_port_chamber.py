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
import pytest

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


@pytest.fixture(scope="module")
def three_port_solution() -> tuple[NetworkSolver, dict]:
    """One shared solve of the default 3-port net for the test_three_port_*
    assertions below (each previously re-ran the identical ~3 min solve).

    The primary attempt on this net is known-doomed (hybr stalls, the LM
    fallback stalls too) and the outlet-referenced warm-start auto-retry
    rescues it in seconds; the finite timeout bounds the doomed phase and
    pins that rescue-within-budget behavior.
    """
    net = _three_port_net()
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=120.0)
    return solver, sol


def test_three_port_network_converges(three_port_solution):
    """Solver assembles and converges on the basic 3-port momentum-CV network."""
    solver, sol = three_port_solution
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    # Residual norm at the converged state is small.
    x = np.array([sol.get(n, 0.0) for n in solver.unknown_names])
    res, _ = solver._residuals_and_jacobian(x, compute_jacobian=False)
    assert max(abs(r) for r in res) < 1e-3, f"residuals not small: {res}"

    # Document the rescue path: this net converges via the outlet-ref
    # warm-start auto-retry after the doomed primary fails fast.
    ssu = solver._diagnostic_data["solver_settings_used"]
    assert ssu["init_strategy"] == "outletref_warmstart"
    assert ssu.get("auto_retry") is True


def test_three_port_doomed_primary_fails_fast_via_lm_stall():
    """Without the auto-retry, the doomed primary must fail FAST: hybr
    stalls (handover), the LM fallback then plateaus far from tolerance
    and the LM-phase stall detector aborts it. Pre-detector this burned
    ~187 s flat at |F| ~ 4e3 (timeout=None leaves LM unbounded); with it
    the attempt ends ~10 s after LM flatlines.
    """
    import time

    net = _three_port_net()
    solver = NetworkSolver(net)
    t0 = time.perf_counter()
    with pytest.warns(UserWarning, match="did not converge"):
        sol = solver.solve(timeout=None, auto_retry=False)
    wall = time.perf_counter() - t0
    assert not sol["__success__"]
    diag = solver._diagnostic_data
    assert diag["stall_handover"] is True
    assert diag["lm_stall"] is True
    assert "LM fallback also stalled" in str(sol.get("__message__"))
    assert wall < 60.0, f"doomed primary took {wall:.1f} s; LM-stall abort broken"


def test_three_port_mass_conservation(three_port_solution):
    """All three port flows balance at the junction (sum-mass residual = 0).

    Junction convention: ch_in is the declared inlet (sign -1: +ch_in is flow
    INTO junction), ch_str and ch_bra are outlets (sign +1: +ch_str / +ch_bra
    are flow OUT of junction). Conservation:
        -ch_in.m_dot + ch_str.m_dot + ch_bra.m_dot = 0
    """
    _, sol = three_port_solution
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    m_in = sol["ch_in.m_dot"]
    m_str = sol["ch_str.m_dot"]
    m_bra = sol["ch_bra.m_dot"]
    imbalance = -m_in + m_str + m_bra
    assert abs(imbalance) < 1e-6, (
        f"mass imbalance {imbalance} kg/s: m_in={m_in}, m_str={m_str}, m_bra={m_bra}"
    )


def test_three_port_converges_with_conservation(three_port_solution):
    """Network converges and conserves mass.

    Directional behaviour at asymmetric pressure BCs depends on the cross-
    coupling: at sufficiently asymmetric collector Pt the junction acts as
    an ejector and the lower-Pt branch may source flow INTO the junction
    (water-jet-pump effect). This is correct physics per Bassett 2001
    Section 3 (and the basis of the K6 < 0 region at low q). We assert
    convergence + mass conservation only; direction is BC-dependent.
    """
    _, sol = three_port_solution
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    m_in = sol["ch_in.m_dot"]
    m_str = sol["ch_str.m_dot"]
    m_bra = sol["ch_bra.m_dot"]
    # Common inlet always forward at these BCs (common Pt is highest).
    assert m_in > 0.0, f"inlet reversed: m_in={m_in}"
    # Mass conservation at the junction.
    imbalance = m_in - m_str - m_bra
    assert abs(imbalance) < 1e-3, (
        f"mass imbalance {imbalance}: m_in={m_in}, m_str={m_str}, m_bra={m_bra}"
    )


def test_three_port_p_jct_in_physical_range(three_port_solution):
    """P_jct should lie between the lowest collector Pt and the source Pt."""
    # Fixture defaults (see _three_port_net signature).
    Pt_in = 2.1e5
    Pt_str = 2.05e5
    Pt_bra = 2.00e5
    _, sol = three_port_solution
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
