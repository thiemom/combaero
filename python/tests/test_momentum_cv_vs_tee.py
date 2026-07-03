"""
Tier 4 cross-check: MultiPortChamberElement vs TeeJunctionElement at low Mach.

Status: scoped down from "literal K-coefficient agreement" to "qualitative
agreement in physical behaviour where both models share an unambiguous
answer". The two models cannot be compared numerically point-for-point:

- TeeJunctionElement's K-closure uses a stagnation-pressure reference and
  encodes turning loss directly in K_branch. It does NOT support bidirectional
  flow at the lateral port -- the K-correlation forces direction.
- MultiPortChamberElement's impulse-equality residuals are sign-free and
  support flow reversal natively (the architectural win, PDF Section 2 and
  addendum Finding 6). Turning loss is on a separate BorderCarnotLossElement.

At near-symmetric BCs the two models give qualitatively DIFFERENT answers
(MPCE happily routes flow between collectors; v3 forces all-forward). That
divergence is the feature, not a bug. So a strict numerical cross-check
would either be vacuous (uniform-Pt zero-flow) or test a regime where the
models genuinely disagree.

What's worth verifying here:
- Both models recover the trivial zero-flow root at uniform Pt.
- At a strongly asymmetric one-source / two-sink setup, both models produce
  an all-forward dominant-inlet flow (this is the regime where the K-closure
  reproduces Hager / Bassett incompressible-flow results that the PDF promises
  the momentum-CV matches at M -> 0).

Tier-1 (Hager / Bassett digitised CSVs under docs/bassett/) and Tier-2
(Perez-Garcia / Wang compressible reference data under docs/junction/) are
the strict numerical validations; they live in separate test files (TBD).
This file is the qualitative agreement check.
"""

import math

import combaero as cb
from combaero.network import (
    BorderCarnotLossElement,
    FlowNetwork,
    LosslessConnectionElement,
    MomentumChamberNode,
    MultiPortChamberElement,
    NetworkSolver,
    PressureBoundary,
    TeeJunctionElement,
)

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_THETA_90 = math.pi / 2.0
_F_C = 0.01  # 100 cm^2


def _v3_branching_net_direct(
    Pt_common: float = 2.1e5,
    Pt_straight: float = 2.098e5,
    Pt_branch: float = 2.06e5,
) -> FlowNetwork:
    """v3 K-closure tee, native topology (direct PB-to-PB).

    Matches the test_tee_network._branching_net default fixture. Uses BCs
    inside the K-closure's K_run ~ 0.30 envelope (addendum Finding 6).
    """
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(PressureBoundary("pb_common", Pt=Pt_common, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_straight", Pt=Pt_straight, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_branch", Pt=Pt_branch, Tt=300.0, Y=Y))
    net.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pb_common",
            straight_node="pb_straight",
            branch_node="pb_branch",
            theta=_THETA_90,
            F_C=_F_C,
            psi=1.0,
            tee_type="branching",
        )
    )
    return net


def _mpce_branching_net(
    Pt_common: float = 2.1e5,
    Pt_straight: float = 2.098e5,
    Pt_branch: float = 2.06e5,
) -> FlowNetwork:
    """MPCE with port-MCNs and a BorderCarnotLoss on the branch.

    MPCE structurally requires MomentumChamberNode port faces (PDF Section 2);
    direct PB-to-PB wiring is not supported. We use lossless connections from
    the boundary PBs to the port-MCNs, mirroring _v3_branching_net's BCs.
    """
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(PressureBoundary("pb_common", Pt=Pt_common, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_straight", Pt=Pt_straight, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_branch", Pt=Pt_branch, Tt=300.0, Y=Y))
    net.add_node(MomentumChamberNode("port_com", area=_F_C))
    net.add_node(MomentumChamberNode("port_str", area=_F_C))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C))
    net.add_node(MomentumChamberNode("port_bra_post", area=_F_C))
    net.add_element(LosslessConnectionElement("lc_com", "pb_common", "port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_straight"))
    net.add_element(
        BorderCarnotLossElement(
            "loss_bra",
            from_node="port_bra",
            to_node="port_bra_post",
            delta_geom_deg=90.0,
            area=_F_C,
        )
    )
    net.add_element(LosslessConnectionElement("lc_bra", "port_bra_post", "pb_branch"))
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["port_com"],
            outlet_nodes=["port_str", "port_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[_F_C, _F_C, _F_C],
        )
    )
    return net


def test_both_models_recover_zero_flow_at_uniform_pt():
    """At uniform Pt (no drive), both models give zero flow."""
    Pt = 1.5e5
    v3_sol = NetworkSolver(
        _v3_branching_net_direct(Pt_common=Pt, Pt_straight=Pt, Pt_branch=Pt)
    ).solve()
    mp_sol = NetworkSolver(_mpce_branching_net(Pt_common=Pt, Pt_straight=Pt, Pt_branch=Pt)).solve()

    assert v3_sol["__success__"]
    assert mp_sol["__success__"]
    assert abs(v3_sol["tee.m_dot_com"]) < 1e-3, f"v3 nonzero flow: {v3_sol['tee.m_dot_com']}"
    assert abs(mp_sol["lc_com.m_dot"]) < 1e-3, f"MPCE nonzero flow: {mp_sol['lc_com.m_dot']}"


def test_both_models_give_forward_dominant_inlet_flow():
    """At a one-source / two-sink BC inside the K-closure envelope, both
    models converge with positive common inflow.

    This is the regime where the K-closure works correctly (per PDF addendum
    Finding 6: pb_str at 2.098e5 puts the straight drop inside K_run's 0.30
    envelope) and where the momentum-CV with the (3/4)-correction lateral
    loss should also reproduce Hager-like behaviour at low Mach.

    We do NOT assert the m_str signs agree: at these BCs the K-closure forces
    all-forward, while MPCE may produce bidirectional flow at the straight
    collector if the impulse-equality dynamics favour it (and that is
    physically valid).
    """
    v3_sol = NetworkSolver(_v3_branching_net_direct()).solve()

    # MPCE's residual set is even in the mass flows apart from the linear
    # mass row, so the network is symmetric under a global sign flip: the
    # all-forward root and its mirrored all-reverse image are equally exact,
    # and the default cold start may land either. Seed the forward basin
    # explicitly -- this test compares the models' forward solutions, it
    # does not test cold-start basin selection.
    mp_net = _mpce_branching_net()
    mp_net.elements["lc_com"].initial_guess = {"lc_com.m_dot": 0.4}
    mp_net.elements["lc_str"].initial_guess = {"lc_str.m_dot": 0.25}
    mp_net.elements["lc_bra"].initial_guess = {"lc_bra.m_dot": 0.15}
    mp_sol = NetworkSolver(mp_net).solve()

    assert v3_sol["__success__"], f"v3 did not converge: {v3_sol.get('__message__')}"
    assert mp_sol["__success__"], f"MPCE did not converge: {mp_sol.get('__message__')}"

    # Common inflow is forward in both.
    m_com_v3 = v3_sol["tee.m_dot_com"]
    m_com_mp = mp_sol["lc_com.m_dot"]
    assert m_com_v3 > 0.0, f"v3 common reversed: {m_com_v3}"
    assert m_com_mp > 0.0, f"MPCE common reversed: {m_com_mp}"

    # Magnitudes should be of the same order. Loose order-of-magnitude check.
    # The K-closure and the momentum-CV with Hager loss can differ by O(1) at
    # finite Mach (different reference dynamic heads); we only check they are
    # within a factor of 5 of each other.
    ratio = m_com_v3 / m_com_mp
    assert 0.2 < ratio < 5.0, (
        f"common-flow magnitudes disagree by > 5x: v3={m_com_v3:.4f}, MPCE={m_com_mp:.4f}"
    )
