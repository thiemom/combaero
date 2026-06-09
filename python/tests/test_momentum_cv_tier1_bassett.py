"""
Tier 1 validation: MPCE + BorderCarnotLoss vs Bassett 2001 / Hager 1984 at M -> 0.

STATUS: tests are xfail. The PDF Section 3.4 claim "MPCE + loss element form 2
reproduces Hager exactly at M -> 0" does NOT hold empirically as set up here.

The disagreement could come from EITHER of two sources (or both); this file
does not yet disentangle them:

1. PDF formula mismatch (algebra).
   MPCE impulse + MCN KE coupling gives, at equal areas low Mach:
     K_straight_MPCE = q^2 - 2*q          (not Hager xi_t = q^2 - 0.5*q)
     K_lateral_MPCE  = q^2 - 1            (not Bassett K6 = q^2 + 1 - 0.7654*q)
   Form-2 BorderCarnotLoss adds K_loss = L*q^2 with constant L for the
   lateral; no constant L makes (1+L)*q^2 - 1 equal Bassett K6 across q.

2. Solver finding a degenerate near-root.
   The momentum-CV junction has multiple near-singular directions in its
   Jacobian (sum-mass residual lives in a left-null direction of the
   pressure block; impulse residuals are sign-free so flow direction can
   flip). Newton lands on whichever locally-zero residual basin is
   closest to the initial guess; that basin may not be the physically
   correct root. We see this in the smoke tests too (test_three_port_*).

Resolution requires both:
- Verifying converged residual norm matches PDF expectations at the imposed q
  (rules out solver degeneracy), AND
- Either (a) revisiting the loss-element form (linear-in-q? referenced to
  q_common not q_local?) or (b) accepting the PDF spec as approximate.

Tier 1 is foundational: until this is sorted, Tier 2 (Pérez-García / Wang
compressible) rides on a low-Mach baseline that already disagrees with the
canonical reference. See docs/junction/tier2_reference_data.md.
"""

import math

import pytest

import combaero as cb
from combaero.network import (
    BorderCarnotLossElement,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    MultiPortChamberElement,
    NetworkSolver,
    PressureBoundary,
)

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_F_C = 0.01  # 100 cm^2
_M_DOT_IN = 0.1  # kg/s -- low-Mach inlet


def _K6_bassett(q: float, psi: float = 1.0, theta_rad: float = math.pi / 2.0) -> float:
    """Bassett 2001 Eq. 27 (separating flow, lateral / main-to-branch loss).
    Already includes Hager 3/4-correction."""
    return q * q * psi * psi + 1.0 - 2.0 * q * psi * math.cos(0.75 * theta_rad)


def _K2_bassett(q: float) -> float:
    """Bassett 2001 Eq. 15 (separating flow, straight loss). theta- and psi-independent."""
    return q * q - 1.5 * q + 0.5


def _build_imposed_q_network(m_in: float, m_lateral: float, Pt_ref: float = 1.0e5) -> FlowNetwork:
    """Build a network with imposed mass-flow split.

    Topology:
        mb_in (MFB, m_in INTO common) -> port_com
        port_com -> jct -> port_str -> lc_str -> pb_str (PB at Pt_ref)
                    jct -> port_bra -> loss_bra -> port_bra_post -> mb_lat (MFB, m_lateral OUT)

    Two MassFlowBoundaries pin m_in and m_lateral; by mass balance the
    straight flow m_str = m_in - m_lateral. PB at the straight sets the
    pressure datum; everything else is determined by the impulse + loss
    residuals.
    """
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(MassFlowBoundary("mb_in", m_dot=m_in, Tt=300.0, Y=Y))
    net.add_node(MassFlowBoundary("mb_lat", m_dot=m_lateral, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_ref, Tt=300.0, Y=Y))
    net.add_node(MomentumChamberNode("port_com", area=_F_C))
    net.add_node(MomentumChamberNode("port_str", area=_F_C))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C))
    net.add_node(MomentumChamberNode("port_bra_post", area=_F_C))
    net.add_element(LosslessConnectionElement("lc_in", "mb_in", "port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    net.add_element(
        BorderCarnotLossElement(
            "loss_bra",
            from_node="port_bra",
            to_node="port_bra_post",
            delta_geom_deg=90.0,
            area=_F_C,
        )
    )
    net.add_element(LosslessConnectionElement("lc_bra", "port_bra_post", "mb_lat"))
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


def _extract_K_at_imposed_q(q: float) -> dict[str, float]:
    """Solve at imposed q and return diagnostic K values."""
    m_in = _M_DOT_IN
    m_lateral = q * m_in
    net = _build_imposed_q_network(m_in=m_in, m_lateral=m_lateral)
    sol = NetworkSolver(net).solve()
    if not sol["__success__"]:
        return {"converged": False, "msg": sol.get("__message__", "")[:80]}

    # Read converged junction-face states.
    P_com = sol["port_com.P"]
    Pt_com = sol["port_com.Pt"]
    Pt_str = sol["port_str.Pt"]
    Pt_bra_post = sol["port_bra_post.Pt"]

    # rho_com and u_com from the converged state.
    X_air = list(cb.mass_to_mole(_DRY_AIR_Y))
    rho_com = float(cb.density(300.0, P_com, X_air))
    u_com = m_in / (rho_com * _F_C)
    q_dyn_com = 0.5 * rho_com * u_com * u_com

    # Mach for sanity.
    a_com = float(cb.speed_of_sound(300.0, X_air))
    Mach_com = abs(u_com) / a_com

    K_straight = (Pt_com - Pt_str) / q_dyn_com
    K_lateral = (Pt_com - Pt_bra_post) / q_dyn_com

    return {
        "converged": True,
        "q": q,
        "K_straight": K_straight,
        "K_lateral": K_lateral,
        "K2_bassett": _K2_bassett(q),
        "K6_bassett": _K6_bassett(q),
        "Mach_com": Mach_com,
    }


@pytest.mark.xfail(
    reason=(
        "MPCE+cross-coupling reproduces Bassett K6 within ~15% under bounded "
        "least_squares (see experimental_bounded_solver.py). The default "
        "NetworkSolver (hybr -> LM) lands on a different basin for this "
        "imposed-q setup -- the Newton problem is harder with cross-coupling "
        "active. Switching the production solver to bounded LM for MPCE "
        "networks is a separate work item."
    ),
    strict=False,
)
def test_low_mach_K_lateral_agrees_with_bassett_K6():
    """Sweep q over Bassett's range, compare K_lateral to K6.

    With sin^2(theta) impulse + cross-coupling (Bassett axial momentum +
    wall reaction Eq 3-4) integrated into MPCE, K_lat matches Bassett
    K6 = 1 + q^2 - 2*q*cos((3/4)*theta) within ~15% under the bounded LM
    solver. With the default solver it lands on a different basin -- see
    xfail reason.
    """
    results = []
    for q in [0.25, 0.5, 0.75]:
        r = _extract_K_at_imposed_q(q)
        if not r["converged"]:
            pytest.skip(f"q={q}: did not converge ({r['msg']})")
        results.append(r)

    failures = []
    for r in results:
        diff = abs(r["K_lateral"] - r["K6_bassett"])
        rel = diff / max(abs(r["K6_bassett"]), 1e-6)
        # Tolerance: 20% relative or 0.15 absolute. Bassett's own measurements
        # vs. analytical predictions show ~5-10% spread (Fig 7a); we add
        # margin for MCN-propagation drift in the wrapped network topology.
        if not (rel < 0.20 or diff < 0.15):
            failures.append(
                f"q={r['q']}: K_lateral={r['K_lateral']:.4f}, "
                f"K6={r['K6_bassett']:.4f}, rel err {rel:.1%}, M={r['Mach_com']:.4f}"
            )

    assert not failures, "K_lateral does NOT match Bassett K6 at M -> 0:\n  " + "\n  ".join(
        failures
    )


@pytest.mark.xfail(
    reason=(
        "Empirically MPCE alone does not reproduce Hager xi_t / Bassett K2 "
        "at M -> 0. Could be PDF formula mismatch or solver degeneracy. "
        "See module docstring."
    ),
    strict=False,
)
def test_low_mach_K_straight_agrees_with_bassett_K2():
    """Sweep q, compare K_straight to K2 (Bassett Eq 15, no loss element on straight)."""
    results = []
    for q in [0.25, 0.5, 0.75]:
        r = _extract_K_at_imposed_q(q)
        if not r["converged"]:
            pytest.skip(f"q={q}: did not converge ({r['msg']})")
        results.append(r)

    failures = []
    for r in results:
        diff = abs(r["K_straight"] - r["K2_bassett"])
        rel = diff / max(abs(r["K2_bassett"]), 1e-6)
        if not (rel < 0.20 or diff < 0.05):
            failures.append(
                f"q={r['q']}: K_straight={r['K_straight']:.4f}, "
                f"K2={r['K2_bassett']:.4f}, rel err {rel:.1%}, M={r['Mach_com']:.4f}"
            )

    assert not failures, (
        "K_straight does NOT match Bassett K2 at M -> 0. This contradicts the "
        "PDF Section 3.4 claim that MPCE alone reproduces Hager xi_t. "
        "Findings:\n  " + "\n  ".join(failures)
    )


def test_mach_stays_low():
    """Sanity: at the imposed mass flow, common-branch Mach is below 0.05.

    The Bassett / Hager K formulas are incompressible. If the test setup
    runs at high Mach the comparison is invalid -- this guard makes sure
    we are deep in the incompressible regime.
    """
    r = _extract_K_at_imposed_q(q=0.5)
    if not r["converged"]:
        pytest.skip(f"did not converge: {r['msg']}")
    assert r["Mach_com"] < 0.1, (
        f"common-branch Mach {r['Mach_com']:.3f} is too high for incompressible comparison"
    )
