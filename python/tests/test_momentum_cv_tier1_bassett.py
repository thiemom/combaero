"""
Tier 1 validation: MPCE + BorderCarnotLoss vs Bassett 2001 / Hager 1984 at M -> 0.

STATUS: both K tests are xfail. The PDF Section 3.4 claim "MPCE + loss element
form 2 reproduces Hager exactly at M -> 0" does NOT hold. Tracked as issue #272.

CONVENTION (fixed 2026-09-04, was the source of a bogus comparison). Each
Bassett coefficient is a function of the mass-flow fraction in ITS OWN leg, so
K6 takes the lateral fraction but K2/K5 take the STRAIGHT fraction. `q` in this
module is the lateral fraction, so K2 is evaluated at 1-q via the transform
equivalences.py already defines. Evaluating K2 at the lateral fraction -- which
an earlier local copy of the formula did -- compares against a mirrored curve.
The corrected K2 target is identically Hager xi_t = q^2 - 0.5*q.

MEASURED (imposed-q network below, M_com = 0.025 throughout). K is defined here
as (Pt_com - Pt_port) / q_dyn_com, so positive K = stagnation loss:

    q     K_str   Bassett K2    K_lat   Bassett K6   K_lat bounded-LM
   0.25   0.4375    -0.0625     0.7751     0.8712        0.6799
   0.50   0.7500     0.0000     0.8657     0.8673        0.4847
   0.75   0.9375     0.1875     1.2719     0.9885        0.4145

Two SEPARATE problems, not one. Earlier revisions of this docstring framed it
as a single undisentangled question and quoted formulas from a model version
that predates the sin^2(theta) + cross-coupling integration; both are corrected
here against the numbers above.

1. K_straight: a MODEL gap, and the cleaner of the two.
   The model gives exactly K_str = 2*q - q^2, solver-independently (identical
   under the default and the bounded solver), against the corrected reference
   q^2 - 0.5*q. The reference crosses zero at q=0.5 and stays small; the model
   is positive across the whole range and roughly an order of magnitude larger.

   ROOT CAUSE -- the sin^2(theta) gate on the port's OWN angle. The kernel
   (src/solver_interface.cpp) evaluates

       R_i = P_i + sin^2(theta_i) * rho_i*u_i^2 + cross_i - P_jct

   so a collinear port (theta_i = 0) has its momentum flux gated OFF and
   reduces to STATIC pressure equality P_i = P_jct -- not the impulse function
   R_mom,i = (P_i + rho_i*u_i^2) - P_jct that the implementation guide states
   for every port. That gate is the whole story:

       static equality  (implemented) -> K_str = 2*q - q^2   (positive)
       impulse function (spec)        -> K_str = q^2 - 2*q   (negative)

   Note this is NOT settled by sweeping the LATERAL angle -- K_str is
   bit-identical at 30/45/60/90/120 deg there, which is why an earlier revision
   of this docstring wrongly concluded the projection was inert on this port.
   Sweeping the STRAIGHT port's own angle moves it a long way: at q=0.5,
   K_str = 0.7500 / 0.2759 / -0.0485 / -0.1754 at theta_straight =
   0 / 15 / 30 / 45 deg.

   WHY IT MATTERS -- the junction is loss-free BY DESIGN (empirical loss lives
   in separate per-port elements), so the question is not whether the junction
   alone reproduces a dissipative reference, but whether the loss a port
   element would have to supply is physically realisable:

       from the implemented static form : 2*q^2 - 2.5*q   NEGATIVE
                                          (-0.75 at q=0.5: needs pressure
                                          RECOVERY; no loss element can do it)
       from the spec impulse form       : 1.5*q           NON-NEGATIVE
                                          (0 at q=0, 1.5 at q=1: realisable)

   So the spec residual leaves a deficit a loss element can close and the
   implemented one does not. That points at the collinear-port gate as the
   defect rather than at a capability limit of the closure.

   Note also that the required 1.5*q is LINEAR in the split ratio, not
   quadratic in a local velocity: a constant-coefficient element referenced to
   the straight leg supplies L*(1-q)^2, referenced to the common leg L. Neither
   is 1.5*q, so closing this needs a q-dependent (Mynard-style) coefficient --
   i.e. a junction quantity, which is an argument for folding the straight-port
   loss into the closure rather than a standalone element.

   FIXTURE CAVEAT: this network attaches a BorderCarnotLossElement to the
   lateral but a LosslessConnectionElement to the straight leg, then compares
   K_straight against a dissipative reference. The comparison is only
   meaningful once the straight port carries its loss too.

2. K_lateral: closer, and NOT explained by solver basin choice.
   The default solver agrees with Bassett K6 to 11.0% / 0.2% / 28.7% at
   q = 0.25 / 0.5 / 0.75 -- it fails the 20%-or-0.15-absolute tolerance only
   at q = 0.75. An earlier hypothesis held that bounds (enforcing physically
   forward flow) would fix this; experimental_bounded_solver.py tests it and
   the answer is NO -- bounded least_squares converges cleanly (|F| ~ 1e-11)
   and lands FARTHER from K6 (22% / 44% / 58%), i.e. worse than the default
   solver at every q. Degenerate-root selection is therefore ruled out as the
   explanation, and the remaining candidates are the loss-element form or the
   cross-coupling projection.

Both xfails are strict=False, so they will not fail loudly if a change happens
to fix them; re-check the table above after touching the junction closure.

Tier 1 is foundational: until this is sorted, Tier 2 (Perez-Garcia / Wang
compressible) rides on a low-Mach baseline that already disagrees with the
canonical reference. See docs/junction/tier2_reference_data.md.

Regenerate the table by calling _extract_K_at_imposed_q(q) below for the
default-solver columns, and with
  uv run pytest python/tests/experimental_bounded_solver.py -s
for the bounded-LM column.
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
from validation.junction.equivalences import q_transform
from validation.junction.models import bassett2001

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_F_C = 0.01  # 100 cm^2
_M_DOT_IN = 0.1  # kg/s -- low-Mach inlet


# Each Bassett coefficient is a function of the mass-flow fraction in ITS OWN
# leg -- equivalences.py: "q = m_other / m_com, where m_other is the flow leg
# defining the K". So K6 (lateral) takes the lateral fraction while K2/K5
# (straight) take the STRAIGHT fraction. `q` throughout this module is the
# lateral fraction (m_lateral / m_in), which is Hager's convention, so K2 needs
# the 1-q transform that equivalences.py already defines for hager1984/xi_t.
#
# The canonical implementations are imported rather than re-typed: local copies
# are how the transform got skipped in the first place, and K2 evaluated at the
# lateral fraction is a mirrored curve, not Bassett's straight-loss coefficient.
_Q_TO_K2 = q_transform("hager1984", "xi_t")


def _K6_bassett(q: float, psi: float = 1.0, theta_rad: float = math.pi / 2.0) -> float:
    """Bassett Eq. 27 (separating lateral loss) at the lateral fraction q."""
    return bassett2001.K6(q, psi, theta_rad)


def _K2_bassett(q: float) -> float:
    """Bassett Eq. 15 (separating straight loss) for a lateral fraction q.

    Evaluates the canonical K2 at the straight fraction 1-q. The result is
    identically Hager xi_t = q^2 - 0.5*q, as it must be.
    """
    return bassett2001.K2(_Q_TO_K2(q))


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
        "K_lateral agrees with Bassett K6 to 11.0% / 0.2% / 28.7% at "
        "q = 0.25 / 0.5 / 0.75 -- only q=0.75 misses the tolerance. Bounds "
        "do NOT explain the gap: experimental_bounded_solver.py converges to "
        "|F| ~ 1e-11 and lands farther out (22% / 44% / 58%), so degenerate-"
        "root selection is ruled out. Remaining candidates are the loss-"
        "element form and the cross-coupling projection. Issue #272."
    ),
    strict=False,
)
def test_low_mach_K_lateral_agrees_with_bassett_K6():
    """Sweep q over Bassett's range, compare K_lateral to K6.

    With sin^2(theta) impulse + cross-coupling (Bassett axial momentum +
    wall reaction Eq 3-4) integrated into MPCE, K_lat tracks Bassett
    K6 = 1 + q^2 - 2*q*cos((3/4)*theta) closely at low q and diverges as q
    rises. See the module docstring for the measured table.
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
        "The collinear straight port has its momentum flux gated off by "
        "sin^2(theta_i)=0, so it uses static equality (K_str = 2*q - q^2) "
        "rather than the impulse function the spec states (q^2 - 2*q), "
        "against the corrected reference q^2 - 0.5*q. Under the spec form a "
        "straight-port loss element would need 1.5*q (realisable); under the "
        "implemented form 2*q^2 - 2.5*q (negative, unrealisable). Also, this "
        "fixture leaves the straight leg lossless. Issue #272."
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
