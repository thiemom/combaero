"""The expansion factor Y, checked against combaero's own isentropic solve.

Eq. (5)'s Y_n IS the isentropic nozzle expansion factor. combaero already
computes that mass flow exactly on the regime='compressible' path, through
``nozzle_flow`` with real-gas properties -- a completely independent code
path from this closed form. Agreement between them does two jobs:

  - it validates the Eq. (5) transcription against something that was not
    derived from the paper, and
  - it is the evidence that applying Y on the compressible path would
    correct for compressibility twice.

See validation/cooling/extractions/orifice_discharge_coefficient.md, decision
D10.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb

GAMMA = 1.4
T0 = 300.0
P0 = 4.0e5
AREA = 1.0e-5
CD = 0.80


@pytest.fixture(scope="module")
def air():
    return cb.standard_dry_air_composition()


def test_nozzle_form_reproduces_the_isentropic_solve(air) -> None:
    """Eq. (5) against nozzle_flow: different code path, real gas properties.

    The comparison stops at the critical pressure ratio, below which the
    isentropic form has no valid branch and nozzle_flow chokes.
    """
    rho1 = cb.density(T0, P0, air)
    worst = 0.0
    for s in (0.99, 0.97, 0.95, 0.92, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60):
        exact = cb.nozzle_flow(T0, P0, s * P0, CD * AREA, air).mdot
        incompressible = CD * AREA * math.sqrt(2.0 * rho1 * (P0 - s * P0))
        y_implied = exact / incompressible
        y_eq5 = cb.mcgreehan_schotsch_1988_expansion_nozzle(s, GAMMA)
        worst = max(worst, abs(y_eq5 - y_implied) / y_implied)
    assert worst < 5e-4, f"Eq. (5) vs nozzle_flow drifted to {worst:.2e}"


def test_orifice_form_is_deliberately_not_the_nozzle_form(air) -> None:
    """An orifice is not a nozzle, which is why the paper carries both.

    If these two ever agreed, the Cd-dependent blend between them would be
    doing nothing and one of the transcriptions would be wrong.
    """
    for s in (0.9, 0.8, 0.7, 0.6):
        y_o = cb.mcgreehan_schotsch_1988_expansion_orifice(s, GAMMA)
        y_n = cb.mcgreehan_schotsch_1988_expansion_nozzle(s, GAMMA)
        assert y_o > y_n, f"ordering wrong at S={s}"
    # They differ by up to ~17% at S = 0.6.
    y_o = cb.mcgreehan_schotsch_1988_expansion_orifice(0.6, GAMMA)
    y_n = cb.mcgreehan_schotsch_1988_expansion_nozzle(0.6, GAMMA)
    assert (y_o - y_n) / y_n == pytest.approx(0.166, abs=0.01)


def test_blend_endpoints_select_the_right_form() -> None:
    s = 0.7
    sharp = cb.mcgreehan_schotsch_1988_expansion_factor(0.70, s, GAMMA)
    rounded = cb.mcgreehan_schotsch_1988_expansion_factor(0.995, s, GAMMA)
    assert sharp == pytest.approx(cb.mcgreehan_schotsch_1988_expansion_orifice(s, GAMMA), abs=5e-3)
    assert rounded == pytest.approx(cb.mcgreehan_schotsch_1988_expansion_nozzle(s, GAMMA), abs=5e-3)


def test_hard_clamp_has_a_dead_derivative_and_the_default_does_not() -> None:
    """The reviewer's concern, pinned.

    eps=0 is the paper's exact ramp clamped hard: dY/dCd is exactly zero
    outside [0.82, 0.94]. The default saturates smoothly instead.
    """
    s, h = 0.7, 1e-6

    def d_dcd(cd, eps):
        return (
            cb.mcgreehan_schotsch_1988_expansion_factor(cd + h, s, GAMMA, eps)
            - cb.mcgreehan_schotsch_1988_expansion_factor(cd - h, s, GAMMA, eps)
        ) / (2 * h)

    assert d_dcd(0.78, 0.0) == 0.0
    assert d_dcd(0.99, 0.0) == 0.0

    for cd in (0.70, 0.78, 0.818, 0.822, 0.90, 0.936, 0.944, 0.97, 0.99):
        assert d_dcd(cd, 0.05) != 0.0, f"dead derivative at Cd = {cd}"


def test_y_saturates_at_choking() -> None:
    """Below the critical pressure ratio the isentropic form turns over and
    predicts decreasing flow. The mass-flux proxy must not follow it down."""
    s_crit = (2.0 / (GAMMA + 1.0)) ** (GAMMA / (GAMMA - 1.0))

    def flux(s):
        y = cb.mcgreehan_schotsch_1988_expansion_factor(0.99, s, GAMMA)
        return y * math.sqrt(max(1.0 - s, 0.0))

    at_crit = flux(s_crit)
    for s in (0.45, 0.35, 0.20, 0.05):
        assert flux(s) > 0.97 * at_crit, f"flow collapses below S* at {s}"


def test_applying_y_on_top_of_nozzle_flow_would_double_count(air) -> None:
    """Guards the architecture decision, not just the number.

    nozzle_flow already contains this factor. Multiplying it in again is a
    large, silent error -- this records how large, so the decision is not
    quietly reversed later.
    """
    s = 0.7
    exact = cb.nozzle_flow(T0, P0, s * P0, CD * AREA, air).mdot
    y = cb.mcgreehan_schotsch_1988_expansion_factor(CD, s, GAMMA)
    doubled = exact * y
    assert (exact - doubled) / exact > 0.08, (
        "double-counting Y should cost well over 8% here; if this shrinks, the "
        "compressible path has changed and D10 needs revisiting"
    )


# ---------------------------------------------------------------------------
# Standing smoothness guard over the whole surface
# ---------------------------------------------------------------------------
#
# Point assertions catch the hazards someone thought to look for. This scans
# instead, which is how the hard min(S, 1) was found after two other bounds
# had been smoothed deliberately.


def test_expansion_factor_surface_is_solver_smooth() -> None:
    """No floor, kink or divergence anywhere in Y, on either axis."""
    from validation.solver_smoothness import assert_smooth

    for cd in (0.70, 0.85, 0.90, 0.95, 0.99):
        assert_smooth(
            lambda s, c=cd: cb.mcgreehan_schotsch_1988_expansion_factor(c, s, GAMMA),
            0.05,
            1.20,
            label=f"Y vs pressure ratio at Cd={cd}",
        )
    for s in (0.6, 0.7, 0.85, 0.95):
        assert_smooth(
            lambda c, sv=s: cb.mcgreehan_schotsch_1988_expansion_factor(c, sv, GAMMA),
            0.60,
            1.05,
            label=f"Y vs Cd at S={s}",
        )


def test_the_scan_catches_the_papers_hard_clamp() -> None:
    """Falsifies the guard above: with eps=0 the same scan must go red, and
    report both the flat regions and the corners."""
    from validation.solver_smoothness import Hazard, scan

    found = scan(
        lambda c: cb.mcgreehan_schotsch_1988_expansion_factor(c, 0.7, GAMMA, 0.0),
        0.60,
        1.05,
    )
    hazards = {g.hazard for g in found}
    assert Hazard.FLOOR in hazards
    assert Hazard.KINK in hazards
    corners = sorted(g.x for g in found if g.hazard is Hazard.KINK)
    assert corners[0] == pytest.approx(0.82, abs=0.01)
    assert corners[-1] == pytest.approx(0.94, abs=0.01)


def test_cd_chain_has_exactly_the_two_hazards_tracked_in_383() -> None:
    """The Cd chain shipped in #381 has two known hazards and no others.

    They are named in ``allow`` rather than skipping the scan, so the rest of
    the surface stays guarded and a NEW hazard would still fail here. Remove
    the entries as #383 resolves them.
    """
    from validation.solver_smoothness import Hazard, assert_smooth, scan

    # Known hazard 1: dCd/dRe is exactly zero below the Re = 1e4 validity
    # floor, and the floor's edge is a kink. Two faces of one clamp.
    re_found = scan(
        lambda r: cb.mcgreehan_schotsch_1988_cd(r, 0.0, 1.0, 0.0),
        1.0e3,
        1.0e6,
        log=True,
    )
    assert {g.hazard for g in re_found} == {Hazard.FLOOR, Hazard.KINK}
    assert all(g.x < 1.1e4 for g in re_found), "a hazard away from the Re floor"

    # Above the floor the Reynolds dependence is clean.
    assert_smooth(
        lambda r: cb.mcgreehan_schotsch_1988_cd(r, 0.0, 1.0, 0.0),
        1.2e4,
        1.0e6,
        log=True,
        label="Cd vs Re above the validity floor",
    )

    # Known hazard 2: dCd/d(U1/Vi) is unbounded at zero crossflow -- a
    # DIVERGENCE, not a kink, which is why #383 cannot fix it by smoothing a
    # corner. It sits exactly on the domain edge.
    u_found = scan(lambda u: cb.mcgreehan_schotsch_1988_cd(3.2e4, 0.0, 1.0, u), 0.0, 3.0)
    assert {g.hazard for g in u_found} == {Hazard.DIVERGENCE}
    assert all(g.x == 0.0 for g in u_found), "a new hazard away from the origin"

    # Above the singular point the crossflow term is clean.
    assert_smooth(
        lambda u: cb.mcgreehan_schotsch_1988_cd(3.2e4, 0.0, 1.0, u),
        0.05,
        3.0,
        label="Cd vs U1/Vi away from the origin",
    )
