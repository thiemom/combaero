"""Validate the Phase 1A subcritical-droop closed forms of the ejector model.

Covers the operating-regime extension's dead-head anchor P_b0
(`dead_head_back_pressure`), the Tier-1 linear droop
(`subcritical_entrainment_ratio`), and the C1 smooth-min blend
(`blended_entrainment_ratio`) -- see
validation/ejector/OPERATING_REGIMES_DESIGN.md sec 2.2, 3.2.

Ground truth is physics/literature, not the code's own output:
  * P_b0 is Kracik & Dvorak's SAME mixing closure at omega = 0, so it must sit
    above P_c* (dead-head holds more back pressure than the critical point) and
    below the primary supply, and be independent of A_3/A_t (which only sets
    omega, and omega = 0 here).
  * omega -> 0 continuity: at the geometry where A_3/A_t = A_py/A_t the critical
    entrainment is exactly zero, so P_c* must equal P_b0 exactly.
  * Tier-1 droop is exact at both anchors and linear between them.
  * the smooth-min rounds min(omega_crit, omega_sub) by at most eps/2 and is C1
    (continuous derivative) across P_c*, unlike a bare min().
  * the near-linear droop assumption is cross-checked against measured AIR data,
    Henry et al. (2007) HEFAT2007 Fig. 5 (validation/ejector/data/henry2007_fig5.py).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import validation.ejector.data as _data_pkg
from validation.ejector.data.henry2007_fig5 import FIG5_ANCHORS, FIG5_OMEGA_VS_PB
from validation.ejector.models.huang1999 import (
    EjectorGeometry,
    blended_entrainment_ratio,
    critical_back_pressure,
    dead_head_back_pressure,
    entrainment_ratio,
    subcritical_entrainment_ratio,
)

_FIXTURE = json.loads((Path(_data_pkg.__file__).parent / "huang1999_reference.json").read_text())
_CASES = _FIXTURE["cases"]
_CASE = _CASES[0]  # Huang T3-EH, a healthy motive ejector (M_py > 1)
_RTOL = 1.0e-9


def _inputs(case: dict) -> tuple[float, float, float, float, EjectorGeometry, float, float]:
    i = case["inputs"]
    geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
    return i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom, i["gamma"], i["r_gas"]


def _anchors(case: dict) -> tuple[float, float, float]:
    """(omega_crit, P_c*, P_b0) for a fixture case."""
    p_g, t_g, p_e, t_e, geom, gamma, r = _inputs(case)
    omega = entrainment_ratio(p_g, t_g, p_e, t_e, geom, gamma).omega
    p_c = critical_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_c_pa
    p_b0 = dead_head_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_b0_pa
    return omega, p_c, p_b0


@pytest.mark.parametrize("case", _CASES, ids=lambda c: c["id"])
def test_dead_head_above_critical_and_bounded(case: dict) -> None:
    """P_b0 is the upper (omega -> 0) anchor: P_e < P_c* < P_b0 < P_g. The
    dead-head holds strictly more back pressure than the critical point, but
    still cannot recover the full primary supply."""
    p_g, _, p_e, _, _, _, _ = _inputs(case)
    _omega, p_c, p_b0 = _anchors(case)
    assert p_e < p_c < p_b0 < p_g


def test_dead_head_independent_of_area_ratio_mix() -> None:
    """P_b0 depends only on the y-y state (mach_py, set by A_p1/A_t and the
    pressures), NOT on A_3/A_t -- which enters only through omega, and omega
    is zero at dead-head. Changing A_3/A_t must leave P_b0 unchanged."""
    p_g, t_g, p_e, t_e, geom, gamma, r = _inputs(_CASE)
    base = dead_head_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_b0_pa
    for d in (-3.0, +5.0):
        other = EjectorGeometry(geom.area_ratio_nozzle, geom.area_ratio_mix + d)
        assert dead_head_back_pressure(
            p_g, t_g, p_e, t_e, other, gamma, r
        ).p_b0_pa == pytest.approx(base, rel=_RTOL)


def test_critical_equals_dead_head_at_zero_entrainment() -> None:
    """omega -> 0 continuity, exact. At A_3/A_t = A_py/A_t the entrained area is
    zero so omega_crit = 0, and the critical closure (gam = omega = 0) must
    land exactly on the dead-head closure (gam = 0)."""
    p_g, t_g, p_e, t_e, geom, gamma, r = _inputs(_CASE)
    area_py = entrainment_ratio(p_g, t_g, p_e, t_e, geom, gamma).area_ratio_primary_core
    zero_geom = EjectorGeometry(geom.area_ratio_nozzle, area_py)
    omega0 = entrainment_ratio(p_g, t_g, p_e, t_e, zero_geom, gamma).omega
    assert omega0 == pytest.approx(0.0, abs=1e-12)
    p_c = critical_back_pressure(p_g, t_g, p_e, t_e, zero_geom, gamma, r).p_c_pa
    p_b0 = dead_head_back_pressure(p_g, t_g, p_e, t_e, zero_geom, gamma, r).p_b0_pa
    assert p_c == pytest.approx(p_b0, rel=_RTOL)


def test_tier1_droop_anchors_midpoint_and_linear() -> None:
    """Tier-1 droop: exact at both anchors, half-value at the midpoint, and
    linear (zero second difference) between them."""
    omega, p_c, p_b0 = _anchors(_CASE)
    assert subcritical_entrainment_ratio(omega, p_c, p_b0, p_c) == pytest.approx(omega, rel=_RTOL)
    assert subcritical_entrainment_ratio(omega, p_c, p_b0, p_b0) == pytest.approx(0.0, abs=1e-12)
    mid = 0.5 * (p_c + p_b0)
    assert subcritical_entrainment_ratio(omega, p_c, p_b0, mid) == pytest.approx(
        0.5 * omega, rel=_RTOL
    )
    xs = [p_c + (p_b0 - p_c) * k / 10.0 for k in range(11)]
    ys = [subcritical_entrainment_ratio(omega, p_c, p_b0, x) for x in xs]
    assert all(ys[k + 1] <= ys[k] + 1e-15 for k in range(len(ys) - 1))  # monotone
    second = [ys[k + 1] - 2 * ys[k] + ys[k - 1] for k in range(1, len(ys) - 1)]
    assert max(abs(s) for s in second) < 1e-12 * omega  # linear


def test_tier1_raises_on_collapsed_window() -> None:
    """A collapsed droop window (P_b0 <= P_c*) is the degenerate unchoked-
    primary corner the jet-pump regime owns; the closed form must refuse it."""
    with pytest.raises(ValueError, match="P_b0 > P_c"):
        subcritical_entrainment_ratio(0.4, 120000.0, 100000.0, 110000.0)


def test_blended_reduces_to_plateau_and_droop() -> None:
    """Below P_c* the blend returns the plateau omega_crit; above it, the
    droop omega_sub -- each recovered to within the eps rounding away from the
    corner."""
    omega, p_c, p_b0 = _anchors(_CASE)
    eps = 1.0e-3 * omega
    # Well below P_c* (plateau side): omega_eff ~ omega_crit.
    below = blended_entrainment_ratio(omega, p_c, p_b0, p_c - 0.3 * (p_b0 - p_c))
    assert below == pytest.approx(omega, abs=eps)
    # Well above P_c* (droop side): omega_eff ~ omega_sub(P_b).
    pb = 0.5 * (p_c + p_b0)
    droop = blended_entrainment_ratio(omega, p_c, p_b0, pb)
    assert droop == pytest.approx(subcritical_entrainment_ratio(omega, p_c, p_b0, pb), abs=eps)


def test_blended_smoothmin_bound_and_kink() -> None:
    """|omega_eff - min(omega_crit, omega_sub)| <= eps/2 everywhere, with
    equality exactly at the corner P_b = P_c* (where omega_sub = omega_crit)."""
    omega, p_c, p_b0 = _anchors(_CASE)
    eps_frac = 1.0e-3
    eps = eps_frac * omega
    xs = [p_c + (p_b0 - p_c) * k / 200.0 for k in range(-40, 201)]
    for x in xs:
        raw = min(omega, subcritical_entrainment_ratio(omega, p_c, p_b0, x))
        blended = blended_entrainment_ratio(omega, p_c, p_b0, x, eps_frac)
        # smooth-min sits at or just below min, by at most eps/2.
        assert 0.0 <= raw - blended <= eps / 2 + 1e-12
    at_kink = blended_entrainment_ratio(omega, p_c, p_b0, p_c, eps_frac)
    assert at_kink == pytest.approx(omega - eps / 2, rel=1e-9)


def test_blended_is_c1_across_critical() -> None:
    """The blend is C1: its derivative is continuous across P_c*. A bare min()
    would jump from slope 0 (plateau) to -omega_crit/(P_b0-P_c*) (droop); the
    smooth-min's local slope just left and just right of the corner differ by
    only a small fraction of that jump."""
    omega, p_c, p_b0 = _anchors(_CASE)
    eps_frac = 1.0e-3
    eps_pb = eps_frac * omega * (p_b0 - p_c) / omega  # eps expressed in P_b units
    slope_jump = omega / (p_b0 - p_c)  # magnitude of the raw min() kink

    h = 1.0e-3 * eps_pb  # FD step, << the sampling offset

    def deriv(x: float) -> float:
        hi = blended_entrainment_ratio(omega, p_c, p_b0, x + h, eps_frac)
        lo = blended_entrainment_ratio(omega, p_c, p_b0, x - h, eps_frac)
        return (hi - lo) / (2.0 * h)

    # Sample just inside the rounded corner (<< eps): both one-sided slopes
    # approach the smooth mean -slope_jump/2, so their difference is a tiny
    # fraction of slope_jump. A bare min() would give 0 and -slope_jump here --
    # a full slope_jump discontinuity.
    small = 1.0e-2 * eps_pb
    assert abs(deriv(p_c + small) - deriv(p_c - small)) < 0.05 * slope_jump
    # At the corner the slope is the smooth mean of the two one-sided slopes.
    assert deriv(p_c) == pytest.approx(-0.5 * slope_jump, rel=0.05)


# --- Literature cross-check: Henry et al. (2007) HEFAT2007 Fig. 5 (air) ---


@pytest.mark.parametrize("p1_bar", (4, 5, 6))
def test_henry_fig5_droop_is_near_linear(p1_bar: int) -> None:
    """The measured air-ejector subcritical droop is NEAR-linear -- the core
    justification for the Tier-1 linear model. A least-squares line through the
    droop points (between the plateau knee and breakdown) fits with R^2 >= 0.93
    and residuals within 10% of the plateau omega:
        4 bar: R^2 0.992, max residual 3.8% | 5 bar: 0.989, 3.3% |
        6 bar: 0.939, 8.5% (mildly concave -- the documented Tier-2 trigger,
        OPERATING_REGIMES_DESIGN.md sec 3.2).
    This validates the linear SHAPE only; the chord's slope/offset (anchored at
    the model's P_c* and P_b0) needs absolute-pressure data (Akbarnejad 2025) to
    validate separately -- eyeballing the knee here would conflate the two."""
    pts = FIG5_OMEGA_VS_PB[p1_bar]
    anch = FIG5_ANCHORS[p1_bar]
    knee_x, plateau = anch["pb_knee"], anch["omega_plateau"]
    brk_x = pts[-1][0]  # last measured point (near breakdown)
    droop = [(x, y) for (x, y) in pts if knee_x < x < brk_x]
    assert len(droop) >= 3

    xs = [x for x, _ in droop]
    ys = [y for _, y in droop]
    n = len(droop)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in droop)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    intercept = (sy - slope * sx) / n
    residuals = [y - (slope * x + intercept) for x, y in droop]
    ybar = sy / n
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    ss_res = sum(r * r for r in residuals)
    r_squared = 1.0 - ss_res / ss_tot

    assert r_squared >= 0.93
    assert max(abs(r) for r in residuals) < 0.10 * plateau


@pytest.mark.parametrize("case", _CASES, ids=lambda c: c["id"])
def test_fixture_p_b0_reproduces_model(case: dict) -> None:
    """The golden P_b0 stored on every case still matches the live model --
    the same regression lock the critical-mode fields already have."""
    p_g, t_g, p_e, t_e, geom, gamma, r = _inputs(case)
    live = dead_head_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_b0_pa
    assert live == pytest.approx(case["reference"]["p_b0_pa"], rel=_RTOL)


def test_fixture_subcritical_sweep_reproduces_model() -> None:
    """The stored subcritical sweep (omega_eff and its d/dP_b) reproduces the
    live blended closure -- keeps the Phase 2/3 C++ golden target honest."""
    sub = _FIXTURE["subcritical_cases"][0]
    a = sub["anchors"]
    for pt in sub["sweep"]:
        live = blended_entrainment_ratio(
            a["omega_crit"], a["p_c_pa"], a["p_b0_pa"], pt["p_b_pa"], a["eps_frac"]
        )
        assert live == pytest.approx(pt["omega_eff"], rel=_RTOL)


def test_henry_fig5_physical_trends() -> None:
    """Cross-curve physics the model must also honour: raising the primary
    pressure LOWERS the plateau entrainment (stronger motive, less entrained
    per unit primary) but RAISES the breakdown back pressure (stronger motive
    holds more)."""
    plateaus = [FIG5_ANCHORS[p]["omega_plateau"] for p in (4, 5, 6)]
    breakdowns = [FIG5_OMEGA_VS_PB[p][-1][0] for p in (4, 5, 6)]
    assert plateaus[0] > plateaus[1] > plateaus[2]
    assert breakdowns[0] < breakdowns[1] < breakdowns[2]
