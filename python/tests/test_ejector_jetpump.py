"""Validate the Phase 1B unchoked-primary jet-pump closure of the ejector model.

`jet_pump_entrainment_ratio` is the Keenan constant-area mixing solve for the
unchoked-primary (subsonic jet-pump) regime -- see
validation/ejector/OPERATING_REGIMES_DESIGN.md sec 3.3.

Ground truth is physics/literature, not the code's own output:
  * as Mach -> 0 the closure must reduce to the LOSSLESS incompressible
    jet-pump head-flow relation N(M, R) -- the same conservation laws in the
    constant-density limit (Sanger, NASA TN D-4445, "Noncavitating Performance
    of Two Low-Area-Ratio Water Jet Pumps", 1968, Appendix B). `N_lossless`
    below is that relation, derived from Bernoulli inlets + constant-area throat
    momentum + ideal diffuser; it is verified against an independent numeric
    integration to 1e-4 in the docstring's provenance record.
  * that lossless relation OVERPREDICTS Sanger's measured (with-loss) head by
    ~1.4x -- the data-anchored basis for recovery_efficiency ~= 0.7 (kept
    unfitted at 1.0 here, per the calibration policy).
  * on the reported degenerate network the closure must give a BOUNDED omega
    O(few), not the ~33 the double-choked construction wrongly returns.

R (Sanger area ratio) = An/At = A_p1/A_3 = area_ratio_nozzle / area_ratio_mix.
M (flow ratio) = omega at low Mach. N (head ratio) = (P_b - P_e)/(P_g - P_b).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import validation.ejector.data as _data_pkg
from combaero.network._ejector_huang1999 import (
    cd_nozzle_exit_static,
    cd_nozzle_mass_flow,
    choked_mass_flow,
)
from validation.ejector.models.huang1999 import (
    EjectorGeometry,
    jet_pump_entrainment_ratio,
)

_GAMMA = 1.40
_FIXTURE = json.loads((Path(_data_pkg.__file__).parent / "huang1999_reference.json").read_text())


def N_lossless(R: float, M: float) -> float:
    """Lossless incompressible jet-pump head ratio N(M, R) (Sanger D-4445).

    Derived from Bernoulli inlets + constant-area throat momentum + ideal
    diffuser; the head ratio is independent of the absolute driving pressure
    (incompressible), so it is a pure function of the area ratio R = A_p1/A_3
    and the flow ratio M."""
    g = 1.0 / (1.0 + M) ** 2
    A = g * (1.0 / R + M * M / (1.0 - R))
    num = A - 0.5 - 0.5 * g * M * M / (1.0 - R) ** 2
    den = 0.5 * g / R**2 - A + 0.5
    return num / den


def _geom_for_R(R: float) -> EjectorGeometry:
    # Only A_p1/A_3 = R matters to omega; pick A_p1/A_t = 1 so A_3/A_t = 1/R.
    return EjectorGeometry(area_ratio_nozzle=1.0, area_ratio_mix=1.0 / R)


def _feasible_window(p_g: float, p_e: float, geom: EjectorGeometry) -> tuple[float, float]:
    """Discharge back-pressure window [lower, upper] of interior solutions,
    found via the closure's own boundary flags (public API only)."""

    def boundary(p_b: float) -> str:
        return jet_pump_entrainment_ratio(p_g, 300.0, p_e, 300.0, geom, _GAMMA, p_b).boundary

    lo, hi = p_e, p_g  # upper edge: largest p_b that is not "back_pressure"
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if boundary(mid) == "back_pressure":
            hi = mid
        else:
            lo = mid
    upper = 0.5 * (lo + hi)
    lo, hi = p_e, p_g  # lower edge: smallest p_b that is not "primary_choke"
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if boundary(mid) == "primary_choke":
            lo = mid
        else:
            hi = mid
    lower = 0.5 * (lo + hi)
    return lower, upper


@pytest.mark.parametrize("drive", (0.02, 0.006, 0.002, 0.0006))
def test_reduces_to_lossless_closed_form_as_mach_falls(drive: float) -> None:
    """The compressible closure's (M, N) point converges to the lossless
    incompressible relation N(M, R) as the driving pressure (Mach) falls.
    Error must shrink with drive and be < 1.5% even at the largest drive."""
    R = 0.125
    geom = _geom_for_R(R)
    p_g, p_e = 1e5 * (1.0 + drive), 1e5 * (1.0 - 0.2 * drive)
    lower, upper = _feasible_window(p_g, p_e, geom)
    p_b = lower + 0.4 * (upper - lower)  # interior
    r = jet_pump_entrainment_ratio(p_g, 300.0, p_e, 300.0, geom, _GAMMA, p_b)
    assert r.boundary == ""
    M = r.omega
    N = (p_b - p_e) / (p_g - p_b)
    rel_err = abs(N - N_lossless(R, M)) / N_lossless(R, M)
    assert rel_err < 0.015
    # tighter bound at low Mach (the reduction is exact in the limit)
    if drive <= 0.002:
        assert rel_err < 0.002


def test_low_mach_error_decreases_monotonically() -> None:
    """The lossless-reduction error must shrink monotonically toward zero as
    the driving pressure falls -- the signature of an exact limit."""
    R = 0.125
    geom = _geom_for_R(R)
    errs = []
    for drive in (0.02, 0.006, 0.002, 0.0006, 0.0002):
        p_g, p_e = 1e5 * (1.0 + drive), 1e5 * (1.0 - 0.2 * drive)
        lower, upper = _feasible_window(p_g, p_e, geom)
        p_b = lower + 0.4 * (upper - lower)
        r = jet_pump_entrainment_ratio(p_g, 300.0, p_e, 300.0, geom, _GAMMA, p_b)
        N = (p_b - p_e) / (p_g - p_b)
        errs.append(abs(N - N_lossless(R, r.omega)) / N_lossless(R, r.omega))
    assert all(errs[i + 1] < errs[i] for i in range(len(errs) - 1))
    assert errs[-1] < 1e-3


def test_lossless_overpredicts_sanger_measured_head() -> None:
    """The lossless relation overpredicts Sanger's MEASURED water-jet-pump head
    by ~1.4x (real pump ~70% of ideal) -- the data-anchored evidence that
    recovery ~= 0.7, kept unfitted at 1.0 here. Guards the audit finding."""
    for R, M, N_measured in ((0.066, 3.5, 0.084), (0.197, 1.4, 0.255)):
        ideal = N_lossless(R, M)
        assert ideal > N_measured  # lossless overpredicts
        assert 1.3 < ideal / N_measured < 1.6  # by ~1.4x


def test_reported_case_omega_is_bounded() -> None:
    """The reported degenerate network (secondary ~= outlet ~= primary) must
    yield a bounded omega O(few) across the whole feasible window -- NOT the
    ~33 the double-choked construction returns (OPERATING_REGIMES_DESIGN 2.3)."""
    geom = EjectorGeometry(1.0e-4 / 3.14e-5, 8.0e-4 / 3.14e-5)  # A_p1/A_3 = 0.125
    p_g, p_e = 102435.0, 101325.0
    lower, upper = _feasible_window(p_g, p_e, geom)
    for frac in (0.1, 0.3, 0.5, 0.7, 0.9):
        p_b = lower + frac * (upper - lower)
        r = jet_pump_entrainment_ratio(p_g, 300.0, p_e, 288.15, geom, _GAMMA, p_b)
        assert r.boundary == ""
        assert 0.0 < r.omega < 10.0  # bounded and forward, not 33


def test_omega_increases_as_back_pressure_falls() -> None:
    """More suction (lower discharge pressure) entrains more flow: omega must
    increase monotonically as p_b decreases through the feasible window."""
    geom = _geom_for_R(0.125)
    p_g, p_e = 1.02e5, 0.997e5
    lower, upper = _feasible_window(p_g, p_e, geom)
    omegas = []
    for frac in (0.9, 0.7, 0.5, 0.3, 0.1):  # decreasing p_b
        p_b = lower + frac * (upper - lower)
        omegas.append(jet_pump_entrainment_ratio(p_g, 300.0, p_e, 300.0, geom, _GAMMA, p_b).omega)
    assert all(omegas[i + 1] > omegas[i] for i in range(len(omegas) - 1))


def test_boundary_flags_and_forward_direction() -> None:
    """Above the window: no forward entrainment (back_pressure, omega=0). Below:
    primary chokes (primary_choke). Interior: forward pumping, p_e < p_b."""
    geom = _geom_for_R(0.125)
    p_g, p_e = 1.02e5, 0.997e5
    lower, upper = _feasible_window(p_g, p_e, geom)
    high = jet_pump_entrainment_ratio(
        p_g, 300.0, p_e, 300.0, geom, _GAMMA, upper + 0.5 * (p_g - upper)
    )
    assert high.boundary == "back_pressure" and high.omega == 0.0
    low = jet_pump_entrainment_ratio(
        p_g, 300.0, p_e, 300.0, geom, _GAMMA, lower - 0.5 * (lower - p_e)
    )
    assert low.boundary == "primary_choke"
    mid = jet_pump_entrainment_ratio(p_g, 300.0, p_e, 300.0, geom, _GAMMA, 0.5 * (lower + upper))
    assert mid.boundary == "" and mid.omega > 0.0
    assert p_e < 0.5 * (lower + upper)  # forward pumping raises the secondary


def test_fixture_jetpump_sweep_reproduces_model() -> None:
    """The stored air jet-pump sweep (omega at each P_b) reproduces the live
    closure -- the regression lock on the Phase 2/3 C++ golden target."""
    jp = _FIXTURE["jetpump_cases"][0]
    i = jp["inputs"]
    geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
    for pt in jp["sweep"]:
        live = jet_pump_entrainment_ratio(
            i["p_g_pa"],
            i["t_g_k"],
            i["p_e_pa"],
            i["t_e_k"],
            geom,
            i["gamma"],
            pt["p_b_pa"],
            i["recovery_efficiency"],
        ).omega
        assert live == pytest.approx(pt["omega"], rel=1e-9)


def test_cd_nozzle_saturates_at_choked_cap() -> None:
    """The C-D nozzle (R0 primary residual) saturates at the sonic-throat mass
    flow once the exit pressure is low enough to choke -- deep expansion gives
    the choked_mass_flow value (to within the smooth-min rounding)."""
    A_t, A_e, g, R, Tt, p0 = 3.14e-5, 1.0e-4, 1.40, 287.0, 300.0, 101325.0
    cap = choked_mass_flow(p0, Tt, A_t, g, R)
    deep = cd_nozzle_mass_flow(p0, Tt, 0.3 * p0, A_t, A_e, g, R)
    assert deep == pytest.approx(cap, rel=1e-3)  # within eps_frac rounding


def test_cd_nozzle_monotone_rises_as_exit_pressure_falls() -> None:
    """Unchoked: lower exit pressure -> more expansion -> more mass flow, up to
    the choke cap. Strictly increasing as p_static_down falls, then flat."""
    A_t, A_e, g, R, Tt, p0 = 3.14e-5, 1.0e-4, 1.40, 287.0, 300.0, 101325.0
    ps = [p0 * f for f in (0.999, 0.99, 0.97, 0.95, 0.9)]
    mdots = [cd_nozzle_mass_flow(p0, Tt, p, A_t, A_e, g, R) for p in ps]
    assert all(mdots[i + 1] >= mdots[i] for i in range(len(mdots) - 1))
    cap = choked_mass_flow(p0, Tt, A_t, g, R)
    assert all(m <= cap * (1.0 + 1e-9) for m in mdots)  # never exceeds the cap


def test_cd_nozzle_is_c1_across_the_choke_threshold() -> None:
    """The smooth-min rounds the choke corner: the finite-difference slope is
    continuous across the threshold (a bare min() would jump to zero)."""
    A_t, A_e, g, R, Tt, p0 = 3.14e-5, 1.0e-4, 1.40, 287.0, 300.0, 101325.0

    def slope(p):
        h = 1.0
        return (
            cd_nozzle_mass_flow(p0, Tt, p + h, A_t, A_e, g, R)
            - cd_nozzle_mass_flow(p0, Tt, p - h, A_t, A_e, g, R)
        ) / (2.0 * h)

    # scan the corner; consecutive slopes must not jump discontinuously
    ps = [p0 * f for f in (0.975, 0.973, 0.971, 0.969, 0.967)]
    slopes = [slope(p) for p in ps]
    for i in range(len(slopes) - 1):
        assert abs(slopes[i + 1] - slopes[i]) < 5.0e-7  # bounded curvature, no kink


def test_cd_nozzle_exit_static_round_trips() -> None:
    """cd_nozzle_exit_static inverts cd_nozzle_mass_flow: feeding the exit
    static back through the nozzle must recover the input mass flow (the
    derived mixing pressure P_py of the unchoked-primary element branch)."""
    A_t, A_e, g, R, Tt, p0 = 3.14e-5, 1.0e-4, 1.40, 287.0, 300.0, 101325.0
    cap = choked_mass_flow(p0, Tt, A_t, g, R)
    for mp in (0.2 * cap, 0.5 * cap, 0.8 * cap, 0.95 * cap):
        p_py = cd_nozzle_exit_static(mp, p0, Tt, A_t, A_e, g, R)
        assert cd_nozzle_mass_flow(p0, Tt, p_py, A_t, A_e, g, R) == pytest.approx(mp, rel=1e-9)
        assert p0 * (2.0 / 2.4) ** (1.4 / 0.4) < p_py < p0  # within (p0*r_crit, p0)


def test_cd_nozzle_exit_static_decreases_with_mass_flow() -> None:
    """More mass flow needs more expansion -> a LOWER exit static. Monotone
    decreasing in m_dot."""
    A_t, A_e, g, R, Tt, p0 = 3.14e-5, 1.0e-4, 1.40, 287.0, 300.0, 101325.0
    cap = choked_mass_flow(p0, Tt, A_t, g, R)
    ps = [cd_nozzle_exit_static(f * cap, p0, Tt, A_t, A_e, g, R) for f in (0.2, 0.4, 0.6, 0.8)]
    assert all(ps[i + 1] < ps[i] for i in range(len(ps) - 1))


def test_raises_when_primary_cannot_be_unchoked() -> None:
    """If the primary supply is strong enough to choke before the secondary can
    flow, this is an ejector, not a jet pump -- the closure must refuse it."""
    geom = _geom_for_R(0.125)
    with pytest.raises(ValueError, match="unchoked-primary window"):
        # p_g >> p_e: primary chokes well above p_e -> no jet-pump window.
        jet_pump_entrainment_ratio(5.0e5, 300.0, 1.0e5, 300.0, geom, _GAMMA, 1.5e5)
