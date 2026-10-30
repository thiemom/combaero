"""Jobson (1955) against the sharp-orifice data of his Fig. 5 (#854, #861).

Fig. 5a is Perry's (1949) air (n = 1.4), Fig. 5b Schiller's (1933)
superheated steam (n = 1.3), both digitised and user-confirmed on a scan
overlay (validation/cooling/data/jobson1955/metadata.yaml). Jobson fitted
nothing to either: his Cd 0.600 is slit theory's Ci = pi/(pi + 2) times a
hydraulics velocity coefficient Cv = 0.982. So this is CROSS-SOURCE ACCURACY
of Jobson's theory (docs/VALIDATION_POLICY.md), never fidelity.

The bands below are the MODEL'S ERROR ON EACH SERIES, widened so a harmless
numerical change does not trip them. When written: Perry bias -1.43 %, RMS
1.57 %, worst -2.52 %; Schiller bias +1.13 %, RMS 1.91 %, worst +5.07 %.
They are NOT a measurement band: Jobson states none for either test series.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero import _solver_tools
from validation.cooling import orifice_runner, scorecard
from validation.cooling.schema import DATA_ROOT, load_dataset, load_points

CI_SLIT = math.pi / (math.pi + 2.0)
CV = 0.982

# figure -> (n, points, |bias| <, RMS <, |worst| <)
SERIES = {
    "5a": (1.4, 15, 0.025, 0.025, 0.035),
    "5b": (1.3, 9, 0.025, 0.030, 0.065),
}


@pytest.fixture(scope="module")
def jobson_series():
    series = [s for s in load_dataset(DATA_ROOT) if s.source.name == "jobson1955"]
    assert sorted(s.figure for s in series) == sorted(SERIES)
    return {s.figure: s for s in series}


def _scores(pred, points) -> tuple[float, float, float]:
    errs = [pred(p.x) / p.y - 1.0 for p in points]
    bias = sum(errs) / len(errs)
    rms = math.sqrt(sum(e * e for e in errs) / len(errs))
    return bias, rms, max(errs, key=abs)


def _K_jobson(r: float, n: float) -> float:
    return (
        CV
        * cb.jobson_contraction_coefficient(r, CI_SLIT, n)
        * (orifice_runner.nozzle_flux_coefficient(r, n))
    )


def _K_nozzle(r: float, n: float) -> float:
    # The default compressible path's shape: Cd K_N(r), flat below r_c.
    return CV * CI_SLIT * orifice_runner.nozzle_flux_coefficient(r, n)


@pytest.mark.parametrize("figure", sorted(SERIES))
def test_series_is_cross_source_and_owned_by_the_orifice_runner(jobson_series, figure) -> None:
    s = jobson_series[figure]
    n, count, *_ = SERIES[figure]
    assert len(load_points(s)) == count
    assert float(s.geometry["n"]) == n
    assert orifice_runner.owns(s)
    assert not scorecard.rib_owned(s)
    assert scorecard.basis_of(s, s.scores) == "accuracy"


@pytest.mark.parametrize("figure", sorted(SERIES))
def test_jobson_scores_inside_its_band(jobson_series, figure) -> None:
    """Cross-source accuracy; the band is the model's error, not a
    measurement band (see the module docstring)."""
    s = jobson_series[figure]
    n, count, max_bias, max_rms, max_worst = SERIES[figure]
    bias, rms, worst = _scores(lambda r: _K_jobson(r, n), load_points(s))
    assert abs(bias) < max_bias
    assert rms < max_rms
    assert abs(worst) < max_worst
    # The runner scores the same thing the scorecard reports.
    sc = orifice_runner.score(orifice_runner.run_series(s))
    assert sc.n == count
    assert sc.rms == pytest.approx(rms, rel=1e-12)


@pytest.mark.parametrize("figure", sorted(SERIES))
def test_jobson_beats_the_default_nozzle(jobson_series, figure) -> None:
    """The nozzle's flow is flat below r_c where a sharp orifice keeps
    passing more, so no Cd fixes its shape (#854). Same Cd 0.600 for both."""
    n = SERIES[figure][0]
    pts = load_points(jobson_series[figure])
    _, rms_j, _ = _scores(lambda r: _K_jobson(r, n), pts)
    _, rms_n, _ = _scores(lambda r: _K_nozzle(r, n), pts)
    assert rms_n > 0.10
    assert rms_n > 5.0 * rms_j


def test_network_kernel_agrees_with_the_paper_form_on_perry(jobson_series) -> None:
    """The solver's kernel (real-air gamma at 300 K, not exactly 1.4) gives
    the same K as the paper's form at every Perry pressure ratio (3e-5 when
    written). Schiller's n = 1.3 is Jobson's polytropic index for steam, not
    a real-gas gamma the kernel would compute, so it is not checked here."""
    X = cb.species.dry_air()
    Y = list(cb.mole_to_mass(X))
    Pt, Tt, A = 2.0e5, 300.0, 1e-4
    rho0 = cb.density(Tt, Pt, X)
    for p in load_points(jobson_series["5a"]):
        r = _solver_tools.orifice_jobson_residuals_and_jacobian(
            0.0, Pt, Pt, Tt, Y, p.x * Pt, CV * CI_SLIT, A, 0.0, CV
        )
        K = r.m_dot_calc / (A * math.sqrt(Pt * rho0))
        assert pytest.approx(_K_jobson(p.x, 1.4), rel=1e-3) == K, p.x
