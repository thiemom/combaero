"""Jobson (1955) against Perry's (1949) sharp-orifice air data (#854).

Perry's points are Jobson's Fig. 5a, digitised and user-confirmed on a scan
overlay (validation/cooling/data/jobson1955/metadata.yaml). Jobson fitted
nothing to them: his Cd 0.600 is slit theory's Ci = pi/(pi + 2) times a
hydraulics velocity coefficient Cv = 0.982. So this is CROSS-SOURCE ACCURACY
of Jobson's theory (docs/VALIDATION_POLICY.md), never fidelity.

The bands below are the MODEL'S ERROR ON THIS SERIES (bias -1.43 %, RMS
1.57 %, worst -2.52 % when written), widened so a harmless numerical change
does not trip them. They are NOT a measurement band: Jobson states none for
Perry's tests.
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
N = 1.4


@pytest.fixture(scope="module")
def perry():
    series = [s for s in load_dataset(DATA_ROOT) if s.source.name == "jobson1955"]
    assert len(series) == 1
    return series[0]


def _scores(pred, points) -> tuple[float, float, float]:
    errs = [pred(p.x) / p.y - 1.0 for p in points]
    bias = sum(errs) / len(errs)
    rms = math.sqrt(sum(e * e for e in errs) / len(errs))
    return bias, rms, max(errs, key=abs)


def _K_jobson(r: float) -> float:
    return (
        CV
        * cb.jobson_contraction_coefficient(r, CI_SLIT, N)
        * (orifice_runner.nozzle_flux_coefficient(r, N))
    )


def _K_nozzle(r: float) -> float:
    # The default compressible path's shape: Cd K_N(r), flat below r_c.
    return CV * CI_SLIT * orifice_runner.nozzle_flux_coefficient(r, N)


def test_series_is_cross_source_and_owned_by_the_orifice_runner(perry) -> None:
    assert len(load_points(perry)) == 15
    assert orifice_runner.owns(perry)
    assert not scorecard.rib_owned(perry)
    assert scorecard.basis_of(perry, perry.scores) == "accuracy"


def test_jobson_scores_inside_its_band_on_perry(perry) -> None:
    """Cross-source accuracy; the band is the model's error, not a
    measurement band (see the module docstring)."""
    bias, rms, worst = _scores(_K_jobson, load_points(perry))
    assert abs(bias) < 0.025
    assert rms < 0.025
    assert abs(worst) < 0.035
    # The runner scores the same thing the scorecard reports.
    s = orifice_runner.score(orifice_runner.run_series(perry))
    assert s.n == 15
    assert s.rms == pytest.approx(rms, rel=1e-12)


def test_jobson_beats_the_default_nozzle_on_perry(perry) -> None:
    """The nozzle's flow is flat below r_c where Perry's orifice keeps
    passing more, so no Cd fixes its shape (#854). Same Cd 0.600 for both."""
    pts = load_points(perry)
    _, rms_j, _ = _scores(_K_jobson, pts)
    _, rms_n, _ = _scores(_K_nozzle, pts)
    assert rms_n > 0.10
    assert rms_n > 5.0 * rms_j


def test_network_kernel_agrees_with_the_paper_form_on_perry(perry) -> None:
    """The solver's kernel (real-air gamma at 300 K, not exactly 1.4) gives
    the same K as the paper's form at every Perry pressure ratio (3e-5 when
    written)."""
    X = cb.species.dry_air()
    Y = list(cb.mole_to_mass(X))
    Pt, Tt, A = 2.0e5, 300.0, 1e-4
    rho0 = cb.density(Tt, Pt, X)
    for p in load_points(perry):
        r = _solver_tools.orifice_jobson_residuals_and_jacobian(
            0.0, Pt, Pt, Tt, Y, p.x * Pt, CV * CI_SLIT, A, 0.0, CV
        )
        K = r.m_dot_calc / (A * math.sqrt(Pt * rho0))
        assert pytest.approx(_K_jobson(p.x), rel=1e-3) == K, p.x
