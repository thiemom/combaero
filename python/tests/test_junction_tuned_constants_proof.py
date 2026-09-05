"""The on/off proofs for the tuned constants, pinned against the data.

``test_junction_tuned_constants.py`` pins the *values*. This file pins the
*reason they are allowed to have those values*: each term must improve
agreement with the digitised data on the cells it acts on, and must not act
anywhere else. If a change to the closure or the data flips one of these,
the constant's label on issue #271 is no longer true and it needs a new
table before it can stay.

Each assertion is a direction with a margin, not a stored number, so a small
re-digitisation does not trip it while a real regression does. Runs the
imposed_q scorecard, ~4 s per configuration.
"""

from __future__ import annotations

import statistics

import pytest

from validation.junction.models.mpce_v2_network import MPCEv2Network
from validation.junction.network_runner import run_network


def _mae_bias(alpha: float, eta: float, paper: str, K_ids: set[str]) -> tuple[float, float, int]:
    recs = [
        r
        for r in run_network(
            MPCEv2Network(joining_etransfer_alpha=alpha, eta_scale=eta),
            topologies=("imposed_q",),
        )
        if r.paper == paper and r.K_id in K_ids and r.converged and r.K_extracted is not None
    ]
    errs = [r.K_extracted - r.K_measured for r in recs]
    return statistics.fmean(abs(e) for e in errs), statistics.fmean(errs), len(errs)


@pytest.fixture(scope="module")
def hager():
    return {eta: _mae_bias(0.2, eta, "hager1984", {"xi_t"}) for eta in (0.0, 1.0)}


@pytest.fixture(scope="module")
def bassett_dividing():
    return {eta: _mae_bias(0.2, eta, "bassett2001", {"K5", "K6"}) for eta in (0.0, 1.0)}


@pytest.fixture(scope="module")
def idelchik_by_alpha():
    return {a: _mae_bias(a, 1.0, "idelchik1966", {"K11", "K12"}) for a in (0.0, 0.2, 0.3)}


def test_eta_improves_the_independent_straight_flow_data(hager):
    """Hager xi_t is the one K_straight source that is not Bassett.

    With eta off the closure over-predicts by a third of a dynamic head;
    with it on the bias falls by ~0.28 and MAE by ~0.05.
    """
    mae_off, bias_off, n_off = hager[0.0]
    mae_on, bias_on, n_on = hager[1.0]

    assert n_on == n_off == 45
    assert mae_on < mae_off - 0.03
    assert abs(bias_on) < abs(bias_off) - 0.2


def test_eta_improves_bassett_dividing_flow(bassett_dividing):
    mae_off, bias_off, _ = bassett_dividing[0.0]
    mae_on, bias_on, _ = bassett_dividing[1.0]

    assert mae_on < mae_off - 0.01
    assert abs(bias_on) < abs(bias_off)


def test_eta_does_not_touch_joining_cells():
    """Its (1 - lambda) factor is zero for a single collector, so joining
    coefficients must be identical with the term on and off -- a difference
    here is a plumbing defect, not a physics result."""
    off = _mae_bias(0.2, 0.0, "idelchik1966", {"K11", "K12"})
    on = _mae_bias(0.2, 1.0, "idelchik1966", {"K11", "K12"})

    assert on[2] == off[2]
    assert on[0] == pytest.approx(off[0], abs=2e-3)


def test_alpha_improves_joining_flow_over_off(idelchik_by_alpha):
    """alpha = 0 leaves Idelchik K11 under-predicted by ~0.3; 0.2 removes it."""
    mae0, bias0, n0 = idelchik_by_alpha[0.0]
    mae2, bias2, n2 = idelchik_by_alpha[0.2]

    assert n0 == n2
    assert mae2 < mae0 - 0.1
    assert abs(bias2) < abs(bias0) - 0.2


def test_alpha_default_beats_the_anchor_refit_in_network(idelchik_by_alpha):
    """The corrected-axis anchor refit (#283) proposed ~0.3; the measured and
    tabulated data prefer 0.2. Pins that the in-network table, not the
    analytical anchors, decides the value."""
    mae2, _, _ = idelchik_by_alpha[0.2]
    mae3, _, _ = idelchik_by_alpha[0.3]

    assert mae2 < mae3 - 0.02


def test_alpha_does_not_touch_dividing_cells():
    """It needs two suppliers; a dividing tee has one."""
    a0 = _mae_bias(0.0, 1.0, "hager1984", {"xi_t"})
    a2 = _mae_bias(0.2, 1.0, "hager1984", {"xi_t"})

    assert a0 == pytest.approx(a2)
