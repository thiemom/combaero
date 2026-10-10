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

from validation.junction.models.mpce_network import MPCENetwork
from validation.junction.network_runner import run_network


def _errors(alpha: float, eta: float, paper: str, K_ids: set[str]) -> dict[tuple, float]:
    """K_extracted - K_measured per scored record, keyed by its cell and point.

    K_measured is in the key because a source can carry two points at one
    (psi, theta, q) -- Hager does."""
    return {
        (r.K_id, r.psi, r.theta_deg, r.q, r.K_measured): r.K_extracted - r.K_measured
        for r in run_network(
            MPCENetwork(joining_etransfer_alpha=alpha, eta_scale=eta),
            topologies=("imposed_q",),
        )
        if r.paper == paper and r.K_id in K_ids and r.converged and r.K_extracted is not None
    }


def _mae_bias(alpha: float, eta: float, paper: str, K_ids: set[str]) -> tuple[float, float, int]:
    errs = list(_errors(alpha, eta, paper, K_ids).values())
    return statistics.fmean(abs(e) for e in errs), statistics.fmean(errs), len(errs)


def _common(a: dict[tuple, float], b: dict[tuple, float]) -> tuple[dict, dict]:
    """The two scored sets on the records both converged.

    A record converging under one setting and not the other is allowed only at
    a curve endpoint (q = 0 or 1). There one port carries zero flow, its
    residual has a kink where the flow sign flips, and whether a solve that
    stalls next to the kink clears the verifier's m_dot > 1e-6 threshold is
    incidental (#272: Idelchik K11/K12, psi = 10, theta = 30, q = 0 accepted
    at eta = 1 with 1.3e-6 kg/s on the dead port, rejected at eta = 0). Every
    interior record must converge under both, or the comparison is not like
    for like.
    """
    differing = set(a) ^ set(b)
    assert all(k[3] in (0.0, 1.0) for k in differing), sorted(differing)
    assert len(differing) <= 2
    keys = set(a) & set(b)
    return {k: a[k] for k in keys}, {k: b[k] for k in keys}


def _stats(errs: dict) -> tuple[float, float, int]:
    e = list(errs.values())
    return statistics.fmean(abs(x) for x in e), statistics.fmean(e), len(e)


@pytest.fixture(scope="module")
def hager():
    return {eta: _mae_bias(0.2, eta, "hager1984", {"xi_t"}) for eta in (0.0, 1.0)}


@pytest.fixture(scope="module")
def bassett_dividing():
    return {eta: _mae_bias(0.2, eta, "bassett2001", {"K5", "K6"}) for eta in (0.0, 1.0)}


@pytest.fixture(scope="module")
def idelchik_by_alpha():
    return {a: _errors(a, 1.0, "idelchik1966", {"K11", "K12"}) for a in (0.0, 0.2, 0.3)}


def test_eta_no_longer_earns_its_place_on_the_independent_straight_data(hager):
    """Hager xi_t is the one K_straight source that is not Bassett.

    THIS ASSERTION IS INVERTED FROM ITS ORIGINAL FORM, and the reason is the
    matched pair. Mynard's eta was fitted to CFD in a formulation that had lost
    the dividing-streamline recovery on the continuing collector. With that
    recovery restored (`_mynard2010.DIVIDING_STREAMLINE_RECOVERY`) the two do
    the same job and eta is now a duplicate. Measured on the digitised
    dividing data, all four combinations:

        configuration          Hager xi_t MAE    bias   Bassett K5+K6 MAE   bias
        no term, eta=1 (was)          0.2859  +0.0546              0.1152 +0.0366
        term,    eta=1                0.3007  -0.2379              0.1205 -0.0868
        term,    eta=0 (now)          0.0859  +0.0358              0.0564 -0.0207
        no term, eta=0                0.3385  +0.3385              0.1569 +0.0874

    Neither half alone is good and both together are worse than the derived
    term alone, which is what a matched pair looks like when one half is
    standing in for missing physics.
    """
    mae_off, bias_off, n_off = hager[0.0]
    mae_on, _bias_on, _n_on = hager[1.0]

    assert n_off == 45
    assert mae_off < mae_on - 0.15, (
        "eta is supposed to be the worse option now; if this fails the term "
        "and the transfer are no longer duplicates and the default is worth "
        "re-deciding, not silently flipping"
    )
    assert abs(bias_off) < 0.05


def test_eta_no_longer_earns_its_place_on_bassett_dividing_flow(bassett_dividing):
    """Same inversion, on Bassett's own dividing pair. See the table above."""
    mae_off, bias_off, _ = bassett_dividing[0.0]
    mae_on, _bias_on, _ = bassett_dividing[1.0]

    assert mae_off < mae_on - 0.04
    assert abs(bias_off) < 0.05


def test_eta_does_not_touch_joining_cells():
    """Its (1 - lambda) factor is zero for a single collector, so joining
    coefficients must be identical with the term on and off -- a difference
    here is a plumbing defect, not a physics result."""
    off, on = _common(
        _errors(0.2, 0.0, "idelchik1966", {"K11", "K12"}),
        _errors(0.2, 1.0, "idelchik1966", {"K11", "K12"}),
    )

    # The MAE, not each record: the closure's joining K is eta-independent to
    # the last digit, but where each solve stops inside its tolerance still
    # moves single extracted K by up to ~0.015 at psi = 10 (also on main
    # before #272), and that averages out.
    assert _stats(on)[0] == pytest.approx(_stats(off)[0], abs=2e-3)


def test_alpha_improves_joining_flow_over_off(idelchik_by_alpha):
    """alpha = 0 leaves Idelchik K11 under-predicted by ~0.3; 0.2 removes it."""
    e0, e2 = _common(idelchik_by_alpha[0.0], idelchik_by_alpha[0.2])
    mae0, bias0, _ = _stats(e0)
    mae2, bias2, _ = _stats(e2)

    assert mae2 < mae0 - 0.1
    assert abs(bias2) < abs(bias0) - 0.2


def test_alpha_default_beats_the_anchor_refit_in_network(idelchik_by_alpha):
    """The corrected-axis anchor refit (#283) proposed ~0.3; the measured and
    tabulated data prefer 0.2. Pins that the in-network table, not the
    analytical anchors, decides the value."""
    e2, e3 = _common(idelchik_by_alpha[0.2], idelchik_by_alpha[0.3])
    mae2, _, _ = _stats(e2)
    mae3, _, _ = _stats(e3)

    assert mae2 < mae3 - 0.02


def test_alpha_does_not_touch_dividing_cells():
    """It needs two suppliers; a dividing tee has one."""
    a0 = _mae_bias(0.0, 1.0, "hager1984", {"xi_t"})
    a2 = _mae_bias(0.2, 1.0, "hager1984", {"xi_t"})

    assert a0 == pytest.approx(a2)
