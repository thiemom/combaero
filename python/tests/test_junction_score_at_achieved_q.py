"""The junction scorecard must say WHERE it measured, and refuse to score a fixture.

Two defects behind these tests (issue #271, record in
docs/archive/JUNCTION_OPERATING_POINT_271.md):

* Only ``imposed_q`` imposes the split. In the pressure-driven topologies the
  solve settles wherever the model's coefficients reproduce the imposed
  pressure differences, which is generally NOT the target q -- and the
  scorecard compared the extracted K against the paper at the target anyway.
* ``three_pb`` imposes every total pressure through lossless connections, and
  ``_extract_K`` normalises by the fixed reference flow, so its extracted K is
  the imposed target by construction. Its K column measured the fixture.

The falsification test below is the load-bearing one: it perturbs the MODEL and
shows the ``three_pb`` K does not move, which is what makes the exclusion a
finding rather than an opinion.
"""

from __future__ import annotations

import math

import pytest

from validation.junction.models import bassett2001
from validation.junction.models.mpce_v2_network import MPCEv2Network
from validation.junction.network_runner import (
    _K_SCORING_TOPOLOGIES,
    _Q_OFF_POINT,
    NetworkRecord,
    _curve_value_at,
    build_network_cells,
    format_network_scorecard,
)

_THETA = math.radians(45.0)
_PSI = 3.0


@pytest.fixture
def model():
    return MPCEv2Network(strict=False)


# ---------------------------------------------------------------------------
# The tautology, falsified
# ---------------------------------------------------------------------------


def test_three_pb_K_is_the_imposed_target_whatever_the_model_does():
    """A deliberately wrong model must still reproduce three_pb's target.

    eta_scale=3.0 triples Mynard's energy-transfer factor. It moves the model's
    own answer (imposed_q) by more than 6%. If three_pb's K moved with it, the
    topology would be scoring the model and the exclusion below would be wrong.
    """
    honest = MPCEv2Network(strict=False)
    wrong = MPCEv2Network(strict=False, eta_scale=3.0)
    q = 0.8
    target = bassett2001.K6(q, _PSI, _THETA)

    honest_imposed = honest.evaluate_network("bassett2001", "K6", q, _PSI, _THETA)
    wrong_imposed = wrong.evaluate_network("bassett2001", "K6", q, _PSI, _THETA)
    assert honest_imposed.converged and wrong_imposed.converged
    moved = abs(wrong_imposed.K_lateral - honest_imposed.K_lateral)
    assert moved > 0.15, f"the perturbation must actually change the model: moved only {moved:.4f}"

    wrong_three_pb = wrong.evaluate_network(
        "bassett2001", "K6", q, _PSI, _THETA, topology="three_pb"
    )
    assert wrong_three_pb.converged
    assert wrong_three_pb.K_lateral == pytest.approx(target, abs=0.01), (
        "three_pb no longer reproduces the target from a wrong model -- if this "
        "fails, the topology has become informative and the exclusion in "
        "_K_SCORING_TOPOLOGIES should be revisited, not deleted"
    )


def test_three_pb_is_excluded_from_K_scoring():
    assert "three_pb" not in _K_SCORING_TOPOLOGIES
    assert {"imposed_q", "mfb_two_pb"} <= _K_SCORING_TOPOLOGIES


def _record(
    topology: str,
    q: float,
    q_conv: float | None,
    K_ext: float,
    K_at_conv: float | None = 1.0,
) -> NetworkRecord:
    return NetworkRecord(
        paper="bassett2001",
        K_id="K6",
        canonical_K="K_lateral_sep",
        psi=1.0,
        theta_deg=45.0,
        q=q,
        K_measured=1.0,
        K_extracted=K_ext,
        converged=True,
        residual_norm=1e-12,
        wall_time_s=0.01,
        which="lateral",
        topology=topology,
        q_converged=q_conv,
        K_measured_at_q_converged=K_at_conv,
    )


def test_excluded_topology_reports_no_rmse_and_says_why():
    cells = build_network_cells([_record("three_pb", 0.5, 0.9, 5.0)])

    assert math.isnan(cells[0].rmse_meas)
    assert math.isnan(cells[0].bias_meas)
    text = format_network_scorecard("m", cells)
    assert "taut" in text, "an excluded K must not render as a blank or a dash"


def test_scored_topology_still_reports_rmse():
    cells = build_network_cells([_record("mfb_two_pb", 0.5, 0.5, 1.5)])

    assert cells[0].rmse_meas == pytest.approx(0.5)
    assert cells[0].bias_meas == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# The achieved operating point
# ---------------------------------------------------------------------------


def test_drift_is_reported_per_cell():
    cells = build_network_cells(
        [_record("mfb_two_pb", 0.5, 0.9, 1.5), _record("mfb_two_pb", 0.5, 0.7, 1.5)]
    )

    assert cells[0].median_q_drift == pytest.approx(0.3)


def test_a_record_that_did_not_converge_has_no_drift():
    r = _record("mfb_two_pb", 0.5, None, 1.0)
    assert r.q_drift is None


@pytest.mark.parametrize("q", [0.4, 0.6, 0.8])
def test_imposed_q_reaches_the_operating_point_it_was_given(model, q):
    """The split is a boundary condition there, so the instrument must read
    back exactly what was imposed. This is the calibration of the measurement."""
    r = model.evaluate_network("bassett2001", "K6", q, _PSI, _THETA, topology="imposed_q")

    assert r.converged, r.message
    assert r.q_converged == pytest.approx(q, abs=1e-6)


def test_pressure_driven_solve_reports_where_it_actually_landed(model):
    """Bassett Fig 7b mfb_two_pb q=0.8 settles at q=0.857, not 0.8.

    Reported as a mirror root before this field existed, because the K was
    compared against the paper at 0.8 while the solve sat elsewhere.
    """
    r = model.evaluate_network("bassett2001", "K6", 0.8, _PSI, _THETA, topology="mfb_two_pb")

    assert r.converged, r.message
    assert r.q_converged == pytest.approx(0.857, abs=0.005)
    assert r.q_converged != pytest.approx(0.8, abs=0.01)


def test_bassett_K11_drift_is_reported_on_the_papers_own_axis(model):
    """K11 is indexed on the STRAIGHT inlet, so the adapter feeds the network
    1 - q. The achieved q has to come back through the same transform or the
    drift is measured on a mirrored axis -- the defect class of #276/#277/#279.
    """
    q = 0.3
    r = model.evaluate_network("bassett2001", "K11", q, 1.0, _THETA, topology="imposed_q")

    assert r.converged, r.message
    # imposed_q pins the split, so on the file's axis the achieved q is q.
    assert r.q_converged == pytest.approx(q, abs=1e-6)


def test_idelchik_K11_is_not_axis_flipped(model):
    """Idelchik tabulates every diagram on the lateral fraction, so K11 there
    passes q through unchanged. Applying Bassett's 1 - q would mirror it."""
    q = 0.3
    r = model.evaluate_network("idelchik1966", "K11", q, 3.333, _THETA, topology="imposed_q")

    assert r.converged, r.message
    assert r.q_converged == pytest.approx(q, abs=1e-6)


# ---------------------------------------------------------------------------
# Item 9b: the comparison happens ON THE CURVE, where the solve landed
# ---------------------------------------------------------------------------

_CURVE = [(0.0, 0.0), (0.5, 1.0), (1.0, 3.0)]


@pytest.mark.parametrize(
    "q, expected",
    [(0.0, 0.0), (0.25, 0.5), (0.5, 1.0), (0.75, 2.0), (1.0, 3.0)],
)
def test_curve_is_interpolated_between_its_digitised_points(q, expected):
    assert _curve_value_at(_CURVE, q) == pytest.approx(expected)


@pytest.mark.parametrize("q", [-0.01, 1.01, 5.0])
def test_outside_the_digitised_range_there_is_no_ground_truth(q):
    """Extrapolating a digitised curve is not measured data (policy sec 0
    item 1). Such a record is reported as off-curve, never scored."""
    assert _curve_value_at(_CURVE, q) is None


def test_an_empty_curve_yields_nothing():
    assert _curve_value_at([], 0.5) is None


def test_error_is_measured_where_the_solve_landed_not_where_it_was_aimed():
    r = _record("mfb_two_pb", q=0.2, q_conv=0.9, K_ext=3.0, K_at_conv=2.5)

    assert r.error == pytest.approx(0.5)
    assert r.K_measured == 1.0  # the value at the target q, deliberately unused


def test_a_record_off_the_curve_is_not_scored():
    r = _record("mfb_two_pb", q=0.2, q_conv=1.5, K_ext=3.0, K_at_conv=None)

    assert r.off_curve
    assert r.error is None
    cells = build_network_cells([r])
    assert cells[0].n_off_curve == 1
    assert math.isnan(cells[0].rmse_meas), "an unscorable record must not enter RMSE"


def test_a_record_that_landed_elsewhere_is_scored_but_counted_as_lost_coverage():
    """It is a correct measurement of the model at 0.9, and no evidence at all
    about 0.2. Both facts have to survive into the scorecard."""
    r = _record("mfb_two_pb", q=0.2, q_conv=0.9, K_ext=3.0, K_at_conv=2.5)

    assert r.off_point
    cells = build_network_cells([r])
    assert cells[0].n_off_point == 1
    assert cells[0].rmse_meas == pytest.approx(0.5), "it must still be scored"
    assert "off-point" in format_network_scorecard("m", cells)


def test_a_record_within_the_off_point_threshold_keeps_its_coverage():
    r = _record("mfb_two_pb", q=0.5, q_conv=0.5 + 0.5 * _Q_OFF_POINT, K_ext=1.0)
    assert not r.off_point


def test_imposed_q_is_scored_against_the_measured_value_itself(model):
    """The split is imposed, so nothing drifted and no interpolation should
    happen. An earlier cut compared floats exactly, and float noise pushed 9
    imposed_q records off the end of their curves and out of the score.
    """
    from validation.junction.network_runner import iter_network_records
    from validation.junction.schema import load_dataset

    records = [
        r
        for r in iter_network_records(model, load_dataset(), topologies=("imposed_q",))
        if r.converged and r.K_extracted is not None
    ]

    assert records, "no converged imposed_q records to check"
    assert all(r.K_measured_at_q_converged == r.K_measured for r in records)
    assert not any(r.off_curve for r in records)
    assert not any(r.off_point for r in records)
