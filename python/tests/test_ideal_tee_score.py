"""The ideal tee's scoring against its pre-registration (#517), pinned.

``validation/junction/ideal_tee_score.py`` is the report; this holds its
numbers so a change to the element, the harness or the data shows up as a
diff rather than as a silently different report. Nothing here is a target:
the values are what the parameter-free model measured on Bassett 2001
(role: SELECTED), and a closure change that moves them should move these.

The metric is falsified alongside: the same scorer, handed the shipped tee,
must report a different bias -- a scorer that reports the same number for
both models is measuring the data, not the model.
"""

from __future__ import annotations

import pytest

from validation.junction import ideal_tee_score as its
from validation.junction.models.ideal_network import IdealMomentumNetwork
from validation.junction.models.mpce_network import MPCENetwork
from validation.junction.network_runner import run_network


@pytest.fixture(scope="module")
def records():
    return (
        run_network(IdealMomentumNetwork(), topologies=("imposed_q",)),
        run_network(MPCENetwork(), topologies=("imposed_q",)),
    )


@pytest.fixture(scope="module")
def summary(records):
    return its.dividing_report(*records)


def test_s1_straight_run_overpredicts_as_preregistered(summary) -> None:
    """S1: the kappa = 1 straight run sits above Bassett's measured K5 on
    every slice; +0.26 over the 43 points."""
    s1 = summary["K5"]
    assert s1["n"] == 43
    assert s1["bias_ideal"] == pytest.approx(0.264, abs=0.005)


def test_s2_branch_bias_is_reported(summary) -> None:
    """S2: size not predicted; measured +0.23 over 62 points."""
    s2 = summary["K6"]
    assert s2["n"] == 62
    assert s2["bias_ideal"] == pytest.approx(0.233, abs=0.005)


def test_the_scorer_tells_the_two_models_apart(summary) -> None:
    for K_id in ("K5", "K6"):
        assert abs(summary[K_id]["bias_ideal"] - summary[K_id]["bias_shipped"]) > 0.2


def test_m1_network_reproduces_the_scored_bm_theory_variant(records) -> None:
    assert its.merge_check(records[0]) <= its.M1_TOLERANCE
