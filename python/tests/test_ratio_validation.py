"""Taslim & Spring (1987) ratio sets scored on the paper's own data (#444).

Each set's coefficients come from figure 9 (Nu) and figure 11 (f); the
scored series are figures 4 and 5. Same paper, different figures, so this is
fidelity: a miss is a disagreement between the paper's own figures or our
digitisation of them, never a model limitation.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

import combaero as cb
from validation.cooling import jet_array_runner, orifice_runner, ratio_runner
from validation.cooling.schema import load_dataset


@pytest.fixture(scope="module")
def dataset():
    return load_dataset()


@pytest.fixture(scope="module")
def taslim(dataset):
    return [s for s in dataset if ratio_runner.owns(s)]


def _mae(records) -> float:
    errs = [abs(r.rel_error) for r in records if r.rel_error is not None]
    return sum(errs) / len(errs)


def test_every_shipped_set_scores_one_series(taslim) -> None:
    assert len(taslim) == len(ratio_runner.SETS) == 7
    assert {s.scores for s in taslim} == set(ratio_runner.SETS)


def test_ownership_is_disjoint_from_the_other_specialists(dataset) -> None:
    for s in dataset:
        owners = [m.owns(s) for m in (ratio_runner, jet_array_runner, orifice_runner)]
        assert sum(owners) <= 1, s.label


def test_every_point_is_scored_inside_its_fitted_range(taslim) -> None:
    for s in taslim:
        records = ratio_runner.run_series(s)
        assert records, s.label
        for r in records:
            assert r.reason is None, (s.label, r.reason)
            assert r.predicted is not None
            assert not r.extrapolated, (s.label, r.x)


def test_fidelity_per_configuration(taslim) -> None:
    """Measured 2.6-10.7%; the two worst are figure 9 vs figure 4/5 offsets.

    The bound is the measurement plus margin, not a target: a regression that
    moves any configuration past it is a change to look at.
    """
    for s in taslim:
        assert _mae(ratio_runner.run_series(s)) < 0.12, s.label
    pooled = [r for s in taslim for r in ratio_runner.run_series(s)]
    assert _mae(pooled) < 0.08


def test_a_wrong_reynolds_exponent_scores_worse(taslim, monkeypatch) -> None:
    """Falsify the metric: drop the paper's Re^-0.2 and the fit must degrade.

    With Nu/Nu0 held at C, Nu goes as Re^0.8 instead of the paper's Re^0.6;
    over a factor of 3-6 in Re that is 25-40% at the ends.
    """
    pooled = [r for s in taslim for r in ratio_runner.run_series(s)]
    baseline = _mae(pooled)

    for name, factory in list(ratio_runner.SETS.items()):

        def mutated(factory=factory) -> cb.RibRatioSet:
            s = factory()
            s.Nu_Re = cb.RibTerm(0.0, s.Nu_Re.reference)
            return s

        monkeypatch.setitem(ratio_runner.SETS, name, mutated)

    worse = [r for s in taslim for r in ratio_runner.run_series(s)]
    assert _mae(worse) > baseline + 0.05


def test_one_side_series_are_refused_not_scored(dataset) -> None:
    s = next(s for s in dataset if ratio_runner.owns(s))
    one_side = replace(s, geometry={**s.geometry, "turbulated_walls": 1})
    records = ratio_runner.run_series(one_side)
    assert all(r.predicted is None and r.reason for r in records)
