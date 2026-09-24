"""The aggregate scorecard: one runner per series, one view of all of them.

Three runners exist because three physics families score differently. Before
this, `scorecard.py` knew only the rib runner, so the jet-array and orifice
series read as "not scored by any set" and there was no single place to see
whether cooling as a whole was healthy.
"""

from __future__ import annotations

import math

import pytest

from validation.cooling import jet_array_runner, orifice_runner
from validation.cooling.schema import load_dataset
from validation.cooling.scorecard import build, rollup, run_dataset


@pytest.fixture(scope="module")
def dataset():
    return load_dataset()


@pytest.fixture(scope="module")
def cells(dataset):
    return build(run_dataset(dataset))


def test_every_series_is_claimed_exactly_once(dataset) -> None:
    """The bug this guards against is double counting.

    Every runner returns reason-carrying records for series it cannot score,
    so "the runner returned something" is not an ownership test -- the first
    version of run_dataset used it and the jet-array runner swallowed the
    orifice series, reporting them unscored.
    """
    records = run_dataset(dataset)
    per_series: dict[str, set[str]] = {}
    for r in records:
        per_series.setdefault(r.series.label, set()).add(type(r).__module__)

    assert len(per_series) == len(dataset), "a series was dropped or duplicated"
    for label, modules in per_series.items():
        assert len(modules) == 1, f"{label} scored by {modules}"


def test_ownership_predicates_are_disjoint(dataset) -> None:
    """No series may be owned by two specialists."""
    for s in dataset:
        assert not (jet_array_runner.owns(s) and orifice_runner.owns(s)), s.label


def test_the_specialist_series_are_now_scored(cells) -> None:
    """Florschuetz and the orifice sources used to report 0 scored points."""
    scored = {c.label: c for c in cells if c.n_scored > 0}

    assert any(c.startswith("florschuetz1981/") for c in scored)
    assert any(c.startswith("rohde1969/") for c in scored)
    assert any(c.startswith("mcgreehan_schotsch1988/") for c in scored)


def test_rohde_is_reported_per_velocity_head_band(cells) -> None:
    """Pooling Rohde's bands would mix a near-exact comparison with one
    dominated by the total-to-static conversion (factor 1.01 vs 1.33), so the
    runner splits them and the scorecard must keep them split."""
    rohde = [c.label for c in cells if c.label.startswith("rohde1969/fig10_rd0_")]
    assert len(rohde) == 3, rohde
    assert all("[VHR " in lbl for lbl in rohde)

    # And the split is not cosmetic: the bands disagree by a wide margin.
    by_band = {
        c.label.split("[")[1]: c for c in cells if c.label.startswith("rohde1969/fig10_rd0_")
    }
    low = by_band["VHR 1-10]"].bias
    high = by_band["VHR 25-60]"].bias
    assert low > 3.0 * high, f"bands no longer separate: {low:.3f} vs {high:.3f}"


def test_rollup_covers_every_implemented_set(cells) -> None:
    sets = {c.scored_by for c in rollup(cells)}
    for expected in (
        "han_1988_orthogonal",
        "han_park_1988_angled",
        "rallabandi_2009_high_re",
        "florschuetz_1981_inline",
        "mcgreehan_schotsch_1988_cd",
        "mcgreehan_schotsch_1988_crossflow_cd",
    ):
        assert expected in sets, f"{expected} missing from the rollup"


def test_rollup_matches_the_runners_own_numbers(dataset) -> None:
    """The aggregate must not quietly disagree with the runner it aggregates.

    Scored independently here, through each runner's own API, and compared
    against what the rollup reports.
    """
    summary = {c.label: c for c in rollup(build(run_dataset(dataset)))}

    flor = [s for s in dataset if s.source.name == "florschuetz1981"]
    recs = [r for r in jet_array_runner.run_all(flor) if r.predicted is not None]
    errs = [r.rel_error for r in recs]
    direct_bias = sum(errs) / len(errs)
    assert summary["florschuetz_1981_inline"].bias == pytest.approx(direct_bias, abs=1e-9)
    assert summary["florschuetz_1981_inline"].n_scored == len(errs)

    ms = [s for s in dataset if s.source.name == "mcgreehan_schotsch1988"]
    sc = orifice_runner.score(orifice_runner.run_all(ms))
    assert summary["mcgreehan_schotsch_1988_crossflow_cd"].bias == pytest.approx(sc.bias, abs=1e-9)


def test_a_set_that_answered_nothing_cannot_look_perfect(cells) -> None:
    """The shape of mistake this harness exists to catch: an unscored series
    reporting 0.0% error rather than '-'."""
    for c in cells:
        if c.n_scored == 0:
            assert math.isnan(c.mae) and math.isnan(c.rmse) and math.isnan(c.bias)


def test_within_is_blank_where_no_band_is_stated(cells) -> None:
    """`uncertainty` is the band agreement is judged against. Han and
    Florschuetz state one (3-8%); Rohde and McGreehan-Schotsch state none, so
    their `within` must read '-' rather than a structural 0%.
    """
    for c in cells:
        if c.label.startswith(("rohde1969/", "mcgreehan_schotsch1988/")):
            assert math.isnan(c.within), f"{c.label} reports a within fraction"
        if c.label.startswith("han2012/fig4.46_R_eD") and c.n_scored:
            assert not math.isnan(c.within)
