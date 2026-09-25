"""The aggregate scorecard: one runner per series, one view of all of them.

Three runners exist because three physics families score differently. Before
this, `scorecard.py` knew only the rib runner, so the jet-array and orifice
series read as "not scored by any set" and there was no single place to see
whether cooling as a whole was healthy.
"""

from __future__ import annotations

import math

import pytest

from validation.cooling import jet_array_runner, orifice_runner, runner
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


# ---------------------------------------------------------------------------
# Lau 1990 -- the first rib source outside Han (#385)
# ---------------------------------------------------------------------------


def test_lau_independently_confirms_han_G(cells) -> None:
    """A different lab, rig and decade reproducing han_1988_orthogonal's
    heat-transfer roughness function.

    Han prints G = 3.7 (e+)^0.28, Lau 4.218 (e+)^0.257 -- different
    coefficient AND exponent, agreeing because the forms cross near e+ = 300.
    """
    c = next(c for c in cells if c.label == "lau1990/fig_table2_90deg_G")
    assert c.n_scored == 9
    assert c.rmse < 0.03, f"agreement degraded to {c.rmse:.4f}"
    assert abs(c.bias) < 0.02
    # Every point inside Lau's own stated +/-5.8% Stanton uncertainty -- a
    # MEASUREMENT band from the source, not a model-derived one.
    assert c.within == 1.0


def test_lau_R_is_committed_but_not_scored(cells) -> None:
    """runner.py's e+ path has no absolute-R branch. The series is committed
    so the ~13% disagreement is visible rather than absent, and reports '-'
    rather than a zero that would read as agreement."""
    import math

    c = next(c for c in cells if c.label == "lau1990/fig_table2_90deg_R")
    assert c.n == 9
    assert c.n_scored == 0
    assert math.isnan(c.rmse)


def test_lau_Gbar_independently_confirms_hans_1p2_factor(cells) -> None:
    """Han's `G_bar = 1.2 G` rests on a single printed label on figure 4.47.
    Lau's Eq. (8) Gbar is the same construction -- St_avg over two ribbed and
    two smooth walls -- so it is a like-for-like check by another lab.

    That Han's Nu(AV) really is that average is verified, not assumed: NASA
    CR-3837's appendix prints Nu(R), Nu(S) and Nu(AV) per run, and Nu(AV) is
    their two-way mean to 0.084% over 33 rows.

    An earlier version of this test asserted the opposite -- that tagging
    Lau's Gbar would apply a "Prandtl factor" -- by quoting the WITHDRAWN
    form of han_ribbed.md item 10. See #392.
    """
    c = next(c for c in cells if c.label == "lau1990/fig_table2_90deg_Gbar")
    assert c.n == 9
    assert c.n_scored == 9
    # Two labs' wall-averaging ratios differ by 3-4% (Lau 1.235-1.251 against
    # Han's constant 1.200), so this confirms the relationship without
    # reproducing it exactly. Pinned as a band, not a point.
    assert 0.02 < c.rmse < 0.06, f"RMSE {c.rmse:.3f} outside the recorded band"
    assert c.bias < 0.0, "Han's 1.2 should sit below Lau's ratio, not above"
    assert c.within >= 0.85


def test_gbar_path_stays_restricted_to_90_deg(dataset, cells) -> None:
    """The real guard on the G_bar path is the rib angle, not the source.
    G_BAR_OVER_G = 1.2 was established at 90 deg; applying it to an angled
    rib is wrong by about 4.5%, larger than the measurement scatter, so
    runner.py refuses off-90 series rather than scaling them.

    This pins the refusal behaviourally: an off-90 G_bar series may carry a
    `scores:` target, but it must come back unscored.
    """
    by_label = {c.label: c for c in cells}
    seen_on, seen_off = False, False
    for series in dataset:
        if series.y_axis != "G_bar" or series.scores is None:
            continue
        cell = by_label.get(series.label)
        if cell is None:
            continue
        if series.alpha_deg == runner.G_BAR_VALID_ALPHA:
            seen_on = True
        else:
            seen_off = True
            assert cell.n_scored == 0, (
                f"{series.label} at alpha={series.alpha_deg} was scored "
                "through the G_bar path; G_BAR_OVER_G only holds at 90 deg"
            )
    assert seen_on and seen_off, "fixture no longer covers both sides of the guard"
