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
    return build(run_dataset(dataset), dataset)


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
    rows = rollup(build(run_dataset(dataset), dataset))

    def only(set_name: str):
        """The one rollup row for a set.

        Rollup rows are segregated by sampling completeness, so a set with
        both completely- and partially-sampled series has several. The two
        cross-checked here have one apiece; assert that rather than assume
        it, so this stops being a valid comparison loudly if it changes.
        """
        matching = [c for c in rows if c.scored_by == set_name]
        assert len(matching) == 1, (
            f"{set_name} now spans {len(matching)} sampling classes; comparing "
            "it against a single whole-set number from the runner is no "
            "longer like-for-like"
        )
        return matching[0]

    flor = [s for s in dataset if s.source.name == "florschuetz1981"]
    recs = [r for r in jet_array_runner.run_all(flor) if r.predicted is not None]
    errs = [r.rel_error for r in recs]
    direct_bias = sum(errs) / len(errs)
    assert only("florschuetz_1981_inline").bias == pytest.approx(direct_bias, abs=1e-9)
    assert only("florschuetz_1981_inline").n_scored == len(errs)

    ms = [s for s in dataset if s.source.name == "mcgreehan_schotsch1988"]
    sc = orifice_runner.score(orifice_runner.run_all(ms))
    assert only("mcgreehan_schotsch_1988_crossflow_cd").bias == pytest.approx(sc.bias, abs=1e-9)


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


def test_lau_R_is_scored_now_that_an_absolute_R_branch_exists(cells) -> None:
    """This asserted the opposite until the e+ path gained an absolute-R
    branch: the series was committed `scores: null` so the ~13%
    disagreement stayed visible rather than absent, because scoring it
    would have compared R against G.

    The disagreement itself is pinned in
    test_lau_R_disagreement_is_now_scored_not_just_narrated."""
    c = next(c for c in cells if c.label == "lau1990/fig_table2_90deg_R")
    assert c.n == 9
    assert c.n_scored == 9
    assert not math.isnan(c.rmse)


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


def test_rollup_never_pools_across_sampling_completeness(dataset) -> None:
    """A tabulated series and an overplotted one do not measure the same
    thing, so the rollup must not average them into one row.

    Where a figure overplots symbol classes, only the spatially isolated
    marks can be digitised, and those are the ones furthest from the
    cluster centre. So a `partial` row's error is an upper bound, and a
    `complete` row's is an estimate. One number covering both is neither.
    """
    summary = rollup(build(run_dataset(dataset), dataset))
    han = {c.sampling: c for c in summary if c.scored_by == "han_1988_orthogonal"}
    assert "complete" in han and "partial" in han, (
        "han_1988_orthogonal must report completely- and partially-sampled "
        f"series separately; got {sorted(han)}"
    )
    assert han["complete"].label != han["partial"].label


def test_segregation_recovers_the_sources_own_stated_agreement(dataset) -> None:
    """The reason segregation matters, pinned as a number.

    Han prints that 95% of his data lies within 8% on G and 6% on R.
    Pooled, `han_1988_orthogonal` reported 73% within, which reads as a
    deficiency in the correlation. Segregated, the completely-sampled
    series reach about 91% -- close to the source's own claim -- and the
    partially-sampled ones about 65%.

    The gap is the pick, not the model. If this test fails because the
    complete subset fell, something real changed in the correlation; if
    it fails because the two converged, the sampling classification has
    stopped discriminating and is no longer worth its complexity.
    """
    # Restricted to G. The comparison has to be like-for-like in QUANTITY:
    # pooling R in drags the completely-sampled figure from 84% to 76%, not
    # because sampling got worse but because lau1990's R carries a genuine
    # ~11% disagreement between two labs. A real model disagreement in one
    # quantity would otherwise masquerade as a sampling effect -- the same
    # comparability trap #389 is about.
    cells_by_label = {c.label.split("  [")[0]: c for c in build(run_dataset(dataset), dataset)}
    buckets: dict[str, list[tuple[float, float, int]]] = {"complete": [], "partial": []}
    for series in dataset:
        cell = cells_by_label.get(series.label)
        if cell is None or cell.scored_by != "han_1988_orthogonal":
            continue
        if series.y_axis != "G" or not cell.n_scored:
            continue
        if cell.sampling in buckets:
            buckets[cell.sampling].append((cell.within, cell.mae, cell.n_scored))

    def weighted(rows, index):
        n = sum(k for *_, k in rows)
        return sum(row[index] * row[2] for row in rows) / n

    assert buckets["complete"] and buckets["partial"]
    complete_within = weighted(buckets["complete"], 0)
    partial_within = weighted(buckets["partial"], 0)

    assert complete_within > 0.80, (
        f"completely-sampled G series now agree only {complete_within:.0%} of "
        "the time, against the source's stated 95%"
    )
    assert partial_within < 0.70
    assert complete_within - partial_within > 0.15, (
        "sampling completeness no longer separates the two populations "
        f"({complete_within:.0%} against {partial_within:.0%})"
    )
    # And the upper-bound reading must hold: partial error is the larger.
    assert weighted(buckets["partial"], 1) > weighted(buckets["complete"], 1)


def test_recovered_bounds_are_model_free_and_discarded_when_loose(dataset) -> None:
    """Recovery bounds a missing mark by the OTHER marks near that
    abscissa -- never by the correlation plus its stated band, which
    would bound the observation with the model being scored.

    The price is that the envelope is only sometimes tight enough to
    constrain, and that is a property of the figure. Figure 4.46's
    classes hug one correlation, so its envelopes are usable; figure
    4.51's classes genuinely separate, so most of its are not and are
    discarded rather than counted as weak evidence.
    """
    from validation.cooling.recovery import MAX_WIDTH, recover

    by_figure: dict[str, int] = {}
    for series in dataset:
        bounds = recover(series, dataset)
        for bound in bounds:
            assert bound.lo < bound.hi
            assert bound.width <= MAX_WIDTH, (
                f"{series.label} kept a {bound.width:.0%} envelope; anything "
                f"over {MAX_WIDTH:.0%} admits nearly any prediction"
            )
            assert bound.n_marks >= 2
        if bounds:
            by_figure[series.figure] = by_figure.get(series.figure, 0) + len(bounds)

    assert by_figure.get("4.46", 0) >= 20, (
        f"figure 4.46 should recover a useful number of runs; got {by_figure.get('4.46', 0)}"
    )
    assert by_figure.get("4.51", 0) <= 3, (
        "figure 4.51's classes separate, so its envelopes should almost all "
        f"be discarded; kept {by_figure.get('4.51', 0)}"
    )


def test_absolute_R_series_are_scored_against_R_not_G(cells) -> None:
    """Sources that TABULATE the roughness functions print absolute `R`,
    not the normalised ordinate figure 4.46 plots. The e+ path had no
    branch for it, so six series sat committed with `scores: null` purely
    for want of one -- the measurements were never in doubt.

    Pinned as a band rather than a point: these are measurements, and the
    numbers should move if the correlation does.
    """
    for label in (
        "han_park_lei1984/table_a90_R",
        "han_park_lei1984/table_a45_R",
        "lau1990/fig_table2_90deg_R",
    ):
        cell = next(c for c in cells if c.label == label)
        assert cell.n_scored == cell.n > 0, f"{label} is still unscored"
        assert not math.isnan(cell.rmse)


def test_eq_4_17_is_confirmed_out_of_sample_on_R(cells) -> None:
    """CR-3837 is e/D = 0.063, a blockage Han and Park never measured, so
    scoring its R against Eq. 4.17 is a genuine out-of-sample test of the
    friction roughness function.

    It passes: inside the source's own 6.6% friction uncertainty at 90 and
    45 deg, and the 90 deg set has every point within band.
    """
    ninety = next(c for c in cells if c.label == "han_park_lei1984/table_a90_R")
    assert ninety.mae < 0.066, (
        f"Eq. 4.17 now misses CR-3837's 90 deg R by {ninety.mae:.1%}, outside "
        "the source's stated 6.6% friction uncertainty"
    )
    assert ninety.within == 1.0


def test_lau_R_disagreement_is_now_scored_not_just_narrated(cells) -> None:
    """Reviewer item L2. Lau's R runs well above `han_1988_orthogonal`'s
    constant 3.2, which converts to two labs ~12% apart in ribbed-wall
    friction for nominally the same 90 deg configuration.

    Until the absolute-R branch existed this lived only in prose. It is now
    a harness number, and the point of pinning it is that it must stay
    VISIBLE -- a silent drift to agreement would mean the correlation or the
    conversion moved, not that two labs reconciled.
    """
    cell = next(c for c in cells if c.label == "lau1990/fig_table2_90deg_R")
    assert cell.bias < -0.08, (
        f"Lau's R bias is now {cell.bias:.1%}; the recorded disagreement is "
        "about -11%, so something moved"
    )
    # Lau states +/-10.9% on friction, and the gap is larger than his own band.
    assert cell.within == 0.0


def test_unknown_y_axis_on_the_eplus_path_is_refused(dataset) -> None:
    """The branch used to end in `else: predicted = g`, so any y_axis that
    was not G_bar or R_normalised -- absolute R among them -- was silently
    scored against G. Nothing hit it only because every such series carried
    `scores: null`.

    A quantity the path cannot predict must come back unscored with a
    reason, never compared against a different quantity.
    """
    import dataclasses

    from validation.cooling.runner import run_series

    series = next(s for s in dataset if s.label == "han_park_lei1984/table_a90_R" and s.scores)
    impostor = dataclasses.replace(series, y_axis="Nu_ratio")
    records = run_series(impostor)
    assert records
    for record in records:
        assert record.predicted is None, (
            "an unpredictable y_axis was given a number; the catch-all is back"
        )
        assert "Nu_ratio" in (record.reason or "")


def test_fidelity_and_accuracy_are_never_pooled(dataset) -> None:
    """Fidelity (the correlation's own paper) and accuracy (cross-source)
    answer different questions and fail differently: a fidelity miss is our
    transcription, an accuracy miss is the model's limitation. One number
    covering both is neither.

    See docs/VALIDATION_POLICY.md.
    """
    rows = rollup(build(run_dataset(dataset), dataset))
    for set_name in ("han_1988_orthogonal", "han_park_1988_angled"):
        bases = {r.basis for r in rows if r.scored_by == set_name}
        assert {"fidelity", "accuracy"} <= bases, (
            f"{set_name} reports {sorted(bases)}; both must appear as separate rows"
        )
    for row in rows:
        assert row.basis in ("fidelity", "accuracy", "unknown")


def test_every_scored_set_declares_which_paper_it_is(dataset) -> None:
    """`basis_of` returns "unknown" for a set with no recorded origin rather
    than defaulting to fidelity, so a new set cannot quietly claim the
    stronger of the two. This asserts none is currently unknown -- i.e. the
    map has kept up with the sets."""
    from validation.cooling.scorecard import SET_ORIGIN

    scored = {c.scored_by for c in build(run_dataset(dataset), dataset) if c.n_scored}
    missing = {s for s in scored if s not in SET_ORIGIN}
    assert not missing, (
        f"these sets score series but declare no origin paper: {sorted(missing)}. "
        "Add them to SET_ORIGIN or their results cannot be read as fidelity "
        "or accuracy."
    )


def test_han_1988_fidelity_rests_on_very_few_clean_points(dataset) -> None:
    """The interaction the pooled number hid, pinned so it stays visible.

    `han_1988_orthogonal`'s own data survives only as an overplotted
    scatter, so almost all its FIDELITY evidence is sampling-biased and its
    cross-source ACCURACY looks better than its fidelity purely because the
    cross-source data is tabulated.

    If this starts failing because the clean fidelity set grew, that is
    good news and the number should be updated.
    """
    rows = rollup(build(run_dataset(dataset), dataset))
    clean = [
        r
        for r in rows
        if r.scored_by == "han_1988_orthogonal"
        and r.basis == "fidelity"
        and r.sampling == "complete"
    ]
    assert clean, "han_1988_orthogonal reports no cleanly-sampled own-data row"
    assert sum(r.n_scored for r in clean) < 30, (
        "the clean fidelity set has grown; update docs/VALIDATION_POLICY.md's "
        "worked example, which says it rests on eight points"
    )
