"""Scored validation of the cooling correlations against digitised data.

These are not self-consistency checks. Every expectation here is a number
read off a published figure, and the bands come from the source's own
stated accuracy -- not from what the implementation currently returns.
See validation/cooling/ and issue #333.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from validation.cooling.runner import Record, run_all, run_series
from validation.cooling.schema import load_dataset, load_points


@pytest.fixture(scope="module")
def dataset():
    return load_dataset()


@pytest.fixture(scope="module")
def records(dataset) -> list[Record]:
    return run_all(dataset)


def _series(dataset, stem: str):
    for s in dataset:
        if s.path.stem == stem:
            return s
    raise AssertionError(f"no series {stem}")


def test_every_series_declares_a_cross_check(dataset) -> None:
    """The harness rule: a set without a check that could fail is not data."""
    for s in dataset:
        assert s.cross_check.strip(), f"{s.label} has no cross_check"
        assert s.item and s.series, f"{s.label} is missing provenance"


def test_han_1988_reproduces_its_own_figure_454_line(dataset) -> None:
    """The implemented set against the drawn 90 deg correlation line.

    Han and Zhang (1992) state the 90 deg reference line IS Han (1988), so
    this is the set meeting a redrawing of itself in a later paper. The 6%
    band is Han's own stated accuracy for 95% of the data, not a tolerance
    fitted to what the code returns.
    """
    recs = run_series(_series(dataset, "fig4.54_90deg_continuous_Gbar_curve"))
    errs = [r.rel_error for r in recs if r.rel_error is not None]
    assert len(errs) == 7
    rmse = math.sqrt(sum(e * e for e in errs) / len(errs))
    assert rmse < 0.06, f"RMSE {rmse:.1%} exceeds Han's stated 6%"


def test_han_1988_against_figure_454_measured_points(dataset) -> None:
    """Measured scatter must be worse than the drawn line, and still in band.

    If the points scored BETTER than the line they were drawn through, the
    two series would have been swapped in the metadata.
    """
    curve = run_series(_series(dataset, "fig4.54_90deg_continuous_Gbar_curve"))
    data = run_series(_series(dataset, "fig4.54_90deg_continuous_Gbar_data"))

    def rmse(recs: list[Record]) -> float:
        errs = [r.rel_error for r in recs if r.rel_error is not None]
        return math.sqrt(sum(e * e for e in errs) / len(errs))

    assert rmse(data) > rmse(curve)
    assert rmse(data) < 0.08


def test_45deg_series_is_refused_not_scored(dataset) -> None:
    """A 90 deg-only set must decline 45 deg data rather than score it.

    Scoring it produced a 26.7% bias that reads as model error when it is
    really the wrong correlation entirely. The refusal is the feature.

    Built from a synthetic series (dataclasses.replace on a real one) rather
    than pinned to a specific filed CSV: which series happen to carry
    scores=han_1988_orthogonal at alpha=45 is metadata, and changes as
    understanding improves -- fig4.193c_old_correlation_G_curve was scored
    against han_1988_orthogonal via this exact refusal path once, then
    re-labelled scores=null once a better reason (drawn line, no
    angled-correlation set yet) was found. The mechanism under test should
    not depend on any one series still holding that combination.
    """
    import dataclasses

    series_45 = dataclasses.replace(
        next(s for s in dataset if s.y_axis == "G"),
        alpha_deg=45.0,
        scores="han_1988_orthogonal",
        geometry=None,
    )
    recs = run_series(series_45)
    assert all(r.predicted is None for r in recs)
    assert all("90-90 deg only" in (r.reason or "") for r in recs)


def test_figure_454_x_scale_is_the_only_consistent_decade(dataset) -> None:
    """Falsify the recorded x_scale: neighbouring decades must be far worse.

    The digitiser dropped the axis's `x 10^-2` multiplier. This pins the
    recovery so a future re-digitisation cannot silently reintroduce it.
    """
    series = _series(dataset, "fig4.54_45deg_vbroken_Gbar_data")
    assert series.x_scale == 1.0e4
    raw = [(p.x / series.x_scale, p.y) for p in load_points(series)]

    def rms_against_fit(scale: float) -> float:
        # The extraction's free fit for this series, G_bar = 3.242 e+^0.2882.
        errs = [3.242 * (x * scale) ** 0.2882 / y - 1.0 for x, y in raw]
        return math.sqrt(sum(e * e for e in errs) / len(errs))

    assert rms_against_fit(1.0e4) < 0.02
    for wrong in (1.0, 1.0e3, 1.0e5):
        assert rms_against_fit(wrong) > 0.40


def test_45deg_vbroken_outperforms_90deg_line(dataset) -> None:
    """Lower G is better. The text calls broken V ribs the best performers.

    Reads the two series against each other at a shared e+, which a sign or
    axis error in either digitisation would break.
    """
    v_broken = load_points(_series(dataset, "fig4.54_45deg_vbroken_Gbar_data"))
    line = load_points(_series(dataset, "fig4.54_90deg_continuous_Gbar_curve"))

    def fit(points):
        n = len(points)
        sx = sum(math.log(p.x) for p in points)
        sy = sum(math.log(p.y) for p in points)
        sxx = sum(math.log(p.x) ** 2 for p in points)
        sxy = sum(math.log(p.x) * math.log(p.y) for p in points)
        b = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        return math.exp((sy - b * sx) / n), b

    c_v, n_v = fit(v_broken)
    c_l, n_l = fit(line)
    at = 500.0
    ratio = (c_v * at**n_v) / (c_l * at**n_l)
    assert 0.70 < ratio < 0.78, f"expected ~0.739 at e+=500, got {ratio:.3f}"


def test_every_digitised_series_passes_its_figure_card() -> None:
    """The figure card and the points are two independent channels.

    The card records what the printed axes and equations say, read off the
    page without reference to where the digitiser put anything. A
    disagreement means one of the two is wrong, and which one is not
    decided here -- it is raised for a human.

    Verified against injected bugs: a dropped axis multiplier (the real
    figure 4.54 defect) fails x-span, two curves swapped fails
    printed-exponent, a double-picked mark fails count and distinct.
    """
    from validation.cooling.verify import check_all

    failures = [f for f in check_all() if not f.ok]
    assert not failures, "\n".join(f"{f.series}: {f.check}: {f.detail}" for f in failures)


def test_pooled_figure_446_validates_the_implemented_set(records) -> None:
    """han_1988_orthogonal against the figure it was extracted from.

    Pooled rather than per-class: symbols on figure 4.46 overlap and a
    mark cannot always be assigned to its class, but the set carries zero
    geometry exponents in both R and G, so every class must land on one
    curve and pooling removes a dependence on labels that cannot be fully
    trusted. See the note at the head of the han2012 metadata.

    The bands are Han's own stated accuracy, not fitted to the code.
    """
    from validation.cooling.scorecard import pool

    n_r, mae_r, _, bias_r, _ = pool(records, "han2012/fig4.46_R_eD")
    assert n_r >= 60, f"only {n_r} R points pooled"
    assert abs(bias_r) < 0.02, f"R bias {bias_r:.1%} against the printed 3.2"
    assert mae_r < 0.06

    n_g, mae_g, _, bias_g, _ = pool(records, "han2012/fig4.46_G_eD")
    assert n_g >= 40, f"only {n_g} G points pooled"
    assert abs(bias_g) < 0.06, f"G bias {bias_g:.1%}"


def test_stated_6pct_behaves_as_one_sigma_not_a_95pct_bound(records) -> None:
    """Han states 'within 6% for 95% of the data'. The data says otherwise.

    About 70% of the lower-panel cloud falls inside 6%, with a measured
    standard deviation near 6%. That is 1 sigma, not a 95% bound -- which
    is what harness tolerances here are set from. Pinned because the
    alternative reading would justify a band roughly twice as wide, and a
    band twice as wide is how a real error hides.
    """
    from validation.cooling.scorecard import pool

    _, _, rmse, _, within = pool(records, "han2012/fig4.46_R_eD")
    assert 0.04 < rmse < 0.08, f"spread {rmse:.1%} is not near Han's 6%"
    assert 0.55 < within < 0.85, (
        f"{within:.0%} inside 6% -- consistent with 1 sigma (~68%), not with the stated 95%"
    )


def test_disputed_class_labels_are_declared_and_surfaced(dataset) -> None:
    """A disputed class label must stay visible, not decay into a comment.

    The four fig4.51 45-deg series return G_bar/G below 1, which is not
    physically attainable; which panel duplicated the other cannot be
    settled from the digitised data and needs the primary paper, Han,
    Zhang and Lee (1991) ASME JHT 113, 590-598.

    fig4.46_G_eD0.047_pe10_wh2 was the fifth and is now resolved -- see
    test_fig446_wh2_class_is_pinned_to_its_source_runs.

    This pins three things: the series stays marked disputed, its
    cross_check still carries the specific question to ask of the book,
    and the verifier reports it on every run. Deleting any of those makes
    the uncertainty invisible, which is worse than the uncertainty.
    """
    from validation.cooling.verify import check_all

    disputed = [s for s in dataset if s.class_confidence == "disputed"]
    assert disputed, "the known disputed series has lost its marking"

    for s in disputed:
        # The route to resolution differs -- a figure may settle one, a
        # primary paper another -- but a dispute with no route recorded is
        # just an unexplained flag.
        assert "NEEDS THE" in s.cross_check, (
            f"{s.label} is disputed but records no route to resolve it"
        )

    reported = {f.series for f in check_all() if f.check == "class-label"}
    for s in disputed:
        assert s.label in reported, f"{s.label} is disputed but not surfaced"


def test_fig446_wh2_class_is_pinned_to_its_source_runs(dataset) -> None:
    """The e/D 0.047, P/e 10, W/H 2 class was disputed on the grounds that
    its two panels pair on only 2 of 4 marks. Resolved 2026-09-26 by going
    behind the figure to NASA CR-4015 appendix 7.3 (pages 141-145), which
    tabulates the five source runs with self-identifying headers, so the
    class is found by label and the overplotted cluster is never resolved.

    This pins the finding rather than the prose: every digitised mark in
    both panels must land on one of the five documented runs. If someone
    re-picks these series and a mark drifts off the run grid, that is a
    digitisation error and should fail here.
    """
    # CR-4015 p.141-145, E/D=0.047 P/E=10.0 ALPHA=90 HYD DIA=2.667 IN.
    # e+ from the printed R = 3.2 via f; the paper states R holds 95% of
    # its data within 6%, which is about +/-3% on e+, so the tolerance is
    # that plus digitisation error at the crowded right-hand edge.
    RUN_EPLUS = [80.7, 150.5, 256.2, 485.8, 512.1]
    TOL = 0.10

    seen = 0
    for series in dataset:
        if "fig4.46" not in series.label or "eD0.047_pe10_wh2" not in series.label:
            continue
        assert series.class_confidence == "confirmed", (
            f"{series.label} lost its resolution against CR-4015"
        )
        assert "CR-4015" in series.cross_check, (
            f"{series.label} no longer records how the dispute was closed"
        )
        seen += 1
        for pt in load_points(series):
            x = pt.x
            nearest = min(RUN_EPLUS, key=lambda r: abs(x / r - 1.0))
            assert abs(x / nearest - 1.0) <= TOL, (
                f"{series.label}: mark at e+={x:.1f} is "
                f"{100 * abs(x / nearest - 1):.1f}% from the nearest "
                f"CR-4015 run ({nearest}); no source run explains it"
            )
    assert seen == 2, f"expected the G and R pair, found {seen}"


def test_disputed_labels_are_refused_by_geometry_binding_sets(dataset) -> None:
    """Harmless today, refused the moment it would matter.

    han_1988_orthogonal carries zero geometry exponents, so a wrong class
    label changes nothing and the series scores normally. A set that binds
    geometry uses the label AS an input, and the runner must decline
    rather than feed a disputed one in. This guards the 4.47/4.48 work
    before it exists.
    """
    import combaero as cb
    from validation.cooling.runner import _binds_geometry

    assert not _binds_geometry(cb.han_1988_orthogonal()), (
        "han_1988_orthogonal now binds geometry; the disputed series must "
        "be resolved against the page before it can be scored again"
    )

    # A dispute must not by itself stop a series scoring. Check one that
    # is otherwise scorable: 90 degrees, so the set's valid_alpha admits
    # it, and G rather than G_bar so no ratio conversion is involved.
    from validation.cooling.runner import run_series

    scorable = [
        s
        for s in dataset
        if s.class_confidence == "disputed"
        and s.scores
        and s.alpha_deg in (None, 90.0)
        and s.y_axis == "G"
    ]
    for s in scorable:
        recs = run_series(s)
        assert any(r.predicted is not None for r in recs), (
            f"{s.label} is disputed but otherwise scorable, and was refused"
        )


def test_every_figure_redraws(tmp_path, dataset) -> None:
    """The plotter must run for every figure in the dataset.

    A visual check only helps if it still works, and it is the kind of
    tool that rots silently -- nobody notices until they need it. This
    also pins that every series carries the axis spec the plot needs, so
    a new series cannot be added without one.
    """
    pytest.importorskip("matplotlib", reason="needs the 'examples' extra")

    from validation.cooling.plot import plot_figure

    figures = {s.figure for s in dataset if s.figure}
    assert figures, "no series declares a figure"

    for figure in sorted(figures):
        out = plot_figure(figure, tmp_path / f"{figure}.png", dataset)
        assert out.exists() and out.stat().st_size > 5000, f"{figure} drew nothing"


def test_every_figure_has_a_committed_plot(dataset) -> None:
    """A figure in the dataset must have its redrawn image committed.

    The images are committed so a reviewer sees the cloud's shape in the
    diff. That only works if they exist, and the way it fails is filing a
    new figure and forgetting to regenerate -- so that case is pinned.

    This does not check the image is CURRENT; PNG bytes are not stable
    across matplotlib versions, so comparing them would fail for reasons
    that have nothing to do with the data. Regenerate with
    `uv run python -m validation.cooling.plot --all` after changing data
    or a card.
    """
    from validation.cooling.plot import plot_path

    pairs = sorted({(s.source.name, s.figure) for s in dataset if s.figure})
    for source, figure in pairs:
        img = plot_path(source, figure)
        assert img.exists(), (
            f"{source}/{figure} has no committed plot; run "
            "`uv run python -m validation.cooling.plot --all`"
        )
        assert img.stat().st_size > 5000, f"{img.name} looks empty"


def test_figure_references_are_namespaced_by_source(dataset) -> None:
    """A figure number alone is not a unique identifier.

    Two books can each have a figure 4.46. Data, cards and plots are all
    namespaced by source; a bare number is accepted only while it happens
    to be unique, and becomes an error naming the alternatives rather than
    a silent pick once it is not.
    """
    from validation.cooling.plot import resolve_figure

    source, figure = resolve_figure(dataset, "4.46")
    assert (source, figure) == ("han2012", "4.46")
    assert resolve_figure(dataset, "han2012/4.46") == ("han2012", "4.46")

    for bad in ("nope/4.46", "han2012/9.99", "9.99"):
        with pytest.raises(SystemExit):
            resolve_figure(dataset, bad)

    # Ambiguity must be refused, not resolved by luck of ordering.
    import dataclasses

    other = dataclasses.replace(
        next(s for s in dataset if s.figure == "4.46"),
        source=dataclasses.replace(
            next(s for s in dataset if s.figure == "4.46").source,
            name="someoneelse2030",
        ),
    )
    with pytest.raises(SystemExit, match="ambiguous"):
        resolve_figure([*dataset, other], "4.46")


def test_axis_specs_are_declared_for_plotting(dataset) -> None:
    """Every series needs a log/linear declaration and plot limits.

    Defaulting these would silently draw a log figure on linear axes,
    which looks plausible and is wrong -- figure 4.53 is linear while the
    rest are log-log.
    """
    for s in dataset:
        card = s.verification or {}
        assert card.get("x_axis_type") in ("log", "linear"), f"{s.label} declares no x_axis_type"
        assert card.get("y_axis_type") in ("log", "linear"), f"{s.label} declares no y_axis_type"
        assert s.figure and s.panel, f"{s.label} is not assigned to a panel"


def test_gbar_ratio_is_not_applied_off_90_degrees(dataset) -> None:
    """G_bar/G = 1.2 is a 90 degree result, not a universal one.

    Digitising both panels of figure 4.51 measures G_bar/G directly for
    seven rib configurations: 90 deg gives 1.2193, matching the printed
    1.2162 to 0.3%, while the angled ones average 1.155 -- 4.5% below it,
    against a 2.6% measurement scatter. The ratio is configuration
    dependent, so the 90 deg value is not a universal converter.

    Applying 1.2 off 90 degrees would manufacture a number that looks like
    a measurement. The runner must refuse and say why.
    """
    from validation.cooling.runner import _gbar_reason, run_series

    # The guard itself, at every angle the figures carry. It is preventive:
    # figure 4.51 brings eight angled configurations, and none of the G_bar
    # series filed today is both off 90 degrees AND scored.
    for alpha in (30.0, 45.0, 60.0):
        series = dataclasses.replace(
            next(s for s in dataset if s.y_axis == "G_bar"),
            alpha_deg=alpha,
        )
        reason = _gbar_reason(series)
        assert reason and "90 deg result" in reason, f"G_bar at {alpha} deg was not refused"
    assert (
        _gbar_reason(
            dataclasses.replace(next(s for s in dataset if s.y_axis == "G_bar"), alpha_deg=90.0)
        )
        is None
    ), "90 degrees must still convert"

    # And no off-90 G_bar series may come back with a prediction.
    for s in dataset:
        if s.y_axis != "G_bar" or s.alpha_deg in (None, 90.0):
            continue
        recs = run_series(s)
        assert all(r.predicted is None for r in recs), (
            f"{s.label} is G_bar at {s.alpha_deg} deg and was scaled by 1.2"
        )


def test_90_degree_gbar_still_scores(dataset) -> None:
    """The guard must not block where the ratio was actually established.

    Figure 4.46 prints both lines, 3.7 and 4.5, a ratio of 1.216, and
    figure 4.47 labels its dashed line G_bar = 1.2 G. At 90 degrees the
    conversion is evidenced and must keep working.
    """
    from validation.cooling.runner import run_series

    scored = 0
    for s in dataset:
        if s.y_axis != "G_bar" or s.alpha_deg != 90.0 or not s.scores:
            continue
        recs = run_series(s)
        if any(r.predicted is not None for r in recs):
            scored += 1

    assert scored, "the 90 degree G_bar series stopped scoring"


def test_candidate_mode_runs(tmp_path, dataset) -> None:
    """The digitising path must keep working, not just the filed path.

    check_candidate builds a SeriesMetadata by hand, so every field added
    to that dataclass has to be added here too. It broke silently when
    `figure`, `panel` and `class_confidence` arrived, because 262 checks
    covered the filed data and none covered the path actually used while
    picking points.
    """
    import yaml

    from validation.cooling.verify import check_candidate

    csv_path = tmp_path / "candidate.csv"
    csv_path.write_text("x, y\n100, 12.0\n300, 17.0\n900, 22.0\n")
    card_path = tmp_path / "card.yaml"
    card_path.write_text(
        yaml.safe_dump(
            {
                "verification": {
                    "x_axis_type": "log",
                    "y_axis_type": "log",
                    "x_ticks": [50, 100, 1000],
                    "x_multiplier": 1.0,
                    "y_ticks": [8, 40],
                    "printed_curve": {"power_law": {"C": 3.7, "n": 0.28}},
                    "curve_tolerance": 0.2,
                }
            }
        )
    )
    findings = check_candidate(csv_path, card_path)
    assert findings, "candidate mode produced no findings at all"
    assert all(f.ok for f in findings), [f"{f.check}: {f.detail}" for f in findings if not f.ok]


def test_candidate_mode_catches_a_dropped_multiplier(tmp_path) -> None:
    """The bug candidate mode exists to catch, end to end."""
    import yaml

    from validation.cooling.verify import check_candidate

    csv_path = tmp_path / "candidate.csv"
    # Two decades low: the real figure 4.54 and 4.51 defect.
    csv_path.write_text("x, y\n1.0, 12.0\n3.0, 17.0\n9.0, 22.0\n")
    card_path = tmp_path / "card.yaml"
    card_path.write_text(
        yaml.safe_dump(
            {
                "verification": {
                    "x_axis_type": "log",
                    "y_axis_type": "log",
                    "x_ticks": [1, 10],
                    "x_multiplier": 1.0e2,
                    "y_ticks": [8, 40],
                }
            }
        )
    )
    findings = check_candidate(csv_path, card_path)
    failed = [f.check for f in findings if not f.ok]
    assert "x-span" in failed, f"dropped multiplier not caught; got {failed}"


def test_figure_448_resolves_extraction_item_41(dataset) -> None:
    """The W/H=1/2 branch boundary: which side does Han's own figure use?

    Extraction item 41 found Eq. 4.19's two branches disagreeing by ~10%
    at W/H=0.5, unresolved. Figure 4.48's own drawn correlation at that
    exact W/H settles it: its exponent must sit far closer to the wide
    branch (n=0.35) than the narrow one (n=0.258).
    """
    series = next(s for s in dataset if s.path.name == "fig4.48_G_correlation_WH0.5.csv")
    pts = load_points(series)
    fitted = _power_slope_test(pts)
    assert abs(fitted - 0.35) < abs(fitted - 0.258), (
        f"fitted exponent {fitted:.4f} is closer to the narrow branch "
        "(0.258) than the wide one (0.35) -- item 41's resolution direction "
        "has reversed"
    )
    assert 0.30 < fitted < 0.36


def _power_slope_test(points) -> float:
    xs = [math.log(p.x) for p in points]
    ys = [math.log(p.y) for p in points]
    n = len(points)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(v * v for v in xs)
    sxy = sum(a * b for a, b in zip(xs, ys, strict=False))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def test_figure_448_R_axis_is_linear_not_log(dataset) -> None:
    """Panel (a)'s y-axis was initially misread as log; it is linear.

    The tell is the digitised frame's own picks: 1.5, 2.5, 3.5, 4.5 are
    half-integer minor ticks, which a log axis cannot produce -- its
    conventional minor ticks sit at the same mantissa (2..9) in every
    decade, never at a half-integer. Confirmed on the printed page too:
    one minor tick sits exactly centred between each major label.
    """
    for s in dataset:
        if s.figure == "4.48" and s.panel == "top":
            card = s.verification or {}
            assert card.get("y_axis_type") == "linear", (
                f"{s.label} declares y_axis_type={card.get('y_axis_type')}, "
                "but panel (a)'s R axis is linear"
            )


def test_figure_448_nothing_is_scored_against_the_wrong_paper(dataset) -> None:
    """Overlapping data must not silently score against a cross-paper gap.

    The '[Ref. 3]' W/H=1 curve on figure 4.48 sits inside
    han_1988_orthogonal's valid range, but it reproduces Eq. 4.18
    (C=2.24, n=0.35), a DIFFERENT published form from the one
    han_1988_orthogonal encodes (Eq. 4.15/16, n=0.28). Scoring it would
    report that ~15% cross-paper gap as model error. Every figure 4.48
    series is therefore unscored, regardless of geometry.
    """
    fig448 = [s for s in dataset if s.figure == "4.48"]
    assert fig448, "figure 4.48 is not filed"
    assert all(s.scores is None for s in fig448), [s.label for s in fig448 if s.scores is not None]


def test_han_park_1988_angled_scores_its_own_source(records) -> None:
    """han_park_1988_angled must reproduce figure 4.47, its own source.

    Both R and G panels scored through evaluate_rib with a representative
    geometry (the cloud carries no legend). Bands are what this project
    measured, not an author claim -- see the C++ set's accuracy_R/accuracy_G
    comment for why that distinction matters here specifically.
    """
    from validation.cooling.scorecard import build

    cells = {c.label: c for c in build(records)}
    r_cell = cells["han2012/fig4.47_R_vs_alpha"]
    g_cell = cells["han2012/fig4.47_G_vs_eplus"]

    assert r_cell.n_scored == 39
    assert r_cell.rmse < 0.12
    assert g_cell.n_scored == 115
    assert g_cell.rmse < 0.10


def test_han_park_1988_angled_confirmed_independently_by_figure_451(
    records,
) -> None:
    """A different paper's parallel-rib data, scored against this equation.

    Figure 4.51 (Han et al. 1991) is not han_park_1988_angled's own source
    (Han and Park 1988) -- an independent cross-paper check, on the two
    classes (45/60 deg PARALLEL) that match the rib shape Eq. 4.17/4.18
    was fitted to. The other seven shapes on figure 4.51 (crossed, V,
    lambda) are deliberately not scored against this set: it has no way to
    represent rib shape beyond angle, and scoring a shape it cannot model
    would repeat the mistake extraction item 23 already recorded once.
    """
    from validation.cooling.scorecard import build

    cells = {c.label: c for c in build(records)}
    scored_parallel = [
        c for name, c in cells.items() if "par" in name and "4.51" in name and not c.unsupported
    ]
    assert scored_parallel, "no figure 4.51 parallel class scored"
    for c in scored_parallel:
        assert c.rmse < 0.15, f"{c.label}: RMSE {c.rmse:.1%}"

    # The other seven shapes must NOT be silently scored against this set.
    unscored_shapes = [
        c
        for name, c in cells.items()
        if "4.51" in name and any(t in name for t in ("crs", "vee", "lam")) and c.n > 0
    ]
    assert unscored_shapes
    assert all(c.unsupported for c in unscored_shapes), [
        c.label for c in unscored_shapes if not c.unsupported
    ]


def test_han_park_r_alpha_normalisation_is_not_raw_R(records) -> None:
    """Figure 4.47's R panel plots R divided by its own normalisers.

    Falsifiable directly: if the runner scored raw R instead of the
    normalised quantity, the bias would be enormous (R_pe and the W/H^m
    term are both > 1 across the validity box), not a few percent.
    """
    from validation.cooling.scorecard import build

    cells = {c.label: c for c in build(records)}
    r_cell = cells["han2012/fig4.47_R_vs_alpha"]
    assert abs(r_cell.bias) < 0.15, (
        f"bias {r_cell.bias:.1%} -- looks like raw R is being compared "
        "against the normalised digitised values"
    )


def test_fig446_span_check_is_one_sided_and_its_exceedances_are_known(dataset) -> None:
    """The fig4.46 cross_checks record that their span argument is ONE-SIDED:
    only isolated marks can be picked, so the picked span is a subset of the
    true span and fitting inside the predicted band is close to vacuous.
    Exceeding it is the informative direction.

    This pins both halves of that claim so neither can rot:

      1. the check really is weak where it passes -- mean coverage of the
         predicted range is well under 100%;
      2. the series that DO exceed the band are exactly the four recorded.

    A re-pick that moves a series across the boundary changes this list and
    should make someone look, rather than silently joining a set of series
    whose prose says they agree.
    """
    import csv as _csv
    import math

    # Figure 4.46's own legend: "This study" 10,000-60,000 (e/D 0.047, 0.078);
    # "Han (1984)" 8,000-80,000 (e/D 0.063, 0.042, 0.021).
    RE_RANGE = {0.047: (10000, 60000), 0.078: (10000, 60000)}
    DEFAULT_RE = (8000, 80000)
    MARGIN = 0.15  # f comes from the printed R = 3.2, worth about +/-3% on e+
    KNOWN_EXCEEDANCES = {
        "fig4.46_G_eD0.063_pe10_wh1",
        "fig4.46_R_eD0.042_pe10_wh1",
        "fig4.46_R_eD0.063_pe10_wh1",
        "fig4.46_R_eD0.063_pe20_wh1",
    }

    def e_plus(re_d: float, e_over_d: float, w_over_h: float) -> float:
        geom = 2 * e_over_d * (2 * w_over_h / (w_over_h + 1.0))
        f = 2.0 / (3.2 - 2.5 * math.log(geom) - 2.5) ** 2
        return e_over_d * re_d * math.sqrt(f / 2.0)

    coverage: list[float] = []
    exceed: set[str] = set()
    for series in dataset:
        stem = series.path.stem
        if not stem.startswith("fig4.46_") or series.kind != "measured":
            continue
        geom = series.geometry or {}
        e_over_d, w_over_h = geom.get("e_D"), geom.get("W_H")
        if e_over_d is None or w_over_h is None:
            continue
        with open(series.path) as fh:
            xs = [float(r[0]) for r in list(_csv.reader(fh))[1:] if r]
        if not xs:
            continue
        lo_re, hi_re = RE_RANGE.get(e_over_d, DEFAULT_RE)
        lo, hi = e_plus(lo_re, e_over_d, w_over_h), e_plus(hi_re, e_over_d, w_over_h)
        coverage.append((max(xs) / min(xs)) / (hi / lo))
        if any(x < lo * (1 - MARGIN) or x > hi * (1 + MARGIN) for x in xs):
            exceed.add(stem)

    assert coverage, "no fig4.46 measured series found"
    mean_cov = sum(coverage) / len(coverage)
    assert mean_cov < 0.95, (
        f"picked marks now cover {mean_cov:.0%} of the predicted range; the "
        "cross_checks describe this check as weak because coverage is partial"
    )
    assert exceed == KNOWN_EXCEEDANCES, (
        f"the set of series exceeding the predicted band changed: "
        f"added {sorted(exceed - KNOWN_EXCEEDANCES)}, "
        f"gone {sorted(KNOWN_EXCEEDANCES - exceed)}"
    )


def test_rib_angle_reaches_the_correlation_on_the_eplus_path() -> None:
    """A series' rib angle must reach `evaluate_rib`, not be replaced by the
    probe geometry's default of 90 degrees.

    It was. `run_series` copied `e_D`, `p_e` and `W_H` out of
    `series.geometry` but never `alpha_deg`, which is a top-level field
    rather than part of that mapping. Every angled series on the e+ path was
    therefore scored as though its ribs were transverse -- against
    `han_park_1988_angled`, the one set whose entire subject is rib angle.

    Checked against the library directly rather than against a stored
    number, so the test cannot drift with the data.
    """
    import combaero as cb
    from validation.cooling.runner import run_series

    dataset = load_dataset()
    series = next(s for s in dataset if s.path.stem == "fig4.51_G_60par")
    assert series.alpha_deg == 60.0

    record = next(r for r in run_series(series) if r.predicted is not None)
    geometry = cb.RibGeometry(e_D=0.0625, p_e=10.0, W_H=1.0, alpha_deg=series.alpha_deg)
    at_declared = cb.evaluate_rib(cb.han_park_1988_angled(), geometry, record.re_used)
    geometry.alpha_deg = 90.0
    at_ninety = cb.evaluate_rib(cb.han_park_1988_angled(), geometry, record.re_used)
    # The two must be far enough apart for this to be a real check.
    assert abs(at_declared.G / at_ninety.G - 1.0) > 0.05

    assert abs(record.predicted / at_declared.G - 1.0) < 0.05, (
        f"prediction {record.predicted:.3f} is not the declared "
        f"alpha={series.alpha_deg} value ({at_declared.G:.3f}); it matches "
        f"alpha=90 ({at_ninety.G:.3f}), so the rib angle is being dropped"
    )
