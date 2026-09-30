"""The multi-row film/effusion scoring runner.

Ground truth here is the runner's OWN contract -- which rows are upstream,
how the two sides are averaged, where the row comb sits -- not the physics,
which `tests/test_film_effectiveness.cpp` and `test_film_superposition.cpp`
already pin against the papers. What this file guards is the bookkeeping
between the data and those correlations, because that is where a scoring
harness quietly lies.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import combaero as cb  # noqa: E402
from validation.cooling import film_runner as fr  # noqa: E402
from validation.cooling.schema import load_dataset, load_points  # noqa: E402


@pytest.fixture(scope="module")
def andrei():
    series = [s for s in load_dataset() if s.label.startswith("andrei2014/")]
    assert series, "andrei2014 is not in the dataset"
    return series


@pytest.fixture(scope="module")
def one(andrei):
    return next(s for s in andrei if s.label.endswith("BR1_DR1"))


def _contributing_rows(series, x: float) -> int:
    """How many rows the runner counts as upstream of ``x``."""
    return sum(1 for x_row in fr.row_positions(series) if x_row < x)


def test_owns_only_the_effusion_effectiveness_series(andrei) -> None:
    """Dispatch is an explicit predicate, so it must be exactly right.

    The first version of `run_dataset` dispatched on truthiness and the
    jet-array runner swallowed the orifice series; every runner returns
    reason-carrying records, so "it returned something" means nothing.
    """
    assert all(fr.owns(s) for s in andrei)

    # murray2018 is the other family this runner owns: same quantity, but
    # its abscissa is x/D where andrei2014's is x/s_x.
    murray = [s for s in load_dataset() if s.label.startswith("murray2018/")]
    assert murray, "murray2018 is not in the dataset"
    curves = [s for s in murray if s.kind != "frame"]
    assert all(fr.owns(s) for s in curves), "a murray2018 curve went unclaimed"

    others = [s for s in load_dataset() if not s.label.startswith(("andrei2014/", "murray2018/"))]
    claimed = [s.label for s in others if fr.owns(s)]
    assert not claimed, f"film runner claims series it does not score: {claimed}"


def test_row_comb_comes_from_metadata_and_closes_on_the_data(one) -> None:
    """18 rows, unit pitch, and the comb must fit the committed abscissae.

    The phase is an extraction fact recovered from the data's own injection
    points (metadata header, cross-checked three ways). This asserts the
    consequence that makes it checkable: the comb has to leave room before
    row 1 and after row 18 and no room for a 19th.
    """
    rows = fr.row_positions(one)
    assert len(rows) == 18
    gaps = [rows[i + 1] - rows[i] for i in range(len(rows) - 1)]
    assert all(g == pytest.approx(1.0) for g in gaps)

    xs = [p.x for p in load_points(one)]
    assert min(xs) < rows[0], "no samples ahead of row 1"
    assert rows[-1] < max(xs), "row 18 sits past the end of the data"
    assert max(xs) - rows[-1] < 1.0, "there is room for a 19th row"
    assert rows[0] - min(xs) < 0.1, "an implausible run-up before row 1"


def test_effectiveness_is_negligible_before_the_first_row(one) -> None:
    """Physics the comb has to respect: no coolant upstream of row 1.

    Measured, not exactly zero: the two committed samples ahead of row 1 read
    0.0 and 0.00065, the latter 0.18% of this curve's final effectiveness.
    That is the measurement, not a comb error -- so the bound is "negligible
    against the plate's own scale", not "identically zero". A phase placed a
    pitch too late would put a fully developed eta of ~0.06 here and blow
    straight through it.
    """
    ahead = [p.y for p in load_points(one) if p.x < fr.row_positions(one)[0]]
    assert ahead, "the phase leaves no samples ahead of row 1 to check"
    final = max(p.y for p in load_points(one))
    assert max(ahead) < 0.01 * final

    # The runner itself must offer nothing at all up there.
    value, _ = fr.predict(one, fr.row_positions(one)[0] - 0.01)
    assert value is None, "predicted cooling upstream of every hole"


def test_only_upstream_rows_contribute(one) -> None:
    """A row downstream of the point must not cool it.

    Off-by-one here would be invisible in the aggregate -- it shifts every
    prediction by one pitch, which looks like a modest bias rather than a
    bug.
    """
    rows = fr.row_positions(one)
    with cb.suppress_warnings():
        just_after_first, _ = fr.predict(one, rows[0] + 0.01)
        just_after_second, _ = fr.predict(one, rows[1] + 0.01)
    assert just_after_first is not None and just_after_second is not None
    assert just_after_second > just_after_first, "the second row added nothing"

    # The row COUNT is what the boundary controls, so assert that rather than
    # the value. Straddling a row by 1e-6 does NOT raise eta: the row just
    # crossed sits at x/D ~ 1e-5 and contributes essentially nothing, while
    # every established row has decayed over the same step, so the total dips
    # very slightly. That is correct, and it is why a value comparison at the
    # boundary would be the wrong assertion.
    for k in (1, 5, 12):
        assert _contributing_rows(one, rows[k] - 1e-6) == k
        assert _contributing_rows(one, rows[k] + 1e-6) == k + 1


def test_both_sides_are_averaged_over_the_same_abscissae(one) -> None:
    """The comparison must be mean-to-mean, not mean-to-point-value.

    An earlier version averaged the measurement across each pitch but
    evaluated the prediction once at the interval centre. eta climbs steeply
    after each injection, so over row 1 -- where it rises from zero -- the
    centre value overstates the interval mean substantially. Averaging both
    sides over the SAME sample abscissae also means an uneven digitisation
    cannot bias one side alone.
    """
    with cb.suppress_warnings():
        records = fr.run_series(one)
    first = records[0]

    intervals = fr.row_intervals(one, load_points(one))
    x_mid, _end, xs, _ys = intervals[0]

    with cb.suppress_warnings():
        at_centre, _ = fr.predict(one, x_mid)
        over_interval = [v for v, _ in (fr.predict(one, x) for x in xs) if v is not None]
    mean_over_interval = sum(over_interval) / len(over_interval)

    assert first.predicted == pytest.approx(mean_over_interval)
    # And the two genuinely differ over the first row, which is why it matters.
    assert abs(at_centre - mean_over_interval) > 0.01 * mean_over_interval


def test_every_record_is_flagged_extrapolated(andrei) -> None:
    """andrei2014 is outside Baldauf's envelope on every curve.

    s/D 7.37 against a published maximum of 5 does it alone, on all six.
    This is what makes the whole dataset a cross-source ACCURACY measurement
    rather than a fidelity check, so it must never be reported unflagged.
    The flag comes from the correlation, not from a judgement in the runner.
    """
    for series in andrei:
        with cb.suppress_warnings():
            records = [r for r in fr.run_series(series) if r.predicted is not None]
        assert records, series.label
        assert all(r.extrapolated for r in records), series.label


def test_the_runner_does_not_reimplement_the_correlation(one) -> None:
    """Predictions must move when the real correlation moves.

    A runner that inlined the algebra would keep answering after the library
    changed underneath it, and the scorecard would report a stale model.
    Superposition is monotone in each row's eta, so raising the turbulence
    intensity -- a real Baldauf input -- must move the prediction.
    """
    rows = fr.row_positions(one)
    x = rows[6] + 0.5
    with cb.suppress_warnings():
        low, _ = fr.predict(one, x, Tu=0.01)
        high, _ = fr.predict(one, x, Tu=0.07)
    assert low != high, "prediction ignores an input the correlation uses"


def test_turbulence_sensitivity_is_reported_not_minimised(one) -> None:
    """Andrei state no Tu, so the assumption's cost must be visible.

    Reporting it is the point. Choosing the Tu that scores best would be
    fitting an unmeasured input against the metric that judges it.
    """
    with cb.suppress_warnings():
        sens = fr.turbulence_sensitivity(one)
    assert set(sens) == {0.01, 0.05, 0.075}
    assert all(v > 0.0 for v in sens.values())
    assert fr.ASSUMED_TU in sens, "the assumed value is not among those reported"


def test_gao_alpha_cannot_be_fitted_on_this_family(andrei) -> None:
    """The knob is one-sided, and this family needs both directions.

    Gao's alpha is bounded [0, 1] and only ever REDUCES the superposed value
    (pinned in test_film_superposition.cpp). Andrei's six curves split three
    over-predicted and three under-predicted, so no single (a, b) can serve
    both halves -- fitting it here would trade one against the other rather
    than measure anything.

    Recorded as a test because the roadmap said to fit (a, b) on this family.
    If a future change makes every curve fall the same side, this test goes
    red and the plan becomes viable again -- which is exactly when someone
    should be told.

    NOT a verdict on alpha itself, and the distinction matters. Murray &
    Ireland (2018) applied Sellers to a single-hole CFD result -- no
    correlation involved -- and measured a ~2x OVER-prediction at M ~ 1 for a
    5.75D pitch, growing with blowing, caused by streamwise jet interaction
    and largely cured by tripling the streamwise pitch. So the error alpha
    exists to correct is real, one-signed, and exactly alpha's shape. What
    disqualifies THIS family is that Baldauf's lift-off collapse above
    M ~ 1 drives BR 2-3 the other way and masks it. Murray's blowing range
    of 0.1-1.2 sits entirely below that peak, which is where alpha can be
    measured cleanly. See
    validation/cooling/extractions/murray_ireland_2018_effusion_superposition.md.
    """
    biases = []
    for series in andrei:
        with cb.suppress_warnings():
            errs = [r.rel_error for r in fr.run_series(series) if r.rel_error is not None]
        assert errs, series.label
        biases.append(sum(errs) / len(errs))

    over = [b for b in biases if b > 0]
    under = [b for b in biases if b < 0]
    assert over and under, (
        "all six curves now fall the same side of the data; Gao's alpha could "
        f"be fitted after all -- biases {[round(b, 3) for b in biases]}"
    )


@pytest.fixture(scope="module")
def murray():
    series = [s for s in load_dataset() if s.label.startswith("murray2018/")]
    assert series, "murray2018 is not in the dataset"
    return series


def _interp(xs, ys, x):
    if x < xs[0] or x > xs[-1]:
        return None
    for i in range(1, len(xs)):
        if xs[i] >= x:
            span = xs[i] - xs[i - 1]
            t = (x - xs[i - 1]) / span if span else 0.0
            return ys[i - 1] + t * (ys[i] - ys[i - 1])
    return None


def _points(series):
    pts = sorted((p.x, p.y) for p in load_points(series))
    return [a for a, _ in pts], [b for _, b in pts]


def test_murray_row_comb_is_staggered_half_pitch(murray) -> None:
    """10 rows at half the primary pitch, not 5 at the full pitch.

    The plate is staggered, so the film meets a row every 5.75/2 = 2.875 D.
    Measured rather than assumed: clustering the steepest rises on all five
    curves gives a dominant gap of 2.7-3.1. It also closes the hole count --
    5 primary rows plus 5 staggered is 10, and the paper's Geometry 2 is a
    4 x 5 primary array with 40 holes including the staggered ones.
    """
    curve = next(s for s in murray if s.label.endswith("exp_M0p19"))
    rows = fr.row_positions(curve)
    assert len(rows) == 10
    gaps = [rows[i + 1] - rows[i] for i in range(len(rows) - 1)]
    assert all(g == pytest.approx(5.75 / 2.0) for g in gaps)

    # The comb must end before the data does -- the plate stops at row 10
    # and the film then decays with no further injection.
    xs, ys = _points(curve)
    assert rows[-1] < max(xs), "the last row sits past the data"
    tail = [y for x, y in zip(xs, ys, strict=False) if x > rows[-1] + 0.5]
    assert len(tail) >= 3, "no post-plate tail to check"
    assert tail[-1] < tail[0], "eta does not decay after the last row"


def test_murray_abscissa_is_diameters_not_pitches(murray) -> None:
    """The runner must read x/D and x/s_x differently.

    `andrei2014` plots streamwise distance in PITCHES and `murray2018` in
    DIAMETERS. Treating one as the other misplaces every row by the pitch
    ratio -- here a factor of 2.875 -- which would look like a modest bias
    rather than a bug.
    """
    m = next(s for s in murray if s.label.endswith("exp_M0p19"))
    a = next(s for s in load_dataset() if s.label.startswith("andrei2014/"))
    assert m.x_axis == "x_over_D" and a.x_axis == "x_over_sx"
    assert fr._diameters_per_x(m) == 1.0
    assert fr._diameters_per_x(a) == pytest.approx(float(a.geometry["sx_over_d"]))


def test_murray_is_uniformly_over_predicted(murray) -> None:
    """Every Murray curve falls the SAME side, unlike andrei2014.

    This is what makes the family usable for fitting Gao's alpha, which is
    bounded [0, 1] and can only ever reduce. `andrei2014` split three over
    and three under -- see the test above -- because Baldauf's lift-off
    collapse drives its BR 2-3 curves the other way. Murray runs M 0.19 to
    0.96, entirely below that peak, so the closure does not change sign and
    the superposition error is measurable on its own.

    If this ever goes red, the family has stopped being a valid fitting set
    and any alpha fitted on it must be revisited.
    """
    biases = {}
    for series in (s for s in murray if s.kind == "measured"):
        with cb.suppress_warnings():
            errs = [r.rel_error for r in fr.run_series(series) if r.rel_error is not None]
        assert errs, series.label
        biases[series.label.split("/")[1]] = sum(errs) / len(errs)

    assert len(biases) == 3
    assert all(b > 0 for b in biases.values()), (
        f"not all Murray curves over-predict, so alpha cannot serve them all: {biases}"
    )
    # And the over-prediction must grow with blowing, which is the trend
    # alpha has to reproduce.
    assert biases["fig6_eta_exp_M0p96"] > biases["fig6_eta_exp_M0p19"]


def test_superposition_is_the_larger_error_not_the_closure(murray) -> None:
    """The decomposition #420 could not make.

    Figure 6 carries the paper's OWN Sellers superposition over a
    single-hole CFD field, so three curves share one axis:

        their/exp   isolates the SUPERPOSITION -- their per-row input is CFD
        ours/their  isolates the CLOSURE       -- both use Sellers
        ours/exp    the total

    Measured over x/D 6-27: superposition 1.18x at M = 0.19 rising to 1.74x
    at M = 0.96, while the closure only moves 1.13x to 1.28x. The
    superposition is the dominant term and the one that responds to blowing,
    which is what justifies spending a correction on it.
    """
    ratios = {}
    for tag in ("M0p19", "M0p96"):
        exp = next(s for s in murray if s.label.endswith(f"exp_{tag}"))
        sup = next(s for s in murray if s.label.endswith(f"sup_{tag}"))
        ex, ey = _points(exp)
        sx, sy = _points(sup)
        sup_r, clo_r = [], []
        with cb.suppress_warnings():
            for st in range(6, 28):
                e = _interp(ex, ey, st)
                t = _interp(sx, sy, st)
                o, _ = fr.predict(exp, float(st))
                if not (e and t and o) or e < 0.02:
                    continue
                sup_r.append(t / e)
                clo_r.append(o / t)
        assert sup_r, tag
        ratios[tag] = (sum(sup_r) / len(sup_r), sum(clo_r) / len(clo_r))

    for tag, (sup_err, clo_err) in ratios.items():
        assert sup_err > 1.0, f"{tag}: superposition does not over-predict"
        assert sup_err > clo_err, (
            f"{tag}: the closure ({clo_err:.2f}x) now exceeds the "
            f"superposition ({sup_err:.2f}x); the decomposition has changed"
        )

    # The superposition error responds to blowing far more than the closure.
    d_sup = ratios["M0p96"][0] - ratios["M0p19"][0]
    d_clo = ratios["M0p96"][1] - ratios["M0p19"][1]
    assert d_sup > 2 * d_clo, f"superposition grew {d_sup:.2f} against the closure's {d_clo:.2f}"


def test_murray_calibration_ticks_land_on_round_numbers() -> None:
    """The digitisation's own calibration evidence, checked.

    `fig6_calibration_M0p*.csv` were digitised in DATA coordinates, so the
    tick marks must fall on the printed round numbers. They are not declared
    as series: `kind: frame` means a straight panel EDGE in this harness,
    fitted as a power law to measure skew, and these are a traversal of
    ticks and corners that includes x = 0 -- which makes `log(x)` blow up.
    Checking them here is both correct and stronger than a slope fit.

    Traversal order, stated by the digitiser and load-bearing here: from the
    origin counter-clockwise -- seven x ticks along the bottom, then the
    bottom-right, top-right and top-left corners, then six y ticks down the
    left side. These files are therefore NOT sorted.

    Three panels digitised independently agreeing to ~0.4% is the real
    cross-check; a log/linear slip or a misread tick would not reproduce
    across all three.
    """
    import csv

    x_ticks = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    y_ticks = [0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    x_span, y_span = 32.0, 0.7

    folder = Path(__file__).resolve().parents[2] / "validation" / "cooling" / "data" / "murray2018"
    files = sorted(folder.glob("fig6_calibration_*.csv"))
    assert len(files) == 3, f"expected three calibration panels, got {len(files)}"

    for path in files:
        rows = [
            (float(r[0]), float(r[1]))
            for r in csv.reader(path.open())
            if not r[0].strip().startswith("x")
        ]
        assert len(rows) == 16, f"{path.name}: {len(rows)} points, expected 16"

        xs, corners, ys = rows[0:7], rows[7:10], rows[10:16]

        worst_x = max(abs(p[0] - t) for p, t in zip(xs, x_ticks, strict=True))
        worst_y = max(abs(p[1] - t) for p, t in zip(ys, y_ticks, strict=True))
        assert worst_x / x_span < 0.005, (
            f"{path.name}: x ticks off by {worst_x:.3f} ({worst_x / x_span:.2%} of span)"
        )
        assert worst_y / y_span < 0.01, (
            f"{path.name}: y ticks off by {worst_y:.4f} ({worst_y / y_span:.2%} of span)"
        )

        # x ticks lie on y = 0 and y ticks on x = 0; that is what makes them
        # ticks rather than arbitrary points, and it catches a traversal
        # read in the wrong order.
        assert max(abs(p[1]) for p in xs) < 0.01, f"{path.name}: x ticks off the axis"
        assert max(abs(p[0]) for p in ys) < 0.2, f"{path.name}: y ticks off the axis"

        # The three corners close the plot area at (32, 0), (32, 0.7), (0, 0.7).
        br, tr, tl = corners
        assert abs(br[0] - x_span) < 0.2 and abs(br[1]) < 0.01, f"{path.name}: BR"
        assert abs(tr[0] - x_span) < 0.2 and abs(tr[1] - y_span) < 0.01, f"{path.name}: TR"
        assert abs(tl[0]) < 0.2 and abs(tl[1] - y_span) < 0.01, f"{path.name}: TL"


def _predict_with_alpha(series, x, alpha):
    """Superposed eta at x with a uniform per-row alpha."""
    M, P = fr._blowing_and_density(series)
    g = series.geometry
    rows = []
    for x_row in fr.row_positions(series):
        if x_row >= x:
            break
        rows.append(
            cb.film_effectiveness_baldauf_2002(
                (x - x_row) * fr._diameters_per_x(series),
                M,
                P,
                float(g["angle_deg"]),
                float(g["sz_over_d"]),
                fr.ASSUMED_TU,
            )
        )
    if not rows:
        return None
    if len(rows) == 1:
        return rows[0]
    return cb.film_superposition_corrected(rows, [alpha] * (len(rows) - 1))


def _mae_at_alpha(series, alpha):
    pts = sorted((p.x, p.y) for p in load_points(series))
    errs = []
    for x, y in pts:
        v = _predict_with_alpha(series, x, alpha)
        if v is not None and y > 0.02:
            errs.append(abs(v / y - 1.0))
    return sum(errs) / len(errs)


def test_a_constant_alpha_helps_and_helps_most_where_the_error_is_worst(murray) -> None:
    """One knob is worth having, which is why the fit was attempted at all.

    Sellers uncorrected gives 33%, 44% and 113% MAE at M = 0.19, 0.48 and
    0.96. A single constant alpha takes those to roughly 16%, 21% and 24% --
    at M = 0.96 that is 113% down to 24%, and the residual is close to the
    paper's own stated 15% experimental uncertainty.
    """
    with cb.suppress_warnings():
        for tag, alpha, ceiling in (
            ("M0p19", 0.85, 0.20),
            ("M0p48", 0.85, 0.25),
            ("M0p96", 0.69, 0.28),
        ):
            s = next(x for x in murray if x.label.endswith(f"exp_{tag}"))
            corrected = _mae_at_alpha(s, alpha)
            uncorrected = _mae_at_alpha(s, 1.0)
            assert corrected < ceiling, f"{tag}: {corrected:.1%} at alpha {alpha}"
            assert corrected < uncorrected, f"{tag}: alpha made it worse"

        # It earns its keep at high blowing in particular.
        s96 = next(x for x in murray if x.label.endswith("exp_M0p96"))
        assert _mae_at_alpha(s96, 1.0) > 1.0, "M = 0.96 should be over 100% uncorrected"
        assert _mae_at_alpha(s96, 0.69) < 0.30


def test_gao_alpha_form_cannot_be_fitted_on_murray(murray) -> None:
    """The negative result, pinned so it is not quietly re-attempted.

    Gao Eq. (5) is `alpha = a r/(a r + 1) + b`, strictly MONOTONE in r --
    `d/dr = a/(a r + 1)^2` never changes sign. The best constant alpha per
    blowing ratio runs 0.85, 0.85, 0.69: flat, then falling. No monotone
    function passes through that, and for `a > 0` the form is *increasing*,
    meaning more coolant needs less correction while Murray needs more.

    The assertion is on the ORDERING, not on the fitted numbers, because
    that is what rules the form out. It is also invariant to how `r` is
    defined: any `m_coolant/m_mainstream` differs only by a positive scale
    for one geometry, and `r -> k r` is the same curve with `a -> a/k`.

    Full attempt, including the hold-one-out and the cumulative-r variant,
    in validation/cooling/extractions/gao_alpha_fit_on_murray.md.
    """
    best = {}
    with cb.suppress_warnings():
        for tag in ("M0p19", "M0p48", "M0p96"):
            s = next(x for x in murray if x.label.endswith(f"exp_{tag}"))
            best[tag] = min((0.30 + 0.01 * i for i in range(71)), key=lambda a: _mae_at_alpha(s, a))

    lo, mid, hi = best["M0p19"], best["M0p48"], best["M0p96"]
    # The high-blowing case needs markedly MORE correction than either
    # lower one -- that is the direction Gao's form gets wrong.
    assert hi < lo - 0.10, f"M=0.96 alpha {hi:.2f} not well below M=0.19 {lo:.2f}"
    assert hi < mid - 0.10, f"M=0.96 alpha {hi:.2f} not well below M=0.48 {mid:.2f}"
    # And the two low-blowing cases sit together, so the sequence is not
    # monotone decreasing either -- it is flat, then falls.
    assert abs(mid - lo) < 0.05, (
        f"M=0.19 and M=0.48 alphas ({lo:.2f}, {mid:.2f}) are no longer "
        "together; the sequence may have become monotone and Gao's form "
        "should be re-tested"
    )
