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

    others = [s for s in load_dataset() if not s.label.startswith("andrei2014/")]
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
