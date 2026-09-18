"""Scored validation of the cooling correlations against digitised data.

These are not self-consistency checks. Every expectation here is a number
read off a published figure, and the bands come from the source's own
stated accuracy -- not from what the implementation currently returns.
See validation/cooling/ and issue #333.
"""

from __future__ import annotations

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
    """
    recs = run_series(_series(dataset, "fig4.193c_old_correlation_G_curve"))
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
