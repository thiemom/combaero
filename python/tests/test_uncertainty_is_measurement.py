"""`uncertainty` must be a MEASUREMENT band, never the model's own error.

The defect this guards against (#389) lived in the data, not in a test.
Every `florschuetz1981` series carried a per-series `uncertainty` that had
been set from the model's measured error on that same series -- the series
whose `cross_check` reads "bias +31%, RMS 32%" carried a band of 0.324. So
`within` asked whether the model sits inside its own error, which it largely
must, and reported a meaningless 58.6%.

Measured before the fix: the declared band correlated with the measured model
RMS at r = +0.9987 over 27 series.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.cooling.schema import load_dataset  # noqa: E402
from validation.cooling.scorecard import build, run_dataset  # noqa: E402


@pytest.fixture(scope="module")
def scored():
    """Series that are actually scored, with their cell."""
    dataset = load_dataset()
    cells = {c.label: c for c in build(run_dataset(dataset), dataset)}
    return [
        (s, cells[s.label]) for s in dataset if s.label in cells and cells[s.label].n_scored > 0
    ]


def _pearson(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys, strict=False))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def test_band_is_constant_within_a_figure_and_quantity(scored) -> None:
    """The structural guard, and the one that actually binds.

    A measurement band has exactly two components and each has a scope:

      experimental uncertainty -> a property of the INSTRUMENT, so constant
                                  across a source
      digitisation precision   -> a property of the FIGURE, so constant
                                  across a figure

    Nothing about a measurement varies series-to-series within one figure and
    one quantity. A band that does is a property of the MODEL, which is the
    defect.

    Restricted to SCORED series on purpose: an unscored band judges nothing,
    so it cannot have been derived from the error it judges. `han2012`'s
    error-bar series are the case in point -- their band is the measured
    width of the printed bars (8.8%, 9.6%, 8.6%), a real measurement
    statement, and they carry `scores: null`.

    This is strictly stronger than the correlation check below: it caught two
    `han2012` groups that r = +0.235 did not.
    """
    groups: dict[tuple[str, str, str], list[tuple[float, str]]] = {}
    for series, _cell in scored:
        if series.uncertainty is None:
            continue
        key = (series.label.split("/")[0], str(series.figure), series.y_axis)
        groups.setdefault(key, []).append((series.uncertainty, series.label.split("/")[1]))

    offenders = {
        key: sorted({u for u, _ in v}) for key, v in groups.items() if len({u for u, _ in v}) > 1
    }
    assert not offenders, (
        "uncertainty varies between scored series of the same figure and "
        f"quantity, so it is not a measurement property: {offenders}"
    )
    assert groups, "no scored series carried a band -- the guard would be vacuous"


def test_band_does_not_track_the_model_error_it_judges(scored) -> None:
    """The check #389 asked for, kept as a second, independent signal.

    Weaker than the structural test above and it needs a caveat, so it is not
    the primary guard: with few series and few distinct bands a legitimate
    per-quantity assignment can correlate by coincidence. `lau1990` does
    exactly that -- three series, two bands (+/-5.8% Stanton for G and G_bar,
    +/-10.9% friction for R, both from Kline and McClintock in the paper),
    giving r = +0.97 purely because friction happens to have both the larger
    band and the larger error. Requiring several distinct bands keeps that
    correct assignment from being flagged.
    """
    by_source: dict[str, list[tuple[float, float]]] = {}
    for series, cell in scored:
        if series.uncertainty is None or not math.isfinite(cell.rmse):
            continue
        by_source.setdefault(series.label.split("/")[0], []).append((series.uncertainty, cell.rmse))

    suspicious = {}
    for source, pairs in by_source.items():
        bands = [u for u, _ in pairs]
        # Fewer than four distinct bands cannot distinguish a derived band
        # from a per-quantity one; see the lau1990 note above.
        if len(set(bands)) < 4:
            continue
        r = _pearson(bands, [e for _, e in pairs])
        if r is not None and abs(r) > 0.8:
            suspicious[source] = round(r, 4)

    assert not suspicious, (
        "declared uncertainty tracks the measured model error, so `within` "
        f"is asking whether the model is inside its own error: {suspicious}"
    )


def test_florschuetz_carries_the_measurement_band_it_should(scored) -> None:
    """The specific repair, pinned so it cannot silently regress.

    experimental +/-5% at 95% confidence (the primary paper's own Nu
    measurement uncertainty, extraction item 12) combined in root-sum-square
    with 0.97% digitisation (worst y-tick deviation over all nine panels):

        sqrt(0.05^2 + 0.0097^2) = 0.0509 -> 0.051

    Deliberately NOT item 23's 5.6% standard error, which is how well the
    CORRELATION fits its own data. That is a property of the model and so
    cannot live in this field.
    """
    bands = {s.uncertainty for s, _ in scored if s.label.startswith("florschuetz1981/")}
    assert bands == {0.051}, f"expected a single measurement band, got {bands}"

    expected = math.sqrt(0.05**2 + 0.0097**2)
    assert pytest.approx(expected, abs=5e-4) == 0.051


def test_han2012_scored_bands_are_the_sources_stated_value(scored) -> None:
    """The second half of the repair, pinned.

    The structural test above cannot see these three: once unscored series
    are dropped, each is alone in its (figure, quantity) group, so there is
    nothing for its band to be inconsistent WITH. They were found by reading
    provenance, and a pin is what holds them -- a numeric "band equals this
    series' own RMSE" rule fires on coincidence too (with a constant 0.051
    over 27 florschuetz series, two land on their own RMSE by chance), so it
    would be a guard that cries wolf rather than one that binds.

    What was wrong, and what each now carries:

      fig4.47_R_vs_alpha    0.105 -> 0.06    its own RMS against
                                             han_park_1988_angled; the
                                             cross_check beside it already
                                             said "Han's stated band is 6%"
      fig4.47_G_vs_eplus    0.088 -> 0.06    its own RMS
      fig4.193c_G_scatter   0.069 -> null    han_ribbed_high_re.md's own
                                             "Eq. (18) against the cloud,
                                             RMS 6.9%"

    0.06 is Han and Park's published band (han_ribbed.md item 7) and what
    every other scored series in this source carries. fig4.193c gets null
    rather than 0.06 because it is a different primary paper -- Rallabandi,
    Yang and Han (2009) -- whose own band this dataset has not read, and a
    band borrowed from another rig is not a measurement either. `within`
    renders "-" for it, which is the honest report.
    """
    bands = {s.label: s.uncertainty for s, _ in scored if s.label.startswith("han2012/")}
    assert bands, "no han2012 series scored"

    assert bands.get("han2012/fig4.47_R_vs_alpha") == 0.06
    assert bands.get("han2012/fig4.47_G_vs_eplus") == 0.06
    assert "han2012/fig4.193c_G_scatter" not in bands or (
        bands["han2012/fig4.193c_G_scatter"] is None
    )

    # And nothing else in this source invented a band of its own.
    stray = {k: v for k, v in bands.items() if v is not None and v != 0.06}
    assert not stray, (
        f"scored han2012 series with a band other than the source's stated 6%: {stray}"
    )


def test_within_reports_a_real_fraction_after_the_fix(scored) -> None:
    """The repair's effect, pinned.

    Against a fixed band the florschuetz pool fell from 58.7% to 47.5% and
    individual series now reach both 0.0% and 100.0%. The spread is the
    point, and it is what falsifies: with the self-derived band the range
    was 33.3% to 88.9%, because a series always sits near its own RMS and so
    can neither fall wholly outside a band built from it nor wholly inside.
    """
    flor = [
        c.within
        for s, c in scored
        if s.label.startswith("florschuetz1981/") and math.isfinite(c.within)
    ]
    assert len(flor) == 27, f"expected 27 scored series, got {len(flor)}"
    assert min(flor) == 0.0, "no series sits wholly outside a fixed band"
    assert max(flor) == 1.0, "no series sits wholly inside a fixed band"
