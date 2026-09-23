"""Score mcgreehan_schotsch_1988_cd against Rohde's data, as the source plots it.

See validation/cooling/extractions/orifice_discharge_coefficient.md. The data
is McGreehan and Schotsch Fig. 6's open-square series -- Rohde's measurements
converted to static parameters by the paper's own authors, so no conversion of
ours sits between the measurement and the comparison.
"""

from __future__ import annotations

import pytest

from validation.cooling.orifice_runner import run_all, score
from validation.cooling.schema import load_dataset


@pytest.fixture(scope="module")
def dataset():
    return [s for s in load_dataset() if s.source.name == "mcgreehan_schotsch1988"]


def test_the_series_loads(dataset) -> None:
    assert len(dataset) == 1
    s = dataset[0]
    assert s.geometry == {"r_over_d": 0.0, "L_over_d": 0.51, "cd_baseline": 0.64}
    assert s.kind == "measured"


def test_eq17_against_rohde(dataset) -> None:
    """Eq. (17) anchored at the figure's own baseline, which is what Fig. 6 shows.

    The correlation reads HIGH against Rohde and the band below is not
    centred on zero: the paper itself notes Rohde's values run low. The band
    records the measured disagreement; it is not a target that was widened to
    admit the result.
    """
    s = score(run_all(dataset, mode="eq17"))
    assert s.n == 11
    assert 0.03 <= s.bias <= 0.09, f"bias drifted: {s.bias:+.4f}"
    assert s.rms <= 0.09, f"RMS drifted: {s.rms:.4f}"


def test_the_metric_responds_to_the_crossflow_term(dataset) -> None:
    """Falsify the metric, not just the change.

    If perturbing the crossflow input leaves the score alone, the score is
    measuring an imposed quantity and proves nothing. Perturb U1/Vi and the
    RMS must degrade in both directions.
    """
    from validation.cooling.orifice_runner import Record, predict

    base = run_all(dataset, mode="eq17")
    base_rms = score(base).rms

    for scale in (0.5, 2.0):
        perturbed = [
            Record(
                series=r.series,
                x=r.x,
                measured=r.measured,
                predicted=predict(r.series, r.x * scale, "eq17"),
            )
            for r in base
        ]
        assert score(perturbed).rms > 2.0 * base_rms, (
            f"scaling U1/Vi by {scale} barely moved the score: "
            f"{score(perturbed).rms:.4f} vs {base_rms:.4f}"
        )


def test_chain_mode_carries_the_baseline_gap(dataset) -> None:
    """The unanchored chain is worse, and that is the documented reason Fig. 6
    anchors its baseline rather than predicting it."""
    eq17 = score(run_all(dataset, mode="eq17"))
    chain = score(run_all(dataset, mode="chain"))
    assert chain.bias > eq17.bias
    assert chain.rms > eq17.rms
    # Rohde's own basic value for this geometry is 0.64; the chain gives ~0.669.
    assert 0.10 <= chain.bias <= 0.18, f"baseline gap drifted: {chain.bias:+.4f}"


@pytest.fixture(scope="module")
def rohde():
    return [s for s in load_dataset() if s.source.name == "rohde1969"]


def test_rohde_fig10_loads(rohde) -> None:
    assert len(rohde) == 3
    by_rd = {s.geometry["r_over_d"]: s for s in rohde}
    assert sorted(by_rd) == pytest.approx([0.0, 0.1951, 0.4878])
    for s in rohde:
        assert s.geometry["L_over_d"] == 1.06
        assert s.x_axis == "velocity_head_ratio"


def test_rohde_ordering_invariant(rohde) -> None:
    """More inlet rounding must give a higher Cd at every velocity head ratio.

    Physically required, and independent of the correlation entirely. A
    swapped series label or a calibration slip in one series breaks it.
    """
    import math

    from validation.cooling.schema import load_points

    curves = {
        s.geometry["r_over_d"]: sorted((math.log(p.x), p.y) for p in load_points(s)) for s in rohde
    }

    def at(curve, lx):
        if lx < curve[0][0] or lx > curve[-1][0]:
            return None
        for (x0, y0), (x1, y1) in zip(curve, curve[1:], strict=False):
            if x0 <= lx <= x1:
                return y0 + (lx - x0) / (x1 - x0) * (y1 - y0)
        return None

    checked = 0
    for vhr in (4, 6, 8, 11, 15, 20, 27, 35, 45):
        lx = math.log(vhr)
        sharp, mid, round_ = (at(curves[0.0], lx), at(curves[0.1951], lx), at(curves[0.4878], lx))
        if None in (sharp, mid, round_):
            continue
        checked += 1
        assert round_ > mid > sharp, f"ordering broken at VHR {vhr}"
    assert checked >= 8


def test_rohde_reveals_the_crossflow_limit_of_the_rd_term(rohde) -> None:
    """The measured disagreement, recorded rather than tuned away.

    The correlation reads HIGH against Rohde and the error grows sharply as
    crossflow does. Critically, it grows MOST for the sharp orifice and least
    for the most-rounded one -- an ordering the total-to-static conversion
    cannot produce, since that factor is identical for all three series at a
    given velocity head ratio. So this is a real limit of the model, not an
    artefact of getting Rohde into the correlation's coordinates.
    """
    from validation.cooling.orifice_runner import run_series, score

    low = {}
    high = {}
    for s in rohde:
        recs = run_series(s)
        rd = s.geometry["r_over_d"]
        low[rd] = score([r for r in recs if r.vhr < 10.0])
        high[rd] = score([r for r in recs if r.vhr >= 25.0])

    # At low crossflow all three converge on the same modest overprediction,
    # matching the independent Fig. 6 result (+5.95%).
    for rd, sc in high.items():
        assert 0.0 < sc.bias < 0.12, f"r/d={rd} high-VHR bias {sc.bias:+.3f}"

    # At high crossflow the error explodes, and monotonically in r/d.
    assert low[0.0].bias > low[0.1951].bias > low[0.4878].bias
    assert low[0.0].bias > 0.30, f"sharp low-VHR bias {low[0.0].bias:+.3f}"


def test_rohde_deskew_is_derived_from_committed_calibration(rohde) -> None:
    """The y-skew correction must come from the corners file, not a constant.

    Confirmed against a third channel: Rohde's prose says the largest radius
    reaches "as high as 0.94". The raw circle maximum is 0.957; deskewed it is
    0.936.
    """
    from validation.cooling.orifice_runner import apply_deskew, deskew
    from validation.cooling.schema import load_points

    circles = next(s for s in rohde if s.geometry["r_over_d"] == 0.4878)
    dy_lo, dy_hi = deskew(circles)
    assert dy_hi - dy_lo == pytest.approx(0.0205, abs=0.002)
    assert abs(dy_lo) < 0.003, "left edge should be nearly true"

    pts = load_points(circles)
    raw_max = max(p.y for p in pts)
    x_at = next(p.x for p in pts if p.y == raw_max)
    assert raw_max == pytest.approx(0.957, abs=0.002)
    assert apply_deskew(x_at, raw_max, dy_lo, dy_hi) == pytest.approx(0.94, abs=0.01)
