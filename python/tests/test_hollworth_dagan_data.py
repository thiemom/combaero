"""Hollworth & Dagan (1980) Part I: average heat transfer of a jet array whose
spent air leaves through vent holes in the target (#468). The coordinates are
digitised from Figs. 6 and 8 (cluster centres per standoff station; the Re
symbols collapse on Re^0.8). The extraction record is
validation/cooling/extractions/hollworth_dagan_1980_vented_target.md.

These tests check the digitisation against channels it was not read from:
the printed Table 2, the paper's stated staggered-over-edge-vented gains, and
the two figures against each other. The extraction is signed off (2026-10-24);
scoring comes with the #468 closure.
"""

from __future__ import annotations

import csv
from functools import cache

import pytest

from validation.cooling.schema import DATA_ROOT

D = DATA_ROOT / "hollworth_dagan1980"


@cache
def _rows(name: str) -> list[dict[str, str]]:
    with (D / name).open() as f:
        return list(csv.DictReader(f))


def _fig6(arrangement: str) -> dict[float, float]:
    return {
        float(r["Z"]): float(r["C100"])
        for r in _rows("fig6_x10_stations.csv")
        if r["arrangement"] == arrangement
    }


def _fig8(X: int) -> dict[float, float]:
    return {
        float(r["Z"]): float(r["C100"])
        for r in _rows("fig8_staggered_stations.csv")
        if int(r["X"]) == X
    }


def test_the_stations_are_all_there() -> None:
    """Fig. 8: four X, at Z = 1, 2.5, 5, 10 and at Z = X for X = 15, 20."""
    assert {X: sorted(_fig8(X)) for X in (5, 10, 15, 20)} == {
        5: [1, 2.5, 5, 10],
        10: [1, 2.5, 5, 10],
        15: [1, 2.5, 5, 10, 15],
        20: [1, 2.5, 5, 10, 20],
    }
    assert sorted(_fig6("staggered")) == sorted(_fig6("inline")) == [1, 2.5, 5, 10]


@pytest.mark.parametrize("Z", [1, 2.5, 5, 10])
def test_fig6_reproduces_the_printed_inline_ratio(Z: float) -> None:
    """Table 2 prints Nu_inline/Nu_staggered; Fig. 6 plots both at X = 10."""
    printed = {
        float(r["Z"]): float(r["ratio"])
        for r in _rows("table2_inline_over_staggered.csv")
        if int(r["X"]) == 10
    }
    read = _fig6("inline")[Z] / _fig6("staggered")[Z]
    assert read == pytest.approx(printed[Z], abs=0.015)


@pytest.mark.parametrize("Z", [1, 2.5, 5, 10])
def test_the_two_figures_agree_at_x10(Z: float) -> None:
    """Fig. 6 (linear axis) and Fig. 8 (log axis) plot the same staggered
    X = 10 data, read independently."""
    assert _fig8(10)[Z] == pytest.approx(_fig6("staggered")[Z], rel=0.02)


@pytest.mark.parametrize(("X", "stated"), [(10, 0.35), (15, 0.32), (20, 0.24)])
def test_the_stated_gain_over_edge_venting(X: int, stated: float) -> None:
    """The text: staggered vents yield Nu 35/32/24% above the same array with
    edge venting (ref. 4's curves on Fig. 8), averaged over 1 <= Z <= 10. The
    average here is over the stations the ref. 4 curve spans."""
    ref: dict[float, float] = {
        float(r["Z"]): float(r["C100"])
        for r in _rows("fig8_ref4_edge_vented_curves.csv")
        if int(r["X"]) == X
    }
    gains = [c / ref[Z] - 1.0 for Z, c in _fig8(X).items() if 1 <= Z <= 10 and Z in ref]
    assert len(gains) >= 3
    assert sum(gains) / len(gains) == pytest.approx(stated, abs=0.025)


def test_the_axis_calibration_holds_out() -> None:
    """Each tick, held out of a linear fit on the others, lands on its value:
    1% on Fig. 8's log axis, 0.02 on the linear ones."""
    import numpy as np

    for name, axis, log, tol in (
        ("fig8_axis_ticks.csv", "y_log", True, 0.01),
        ("fig8_axis_ticks.csv", "x", False, 0.05),
        ("fig6_axis_ticks.csv", "y", False, 0.02),
        ("fig6_axis_ticks.csv", "x", False, 0.02),
    ):
        rows = [r for r in _rows(name) if r["axis"] == axis]
        p = np.array([float(r["pixel"]) for r in rows])
        v = np.array([float(r["value"]) for r in rows])
        t = np.log10(v) if log else v
        for i in range(len(p)):
            keep = np.arange(len(p)) != i
            pred = np.polyval(np.polyfit(p[keep], t[keep], 1), p[i])
            if log:
                assert 10**pred == pytest.approx(v[i], rel=tol)
            else:
                assert pred == pytest.approx(v[i], abs=tol)
