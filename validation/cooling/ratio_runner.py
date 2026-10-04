"""Score ratio-form rib sets (#444) against digitised series.

A ratio set (``cb.RibRatioSet``) predicts Nu and f directly as multipliers on
its declared smooth-duct baseline, so the scored quantity is the paper's own:
no law-of-the-wall conversion on the measured side. Like the other runners it
never recomputes a correlation itself -- every prediction comes from
``cb.evaluate_rib_ratio``.
"""

from __future__ import annotations

from dataclasses import dataclass

import combaero as cb
from validation.cooling.schema import Point, SeriesMetadata, load_points

# Taslim & Spring tested air; their Dittus-Boelter smooth-duct check
# reproduces at Pr 0.70 (extractions/taslim_spring_1987_aspect_ratio.md), so
# the ratio's baseline is evaluated there too.
TASLIM_PR = 0.70

# Two-side configurations: (Taslim aspect ratio, e/D_H).
TASLIM_CONFIGS = (
    (0.5, 0.125), (0.5, 0.250), (1.0, 0.083), (1.0, 0.167),
    (3.5, 0.053), (3.5, 0.107), (3.5, 0.161),
)

SETS = {
    f"taslim_spring_1987_ar{ar:.1f}_eD{ed:.3f}": (
        lambda ar=ar, ed=ed: cb.taslim_spring_1987(ar, ed)
    )
    for ar, ed in TASLIM_CONFIGS
}

# What each y_axis is predicted by, and at which Prandtl number.
PREDICTS = {
    "Nu_turbulated": ("Nu", TASLIM_PR),
    "f_fanning_passage": ("f", TASLIM_PR),
}


@dataclass(frozen=True)
class Record:
    """One digitised point, evaluated."""

    series: SeriesMetadata
    x: float  # Re
    measured: float
    predicted: float | None
    extrapolated: bool
    reason: str | None = None

    @property
    def rel_error(self) -> float | None:
        if self.predicted is None or self.measured == 0.0:
            return None
        return self.predicted / self.measured - 1.0

    @property
    def within_uncertainty(self) -> bool:
        band = self.series.uncertainty
        err = self.rel_error
        if band is None or err is None:
            return False
        return abs(err) <= band


def owns(series: SeriesMetadata) -> bool:
    """A series is this runner's exactly when a ratio set scores it."""
    return series.scores in SETS


def run_series(series: SeriesMetadata) -> list[Record]:
    points: list[Point] = load_points(series)

    def refuse(reason: str) -> list[Record]:
        return [Record(series, p.x, p.y, None, False, reason) for p in points]

    if series.scores not in SETS:
        return refuse("not scored by a ratio set")
    if series.x_axis != "Re_Dh" or series.y_axis not in PREDICTS:
        return refuse(f"x_axis={series.x_axis}, y_axis={series.y_axis} not implemented")
    geom_in = series.geometry or {}
    if not all(k in geom_in for k in ("e_D", "p_e", "W_H")):
        return refuse("series carries no e_D/p_e/W_H geometry")
    if geom_in.get("turbulated_walls", 2) != 2:
        return refuse("one-side-turbulated; the set is for two turbulated walls")

    ratio_set = SETS[series.scores]()
    geom = cb.RibGeometry(
        e_D=float(geom_in["e_D"]),
        p_e=float(geom_in["p_e"]),
        W_H=float(geom_in["W_H"]),
        alpha_deg=float(series.alpha_deg if series.alpha_deg is not None else 90.0),
    )
    quantity, pr = PREDICTS[series.y_axis]
    records = []
    for p in points:
        res = cb.evaluate_rib_ratio(ratio_set, geom, p.x, pr)
        predicted = res.Nu if quantity == "Nu" else res.f
        records.append(Record(series, p.x, p.y, predicted, res.extrapolated))
    return records


def run_all(dataset: list[SeriesMetadata]) -> list[Record]:
    out: list[Record] = []
    for series in dataset:
        out.extend(run_series(series))
    return out
