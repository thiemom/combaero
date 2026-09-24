"""Score the McGreehan and Schotsch Cd chain against Rohde's data.

Kept separate from runner.py and jet_array_runner.py for the same reason
those are separate from each other: there is no Reynolds-number bisection
and no correlation-set dispatch here. A point is (U1/Vi, Cd) and the
prediction is one call.

The baseline subtlety this runner exists to respect. McGreehan and Schotsch
Fig. 6 does NOT plot chain predictions: it anchors each curve to "a set
baseline point at U1/Vi = 0" taken from Rohde's own measurement, and applies
only Eq. (17) to it. The chain predicts a baseline 3-12% higher for the same
geometry, which is the paper's own observation that Rohde's "basic values are
lower". Scoring the full chain against this data would therefore measure the
sum of two disagreements and attribute both to the crossflow term.

So EQ17 is the scoring mode that reproduces the source's own validation, and
CHAIN is offered alongside it to show what the baseline gap costs -- never
as a substitute.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from typing import Literal

import combaero as cb

from validation.cooling.schema import SeriesMetadata, load_points

Mode = Literal["eq17", "chain"]

# The paper's own reference Reynolds number (p.213), used only by CHAIN mode.
# Eq. (8) is flat to 0.4% between Re = 1e4 and 1e6, so this choice is not
# what drives the CHAIN-mode disagreement -- the baseline gap is.
REFERENCE_RE = 3.2e4


# ---------------------------------------------------------------------------
# Rohde Fig. 10 lives in the REPORT's coordinates, not the correlation's.
# ---------------------------------------------------------------------------


def deskew(series: SeriesMetadata) -> "tuple[float, float]":
    """Measure the canvas y-skew from the committed corners file.

    Derived from data rather than hard-coded, so the correction is
    reproducible and moves if the calibration picks are ever redone. Returns
    (dy at the left edge, dy at the right edge), both in Cd.
    """
    path = series.path.parent / "fig10_corners.csv"
    pts: list[tuple[float, float]] = []
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if not row or row[0].strip() in ("", "x"):
                continue
            pts.append((float(row[0]), float(row[1])))
    # Corner order as digitised: TL, TR, BR, BL; nominal box y = [0.2, 1.0].
    (_, tl), (_, tr), (_, br), (_, bl) = pts
    return ((tl - 1.0) + (bl - 0.2)) / 2.0, ((tr - 1.0) + (br - 0.2)) / 2.0


VHR_AXIS_LO, VHR_AXIS_HI = 1.0, 60.0


def apply_deskew(vhr: float, cd: float, dy_lo: float, dy_hi: float) -> float:
    frac = math.log(vhr / VHR_AXIS_LO) / math.log(VHR_AXIS_HI / VHR_AXIS_LO)
    return cd - (dy_lo + (dy_hi - dy_lo) * frac)


def vhr_to_crossflow_ratio(vhr: float) -> float:
    """U1/Vi = 1/sqrt(VHR - 1). Derived in the extraction."""
    return 1.0 / math.sqrt(vhr - 1.0)


def vhr_to_static_cd(vhr: float, cd_total_referenced: float) -> float:
    """Rohde references his ideal flow to duct TOTAL pressure; Eq. (17) wants
    it referenced to duct STATIC. The factor diverges as VHR -> 1, which is
    why scores are reported per VHR band and never pooled."""
    return cd_total_referenced * math.sqrt(vhr / (vhr - 1.0))


# Rohde scores are reported per velocity-head-ratio band, never pooled: the
# static-referencing factor sqrt(VHR/(VHR-1)) is 1.33 at VHR 2 and 1.01 at
# VHR 50, so one number would mix a near-exact comparison with one dominated
# by the conversion. The bands travel with the records so the scorecard can
# keep them apart without knowing why.
VHR_BANDS = ((1.0, 10.0), (10.0, 25.0), (25.0, 60.0))


def _band_label(vhr: float) -> str | None:
    for lo, hi in VHR_BANDS:
        if lo <= vhr < hi:
            return f"VHR {lo:g}-{hi:g}"
    return None


@dataclass(frozen=True)
class Record:
    """One digitised point, evaluated.

    Carries the same reporting surface as the rib and jet-array runners
    (``rel_error``, ``within_uncertainty``, ``extrapolated``, ``reason``) so
    one scorecard can aggregate all three without special-casing.
    """

    series: SeriesMetadata
    x: float  # U1/Vi
    measured: float  # Cd
    predicted: float | None
    vhr: float | None = None  # set for Rohde series, else None
    extrapolated: bool = False
    reason: str | None = None
    #: Sub-partition for reporting. Rohde series are split by VHR band.
    group: str | None = None
    #: Which correlation actually produced the prediction.
    scored_by: str | None = None

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


def predict(series: SeriesMetadata, x: float, mode: Mode) -> float:
    geom = series.geometry or {}
    if mode == "eq17" and "cd_baseline" in geom:
        return cb.mcgreehan_schotsch_1988_crossflow_cd(geom["cd_baseline"], x)
    return cb.mcgreehan_schotsch_1988_cd(
        REFERENCE_RE, geom["r_over_d"], geom["L_over_d"], x
    )


def score_by_vhr_band(records: list[Record], bands) -> "dict":
    """Rohde scores are reported per VHR band, never pooled.

    The static-referencing factor sqrt(VHR/(VHR-1)) is 1.33 at VHR 2 and 1.01
    at VHR 50, so a pooled number would mix a near-exact comparison with one
    dominated by the conversion.
    """
    out = {}
    for lo, hi in bands:
        sel = [r for r in records if r.vhr is not None and lo <= r.vhr < hi]
        out[(lo, hi)] = score(sel)
    return out


def owns(series: SeriesMetadata) -> bool:
    """Whether this runner is the one that should score ``series``."""
    return series.y_axis == "Cd" and series.x_axis in (
        "U1_over_Vi",
        "velocity_head_ratio",
    )


def run_series(series: SeriesMetadata, mode: Mode = "eq17") -> list[Record]:
    if not owns(series):
        return []

    if series.x_axis == "U1_over_Vi":
        return [
            Record(series=series, x=p.x, measured=p.y,
                   predicted=predict(series, p.x, mode),
                   scored_by=series.scores)
            for p in load_points(series)
        ]

    if series.x_axis == "velocity_head_ratio":
        # Rohde's own coordinates. Deskew, then convert both axes.
        dy_lo, dy_hi = deskew(series)
        out = []
        for p in load_points(series):
            if p.x <= 1.02:  # the conversion is meaningless at VHR -> 1
                continue
            cd = vhr_to_static_cd(p.x, apply_deskew(p.x, p.y, dy_lo, dy_hi))
            u = vhr_to_crossflow_ratio(p.x)
            out.append(
                Record(series=series, x=u, measured=cd,
                       predicted=predict(series, u, "chain"), vhr=p.x,
                       group=_band_label(p.x), scored_by=series.scores)
            )
        return out

    return []


def run_all(dataset: list[SeriesMetadata], mode: Mode = "eq17") -> list[Record]:
    out: list[Record] = []
    for s in dataset:
        out.extend(run_series(s, mode))
    return out


@dataclass(frozen=True)
class Score:
    n: int
    bias: float  # fractional, not percent
    rms: float
    mae: float


def score(records: list[Record]) -> Score:
    errs = [(r.predicted - r.measured) / r.measured for r in records]
    n = len(errs)
    if n == 0:
        return Score(0, math.nan, math.nan, math.nan)
    return Score(
        n=n,
        bias=sum(errs) / n,
        rms=math.sqrt(sum(e * e for e in errs) / n),
        mae=sum(abs(e) for e in errs) / n,
    )
