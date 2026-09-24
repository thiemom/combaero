"""Evaluate jet_array_impingement_nu against Florschuetz's own Figure 6.

Unlike every rib series in runner.py, this needs no Reynolds-number
bisection at all: Nu/Nu1 = 1 - B[(z/d)(Gc/Gj)]^n exactly, because Nu1 is
the same correlation's own value at Gc/Gj = 0 (bracket = 1), so the
A*Re_j^m*Pr^(1/3) factor cancels identically between Nu and Nu1. Calling
jet_array_impingement_nu twice at the same (arbitrary) Re_j and Pr -- once
at the point's own Gc/Gj, once at Gc/Gj = 0 -- and taking the ratio scores
the real implementation exactly, the same discipline runner.py uses
(never reimplement the formula here to compute an "expected" value).
"""

from __future__ import annotations

from dataclasses import dataclass

import combaero as cb
from validation.cooling.schema import Point, SeriesMetadata, load_points

SETS = {
    "florschuetz_1981_inline": cb.florschuetz_1981_inline,
    "florschuetz_1981_staggered": cb.florschuetz_1981_staggered,
}

# Re_j and Pr cancel in the Nu/Nu1 ratio (see module docstring), so any
# value works -- these just need to be finite and inside a plausible range
# for the extrapolation check below to mean anything on the Re_j axis.
ARBITRARY_RE_J = 20000.0
ARBITRARY_PR = 0.7


@dataclass(frozen=True)
class Record:
    """One digitised point, evaluated."""

    series: SeriesMetadata
    x: float  # Gc/Gj
    measured: float  # Nu/Nu1
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


def _nu_ratio(
    jet_set: cb.JetArrayCorrelationSet, Gc_Gj: float, xn_d: float, yn_d: float, z_d: float
) -> tuple[float, bool]:
    at_point = cb.jet_array_impingement_nu(
        jet_set, ARBITRARY_RE_J, Gc_Gj, ARBITRARY_PR, xn_d, yn_d, z_d
    )
    at_zero = cb.jet_array_impingement_nu(
        jet_set, ARBITRARY_RE_J, 0.0, ARBITRARY_PR, xn_d, yn_d, z_d
    )
    ratio = at_point.Nu / at_zero.Nu if at_zero.Nu else float("nan")
    # Re_j's own extrapolation flag is meaningless here (ARBITRARY_RE_J is
    # not a real measurement) -- only the geometry/Gc_Gj flags from either
    # call carry information, and both calls share the same geometry.
    return ratio, at_point.extrapolated


def owns(series: SeriesMetadata) -> bool:
    """Whether this runner is the one that should score ``series``.

    Explicit, because ``run_series`` deliberately returns reason-carrying
    records rather than nothing for series it cannot score -- useful when it
    is driven directly, but useless as an ownership test.
    """
    return series.x_axis == "Gc_Gj" and series.y_axis == "Nu_over_Nu1"


def run_series(series: SeriesMetadata) -> list[Record]:
    """Evaluate one series. Unscored series yield records with no prediction."""
    points: list[Point] = load_points(series)

    if series.scores is None:
        return [Record(series, p.x, p.y, None, False, "not scored by any set") for p in points]
    if series.x_axis != "Gc_Gj" or series.y_axis != "Nu_over_Nu1":
        return [
            Record(
                series,
                p.x,
                p.y,
                None,
                False,
                f"x_axis={series.x_axis}, y_axis={series.y_axis} not implemented",
            )
            for p in points
        ]
    if not series.geometry or not all(k in series.geometry for k in ("xn_d", "yn_d", "z_d")):
        return [
            Record(
                series,
                p.x,
                p.y,
                None,
                False,
                "series carries no xn_d/yn_d/z_d geometry",
            )
            for p in points
        ]

    jet_set = SETS[series.scores]()
    xn_d = float(series.geometry["xn_d"])
    yn_d = float(series.geometry["yn_d"])
    z_d = float(series.geometry["z_d"])

    records = []
    for p in points:
        predicted, extrapolated = _nu_ratio(jet_set, p.x, xn_d, yn_d, z_d)
        records.append(Record(series, p.x, p.y, predicted, extrapolated))
    return records


def run_all(dataset: list[SeriesMetadata]) -> list[Record]:
    out: list[Record] = []
    for series in dataset:
        out.extend(run_series(series))
    return out
