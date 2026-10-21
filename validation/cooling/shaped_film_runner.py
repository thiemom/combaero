"""Score shaped-hole film effectiveness: Colban, Thole and Bogard (2011).

One row of shaped holes, laterally averaged adiabatic eta against x/D,
through `film_effectiveness_colban_2011` -- the real implementation, never
restated here. Point by point: unlike the multi-row series `film_runner`
reduces to per-row means, a single row has no sawtooth to phase.

The colban2011 series are Gritsch et al. (2005) data the correlation was
FITTED on, so their scores are fidelity on the fitting data -- an optimistic
basis, and labelled as such by the scorecard's provenance match.

OWNS EVERY colban2011 SERIES, and runs before `film_runner`, which would
otherwise claim any `x_over_D` / `eta_adiabatic` series and score it as
cylindrical-hole Baldauf + Sellers.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import combaero as cb

from validation.cooling.schema import SeriesMetadata, load_points

SET = "colban_2011_shaped"
NEEDED = ("M", "area_ratio", "P_over_D", "t_over_P", "angle_deg", "density_ratio")


@dataclass(frozen=True)
class Record:
    """One digitised point, evaluated."""

    series: SeriesMetadata
    x: float  # x/D
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
    """Every series of the colban2011 source, scored or not."""
    return series.source.name == "colban2011" or series.scores == SET


@lru_cache(maxsize=None)
def _is_extrapolated(
    M: float, AR: float, PD: float, tP: float, alpha_deg: float, DR: float
) -> bool:
    """Ask the correlation, via its warnings, rather than restating Table 4.

    Cached for the reason `film_runner._is_extrapolated` gives: each probe
    swaps the global warning handler, and restoring nests a wrapper.
    """
    seen: list[str] = []
    previous = cb.get_warning_handler()
    try:
        cb.set_warning_handler(seen.append)
        cb.film_effectiveness_colban_2011(10.0, M, AR, PD, tP, alpha_deg, DR)
    finally:
        cb.set_warning_handler(previous)
    return any("outside the fitted range" in m for m in seen)


def predict(series: SeriesMetadata, x_over_D: float) -> tuple[float, bool]:
    g = series.geometry or {}
    args = tuple(float(g[k]) for k in NEEDED)
    M, AR, PD, tP, alpha_deg, DR = args
    eta = cb.film_effectiveness_colban_2011(x_over_D, M, AR, PD, tP, alpha_deg, DR)
    return eta, _is_extrapolated(*args)


def run_series(series: SeriesMetadata) -> list[Record]:
    points = load_points(series)
    if series.scores is None:
        return [Record(series, p.x, p.y, None, False, "not scored by any set") for p in points]
    if series.scores != SET or series.x_axis != "x_over_D":
        return [
            Record(series, p.x, p.y, None, False, f"scores={series.scores} not implemented")
            for p in points
        ]
    missing = [k for k in NEEDED if k not in (series.geometry or {})]
    if missing:
        return [
            Record(series, p.x, p.y, None, False, f"geometry lacks {', '.join(missing)}")
            for p in points
        ]
    out: list[Record] = []
    for p in points:
        eta, extrapolated = predict(series, p.x)
        out.append(Record(series, p.x, p.y, eta, extrapolated))
    return out
