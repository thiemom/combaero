"""Score effusion plate INTERNAL heat transfer against Andrews' Figure 8.

The chain under test is Andrews 86-GT-225 Eq. (19) -- Sparrow's hole
approach summed with Mills' short-hole throat -- called through the real
implementation, never reimplemented here. The only arithmetic this module
owns is the two conversions the source's own axes require.

WHAT THE SOURCE PLOTS. Andrews 88-GT-290 Fig. 8 is `h` against `G`, where
`G` is coolant mass flow per unit plate area and `h` is referenced to the
plate approach area `A = X^2 - pi D^2/4`. Eq. (19) gives a Nusselt number
on the hole diameter and the HOLE INTERNAL area, so two conversions stand
between them:

    Re = 4 (G X^2) / (pi D mu)          coolant per hole from G
    h  = Nu k / D * A_h / A             hole-area Nu to plate-area h

`A_h/A` is 1/3.46 for plate C, so omitting it would overstate `h` by that
factor -- a definitional error of the #389 class, not a modelling one.

WHAT A MISS MEANS. Tempting to call this fidelity -- same author, same
group -- but it is ACCURACY. The correlations are Andrews 86-GT-225 and the
data is Andrews 88-GT-290: same lab, different study, which the validation
policy counts as cross-source. The scorecard label follows from
`SET_ORIGIN` naming the 1986 paper rather than the bare author name, which
would have matched the 1988 source and overstated the claim.

The paper states no coolant temperature, so `ASSUMED_T` is an assumption;
`temperature_sensitivity` reports what it costs rather than hiding it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import combaero as cb

from validation.cooling.schema import SeriesMetadata, load_points

# The rig is a transient technique with laboratory-temperature coolant and
# the paper states no value, so this is an ASSUMPTION. 280-320 K brackets
# any plausible reading and moves the result by about 4% end to end --
# reported by `temperature_sensitivity`, never tuned.
ASSUMED_T = 300.0
ASSUMED_P = 101325.0


@dataclass(frozen=True)
class Record:
    """One digitised point, evaluated."""

    series: SeriesMetadata
    x: float  # G [kg/s/m^2]
    measured: float  # h [W/m^2 K]
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
    """Whether this runner should score ``series``.

    Explicit predicate, because `run_series` returns reason-carrying records
    for series it cannot score, which makes truthiness useless for dispatch.
    """
    return series.x_axis == "G_coolant" and series.y_axis == "h_internal"


def _geometry(series: SeriesMetadata) -> tuple[float, float, float] | None:
    """(D, X, L) in metres, or None if the series does not carry them.

    The hole LENGTH is not the plate thickness in general -- an inclined
    hole is longer by 1/sin(alpha). Andrews' plates are drilled normal to
    the surface (alpha = 90 deg) so the two coincide here, but computing it
    rather than assuming it keeps the runner honest for an angled plate.
    """
    g = series.geometry or {}
    if not all(k in g for k in ("D_mm", "X_mm", "thickness_mm")):
        return None
    alpha = math.radians(float(series.alpha_deg or 90.0))
    if not (0.0 < alpha <= math.pi / 2 + 1e-12):
        return None
    length = float(g["thickness_mm"]) * 1e-3 / math.sin(alpha)
    return float(g["D_mm"]) * 1e-3, float(g["X_mm"]) * 1e-3, length


def predict(
    series: SeriesMetadata, G: float, T: float = ASSUMED_T
) -> float | None:
    """Plate-area internal heat transfer coefficient at coolant flux G."""
    geom = _geometry(series)
    if geom is None or G <= 0.0:
        return None
    D, X, L = geom

    air = cb.standard_dry_air_composition()
    k = cb.thermal_conductivity(T, ASSUMED_P, air)
    mu = cb.viscosity(T, ASSUMED_P, air)
    pr = cb.prandtl(T, ASSUMED_P, air)

    # G is per unit GROSS plate area, so the flow through one hole is G
    # times the cell area X^2. Confirmed against the source's own hole
    # count: N = 4306 per m^2 and 1/X^2 = 4305.8.
    m_hole = G * X * X
    Re = 4.0 * m_hole / (math.pi * D * mu)

    nu = cb.effusion_internal_nusselt(Re, pr, X / L, L / D)

    # Nu is on the hole internal area; Fig. 8's h is on the plate approach
    # area. Both are geometry, no correlation involved.
    area_plate = X * X - math.pi * D * D / 4.0
    area_hole = math.pi * D * L
    return nu * k / D * area_hole / area_plate


def run_series(series: SeriesMetadata) -> list[Record]:
    points = load_points(series)

    if series.scores is None:
        return [
            Record(series, p.x, p.y, None, False, "not scored by any set")
            for p in points
        ]
    if not owns(series):
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
    if _geometry(series) is None:
        return [
            Record(series, p.x, p.y, None, False, "geometry lacks D_mm/X_mm/thickness_mm")
            for p in points
        ]

    return [
        Record(series, p.x, p.y, predict(series, p.x), False) for p in points
    ]


def temperature_sensitivity(
    series: SeriesMetadata, levels: tuple[float, ...] = (280.0, 300.0, 320.0)
) -> dict[float, float]:
    """Mean relative error against assumed coolant temperature.

    The paper states none. Reported, never minimised: choosing the
    temperature that scores best would be fitting an unmeasured input to
    the metric that judges it.
    """
    points = load_points(series)
    out: dict[float, float] = {}
    for T in levels:
        errs = []
        for p in points:
            v = predict(series, p.x, T)
            if v is not None and p.y:
                errs.append(v / p.y - 1.0)
        if errs:
            out[T] = sum(errs) / len(errs)
    return out


def run_all(dataset: list[SeriesMetadata]) -> list[Record]:
    out: list[Record] = []
    for series in dataset:
        out.extend(run_series(series))
    return out
