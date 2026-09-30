"""Score effusion plate INTERNAL heat transfer against Andrews' Figure 8.

The chain under test is Andrews 86-GT-225 Eq. (19) -- Sparrow's hole
approach summed with Mills' short-hole throat -- called through the real
implementation, never reimplemented here. The only arithmetic this module
owns is the conversions the sources' own axes require.

TWO PAPERS PLOT "h" AGAINST "G" AND THEY DO NOT MEAN THE SAME h.
This is the single thing to get right here, and it is worth a factor of
`A/A_h` -- between 1.4 and 9.8 across the plates involved.

  * **88-GT-290 Fig. 8** (`y_axis: h_internal`) says "h Convective heat
    transfer coefficient based on the surface area, A", with "A Total hole
    approach surface area (A = X^2 - pi/4 D^2)". PLATE AREA.
  * **86-GT-225 Fig. 8** (`y_axis: h_hole_length`) says "h_m Average heat
    transfer coefficient, W/m2K, over the hole length". HOLE INTERNAL
    AREA, pi D L -- the same basis Eq. (19)'s Nusselt number already uses,
    so no area conversion at all.

Eq. (19) gives a Nusselt number on the hole diameter and the hole internal
area, so the chain is

    Re = 4 (G X^2) / (pi D mu)          coolant per hole from G
    h  = Nu k / D                       86-GT-225's h_m, directly
    h  = Nu k / D * A_h / A             88-GT-290's h, rebased on A

Using the plate-area form on 86-GT-225's points scores -60%; the correct
one scores -10.4%. That is a definitional error of the #389 class, not a
modelling one, and the two conventions are carried as distinct `y_axis`
values precisely so the runner cannot silently pick the wrong one.

WHAT A MISS MEANS -- and it differs between the two figures:

  * **86-GT-225 Fig. 8 is FIDELITY.** Same paper as the correlations, so a
    miss would normally be our transcription bug. It is not, here: Fig. 10
    of that paper pins our Eq. (19) against the author's own evaluation of
    Eq. (19) at 1.5%, so the -10.4% is the correlation missing its
    author's measurements.
  * **88-GT-290 Fig. 8 is ACCURACY.** Tempting to call it fidelity -- same
    author, same group -- but the correlations are 1986 and the data is
    1988: same lab, different study, which the validation policy counts as
    cross-source. The scorecard label follows from `SET_ORIGIN` naming the
    1986 paper rather than the bare author name, which would have matched
    the 1988 source and overstated the claim.

Neither paper states a coolant temperature, so `ASSUMED_T` is an
assumption; `temperature_sensitivity` reports what it costs rather than
hiding it.
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


# The two h conventions, named rather than spelled inline, because the
# whole point is that a reader must never confuse them. See the module
# docstring for the nomenclature each paper prints.
H_PLATE_AREA = "h_internal"  # 88-GT-290: h on A = X^2 - pi D^2/4
H_HOLE_LENGTH = "h_hole_length"  # 86-GT-225: h_m over the hole length, pi D L

SCORED_Y_AXES = (H_PLATE_AREA, H_HOLE_LENGTH)


def owns(series: SeriesMetadata) -> bool:
    """Whether this runner should score ``series``.

    Explicit predicate, because `run_series` returns reason-carrying records
    for series it cannot score, which makes truthiness useless for dispatch.
    """
    return series.x_axis == "G_coolant" and series.y_axis in SCORED_Y_AXES


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

    # Nu is already on the hole internal area, which is exactly
    # 86-GT-225's own h_m basis -- so that convention needs no conversion.
    h_hole = nu * k / D
    if series.y_axis == H_HOLE_LENGTH:
        return h_hole

    # 88-GT-290 references h to the plate approach area instead. Pure
    # geometry, no correlation involved.
    area_plate = X * X - math.pi * D * D / 4.0
    area_hole = math.pi * D * L
    return h_hole * area_hole / area_plate


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
