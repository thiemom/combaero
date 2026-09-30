"""Score OVERALL effusion cooling effectiveness against Andrews' Figure 10.

`eta = (T_g - T_w)/(T_g - T_c)`, the paper's Eq. (3). Unlike the internal
`h` runner beside this one, eta is not a correlation output -- it is a
NETWORK result, and the point of scoring it is to find out whether the
network is the right one.

WHAT IS SHIPPED, AND WHAT IT LUMPS

    eta = h_i / (h_i + h_gas)

`h_gas` is the caller's gas side. The closure carries NO separate film
term, and the reason is not that there is no film. It is that the film
and the gas-side augmentation that comes with it are two halves of one
pair, and only the pair is identifiable from this data -- so both are
lumped into the one `h_gas` the caller supplies.

WHY HALF A PAIR IS NOT ENOUGH. The two-temperature form is

    q = h_f (T_aw - T_w),   T_aw = T_g - eta_f (T_g - T_c)
    eta = (h_i + h_f eta_f) / (h_i + h_f)

and it needs both `eta_f` and `h_f`. Baldauf supplies `eta_f` only, and
structurally cannot supply `h_f`: an adiabatic wall has `q = 0` by
construction, so the measurement determines `T_aw` and says nothing about
the coefficient.

WHAT BALDAUF DOES AND DOES NOT CONTAIN -- worth being exact about,
because an earlier reading of this had it backwards. The counter-rotating
vortex pair entraining hot gas under a lifted jet changes `T_aw`, which
is precisely what an IR-on-insulated-wall measurement records. Baldauf
HAS that, in Eq. (38)'s decay branch, and steeply: the falling exponent
`b_pk = exp(1.92 - 7.5 (s/D)^-1.5)` is 4.57 at plate B's spacing and 3.24
at plate C's. Andrews' own sentence about "jet stirring of the cooling
film" describes that same film destruction -- so it names the mechanism
Baldauf models, not a missing one.

SO WHAT THE MEASUREMENT SAYS, BOTH WAYS:

  * Offered at the smooth-duct `h_g`, the film is refused. The admissible
    `eta_f` is negative at all 66 of plate B's points and reaches at most
    +0.105 on plate C, against 0.27 to 0.58 from Baldauf over ten rows.
  * Admitted with a matching augmentation, it is consistent: the same
    measurement then demands `F = h_f/h_0` of 1.6 to 4.3.

THE ABSOLUTE LEVEL IS CONFOUNDED AND THE RATIO IS NOT. The test plate
sits only `x/D_h = 1.5` into the duct, deep in the entry region where
Dittus-Boelter's fully developed value understates `h_0`; a standard
entry correction raises it 1.75x and brings plate C's F to 0.9-1.6,
ordinary film-cooling augmentation. That factor divides out of every B/C
ratio, which is why the finding is stated as a ratio:

    B/C augmentation ratio     no film    with Baldauf's film
      G = 0.6                    1.82           1.75
      G = 1.0                    1.86           1.67
      G = 1.4                    1.96           1.79

Plate B needs about 1.7 times plate C's gas-side coefficient either way.
ADMITTING THE FILM DOES NOT EXPLAIN THE SPLIT, which is the result: the
residual is a geometry dependence that two plates cannot separate from a
jet one. Three candidates survive and none is chosen -- Baldauf
extrapolating hard (M to 7.0 against 2.5, s/D 7.06 against 5), a genuine
difference in `h_f/h_0`, and Sellers over ten rows, which Gao et al.
(2025) measure as always overestimating and accumulating with row count.

NO THRESHOLD IS OFFERED. The enhancement collapses on none of the
velocity ratio, the blowing ratio or the momentum flux ratio: at equal
values the plates still differ by well over 25%. So `jet_regime()`
REPORTS all three beside the result and fits nothing -- the same
report-don't-enforce idiom the correlation status codes use. A user CAN
reach plate B's regime without meaning to, and there the closure misses
by 24%.

THE GAS SIDE IS THE RIG, NOT THE MODEL. `h_g` here is Dittus-Boelter on
the duct the paper describes: 76 x 152 mm, product gases at 750 K and
M = 0.05. That is a boundary condition of Andrews' facility, so it lives
in this runner rather than in the library, and `gas_side_sensitivity()`
reports what it is worth. It is never chosen to improve the score.

WHAT A MISS MEANS: ACCURACY. The correlation is 86-GT-225 and the data is
88-GT-290 -- same lab, different study, which the validation policy counts
as cross-source.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import combaero as cb

from validation.cooling.schema import SeriesMetadata, load_points

# The facility, page 3: "the wall of a 76 mm by 152 mm wide air cooled duct
# through which the product gases from a propane preheater flowed at 750K
# and a Mach number of 0.05".
DUCT_H_M = 0.076
DUCT_W_M = 0.152
GAS_MACH = 0.05
# Combustion products, so not air; 1.33 is the usual hot-products value and
# only enters the speed of sound, hence U_g, hence h_g to the 0.8 power.
GAS_GAMMA = 1.33
R_AIR = 287.0

# Page 7: "The measurements were made at a Tg of 750K and Tc of 295K". The
# PLOT prints 744 K and 293 K. The text is used and the figure's values are
# reported by `temperature_sensitivity`; eta is a ratio, so it is worth
# under 0.5% either way.
T_GAS = 750.0
T_COOLANT = 295.0
T_GAS_PLOTTED = 744.0
T_COOLANT_PLOTTED = 293.0
P_AMBIENT = 101325.0


@dataclass(frozen=True)
class Record:
    """One digitised point, evaluated."""

    series: SeriesMetadata
    x: float  # G [kg/s/m^2]
    measured: float  # eta [-]
    predicted: float | None
    extrapolated: bool
    reason: str | None = None
    # Reported PER PLATE, so the rollup cannot average a 3% fit with a 24%
    # miss into a 13% number that describes neither. Per plate and NOT per
    # jet regime, deliberately -- see `_group_of`.
    group: str | None = None

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
    """Whether this runner should score ``series``."""
    return series.x_axis == "G_coolant" and series.y_axis == "eta_overall"


def _group_of(series: SeriesMetadata) -> str:
    """The rollup reporting group: one per plate.

    PER PLATE AND NOT PER JET REGIME, and the difference is not cosmetic.
    Splitting the points by velocity ratio instead looks informative --
    +6.9% below VR = 1 against +19.1% above it -- but it is an artefact of
    which plate supplies the points, because both plates cross VR = 1
    inside their own G range. Holding the regime fixed and varying the
    plate:

        VR < 1    plate B +21.2% (n=17)   plate C +2.2% (n=52)
        VR >= 1   plate B +24.6% (n=49)   plate C +5.6% (n=20)

    The plate moves the error by a factor of ten; the regime label moves it
    by a few points. A VR-grouped row would report the plate mix and call
    it physics. Pinned by
    `test_a_velocity_ratio_split_would_misreport_the_plate_effect`.
    """
    plate = (series.geometry or {}).get("plate")
    return f"overall eta, plate {plate}" if plate else "overall eta"


def _geometry(series: SeriesMetadata) -> tuple[float, float, float] | None:
    """(D, X, L) in metres, or None if the series does not carry them."""
    g = series.geometry or {}
    if not all(k in g for k in ("D_mm", "X_mm", "thickness_mm")):
        return None
    alpha = math.radians(float(series.alpha_deg or 90.0))
    if not (0.0 < alpha <= math.pi / 2 + 1e-12):
        return None
    length = float(g["thickness_mm"]) * 1e-3 / math.sin(alpha)
    return float(g["D_mm"]) * 1e-3, float(g["X_mm"]) * 1e-3, length


def gas_side() -> dict[str, float]:
    """The rig's gas side: duct velocity, density and heat transfer.

    Dittus-Boelter on the duct hydraulic diameter. An ASSUMPTION about the
    facility, not a combaero model -- reported so a reader can substitute
    their own, and never chosen to improve the score.
    """
    air = cb.standard_dry_air_composition()
    area = DUCT_H_M * DUCT_W_M
    perimeter = 2.0 * (DUCT_H_M + DUCT_W_M)
    d_h = 4.0 * area / perimeter
    u_g = GAS_MACH * math.sqrt(GAS_GAMMA * R_AIR * T_GAS)
    rho_g = P_AMBIENT / (R_AIR * T_GAS)
    k = cb.thermal_conductivity(T_GAS, P_AMBIENT, air)
    mu = cb.viscosity(T_GAS, P_AMBIENT, air)
    pr = cb.prandtl(T_GAS, P_AMBIENT, air)
    re = rho_g * u_g * d_h / mu
    return {
        "D_h": d_h,
        "U_g": u_g,
        "rho_g": rho_g,
        "Re_duct": re,
        "h_g": 0.023 * re**0.8 * pr**0.4 * k / d_h,
    }


def internal_h(series: SeriesMetadata, G: float) -> float | None:
    """Internal coefficient on the PLATE area, from 86-GT-225.

    The same chain `effusion_internal_runner` scores directly against
    Figure 8, called through the library rather than reimplemented.
    """
    geom = _geometry(series)
    if geom is None or G <= 0.0:
        return None
    D, X, L = geom
    air = cb.standard_dry_air_composition()
    k = cb.thermal_conductivity(T_COOLANT, P_AMBIENT, air)
    mu = cb.viscosity(T_COOLANT, P_AMBIENT, air)
    pr = cb.prandtl(T_COOLANT, P_AMBIENT, air)
    re = 4.0 * (G * X * X) / (math.pi * D * mu)
    nu = cb.effusion_internal_nusselt(re, pr, X / L, L / D)
    return nu * k / D * (math.pi * D * L) / (X * X - math.pi * D * D / 4.0)


def predict(series: SeriesMetadata, G: float, h_g: float | None = None) -> float | None:
    """Overall effectiveness from the two-resistance closure.

        eta = h_i / (h_i + h_g)

    NO separate film term and NO coolant heat-up term, both deliberately:

    * The film -- LUMPED INTO `h_g`, not denied. `eta_f` and `h_f` are
      two halves of one pair and Baldauf supplies only the first; offered
      alone at the smooth-duct `h_g` the data refuses it, and admitted
      with its matching augmentation the same data demands F = 1.6 to
      4.3. Only the pair is identifiable here, so the caller supplies it
      as one number. Zero free parameters, which is what makes this form
      falsifiable where the split form can fit anything.
    * The heat-up -- the coolant gains 12 K at G = 1.4 and 60 K at
      G = 0.2, so it looks obviously missing. It is not. 86-GT-225's h is
      fitted from the plate's transient cooling rate against the SUPPLY
      temperature, so the heat-up is inside it already. Adding it tilts
      the error rather than flattening it, which is the evidence:
      `test_the_coolant_heat_up_is_already_inside_the_measured_h`.

    Wall conduction is omitted for the same reason -- the transient
    technique lumps the plate -- and is 3.5% of `1/h_i` in any case.
    """
    h_i = internal_h(series, G)
    if h_i is None:
        return None
    h_gas = gas_side()["h_g"] if h_g is None else h_g
    return h_i / (h_i + h_gas)


def jet_regime(series: SeriesMetadata, G: float) -> dict[str, float] | None:
    """The three jet parameters, REPORTED and never used to correct.

    A user can reach plate B's regime -- jets faster than the mainstream --
    without meaning to, and there the closure misses by 24%. So the numbers
    that describe it travel with the result.

    None of the three explains the B/C split on its own: at equal VR the
    two plates differ by 1.5-1.6x in implied gas-side h, at equal M by
    1.35-1.45x, at equal I by 1.3-1.5x. Nor does admitting Baldauf's film,
    which leaves the ratio at 1.67-1.79. Two plates cannot separate a jet
    parameter from a geometry one, so nothing here is a threshold.
    """
    geom = _geometry(series)
    if geom is None or G <= 0.0:
        return None
    D, X, _ = geom
    gas = gas_side()
    rho_c = P_AMBIENT / (R_AIR * T_COOLANT)
    mass_flux_c = G * X * X / (math.pi * D * D / 4.0)
    u_jet = mass_flux_c / rho_c
    density_ratio = rho_c / gas["rho_g"]
    blowing = mass_flux_c / (gas["rho_g"] * gas["U_g"])
    return {
        "u_jet": u_jet,
        "U_gas": gas["U_g"],
        "velocity_ratio": u_jet / gas["U_g"],
        "blowing_ratio": blowing,
        "momentum_flux_ratio": blowing * blowing / density_ratio,
        "density_ratio": density_ratio,
    }


# The 152 mm plate carries ten hole rows at the 15.24 mm pitch, and a
# point at row n sees n of them upstream.
FILM_ROWS = 10


def film_effectiveness(series: SeriesMetadata, G: float) -> float | None:
    """Baldauf per row, Sellers-superposed, averaged over the plate.

    NOT used by `predict` -- see its docstring. This exists so that the
    alternative closure can be evaluated and reported rather than merely
    asserted about, and so `implied_gas_side_h` can be given a film to
    work against.

    Far outside Baldauf's envelope here: M reaches 7.0 against a 2.5
    limit, the density ratio is 2.54 against 1.8, and plate B's s/D is
    7.06 against 5. Sellers over ten rows is its own overestimate --
    Gao et al. (2025) measure that the error accumulates with row count
    and always overestimates.
    """
    geom = _geometry(series)
    if geom is None or G <= 0.0:
        return None
    D, X, _ = geom
    gas = gas_side()
    rho_c = P_AMBIENT / (R_AIR * T_COOLANT)
    blowing = (G * X * X / (math.pi * D * D / 4.0)) / (gas["rho_g"] * gas["U_g"])
    density_ratio = rho_c / gas["rho_g"]

    total = 0.0
    for n in range(1, FILM_ROWS + 1):
        rows = [
            cb.film_effectiveness_baldauf_2002(
                j * X / D, blowing, density_ratio, 90.0, X / D, 0.05
            )
            for j in range(1, n + 1)
        ]
        total += cb.film_superposition_sellers(rows)
    return total / FILM_ROWS


def implied_gas_side_h(
    series: SeriesMetadata, G: float, eta: float, eta_film: float = 0.0
) -> float | None:
    """The gas-side h this measurement implies, given combaero's internal h.

    Inverts the two-temperature form,
    `eta = (h_i + h_f eta_film)/(h_i + h_f)`, so

        h_f = h_i (1 - eta) / (eta - eta_film)

    With `eta_film = 0` this is the inversion of `predict` and depends on
    NO assumption about the rig's gas side -- which is why the B/C
    comparison is stated in terms of it. With Baldauf's own film it
    returns the augmented coefficient the measurement then requires,
    F = 1.6 to 4.3 times the smooth duct.

    Returns None where the measurement is at or below the film's own
    effectiveness, which the two-temperature form cannot represent.
    """
    h_i = internal_h(series, G)
    if h_i is None or not 0.0 < eta < 1.0 or eta <= eta_film:
        return None
    return h_i * (1.0 - eta) / (eta - eta_film)


def run_series(series: SeriesMetadata) -> list[Record]:
    points = load_points(series)

    if series.scores is None:
        return [
            Record(series, p.x, p.y, None, False, "not scored by any set", None)
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
                None,
            )
            for p in points
        ]
    if _geometry(series) is None:
        return [
            Record(series, p.x, p.y, None, False, "geometry lacks D_mm/X_mm/thickness_mm", None)
            for p in points
        ]

    group = _group_of(series)
    return [
        Record(series, p.x, p.y, predict(series, p.x), False, None, group)
        for p in points
    ]


def gas_side_sensitivity(
    series: SeriesMetadata, factors: tuple[float, ...] = (0.75, 1.0, 1.5, 2.0)
) -> dict[float, float]:
    """Mean relative error against a multiplier on the duct's gas-side h.

    The gas side is an assumption about Andrews' facility, so its cost is
    reported. Reported, never minimised: picking the multiplier that scores
    best would be fitting an unmeasured input to the metric judging it --
    and on plate B that multiplier is exactly the jet enhancement the model
    is failing to predict, which would hide the finding.
    """
    base = gas_side()["h_g"]
    out: dict[float, float] = {}
    for f in factors:
        errs = []
        for p in load_points(series):
            v = predict(series, p.x, base * f)
            if v is not None and p.y:
                errs.append(v / p.y - 1.0)
        if errs:
            out[f] = sum(errs) / len(errs)
    return out


def run_all(dataset: list[SeriesMetadata]) -> list[Record]:
    out: list[Record] = []
    for series in dataset:
        out.extend(run_series(series))
    return out
