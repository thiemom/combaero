"""What each measured quantity is REFERENCED TO, declared rather than implied.

THE FAILURE THIS EXISTS TO CATCH. Two numbers for the same physical thing
under different conventions are different numbers, and every coordinate is
in range, so no bounds check can see it. #389 calls it a category error
rather than a disagreement, and this project has now hit it three times:

  * **Rohde's `Cd` is referenced to duct TOTAL pressure; McGreehan and
    Schotsch's to STATIC.** Related by `sqrt(VHR/(VHR-1))` -- 1.33 at
    VHR 2. `orifice_runner` applies it correctly and always has, but
    nothing DECLARED that a conversion was needed, so nothing would have
    noticed a future series arriving the other way round.
  * **Andrews 86-GT-225's `h_m` is over the HOLE LENGTH; 88-GT-290's `h`
    is on the PLATE AREA.** Same author, same group, same figure number.
    Taking one for the other scores -60% instead of -10%; it was caught
    by reading both nomenclatures by hand (#430).
  * **A gas-side `h` is meaningless without its unblown baseline.** The
    augmentation Andrews' Fig. 10 requires is 3.1-4.3 against a fully
    developed Dittus-Boelter and 1.6-2.5 against the same duct
    entry-corrected -- a factor of 1.75 from the baseline alone (#431).

THE SHAPE OF THE FIX. A measured series declares what its y quantity is
referenced to; a correlation set declares what it produces; scoring either
matches them, applies a REGISTERED conversion, or refuses. It never
silently compares two different quantities.

WHY THE DECLARATION IS A TABLE AND NOT 103 YAML EDITS. The convention is a
property of (source, quantity), not of the individual series -- every
Rohde `Cd` is total-referenced because that is how Rohde's rig was
instrumented. One auditable table beats the same string copied 103 times
and drifting. `SERIES_OVERRIDE` exists for the exception, and
`test_every_scored_series_declares_a_convention` makes a new source fail
loudly rather than default.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

# ---------------------------------------------------------------------------
# The conventions themselves. The string is the declaration; the comment is
# what a reader needs to know to tell two of them apart.
# ---------------------------------------------------------------------------

CONVENTION_NOTES: dict[str, str] = {
    # Discharge coefficient
    "Cd_total": "Cd referenced to duct TOTAL pressure (Rohde 1969)",
    "Cd_static": "Cd referenced to duct STATIC pressure (McGreehan-Schotsch)",
    # Rib roughness functions
    "G_ribbed_wall": "heat-transfer roughness function on the RIBBED wall only",
    "G_four_wall": "G_bar, averaged over two ribbed and two smooth walls",
    "R_ribbed_wall": "friction roughness function on the ribbed wall",
    "R_normalised_pe": "R divided by (P/e/10)^0.35, Han's figure 4.46 ordinate",
    "f_channel_two_ribbed": (
        "passage-average FANNING friction of a channel with two opposite "
        "ribbed walls -- Han's measured fbar, NOT the four-sided f_r = "
        "fbar + (H/W)(fbar - f_s) his R is built on"
    ),
    "R_normalised_angled": "R normalised for the angled-rib quadratic form",
    "Nu_ratio_ribbed_DB_vs_fbar": (
        "RIBBED-side Nu over Dittus-Boelter 0.023 Re^0.8 Pr^0.4, plotted "
        "against CHANNEL-AVERAGE f over 0.046 Re^-0.2 (Han, Zhang & Lee 1991, "
        "Eqs. 2/4). Not the four-sided f: 1.8x apart at e/D 0.0625. Nu is "
        "averaged over X/D 0-20 including the entrance region"
    ),
    # Jet impingement
    "Nu_over_Nu1": "row Nusselt divided by the first-row value, same correlation",
    # Film and effusion
    "eta_adiabatic": (
        "ADIABATIC film effectiveness, (Tg-Taw)/(Tg-Tc). Measured on an "
        "insulated wall, so it carries no heat transfer coefficient"
    ),
    "eta_overall": (
        "OVERALL cooling effectiveness, (Tg-Tw)/(Tg-Tc). Contains the "
        "internal convection and the gas-side coefficient as well as the "
        "film -- NOT the same quantity as eta_adiabatic and not convertible "
        "to it without the whole resistance network"
    ),
    "h_hole_internal_area": "h referenced to the hole internal surface, pi D L",
    "h_plate_approach_area": "h referenced to the plate area, X^2 - pi D^2/4",
}

# ---------------------------------------------------------------------------
# What each SOURCE publishes, per quantity. Keyed (source, y_axis).
# ---------------------------------------------------------------------------

SOURCE_PUBLISHES: dict[tuple[str, str], str] = {
    ("rohde1969", "Cd"): "Cd_total",
    ("mcgreehan_schotsch1988", "Cd"): "Cd_static",
    ("han2012", "G"): "G_ribbed_wall",
    ("han_park_lei1984", "G"): "G_ribbed_wall",
    ("lau1990", "G"): "G_ribbed_wall",
    ("han2012", "G_bar"): "G_four_wall",
    ("lau1990", "G_bar"): "G_four_wall",
    ("han_park_lei1984", "R"): "R_ribbed_wall",
    # Taslim & Spring (1987) Fig. 11, f = dP gc / (2 (L/D_H) rho V^2) over
    # the passage -- the same definition as Han's Eq. (1) (#403 item 2).
    ("taslim_spring1987", "f_fanning_passage"): "f_channel_two_ribbed",
    ("lau1990", "R"): "R_ribbed_wall",
    # Figure 4.48's R panel (#402) plots raw R, unlike figures 4.46/4.47's
    # normalised ordinate -- same quantity as han_park_lei1984/lau1990's "R".
    ("han2012", "R"): "R_ribbed_wall",
    ("han2012", "R_normalised"): "R_normalised_pe",
    # Figure 4.53, Han and Zhang (1992). Their own normalisation is not on
    # disk; the 1991 paper it reprints from prints it (#435).
    ("han2012", "Nu_ratio"): "Nu_ratio_ribbed_DB_vs_fbar",
    ("han2012", "R_normalised_angled"): "R_normalised_angled",
    ("florschuetz1981", "Nu_over_Nu1"): "Nu_over_Nu1",
    ("andrei2014", "eta_adiabatic"): "eta_adiabatic",
    ("murray2018", "eta_adiabatic"): "eta_adiabatic",
    ("andrews1988", "eta_overall"): "eta_overall",
    # The pair that cost -60% when confused. Same lab, same figure number.
    ("andrews1986", "h_hole_length"): "h_hole_internal_area",
    ("andrews1988", "h_internal"): "h_plate_approach_area",
}

# A series whose convention differs from the rest of its source's. None yet;
# the hook exists so the exception does not force the table to be abandoned.
SERIES_OVERRIDE: dict[str, str] = {}

# ---------------------------------------------------------------------------
# What each CORRELATION SET produces.
# ---------------------------------------------------------------------------

SET_PRODUCES: dict[str, str] = {
    "han_1988_orthogonal": "G_ribbed_wall",
    "han_park_1988_angled": "G_ribbed_wall",
    "han_1989_narrow_channel": "G_ribbed_wall",
    "rallabandi_2009_high_re": "G_ribbed_wall",
    "florschuetz_1981_inline": "Nu_over_Nu1",
    "mcgreehan_schotsch_1988_cd": "Cd_static",
    "mcgreehan_schotsch_1988_crossflow_cd": "Cd_static",
    "baldauf_2002_sellers": "eta_adiabatic",
    # One set, two quantities: it produces an internal h, and the overall
    # eta runner builds a network result from it. Declared on the quantity
    # the runner actually compares, which `resolve` is told per series.
    "andrews_1986_effusion_internal": "h_hole_internal_area",
}

# A set that legitimately produces more than one quantity declares the
# others here, so `resolve` can accept any of them rather than needing a
# conversion that is really a different output.
SET_ALSO_PRODUCES: dict[str, tuple[str, ...]] = {
    # The rib sets produce R as well as G, in several normalisations the
    # runner selects between; and G_bar through Han's published 1.2.
    "han_1988_orthogonal": (
        "R_ribbed_wall", "R_normalised_pe", "G_four_wall", "f_channel_two_ribbed",
        "Nu_ratio_ribbed_DB_vs_fbar",
    ),
    "han_park_1988_angled": (
        "R_ribbed_wall", "R_normalised_pe", "R_normalised_angled",
        "G_four_wall",
    ),
    "rallabandi_2009_high_re": ("R_ribbed_wall", "R_normalised_pe"),
    # Fig. 4.48's R panel plots raw R (no (p/e/10)^0.35 normalisation unlike
    # figures 4.46/4.47), so this set needs no R_normalised_* entry.
    "han_1989_narrow_channel": ("R_ribbed_wall",),
    # The effusion set's overall-eta closure is a network result built from
    # the same internal h -- an OUTPUT of the set, not a conversion of it.
    "andrews_1986_effusion_internal": (
        "h_plate_approach_area", "eta_overall",
    ),
}

# Han, Zhang and Lee (1991) Table 2: R and G on Han's basis (four-sided f_r,
# ribbed-side St), plus the source's OWN printed G_bar; the 60 deg V set also
# answers figure 4.53's performance curve through the f_ratio path.
for _tag in ("90", "60par", "60crs", "60vee", "60lam", "45par", "45crs", "45vee", "45lam"):
    SET_PRODUCES[f"han_zhang_lee_1991_{_tag}"] = "G_ribbed_wall"
    SET_ALSO_PRODUCES[f"han_zhang_lee_1991_{_tag}"] = (
        "R_ribbed_wall", "G_four_wall", "Nu_ratio_ribbed_DB_vs_fbar",
    )


@dataclass(frozen=True)
class Conversion:
    """A registered way to turn one convention into another.

    `factor` takes whatever the runner can supply as context and returns
    the multiplier. Registered conversions are PURE DEFINITION -- geometry
    or an algebraic identity -- never a fitted correction. A conversion
    that needed fitting would be a model, and would belong in the library
    with its own evidence.
    """

    derivation: str
    applied_by: str
    factor: Callable[..., float] | None = None


CONVERSIONS: dict[tuple[str, str], Conversion] = {
    # Rohde's Cd on total -> the static basis McGreehan-Schotsch produces.
    # The factor diverges as VHR -> 1, which is why those scores are
    # reported per VHR band and never pooled.
    ("Cd_total", "Cd_static"): Conversion(
        derivation="Cd_static = Cd_total * sqrt(VHR/(VHR-1)); 1.33 at VHR 2, "
        "1.01 at VHR 50. See extractions/rohde_1969_orifice.md.",
        applied_by="validation.cooling.orifice_runner.vhr_to_static_cd",
        factor=lambda vhr: math.sqrt(vhr / (vhr - 1.0)),
    ),
    # Andrews' two h bases. Pure geometry: A_h/A with A_h = pi D L and
    # A = X^2 - pi D^2/4. Between 1.4 and 9.8 over his plates.
    ("h_hole_internal_area", "h_plate_approach_area"): Conversion(
        derivation="h_plate = h_hole * (pi D L)/(X^2 - pi D^2/4). Pure "
        "geometry, no correlation. See extractions/"
        "andrews_effusion_internal_h.md.",
        applied_by="validation.cooling.effusion_internal_runner.predict",
        factor=None,  # needs the series geometry; the runner owns it
    ),
}

# Pairs that must NEVER be converted, with why. Listing them is the point:
# an unregistered pair is refused by default, but these are the ones a
# reader might assume are convertible.
INCOMPATIBLE: dict[tuple[str, str], str] = {
    ("eta_adiabatic", "eta_overall"): (
        "An adiabatic effectiveness is measured on a wall passing no heat, "
        "so it contains no coefficient. An overall effectiveness contains "
        "the internal convection AND the gas-side coefficient. Going from "
        "one to the other needs the whole resistance network, which is a "
        "model and not a conversion -- and an overall-eta correlation used "
        "as a closure would double-count the convection the network "
        "computes (#387)."
    ),
    ("G_ribbed_wall", "G_four_wall"): (
        "Han's G_bar = 1.2 G is a published CONSTANT fitted at 90 degrees, "
        "not an identity: the measured ratio spans 1.096 to 1.413 across "
        "16 configurations, with 69% of the variance between rigs (#401). "
        "It is applied as the source's own correlation, with its cost "
        "measured -- not registered here as though it were definitional."
    ),
}

STATUS_DIRECT = "direct"
STATUS_CONVERTED = "converted"
STATUS_UNDECLARED = "undeclared"
STATUS_INCOMPATIBLE = "incompatible"
STATUS_UNREGISTERED = "unregistered"
# A series no set scores. Distinct from `undeclared`, which means a
# scored series whose convention the table does not yet carry -- the
# count that must not grow.
STATUS_NOT_SCORED = "not-scored"


@dataclass(frozen=True)
class Resolution:
    status: str
    detail: str
    conversion: Conversion | None = None

    @property
    def scorable(self) -> bool:
        """Whether a comparison may proceed.

        `undeclared` scores: the table is incomplete by construction while
        sources are being added, and refusing would hide working results
        behind bookkeeping. It is COUNTED instead, and
        `test_the_undeclared_count_does_not_grow` is what keeps that from
        becoming permanent.
        """
        return self.status in (
            STATUS_DIRECT,
            STATUS_CONVERTED,
            STATUS_UNDECLARED,
            STATUS_NOT_SCORED,
        )


def series_convention(series) -> str | None:
    """What this series' y quantity is referenced to, or None if undeclared."""
    override = SERIES_OVERRIDE.get(series.label)
    if override:
        return override
    return SOURCE_PUBLISHES.get((series.source.name, series.y_axis))


def resolve(series, set_name: str | None) -> Resolution:
    """Whether `series` and `set_name` measure the same quantity.

    Returns a Resolution rather than raising, so the scorecard can REPORT
    a refusal beside the rows it refused -- a silent drop is the failure
    mode this is replacing.
    """
    if not set_name:
        return Resolution(STATUS_NOT_SCORED, "not scored by any set")

    measured = series_convention(series)
    produced = SET_PRODUCES.get(set_name)
    if measured is None or produced is None:
        missing = []
        if measured is None:
            missing.append(f"series ({series.source.name}, {series.y_axis})")
        if produced is None:
            missing.append(f"set ({set_name})")
        return Resolution(
            STATUS_UNDECLARED, "no convention declared for " + " and ".join(missing)
        )

    if measured == produced or measured in SET_ALSO_PRODUCES.get(set_name or "", ()):
        return Resolution(STATUS_DIRECT, measured)

    why = INCOMPATIBLE.get((produced, measured)) or INCOMPATIBLE.get(
        (measured, produced)
    )
    if why:
        return Resolution(
            STATUS_INCOMPATIBLE, f"{produced} vs {measured}: {why}"
        )

    conversion = CONVERSIONS.get((measured, produced))
    if conversion is not None:
        return Resolution(
            STATUS_CONVERTED, f"{measured} -> {produced}", conversion
        )

    return Resolution(
        STATUS_UNREGISTERED,
        f"{produced} vs {measured}: no registered conversion. Scoring two "
        "different quantities is a category error, not a disagreement -- "
        "register the conversion with its derivation, or declare the pair "
        "INCOMPATIBLE.",
    )
