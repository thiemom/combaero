"""Public Python API for the combaero package.

This module re-exports selected functions from the compiled _core extension.
It is written to work both when running from the source tree and when
combaero is installed as a wheel.
"""

from __future__ import annotations

from importlib import import_module

try:  # Python 3.8+
    from importlib.metadata import PackageNotFoundError, version as _pkg_version
except ModuleNotFoundError:  # pragma: no cover - very old Python
    from importlib_metadata import (  # type: ignore[import-not-found]
        PackageNotFoundError,
        version as _pkg_version,
    )


def _load_version() -> str:
    """Return the installed distribution version, or a safe fallback.

    When running from a source tree without an installed wheel, the
    distribution metadata may not be available; in that case we expose a
    placeholder string instead of raising.
    """

    try:
        return _pkg_version("combaero")
    except PackageNotFoundError:
        return "0.0.0+local"


__version__: str = _load_version()


try:
    # Preferred: local extension in the same package (installed wheel or
    # in-tree build where _core was successfully built next to this file).
    from ._core import (  # type: ignore[attr-defined]
        mixture_h,
        adiabatic_T_wgs,
        cp,
        h,
        s,
        cv,
        density,
        specific_gas_constant,
        isentropic_expansion_coefficient,
        speed_of_sound,
        viscosity,
        thermal_conductivity,
        kinematic_viscosity,
        prandtl,
        mole_to_mass,
        mass_to_mole,
        normalize_fractions,
        convert_to_dry_fractions,
        equivalence_ratio_mole,
        set_equivalence_ratio_mole,
        equivalence_ratio_mass,
        set_equivalence_ratio_mass,
        bilger_stoich_mixture_fraction_mass,
        bilger_Z_from_equivalence_ratio_mass,
        equivalence_ratio_from_bilger_Z_mass,
        standard_dry_air_composition,
        humid_air_composition,
        dewpoint,
        relative_humidity_from_dewpoint,
        species_common_names,
        calc_T_from_h,
        calc_T_from_s,
        calc_T_from_cp,
        oxygen_required_per_mol_fuel,
        oxygen_required_per_kg_fuel,
        oxygen_required_per_mol_mixture,
        oxygen_required_per_kg_mixture,
        complete_combustion_to_CO2_H2O,
        complete_combustion_to_CO2_H2O_with_fraction,
    )
except ModuleNotFoundError:
    # Fallback: attempt to import from an installed combaero package that
    # already has _core available, then re-export the symbols.
    _core = import_module("combaero._core")
    mixture_h = _core.mixture_h
    adiabatic_T_wgs = _core.adiabatic_T_wgs
    cp = _core.cp
    h = _core.h
    s = _core.s
    cv = _core.cv
    density = _core.density
    specific_gas_constant = _core.specific_gas_constant
    isentropic_expansion_coefficient = _core.isentropic_expansion_coefficient
    speed_of_sound = _core.speed_of_sound
    viscosity = _core.viscosity
    thermal_conductivity = _core.thermal_conductivity
    kinematic_viscosity = _core.kinematic_viscosity
    prandtl = _core.prandtl
    mole_to_mass = _core.mole_to_mass
    mass_to_mole = _core.mass_to_mole
    normalize_fractions = _core.normalize_fractions
    convert_to_dry_fractions = _core.convert_to_dry_fractions
    equivalence_ratio_mole = _core.equivalence_ratio_mole
    set_equivalence_ratio_mole = _core.set_equivalence_ratio_mole
    equivalence_ratio_mass = _core.equivalence_ratio_mass
    set_equivalence_ratio_mass = _core.set_equivalence_ratio_mass
    bilger_stoich_mixture_fraction_mass = _core.bilger_stoich_mixture_fraction_mass
    bilger_Z_from_equivalence_ratio_mass = _core.bilger_Z_from_equivalence_ratio_mass
    equivalence_ratio_from_bilger_Z_mass = _core.equivalence_ratio_from_bilger_Z_mass
    standard_dry_air_composition = _core.standard_dry_air_composition
    humid_air_composition = _core.humid_air_composition
    dewpoint = _core.dewpoint
    relative_humidity_from_dewpoint = _core.relative_humidity_from_dewpoint
    species_common_names = _core.species_common_names
    calc_T_from_h = _core.calc_T_from_h
    calc_T_from_s = _core.calc_T_from_s
    calc_T_from_cp = _core.calc_T_from_cp
    oxygen_required_per_mol_fuel = _core.oxygen_required_per_mol_fuel
    oxygen_required_per_kg_fuel = _core.oxygen_required_per_kg_fuel
    oxygen_required_per_mol_mixture = _core.oxygen_required_per_mol_mixture
    oxygen_required_per_kg_mixture = _core.oxygen_required_per_kg_mixture
    complete_combustion_to_CO2_H2O = _core.complete_combustion_to_CO2_H2O
    complete_combustion_to_CO2_H2O_with_fraction = _core.complete_combustion_to_CO2_H2O_with_fraction


__all__ = [
    "mixture_h",
    "adiabatic_T_wgs",
    "cp",
    "h",
    "s",
    "cv",
    "density",
    "specific_gas_constant",
    "isentropic_expansion_coefficient",
    "speed_of_sound",
    "viscosity",
    "thermal_conductivity",
    "kinematic_viscosity",
    "prandtl",
    "mole_to_mass",
    "mass_to_mole",
    "normalize_fractions",
    "convert_to_dry_fractions",
    "equivalence_ratio_mole",
    "set_equivalence_ratio_mole",
    "equivalence_ratio_mass",
    "set_equivalence_ratio_mass",
    "bilger_stoich_mixture_fraction_mass",
    "bilger_Z_from_equivalence_ratio_mass",
    "equivalence_ratio_from_bilger_Z_mass",
    "standard_dry_air_composition",
    "humid_air_composition",
    "dewpoint",
    "relative_humidity_from_dewpoint",
    "species_common_names",
    "calc_T_from_h",
    "calc_T_from_s",
    "calc_T_from_cp",
    "oxygen_required_per_mol_fuel",
    "oxygen_required_per_kg_fuel",
    "oxygen_required_per_mol_mixture",
    "oxygen_required_per_kg_mixture",
    "complete_combustion_to_CO2_H2O",
    "complete_combustion_to_CO2_H2O_with_fraction",
    "__version__",
]
