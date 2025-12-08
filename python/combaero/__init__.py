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
        molar_volume,
        viscosity,
        thermal_conductivity,
        kinematic_viscosity,
        prandtl,
        thermal_diffusivity,
        reynolds,
        peclet,
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
        formula_to_name,
        name_to_formula,
        common_name,
        formula,
        calc_T_from_h,
        calc_T_from_s,
        calc_T_from_cp,
        oxygen_required_per_mol_fuel,
        oxygen_required_per_kg_fuel,
        oxygen_required_per_mol_mixture,
        oxygen_required_per_kg_mixture,
        complete_combustion_to_CO2_H2O,
        complete_combustion_to_CO2_H2O_with_fraction,
        # State-based API
        State,
        Stream,
        mix,
        # Inverse solvers - find fuel stream
        set_fuel_stream_for_Tad,
        set_fuel_stream_for_O2,
        set_fuel_stream_for_O2_dry,
        set_fuel_stream_for_CO2,
        set_fuel_stream_for_CO2_dry,
        # Inverse solvers - find oxidizer stream
        set_oxidizer_stream_for_Tad,
        set_oxidizer_stream_for_O2,
        set_oxidizer_stream_for_O2_dry,
        set_oxidizer_stream_for_CO2,
        set_oxidizer_stream_for_CO2_dry,
        complete_combustion,
        complete_combustion_isothermal,
        wgs_equilibrium,
        wgs_equilibrium_adiabatic,
        smr_wgs_equilibrium,
        smr_wgs_equilibrium_adiabatic,
        reforming_equilibrium,
        reforming_equilibrium_adiabatic,
        combustion_equilibrium,
        # Compressible flow
        CompressibleFlowSolution,
        NozzleStation,
        NozzleSolution,
        FannoStation,
        FannoSolution,
        nozzle_flow,
        solve_A_eff_from_mdot,
        solve_P_back_from_mdot,
        solve_P0_from_mdot,
        critical_pressure_ratio,
        mach_from_pressure_ratio,
        mass_flux_isentropic,
        nozzle_cd,
        fanno_pipe,
        fanno_max_length,
        # Thrust
        ThrustResult,
        nozzle_thrust,
        nozzle_thrust_cd,
        # Friction
        friction_haaland,
        friction_serghides,
        friction_colebrook,
        friction_petukhov,
        # Heat transfer
        nusselt_dittus_boelter,
        nusselt_gnielinski,
        nusselt_sieder_tate,
        nusselt_petukhov,
        htc_from_nusselt,
        overall_htc,
        overall_htc_wall,
        overall_htc_wall_multilayer,
        thermal_resistance,
        thermal_resistance_wall,
        lmtd,
        lmtd_counterflow,
        lmtd_parallelflow,
        # Incompressible flow
        bernoulli_P2,
        bernoulli_v2,
        orifice_mdot,
        orifice_Q,
        orifice_velocity,
        orifice_area,
        orifice_dP,
        pipe_dP,
        pipe_dP_mdot,
        pipe_velocity,
        pipe_mdot,
        dynamic_pressure,
        velocity_from_q,
        hydraulic_diameter,
        hydraulic_diameter_rect,
        hydraulic_diameter_annulus,
        # Orifice Cd correlations
        OrificeGeometry,
        OrificeState,
        Cd_sharp_thin_plate,
        Cd_thick_plate,
        Cd_rounded_entry,
        Cd_orifice,
        orifice_mdot_Cd,
        orifice_dP_Cd,
        orifice_Cd_from_measurement,
        orifice_K_from_Cd,
        orifice_Cd_from_K,
        orifice_thickness_correction,
        # Units query API
        UnitInfo,
        get_units,
        input_units,
        output_units,
        has_units,
        list_functions_with_units,
        all_units,
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
    molar_volume = _core.molar_volume
    viscosity = _core.viscosity
    thermal_conductivity = _core.thermal_conductivity
    kinematic_viscosity = _core.kinematic_viscosity
    prandtl = _core.prandtl
    thermal_diffusivity = _core.thermal_diffusivity
    reynolds = _core.reynolds
    peclet = _core.peclet
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
    formula_to_name = _core.formula_to_name
    name_to_formula = _core.name_to_formula
    common_name = _core.common_name
    formula = _core.formula
    calc_T_from_h = _core.calc_T_from_h
    calc_T_from_s = _core.calc_T_from_s
    calc_T_from_cp = _core.calc_T_from_cp
    oxygen_required_per_mol_fuel = _core.oxygen_required_per_mol_fuel
    oxygen_required_per_kg_fuel = _core.oxygen_required_per_kg_fuel
    oxygen_required_per_mol_mixture = _core.oxygen_required_per_mol_mixture
    oxygen_required_per_kg_mixture = _core.oxygen_required_per_kg_mixture
    complete_combustion_to_CO2_H2O = _core.complete_combustion_to_CO2_H2O
    complete_combustion_to_CO2_H2O_with_fraction = _core.complete_combustion_to_CO2_H2O_with_fraction
    # State-based API
    State = _core.State
    Stream = _core.Stream
    mix = _core.mix
    # Inverse solvers - find fuel stream
    set_fuel_stream_for_Tad = _core.set_fuel_stream_for_Tad
    set_fuel_stream_for_O2 = _core.set_fuel_stream_for_O2
    set_fuel_stream_for_O2_dry = _core.set_fuel_stream_for_O2_dry
    set_fuel_stream_for_CO2 = _core.set_fuel_stream_for_CO2
    set_fuel_stream_for_CO2_dry = _core.set_fuel_stream_for_CO2_dry
    # Inverse solvers - find oxidizer stream
    set_oxidizer_stream_for_Tad = _core.set_oxidizer_stream_for_Tad
    set_oxidizer_stream_for_O2 = _core.set_oxidizer_stream_for_O2
    set_oxidizer_stream_for_O2_dry = _core.set_oxidizer_stream_for_O2_dry
    set_oxidizer_stream_for_CO2 = _core.set_oxidizer_stream_for_CO2
    set_oxidizer_stream_for_CO2_dry = _core.set_oxidizer_stream_for_CO2_dry
    complete_combustion = _core.complete_combustion
    complete_combustion_isothermal = _core.complete_combustion_isothermal
    wgs_equilibrium = _core.wgs_equilibrium
    wgs_equilibrium_adiabatic = _core.wgs_equilibrium_adiabatic
    smr_wgs_equilibrium = _core.smr_wgs_equilibrium
    smr_wgs_equilibrium_adiabatic = _core.smr_wgs_equilibrium_adiabatic
    reforming_equilibrium = _core.reforming_equilibrium
    reforming_equilibrium_adiabatic = _core.reforming_equilibrium_adiabatic
    combustion_equilibrium = _core.combustion_equilibrium
    # Compressible flow
    CompressibleFlowSolution = _core.CompressibleFlowSolution
    NozzleStation = _core.NozzleStation
    NozzleSolution = _core.NozzleSolution
    FannoStation = _core.FannoStation
    FannoSolution = _core.FannoSolution
    nozzle_flow = _core.nozzle_flow
    solve_A_eff_from_mdot = _core.solve_A_eff_from_mdot
    solve_P_back_from_mdot = _core.solve_P_back_from_mdot
    solve_P0_from_mdot = _core.solve_P0_from_mdot
    critical_pressure_ratio = _core.critical_pressure_ratio
    mach_from_pressure_ratio = _core.mach_from_pressure_ratio
    mass_flux_isentropic = _core.mass_flux_isentropic
    nozzle_cd = _core.nozzle_cd
    fanno_pipe = _core.fanno_pipe
    fanno_max_length = _core.fanno_max_length
    # Thrust
    ThrustResult = _core.ThrustResult
    nozzle_thrust = _core.nozzle_thrust
    nozzle_thrust_cd = _core.nozzle_thrust_cd
    # Friction
    friction_haaland = _core.friction_haaland
    friction_serghides = _core.friction_serghides
    friction_colebrook = _core.friction_colebrook
    friction_petukhov = _core.friction_petukhov
    # Heat transfer
    nusselt_dittus_boelter = _core.nusselt_dittus_boelter
    nusselt_gnielinski = _core.nusselt_gnielinski
    nusselt_sieder_tate = _core.nusselt_sieder_tate
    nusselt_petukhov = _core.nusselt_petukhov
    htc_from_nusselt = _core.htc_from_nusselt
    overall_htc = _core.overall_htc
    overall_htc_wall = _core.overall_htc_wall
    overall_htc_wall_multilayer = _core.overall_htc_wall_multilayer
    thermal_resistance = _core.thermal_resistance
    thermal_resistance_wall = _core.thermal_resistance_wall
    lmtd = _core.lmtd
    lmtd_counterflow = _core.lmtd_counterflow
    lmtd_parallelflow = _core.lmtd_parallelflow
    # Incompressible flow
    bernoulli_P2 = _core.bernoulli_P2
    bernoulli_v2 = _core.bernoulli_v2
    orifice_mdot = _core.orifice_mdot
    orifice_Q = _core.orifice_Q
    orifice_velocity = _core.orifice_velocity
    orifice_area = _core.orifice_area
    orifice_dP = _core.orifice_dP
    pipe_dP = _core.pipe_dP
    pipe_dP_mdot = _core.pipe_dP_mdot
    pipe_velocity = _core.pipe_velocity
    pipe_mdot = _core.pipe_mdot
    dynamic_pressure = _core.dynamic_pressure
    velocity_from_q = _core.velocity_from_q
    hydraulic_diameter = _core.hydraulic_diameter
    hydraulic_diameter_rect = _core.hydraulic_diameter_rect
    hydraulic_diameter_annulus = _core.hydraulic_diameter_annulus
    # Orifice Cd correlations
    OrificeGeometry = _core.OrificeGeometry
    OrificeState = _core.OrificeState
    Cd_sharp_thin_plate = _core.Cd_sharp_thin_plate
    Cd_thick_plate = _core.Cd_thick_plate
    Cd_rounded_entry = _core.Cd_rounded_entry
    Cd_orifice = _core.Cd_orifice
    orifice_mdot_Cd = _core.orifice_mdot_Cd
    orifice_dP_Cd = _core.orifice_dP_Cd
    orifice_Cd_from_measurement = _core.orifice_Cd_from_measurement
    orifice_K_from_Cd = _core.orifice_K_from_Cd
    orifice_Cd_from_K = _core.orifice_Cd_from_K
    orifice_thickness_correction = _core.orifice_thickness_correction
    # Units query API
    UnitInfo = _core.UnitInfo
    get_units = _core.get_units
    input_units = _core.input_units
    output_units = _core.output_units
    has_units = _core.has_units
    list_functions_with_units = _core.list_functions_with_units
    all_units = _core.all_units


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
    "molar_volume",
    "viscosity",
    "thermal_conductivity",
    "kinematic_viscosity",
    "prandtl",
    "thermal_diffusivity",
    "reynolds",
    "peclet",
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
    "formula_to_name",
    "name_to_formula",
    "common_name",
    "formula",
    "calc_T_from_h",
    "calc_T_from_s",
    "calc_T_from_cp",
    "oxygen_required_per_mol_fuel",
    "oxygen_required_per_kg_fuel",
    "oxygen_required_per_mol_mixture",
    "oxygen_required_per_kg_mixture",
    "complete_combustion_to_CO2_H2O",
    "complete_combustion_to_CO2_H2O_with_fraction",
    # State-based API
    "State",
    "Stream",
    "mix",
    # Inverse solvers - find fuel stream
    "set_fuel_stream_for_Tad",
    "set_fuel_stream_for_O2",
    "set_fuel_stream_for_O2_dry",
    "set_fuel_stream_for_CO2",
    "set_fuel_stream_for_CO2_dry",
    # Inverse solvers - find oxidizer stream
    "set_oxidizer_stream_for_Tad",
    "set_oxidizer_stream_for_O2",
    "set_oxidizer_stream_for_O2_dry",
    "set_oxidizer_stream_for_CO2",
    "set_oxidizer_stream_for_CO2_dry",
    "complete_combustion",
    "complete_combustion_isothermal",
    "wgs_equilibrium",
    "wgs_equilibrium_adiabatic",
    "smr_wgs_equilibrium",
    "smr_wgs_equilibrium_adiabatic",
    "reforming_equilibrium",
    "reforming_equilibrium_adiabatic",
    "combustion_equilibrium",
    # Compressible flow
    "CompressibleFlowSolution",
    "NozzleStation",
    "NozzleSolution",
    "FannoStation",
    "FannoSolution",
    "nozzle_flow",
    "solve_A_eff_from_mdot",
    "solve_P_back_from_mdot",
    "solve_P0_from_mdot",
    "critical_pressure_ratio",
    "mach_from_pressure_ratio",
    "mass_flux_isentropic",
    "nozzle_cd",
    "fanno_pipe",
    "fanno_max_length",
    # Thrust
    "ThrustResult",
    "nozzle_thrust",
    "nozzle_thrust_cd",
    # Friction
    "friction_haaland",
    "friction_serghides",
    "friction_colebrook",
    "friction_petukhov",
    # Heat transfer
    "nusselt_dittus_boelter",
    "nusselt_gnielinski",
    "nusselt_sieder_tate",
    "nusselt_petukhov",
    "htc_from_nusselt",
    "overall_htc",
    "overall_htc_wall",
    "overall_htc_wall_multilayer",
    "thermal_resistance",
    "thermal_resistance_wall",
    "lmtd",
    "lmtd_counterflow",
    "lmtd_parallelflow",
    # Incompressible flow
    "bernoulli_P2",
    "bernoulli_v2",
    "orifice_mdot",
    "orifice_Q",
    "orifice_velocity",
    "orifice_area",
    "orifice_dP",
    "pipe_dP",
    "pipe_dP_mdot",
    "pipe_velocity",
    "pipe_mdot",
    "dynamic_pressure",
    "velocity_from_q",
    "hydraulic_diameter",
    "hydraulic_diameter_rect",
    "hydraulic_diameter_annulus",
    # Orifice Cd correlations
    "OrificeGeometry",
    "OrificeState",
    "Cd_sharp_thin_plate",
    "Cd_thick_plate",
    "Cd_rounded_entry",
    "Cd_orifice",
    "orifice_mdot_Cd",
    "orifice_dP_Cd",
    "orifice_Cd_from_measurement",
    "orifice_K_from_Cd",
    "orifice_Cd_from_K",
    "orifice_thickness_correction",
    # Units query API
    "UnitInfo",
    "get_units",
    "input_units",
    "output_units",
    "has_units",
    "list_functions_with_units",
    "all_units",
    "__version__",
]
