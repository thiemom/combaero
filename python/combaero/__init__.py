"""Public Python API for the combaero package.

This module re-exports selected functions from the compiled _core extension.
It is written to work both when running from the source tree and when
combaero is installed as a wheel.
"""

from __future__ import annotations

from importlib import import_module

try:  # Python 3.8+
    from importlib.metadata import PackageNotFoundError
    from importlib.metadata import version as _pkg_version
except ModuleNotFoundError:  # pragma: no cover - very old Python
    from importlib_metadata import PackageNotFoundError
    from importlib_metadata import version as _pkg_version  # type: ignore[import-not-found]


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
    from ._core import (
        AcousticMedium,
        AcousticMode,
        AcousticProperties,
        AirProperties,
        AnnularDuctGeometry,
        AnnularMode,
        Annulus,  # type: ignore[attr-defined]
        BlochMode,
        BoundaryCondition,
        CanAnnularFlowGeometry,
        CanAnnularGeometry,
        Cd_orifice,
        Cd_rounded_entry,
        Cd_sharp_thin_plate,
        Cd_thick_plate,
        CdCorrelation,
        CombustionMethod,
        CombustionState,
        CompleteState,
        CompressibleFlowSolution,
        EquilibriumResult,
        FannoSolution,
        FannoStation,
        LinerCavity,
        LinerFlowState,
        LinerOrificeGeometry,
        NozzleSolution,
        NozzleStation,
        OrificeFlowResult,
        OrificeGeometry,
        OrificeState,
        State,
        Stream,
        ThermoState,
        ThrustResult,
        TransferMatrix,
        TransportState,
        Tube,
        UnitInfo,
        absorption_from_impedance_norm,
        acoustic_impedance,
        acoustic_properties,
        adiabatic_T_wgs,
        adiabatic_wall_temperature,
        air_properties,
        all_units,
        annular_area,
        annular_duct_eigenmodes,
        annular_duct_modes_analytical,
        annulus_axial_modes,
        annulus_azimuthal_modes,
        annulus_modes,
        axial_mode_downstream,
        axial_mode_split,
        axial_mode_upstream,
        bandwidth,
        bernoulli_P2,
        bernoulli_v2,
        beta_from_diameters,
        bilger_stoich_mixture_fraction_mass,
        bilger_Z_from_equivalence_ratio_mass,
        bulk_T_from_edge_T_and_q,
        calc_T_from_cp,
        calc_T_from_h,
        calc_T_from_s,
        can_annular_eigenmodes,
        capacity_ratio,
        closest_mode,
        combustion_equilibrium,
        combustion_state,
        combustion_state_from_streams,
        common_name,
        complete_combustion,
        complete_combustion_isothermal,
        complete_combustion_to_CO2_H2O,
        complete_combustion_to_CO2_H2O_with_fraction,
        complete_state,
        convert_to_dry_fractions,
        cooled_wall_heat_flux,
        cp,
        cp_mass,
        critical_pressure_ratio,
        cv,
        cv_mass,
        damping_ratio,
        density,
        dewpoint,
        dimple_friction_multiplier,
        dimple_nusselt_enhancement,
        dT_edge_dq,
        dT_edge_dT_bulk,
        dT_edge_dT_cold,
        dT_edge_dT_hot,
        dynamic_pressure,
        effective_viscothermal_layer,
        effectiveness_counterflow,
        effectiveness_parallelflow,
        effusion_discharge_coefficient,
        effusion_effectiveness,
        equivalence_ratio_from_bilger_Z_mass,
        equivalence_ratio_mass,
        equivalence_ratio_mole,
        expansibility_factor,
        fanno_max_length,
        fanno_pipe,
        film_cooling_effectiveness,
        film_cooling_effectiveness_avg,
        film_cooling_multirow_sellers,
        formula,
        formula_to_name,
        frequency_from_strouhal,
        frequency_from_wavelength,
        friction_colebrook,
        friction_haaland,
        friction_petukhov,
        friction_serghides,
        fuel_lhv_mass,
        fuel_lhv_molar,
        get_units,
        h,
        h_mass,
        half_wave_frequency,
        has_units,
        heat_flux,
        heat_flux_from_T_at_depth,
        heat_flux_from_T_at_edge,
        heat_rate,
        heat_rate_from_effectiveness,
        heat_transfer_area,
        heat_transfer_dT,
        helmholtz_frequency,
        helmholtz_Q,
        htc_from_nusselt,
        htc_pipe,
        humid_air_composition,
        hydraulic_diameter,
        hydraulic_diameter_annulus,
        hydraulic_diameter_rect,
        impingement_nusselt,
        input_units,
        is_whistling_risk,
        isentropic_expansion_coefficient,
        k_aluminum_6061,
        k_haynes230,
        k_inconel718,
        k_stainless_steel_316,
        k_tbc_ysz,
        kinematic_viscosity,
        liner_2dof_serial_absorption,
        liner_2dof_serial_impedance_norm,
        liner_sdof_absorption,
        liner_sdof_impedance_norm,
        list_functions_with_units,
        list_materials,
        lmtd,
        lmtd_counterflow,
        lmtd_parallelflow,
        mach_from_pressure_ratio,
        mass_flux_isentropic,
        mass_to_mole,
        min_mode_separation,
        mix,
        mixture_h,
        modes_in_range,
        molar_volume,
        mole_to_mass,
        mwmix,
        name_to_formula,
        normalize_fractions,
        nozzle_cd,
        nozzle_flow,
        nozzle_thrust,
        nozzle_thrust_cd,
        ntu,
        nusselt_dittus_boelter,
        nusselt_gnielinski,
        nusselt_petukhov,
        nusselt_sieder_tate,
        orifice_area,
        orifice_area_from_beta,
        orifice_Cd_from_K,
        orifice_Cd_from_measurement,
        orifice_dP,
        orifice_dP_Cd,
        orifice_flow,
        orifice_flow_state,
        orifice_impedance_with_flow,
        orifice_K_from_Cd,
        orifice_mdot,
        orifice_mdot_Cd,
        orifice_Q,
        orifice_Re_d_from_mdot,
        orifice_thickness_correction,
        orifice_velocity,
        orifice_velocity_from_mdot,
        output_units,
        overall_htc,
        overall_htc_wall,
        overall_htc_wall_multilayer,
        oxygen_required_per_kg_fuel,
        oxygen_required_per_kg_mixture,
        oxygen_required_per_mol_fuel,
        oxygen_required_per_mol_mixture,
        particle_velocity,
        peclet,
        pin_fin_nusselt,
        pipe_area,
        pipe_dP,
        pipe_dP_mdot,
        pipe_mdot,
        pipe_roughness,
        pipe_velocity,
        pipe_volume,
        prandtl,
        pressure_drop_pipe,
        quarter_wave_frequency,
        quarter_wave_resonator_tmm,
        reforming_equilibrium,
        reforming_equilibrium_adiabatic,
        relative_humidity_from_dewpoint,
        residence_time,
        residence_time_annulus,
        residence_time_can_annular,
        residence_time_mdot,
        residence_time_mdot_can_annular,
        residence_time_tube,
        reynolds,
        rib_enhancement_factor,
        rib_friction_multiplier,
        s,
        s_mass,
        set_equivalence_ratio_mass,
        set_equivalence_ratio_mole,
        set_fuel_stream_for_CO2,
        set_fuel_stream_for_CO2_dry,
        set_fuel_stream_for_O2,
        set_fuel_stream_for_O2_dry,
        set_fuel_stream_for_phi,
        set_fuel_stream_for_Tad,
        set_oxidizer_stream_for_CO2,
        set_oxidizer_stream_for_CO2_dry,
        set_oxidizer_stream_for_O2,
        set_oxidizer_stream_for_O2_dry,
        set_oxidizer_stream_for_Tad,
        smr_wgs_equilibrium,
        smr_wgs_equilibrium_adiabatic,  # deprecated: use reforming_equilibrium
        solve_A_eff_from_mdot,
        solve_orifice_mdot,
        solve_P0_from_mdot,
        solve_P_back_from_mdot,
        sound_pressure_level,
        space_velocity,
        specific_gas_constant,
        speed_of_sound,
        standard_dry_air_composition,
        standard_pipe_roughness,
        stokes_layer,
        strouhal,
        sweep_liner_2dof_serial_absorption,
        sweep_liner_sdof_absorption,
        thermal_conductivity,
        thermal_diffusivity,
        thermal_layer,
        thermal_resistance,
        thermal_resistance_wall,
        thermo_state,
        to_acoustic_geometry,
        transport_state,
        tube_axial_modes,
        tube_Q,
        u,
        u_mass,
        velocity_from_q,
        viscosity,
        wall_temperature_profile,
        wavelength,
        wgs_equilibrium,
        wgs_equilibrium_adiabatic,
    )
except ModuleNotFoundError:
    # Fallback: attempt to import from an installed combaero package that
    # already has _core available, then re-export the symbols.
    _core = import_module("combaero._core")
    mixture_h = _core.mixture_h
    adiabatic_T_wgs = _core.adiabatic_T_wgs
    cp = _core.cp
    cp_mass = _core.cp_mass
    h = _core.h
    h_mass = _core.h_mass
    s = _core.s
    s_mass = _core.s_mass
    cv = _core.cv
    cv_mass = _core.cv_mass
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
    expansibility_factor = _core.expansibility_factor
    solve_orifice_mdot = _core.solve_orifice_mdot
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
    fuel_lhv_molar = _core.fuel_lhv_molar
    fuel_lhv_mass = _core.fuel_lhv_mass
    complete_combustion_to_CO2_H2O = _core.complete_combustion_to_CO2_H2O
    complete_combustion_to_CO2_H2O_with_fraction = (
        _core.complete_combustion_to_CO2_H2O_with_fraction
    )
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
    heat_rate = _core.heat_rate
    heat_flux = _core.heat_flux
    heat_transfer_area = _core.heat_transfer_area
    heat_transfer_dT = _core.heat_transfer_dT
    wall_temperature_profile = _core.wall_temperature_profile
    ntu = _core.ntu
    capacity_ratio = _core.capacity_ratio
    effectiveness_counterflow = _core.effectiveness_counterflow
    effectiveness_parallelflow = _core.effectiveness_parallelflow
    heat_rate_from_effectiveness = _core.heat_rate_from_effectiveness
    heat_flux_from_T_at_edge = _core.heat_flux_from_T_at_edge
    heat_flux_from_T_at_depth = _core.heat_flux_from_T_at_depth
    bulk_T_from_edge_T_and_q = _core.bulk_T_from_edge_T_and_q
    dT_edge_dT_hot = _core.dT_edge_dT_hot
    dT_edge_dT_cold = _core.dT_edge_dT_cold
    dT_edge_dT_bulk = _core.dT_edge_dT_bulk
    dT_edge_dq = _core.dT_edge_dq
    # Acoustics
    Tube = _core.Tube
    Annulus = _core.Annulus
    CanAnnularFlowGeometry = _core.CanAnnularFlowGeometry
    BoundaryCondition = _core.BoundaryCondition
    AcousticMode = _core.AcousticMode
    to_acoustic_geometry = _core.to_acoustic_geometry
    tube_axial_modes = _core.tube_axial_modes
    annulus_axial_modes = _core.annulus_axial_modes
    annulus_azimuthal_modes = _core.annulus_azimuthal_modes
    annulus_modes = _core.annulus_modes
    modes_in_range = _core.modes_in_range
    closest_mode = _core.closest_mode
    min_mode_separation = _core.min_mode_separation
    # Mean flow correction
    axial_mode_upstream = _core.axial_mode_upstream
    axial_mode_downstream = _core.axial_mode_downstream
    axial_mode_split = _core.axial_mode_split
    # Helmholtz resonator
    helmholtz_frequency = _core.helmholtz_frequency
    # Strouhal
    strouhal = _core.strouhal
    frequency_from_strouhal = _core.frequency_from_strouhal
    # Convenience
    quarter_wave_frequency = _core.quarter_wave_frequency
    half_wave_frequency = _core.half_wave_frequency
    # Viscothermal
    stokes_layer = _core.stokes_layer
    thermal_layer = _core.thermal_layer
    effective_viscothermal_layer = _core.effective_viscothermal_layer
    # Quality factor
    helmholtz_Q = _core.helmholtz_Q
    tube_Q = _core.tube_Q
    damping_ratio = _core.damping_ratio
    bandwidth = _core.bandwidth
    # Acoustic properties
    AcousticProperties = _core.AcousticProperties
    acoustic_properties = _core.acoustic_properties
    wavelength = _core.wavelength
    frequency_from_wavelength = _core.frequency_from_wavelength
    acoustic_impedance = _core.acoustic_impedance
    sound_pressure_level = _core.sound_pressure_level
    particle_velocity = _core.particle_velocity
    # Residence time
    residence_time = _core.residence_time
    residence_time_tube = _core.residence_time_tube
    residence_time_annulus = _core.residence_time_annulus
    residence_time_can_annular = _core.residence_time_can_annular
    residence_time_mdot = _core.residence_time_mdot
    residence_time_mdot_can_annular = _core.residence_time_mdot_can_annular
    space_velocity = _core.space_velocity
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
    pressure_drop_pipe = _core.pressure_drop_pipe
    dynamic_pressure = _core.dynamic_pressure
    velocity_from_q = _core.velocity_from_q
    hydraulic_diameter = _core.hydraulic_diameter
    hydraulic_diameter_rect = _core.hydraulic_diameter_rect
    hydraulic_diameter_annulus = _core.hydraulic_diameter_annulus
    pipe_area = _core.pipe_area
    annular_area = _core.annular_area
    pipe_volume = _core.pipe_volume
    pipe_roughness = _core.pipe_roughness
    standard_pipe_roughness = _core.standard_pipe_roughness
    # Orifice Cd correlations
    OrificeGeometry = _core.OrificeGeometry
    OrificeState = _core.OrificeState
    CdCorrelation = _core.CdCorrelation
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
    # Orifice flow utilities
    OrificeFlowResult = _core.OrificeFlowResult
    orifice_flow = _core.orifice_flow
    orifice_flow_state = _core.orifice_flow_state
    orifice_velocity_from_mdot = _core.orifice_velocity_from_mdot
    orifice_area_from_beta = _core.orifice_area_from_beta
    beta_from_diameters = _core.beta_from_diameters
    orifice_Re_d_from_mdot = _core.orifice_Re_d_from_mdot
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
    "cp_mass",
    "h",
    "h_mass",
    "s",
    "s_mass",
    "u_mass",
    "cv",
    "cv_mass",
    "u",
    "u_mass",
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
    "equivalence_ratio_from_bilger_Z_mass",
    "bilger_stoich_mixture_fraction_mass",
    "bilger_Z_from_equivalence_ratio_mass",
    "expansibility_factor",
    "solve_orifice_mdot",
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
    "fuel_lhv_molar",
    "fuel_lhv_mass",
    "complete_combustion_to_CO2_H2O",
    "complete_combustion_to_CO2_H2O_with_fraction",
    # State-based API
    "State",
    "Stream",
    "AirProperties",
    "air_properties",
    "ThermoState",
    "thermo_state",
    "TransportState",
    "transport_state",
    "CompleteState",
    "complete_state",
    "EquilibriumResult",
    "CombustionMethod",
    "CombustionState",
    "combustion_state",
    "combustion_state_from_streams",
    "mix",
    # Materials database
    "k_inconel718",
    "k_haynes230",
    "k_stainless_steel_316",
    "k_aluminum_6061",
    "k_tbc_ysz",
    "list_materials",
    # Advanced cooling correlations
    "rib_enhancement_factor",
    "rib_friction_multiplier",
    "impingement_nusselt",
    "film_cooling_effectiveness",
    "film_cooling_effectiveness_avg",
    "film_cooling_multirow_sellers",
    "effusion_effectiveness",
    "effusion_discharge_coefficient",
    "pin_fin_nusselt",
    "dimple_nusselt_enhancement",
    "dimple_friction_multiplier",
    "adiabatic_wall_temperature",
    "cooled_wall_heat_flux",
    # Inverse solvers - find fuel stream
    "set_fuel_stream_for_Tad",
    "set_fuel_stream_for_phi",
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
    "is_whistling_risk",
    # Heat transfer
    "nusselt_dittus_boelter",
    "nusselt_gnielinski",
    "nusselt_sieder_tate",
    "nusselt_petukhov",
    "htc_from_nusselt",
    "htc_pipe",
    "orifice_impedance_with_flow",
    "quarter_wave_resonator_tmm",
    "TransferMatrix",
    "CanAnnularGeometry",
    "BlochMode",
    "can_annular_eigenmodes",
    "AnnularDuctGeometry",
    "AnnularMode",
    "annular_duct_eigenmodes",
    "annular_duct_modes_analytical",
    "AcousticMedium",
    "LinerOrificeGeometry",
    "LinerFlowState",
    "LinerCavity",
    "absorption_from_impedance_norm",
    "liner_sdof_impedance_norm",
    "liner_sdof_absorption",
    "sweep_liner_sdof_absorption",
    "liner_2dof_serial_impedance_norm",
    "liner_2dof_serial_absorption",
    "sweep_liner_2dof_serial_absorption",
    "overall_htc",
    "overall_htc_wall",
    "overall_htc_wall_multilayer",
    "thermal_resistance",
    "thermal_resistance_wall",
    "lmtd",
    "lmtd_counterflow",
    "lmtd_parallelflow",
    "heat_rate",
    "heat_flux",
    "heat_transfer_area",
    "heat_transfer_dT",
    "wall_temperature_profile",
    "ntu",
    "capacity_ratio",
    "effectiveness_counterflow",
    "effectiveness_parallelflow",
    "heat_rate_from_effectiveness",
    "heat_flux_from_T_at_edge",
    "heat_flux_from_T_at_depth",
    "bulk_T_from_edge_T_and_q",
    "dT_edge_dT_hot",
    "dT_edge_dT_cold",
    "dT_edge_dT_bulk",
    "dT_edge_dq",
    # Acoustics
    "Tube",
    "Annulus",
    "CanAnnularFlowGeometry",
    "BoundaryCondition",
    "AcousticMode",
    "to_acoustic_geometry",
    "tube_axial_modes",
    "annulus_axial_modes",
    "annulus_azimuthal_modes",
    "annulus_modes",
    "modes_in_range",
    "closest_mode",
    "min_mode_separation",
    # Mean flow correction
    "axial_mode_upstream",
    "axial_mode_downstream",
    "axial_mode_split",
    # Helmholtz resonator
    "helmholtz_frequency",
    # Strouhal
    "strouhal",
    "frequency_from_strouhal",
    # Convenience
    "quarter_wave_frequency",
    "half_wave_frequency",
    # Viscothermal
    "stokes_layer",
    "thermal_layer",
    "effective_viscothermal_layer",
    # Quality factor
    "helmholtz_Q",
    "tube_Q",
    "damping_ratio",
    "bandwidth",
    # Acoustic properties
    "AcousticProperties",
    "acoustic_properties",
    "wavelength",
    "frequency_from_wavelength",
    "acoustic_impedance",
    "sound_pressure_level",
    "particle_velocity",
    # Residence time
    "residence_time",
    "residence_time_tube",
    "residence_time_annulus",
    "residence_time_can_annular",
    "residence_time_mdot",
    "residence_time_mdot_can_annular",
    "space_velocity",
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
    "pressure_drop_pipe",
    "dynamic_pressure",
    "velocity_from_q",
    "hydraulic_diameter",
    "hydraulic_diameter_rect",
    "hydraulic_diameter_annulus",
    "pipe_area",
    "annular_area",
    "pipe_volume",
    "pipe_roughness",
    "standard_pipe_roughness",
    # Orifice Cd correlations
    "OrificeGeometry",
    "OrificeState",
    "CdCorrelation",
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
    # Orifice flow utilities
    "OrificeFlowResult",
    "orifice_flow",
    "orifice_flow_state",
    "orifice_velocity_from_mdot",
    "orifice_area_from_beta",
    "beta_from_diameters",
    "orifice_Re_d_from_mdot",
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
