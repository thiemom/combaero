#ifndef UNITS_DATA_H
#define UNITS_DATA_H

#include <cstddef>
#include <string_view>

namespace combaero::units {

struct Entry {
  std::string_view name;
  std::string_view input;
  std::string_view output;
};

// All units defined here - edit this table only
// Format: {function_name, input_units, output_units}
inline constexpr Entry function_units[] = {
    // -------------------------------------------------------------------------
    // composition.h - Composition Utilities
    // -------------------------------------------------------------------------
    {"mwmix", "X: mol/mol", "g/mol"},
    {"mole_to_mass", "X: mol/mol", "kg/kg"},
    {"mass_to_mole", "Y: kg/kg", "mol/mol"},
    {"normalize_fractions", "X: mol/mol", "mol/mol"},
    {"convert_to_dry_fractions", "X: mol/mol", "mol/mol"},
    {"equivalence_ratio", "mole_fractions: mol/mol", "[-]"},

    // -------------------------------------------------------------------------
    // thermo.h - Species Data
    // -------------------------------------------------------------------------
    {"species_name", "index: -", "str"},
    {"species_index_from_name", "name: str", "-"},
    {"species_molar_mass", "index: -", "g/mol"},
    {"species_molar_mass_from_name", "name: str", "g/mol"},
    {"num_species", "-", "-"},

    // -------------------------------------------------------------------------
    // thermo.h - Dimensionless NASA Polynomials
    // -------------------------------------------------------------------------
    {"cp_R", "T: K", "- (Cp/R)"},
    {"h_RT", "T: K", "- (H/(R*T))"},
    {"s_R", "T: K", "- (S/R)"},
    {"g_over_RT", "T: K", "- (G/(R*T))"},

    // -------------------------------------------------------------------------
    // thermo.h - Per-Species Properties
    // -------------------------------------------------------------------------
    {"cp_species", "T: K", "J/(mol*K)"},
    {"h_species", "T: K", "J/mol"},
    {"s_species", "T: K", "J/(mol*K)"},

    // -------------------------------------------------------------------------
    // thermo.h - Mixture Properties (Molar Basis)
    // -------------------------------------------------------------------------
    {"cp", "T: K, X: mol/mol", "J/(mol*K)"},
    {"cv", "T: K, X: mol/mol", "J/(mol*K)"},
    {"h", "T: K, X: mol/mol", "J/mol"},
    {"u", "T: K, X: mol/mol", "J/mol"},
    {"s", "T: K, X: mol/mol, P: Pa, P_ref: Pa", "J/(mol*K)"},
    {"dh_dT", "T: K, X: mol/mol", "J/(mol*K)"},
    {"ds_dT", "T: K, X: mol/mol", "J/(mol*K^2)"},
    {"dcp_dT", "T: K, X: mol/mol", "J/(mol*K^2)"},

    // -------------------------------------------------------------------------
    // thermo.h - Mixture Properties (Mass/Other Basis)
    // -------------------------------------------------------------------------
    {"density", "T: K, P: Pa, X: mol/mol", "kg/m^3"},
    {"molar_volume", "T: K, P: Pa", "m^3/mol"},
    {"specific_gas_constant", "X: mol/mol", "J/(kg*K)"},
    {"isentropic_expansion_coefficient", "T: K, X: mol/mol", "- (gamma)"},
    {"speed_of_sound", "T: K, X: mol/mol", "m/s"},
    {"cp_mass", "T: K, X: mol/mol", "J/(kg*K)"},
    {"cv_mass", "T: K, X: mol/mol", "J/(kg*K)"},
    {"h_mass", "T: K, X: mol/mol", "J/kg"},
    {"s_mass", "T: K, X: mol/mol, P: Pa, P_ref: Pa", "J/(kg*K)"},
    {"u_mass", "T: K, X: mol/mol", "J/kg"},
    {"air_properties", "T: K, P: Pa, humidity: - (0-1)",
     "AirProperties struct"},
    {"thermo_state", "T: K, P: Pa, X: mol/mol, P_ref: Pa (default 101325)",
     "ThermoState struct"},
    {"transport_state", "T: K, P: Pa, X: mol/mol", "TransportState struct"},
    {"complete_state", "T: K, P: Pa, X: mol/mol, P_ref: Pa (default 101325)",
     "CompleteState struct"},
    {"combustion_state",
     "X_fuel: mol/mol, X_ox: mol/mol, phi: -, T_reactants: K, P: Pa, "
     "fuel_name: str (default ''), smooth: bool (default False)",
     "CombustionState struct"},
    {"combustion_state_from_streams",
     "fuel_stream: Stream, ox_stream: Stream, fuel_name: str (default ''), "
     "smooth: bool (default False)",
     "CombustionState struct"},

    // -------------------------------------------------------------------------
    // materials.h - Material Thermal Conductivity
    // -------------------------------------------------------------------------
    {"k_inconel718", "T: K", "W/(m*K)"},
    {"k_haynes230", "T: K", "W/(m*K)"},
    {"k_stainless_steel_316", "T: K", "W/(m*K)"},
    {"k_aluminum_6061", "T: K", "W/(m*K)"},
    {"k_tbc_ysz", "T: K, hours: h (default 0), is_ebpvd: bool (default False)",
     "W/(m*K)"},
    {"list_materials", "-", "list[str]"},
    {"get_material_conductivity", "name: str, T: K", "W/(m*K)"},

    // -------------------------------------------------------------------------
    // cooling_correlations.h - Advanced Cooling Correlations
    // -------------------------------------------------------------------------
    {"rib_enhancement_factor", "e_D: -, pitch_to_height: -, alpha_deg: deg",
     "-"},
    {"rib_enhancement_factor_high_re",
     "e_D: -, pitch_to_height: -, alpha_deg: deg, Re: -", "-"},
    {"rib_friction_multiplier", "e_D: -, pitch_to_height: -", "-"},
    {"rib_friction_multiplier_high_re", "e_D: -, pitch_to_height: -", "-"},
    {"thermal_performance_factor", "Nu_ratio: -, f_ratio: -", "-"},
    {"impingement_nusselt",
     "Re_jet: -, Pr: -, z_D: -, x_D: - (default 0), y_D: - (default 0)", "-"},
    {"film_cooling_effectiveness", "x_D: -, M: -, DR: -, alpha_deg: deg", "-"},
    {"film_cooling_effectiveness_avg",
     "x_D: -, M: -, DR: -, alpha_deg: deg, s_D: - (default 3.0)", "-"},
    {"film_cooling_multirow_sellers",
     "row_positions_xD: list[-], eval_xD: -, M: -, DR: -, alpha_deg: deg", "-"},
    {"effusion_effectiveness",
     "x_D: -, M: -, DR: -, porosity: -, s_D: -, alpha_deg: deg", "-"},
    {"pin_fin_nusselt",
     "Re_d: -, Pr: -, L_D: -, S_D: -, X_D: -, is_staggered: bool (default "
     "True)",
     "-"},
    {"dimple_nusselt_enhancement", "Re_Dh: -, d_Dh: -, h_d: -, S_d: -", "-"},
    {"dimple_friction_multiplier", "Re_Dh: -, d_Dh: -, h_d: -", "-"},
    {"effusion_discharge_coefficient",
     "Re_d: -, P_ratio: -, alpha_deg: deg, L_D: - (default 4.0)", "-"},

    // -------------------------------------------------------------------------
    // acoustics.h - Acoustic Properties
    // -------------------------------------------------------------------------
    {"acoustic_properties",
     "f: Hz, rho: kg/m^3, c: m/s, p_rms: Pa (default 20e-6), p_ref: Pa "
     "(default 20e-6)",
     "AcousticProperties struct"},
    {"wavelength", "f: Hz, c: m/s", "m"},
    {"frequency_from_wavelength", "lambda: m, c: m/s", "Hz"},
    {"acoustic_impedance", "rho: kg/m^3, c: m/s", "Pa*s/m"},
    {"sound_pressure_level", "p_rms: Pa, p_ref: Pa (default 20e-6)", "dB"},
    {"particle_velocity", "p: Pa, rho: kg/m^3, c: m/s", "m/s"},

    // -------------------------------------------------------------------------
    // orifice.h - Orifice Flow Utilities
    // -------------------------------------------------------------------------
    {"orifice_flow",
     "geom: OrificeGeometry, dP: Pa, T: K, P: Pa, mu: Pa*s, Z: - (default 1.0)",
     "OrificeFlowResult struct"},
    {"orifice_velocity_from_mdot",
     "mdot: kg/s, rho: kg/m^3, d: m, Z: - (default 1.0)", "m/s"},
    {"orifice_area_from_beta", "D: m, beta: -", "m^2"},
    {"beta_from_diameters", "d: m, D: m", "-"},
    {"orifice_Re_d_from_mdot", "mdot: kg/s, d: m, mu: Pa*s",
     "- (Reynolds number)"},

    // -------------------------------------------------------------------------
    // thermo.h - Inverse Solvers
    // -------------------------------------------------------------------------
    {"calc_T_from_h", "h_target: J/mol, X: mol/mol", "K"},
    {"calc_T_from_s", "s_target: J/(mol*K), P: Pa, X: mol/mol", "K"},
    {"calc_T_from_cp", "cp_target: J/(mol*K), X: mol/mol", "K"},

    // -------------------------------------------------------------------------
    // transport.h - Transport Properties
    // -------------------------------------------------------------------------
    {"viscosity", "T: K, P: Pa, X: mol/mol", "Pa*s"},
    {"thermal_conductivity", "T: K, P: Pa, X: mol/mol", "W/(m*K)"},
    {"prandtl", "T: K, P: Pa, X: mol/mol", "- (Pr)"},
    {"kinematic_viscosity", "T: K, P: Pa, X: mol/mol", "m^2/s"},
    {"thermal_diffusivity", "T: K, P: Pa, X: mol/mol", "m^2/s"},
    {"reynolds", "T: K, P: Pa, X: mol/mol, V: m/s, L: m", "- (Re)"},
    {"peclet", "T: K, P: Pa, X: mol/mol, V: m/s, L: m", "- (Pe)"},

    // -------------------------------------------------------------------------
    // compressible.h - Nozzle Flow
    // -------------------------------------------------------------------------
    {"nozzle_flow", "T0: K, P0: Pa, P_back: Pa, A_eff: m^2, X: mol/mol",
     "CompressibleFlowSolution"},
    {"nozzle_quasi1d", "T0: K, P0: Pa, P_exit: Pa, x: m, A: m^2, X: mol/mol",
     "NozzleSolution"},
    {"nozzle_cd", "T0: K, P0: Pa, P_exit: Pa, A: m^2, x: m, X: mol/mol",
     "NozzleSolution"},
    {"critical_pressure_ratio", "T0: K, P0: Pa, X: mol/mol", "- (P*/P0)"},
    {"mach_from_pressure_ratio", "T0: K, P0: Pa, P: Pa, X: mol/mol", "- (M)"},
    {"mass_flux_isentropic", "T0: K, P0: Pa, P: Pa, X: mol/mol", "kg/(m^2*s)"},

    // -------------------------------------------------------------------------
    // compressible.h - Fanno Flow
    // -------------------------------------------------------------------------
    {"fanno_channel", "T_in: K, P_in: Pa, u_in: m/s, L: m, D: m, f: -, X: mol/mol",
     "FannoSolution"},
    {"fanno_channel_rough",
     "T_in: K, P_in: Pa, u_in: m/s, L: m, D: m, roughness: m, X: mol/mol, "
     "correlation: str, f_multiplier: - (default 1.0)",
     "FannoSolution"},
    {"fanno_max_length", "T_in: K, P_in: Pa, u_in: m/s, D: m, f: -, X: mol/mol",
     "m"},

    // -------------------------------------------------------------------------
    // compressible.h - Thrust
    // -------------------------------------------------------------------------
    {"nozzle_thrust", "NozzleSolution, P_amb: Pa", "ThrustResult"},

    // -------------------------------------------------------------------------
    // incompressible.h - Bernoulli & Orifice
    // -------------------------------------------------------------------------
    {"bernoulli_P2", "P1: Pa, v1: m/s, v2: m/s, rho: kg/m^3, dz: m, g: m/s^2",
     "Pa"},
    {"bernoulli_v2", "P1: Pa, P2: Pa, v1: m/s, rho: kg/m^3, dz: m, g: m/s^2",
     "m/s"},
    {"orifice_mdot", "P1: Pa, P2: Pa, A: m^2, Cd: -, rho: kg/m^3", "kg/s"},
    {"orifice_Q", "P1: Pa, P2: Pa, A: m^2, Cd: -, rho: kg/m^3", "m^3/s"},
    {"orifice_velocity", "P1: Pa, P2: Pa, rho: kg/m^3", "m/s"},
    {"orifice_area", "mdot: kg/s, P1: Pa, P2: Pa, Cd: -, rho: kg/m^3", "m^2"},
    {"orifice_dP", "mdot: kg/s, A: m^2, Cd: -, rho: kg/m^3", "Pa"},
    {"expansibility_factor", "beta: -, dP: Pa, P_upstream: Pa, kappa: -", "-"},

    // -------------------------------------------------------------------------
    // incompressible.h - Thermo-Aware High-Level Functions
    // -------------------------------------------------------------------------
    {"channel_flow", "T: K, P: Pa, X: mol/mol, u: m/s, L: m, D: m, f: -",
     "IncompressibleFlowSolution"},
    {"channel_flow_rough",
     "T: K, P: Pa, X: mol/mol, u: m/s, L: m, D: m, roughness: m, correlation: "
     "str",
     "IncompressibleFlowSolution"},
    {"orifice_flow_thermo",
     "T: K, P: Pa, X: mol/mol, P_back: Pa, A: m^2, Cd: -",
     "IncompressibleFlowSolution"},
    {"channel_pressure_drop",
     "T: K, P: Pa, X: mol/mol, v: m/s, D: m, L: m, roughness: m, correlation: "
     "str",
     "tuple(Pa, -, -)"},

    // -------------------------------------------------------------------------
    // solver_interface.h - Compressible Flow Elements
    // -------------------------------------------------------------------------
    {"orifice_compressible_mdot_and_jacobian",
     "T0: K, P0: Pa, P_back: Pa, X: mol/mol, Cd: -, area: m^2, beta: -",
     "tuple(kg/s, kg/(s*Pa), kg/(s*Pa), kg/(s*K))"},
    {"orifice_compressible_residuals_and_jacobian",
     "m_dot: kg/s, P_total_up: Pa, T_up: K, Y_up: kg/kg, P_static_down: Pa, "
     "Cd: -, area: m^2, beta: -",
     "OrificeResult"},
    {"channel_compressible_mdot_and_jacobian",
     "T_in: K, P_in: Pa, u_in: m/s, X: mol/mol, L: m, D: m, roughness: m, "
     "friction_model: str, f_multiplier: - (default 1.0)",
     "tuple(Pa, -, Pa/K, Pa/(m/s))"},
    {"channel_compressible_residuals_and_jacobian",
     "m_dot: kg/s, P_total_up: Pa, T_up: K, Y_up: kg/kg, P_static_down: Pa, "
     "L: m, D: m, roughness: m, friction_model: str, f_multiplier: - (default 1.0)",
     "ChannelResult"},

    // -------------------------------------------------------------------------
    // geometry.h - Geometric Utilities
    // -------------------------------------------------------------------------
    {"circular_area", "D: m", "m^2"},
    {"annular_area", "D_outer: m, D_inner: m", "m^2"},
    {"cylinder_volume", "D: m, L: m", "m^3"},
    {"channel_roughness", "material: str", "m"},
    {"standard_channel_roughness", "-", "dict[str, m]"},
    {"hydraulic_diameter", "A: m^2, P_wetted: m", "m"},
    {"hydraulic_diameter_rect", "a: m, b: m", "m"},
    {"hydraulic_diameter_annulus", "D_outer: m, D_inner: m", "m"},

    // -------------------------------------------------------------------------
    // friction.h - Friction Factor Correlations
    // -------------------------------------------------------------------------
    {"friction_haaland", "Re: -, e_D: -", "- (f)"},
    {"friction_serghides", "Re: -, e_D: -", "- (f)"},
    {"friction_colebrook", "Re: -, e_D: -", "- (f)"},
    {"friction_petukhov", "Re: -", "- (f)"},

    // -------------------------------------------------------------------------
    // heat_transfer.h - Heat Transfer Correlations
    // -------------------------------------------------------------------------
    {"nusselt_dittus_boelter", "Re: -, Pr: -, heating: bool", "- (Nu)"},
    {"nusselt_gnielinski", "Re: -, Pr: -, f: - (optional)", "- (Nu)"},
    {"nusselt_sieder_tate", "Re: -, Pr: -, mu_ratio: -", "- (Nu)"},
    {"nusselt_petukhov", "Re: -, Pr: -, f: - (optional)", "- (Nu)"},
    {"htc_from_nusselt", "Nu: -, k: W/(m*K), L: m", "W/(m^2*K)"},
    {"nusselt_circular_channel", "State, V: m/s, D: m", "- (Nu)"},
    {"htc_circular_channel", "State, V: m/s, D: m", "W/(m^2*K)"},
    {"htc_circular_channel",
     "T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, correlation: str, "
     "heating: bool, mu_ratio: -, roughness: m",
     "tuple(W/(m^2*K), -, -)"},
    {"lmtd", "dT1: K, dT2: K", "K"},
    {"lmtd_counterflow",
     "T_hot_in: K, T_hot_out: K, T_cold_in: K, T_cold_out: K", "K"},
    {"lmtd_parallelflow",
     "T_hot_in: K, T_hot_out: K, T_cold_in: K, T_cold_out: K", "K"},
    {"overall_htc", "h_values: W/(m^2*K), t_over_k: m^2*K/W", "W/(m^2*K)"},
    {"overall_htc_wall",
     "h_inner, h_outer: W/(m^2*K), t_over_k_layers: m^2*K/W", "W/(m^2*K)"},
    {"overall_htc_wall",
     "h_inner, h_outer: W/(m^2*K), t_wall: m, k_wall: W/(m*K)", "W/(m^2*K)"},
    {"thermal_resistance", "h: W/(m^2*K), A: m^2", "K/W"},
    {"thermal_resistance_wall", "t: m, k: W/(m*K), A: m^2", "K/W"},
    {"heat_rate", "U: W/(m^2*K), A: m^2, dT: K", "W"},
    {"heat_flux", "U: W/(m^2*K), dT: K", "W/m^2"},
    {"heat_transfer_area", "Q: W, U: W/(m^2*K), dT: K", "m^2"},
    {"heat_transfer_dT", "Q: W, U: W/(m^2*K), A: m^2", "K"},
    {"wall_temperature_profile",
     "T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W",
     "K (vector)"},
    {"ntu", "U: W/(m^2*K), A: m^2, C_min: W/K", "-"},
    {"capacity_ratio", "C_min: W/K, C_max: W/K", "-"},
    {"effectiveness_counterflow", "NTU: -, C_r: -", "-"},
    {"effectiveness_parallelflow", "NTU: -, C_r: -", "-"},
    {"heat_rate_from_effectiveness",
     "epsilon: -, C_min: W/K, T_hot_in: K, T_cold_in: K", "W"},
    {"adiabatic_wall_temperature", "T_hot: K, T_coolant: K, eta: -", "K"},
    {"cooled_wall_heat_flux",
     "T_hot: K, T_coolant: K, h_hot: W/(m^2*K), h_coolant: W/(m^2*K), eta: -, "
     "t_wall: m, k_wall: W/(m*K)",
     "W/m^2"},
    {"heat_flux_from_T_at_edge",
     "T_measured: K, edge_idx, T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), "
     "t_over_k: m^2*K/W",
     "W/m^2"},
    {"heat_flux_from_T_at_depth",
     "T_measured: K, depth: m, T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), "
     "thicknesses: m, conductivities: W/(m*K)",
     "W/m^2"},
    {"bulk_T_from_edge_T_and_q",
     "T_measured: K, edge_idx, q: W/m^2, h_hot, h_cold: W/(m^2*K), t_over_k: "
     "m^2*K/W, solve_for: str",
     "K"},
    {"dT_edge_dT_hot", "edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W",
     "- (dT_edge/dT_hot)"},
    {"dT_edge_dT_cold", "edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W",
     "- (dT_edge/dT_cold)"},
    {"dT_edge_dT_bulk", "edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W",
     "- (dT/dT_hot, dT/dT_cold)"},
    {"dT_edge_dq", "edge_idx, h_hot: W/(m^2*K), t_over_k: m^2*K/W",
     "K*m^2/W (dT_edge/dq)"},

    // -------------------------------------------------------------------------
    // acoustics.h - Acoustic Mode Analysis
    // -------------------------------------------------------------------------
    {"tube_axial_modes",
     "Tube (L: m, D: m), c: m/s, upstream: BC, downstream: BC, n_max",
     "Hz (vector)"},
    {"annulus_axial_modes",
     "Annulus (L: m, D_inner: m, D_outer: m), c: m/s, upstream: BC, "
     "downstream: BC, n_max",
     "Hz (vector)"},
    {"annulus_azimuthal_modes", "Annulus, c: m/s, m_max", "Hz (vector)"},
    {"annulus_modes",
     "Annulus, c: m/s, upstream: BC, downstream: BC, n_max, m_max",
     "Hz (vector)"},
    {"modes_in_range", "modes, f_min: Hz, f_max: Hz", "Hz (vector)"},
    {"closest_mode", "modes, f_target: Hz", "AcousticMode"},
    {"min_mode_separation", "modes", "Hz"},
    {"axial_mode_upstream", "f0: Hz, M: -", "Hz"},
    {"axial_mode_downstream", "f0: Hz, M: -", "Hz"},
    {"axial_mode_split", "f0: Hz, M: -", "(Hz, Hz)"},
    {"helmholtz_frequency",
     "V: m^3, A_neck: m^2, L_neck: m, c: m/s, end_correction: -", "Hz"},
    {"strouhal", "f: Hz, L: m, u: m/s", "- (St)"},
    {"frequency_from_strouhal", "St: -, L: m, u: m/s", "Hz"},
    {"quarter_wave_frequency", "L: m, c: m/s", "Hz"},
    {"orifice_impedance_with_flow",
     "freq: Hz, u_bias: m/s, u_grazing: m/s, d_orifice: m, l_orifice: m, "
     "porosity: -, Cd: -, rho: kg/m^3, c: m/s",
     "complex [-]"},
    {"quarter_wave_resonator_tmm",
     "freq: Hz, L_tube: m, A_duct: m^2, A_tube: m^2, d_orifice: m, l_orifice: "
     "m, porosity: -, Cd: -, u_bias: m/s, u_grazing: m/s, rho: kg/m^3, c: m/s",
     "TransferMatrix"},
    {"absorption_from_impedance_norm", "z_norm: complex [-]", "- (alpha)"},
    {"liner_sdof_impedance_norm",
     "freq: Hz, orifice: LinerOrificeGeometry, cavity: LinerCavity, flow: "
     "LinerFlowState, medium: AcousticMedium",
     "complex [-]"},
    {"liner_sdof_absorption",
     "freq: Hz, orifice: LinerOrificeGeometry, cavity: LinerCavity, flow: "
     "LinerFlowState, medium: AcousticMedium",
     "- (alpha)"},
    {"sweep_liner_sdof_absorption",
     "freqs: Hz (vector), orifice: LinerOrificeGeometry, cavity: LinerCavity, "
     "flow: LinerFlowState, medium: AcousticMedium",
     "- (alpha vector)"},
    {"liner_2dof_serial_impedance_norm",
     "freq: Hz, face_orifice: LinerOrificeGeometry, septum_orifice: "
     "LinerOrificeGeometry, depth_1: m, depth_2: m, face_flow: LinerFlowState, "
     "septum_flow: LinerFlowState, medium: AcousticMedium",
     "complex [-]"},
    {"liner_2dof_serial_absorption",
     "freq: Hz, face_orifice: LinerOrificeGeometry, septum_orifice: "
     "LinerOrificeGeometry, depth_1: m, depth_2: m, face_flow: LinerFlowState, "
     "septum_flow: LinerFlowState, medium: AcousticMedium",
     "- (alpha)"},
    {"sweep_liner_2dof_serial_absorption",
     "freqs: Hz (vector), face_orifice: LinerOrificeGeometry, septum_orifice: "
     "LinerOrificeGeometry, depth_1: m, depth_2: m, face_flow: LinerFlowState, "
     "septum_flow: LinerFlowState, medium: AcousticMedium",
     "- (alpha vector)"},
    {"is_whistling_risk", "freq: Hz, u_bias: m/s, d_orifice: m", "bool"},
    {"annular_duct_modes_analytical",
     "geom: AnnularDuctGeometry, c: m/s, f_max: Hz, bc_ends: BoundaryCondition",
     "AnnularMode vector"},
    {"to_acoustic_geometry",
     "flow_geom: CanAnnularFlowGeometry, n_cans: -, L_can: m, D_can: m",
     "CanAnnularGeometry struct"},
    {"half_wave_frequency", "L: m, c: m/s", "Hz"},
    {"stokes_layer", "nu: m^2/s, f: Hz", "m"},
    {"thermal_layer", "alpha: m^2/s, f: Hz", "m"},
    {"effective_viscothermal_layer", "delta_nu: m, delta_kappa: m, gamma: -",
     "m"},
    {"helmholtz_Q",
     "V: m^3, A_neck: m^2, L_neck: m, nu: m^2/s, alpha: m^2/s, gamma: -, f: Hz",
     "- (Q)"},
    {"tube_Q", "L: m, D: m, nu: m^2/s, alpha: m^2/s, gamma: -, f: Hz", "- (Q)"},
    {"damping_ratio", "Q: -", "- (zeta)"},
    {"bandwidth", "f0: Hz, Q: -", "Hz"},

    // -------------------------------------------------------------------------
    // geometry.h - Residence Time
    // -------------------------------------------------------------------------
    {"residence_time", "V: m^3, Q: m^3/s", "s"},
    {"residence_time_tube", "tube: Tube, Q: m^3/s", "s"},
    {"residence_time_annulus", "annulus: Annulus, Q: m^3/s", "s"},
    {"residence_time_can_annular", "geom: CanAnnularFlowGeometry, Q: m^3/s",
     "s"},
    {"residence_time_mdot", "V: m^3, mdot: kg/s, rho: kg/m^3", "s"},
    {"residence_time_mdot_can_annular",
     "geom: CanAnnularFlowGeometry, mdot: kg/s, rho: kg/m^3", "s"},
    {"space_velocity", "Q: m^3/s, V: m^3", "1/s"},

    // -------------------------------------------------------------------------
    // orifice.h - Discharge Coefficients
    // -------------------------------------------------------------------------
    {"Cd_sharp_thin_plate", "OrificeGeometry, OrificeState", "- (Cd)"},
    {"Cd_thick_plate", "OrificeGeometry, OrificeState", "- (Cd)"},
    {"Cd_rounded_entry", "OrificeGeometry, OrificeState", "- (Cd)"},
    {"Cd", "OrificeGeometry, OrificeState", "- (Cd)"},

    // -------------------------------------------------------------------------
    // humidair.h - Humid Air Properties
    // -------------------------------------------------------------------------
    {"saturation_vapor_pressure", "T: K", "Pa"},
    {"vapor_pressure", "T: K, RH: - (0-1)", "Pa"},
    {"humidity_ratio", "T: K, P: Pa, RH: - (0-1)", "kg/kg"},
    {"water_vapor_mole_fraction", "T: K, P: Pa, RH: - (0-1)", "mol/mol"},
    {"humid_air_composition", "T: K, P: Pa, RH: - (0-1)", "mol/mol"},
    {"dewpoint", "T: K, P: Pa, RH: - (0-1)", "K"},
    {"relative_humidity_from_dewpoint", "T: K, Tdp: K, P: Pa", "- (0-1)"},
    {"wet_bulb_temperature", "T: K, P: Pa, RH: - (0-1)", "K"},
    {"humid_air_enthalpy", "T: K, P: Pa, RH: - (0-1)", "J/kg"},
    {"humid_air_density", "T: K, P: Pa, RH: - (0-1)", "kg/m^3"},
    {"HumidAir::set_TP_RH", "T: K, P: Pa, RH: - (0-1)", "HumidAir&"},
    {"HumidAir::set_TP_dewpoint", "T: K, P: Pa, Tdp: K", "HumidAir&"},
    {"HumidAir::set_TP_omega", "T: K, P: Pa, w: kg/kg", "HumidAir&"},
    {"HumidAir::RH", "-", "- (0-1)"},
    {"HumidAir::dewpoint", "-", "K"},
    {"HumidAir::wet_bulb", "-", "K"},
    {"HumidAir::humidity_ratio", "-", "kg/kg"},
    {"HumidAir::h_mass", "-", "J/kg"},
    {"HumidAir::TP_RH", "-", "(K, Pa, -)"},
    {"HumidAir::TP_dewpoint", "-", "(K, Pa, K)"},

    // -------------------------------------------------------------------------
    // combustion.h - Stoichiometry
    // -------------------------------------------------------------------------
    {"oxygen_required_per_mol_fuel", "-", "mol O2/mol fuel"},
    {"oxygen_required_per_kg_fuel", "-", "mol O2/kg fuel"},
    {"oxygen_required_per_mol_mixture", "X: mol/mol", "mol O2/mol mix"},
    {"oxygen_required_per_kg_mixture", "X: mol/mol", "mol O2/kg mix"},
    {"fuel_lhv_molar", "X_fuel: mol/mol, reference_temperature: K",
     "J/mol fuel"},
    {"fuel_lhv_mass", "X_fuel: mol/mol, reference_temperature: K", "J/kg fuel"},
    {"dryair_required_per_mol_fuel", "-", "mol air/mol fuel"},
    {"dryair_required_per_kg_fuel", "-", "mol air/kg fuel"},
    {"dryair_required_per_mol_mixture", "X: mol/mol", "mol air/mol mix"},
    {"dryair_required_per_kg_mixture", "X: mol/mol", "mol air/kg mix"},

    // -------------------------------------------------------------------------
    // combustion.h - Equivalence Ratio
    // -------------------------------------------------------------------------
    {"equivalence_ratio_mole", "X: mol/mol", "- (phi)"},
    {"set_equivalence_ratio_mole", "phi: -, X: mol/mol", "mol/mol"},
    {"equivalence_ratio_mass", "Y: kg/kg", "- (phi)"},
    {"set_equivalence_ratio_mass", "phi: -, Y: kg/kg", "kg/kg"},

    // -------------------------------------------------------------------------
    // combustion.h - Mixture Fraction (Bilger)
    // -------------------------------------------------------------------------
    {"bilger_beta", "Y: kg/kg", "-"},
    {"bilger_mixture_fraction", "Y: kg/kg", "- (Z)"},
    {"bilger_stoich_mixture_fraction_mass", "Y: kg/kg", "- (Z_st)"},
    {"equivalence_ratio_from_bilger_Z_mass", "Z: -, Y: kg/kg", "- (phi)"},
    {"bilger_Z_from_equivalence_ratio_mass", "phi: -, Y: kg/kg", "- (Z)"},

    // -------------------------------------------------------------------------
    // combustion.h - Complete Combustion
    // -------------------------------------------------------------------------
    {"complete_combustion_to_CO2_H2O", "X: mol/mol", "mol/mol"},
    {"complete_combustion",
     "State (T: K, P: Pa, X: mol/mol), smooth: bool (default False)", "State"},
    {"complete_combustion_isothermal", "State, smooth: bool (default False)",
     "State"},

    // -------------------------------------------------------------------------
    // combustion.h - Stream Solvers
    // -------------------------------------------------------------------------
    {"set_fuel_stream_for_phi", "phi: -", "Stream"},
    {"set_fuel_stream_for_Tad", "T_ad_target: K", "Stream"},
    {"set_fuel_stream_for_O2", "X_O2_target: mol/mol", "Stream"},
    {"set_fuel_stream_for_O2_dry", "X_O2_dry_target: mol/mol", "Stream"},
    {"set_fuel_stream_for_CO2", "X_CO2_target: mol/mol", "Stream"},
    {"set_fuel_stream_for_CO2_dry", "X_CO2_dry_target: mol/mol", "Stream"},
    {"set_oxidizer_stream_for_Tad", "T_ad_target: K", "Stream"},
    {"set_oxidizer_stream_for_O2", "X_O2_target: mol/mol", "Stream"},
    {"set_oxidizer_stream_for_O2_dry", "X_O2_dry_target: mol/mol", "Stream"},
    {"set_oxidizer_stream_for_CO2", "X_CO2_target: mol/mol", "Stream"},
    {"set_oxidizer_stream_for_CO2_dry", "X_CO2_dry_target: mol/mol", "Stream"},

    // -------------------------------------------------------------------------
    // equilibrium.h - Chemical Equilibrium
    // -------------------------------------------------------------------------
    {"wgs_equilibrium", "State (T: K, P: Pa, X: mol/mol)", "State"},
    {"wgs_equilibrium_adiabatic", "State", "State"},
    {"smr_wgs_equilibrium", "State", "State"},
    {"smr_wgs_equilibrium_adiabatic", "State", "State"},
    {"reforming_equilibrium", "State", "State"},
    {"reforming_equilibrium_adiabatic", "State", "State"},
    {"combustion_equilibrium", "State, smooth: bool (default False)", "State"},

    // -------------------------------------------------------------------------
    // state.h - State Properties
    // -------------------------------------------------------------------------
    {"State::T", "-", "K"},
    {"State::P", "-", "Pa"},
    {"State::X", "-", "mol/mol"},
    {"State::Y", "-", "kg/kg"},
    {"State::set_X", "X: mol/mol", "-"},
    {"State::set_Y", "Y: kg/kg", "-"},
    {"State::mw", "-", "g/mol"},
    {"State::cp", "-", "J/(kg*K)"},
    {"State::cv", "-", "J/(kg*K)"},
    {"State::h", "-", "J/kg"},
    {"State::u", "-", "J/kg"},
    {"State::s", "-", "J/(kg*K)"},
    {"State::rho", "-", "kg/m^3"},
    {"State::R", "-", "J/(kg*K)"},
    {"State::gamma", "-", "- (gamma)"},
    {"State::a", "-", "m/s"},
    {"State::mu", "-", "Pa*s"},
    {"State::k", "-", "W/(m*K)"},
    {"State::nu", "-", "m^2/s"},
    {"State::Pr", "-", "- (Pr)"},
    {"State::alpha", "-", "m^2/s"},

    // -------------------------------------------------------------------------
    // state.h - Joint Getters / Setters (C++ & Python tuples)
    // -------------------------------------------------------------------------
    {"State::TPX", "-", "tuple(K, Pa, mol/mol)"},
    {"State::TPY", "-", "tuple(K, Pa, kg/kg)"},
    {"State::TP", "-", "tuple(K, Pa)"},
    {"State::set_TP", "T: K, P: Pa", "-"},
    {"State::set_TPX", "T: K, P: Pa, X: mol/mol", "-"},
    {"State::set_TPY", "T: K, P: Pa, Y: kg/kg", "-"},
    {"State::DP_mass", "-", "tuple(kg/m^3, Pa)"},
    {"State::set_DP_mass", "rho: kg/m^3, P: Pa", "-"},
    {"State::HP_mass", "-", "tuple(J/kg, Pa)"},
    {"State::set_HP_mass", "h: J/kg, P: Pa", "-"},
    {"State::SP_mass", "-", "tuple(J/(kg*K), Pa)"},
    {"State::set_SP_mass", "s: J/(kg*K), P: Pa", "-"},
    {"State::UP_mass", "-", "tuple(J/kg, Pa)"},
    {"State::set_UP_mass", "u: J/kg, Pa", "-"},
    {"State::UV_mass", "-", "tuple(J/kg, m^3/kg)"},
    {"State::set_UV_mass", "u: J/kg, v: m^3/kg", "-"},
    {"State::SV_mass", "-", "tuple(J/(kg*K), m^3/kg)"},
    {"State::set_SV_mass", "s: J/(kg*K), v: m^3/kg", "-"},
    {"State::PV_mass", "-", "tuple(Pa, m^3/kg)"},
    {"State::set_PV_mass", "P: Pa, v: m^3/kg", "-"},
    {"State::VH_mass", "-", "tuple(m^3/kg, J/kg)"},
    {"State::set_VH_mass", "v: m^3/kg, h: J/kg", "-"},
    {"State::SH_mass", "-", "tuple(J/(kg*K), J/kg)"},
    {"State::set_SH_mass", "s: J/(kg*K), h: J/kg", "-"},

    // -------------------------------------------------------------------------
    // state.h - Stream Properties
    // -------------------------------------------------------------------------
    {"Stream::mdot", "-", "kg/s"},
    {"Stream::T", "-", "K"},
    {"Stream::P", "-", "Pa"},
    {"Stream::X", "-", "mol/mol"},
    {"Stream::set_mdot", "mdot: kg/s", "-"},
    {"Stream::set_T", "T: K", "-"},
    {"Stream::set_P", "P: Pa", "-"},
    {"Stream::set_X", "X: mol/mol", "-"},

    // -------------------------------------------------------------------------
    // stagnation.h - Stagnation/Static Conversions
    // -------------------------------------------------------------------------
    {"kinetic_energy", "v: m/s", "J/kg"},
    {"h0_from_static", "h_static: J/kg, v: m/s", "J/kg"},
    {"v_from_h0", "h0: J/kg, h_static: J/kg", "m/s"},
    {"mach_number", "v: m/s, T: K, X: mol/mol", "- (M)"},
    {"T0_from_static", "T: K, M: -, X: mol/mol", "K"},
    {"T0_from_static_v", "T: K, v: m/s, X: mol/mol", "K"},
    {"T_from_stagnation", "T0: K, M: -, X: mol/mol", "K"},
    {"P0_from_static", "P: Pa, T: K, M: -, X: mol/mol", "Pa"},
    {"P_from_stagnation", "P0: Pa, T0: K, M: -, X: mol/mol", "Pa"},
    {"recovery_factor", "Pr: -", "- (r)"},
    {"T_adiabatic_wall", "T_static: K, v: m/s, T: K, P: Pa, X: mol/mol", "K"},
    {"T_adiabatic_wall_mach", "T_static: K, M: -, T: K, P: Pa, X: mol/mol",
     "K"},

    // -------------------------------------------------------------------------
    // heat_transfer.h - ChannelResult fields
    // -------------------------------------------------------------------------
    {"ChannelResult::h", "-", "W/(m^2*K)"},
    {"ChannelResult::Nu", "-", "- (Nu)"},
    {"ChannelResult::Re", "-", "- (Re)"},
    {"ChannelResult::Pr", "-", "- (Pr)"},
    {"ChannelResult::f", "-", "- (friction/loss coeff)"},
    {"ChannelResult::dP", "-", "Pa"},
    {"ChannelResult::M", "-", "- (M)"},
    {"ChannelResult::T_aw", "-", "K"},
    {"ChannelResult::q", "-", "W/m^2"},
    {"ChannelResult::dh_dmdot", "-", "W/(m^2*K*kg/s)"},
    {"ChannelResult::dh_dT", "-", "W/(m^2*K^2)"},
    {"ChannelResult::ddP_dmdot", "-", "Pa*s/kg"},
    {"ChannelResult::ddP_dT", "-", "Pa/K"},
    {"ChannelResult::dT_aw_dmdot", "-", "K*s/kg"},
    {"ChannelResult::dT_aw_dT", "-", "- (dimensionless)"},
    {"ChannelResult::dq_dmdot", "-", "W*s/(m^2*kg)"},
    {"ChannelResult::dq_dT", "-", "W/(m^2*K)"},
    {"ChannelResult::dq_dT_hot", "-", "W/(m^2*K)"},

    // -------------------------------------------------------------------------
    // heat_transfer.h - Channel flow functions (HTC + pressure drop)
    // -------------------------------------------------------------------------
    {"channel_smooth",
     "T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, length: m, T_wall: "
     "K",
     "ChannelResult"},
    {"channel_ribbed",
     "T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, length: m, e_D: -, "
     "pitch_to_height: -, alpha_deg: deg, T_hot: K",
     "ChannelResult"},
    {"channel_dimpled",
     "T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, length: m, d_Dh: -, "
     "h_d: -, S_d: -, T_hot: K",
     "ChannelResult"},
    {"channel_pin_fin",
     "T: K, P: Pa, X: mol/mol, velocity: m/s, channel_height: m, pin_diameter: "
     "m, S_D: -, X_D: -, N_rows: -, T_hot: K",
     "ChannelResult"},
    {"channel_impingement",
     "T: K, P: Pa, X: mol/mol, mdot_jet: kg/s, d_jet: m, z_D: -, x_D: -, y_D: "
     "-, A_target: m^2, T_hot: K",
     "ChannelResult"},

    // -------------------------------------------------------------------------
    // cooling_correlations.h - Pin fin friction (new scalar)
    // -------------------------------------------------------------------------
    {"pin_fin_friction", "Re_d: -", "- (f_pin)"},

    // -------------------------------------------------------------------------
    // correlation_status.h - Extrapolation validity utilities
    // -------------------------------------------------------------------------
    // --- Auto-generated properties ---
    {"mixture_h", "-", "-"},
    {"adiabatic_T_wgs", "-", "-"},
    {"normalize_fractions", "-", "-"},
    {"convert_to_dry_fractions", "-", "-"},
    {"solve_orifice_mdot", "-", "Pa"},
    {"dry_air", "-", "X: mol/mol"},
    {"humid_air_enthalpy_nasa9", "-", "J/kg dry air"},
    {"complete_combustion_to_CO2_H2O_with_fraction", "-", "-"},
    {"HumidAir::rh", "-", "0-1"},
    {"HumidAir::state", "-", "-"},
    {"AirProperties::P", "-", "Pa"},
    {"AirProperties::Pr", "-", "-"},
    {"AirProperties::T", "-", "K"},
    {"AirProperties::a", "-", "m/s"},
    {"AirProperties::alpha", "-", "m²/s"},
    {"AirProperties::cp", "-", "J/(kg·K)"},
    {"AirProperties::cv", "-", "J/(kg·K)"},
    {"AirProperties::gamma", "-", "-"},
    {"AirProperties::k", "-", "W/(m·K)"},
    {"AirProperties::mu", "-", "Pa·s"},
    {"AirProperties::nu", "-", "m²/s"},
    {"AirProperties::rho", "-", "kg/m³"},
    {"ThermoState::P", "-", "Pa"},
    {"ThermoState::T", "-", "K"},
    {"ThermoState::a", "-", "m/s"},
    {"ThermoState::cp", "-", "J/(kg·K)"},
    {"ThermoState::cp_mole", "-", "J/(mol·K)"},
    {"ThermoState::cv", "-", "J/(kg·K)"},
    {"ThermoState::cv_mole", "-", "J/(mol·K)"},
    {"ThermoState::gamma", "-", "-"},
    {"ThermoState::h", "-", "J/kg"},
    {"ThermoState::h_mole", "-", "J/mol"},
    {"ThermoState::mw", "-", "g/mol"},
    {"ThermoState::rho", "-", "kg/m³"},
    {"ThermoState::s", "-", "J/(kg·K)"},
    {"ThermoState::s_mole", "-", "J/(mol·K)"},
    {"ThermoState::u", "-", "J/kg"},
    {"ThermoState::u_mole", "-", "J/mol"},
    {"TransportState::P", "-", "Pa"},
    {"TransportState::Pr", "-", "-"},
    {"TransportState::T", "-", "K"},
    {"TransportState::a", "-", "m/s"},
    {"TransportState::alpha", "-", "m²/s"},
    {"TransportState::cp", "-", "J/(kg·K)"},
    {"TransportState::cv", "-", "J/(kg·K)"},
    {"TransportState::gamma", "-", "-"},
    {"TransportState::k", "-", "W/(m·K)"},
    {"TransportState::mu", "-", "Pa·s"},
    {"TransportState::nu", "-", "m²/s"},
    {"TransportState::rho", "-", "kg/m³"},
    {"CompleteState::thermo", "-", "-"},
    {"CompleteState::transport", "-", "-"},
    {"EquilibriumResult::T_in", "-", "K"},
    {"EquilibriumResult::converged", "-", "-"},
    {"EquilibriumResult::delta_T", "-", "K"},
    {"EquilibriumResult::state", "-", "-"},
    {"CombustionState::fuel_burn_fraction", "-", "0-1"},
    {"CombustionState::fuel_name", "-", "-"},
    {"CombustionState::method", "-", "-"},
    {"CombustionState::mixture_fraction", "-", "-"},
    {"CombustionState::phi", "-", "-"},
    {"CombustionState::products", "-", "-"},
    {"CombustionState::reactants", "-", "-"},
    {"CompressibleFlowSolution::M", "-", "-"},
    {"CompressibleFlowSolution::choked", "-", "-"},
    {"CompressibleFlowSolution::mdot", "-", "kg/s"},
    {"CompressibleFlowSolution::outlet", "-", "-"},
    {"CompressibleFlowSolution::stagnation", "-", "-"},
    {"CompressibleFlowSolution::v", "-", "m/s"},
    {"NozzleStation::A", "-", "m²"},
    {"NozzleStation::M", "-", "-"},
    {"NozzleStation::P", "-", "Pa"},
    {"NozzleStation::T", "-", "K"},
    {"NozzleStation::h", "-", "J/kg"},
    {"NozzleStation::rho", "-", "kg/m³"},
    {"NozzleStation::u", "-", "m/s"},
    {"NozzleStation::x", "-", "m"},
    {"NozzleSolution::A_throat", "-", "m²"},
    {"NozzleSolution::P0", "-", "Pa"},
    {"NozzleSolution::T0", "-", "K"},
    {"NozzleSolution::choked", "-", "-"},
    {"NozzleSolution::h0", "-", "J/kg"},
    {"NozzleSolution::inlet", "-", "-"},
    {"NozzleSolution::mdot", "-", "kg/s"},
    {"NozzleSolution::outlet", "-", "-"},
    {"NozzleSolution::profile", "-", "-"},
    {"NozzleSolution::x_throat", "-", "m"},
    {"FannoStation::M", "-", "-"},
    {"FannoStation::P", "-", "Pa"},
    {"FannoStation::Re", "-", "-"},
    {"FannoStation::T", "-", "K"},
    {"FannoStation::f", "-", "-"},
    {"FannoStation::h", "-", "J/kg"},
    {"FannoStation::rho", "-", "kg/m³"},
    {"FannoStation::s", "-", "J/(kg·K)"},
    {"FannoStation::u", "-", "m/s"},
    {"FannoStation::x", "-", "m"},
    {"FannoSolution::D", "-", "m"},
    {"FannoSolution::L", "-", "m"},
    {"FannoSolution::L_choke", "-", "m"},
    {"FannoSolution::Re_in", "-", "-"},
    {"FannoSolution::choked", "-", "-"},
    {"FannoSolution::f", "-", "-"},
    {"FannoSolution::f_avg", "-", "-"},
    {"FannoSolution::h0", "-", "J/kg"},
    {"FannoSolution::inlet", "-", "-"},
    {"FannoSolution::mdot", "-", "kg/s"},
    {"FannoSolution::outlet", "-", "-"},
    {"FannoSolution::profile", "-", "-"},
    {"solve_A_eff_from_mdot", "-", "-"},
    {"solve_P_back_from_mdot", "-", "-"},
    {"solve_P0_from_mdot", "-", "-"},
    {"ThrustResult::P_exit", "-", "Pa"},
    {"ThrustResult::mdot", "-", "kg/s"},
    {"ThrustResult::specific_impulse", "-", "s"},
    {"ThrustResult::thrust", "-", "N"},
    {"ThrustResult::thrust_coefficient", "-", "-"},
    {"ThrustResult::u_exit", "-", "m/s"},
    {"nozzle_thrust_cd", "-", "-"},
    {"TransferMatrix::T11", "-", "-"},
    {"TransferMatrix::T12", "-", "-"},
    {"TransferMatrix::T21", "-", "-"},
    {"TransferMatrix::T22", "-", "-"},
    {"CanAnnularGeometry::area_can", "-", "-"},
    {"CanAnnularGeometry::area_plenum", "-", "-"},
    {"CanAnnularGeometry::length_can", "-", "-"},
    {"CanAnnularGeometry::n_cans", "-", "-"},
    {"CanAnnularGeometry::radius_plenum", "-", "-"},
    {"BlochMode::frequency", "-", "-"},
    {"BlochMode::m_azimuthal", "-", "-"},
    {"BlochMode::n_cans", "-", "-"},
    {"BlochMode::symmetry_type", "-", "-"},
    {"can_annular_eigenmodes", "-", "combaero._core.BlochMode"},
    {"AnnularDuctGeometry::D_inner", "-", "m"},
    {"AnnularDuctGeometry::D_mean", "-", "m"},
    {"AnnularDuctGeometry::D_outer", "-", "m"},
    {"AnnularDuctGeometry::L", "-", "m"},
    {"AnnularDuctGeometry::area", "-", "m²"},
    {"AnnularDuctGeometry::circumference", "-", "m"},
    {"AnnularDuctGeometry::gap", "-", "m"},
    {"AnnularDuctGeometry::length", "-", "m"},
    {"AnnularDuctGeometry::n_azimuthal_max", "-", "-"},
    {"AnnularDuctGeometry::radius_inner", "-", "m"},
    {"AnnularDuctGeometry::radius_outer", "-", "m"},
    {"AnnularDuctGeometry::volume", "-", "m³"},
    {"AnnularMode::frequency", "-", "-"},
    {"AnnularMode::m_azimuthal", "-", "-"},
    {"AnnularMode::mode_type", "-", "-"},
    {"AnnularMode::n_axial", "-", "-"},
    {"annular_duct_eigenmodes", "-", "combaero._core.AnnularMode"},
    {"AcousticMedium::c", "-", "-"},
    {"AcousticMedium::rho", "-", "-"},
    {"LinerOrificeGeometry::Cd", "-", "-"},
    {"LinerOrificeGeometry::d_orifice", "-", "-"},
    {"LinerOrificeGeometry::l_orifice", "-", "-"},
    {"LinerOrificeGeometry::porosity", "-", "-"},
    {"LinerFlowState::u_bias", "-", "-"},
    {"LinerFlowState::u_grazing", "-", "-"},
    {"LinerCavity::depth", "-", "-"},
    {"overall_htc_wall_multilayer", "-", "-"},
    {"Tube::D", "-", "m"},
    {"Tube::L", "-", "m"},
    {"Tube::area", "-", "m²"},
    {"Tube::perimeter", "-", "m"},
    {"Tube::volume", "-", "m³"},
    {"Annulus::D_inner", "-", "m"},
    {"Annulus::D_mean", "-", "m"},
    {"Annulus::D_outer", "-", "m"},
    {"Annulus::L", "-", "m"},
    {"Annulus::area", "-", "m²"},
    {"Annulus::circumference", "-", "m"},
    {"Annulus::gap", "-", "m"},
    {"Annulus::length", "-", "m"},
    {"Annulus::n_azimuthal_max", "-", "-"},
    {"Annulus::radius_inner", "-", "m"},
    {"Annulus::radius_outer", "-", "m"},
    {"Annulus::volume", "-", "m³"},
    {"CanAnnularFlowGeometry::D_inner", "-", "m"},
    {"CanAnnularFlowGeometry::D_mean", "-", "m"},
    {"CanAnnularFlowGeometry::D_outer", "-", "m"},
    {"CanAnnularFlowGeometry::D_primary", "-", "m"},
    {"CanAnnularFlowGeometry::L", "-", "m"},
    {"CanAnnularFlowGeometry::L_primary", "-", "m"},
    {"CanAnnularFlowGeometry::L_transition", "-", "m"},
    {"CanAnnularFlowGeometry::area", "-", "m²"},
    {"CanAnnularFlowGeometry::area_primary", "-", "m²"},
    {"CanAnnularFlowGeometry::circumference", "-", "m"},
    {"CanAnnularFlowGeometry::gap", "-", "m"},
    {"CanAnnularFlowGeometry::length_total", "-", "m"},
    {"CanAnnularFlowGeometry::volume", "-", "m³"},
    {"CanAnnularFlowGeometry::volume_primary", "-", "m³"},
    {"CanAnnularFlowGeometry::volume_total", "-", "m³"},
    {"CanAnnularFlowGeometry::volume_transition", "-", "m³"},
    {"AcousticMode::frequency", "-", "Hz"},
    {"AcousticMode::label", "-", "-"},
    {"AcousticMode::n_axial", "-", "-"},
    {"AcousticMode::n_azimuthal", "-", "-"},
    {"AcousticProperties::frequency", "-", "Hz"},
    {"AcousticProperties::impedance", "-", "Pa·s/m"},
    {"AcousticProperties::particle_velocity", "-", "m/s"},
    {"AcousticProperties::spl", "-", "dB"},
    {"AcousticProperties::wavelength", "-", "m"},
    {"channel_dP", "-", "Pa"},
    {"channel_dP_mdot", "-", "Pa"},
    {"channel_velocity", "-", "m/s"},
    {"channel_mdot", "-", "kg/s"},
    {"channel_area", "-", "m²"},
    {"channel_volume", "-", "m³"},
    {"channel_roughness", "-", "m"},
    {"pressure_drop_channel", "-", "Pa"},
    {"nusselt_channel", "-", "-"},
    {"htc_channel", "-", "W/(m²·K)"},
    {"dynamic_pressure", "-", "Pa"},
    {"velocity_from_q", "-", "m/s"},
    {"OrificeGeometry::D", "-", "m"},
    {"OrificeGeometry::area", "-", "m²"},
    {"OrificeGeometry::beta", "-", "-"},
    {"OrificeGeometry::bevel", "-", "rad"},
    {"OrificeGeometry::d", "-", "m"},
    {"OrificeGeometry::is_valid", "-", "-"},
    {"OrificeGeometry::r", "-", "m"},
    {"OrificeGeometry::r_over_d", "-", "-"},
    {"OrificeGeometry::t", "-", "m"},
    {"OrificeGeometry::t_over_d", "-", "-"},
    {"OrificeState::Re_D", "-", "-"},
    {"OrificeState::Re_d", "-", "-"},
    {"OrificeState::dP", "-", "Pa"},
    {"OrificeState::mu", "-", "Pa·s"},
    {"OrificeState::rho", "-", "kg/m³"},
    {"CdCorrelation::Miller", "-", "-"},
    {"CdCorrelation::ReaderHarrisGallagher", "-", "-"},
    {"CdCorrelation::Stolz", "-", "-"},
    {"CdCorrelation::name", "-", "-"},
    {"CdCorrelation::value", "-", "-"},
    {"Cd_orifice", "-", "-"},
    {"orifice_mdot_Cd", "-", "kg/s"},
    {"orifice_dP_Cd", "-", "Pa"},
    {"orifice_Cd_from_measurement", "-", "kg/s"},
    {"orifice_K_from_Cd", "-", "-"},
    {"orifice_Cd_from_K", "-", "-"},
    {"orifice_thickness_correction", "-", "-"},
    {"OrificeFlowResult::Cd", "-", "-"},
    {"OrificeFlowResult::Re_D", "-", "-"},
    {"OrificeFlowResult::Re_d", "-", "-"},
    {"OrificeFlowResult::epsilon", "-", "-"},
    {"OrificeFlowResult::mdot", "-", "kg/s"},
    {"OrificeFlowResult::rho_corrected", "-", "kg/m³"},
    {"OrificeFlowResult::v", "-", "m/s"},
    {"orifice_flow_state", "-", "-"},
    // --- Base class / metadata properties ---
    {"FlowSolution::Cd", "-", "-"},
    {"FlowSolution::L_choke", "-", "m"},
    {"FlowSolution::M", "-", "-"},
    {"FlowSolution::P_out", "-", "Pa"},
    {"FlowSolution::Re", "-", "-"},
    {"FlowSolution::T_out", "-", "K"},
    {"FlowSolution::choked", "-", "-"},
    {"FlowSolution::dP", "-", "Pa"},
    {"FlowSolution::f", "-", "-"},
    {"FlowSolution::h0", "-", "J/kg"},
    {"FlowSolution::rho", "-", "kg/m³"},
    {"UnitInfo::input", "-", "-"},
    {"UnitInfo::output", "-", "-"},
    {"is_well_behaved", "v: any, lo: any, hi: any", "bool"},
    {"set_warning_handler", "handler: callable(str)", "-"},
    {"get_warning_handler", "-", "callable(str)"},
    {"suppress_warnings", "-", "context manager"},
};

inline constexpr std::size_t function_count =
    sizeof(function_units) / sizeof(Entry);

} // namespace combaero::units

#endif // UNITS_DATA_H
