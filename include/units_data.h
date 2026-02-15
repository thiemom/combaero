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
    // thermo.h - Species Data
    // -------------------------------------------------------------------------
    {"species_molar_mass",      "-",                                    "g/mol"},
    {"mwmix",                   "X: mol/mol",                           "g/mol"},
    {"mole_to_mass",            "X: mol/mol",                           "kg/kg"},
    {"mass_to_mole",            "Y: kg/kg",                             "mol/mol"},

    // -------------------------------------------------------------------------
    // thermo.h - Dimensionless NASA Polynomials
    // -------------------------------------------------------------------------
    {"cp_R",                    "T: K",                                 "- (Cp/R)"},
    {"h_RT",                    "T: K",                                 "- (H/(R*T))"},
    {"s_R",                     "T: K",                                 "- (S/R)"},
    {"g_over_RT",               "T: K",                                 "- (G/(R*T))"},

    // -------------------------------------------------------------------------
    // thermo.h - Per-Species Properties
    // -------------------------------------------------------------------------
    {"cp_species",              "T: K",                                 "J/(mol*K)"},
    {"h_species",               "T: K",                                 "J/mol"},
    {"s_species",               "T: K",                                 "J/(mol*K)"},

    // -------------------------------------------------------------------------
    // thermo.h - Mixture Properties (Molar Basis)
    // -------------------------------------------------------------------------
    {"cp",                      "T: K, X: mol/mol",                     "J/(mol*K)"},
    {"cv",                      "T: K, X: mol/mol",                     "J/(mol*K)"},
    {"h",                       "T: K, X: mol/mol",                     "J/mol"},
    {"u",                       "T: K, X: mol/mol",                     "J/mol"},
    {"s",                       "T: K, X: mol/mol, P: Pa, P_ref: Pa",   "J/(mol*K)"},
    {"dh_dT",                   "T: K, X: mol/mol",                     "J/(mol*K)"},
    {"ds_dT",                   "T: K, X: mol/mol",                     "J/(mol*K^2)"},
    {"dcp_dT",                  "T: K, X: mol/mol",                     "J/(mol*K^2)"},

    // -------------------------------------------------------------------------
    // thermo.h - Mixture Properties (Mass/Other Basis)
    // -------------------------------------------------------------------------
    {"density",                 "T: K, P: Pa, X: mol/mol",              "kg/m^3"},
    {"molar_volume",            "T: K, P: Pa",                          "m^3/mol"},
    {"specific_gas_constant",   "X: mol/mol",                           "J/(kg*K)"},
    {"isentropic_expansion_coefficient", "T: K, X: mol/mol",            "- (gamma)"},
    {"speed_of_sound",          "T: K, X: mol/mol",                     "m/s"},
    {"cp_mass",                 "T: K, X: mol/mol",                     "J/(kg*K)"},
    {"cv_mass",                 "T: K, X: mol/mol",                     "J/(kg*K)"},
    {"h_mass",                  "T: K, X: mol/mol",                     "J/kg"},
    {"s_mass",                  "T: K, X: mol/mol, P: Pa, P_ref: Pa",   "J/(kg*K)"},
    {"u_mass",                  "T: K, X: mol/mol",                     "J/kg"},
    {"air_properties",          "T: K, P: Pa, humidity: - (0-1)",       "AirProperties struct"},
    {"thermo_state",            "T: K, P: Pa, X: mol/mol, P_ref: Pa (default 101325)", "ThermoState struct"},
    {"transport_state",         "T: K, P: Pa, X: mol/mol",              "TransportState struct"},
    {"complete_state",          "T: K, P: Pa, X: mol/mol, P_ref: Pa (default 101325)", "CompleteState struct"},

    // -------------------------------------------------------------------------
    // acoustics.h - Acoustic Properties
    // -------------------------------------------------------------------------
    {"acoustic_properties",     "f: Hz, rho: kg/m^3, c: m/s, p_rms: Pa (default 20e-6), p_ref: Pa (default 20e-6)", "AcousticProperties struct"},
    {"wavelength",              "f: Hz, c: m/s",                        "m"},
    {"frequency_from_wavelength", "lambda: m, c: m/s",                  "Hz"},
    {"acoustic_impedance",      "rho: kg/m^3, c: m/s",                  "Pa*s/m"},
    {"sound_pressure_level",    "p_rms: Pa, p_ref: Pa (default 20e-6)", "dB"},
    {"particle_velocity",       "p: Pa, rho: kg/m^3, c: m/s",           "m/s"},

    // -------------------------------------------------------------------------
    // orifice.h - Orifice Flow Utilities
    // -------------------------------------------------------------------------
    {"orifice_flow",            "geom: OrificeGeometry, dP: Pa, T: K, P: Pa, mu: Pa*s, Z: - (default 1.0)", "OrificeFlowResult struct"},
    {"orifice_velocity_from_mdot", "mdot: kg/s, rho: kg/m^3, d: m, Z: - (default 1.0)", "m/s"},
    {"orifice_area_from_beta",  "D: m, beta: -",                        "m^2"},
    {"beta_from_diameters",     "d: m, D: m",                           "-"},
    {"orifice_Re_d_from_mdot",  "mdot: kg/s, d: m, mu: Pa*s",           "- (Reynolds number)"},

    // -------------------------------------------------------------------------
    // thermo.h - Inverse Solvers
    // -------------------------------------------------------------------------
    {"calc_T_from_h",           "h_target: J/mol, X: mol/mol",          "K"},
    {"calc_T_from_s",           "s_target: J/(mol*K), P: Pa, X: mol/mol", "K"},
    {"calc_T_from_cp",          "cp_target: J/(mol*K), X: mol/mol",     "K"},

    // -------------------------------------------------------------------------
    // transport.h - Transport Properties
    // -------------------------------------------------------------------------
    {"viscosity",               "T: K, P: Pa, X: mol/mol",              "Pa*s"},
    {"thermal_conductivity",    "T: K, P: Pa, X: mol/mol",              "W/(m*K)"},
    {"prandtl",                 "T: K, P: Pa, X: mol/mol",              "- (Pr)"},
    {"kinematic_viscosity",     "T: K, P: Pa, X: mol/mol",              "m^2/s"},
    {"thermal_diffusivity",     "T: K, P: Pa, X: mol/mol",              "m^2/s"},
    {"reynolds",                "T: K, P: Pa, X: mol/mol, V: m/s, L: m", "- (Re)"},
    {"peclet",                  "T: K, P: Pa, X: mol/mol, V: m/s, L: m", "- (Pe)"},

    // -------------------------------------------------------------------------
    // compressible.h - Nozzle Flow
    // -------------------------------------------------------------------------
    {"nozzle_flow",             "T0: K, P0: Pa, P_back: Pa, A_eff: m^2, X: mol/mol", "CompressibleFlowSolution"},
    {"nozzle_quasi1d",          "T0: K, P0: Pa, P_exit: Pa, x: m, A: m^2, X: mol/mol", "NozzleSolution"},
    {"nozzle_cd",               "T0: K, P0: Pa, P_exit: Pa, A: m^2, x: m, X: mol/mol", "NozzleSolution"},
    {"critical_pressure_ratio", "T0: K, P0: Pa, X: mol/mol",            "- (P*/P0)"},
    {"mach_from_pressure_ratio","T0: K, P0: Pa, P: Pa, X: mol/mol",     "- (M)"},
    {"mass_flux_isentropic",    "T0: K, P0: Pa, P: Pa, X: mol/mol",     "kg/(m^2*s)"},

    // -------------------------------------------------------------------------
    // compressible.h - Fanno Flow
    // -------------------------------------------------------------------------
    {"fanno_pipe",              "T_in: K, P_in: Pa, u_in: m/s, L: m, D: m, f: -, X: mol/mol", "FannoSolution"},
    {"fanno_max_length",        "T_in: K, P_in: Pa, u_in: m/s, D: m, f: -, X: mol/mol", "m"},

    // -------------------------------------------------------------------------
    // compressible.h - Thrust
    // -------------------------------------------------------------------------
    {"nozzle_thrust",           "NozzleSolution, P_amb: Pa",            "ThrustResult"},

    // -------------------------------------------------------------------------
    // incompressible.h - Bernoulli & Orifice
    // -------------------------------------------------------------------------
    {"bernoulli_P2",            "P1: Pa, v1: m/s, v2: m/s, rho: kg/m^3, dz: m, g: m/s^2", "Pa"},
    {"bernoulli_v2",            "P1: Pa, P2: Pa, v1: m/s, rho: kg/m^3, dz: m, g: m/s^2", "m/s"},
    {"orifice_mdot",            "P1: Pa, P2: Pa, A: m^2, Cd: -, rho: kg/m^3", "kg/s"},
    {"orifice_Q",               "P1: Pa, P2: Pa, A: m^2, Cd: -, rho: kg/m^3", "m^3/s"},
    {"orifice_velocity",        "P1: Pa, P2: Pa, rho: kg/m^3",          "m/s"},
    {"orifice_area",            "mdot: kg/s, P1: Pa, P2: Pa, Cd: -, rho: kg/m^3", "m^2"},
    {"orifice_dP",              "mdot: kg/s, A: m^2, Cd: -, rho: kg/m^3", "Pa"},
    {"expansibility_factor",    "beta: -, dP: Pa, P_upstream: Pa, kappa: -", "-"},

    // -------------------------------------------------------------------------
    // pipe_flow.h - Composite Pipe Flow Functions
    // -------------------------------------------------------------------------
    {"pressure_drop_pipe",      "T: K, P: Pa, X: mol/mol, v: m/s, D: m, L: m, roughness: m, correlation: str", "tuple(Pa, -, -)"},
    // -------------------------------------------------------------------------
    // geometry.h - Geometric Utilities
    // -------------------------------------------------------------------------
    {"pipe_area",               "D: m",                                 "m^2"},
    {"annular_area",            "D_outer: m, D_inner: m",               "m^2"},
    {"pipe_volume",             "D: m, L: m",                           "m^3"},
    {"pipe_roughness",          "material: str",                        "m"},
    {"standard_pipe_roughness", "-",                                    "dict[str, m]"},
    {"hydraulic_diameter",      "A: m^2, P_wetted: m",                  "m"},
    {"hydraulic_diameter_rect", "a: m, b: m",                           "m"},
    {"hydraulic_diameter_annulus", "D_outer: m, D_inner: m",            "m"},

    // -------------------------------------------------------------------------
    // friction.h - Friction Factor Correlations
    // -------------------------------------------------------------------------
    {"friction_haaland",        "Re: -, e_D: -",                        "- (f)"},
    {"friction_serghides",      "Re: -, e_D: -",                        "- (f)"},
    {"friction_colebrook",      "Re: -, e_D: -",                        "- (f)"},
    {"friction_petukhov",       "Re: -",                                "- (f)"},

    // -------------------------------------------------------------------------
    // heat_transfer.h - Heat Transfer Correlations
    // -------------------------------------------------------------------------
    {"nusselt_dittus_boelter",  "Re: -, Pr: -, heating: bool",          "- (Nu)"},
    {"nusselt_gnielinski",      "Re: -, Pr: -, f: - (optional)",        "- (Nu)"},
    {"nusselt_sieder_tate",     "Re: -, Pr: -, mu_ratio: -",            "- (Nu)"},
    {"nusselt_petukhov",        "Re: -, Pr: -, f: - (optional)",        "- (Nu)"},
    {"htc_from_nusselt",        "Nu: -, k: W/(m*K), L: m",              "W/(m^2*K)"},
    {"nusselt_pipe",            "State, V: m/s, D: m",                  "- (Nu)"},
    {"htc_pipe",                "State, V: m/s, D: m",                  "W/(m^2*K)"},
    {"htc_pipe",                "T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, correlation: str, heating: bool, mu_ratio: -, roughness: m", "tuple(W/(m^2*K), -, -)"},
    {"lmtd",                    "dT1: K, dT2: K",                       "K"},
    {"lmtd_counterflow",        "T_hot_in: K, T_hot_out: K, T_cold_in: K, T_cold_out: K", "K"},
    {"lmtd_parallelflow",       "T_hot_in: K, T_hot_out: K, T_cold_in: K, T_cold_out: K", "K"},
    {"overall_htc",             "h_values: W/(m^2*K), t_over_k: m^2*K/W", "W/(m^2*K)"},
    {"overall_htc_wall",        "h_inner, h_outer: W/(m^2*K), t_over_k_layers: m^2*K/W", "W/(m^2*K)"},
    {"overall_htc_tube",        "h_inner, h_outer: W/(m^2*K), t_wall: m, k_wall: W/(m*K)", "W/(m^2*K)"},
    {"thermal_resistance",      "h: W/(m^2*K), A: m^2",                 "K/W"},
    {"thermal_resistance_wall", "t: m, k: W/(m*K), A: m^2",             "K/W"},
    {"heat_rate",               "U: W/(m^2*K), A: m^2, dT: K",          "W"},
    {"heat_flux",               "U: W/(m^2*K), dT: K",                  "W/m^2"},
    {"heat_transfer_area",      "Q: W, U: W/(m^2*K), dT: K",            "m^2"},
    {"heat_transfer_dT",        "Q: W, U: W/(m^2*K), A: m^2",           "K"},
    {"wall_temperature_profile", "T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W", "K (vector)"},
    {"ntu",                     "U: W/(m^2*K), A: m^2, C_min: W/K",     "-"},
    {"capacity_ratio",          "C_min: W/K, C_max: W/K",               "-"},
    {"effectiveness_counterflow", "NTU: -, C_r: -",                     "-"},
    {"effectiveness_parallelflow", "NTU: -, C_r: -",                    "-"},
    {"heat_rate_from_effectiveness", "epsilon: -, C_min: W/K, T_hot_in: K, T_cold_in: K", "W"},
    {"heat_flux_from_T_at_edge", "T_measured: K, edge_idx, T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W", "W/m^2"},
    {"heat_flux_from_T_at_depth", "T_measured: K, depth: m, T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), thicknesses: m, conductivities: W/(m*K)", "W/m^2"},
    {"bulk_T_from_edge_T_and_q", "T_measured: K, edge_idx, q: W/m^2, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W, solve_for: str", "K"},
    {"dT_edge_dT_hot",          "edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W", "- (dT_edge/dT_hot)"},
    {"dT_edge_dT_cold",         "edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W", "- (dT_edge/dT_cold)"},
    {"dT_edge_dT_bulk",         "edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W", "- (dT/dT_hot, dT/dT_cold)"},
    {"dT_edge_dq",              "edge_idx, h_hot: W/(m^2*K), t_over_k: m^2*K/W", "K*m^2/W (dT_edge/dq)"},

    // -------------------------------------------------------------------------
    // acoustics.h - Acoustic Mode Analysis
    // -------------------------------------------------------------------------
    {"tube_axial_modes",        "Tube (L: m, D: m), c: m/s, upstream: BC, downstream: BC, n_max", "Hz (vector)"},
    {"annulus_axial_modes",     "Annulus (L: m, D_inner: m, D_outer: m), c: m/s, upstream: BC, downstream: BC, n_max", "Hz (vector)"},
    {"annulus_azimuthal_modes", "Annulus, c: m/s, m_max", "Hz (vector)"},
    {"annulus_modes",           "Annulus, c: m/s, upstream: BC, downstream: BC, n_max, m_max", "Hz (vector)"},
    {"modes_in_range",          "modes, f_min: Hz, f_max: Hz", "Hz (vector)"},
    {"closest_mode",            "modes, f_target: Hz", "AcousticMode"},
    {"min_mode_separation",     "modes", "Hz"},
    {"axial_mode_upstream",     "f0: Hz, M: -", "Hz"},
    {"axial_mode_downstream",   "f0: Hz, M: -", "Hz"},
    {"axial_mode_split",        "f0: Hz, M: -", "(Hz, Hz)"},
    {"helmholtz_frequency",     "V: m^3, A_neck: m^2, L_neck: m, c: m/s, end_correction: -", "Hz"},
    {"strouhal",                "f: Hz, L: m, u: m/s", "- (St)"},
    {"frequency_from_strouhal", "St: -, L: m, u: m/s", "Hz"},
    {"quarter_wave_frequency",  "L: m, c: m/s", "Hz"},
    {"half_wave_frequency",     "L: m, c: m/s", "Hz"},
    {"stokes_layer",            "nu: m^2/s, f: Hz", "m"},
    {"thermal_layer",           "alpha: m^2/s, f: Hz", "m"},
    {"effective_viscothermal_layer", "delta_nu: m, delta_kappa: m, gamma: -", "m"},
    {"helmholtz_Q",             "V: m^3, A_neck: m^2, L_neck: m, nu: m^2/s, alpha: m^2/s, gamma: -, f: Hz", "- (Q)"},
    {"tube_Q",                  "L: m, D: m, nu: m^2/s, alpha: m^2/s, gamma: -, f: Hz", "- (Q)"},
    {"damping_ratio",           "Q: -", "- (zeta)"},
    {"bandwidth",               "f0: Hz, Q: -", "Hz"},

    // -------------------------------------------------------------------------
    // geometry.h - Residence Time
    // -------------------------------------------------------------------------
    {"residence_time",          "V: m^3, Q: m^3/s", "s"},
    {"residence_time_mdot",     "V: m^3, mdot: kg/s, rho: kg/m^3", "s"},
    {"space_velocity",          "Q: m^3/s, V: m^3", "1/s"},

    // -------------------------------------------------------------------------
    // orifice.h - Discharge Coefficients
    // -------------------------------------------------------------------------
    {"Cd_sharp_thin_plate",     "OrificeGeometry, OrificeState",        "- (Cd)"},
    {"Cd_thick_plate",          "OrificeGeometry, OrificeState",        "- (Cd)"},
    {"Cd_rounded_entry",        "OrificeGeometry, OrificeState",        "- (Cd)"},
    {"Cd",                      "OrificeGeometry, OrificeState",        "- (Cd)"},

    // -------------------------------------------------------------------------
    // humidair.h - Humid Air Properties
    // -------------------------------------------------------------------------
    {"saturation_vapor_pressure", "T: K",                               "Pa"},
    {"vapor_pressure",          "T: K, RH: - (0-1)",                    "Pa"},
    {"humidity_ratio",          "T: K, P: Pa, RH: - (0-1)",             "kg/kg"},
    {"water_vapor_mole_fraction", "T: K, P: Pa, RH: - (0-1)",           "mol/mol"},
    {"humid_air_composition",   "T: K, P: Pa, RH: - (0-1)",             "mol/mol"},
    {"dewpoint",                "T: K, P: Pa, RH: - (0-1)",             "K"},
    {"relative_humidity_from_dewpoint", "T: K, Tdp: K, P: Pa",          "- (0-1)"},
    {"wet_bulb_temperature",    "T: K, P: Pa, RH: - (0-1)",             "K"},
    {"humid_air_enthalpy",      "T: K, P: Pa, RH: - (0-1)",             "J/kg"},
    {"humid_air_density",       "T: K, P: Pa, RH: - (0-1)",             "kg/m^3"},

    // -------------------------------------------------------------------------
    // combustion.h - Stoichiometry
    // -------------------------------------------------------------------------
    {"oxygen_required_per_mol_fuel",    "-",                            "mol O2/mol fuel"},
    {"oxygen_required_per_kg_fuel",     "-",                            "mol O2/kg fuel"},
    {"oxygen_required_per_mol_mixture", "X: mol/mol",                   "mol O2/mol mix"},
    {"oxygen_required_per_kg_mixture",  "X: mol/mol",                   "mol O2/kg mix"},
    {"dryair_required_per_mol_fuel",    "-",                            "mol air/mol fuel"},
    {"dryair_required_per_kg_fuel",     "-",                            "mol air/kg fuel"},
    {"dryair_required_per_mol_mixture", "X: mol/mol",                   "mol air/mol mix"},
    {"dryair_required_per_kg_mixture",  "X: mol/mol",                   "mol air/kg mix"},

    // -------------------------------------------------------------------------
    // combustion.h - Equivalence Ratio
    // -------------------------------------------------------------------------
    {"equivalence_ratio_mole",      "X: mol/mol",                       "- (phi)"},
    {"set_equivalence_ratio_mole",  "phi: -, X: mol/mol",               "mol/mol"},
    {"equivalence_ratio_mass",      "Y: kg/kg",                         "- (phi)"},
    {"set_equivalence_ratio_mass",  "phi: -, Y: kg/kg",                 "kg/kg"},

    // -------------------------------------------------------------------------
    // combustion.h - Mixture Fraction (Bilger)
    // -------------------------------------------------------------------------
    {"bilger_beta",                         "Y: kg/kg",                 "-"},
    {"bilger_mixture_fraction",             "Y: kg/kg",                 "- (Z)"},
    {"bilger_stoich_mixture_fraction_mass", "Y: kg/kg",                 "- (Z_st)"},
    {"equivalence_ratio_from_bilger_Z_mass","Z: -, Y: kg/kg",           "- (phi)"},
    {"bilger_Z_from_equivalence_ratio_mass","phi: -, Y: kg/kg",         "- (Z)"},

    // -------------------------------------------------------------------------
    // combustion.h - Complete Combustion
    // -------------------------------------------------------------------------
    {"complete_combustion_to_CO2_H2O",  "X: mol/mol",                   "mol/mol"},
    {"complete_combustion",             "State (T: K, P: Pa, X: mol/mol)", "State"},
    {"complete_combustion_isothermal",  "State",                        "State"},

    // -------------------------------------------------------------------------
    // combustion.h - Stream Solvers
    // -------------------------------------------------------------------------
    {"set_fuel_stream_for_phi",     "phi: -",                           "Stream"},
    {"set_fuel_stream_for_Tad",     "T_ad_target: K",                   "Stream"},
    {"set_fuel_stream_for_O2",      "X_O2_target: mol/mol",             "Stream"},
    {"set_fuel_stream_for_O2_dry",  "X_O2_dry_target: mol/mol",         "Stream"},
    {"set_fuel_stream_for_CO2",     "X_CO2_target: mol/mol",            "Stream"},
    {"set_fuel_stream_for_CO2_dry", "X_CO2_dry_target: mol/mol",        "Stream"},
    {"set_oxidizer_stream_for_Tad", "T_ad_target: K",                   "Stream"},
    {"set_oxidizer_stream_for_O2",  "X_O2_target: mol/mol",             "Stream"},
    {"set_oxidizer_stream_for_O2_dry", "X_O2_dry_target: mol/mol",      "Stream"},
    {"set_oxidizer_stream_for_CO2", "X_CO2_target: mol/mol",            "Stream"},
    {"set_oxidizer_stream_for_CO2_dry", "X_CO2_dry_target: mol/mol",    "Stream"},

    // -------------------------------------------------------------------------
    // equilibrium.h - Chemical Equilibrium
    // -------------------------------------------------------------------------
    {"wgs_equilibrium",                 "State (T: K, P: Pa, X: mol/mol)", "State"},
    {"wgs_equilibrium_adiabatic",       "State",                        "State"},
    {"smr_wgs_equilibrium",             "State",                        "State"},
    {"smr_wgs_equilibrium_adiabatic",   "State",                        "State"},
    {"reforming_equilibrium",           "State",                        "State"},
    {"reforming_equilibrium_adiabatic", "State",                        "State"},
    {"combustion_equilibrium",          "State",                        "State"},

    // -------------------------------------------------------------------------
    // state.h - State Properties
    // -------------------------------------------------------------------------
    {"State::T",        "-",    "K"},
    {"State::P",        "-",    "Pa"},
    {"State::X",        "-",    "mol/mol"},
    {"State::mw",       "-",    "g/mol"},
    {"State::cp",       "-",    "J/(mol*K)"},
    {"State::cv",       "-",    "J/(mol*K)"},
    {"State::h",        "-",    "J/mol"},
    {"State::u",        "-",    "J/mol"},
    {"State::s",        "-",    "J/(mol*K)"},
    {"State::rho",      "-",    "kg/m^3"},
    {"State::R",        "-",    "J/(kg*K)"},
    {"State::gamma",    "-",    "- (gamma)"},
    {"State::a",        "-",    "m/s"},
    {"State::mu",       "-",    "Pa*s"},
    {"State::k",        "-",    "W/(m*K)"},
    {"State::nu",       "-",    "m^2/s"},
    {"State::Pr",       "-",    "- (Pr)"},
    {"State::alpha",    "-",    "m^2/s"},

    // -------------------------------------------------------------------------
    // state.h - Stream Properties
    // -------------------------------------------------------------------------
    {"Stream::mdot",    "-",    "kg/s"},
};

inline constexpr std::size_t function_count = sizeof(function_units) / sizeof(Entry);

} // namespace combaero::units

#endif // UNITS_DATA_H
