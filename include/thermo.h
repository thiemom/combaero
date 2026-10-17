#ifndef THERMO_H
#define THERMO_H

#include "state.h"
#include "composition.h"
#include <cstddef>
#include <string>
#include <vector>

namespace combaero::thermo {

// Universal gas constant [J/(mol·K)]
constexpr double R_GAS = 8.31446261815324;

// Boltzmann constant [J/K]
constexpr double BOLTZMANN = 1.380649e-23;

// Avogadro's number [1/mol]
constexpr double AVOGADRO = 6.02214076e23;

// Species thermo beyond its NASA-9 data (#606). Past each end of a species'
// data (T_lo = 200 K; T_hi = 6000 or 20000 K) cp/R leaves the end interval's
// polynomial C2 (value, slope and curvature matched) and levels off (zero
// slope and curvature) at the end of a fade; h and s are its exact integrals,
// so cp stays positive and h, s monotone at every T for the solvers and the
// inverses.
//   Below the data the gas is ideal down to 0 K: over the whole [0, T_lo] the
//   quintic fades cp/R to the rigid-rotor value, translation plus fully
//   excited rotation (5/2 for an atom, 7/2 linear, 4 nonlinear), as the
//   vibrations freeze out -- translation alone (5/2) for H2, whose rotation
//   freezes inside the fade (Theta_rot = 85 K; H2O's largest, the next, is
//   40 K). h(0 K) is finite, so a stream's enthalpy has a floor: less is more
//   heat out than it holds, at no T. Not modelled: condensation (a real
//   stream liquefies well above 0 K).
//   Above the data the quartic levels off after a fade of TAIL_FRACTION_HIGH
//   T_hi at whatever constant it reaches; not physics (no dissociation).
constexpr double TAIL_FRACTION_HIGH = 0.1;
// cp/R of an ideal gas with translation and rotation fully excited and the
// vibrations frozen: the low tail's value at 0 K (CP_R_RIGID_ATOM, translation
// alone, also for H2).
constexpr double CP_R_RIGID_ATOM = 2.5;
constexpr double CP_R_RIGID_LINEAR = 3.5;
constexpr double CP_R_RIGID_NONLINEAR = 4.0;
// s_R and h_RT hold their value at this T below it (ln T, 1/T); cp_R and the
// dimensional h are defined at every T.
constexpr double T_FLOOR = 10.0;

} // namespace combaero::thermo

namespace combaero {

// Conversion factor from J/mol to J/kg
double J_per_mol_to_J_per_kg(double value, double molar_mass);

// Species metadata and lookup functions
std::string species_name(std::size_t idx);
std::size_t species_index_from_name(const std::string& name);
std::size_t num_species();

// Molar mass of a species [g/mol]
double species_molar_mass(std::size_t idx);
double species_molar_mass_from_name(const std::string& name);

// Which species thermo every function below evaluates (#549).
//   NASA9    : the NASA-9 polynomials as published. C0 in cp at
//              the interval joins (1000 K; 6000 K for N2, O2, AR, CO2, CO,
//              H2) but not C1: dcp/dT jumps 1-4 % at 1000 K (H2 20 %) and up
//              to 93 % at 6000 K, so every Jacobian using it steps there.
//   NASA9_C2 : the same polynomials with each join bridged, inside
//              |T - T_join| < 0.05 T_join, by a degree-7 cp/R that is C2 with
//              them at both window edges and whose h and s rises across the
//              window equal NASA-9's. Identical to NASA9 outside the windows
//              (h and s included); inside, cp within 1e-4 at 1000 K (5e-4 at
//              6000 K), h within 1e-5 cp T, s within 3e-5 R. Bridges are
//              generated (include/thermo_c2_refit_data.h). The default.
// Process-wide; thread-safe to read. Switch before solving, not during.
enum class ThermoModel { NASA9, NASA9_C2 };
void set_thermo_model(ThermoModel model);
ThermoModel thermo_model();

// NASA polynomial evaluations (dimensionless), with the tails beyond the data
// (thermo::TAIL_FRACTION_LOW/HIGH) under either model.
double cp_R(std::size_t species_idx, double T);
// d(cp/R)/dT of one species [1/K]; 0 on the constant tails.
double dcp_R_dT(std::size_t species_idx, double T);
double h_RT(std::size_t species_idx, double T);
double s_R(std::size_t species_idx, double T);
double g_over_RT(std::size_t species_idx, double T);

// Per-species dimensional properties
double cp_species(std::size_t species_idx, double T);  // [J/(mol·K)]
double h_species(std::size_t species_idx, double T);   // [J/mol]
double s_species(std::size_t species_idx, double T);   // [J/(mol·K)] at P_ref

// Thermodynamic properties
double cp(double T, const std::vector<double>& X);
double h(double T, const std::vector<double>& X);
double s(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);
double cv(double T, const std::vector<double>& X);
double u(double T, const std::vector<double>& X);  // Internal energy [J/mol]
double density(double T, double P, const std::vector<double>& X);
double molar_volume(double T, double P);  // Ideal gas molar volume [m³/mol]
double specific_gas_constant(const std::vector<double>& X);
double isentropic_expansion_coefficient(double T, const std::vector<double>& X);
double speed_of_sound(double T, const std::vector<double>& X);

// Mass-specific thermodynamic properties
double cp_mass(double T, const std::vector<double>& X);
double cv_mass(double T, const std::vector<double>& X);
double h_mass(double T, const std::vector<double>& X);
double s_mass(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);
double u_mass(double T, const std::vector<double>& X);

// Compute all air properties at once
AirProperties air_properties(double T, double P, double humidity = 0.0);

// Compute all thermodynamic properties at once
ThermoState thermo_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);

// Compute all thermodynamic and transport properties at once
CompleteState complete_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);

// Derivatives of thermodynamic properties with respect to temperature
double dh_dT(double T, const std::vector<double>& X);
double ds_dT(double T, const std::vector<double>& X);
double dcp_dT(double T, const std::vector<double>& X);
double dg_over_RT_dT(double T, const std::vector<double>& X);

// Inverse calculations (solving for temperature). h, s, u and the sv form are
// strictly increasing in T at every T (the tails keep cp > 0), so each has one
// root: safeguarded Newton (factor-4 steps until bracketed, then bisection
// when a step leaves the bracket). They throw std::runtime_error if they do
// not converge in max_iter, rather than return the last iterate.
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_u(double u_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);

// Inverse calculations for mass-specific targets
double calc_T_from_h_mass(double h_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s_mass(double s_mass_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_u_mass(double u_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_sv_mass(double s_mass_target, double v_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_sh_mass(double s_mass_target, double h_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);

// State-based overloads
double mwmix(const State& s);
double cp(const State& s);
double h(const State& s);
double s(const State& s);
double cv(const State& s);
double u(const State& s);
double density(const State& s);
double specific_gas_constant(const State& s);
double isentropic_expansion_coefficient(const State& s);
double speed_of_sound(const State& s);

} // namespace combaero

#endif // THERMO_H
