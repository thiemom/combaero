#ifndef THERMO_H
#define THERMO_H

#include "state.h"
#include "composition.h"
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

namespace combaero::thermo {

// Universal gas constant [J/(mol*K)]
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
// h_RT (and g_over_RT, the equilibrium's g/RT) holds its value at this T
// below it (1/T); the kernels' stagnation-T inputs are floored here for an
// unphysical solver probe. Not a floor on any solved static T: those go to
// their root, however small (#935).
constexpr double T_FLOOR = 10.0;
// s_R holds its value below this T [K] (ln T), so that a non-positive probe
// T keeps s finite. It lies far below any isentropic root a double can carry
// (T/Tt = (P/Pt)^(R/cp) >= (1e-308)^(1/2.5) ~ 1e-123), so s is exact at every
// T a solve can reach.
constexpr double S_R_T_MIN = 1e-300;
// Below this T [K] the isentropic static-T Newton (solve_T_isentropic) steps
// in ln T, where s is ~ cp ln T and one step nearly lands on the root; above
// it the steps (and so the results) are those of the T Newton (#935). Equal
// to the former floor, so no iterate that stayed above it changes.
constexpr double T_LOG_STEP = 10.0;
// Pressure floor of the duct-state kernels [Pa]: an unphysical solver
// iterate (P <= 0) is evaluated there, finite, instead of returning NaN.
constexpr double P_FLOOR = 1.0;
// calc_T_from_cp searches [CP_INVERSE_T_MIN, CP_INVERSE_T_MAX] K: the span
// every species has NASA-9 data for (from 200 K; to 6000 K for H2O, NH3, CH4
// and the alkanes, the shortest). Beyond it the tails (#606/#607) level cp
// off toward a constant, so cp no longer identifies T. It scans outward from
// T_guess in steps of CP_INVERSE_SCAN_STEP K (both sides alternately) and
// polishes the first sign change of cp - target, the root nearest T_guess to
// within a step. Two roots closer than a step (a cp extremum within it of the
// target) can be missed (#868 R6).
constexpr double CP_INVERSE_T_MIN = 200.0;
constexpr double CP_INVERSE_T_MAX = 6000.0;
constexpr double CP_INVERSE_SCAN_STEP = 10.0;

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
double cp_species(std::size_t species_idx, double T);  // [J/(mol*K)]
double h_species(std::size_t species_idx, double T);   // [J/mol]
double s_species(std::size_t species_idx, double T);   // [J/(mol*K)] at P_ref

// Thermodynamic properties
double cp(double T, const std::vector<double>& X);
double h(double T, const std::vector<double>& X);
double s(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);
double cv(double T, const std::vector<double>& X);
double u(double T, const std::vector<double>& X);  // Internal energy [J/mol]
double density(double T, double P, const std::vector<double>& X);
double molar_volume(double T, double P);  // Ideal gas molar volume [m^3/mol]
double specific_gas_constant(const std::vector<double>& X);
double isentropic_expansion_coefficient(double T, const std::vector<double>& X);
double speed_of_sound(double T, const std::vector<double>& X);

// Mass-specific thermodynamic properties
double cp_mass(double T, const std::vector<double>& X);
// Secant cp over [T1, T2], (h(T1) - h(T2)) / (T1 - T2) [J/(kg K)], and its
// partials in T1 and T2; cp at T1 in the limit T2 -> T1. A stream's capacity
// rate over a span it may cross: the variable-property C of an epsilon-NTU
// exchange, which makes C_min (T_hot - T_cold) its exact maximum heat (#627).
std::tuple<double, double, double> cp_secant_mass_and_jacobian(double T1, double T2,
                                                               const std::vector<double>& X);
// The secant cp's gradient in the MASS fractions [J/(kg K)] per species
// (#717): a wall side's capacity rate |m| cp moves with its composition.
std::vector<double> cp_secant_mass_dY(double T1, double T2, const std::vector<double>& X);
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
// cp(T) is not monotone: the root nearest T_guess in [CP_INVERSE_T_MIN,
// CP_INVERSE_T_MAX] (see there); std::runtime_error when none lies in it,
// std::invalid_argument for a non-finite or non-positive input (#868 R6).
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
