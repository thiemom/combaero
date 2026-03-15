#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "state.h"
// -----------------------------------------------------------------------------
// Network Solver Fast-Path Native Interface
// -----------------------------------------------------------------------------
// Functions in this module execute highly optimized `(f, J)` combined
// evaluations required by the global Python `scipy.optimize.root` solver.
//
// Every function returns a `std::tuple<double, double>`, corresponding to:
//    (Value, Exact Analytical Jacobian / Derivative)
//
// This architecture natively evaluates chain-rule analytical gradients where
// possible, and tightly bounded finite-differences where strictly necessary,
// entirely shielding the Python runtime from grid discontinuity noise and
// locking the C++/Python return overhead to ~0.2 microseconds.

namespace combaero {
namespace solver {

enum class CorrelationValidity : std::uint8_t { VALID, EXTRAPOLATED, INVALID };

template <typename T> struct CorrelationResult {
  T result;
  CorrelationValidity status;
  std::string message;
};

// -----------------------------------------------------------------------------
// Stream Definitions for Network Solvers
// -----------------------------------------------------------------------------

struct Stream {
  double m_dot;
  double T;
  double P_total;
  std::vector<double> Y;
};

struct StreamJacobian {
  double d_mdot;
  double d_T;
  double d_P_total;
  std::vector<double> d_Y;
};

struct MixerResult {
  double T_mix;
  double P_total_mix;
  std::vector<double> Y_mix;
  std::vector<StreamJacobian> dT_mix_d_stream;
  std::vector<StreamJacobian> dP_total_mix_d_stream;
  std::vector<std::vector<StreamJacobian>> dY_mix_d_stream;
};

struct ChamberResult {
  std::vector<double> residuals;
  // Local Jacobian: map from eq_idx to {var_name: derivative}
  // var_name is a string like "P", "P_total", "T", "Y[0]"...
  std::vector<std::map<std::string, double>> local_jacobian;
  // Stream Jacobian: map from eq_idx to list of StreamJacobian (one per
  // upstream stream)
  std::vector<std::vector<StreamJacobian>> stream_jacobian;
};

struct OrificeResult {
  double m_dot_calc;
  double d_mdot_dP_total_up;
  double d_mdot_dP_static_down;
  double d_mdot_dP_static_up;
  double d_mdot_dT_up;
  std::vector<double> d_mdot_dY_up;
};

struct PipeResult {
  double dP_calc;
  double d_dP_d_mdot;
  double d_dP_dP_static_up;
  double d_dP_dT_up;
  std::vector<double> d_dP_dY_up;
};

// -----------------------------------------------------------------------------
// 1. Incompressible Flow Components
// -----------------------------------------------------------------------------

// Calculate Mass Flow and its derivative with respect to differential pressure.
//
// Method: Analytical Chain Rule
// Inputs:
//   dP   : Differential pressure drop across the orifice [Pa]
//   rho  : Fluid density [kg/m³]
//   Cd   : Discharge coefficient [-]
//   beta : Orifice beta ratio (d/D) [-]
//
// Returns:
//   tuple(mdot, d_mdot_d_dP)
//   mdot        : Mass flow rate [kg/s]
//   d_mdot_d_dP : Jacobian gradient [kg/(s*Pa)]
std::tuple<double, double> orifice_mdot_and_jacobian(double dP, double rho,
                                                     double Cd, double area,
                                                     double beta);

// Full orifice evaluation with all derivatives.
OrificeResult orifice_residuals_and_jacobian(double m_dot, double P_total_up,
                                             double P_static_up, double T_up,
                                             const std::vector<double> &Y_up,
                                             double P_static_down, double Cd,
                                             double area, double beta);

// Calculate Pressure Drop and its derivative with respect to velocity.
//
// Method: Analytical Chain Rule
// Inputs:
//   v   : Fluid velocity [m/s]
//   rho : Fluid density [kg/m³]
//   K   : Forward loss coefficient [-]
//
// Returns:
//   tuple(dP, d_dP_d_v)
//   d_dP_d_v : Jacobian gradient [Pa/(m/s)]
std::tuple<double, double> pressure_loss_and_jacobian(double v, double rho,
                                                       double K);

// Full pipe evaluation with all derivatives.
PipeResult pipe_residuals_and_jacobian(double m_dot, double P_total_up,
                                       double P_static_up, double T_up,
                                       const std::vector<double> &Y_up,
                                       double P_static_down, double L, double D,
                                       double roughness,
                                       const std::string &friction_model);


// Calculate Lossless Connection Ideal Pressure and Jacobian
//
// Method: Analytical Difference
// Inputs:
//   P_in  : Total pressure entering the element [Pa]
//   P_out : Total pressure exiting the element [Pa]
//
// Returns:
//   tuple(Residual, d_Residual_d_P_out)
//   Residual           : Difference between P_in and P_out [Pa]
//   d_Residual_d_P_out : Jacobian gradient [-]
std::tuple<double, double> lossless_pressure_and_jacobian(double P_in,
                                                          double P_out);

// -----------------------------------------------------------------------------
// 2. Heat Transfer Components
// -----------------------------------------------------------------------------

// Calculate Nusselt Number and its derivative with respect to Reynolds Number.
// Uses the Dittus-Boelter fundamental relation: Nu = 0.023 * Re^0.8 * Pr^n
//
// Method: Analytical Chain Rule
// Inputs:
//   Re      : Reynolds number [-]
//   Pr      : Prandtl number [-]
//   heating : True if fluid is being heated (n=0.4), False if cooled (n=0.3)
//
// Returns:
//   tuple(Nu, d_Nu_d_Re)
//   Nu        : Nusselt number [-]
//   d_Nu_d_Re : Jacobian gradient [-]
CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_dittus_boelter(double Re, double Pr, bool heating = true);

// Calculate Nusselt Number and its derivative with respect to Reynolds Number
// using the Gnielinski correlation.
//
// Method: Finite Central Difference
CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_gnielinski(double Re, double Pr, double f);

// Calculate Nusselt Number and its derivative with respect to Reynolds Number
// using the Sieder-Tate correlation.
//
// Method: Analytical Chain Rule
CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_sieder_tate(double Re, double Pr, double mu_ratio);

// Calculate Nusselt Number and its derivative with respect to Reynolds Number
// using the Petukhov correlation.
//
// Method: Analytical Chain Rule
CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_petukhov(double Re, double Pr, double f);

// Calculate Friction Factor and its derivative with respect to Reynolds Number.
// Uses the Haaland explicit pipe friction relation.
//
// Method: Analytical Chain Rule
// Inputs:
//   Re  : Reynolds number [-]
//   e_D : Relative roughness [-]
//
// Returns:
//   tuple(f, d_f_d_Re)
//   f        : Darcy friction factor [-]
//   d_f_d_Re : Jacobian gradient [-]
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_haaland(double Re, double e_D);

// Calculate Friction Factor and its derivative using the Serghides explicit
// approximation to Colebrook-White (Steffensen acceleration).
//
// Method: Analytical Chain Rule through A, B, C intermediates
// Inputs:
//   Re  : Reynolds number [-]
//   e_D : Relative roughness [-]
//
// Returns:
//   tuple(f, d_f_d_Re)
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_serghides(double Re, double e_D);

// Calculate Friction Factor and its derivative using the Colebrook-White
// implicit equation (Newton-Raphson converged internally).
//
// Method: Implicit differentiation of converged Colebrook root
// Inputs:
//   Re  : Reynolds number [-]
//   e_D : Relative roughness [-]
//
// Returns:
//   tuple(f, d_f_d_Re)
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_colebrook(double Re, double e_D);

// Calculate Friction Factor and its derivative using the Petukhov smooth-pipe
// correlation: f = (coeff_a * ln(Re) - coeff_b)^(-2)
//
// Method: Analytical Chain Rule
// Inputs:
//   Re : Reynolds number [-]  (smooth pipe: e_D not used)
//
// Returns:
//   tuple(f, d_f_d_Re)
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_petukhov(double Re);

// Convenience dispatcher: selects the correct friction_and_jacobian_* overload
// based on a string tag matching FrictionModelLiteral from Python.
//
// tag must be one of: "haaland", "serghides", "colebrook", "petukhov"
//
// Returns:
//   tuple(f, d_f_d_Re)
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian(const std::string &tag, double Re, double e_D);

// -----------------------------------------------------------------------------
// 3. Cooling Correlations
// -----------------------------------------------------------------------------

// Pin Fin Nusselt Number and Jacobian wrt Re_d
CorrelationResult<std::tuple<double, double>>
pin_fin_nusselt_and_jacobian(double Re_d, double Pr, double L_D, double S_D,
                             double X_D, bool is_staggered = true);

// Pin Fin Friction and Jacobian wrt Re_d
CorrelationResult<std::tuple<double, double>>
pin_fin_friction_and_jacobian(double Re_d, bool is_staggered = true);

// Dimple Nusselt Enhancement and Jacobian wrt Re_Dh
CorrelationResult<std::tuple<double, double>>
dimple_nusselt_enhancement_and_jacobian(double Re_Dh, double d_Dh, double h_d,
                                        double S_d);

// Dimple Friction Multiplier and Jacobian wrt Re_Dh
CorrelationResult<std::tuple<double, double>>
dimple_friction_multiplier_and_jacobian(double Re_Dh, double d_Dh, double h_d);

// Rib Enhancement Factor (High-Re) and Jacobian wrt Re
CorrelationResult<std::tuple<double, double>>
rib_enhancement_factor_high_re_and_jacobian(double e_D, double pitch_to_height,
                                            double alpha_deg, double Re);

// Impingement Nusselt Number and Jacobian wrt Re_jet
CorrelationResult<std::tuple<double, double>>
impingement_nusselt_and_jacobian(double Re_jet, double Pr, double z_D,
                                 double x_D = 0.0, double y_D = 0.0);

// Film Cooling Effectiveness and Jacobian wrt Blowing Ratio M
CorrelationResult<std::tuple<double, double>>
film_cooling_effectiveness_and_jacobian(double x_D, double M, double DR,
                                        double alpha_deg);

// Effusion Effectiveness and Jacobian wrt Blowing Ratio M
CorrelationResult<std::tuple<double, double>>
effusion_effectiveness_and_jacobian(double x_D, double M, double DR,
                                    double porosity, double s_D,
                                    double alpha_deg);

// -----------------------------------------------------------------------------
// 4. Thermodynamic & Transport Components
// -----------------------------------------------------------------------------

// Calculate Mass Density and its analytical derivatives w.r.t. T and P.
// Method: Ideal Gas Law Analytical Chain Rule
// Inputs:
//   T : Temperature [K]
//   P : Pressure [Pa]
//   X : Mole fractions [-]
//
// Returns:
//   tuple(rho, d_rho_d_T, d_rho_d_P)
//   rho       : Density [kg/m³]
//   d_rho_d_T : Jacobian gradient w.r.t Temperature [(kg/m³)/K]
//   d_rho_d_P : Jacobian gradient w.r.t Pressure [(kg/m³)/Pa]
std::tuple<double, double, double>
density_and_jacobians(double T, double P, const std::vector<double> &X);

// Calculate Mass Enthalpy and its analytical derivative w.r.t. T.
// Method: dh/dT is identically the mass-specific heat capacity (cp_mass).
// Inputs:
//   T : Temperature [K]
//   X : Mole fractions [-]
//
// Returns:
//   tuple(h, d_h_d_T)
//   h       : Mass Enthalpy [J/kg]
//   d_h_d_T : Jacobian gradient w.r.t Temperature [J/(kg*K)]
std::tuple<double, double> enthalpy_and_jacobian(double T,
                                                 const std::vector<double> &X);

// Calculate Dynamic Viscosity and its derivatives w.r.t. T and P.
// Method: Central Finite Difference (tightly bracketed)
// Inputs:
//   T : Temperature [K]
//   P : Pressure [Pa]
//   X : Mole fractions [-]
//
// Returns:
//   tuple(mu, d_mu_d_T, d_mu_d_P)
//   mu       : Dynamic viscosity [Pa*s]
//   d_mu_d_T : Jacobian gradient w.r.t Temperature [(Pa*s)/K]
//   d_mu_d_P : Jacobian gradient w.r.t Pressure [(Pa*s)/Pa]
std::tuple<double, double, double>
viscosity_and_jacobians(double T, double P, const std::vector<double> &X);

// -----------------------------------------------------------------------------
// 5. Stagnation State Conversions
// -----------------------------------------------------------------------------

// Calculate Mach number and its derivative wrt Velocity
// Method: Analytical Chain Rule
std::tuple<double, double>
mach_number_and_jacobian_v(double v, double T, const std::vector<double> &X);

// Calculate Adiabatic Wall Temperature and its derivative wrt Velocity
// Method: Analytical Chain Rule
std::tuple<double, double>
T_adiabatic_wall_and_jacobian_v(double T_static, double v, double T, double P,
                                const std::vector<double> &X,
                                bool turbulent = true);

// Calculate Stagnation Temperature from Static and its derivative wrt Mach
// Method: Central Finite Difference
std::tuple<double, double>
T0_from_static_and_jacobian_M(double T, double M, const std::vector<double> &X);

// Calculate Stagnation Pressure from Static and its derivative wrt Mach
// Method: Central Finite Difference
std::tuple<double, double>
P0_from_static_and_jacobian_M(double P, double T, double M,
                              const std::vector<double> &X);

// -----------------------------------------------------------------------------
// 6. Combustion Interfaces
// -----------------------------------------------------------------------------

// Calculate Adiabatic Flame Temperature (Complete Combustion) and its
// derivative wrt inlet temperature.
// Method: Central Finite Difference tightly bracketed.
// Inputs:
//   T_in : Inlet temperature [K]
//   P    : Pressure [Pa]
//   X_in : Inlet mole fractions [-]
// Returns:
//   tuple(T_ad, dT_ad_dT_in, X_products)
//   T_ad        : Adiabatic flame temperature [K]
//   dT_ad_dT_in : Jacobian gradient [-]
//   X_products  : Mole fractions of products [-]
std::tuple<double, double, std::vector<double>>
adiabatic_T_complete_and_jacobian_T(double T_in, double P,
                                    const std::vector<double> &X_in);

// Calculate Adiabatic Flame Temperature (Equilibrium Combustion) and its
// derivatives wrt inlet temperature and pressure.
// Method: Central Finite Difference tightly bracketed.
// Inputs:
//   T_in : Inlet temperature [K]
//   P    : Pressure [Pa]
//   X_in : Inlet mole fractions [-]
// Returns:
//   tuple(T_ad, dT_ad_dT_in, dT_ad_dP, X_products)
//   T_ad        : Adiabatic flame temperature [K]
//   dT_ad_dT_in : Jacobian gradient wrt Inlet Temperature [-]
//   dT_ad_dP    : Jacobian gradient wrt Pressure [K/Pa]
//   X_products  : Mole fractions of products [-]
std::tuple<double, double, double, std::vector<double>>
adiabatic_T_equilibrium_and_jacobians(double T_in, double P,
                                      const std::vector<double> &X_in);

// -----------------------------------------------------------------------------
// 7. Stream-Based Network Solvers
// -----------------------------------------------------------------------------


// Mix generalized non-reacting incoming streams into a single outgoing state,
// and compute all analytical Jacobians w.r.t upstream mass flows, Temperatures,
// and mass fractions (Y_i).
MixerResult
mixer_from_streams_and_jacobians(const std::vector<Stream> &streams);

// Combust multiple generalized incoming streams to exactly complete
// equilibrium, outputting the theoretical adiabatic flame temperature and
// product mass fractions. Computes all exact analytical Jacobians w.r.t
// upstream conditions.
MixerResult adiabatic_T_complete_and_jacobian_T_from_streams(
    const std::vector<Stream> &streams, double P);

// Combust multiple generalized incoming streams to generalized equilibrium,
// outputting the theoretical adiabatic flame temperature and product mass
// fractions. Computes all exact analytical Jacobians w.r.t upstream conditions.
MixerResult adiabatic_T_equilibrium_and_jacobians_from_streams(
    const std::vector<Stream> &streams, double P);

} // namespace solver
} // namespace combaero

#endif // SOLVER_INTERFACE_H
