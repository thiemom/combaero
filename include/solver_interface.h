#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "state.h"
#include "combustion.h"
// -----------------------------------------------------------------------------
// Network Solver Fast-Path Native Interface
// -----------------------------------------------------------------------------
// Functions in this module execute highly optimized `(f, J)` combined
// evaluations required by the global Python `scipy.optimize.root` solver.
// Note: Concepts currently named "Pipe" are scheduled for renaming to "Channel"
// to account for advanced non-circular geometries and aerospace cooling applications.
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
  double dT_mix_d_delta_h = 0.0;
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

struct ChannelResult {
  double dP_calc;
  double d_dP_d_mdot;
  double d_dP_dP_static_up;
  double d_dP_dT_up;
  std::vector<double> d_dP_dY_up;
};

using PipeResult = ChannelResult;

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
//   tuple(mdot, d_mdot_d_dP, d_mdot_d_rho)
//   mdot         : Mass flow rate [kg/s]
//   d_mdot_d_dP  : Jacobian gradient w.r.t. pressure difference [kg/(s*Pa)]
//   d_mdot_d_rho : Jacobian gradient w.r.t. density [kg*m³/(s*kg)] = [m³/s]
std::tuple<double, double, double> orifice_mdot_and_jacobian(double dP, double rho,
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

// Full channel evaluation with all derivatives.
ChannelResult channel_residuals_and_jacobian(double m_dot, double P_total_up,
                                             double P_static_up, double T_up,
                                             const std::vector<double> &Y_up,
                                             double P_static_down, double L,
                                             double D, double roughness,
                                             const std::string &friction_model,
                                             double f_multiplier = 1.0);



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
// 1b. Compressible Flow Components
// -----------------------------------------------------------------------------

// Compressible orifice flow using isentropic nozzle model.
// Applies discharge coefficient and velocity-of-approach factor to nozzle_flow.
//
// Method: Finite Difference Jacobians (nozzle_flow is iterative)
// Inputs:
//   T0   : Stagnation temperature [K]
//   P0   : Stagnation pressure [Pa]
//   P_back : Back pressure [Pa]
//   X    : Mole fractions [-]
//   Cd   : Discharge coefficient [-]
//   area : Geometric throat area [m²]
//   beta : Orifice beta ratio (d/D) [-]
//
// Returns:
//   tuple(mdot, d_mdot_dP0, d_mdot_dP_back, d_mdot_dT0)
//   mdot           : Mass flow rate [kg/s]
//   d_mdot_dP0     : Jacobian w.r.t. stagnation pressure [kg/(s*Pa)]
//   d_mdot_dP_back : Jacobian w.r.t. back pressure [kg/(s*Pa)] (smooth through choked transition)
//   d_mdot_dT0     : Jacobian w.r.t. stagnation temperature [kg/(s*K)]
std::tuple<double, double, double, double> orifice_compressible_mdot_and_jacobian(
    double T0, double P0, double P_back,
    const std::vector<double>& X,
    double Cd, double area, double beta);

// Full compressible orifice evaluation with all derivatives for network solver.
OrificeResult orifice_compressible_residuals_and_jacobian(
    double m_dot, double P_total_up, double T_up,
    const std::vector<double>& Y_up,
    double P_static_down, double Cd, double area, double beta);

// Infeasible-flow barrier constants for the solver-facing compressible
// channel residual.
//
// fanno_channel_rough reports choked flow by truncating the march: at
// M_in >= 1 it returns outlet == inlet (dP = 0, zero slope), and for
// in-channel choking (L_choke < L) the reported friction drop collapses
// toward zero as the requested velocity approaches sonic. A Newton solver
// reads that as "no resistance", which turns the infeasible region into a
// flat attractor with exact spurious network roots (dead branch arms,
// reversed supplies). channel_compressible_mdot_and_jacobian therefore
// replaces the truncated result with a monotone barrier:
//   in-channel choking: dP = dP_truncated + kappa * P_in * (1 - L_choke/L)
//   supersonic inlet:   dP = kappa * P_in + (beta/2) * rho_in * (u^2 - a^2)
// Both branches are continuous at their region boundaries and increase with
// the depth of infeasibility, so the residual always pushes m_dot back
// toward the feasible branch.
constexpr double kChannelChokeBarrierKappa = 1.0;
constexpr double kChannelChokeBarrierBeta = 2.0;

// Compressible channel flow using Fanno model with variable friction.
//
// Method: Finite Difference Jacobians (fanno_channel_rough is iterative)
// Inputs:
//   T_in   : Inlet static temperature [K]
//   P_in   : Inlet static pressure [Pa]
//   u_in   : Inlet velocity [m/s]
//   X      : Mole fractions [-]
//   L      : Pipe length [m]
//   D      : Pipe diameter [m]
//   roughness : Absolute roughness [m]
//   friction_model : Friction correlation name
//
// Returns:
//   tuple(dP, d_dP_dP_in, d_dP_dT_in, d_dP_du_in)
//   dP         : Pressure drop [Pa]
//   d_dP_dP_in : Jacobian w.r.t. inlet pressure [Pa/Pa] = [-]
//   d_dP_dT_in : Jacobian w.r.t. inlet temperature [Pa/K]
//   d_dP_du_in : Jacobian w.r.t. inlet velocity [Pa/(m/s)]
std::tuple<double, double, double, double> channel_compressible_mdot_and_jacobian(
    double T_in, double P_in, double u_in,
    const std::vector<double>& X,
    double L, double D, double roughness,
    const std::string& friction_model,
    double f_multiplier = 1.0,
    bool compute_jacobians = true);


// Full compressible channel evaluation with all derivatives for network solver.
ChannelResult channel_compressible_residuals_and_jacobian(
    double m_dot, double P_total_up, double T_up,
    const std::vector<double>& Y_up,
    double P_static_down, double L, double D, double roughness,
    const std::string& friction_model,
    double f_multiplier = 1.0);


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
// Q [W]: absolute heat transfer rate (positive = heating)
// fraction [-]: relative to total enthalpy flow (positive = heating)
MixerResult
mixer_from_streams_and_jacobians(const std::vector<Stream> &streams,
                                double Q = 0.0, double fraction = 0.0);

// Combust multiple generalized incoming streams to exactly complete
// equilibrium, outputting the theoretical adiabatic flame temperature and
// product mass fractions. Computes all exact analytical Jacobians w.r.t
// upstream conditions.
// Q [W]: absolute heat transfer rate (positive = heating)
// fraction [-]: relative to total enthalpy flow (positive = heating)
MixerResult adiabatic_T_complete_and_jacobian_T_from_streams(
    const std::vector<Stream> &streams, double P, double Q = 0.0, double fraction = 0.0);

// Combust multiple generalized incoming streams to generalized equilibrium,
// outputting the theoretical adiabatic flame temperature and product mass
// fractions. Computes all exact analytical Jacobians w.r.t upstream conditions.
// Q [W]: absolute heat transfer rate (positive = heating)
// fraction [-]: relative to total enthalpy flow (positive = heating)
MixerResult adiabatic_T_equilibrium_and_jacobians_from_streams(
    const std::vector<Stream> &streams, double P, double Q = 0.0, double fraction = 0.0);
// Full combustor evaluation including pressure loss analytical chain rule.
MixerResult combustor_residuals_and_jacobians(
    const std::vector<Stream> &streams, double P, double Q, double fraction,
    const PressureLossCorrelation &pressure_loss, bool use_equilibrium = false);

// Momentum chamber residual: P_total = P_static + 0.5 * rho * v^2
// where v = m_dot / (rho * A)
// Returns residual and Jacobian entries for P, P_total, and m_dot
struct MomentumChamberResult {
  double residual;
  double d_res_dP;
  double d_res_dP_total;
  double d_res_dmdot;
};

MomentumChamberResult momentum_chamber_residual_and_jacobian(
    double P, double P_total, double m_dot, double T,
    const std::vector<double> &Y, double area);

// -----------------------------------------------------------------------------
// Multi-port chamber (momentum-CV junction).
// PDF spec: docs/junction/momentum cv implementation guide.pdf, Section 2.
// Junction owns one scalar P_jct. Emits N per-port impulse-function residuals
//   R_mom,i = (P_i + rho_i * u_i^2) - P_jct = 0
// plus a global mass residual
//   R_mass = sum_i mdot_i = 0
// Sign convention: mdot_i > 0 means flow OUT of the junction through port i.
// u_i^2 is direction-invariant, so the impulse residual is sign-free.
//
// The junction emits 0 empirical loss; per-port turning/contraction losses go
// on companion BorderCarnotLossElement instances bolted onto lateral ports.

// Per-port impulse-residual Jacobian.
// dR_mom_i/dP_jct = -1 (constant, omitted from struct; caller adds).
// Mass-row Jacobian dR_mass/dmdot_i = +1 (constant, omitted).
struct PortImpulseJacobian {
  double dR_dP;     // d(R_mom_i)/d(P_i)
  double dR_dmdot;  // d(R_mom_i)/d(mdot_i)
  double dR_dT;     // d(R_mom_i)/d(T_i)
};

struct MultiPortChamberResult {
  std::vector<double> impulse_residuals;     // size N
  std::vector<PortImpulseJacobian> port_jac; // size N

  // Cross-coupling Jacobian: each non-axial port's impulse residual depends
  // on the AXIAL REFERENCE port (port 0) state through rho_0 and u_0 in the
  // -2*sin(theta_i)*cos((3/4)*theta_i)*rho_0*u_0*u_i term. Zero at axial
  // ports (sin(theta_i)=0) and at port 0 itself.
  std::vector<double> cross_dR_dP_axial;     // size N: d(R_mom_i)/d(P_0)
  std::vector<double> cross_dR_dT_axial;     // size N: d(R_mom_i)/d(T_0)
  std::vector<double> cross_dR_dmdot_axial;  // size N: d(R_mom_i)/d(mdot_0)

  double mass_residual;                      // sum_i mdot_i
};

// Inputs (per-port vectors must all be length N >= 2):
//   P_jct       : junction internal static pressure [Pa]
//   P           : port-face static pressures [Pa]
//   mdot        : port mass flows [kg/s], positive = out of junction
//   T           : port temperatures [K]
//   Y           : per-port mass-fraction vectors (each length n_species)
//   A           : port cross-section areas [m^2]
//   theta_rad   : per-port geometric branch angles [rad]; the impulse term
//                 picks up a cos^2(theta) projection so the residual recovers
//                 1D-duct momentum at theta=0 (straight) and decouples
//                 entirely at theta=pi/2 (perpendicular branch). This is the
//                 angle-aware extension of the PDF Section 2.2 spec; at
//                 theta=0 for every port it reduces to the original
//                 P + rho*u^2 form.
MultiPortChamberResult multi_port_chamber_residuals_and_jacobian(
    double P_jct,
    const std::vector<double> &P,
    const std::vector<double> &mdot,
    const std::vector<double> &T,
    const std::vector<std::vector<double>> &Y,
    const std::vector<double> &A,
    const std::vector<double> &theta_rad);

// -----------------------------------------------------------------------------
// Border-Carnot loss element (companion to multi-port chamber).
// Two-port in-line element: stagnation drop across the port.
//   R = Pt_in - Pt_out - L(delta_geom) * 0.5 * rho_in * u_in^2 = 0
//   L = 4 * (1 - cos((3/4) * delta_geom))^2     (PDF Section 3.1)
// Sign-free (mdot^2 in the dynamic head); matches Hager xi_l and Bassett K_inc
// at M -> 0 on a sharp-edged lateral with the Hager (3/4) correction.
struct BorderCarnotLossResult {
  double residual;
  double d_res_dPt_in;
  double d_res_dPt_out;
  double d_res_dP_in;  // through rho(T_in, P_in, Y_in)
  double d_res_dT_in;  // through rho
  double d_res_dmdot;
};

BorderCarnotLossResult border_carnot_loss_residual_and_jacobian(
    double mdot,
    double Pt_in, double Pt_out,
    double P_in, double T_in,
    const std::vector<double> &Y_in,
    double area, double delta_geom);

// -----------------------------------------------------------------------------
// 4. Area Change Elements
// -----------------------------------------------------------------------------

// Result for area change elements in the network solver.
// All derivatives are w.r.t. the solver state variables (m_dot, P, T, Y).
struct AreaChangeElementResult {
  double dP_calc;            // Total pressure drop [Pa], signed
  double d_dP_d_mdot;        // d(dP)/d(m_dot)       [Pa*s/kg]
  double d_dP_dP_static_up;  // d(dP)/d(P_static_up) [Pa/Pa]
  double d_dP_dT_up;         // d(dP)/d(T_up)        [Pa/K]
  std::vector<double> d_dP_dY_up;  // d(dP)/d(Y[i])  [Pa]
  bool   mach_clamped{false};      // true if |Mach| was reduced to MACH_CLAMP
};

// Network-level wrapper for sharp-edge area change.
// Evaluates density and viscosity from (T_up, P_static_up, Y_up) internally.
//
// Inputs:
//   m_dot          : mass flow rate [kg/s], signed (positive = F0 -> F1)
//   P_total_up     : upstream total pressure [Pa] (unused, for API consistency)
//   P_static_up    : upstream static pressure [Pa] (used for rho, mu)
//   T_up           : upstream temperature [K]
//   Y_up           : upstream mass fractions [-]
//   P_static_down  : downstream static pressure [Pa] (unused)
//   F0             : area at node-0 side [m^2]
//   F1             : area at node-1 side [m^2]
//   m_scale        : sigmoid transition width [kg/s]
//   D_h            : hydraulic diameter [m] (default 0 = circular)
AreaChangeElementResult area_change_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double T_up,
    const std::vector<double>& Y_up, double P_static_down,
    double F0, double F1, double m_scale = 1e-4, double D_h = 0.0);

// Network-level wrapper for conical (gradual) area change.
// Same interface as area_change_residuals_and_jacobian with additional length.
//
// Inputs:
//   length         : axial length of conical section [m]
//   (all others identical to area_change_residuals_and_jacobian)
AreaChangeElementResult conical_area_change_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double T_up,
    const std::vector<double>& Y_up, double P_static_down,
    double F0, double F1, double length, double m_scale = 1e-4);

// -----------------------------------------------------------------------------
// 5. Tee Junction Elements (Bassett 2001)
// -----------------------------------------------------------------------------
// Three-port pressure-loss elements: MergingTee and BranchingTee.
// Both residuals are normalised by 0.5*rho*u_com^2.
//
// Port convention:
//   MergingTee:   MAIN_INLET=A, BRANCH=B (inflow),  MAIN_OUTLET=C (common)
//   BranchingTee: MAIN_OUTLET=A, BRANCH=B (outflow), MAIN_INLET=C (common)
//
// Inputs common to both wrappers:
//   m_dot_com    : total flow at common port C [kg/s] (positive = design direction)
//   m_dot_branch : branch flow at port B [kg/s] (positive = design direction)
//   dP0_straight : p0_upstream - p0_downstream for the straight path [Pa]
//   dP0_branch   : p0_upstream - p0_downstream for the branch path [Pa]
//   P_static_com : static pressure at common port [Pa] (used for density)
//   T_com        : temperature at common port [K]
//   Y_com        : mass fractions at common port [-]
//   theta        : branch angle [rad] (validated range: (0, pi/2])
//   psi          : area ratio F_C/F_B [-] (validated range: >= 0.05)
//   F_C          : main duct cross-section area [m^2]
//   blend_k      : tanh blend sharpness (default 30.0)

struct TeeJunctionResult {
    // Straight path residual and m_dot Jacobians (analytical)
    double R_straight;
    double dR_straight_d_mdot_com;
    double dR_straight_d_mdot_branch;
    // Straight path: density-dependent Jacobians (finite difference)
    double dR_straight_dP_static_com;
    double dR_straight_dT_com;
    std::vector<double> dR_straight_dY_com;

    // Branch path residual and m_dot Jacobians (analytical)
    double R_branch;
    double dR_branch_d_mdot_com;
    double dR_branch_d_mdot_branch;
    // Branch path: density-dependent Jacobians (finite difference)
    double dR_branch_dP_static_com;
    double dR_branch_dT_com;
    std::vector<double> dR_branch_dY_com;

    // Diagnostics
    double K_straight;      // effective straight-path loss coefficient
    double K_branch;        // effective branch-path loss coefficient
    double q;               // flow ratio m_dot_branch / m_dot_com
    double blend_w;         // tanh blend weight (1 = primary topology)
    bool topology_valid;    // true if q in design direction (within 5% margin)
    CorrelationValidity status; // VALID if inputs in validated Bassett range
};

TeeJunctionResult merging_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com,
    const std::vector<double>& Y_com,
    double theta, double psi, double F_C,
    double blend_k = 30.0);

TeeJunctionResult branching_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com,
    const std::vector<double>& Y_com,
    double theta, double psi, double F_C,
    double blend_k = 30.0);

// -----------------------------------------------------------------------------
// 6. Compressible Unified0D Tee Junction (Mynard & Valen-Sendstad 2015)
// -----------------------------------------------------------------------------
// Replaces the Bassett incompressible K tables with a geometry-derived K to
// O(M^2) that is consistent with mass-flow continuity and stagnation pressure.
// See docs/junction/junction_model.tex for derivation.
//
// Per-branch input struct.  mdot sign: >0 = supplier (inflow to junction).
// P_static and Pt are both passed explicitly (available as state_X.P / Pt in
// the Python solver; C++ does not need to invert the stagnation relation).
struct BranchInput {
    double P_static; // static pressure [Pa]
    double Pt;       // stagnation pressure [Pa]
    double T;        // static temperature [K]
    double m_dot;    // signed mass flow [kg/s]  (>0 = supplier)
    double A;        // cross-section area [m^2]
    double theta;    // centreline angle [rad]
    double gamma_eff;
    double R_gas;    // specific gas constant [J/(kg*K)]
};

// Result: two residuals [Pa] and their Jacobians w.r.t. all coupled unknowns.
// Naming: dR{0|1}_d{Pt|P|T|mdot}_{com|str|bra}
// Pt = stagnation pressure, P = static pressure, T = static temperature.
// Not all fields are nonzero for both tee types -- zero fields are explicitly set.
struct CompressibleTeeResult {
    double R_0;              // first residual [Pa]
    double R_1;              // second residual [Pa]
    // R_0 Jacobians
    double dR0_dPt_com;      // d(R_0)/d(Pt_com)
    double dR0_dP_com;       // d(R_0)/d(P_static_com)
    double dR0_dT_com;       // d(R_0)/d(T_com)
    double dR0_dmdot_com;    // d(R_0)/d(m_dot_com)
    double dR0_dPt_str;      // d(R_0)/d(Pt_str)
    double dR0_dP_str;       // d(R_0)/d(P_static_str)
    double dR0_dT_str;
    double dR0_dPt_bra;
    double dR0_dP_bra;       // d(R_0)/d(P_static_bra)
    double dR0_dT_bra;
    double dR0_dmdot_branch; // d(R_0)/d(m_dot_branch element unknown)
    // R_1 Jacobians
    double dR1_dPt_com;
    double dR1_dP_com;
    double dR1_dT_com;
    double dR1_dmdot_com;
    double dR1_dPt_str;
    double dR1_dP_str;
    double dR1_dT_str;
    double dR1_dPt_bra;
    double dR1_dP_bra;
    double dR1_dT_bra;
    double dR1_dmdot_branch;
};

// Branching tee: common is the single supplier; straight and branch are collectors.
//   R_0 = p0_dat - p0_straight - K_straight * q_ref
//   R_1 = p0_dat - p0_branch   - K_branch   * q_ref
// com.m_dot and bra.m_dot carry the element unknowns m_dot_com, m_dot_branch.
// str.m_dot = com.m_dot - bra.m_dot (derived; pass the current iterate value).
CompressibleTeeResult compressible_branching_tee_rj(
    const BranchInput& com,
    const BranchInput& str,
    const BranchInput& bra);

// Merging tee: straight and branch are suppliers; common is the collector.
//   R_0 = P_static_straight - P_static_branch   (supplier static-pressure equality)
//   R_1 = p0_dat - p0_com - K_com * q_ref
// All three m_dot values are passed (str.m_dot + bra.m_dot == com.m_dot by mass balance).
CompressibleTeeResult compressible_merging_tee_rj(
    const BranchInput& com,
    const BranchInput& str,
    const BranchInput& bra);

} // namespace solver
} // namespace combaero

#endif // SOLVER_INTERFACE_H
