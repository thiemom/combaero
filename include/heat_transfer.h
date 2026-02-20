#ifndef HEAT_TRANSFER_H
#define HEAT_TRANSFER_H

#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "correlation_status.h"

// -------------------------------------------------------------
// Heat Transfer Correlations for Forced Convection
// -------------------------------------------------------------
// All Nusselt correlations return Nu [-], from which:
//   h = Nu * k / L   [W/(m²·K)]
// where k is thermal conductivity [W/(m·K)] and L is characteristic length [m].
//
// For pipe flow: L = D (diameter)
// For flat plate: L = x (distance from leading edge) or plate length
//
// References:
// - Dittus-Boelter (1930): Univ. California Publ. Eng., 2, 443
// - Gnielinski (1976): Int. Chem. Eng., 16, 359
// - Sieder-Tate (1936): Ind. Eng. Chem., 28, 1429
// - Petukhov (1970): Advances in Heat Transfer, 6, 503

// -------------------------------------------------------------
// Internal Flow (Pipe/Duct)
// -------------------------------------------------------------

// Dittus-Boelter correlation (1930)
// Nu = 0.023 * Re^0.8 * Pr^n
// where n = 0.4 for heating (fluid being heated)
//       n = 0.3 for cooling (fluid being cooled)
//
// Valid for:
//   - Fully developed turbulent flow
//   - Re > 10,000
//   - 0.6 < Pr < 160
//   - L/D > 10 (entrance effects negligible)
//
// Parameters:
//   Re      : Reynolds number [-]
//   Pr      : Prandtl number [-]
//   heating : true if fluid is being heated, false if cooled
//   status  : if non-null, set to Extrapolated when Re or Pr is outside the
//             validated range; no warning is emitted in that case
double nusselt_dittus_boelter(double Re, double Pr, bool heating = true,
                              combaero::CorrelationStatus* status = nullptr);

// Gnielinski correlation (1976)
// Nu = (f/8) * (Re - 1000) * Pr / (1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1))
//
// More accurate than Dittus-Boelter, especially in transition region.
// Valid for:
//   - 2300 < Re < 5x10^6
//   - 0.5 < Pr < 2000
//
// Parameters:
//   Re     : Reynolds number [-]
//   Pr     : Prandtl number [-]
//   f      : Darcy friction factor [-] (use friction_colebrook or similar)
//   status : if non-null, set to Extrapolated when Re or Pr is outside the
//            validated range; no warning is emitted in that case.
//            Below Re=2300 a C1-smooth Hermite blend to laminar Nu is used
//            (the raw formula goes negative at Re<1000).
double nusselt_gnielinski(double Re, double Pr, double f,
                          combaero::CorrelationStatus* status = nullptr);

// Gnielinski with automatic friction factor (smooth pipe)
// Uses Petukhov friction correlation: f = (0.790*ln(Re) - 1.64)^(-2)
double nusselt_gnielinski(double Re, double Pr,
                          combaero::CorrelationStatus* status = nullptr);

// Sieder-Tate correlation (1936)
// Nu = 0.027 * Re^0.8 * Pr^(1/3) * (μ_bulk / μ_wall)^0.14
//
// Accounts for viscosity variation between bulk and wall temperatures.
// Valid for:
//   - Re > 10,000
//   - 0.7 < Pr < 16,700
//   - L/D > 10
//
// Parameters:
//   Re       : Reynolds number at bulk temperature [-]
//   Pr       : Prandtl number at bulk temperature [-]
//   mu_ratio : μ_bulk / μ_wall [-] (typically 0.5-2.0)
//   status   : if non-null, set to Extrapolated when Re or Pr is outside the
//              validated range; no warning is emitted in that case
double nusselt_sieder_tate(double Re, double Pr, double mu_ratio = 1.0,
                           combaero::CorrelationStatus* status = nullptr);

// Petukhov correlation (1970)
// Nu = (f/8) * Re * Pr / (1.07 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1))
//
// Basis for Gnielinski; valid for fully turbulent flow.
// Valid for:
//   - 10^4 < Re < 5x10^6
//   - 0.5 < Pr < 2000
//
// Parameters:
//   Re     : Reynolds number [-]
//   Pr     : Prandtl number [-]
//   f      : Darcy friction factor [-]
//   status : if non-null, set to Extrapolated when Re or Pr is outside the
//            validated range; no warning is emitted in that case
double nusselt_petukhov(double Re, double Pr, double f,
                        combaero::CorrelationStatus* status = nullptr);

// Petukhov with automatic friction factor (smooth pipe)
double nusselt_petukhov(double Re, double Pr,
                        combaero::CorrelationStatus* status = nullptr);

// -------------------------------------------------------------
// Laminar Flow (for completeness)
// -------------------------------------------------------------

// Fully developed laminar flow in circular pipe
// Nu = 3.66 (constant wall temperature)
// Nu = 4.36 (constant heat flux)
constexpr double NU_LAMINAR_CONST_T = 3.66;
constexpr double NU_LAMINAR_CONST_Q = 4.36;

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

// Heat transfer coefficient from Nusselt number
// h = Nu * k / L  [W/(m²·K)]
//
// Parameters:
//   Nu : Nusselt number [-]
//   k  : thermal conductivity [W/(m·K)]
//   L  : characteristic length [m] (diameter for pipe flow)
double htc_from_nusselt(double Nu, double k, double L);

// -------------------------------------------------------------
// Overall Heat Transfer Coefficient (Thermal Resistance Networks)
// -------------------------------------------------------------

// Overall HTC from individual resistances (per unit area)
// 1/U = 1/h1 + 1/h2 + ... + t1/k1 + t2/k2 + ...
//
// For a typical tube heat exchanger:
//   1/U = 1/h_inner + t_wall/k_wall + 1/h_outer + R_fouling
//
// Parameters:
//   h_values : vector of convective HTCs [W/(m²·K)]
//   t_over_k : vector of thickness/conductivity ratios [m / W/(m·K) = m²·K/W]
// Returns: overall HTC [W/(m²·K)]
double overall_htc(const std::vector<double>& h_values,
                   const std::vector<double>& t_over_k = {});

// Convenience: overall HTC for wall with inner/outer convection and single layer
// 1/U = 1/h_inner + t_wall/k_wall + 1/h_outer
double overall_htc_wall(double h_inner, double h_outer,
                        double t_wall, double k_wall);

// Multi-layer wall: h_in, [t1/k1, t2/k2, ...], h_out
// 1/U = 1/h_inner + Σ(t_i/k_i) + 1/h_outer
//
// Example: insulated pipe with steel + insulation
//   overall_htc_wall(500, 10, {0.003/50, 0.05/0.04})
//   = h_in=500, h_out=10, steel 3mm @ k=50, insulation 50mm @ k=0.04
double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double>& t_over_k_layers);

// With additional fouling resistance (R_f in m²·K/W)
double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double>& t_over_k_layers,
                        double R_fouling);

// Legacy alias for single-layer wall
inline double overall_htc_tube(double h_inner, double h_outer,
                               double t_wall, double k_wall) {
    return overall_htc_wall(h_inner, h_outer, t_wall, k_wall);
}

inline double overall_htc_tube(double h_inner, double h_outer,
                               double t_wall, double k_wall,
                               double R_fouling) {
    std::vector<double> layers = {t_wall / k_wall};
    return overall_htc_wall(h_inner, h_outer, layers, R_fouling);
}

// Thermal resistance from HTC and area
// R = 1 / (h * A)  [K/W]
double thermal_resistance(double h, double A);

// Thermal resistance for conduction through wall
// R = t / (k * A)  [K/W]
double thermal_resistance_wall(double thickness, double k, double A);

// -------------------------------------------------------------
// Log Mean Temperature Difference (LMTD)
// -------------------------------------------------------------

// LMTD for heat exchangers
// LMTD = (ΔT1 - ΔT2) / ln(ΔT1 / ΔT2)
//
// For counter-flow:
//   ΔT1 = T_hot_in - T_cold_out
//   ΔT2 = T_hot_out - T_cold_in
//
// For parallel-flow:
//   ΔT1 = T_hot_in - T_cold_in
//   ΔT2 = T_hot_out - T_cold_out
//
// Parameters:
//   dT1 : temperature difference at one end [K]
//   dT2 : temperature difference at other end [K]
// Returns: log mean temperature difference [K]
//
// Note: If dT1 ≈ dT2, returns arithmetic mean to avoid 0/0.
double lmtd(double dT1, double dT2);

// LMTD for counter-flow heat exchanger (convenience function)
// Parameters:
//   T_hot_in   : hot fluid inlet temperature [K]
//   T_hot_out  : hot fluid outlet temperature [K]
//   T_cold_in  : cold fluid inlet temperature [K]
//   T_cold_out : cold fluid outlet temperature [K]
double lmtd_counterflow(double T_hot_in, double T_hot_out,
                        double T_cold_in, double T_cold_out);

// LMTD for parallel-flow heat exchanger (convenience function)
double lmtd_parallelflow(double T_hot_in, double T_hot_out,
                         double T_cold_in, double T_cold_out);

// -------------------------------------------------------------
// Heat Transfer Rate and Heat Flux
// -------------------------------------------------------------

// Heat transfer rate Q = U * A * ΔT  [W]
// U  : overall heat transfer coefficient [W/(m²·K)]
// A  : heat transfer area [m²]
// dT : temperature difference [K]
inline double heat_rate(double U, double A, double dT) {
    return U * A * dT;
}

// Heat flux q = U * ΔT  [W/m²]
// (heat rate per unit area)
inline double heat_flux(double U, double dT) {
    return U * dT;
}

// Required area for given heat rate
// A = Q / (U * ΔT)  [m²]
inline double heat_transfer_area(double Q, double U, double dT) {
    return Q / (U * dT);
}

// Temperature difference for given heat rate
// ΔT = Q / (U * A)  [K]
inline double heat_transfer_dT(double Q, double U, double A) {
    return Q / (U * A);
}

// -------------------------------------------------------------
// Temperature Profile Through Wall
// -------------------------------------------------------------

// Temperature profile through a multi-layer wall
// Returns N+1 wall-surface temperatures (not bulk fluid temperatures):
//   [T_wall_hot_surface, T_layer1_2, T_layer2_3, ..., T_wall_cold_surface]
// where N = number of conductive layers.
//
// Note: T_wall_hot_surface = T_hot - q/h_hot  (hot wall surface, inside boundary layer)
//       T_wall_cold_surface = last element     (cold wall surface, inside boundary layer)
//       To recover cold bulk: T_cold_bulk = T_wall_cold_surface - q/h_cold
//
// Parameters:
//   T_hot      : hot-side bulk fluid temperature [K]
//   T_cold     : cold-side bulk fluid temperature [K]
//   h_hot      : hot-side convective HTC [W/(m²·K)]
//   h_cold     : cold-side convective HTC [W/(m²·K)]
//   t_over_k   : vector of thickness/conductivity for each layer [m²·K/W]
//   q          : (output) heat flux through wall [W/m²]
//
// The heat flux is constant through all layers (steady state).
std::vector<double> wall_temperature_profile(
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k,
    double& q);

// Simplified version without heat flux output
std::vector<double> wall_temperature_profile(
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k);

// Single-layer wall convenience function
// Returns: {T_hot_surface, T_cold_surface}
inline std::vector<double> wall_temperature_profile(
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    double t_wall, double k_wall) {
    return wall_temperature_profile(T_hot, T_cold, h_hot, h_cold,
                                    std::vector<double>{t_wall / k_wall});
}

// -------------------------------------------------------------
// Heat Flux from Measured Temperature
// -------------------------------------------------------------
//
// Solve for heat flux given a measured temperature at a known location.
// Useful for inferring heat transfer from thermocouple measurements.
//
// Edge indexing (for N conductive layers):
//   Edge 0: hot-side wall surface (after convective boundary layer)
//   Edge 1: interface between layer 0 and layer 1
//   Edge 2: interface between layer 1 and layer 2
//   ...
//   Edge N: cold-side wall surface (before convective boundary layer)
//
// Example: 2-layer wall (steel + insulation)
//   Edge 0: steel hot surface
//   Edge 1: steel-insulation interface
//   Edge 2: insulation cold surface

// Heat flux from temperature measured at a wall edge
// Parameters:
//   T_measured : measured temperature at the edge [K]
//   edge_idx   : edge index (0 = hot surface, N = cold surface)
//   T_hot      : hot-side bulk fluid temperature [K]
//   T_cold     : cold-side bulk fluid temperature [K]
//   h_hot      : hot-side convective HTC [W/(m²·K)]
//   h_cold     : cold-side convective HTC [W/(m²·K)]
//   t_over_k   : vector of thickness/conductivity for each layer [m²·K/W]
//
// Returns: heat flux q [W/m²]
double heat_flux_from_T_at_edge(
    double T_measured, std::size_t edge_idx,
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k);

// Heat flux from temperature measured at a depth from hot surface
// Parameters:
//   T_measured     : measured temperature [K]
//   depth_from_hot : depth from hot surface into wall [m]
//   T_hot, T_cold  : bulk fluid temperatures [K]
//   h_hot, h_cold  : convective HTCs [W/(m²·K)]
//   thicknesses    : vector of layer thicknesses [m]
//   conductivities : vector of layer thermal conductivities [W/(m·K)]
//
// Returns: heat flux q [W/m²]
// Throws if depth is outside the wall or negative
double heat_flux_from_T_at_depth(
    double T_measured, double depth_from_hot,
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& thicknesses,
    const std::vector<double>& conductivities);

// Infer unknown bulk temperature from measured edge temperature and heat flux
// Useful when you know q (e.g., from power input) and measure a surface T
//
// Parameters:
//   T_measured : measured temperature at edge [K]
//   edge_idx   : edge index
//   q          : known heat flux [W/m²] (positive = hot to cold)
//   h_hot      : hot-side HTC [W/(m²·K)]
//   h_cold     : cold-side HTC [W/(m²·K)]
//   t_over_k   : layer resistances [m²·K/W]
//   solve_for  : "hot" or "cold" - which bulk temperature to solve for
//
// Returns: inferred bulk fluid temperature [K]
double bulk_T_from_edge_T_and_q(
    double T_measured, std::size_t edge_idx, double q,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k,
    const std::string& solve_for);

// -------------------------------------------------------------
// Temperature Sensitivity (for risk assessment)
// -------------------------------------------------------------
//
// Partial derivatives of edge temperatures with respect to bulk temperatures.
// Useful for:
//   - Risk assessment: "If T_hot increases by 50K, how much does the
//     steel-ceramic interface temperature increase?"
//   - Uncertainty propagation
//   - Control system design
//
// For steady-state 1D conduction:
//   T_edge = T_hot - q * R_hot_to_edge
//          = T_hot - (T_hot - T_cold) / R_total * R_hot_to_edge
//          = T_hot * (1 - R_hot_to_edge/R_total) + T_cold * (R_hot_to_edge/R_total)
//
// Therefore:
//   ∂T_edge/∂T_hot  = R_edge_to_cold / R_total  (always positive, ≤1)
//   ∂T_edge/∂T_cold = R_hot_to_edge / R_total   (always positive, ≤1)
//   Sum = 1 (as expected for linear interpolation)

// Sensitivity of edge temperature to hot-side bulk temperature change
// Returns: ∂T_edge/∂T_hot [-] (dimensionless, 0 to 1)
double dT_edge_dT_hot(
    std::size_t edge_idx,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k);

// Sensitivity of edge temperature to cold-side bulk temperature change
// Returns: ∂T_edge/∂T_cold [-] (dimensionless, 0 to 1)
double dT_edge_dT_cold(
    std::size_t edge_idx,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k);

// Both sensitivities at once (more efficient)
// Returns: {∂T_edge/∂T_hot, ∂T_edge/∂T_cold}
std::pair<double, double> dT_edge_dT_bulk(
    std::size_t edge_idx,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k);

// Sensitivity of edge temperature to heat flux change
// Returns: ∂T_edge/∂q [K·m²/W] (negative for edges closer to hot side)
// Useful for: "If heat flux increases by 100 W/m², how does edge T change?"
double dT_edge_dq(
    std::size_t edge_idx,
    double h_hot,
    const std::vector<double>& t_over_k);

// -------------------------------------------------------------
// NTU-Effectiveness Method
// -------------------------------------------------------------

// Number of Transfer Units
// NTU = U * A / C_min
// C_min : minimum heat capacity rate [W/K]
inline double ntu(double U, double A, double C_min) {
    return U * A / C_min;
}

// Heat capacity rate ratio
// C_r = C_min / C_max
inline double capacity_ratio(double C_min, double C_max) {
    return C_min / C_max;
}

// Effectiveness for counter-flow heat exchanger
// ε = (1 - exp(-NTU*(1-C_r))) / (1 - C_r*exp(-NTU*(1-C_r)))
// For C_r = 1: ε = NTU / (1 + NTU)
double effectiveness_counterflow(double NTU, double C_r);

// Effectiveness for parallel-flow heat exchanger
// ε = (1 - exp(-NTU*(1+C_r))) / (1 + C_r)
double effectiveness_parallelflow(double NTU, double C_r);

// Actual heat transfer from effectiveness
// Q = ε * C_min * (T_hot_in - T_cold_in)
inline double heat_rate_from_effectiveness(double epsilon, double C_min,
                                           double T_hot_in, double T_cold_in) {
    return epsilon * C_min * (T_hot_in - T_cold_in);
}

// -------------------------------------------------------------
// State-based convenience functions
// -------------------------------------------------------------

struct State;  // Forward declaration

// Nusselt number for pipe flow using State
// Automatically computes Re, Pr from state and velocity/diameter
double nusselt_pipe(const State& s, double velocity, double diameter,
                    bool heating = true, double roughness = 0.0);

// Heat transfer coefficient for pipe flow [W/(m²·K)]
double htc_pipe(const State& s, double velocity, double diameter,
                bool heating = true, double roughness = 0.0);

// Composite function: compute HTC, Nu, and Re from thermodynamic state
// Returns: tuple (h [W/(m²·K)], Nu [-], Re [-])
//
// Parameters:
//   T           : bulk static temperature [K]
//                 Re and Pr are always evaluated at static T (correct).
//                 For the heat flux driving temperature:
//                   M < 0.3 : pass T_static as T_hot to cooled_wall_heat_flux()
//                   M > 0.3 : pass T_adiabatic_wall() from stagnation.h instead
//                             (T_aw is between T_static and T_total; using T_total
//                              overcorrects since recovery factor r < 1)
//   P           : pressure [Pa]
//   X           : mole fractions [mol/mol]
//   velocity    : flow velocity [m/s]
//   diameter    : pipe diameter [m]
//   correlation : "gnielinski" (default), "dittus_boelter", "sieder_tate", "petukhov"
//   heating     : true for heating, false for cooling (affects Dittus-Boelter)
//   mu_ratio    : μ_bulk / μ_wall for Sieder-Tate viscosity correction (default: 1.0)
//   roughness   : absolute roughness [m] (default: 0.0 = smooth pipe)
//
// Automatically computes: ρ, μ, k, Pr, Re from (T, P, X)
// Selects appropriate correlation and handles laminar flow (Re < 2300)
std::tuple<double, double, double> htc_pipe(
    double T, double P, const std::vector<double>& X,
    double velocity, double diameter,
    const std::string& correlation = "gnielinski",
    bool heating = true,
    double mu_ratio = 1.0,
    double roughness = 0.0);

// -------------------------------------------------------------
// Combined convective heat transfer + pressure loss result
// -------------------------------------------------------------
//
// Returned by all channel_* functions below.
// All properties are evaluated from (T, P, X) in a single pass.
//
// T_aw is always computed continuously for all M:
//   T_aw = T_static + r * v² / (2 * cp)
//   r = Pr^(1/3) turbulent, Pr^(1/2) laminar
//   At v=0: T_aw = T_static exactly (no discontinuity, no threshold).
//   This keeps the Jacobian smooth for network solvers.
//
// q is nan when T_wall is not supplied (flow-only or unknown wall T).
//
// f meaning depends on geometry:
//   smooth / ribbed / dimpled : Darcy friction factor [-]
//   pin_fin                   : pin-array friction coefficient [-]
//                               (dP = N_rows * f * rho*v_max²/2)
//   impingement               : jet-plate loss coefficient [-]
//                               (dP = f * rho*v_jet²/2)
struct ChannelResult {
    double h;     // convective HTC [W/(m²·K)]
    double Nu;    // Nusselt number [-]
    double Re;    // Reynolds number [-]  (hydraulic D or pin D depending on geometry)
    double Pr;    // Prandtl number [-]
    double f;     // friction / loss coefficient [-]  (see note above)
    double dP;    // pressure drop [Pa]
    double M;     // Mach number [-]  (computed internally from v/a(T,X))
    double T_aw;  // adiabatic wall temperature [K]  (always computed, see above)
    double q;     // heat flux [W/m²]  (nan if T_wall not supplied)
};

// -------------------------------------------------------------
// High-level channel functions — each returns a ChannelResult
// -------------------------------------------------------------
//
// All functions:
//   - evaluate fluid properties (rho, mu, k, Pr, cp) once from (T, P, X)
//   - compute M = v / a(T, X) internally
//   - compute T_aw = T_static + r*v²/(2*cp) continuously
//   - set q = h*(T_aw - T_wall) if T_wall is finite, else q = nan
//
// User contract: choose the model that matches the geometry.
// Do not mix models (e.g. do not apply rib multipliers to a smooth result).

// Smooth pipe or duct (Gnielinski / Dittus-Boelter / Sieder-Tate / Petukhov)
//
// Parameters:
//   T, P, X    : bulk static thermodynamic state
//   velocity   : bulk velocity [m/s]
//   diameter   : hydraulic diameter [m]
//   length     : channel length [m]  (for dP)
//   T_wall     : wall temperature [K]  (pass NaN to skip q computation)
//   correlation: "gnielinski" (default), "dittus_boelter", "sieder_tate", "petukhov"
//   heating    : true = fluid is heated (affects Dittus-Boelter exponent)
//   mu_ratio   : mu_bulk / mu_wall for Sieder-Tate (default 1.0)
//   roughness  : absolute wall roughness [m]  (default 0.0 = smooth)
ChannelResult channel_smooth(
    double T, double P, const std::vector<double>& X,
    double velocity, double diameter, double length,
    double T_wall = std::numeric_limits<double>::quiet_NaN(),
    const std::string& correlation = "gnielinski",
    bool heating = true,
    double mu_ratio = 1.0,
    double roughness = 0.0);

// Rib-enhanced cooling channel (Han et al. 1988)
//
// Applies rib_enhancement_factor to Nu and rib_friction_multiplier to f
// from a smooth-pipe baseline at the same Re.
//
// Parameters:
//   e_D   : rib height / hydraulic diameter [-]  (valid: 0.02-0.1)
//   P_e   : rib pitch / rib height [-]           (valid: 5-20)
//   alpha : rib angle [degrees]                  (valid: 30-90)
ChannelResult channel_ribbed(
    double T, double P, const std::vector<double>& X,
    double velocity, double diameter, double length,
    double e_D, double P_e, double alpha,
    double T_wall = std::numeric_limits<double>::quiet_NaN(),
    bool heating = true);

// Dimpled surface cooling channel (Chyu et al. 1997)
//
// Applies dimple_nusselt_enhancement to Nu and dimple_friction_multiplier to f
// from a smooth-pipe baseline at the same Re.
//
// Parameters:
//   d_Dh : dimple diameter / channel height [-]  (valid: 0.1-0.3)
//   h_d  : dimple depth / diameter [-]           (valid: 0.1-0.3)
//   S_d  : dimple spacing / diameter [-]         (valid: 1.5-3.0)
ChannelResult channel_dimpled(
    double T, double P, const std::vector<double>& X,
    double velocity, double diameter, double length,
    double d_Dh, double h_d, double S_d,
    double T_wall = std::numeric_limits<double>::quiet_NaN(),
    bool heating = true);

// Pin-fin array cooling channel (Metzger et al. 1982)
//
// Re is based on pin diameter d.
// h = Nu * k / d.
// dP = N_rows * f_pin * (rho * v_max² / 2)
//   where v_max = velocity * S_D / (S_D - 1)  (minimum cross-section velocity)
//   and   f_pin from pin_fin_friction() in cooling_correlations.h
//
// Parameters:
//   velocity       : approach (upstream) velocity [m/s]
//   channel_height : pin length H [m]
//   pin_diameter   : pin diameter d [m]
//   S_D            : spanwise pitch / d [-]   (valid: 1.5-4.0)
//   X_D            : streamwise pitch / d [-] (valid: 1.5-4.0)
//   N_rows         : number of pin rows in streamwise direction
//   is_staggered   : true = staggered array (default), false = inline
ChannelResult channel_pin_fin(
    double T, double P, const std::vector<double>& X,
    double velocity, double channel_height, double pin_diameter,
    double S_D, double X_D, int N_rows,
    double T_wall = std::numeric_limits<double>::quiet_NaN(),
    bool is_staggered = true);

// Impingement jet array cooling (Florschuetz et al. 1981 / Martin 1977)
//
// Re is based on jet diameter d_jet.
// h = Nu * k / d_jet.
// dP = f * rho * v_jet² / 2  where f = 1/Cd_jet²  (jet-plate orifice loss)
//   v_jet = mdot_jet / (rho * A_jet)
//
// Parameters:
//   mdot_jet : jet mass flow rate [kg/s]  (total for all jets)
//   d_jet    : jet hole diameter [m]
//   z_D      : jet-to-target distance / d_jet [-]  (valid: 1-12)
//   x_D      : streamwise jet spacing / d_jet [-]  (0 = single jet; valid array: 4-16)
//   y_D      : spanwise jet spacing / d_jet [-]    (0 = single jet; valid array: 4-16)
//   A_target : target surface area [m²]  (used to compute average h)
//   Cd_jet   : jet hole discharge coefficient [-]  (default 0.65)
ChannelResult channel_impingement(
    double T, double P, const std::vector<double>& X,
    double mdot_jet, double d_jet,
    double z_D, double x_D, double y_D,
    double A_target,
    double T_wall = std::numeric_limits<double>::quiet_NaN(),
    double Cd_jet = 0.65);

#endif // HEAT_TRANSFER_H
