#ifndef HEAT_TRANSFER_H
#define HEAT_TRANSFER_H

#include <string>
#include <utility>
#include <vector>

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
double nusselt_dittus_boelter(double Re, double Pr, bool heating = true);

// Gnielinski correlation (1976)
// Nu = (f/8) * (Re - 1000) * Pr / (1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1))
//
// More accurate than Dittus-Boelter, especially in transition region.
// Valid for:
//   - 2300 < Re < 5×10^6
//   - 0.5 < Pr < 2000
//
// Parameters:
//   Re  : Reynolds number [-]
//   Pr  : Prandtl number [-]
//   f   : Darcy friction factor [-] (use friction_colebrook or similar)
double nusselt_gnielinski(double Re, double Pr, double f);

// Gnielinski with automatic friction factor (smooth pipe)
// Uses Petukhov friction correlation: f = (0.790*ln(Re) - 1.64)^(-2)
double nusselt_gnielinski(double Re, double Pr);

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
double nusselt_sieder_tate(double Re, double Pr, double mu_ratio = 1.0);

// Petukhov correlation (1970)
// Nu = (f/8) * Re * Pr / (1.07 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1))
//
// Basis for Gnielinski; valid for fully turbulent flow.
// Valid for:
//   - 10^4 < Re < 5×10^6
//   - 0.5 < Pr < 2000
//
// Parameters:
//   Re  : Reynolds number [-]
//   Pr  : Prandtl number [-]
//   f   : Darcy friction factor [-]
double nusselt_petukhov(double Re, double Pr, double f);

// Petukhov with automatic friction factor (smooth pipe)
double nusselt_petukhov(double Re, double Pr);

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
// Returns vector of temperatures at each interface:
//   [T_hot_surface, T_layer1_2, T_layer2_3, ..., T_cold_surface]
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

#endif // HEAT_TRANSFER_H
