#ifndef HEAT_TRANSFER_H
#define HEAT_TRANSFER_H

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
