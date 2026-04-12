#include "../include/heat_transfer.h"
#include "../include/cooling_correlations.h" // rib_*, dimple_*, pin_fin_*, impingement_*
#include "../include/correlation_status.h"
#include "../include/friction.h"       // friction_petukhov, friction_colebrook
#include "../include/math_constants.h" // M_PI cross-platform
#include "../include/solver_interface.h" // density_and_jacobians, viscosity_and_jacobians
#include "../include/stagnation.h"     // T_adiabatic_wall
#include "../include/state.h"
#include "../include/thermo.h"    // density, cp_mass, speed_of_sound
#include "../include/transport.h" // viscosity, thermal_conductivity, prandtl
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>

namespace combaero {

// -------------------------------------------------------------
// Internal Flow Correlations
// -------------------------------------------------------------

double nusselt_dittus_boelter(double Re, double Pr, bool heating,
                              combaero::CorrelationStatus *status) {
  bool extrapolated = false;
  if (Re < 10000.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_dittus_boelter: Re = " + std::to_string(Re) +
                     " is below validated range [10000, inf]."
                     " Extrapolating power-law correlation; check results.");
    }
  }
  if (Pr < 0.6 || Pr > 160.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_dittus_boelter: Pr = " + std::to_string(Pr) +
                     " is outside validated range [0.6, 160]."
                     " Extrapolating; check results.");
    }
  }
  if (status) {
    *status = extrapolated ? combaero::CorrelationStatus::Extrapolated
                           : combaero::CorrelationStatus::Valid;
  }
  double n = heating ? dittus_boelter::coeff_prandtl_heating_exp
                     : dittus_boelter::coeff_prandtl_cooling_exp;
  return dittus_boelter::coeff_outer *
         std::pow(Re, dittus_boelter::coeff_reynolds_exp) * std::pow(Pr, n);
}

double nusselt_gnielinski(double Re, double Pr, double f,
                          combaero::CorrelationStatus *status) {
  if (f <= 0.0) {
    throw std::invalid_argument(
        "nusselt_gnielinski: friction factor must be positive");
  }

  bool extrapolated = false;

  if (Re > 5.0e6) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_gnielinski: Re = " + std::to_string(Re) +
                     " is above validated range [2300, 5e6]."
                     " Extrapolating; check results.");
    }
  }
  if (Pr < 0.5 || Pr > 2000.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_gnielinski: Pr = " + std::to_string(Pr) +
                     " is outside validated range [0.5, 2000]."
                     " Extrapolating; check results.");
    }
  }
  if (status) {
    *status = extrapolated ? combaero::CorrelationStatus::Extrapolated
                           : combaero::CorrelationStatus::Valid;
  }

  double f8 = f / 8.0;
  double sqrt_f8 = std::sqrt(f8);
  double Pr_23 = std::pow(Pr, 2.0 / 3.0);
  double denom = 1.0 + gnielinski::coeff_prandtl_factor * sqrt_f8 * (Pr_23 - 1.0);

  if (Re >= 2300.0) {
    return f8 * (Re - gnielinski::coeff_re_offset) * Pr / denom;
  }

  // Below Re=2300 the raw Gnielinski formula goes negative at Re<1000.
  // Use a C1-smooth cubic Hermite blend between Re=1000 (laminar Nu) and
  // Re=2300 (Gnielinski value + slope), so the Jacobian is continuous.
  if (status && *status == combaero::CorrelationStatus::Valid) {
    *status = combaero::CorrelationStatus::Extrapolated;
  }
  if (!status) {
    combaero::warn("nusselt_gnielinski: Re = " + std::to_string(Re) +
                   " is below validated range [2300, 5e6]."
                   " Using C1 Hermite blend to laminar Nu; check results.");
  }

  // Anchor points for the Hermite cubic
  constexpr double Re0 = 1000.0;         // lower anchor: laminar Nu
  constexpr double Re1 = 2300.0;         // upper anchor: Gnielinski
  const double Nu0 = NU_LAMINAR_CONST_T; // 3.66
  const double dNu0 = 0.0;               // laminar: flat w.r.t. Re
  const double Nu1 = f8 * (Re1 - gnielinski::coeff_re_offset) * Pr / denom;
  const double dNu1 = f8 * Pr / denom; // d(Nu)/d(Re) at Re1

  // Normalised coordinate t in [0, 1]
  double h = Re1 - Re0;
  double t = (Re - Re0) / h;
  t = std::max(0.0, std::min(1.0, t)); // clamp outside [Re0, Re1]

  // Cubic Hermite basis
  double t2 = t * t;
  double t3 = t2 * t;
  double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
  double h10 = t3 - 2.0 * t2 + t;
  double h01 = -2.0 * t3 + 3.0 * t2;
  double h11 = t3 - t2;

  return h00 * Nu0 + h10 * h * dNu0 + h01 * Nu1 + h11 * h * dNu1;
}

double nusselt_gnielinski(double Re, double Pr,
                          combaero::CorrelationStatus *status) {
  double f =
      friction_petukhov(std::max(Re, 3000.0)); // avoid log(0) in Petukhov
  return nusselt_gnielinski(Re, Pr, f, status);
}

double nusselt_sieder_tate(double Re, double Pr, double mu_ratio,
                           combaero::CorrelationStatus *status) {
  if (mu_ratio <= 0.0) {
    throw std::invalid_argument(
        "nusselt_sieder_tate: mu_ratio must be positive");
  }
  bool extrapolated = false;
  if (Re < 10000.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_sieder_tate: Re = " + std::to_string(Re) +
                     " is below validated range [10000, inf]."
                     " Extrapolating power-law correlation; check results.");
    }
  }
  if (Pr < 0.7 || Pr > 16700.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_sieder_tate: Pr = " + std::to_string(Pr) +
                     " is outside validated range [0.7, 16700]."
                     " Extrapolating; check results.");
    }
  }
  if (status) {
    *status = extrapolated ? combaero::CorrelationStatus::Extrapolated
                           : combaero::CorrelationStatus::Valid;
  }
  return sieder_tate::coeff_outer
         * std::pow(Re, sieder_tate::coeff_reynolds_exp)
         * std::pow(Pr, sieder_tate::coeff_prandtl_exp)
         * std::pow(mu_ratio, sieder_tate::coeff_viscosity_exp);
}

double nusselt_petukhov(double Re, double Pr, double f,
                        combaero::CorrelationStatus *status) {
  if (f <= 0.0) {
    throw std::invalid_argument(
        "nusselt_petukhov: friction factor must be positive");
  }
  bool extrapolated = false;
  if (Re < 1.0e4 || Re > 5.0e6) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_petukhov: Re = " + std::to_string(Re) +
                     " is outside validated range [1e4, 5e6]."
                     " Extrapolating; check results.");
    }
  }
  if (Pr < 0.5 || Pr > 2000.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_petukhov: Pr = " + std::to_string(Pr) +
                     " is outside validated range [0.5, 2000]."
                     " Extrapolating; check results.");
    }
  }
  if (status) {
    *status = extrapolated ? combaero::CorrelationStatus::Extrapolated
                           : combaero::CorrelationStatus::Valid;
  }
  double f8 = f / 8.0;
  double sqrt_f8 = std::sqrt(f8);
  double Pr_23 = std::pow(Pr, 2.0 / 3.0);
  return f8 * Re * Pr / (petukhov_ht::coeff_offset + petukhov_ht::coeff_prandtl_factor * sqrt_f8 * (Pr_23 - 1.0));
}

double nusselt_petukhov(double Re, double Pr,
                        combaero::CorrelationStatus *status) {
  double f =
      friction_petukhov(std::max(Re, 3000.0)); // avoid log(0) in Petukhov
  return nusselt_petukhov(Re, Pr, f, status);
}

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

double htc_from_nusselt(double Nu, double k, double L) {
  if (L <= 0) {
    throw std::invalid_argument(
        "htc_from_nusselt: characteristic length must be positive");
  }
  if (k <= 0) {
    throw std::invalid_argument(
        "htc_from_nusselt: thermal conductivity must be positive");
  }
  return Nu * k / L;
}

// -------------------------------------------------------------
// Overall Heat Transfer Coefficient
// -------------------------------------------------------------

double overall_htc(const std::vector<double> &h_values,
                   const std::vector<double> &t_over_k) {
  if (h_values.empty()) {
    throw std::invalid_argument("overall_htc: at least one h value required");
  }

  double R_total = 0.0;

  // Add convective resistances: 1/h
  for (double h : h_values) {
    if (h <= 0) {
      throw std::invalid_argument("overall_htc: h values must be positive");
    }
    R_total += 1.0 / h;
  }

  // Add conductive resistances: t/k
  for (double tk : t_over_k) {
    if (tk < 0) {
      throw std::invalid_argument(
          "overall_htc: t/k values must be non-negative");
    }
    R_total += tk;
  }

  return 1.0 / R_total;
}

double overall_htc_wall(double h_inner, double h_outer, double t_wall,
                        double k_wall) {
  if (k_wall <= 0) {
    throw std::invalid_argument("overall_htc_wall: k_wall must be positive");
  }
  return overall_htc({h_inner, h_outer}, {t_wall / k_wall});
}

double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double> &t_over_k_layers) {
  return overall_htc({h_inner, h_outer}, t_over_k_layers);
}

double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double> &t_over_k_layers,
                        double R_fouling) {
  if (R_fouling < 0) {
    throw std::invalid_argument(
        "overall_htc_wall: R_fouling must be non-negative");
  }
  // Append fouling resistance to layers
  std::vector<double> all_resistances = t_over_k_layers;
  all_resistances.push_back(R_fouling);
  return overall_htc({h_inner, h_outer}, all_resistances);
}

double thermal_resistance(double h, double A) {
  if (h <= 0 || A <= 0) {
    throw std::invalid_argument("thermal_resistance: h and A must be positive");
  }
  return 1.0 / (h * A);
}

double thermal_resistance_wall(double thickness, double k, double A) {
  if (thickness < 0 || k <= 0 || A <= 0) {
    throw std::invalid_argument("thermal_resistance_wall: invalid parameters");
  }
  return thickness / (k * A);
}

// -------------------------------------------------------------
// Log Mean Temperature Difference (LMTD)
// -------------------------------------------------------------

double lmtd(double dT1, double dT2) {
  if (dT1 <= 0 || dT2 <= 0) {
    throw std::invalid_argument(
        "lmtd: temperature differences must be positive");
  }

  // If dT1 ≈ dT2, use arithmetic mean to avoid 0/0
  double ratio = dT1 / dT2;
  if (std::abs(ratio - 1.0) < 1e-6) {
    return (dT1 + dT2) / 2.0;
  }

  return (dT1 - dT2) / std::log(ratio);
}

double lmtd_counterflow(double T_hot_in, double T_hot_out, double T_cold_in,
                        double T_cold_out) {
  double dT1 = T_hot_in - T_cold_out;
  double dT2 = T_hot_out - T_cold_in;
  return lmtd(dT1, dT2);
}

double lmtd_parallelflow(double T_hot_in, double T_hot_out, double T_cold_in,
                         double T_cold_out) {
  double dT1 = T_hot_in - T_cold_in;
  double dT2 = T_hot_out - T_cold_out;
  return lmtd(dT1, dT2);
}

// -------------------------------------------------------------
// Temperature Profile Through Wall
// -------------------------------------------------------------

std::vector<double>
wall_temperature_profile(double T_hot, double T_cold, double h_hot,
                         double h_cold, const std::vector<double> &t_over_k,
                         double R_fouling, double &q) {

  if (h_hot <= 0 || h_cold <= 0) {
    throw std::invalid_argument(
        "wall_temperature_profile: HTCs must be positive");
  }

  // Total thermal resistance per unit area [m²·K/W]
  double R_total = 1.0 / h_hot + 1.0 / h_cold + R_fouling;
  for (double r : t_over_k) {
    if (r < 0) {
      throw std::invalid_argument(
          "wall_temperature_profile: t/k values must be non-negative");
    }
    R_total += r;
  }

  // Heat flux (constant through all layers)
  q = (T_hot - T_cold) / R_total;

  // Build temperature profile
  // Start from hot side, subtract temperature drops
  std::vector<double> temps;
  temps.reserve(t_over_k.size() + 1);

  // Hot surface temperature (after convective resistance)
  double T = T_hot - q / h_hot;
  temps.push_back(T);

  // Temperature at each layer interface
  for (double r : t_over_k) {
    T -= q * r;
    temps.push_back(T);
  }

  return temps;
}

std::vector<double>
wall_temperature_profile(double T_hot, double T_cold, double h_hot,
                         double h_cold, const std::vector<double> &t_over_k,
                         double R_fouling) {
  double q_unused;
  return wall_temperature_profile(T_hot, T_cold, h_hot, h_cold, t_over_k,
                                  R_fouling, q_unused);
}

// -------------------------------------------------------------
// Heat Flux from Measured Temperature
// -------------------------------------------------------------

double heat_flux_from_T_at_edge(double T_measured, std::size_t edge_idx,
                                double T_hot, double T_cold, double h_hot,
                                double h_cold,
                                const std::vector<double> &t_over_k) {

  std::size_t n_layers = t_over_k.size();
  if (edge_idx > n_layers) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_edge: edge_idx out of range (max = " +
        std::to_string(n_layers) + ")");
  }

  if (h_hot <= 0 || h_cold <= 0) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_edge: HTCs must be positive");
  }

  // Thermal resistance from hot bulk to edge
  double R_hot_to_edge = 1.0 / h_hot; // convective resistance
  for (std::size_t i = 0; i < edge_idx; ++i) {
    R_hot_to_edge += t_over_k[i];
  }

  // Thermal resistance from edge to cold bulk
  double R_edge_to_cold = 0.0;
  for (std::size_t i = edge_idx; i < n_layers; ++i) {
    R_edge_to_cold += t_over_k[i];
  }
  R_edge_to_cold += 1.0 / h_cold; // convective resistance

  // Heat flux: q = (T_hot - T_measured) / R_hot_to_edge
  //          = (T_measured - T_cold) / R_edge_to_cold
  // Both should give same q in steady state
  // Use the side with larger resistance for better numerical stability
  if (R_hot_to_edge >= R_edge_to_cold) {
    return (T_hot - T_measured) / R_hot_to_edge;
  } else {
    return (T_measured - T_cold) / R_edge_to_cold;
  }
}

double heat_flux_from_T_at_depth(double T_measured, double depth_from_hot,
                                 double T_hot, [[maybe_unused]] double T_cold,
                                 double h_hot, double h_cold,
                                 const std::vector<double> &thicknesses,
                                 const std::vector<double> &conductivities) {

  if (thicknesses.size() != conductivities.size()) {
    throw std::invalid_argument("heat_flux_from_T_at_depth: thicknesses and "
                                "conductivities must have same size");
  }
  if (depth_from_hot < 0) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_depth: depth must be non-negative");
  }
  if (h_hot <= 0 || h_cold <= 0) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_depth: HTCs must be positive");
  }

  // Calculate total wall thickness
  double total_thickness = 0.0;
  for (double t : thicknesses) {
    total_thickness += t;
  }

  if (depth_from_hot > total_thickness) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_depth: depth exceeds total wall thickness (" +
        std::to_string(total_thickness) + " m)");
  }

  // Find which layer the depth is in and compute resistance to that point
  double R_hot_to_point = 1.0 / h_hot; // convective resistance
  double cumulative_depth = 0.0;

  for (std::size_t i = 0; i < thicknesses.size(); ++i) {
    double layer_end = cumulative_depth + thicknesses[i];

    if (depth_from_hot <= layer_end) {
      // Point is in this layer
      double depth_in_layer = depth_from_hot - cumulative_depth;
      R_hot_to_point += depth_in_layer / conductivities[i];
      break;
    } else {
      // Point is past this layer
      R_hot_to_point += thicknesses[i] / conductivities[i];
      cumulative_depth = layer_end;
    }
  }

  // q = (T_hot - T_measured) / R_hot_to_point
  return (T_hot - T_measured) / R_hot_to_point;
}

double bulk_T_from_edge_T_and_q(double T_measured, std::size_t edge_idx,
                                double q, double h_hot, double h_cold,
                                const std::vector<double> &t_over_k,
                                const std::string &solve_for) {

  std::size_t n_layers = t_over_k.size();
  if (edge_idx > n_layers) {
    throw std::invalid_argument(
        "bulk_T_from_edge_T_and_q: edge_idx out of range (max = " +
        std::to_string(n_layers) + ")");
  }

  if (h_hot <= 0 || h_cold <= 0) {
    throw std::invalid_argument(
        "bulk_T_from_edge_T_and_q: HTCs must be positive");
  }

  if (solve_for == "hot") {
    // T_hot = T_measured + q * R_hot_to_edge
    double R_hot_to_edge = 1.0 / h_hot;
    for (std::size_t i = 0; i < edge_idx; ++i) {
      R_hot_to_edge += t_over_k[i];
    }
    return T_measured + q * R_hot_to_edge;
  } else if (solve_for == "cold") {
    // T_cold = T_measured - q * R_edge_to_cold
    double R_edge_to_cold = 0.0;
    for (std::size_t i = edge_idx; i < n_layers; ++i) {
      R_edge_to_cold += t_over_k[i];
    }
    R_edge_to_cold += 1.0 / h_cold;
    return T_measured - q * R_edge_to_cold;
  } else {
    throw std::invalid_argument(
        "bulk_T_from_edge_T_and_q: solve_for must be 'hot' or 'cold'");
  }
}

// -------------------------------------------------------------
// Temperature Sensitivity
// -------------------------------------------------------------

std::pair<double, double> dT_edge_dT_bulk(std::size_t edge_idx, double h_hot,
                                          double h_cold,
                                          const std::vector<double> &t_over_k) {

  std::size_t n_layers = t_over_k.size();
  if (edge_idx > n_layers) {
    throw std::invalid_argument(
        "dT_edge_dT_bulk: edge_idx out of range (max = " +
        std::to_string(n_layers) + ")");
  }

  if (h_hot <= 0 || h_cold <= 0) {
    throw std::invalid_argument("dT_edge_dT_bulk: HTCs must be positive");
  }

  // Total thermal resistance
  double R_total = 1.0 / h_hot + 1.0 / h_cold;
  for (double r : t_over_k) {
    R_total += r;
  }

  // Resistance from hot bulk to edge
  double R_hot_to_edge = 1.0 / h_hot;
  for (std::size_t i = 0; i < edge_idx; ++i) {
    R_hot_to_edge += t_over_k[i];
  }

  // Resistance from edge to cold bulk
  double R_edge_to_cold = R_total - R_hot_to_edge;

  // Sensitivities
  // T_edge = T_hot * (R_edge_to_cold / R_total) + T_cold * (R_hot_to_edge /
  // R_total)
  double dT_dT_hot = R_edge_to_cold / R_total;
  double dT_dT_cold = R_hot_to_edge / R_total;

  return {dT_dT_hot, dT_dT_cold};
}

double dT_edge_dT_hot(std::size_t edge_idx, double h_hot, double h_cold,
                      const std::vector<double> &t_over_k) {
  return dT_edge_dT_bulk(edge_idx, h_hot, h_cold, t_over_k).first;
}

double dT_edge_dT_cold(std::size_t edge_idx, double h_hot, double h_cold,
                       const std::vector<double> &t_over_k) {
  return dT_edge_dT_bulk(edge_idx, h_hot, h_cold, t_over_k).second;
}

double dT_edge_dq(std::size_t edge_idx, double h_hot,
                  const std::vector<double> &t_over_k) {

  std::size_t n_layers = t_over_k.size();
  if (edge_idx > n_layers) {
    throw std::invalid_argument("dT_edge_dq: edge_idx out of range (max = " +
                                std::to_string(n_layers) + ")");
  }

  if (h_hot <= 0) {
    throw std::invalid_argument("dT_edge_dq: h_hot must be positive");
  }

  // T_edge = T_hot - q * R_hot_to_edge
  // ∂T_edge/∂q = -R_hot_to_edge
  double R_hot_to_edge = 1.0 / h_hot;
  for (std::size_t i = 0; i < edge_idx; ++i) {
    R_hot_to_edge += t_over_k[i];
  }

  return -R_hot_to_edge;
}

// -------------------------------------------------------------
// NTU-Effectiveness Method
// -------------------------------------------------------------

double effectiveness_counterflow(double NTU, double C_r) {
  if (NTU < 0) {
    throw std::invalid_argument(
        "effectiveness_counterflow: NTU must be non-negative");
  }
  if (C_r < 0 || C_r > 1) {
    throw std::invalid_argument(
        "effectiveness_counterflow: C_r must be in [0, 1]");
  }

  // Special case: C_r = 1 (balanced heat exchanger)
  if (std::abs(C_r - 1.0) < 1e-10) {
    return NTU / (1.0 + NTU);
  }

  // Special case: C_r = 0 (one fluid has infinite capacity, e.g.,
  // condenser/evaporator)
  if (C_r < 1e-10) {
    return 1.0 - std::exp(-NTU);
  }

  // General case
  double exp_term = std::exp(-NTU * (1.0 - C_r));
  return (1.0 - exp_term) / (1.0 - C_r * exp_term);
}

double effectiveness_parallelflow(double NTU, double C_r) {
  if (NTU < 0) {
    throw std::invalid_argument(
        "effectiveness_parallelflow: NTU must be non-negative");
  }
  if (C_r < 0 || C_r > 1) {
    throw std::invalid_argument(
        "effectiveness_parallelflow: C_r must be in [0, 1]");
  }

  // Special case: C_r = 0
  if (C_r < 1e-10) {
    return 1.0 - std::exp(-NTU);
  }

  // General case
  return (1.0 - std::exp(-NTU * (1.0 + C_r))) / (1.0 + C_r);
}

// -------------------------------------------------------------
// State-based convenience functions
// -------------------------------------------------------------

double nusselt_pipe(const State &s, double velocity, double diameter,
                    bool heating, double roughness) {
  if (velocity <= 0) {
    throw std::invalid_argument("nusselt_pipe: velocity must be positive");
  }
  if (diameter <= 0) {
    throw std::invalid_argument("nusselt_pipe: diameter must be positive");
  }

  // Re = ρ * V * D / μ
  double rho = s.rho();
  double mu = s.mu();
  double Re = rho * velocity * diameter / mu;
  double Pr = s.Pr();

  // Laminar flow
  if (Re < 2300) {
    return heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
  }

  // Transition region (2300 < Re < 10000): use Gnielinski
  // Turbulent (Re > 10000): Gnielinski is still good, or could use
  // Dittus-Boelter

  // Get friction factor
  double e_D = roughness / diameter;
  double f;
  if (e_D > 0 && Re > 4000) {
    f = friction_colebrook(Re, e_D);
  } else {
    f = friction_petukhov(std::max(Re, 3000.0));
  }

  return nusselt_gnielinski(Re, Pr, f);
}

double htc_pipe(const State &s, double velocity, double diameter, bool heating,
                double roughness) {
  double Nu = nusselt_pipe(s, velocity, diameter, heating, roughness);
  double k = s.k();
  return htc_from_nusselt(Nu, k, diameter);
}

// Composite function: compute HTC from thermodynamic state
// Returns: (h [W/(m²·K)], Nu [-], Re [-])
std::tuple<double, double, double>
htc_pipe(double T, double P, const std::vector<double> &X, double velocity,
         double diameter, const std::string &correlation, bool heating,
         double mu_ratio, double roughness) {
  if (velocity <= 0) {
    throw std::invalid_argument("htc_pipe: velocity must be positive");
  }
  if (diameter <= 0) {
    throw std::invalid_argument("htc_pipe: diameter must be positive");
  }
  if (T <= 0) {
    throw std::invalid_argument("htc_pipe: temperature must be positive");
  }
  if (P <= 0) {
    throw std::invalid_argument("htc_pipe: pressure must be positive");
  }

  // Compute thermodynamic and transport properties
  double rho = density(T, P, X);
  double mu = viscosity(T, P, X);
  double k = thermal_conductivity(T, P, X);
  [[maybe_unused]] double cp = cp_mass(T, X);
  double Pr = prandtl(T, P, X);

  // Reynolds number
  double Re = rho * velocity * diameter / mu;

  // Select correlation and compute Nusselt number
  double Nu;

  if (correlation == "gnielinski") {
    // Gnielinski: 2300 < Re < 5e6, 0.5 < Pr < 2000
    if (Re < 2300) {
      // Laminar flow
      Nu = heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
    } else {
      // Get friction factor
      double e_D = roughness / diameter;
      double f;
      if (e_D > 0 && Re > 4000) {
        f = friction_colebrook(Re, e_D);
      } else {
        f = friction_petukhov(std::max(Re, 3000.0));
      }
      Nu = nusselt_gnielinski(Re, Pr, f);
    }
  } else if (correlation == "dittus_boelter") {
    // Dittus-Boelter: Re > 10000, 0.6 < Pr < 160
    if (Re < 2300) {
      Nu = heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
    } else if (Re < 10000) {
      throw std::invalid_argument(
          "htc_pipe: Dittus-Boelter requires Re > 10000 (got Re=" +
          std::to_string(Re) + "). Use 'gnielinski' for transition region.");
    } else {
      Nu = nusselt_dittus_boelter(Re, Pr, heating);
    }
  } else if (correlation == "sieder_tate") {
    // Sieder-Tate: Re > 10000, 0.7 < Pr < 16700
    if (Re < 2300) {
      Nu = heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
    } else if (Re < 10000) {
      throw std::invalid_argument(
          "htc_pipe: Sieder-Tate requires Re > 10000 (got Re=" +
          std::to_string(Re) + "). Use 'gnielinski' for transition region.");
    } else {
      Nu = nusselt_sieder_tate(Re, Pr, mu_ratio);
    }
  } else if (correlation == "petukhov") {
    // Petukhov: 1e4 < Re < 5e6, 0.5 < Pr < 2000
    if (Re < 2300) {
      Nu = heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
    } else if (Re < 1e4) {
      throw std::invalid_argument(
          "htc_pipe: Petukhov requires Re > 10000 (got Re=" +
          std::to_string(Re) + "). Use 'gnielinski' for transition region.");
    } else {
      double f = friction_petukhov(Re);
      Nu = nusselt_petukhov(Re, Pr, f);
    }
  } else {
    throw std::invalid_argument("htc_pipe: unknown correlation '" +
                                correlation + "'. " +
                                "Valid options: 'gnielinski', "
                                "'dittus_boelter', 'sieder_tate', 'petukhov'");
  }

  // Heat transfer coefficient
  double h = htc_from_nusselt(Nu, k, diameter);

  return std::make_tuple(h, Nu, Re);
}

// -------------------------------------------------------------
// Internal helper: build ChannelResult from common fields
// -------------------------------------------------------------

static ChannelResult make_channel_result(double h, double Nu, double Re,
                                         double Pr, double f, double dP,
                                         double M, double T_aw, double T_wall) {
  double q = std::numeric_limits<double>::quiet_NaN();
  if (std::isfinite(T_wall)) {
    q = h * (T_aw - T_wall);
  }
  return ChannelResult{h, Nu, Re, Pr, f, dP, M, T_aw, q};
}


// -------------------------------------------------------------
// channel_smooth
// -------------------------------------------------------------

ChannelResult channel_smooth(double T, double P, const std::vector<double> &X,
                             double velocity, double diameter, double length,
                             double T_wall, const std::string &correlation,
                             bool heating, double mu_ratio, double roughness,
                             double Nu_multiplier, double f_multiplier) {
  if (velocity < 0.0) {
    throw std::invalid_argument(
        "channel_smooth: velocity must be non-negative");
  }
  if (diameter <= 0.0) {
    throw std::invalid_argument("channel_smooth: diameter must be positive");
  }
  if (length < 0.0) {
    throw std::invalid_argument("channel_smooth: length must be non-negative");
  }

  // Fluid properties and their derivatives
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  auto [mu, dmu_dT, dmu_dP] = solver::viscosity_and_jacobians(T, P, X);
  double k = thermal_conductivity(T, P, X);
  double Pr = prandtl(T, P, X);

  // Reynolds number
  double Re = (velocity > 0.0) ? rho * velocity * diameter / mu : 0.0;

  // Friction factor
  double e_D = roughness / diameter;
  double f;
  if (Re < 2300.0) {
    f = (Re > 0.0) ? 64.0 / Re : 0.0;
  } else if (e_D > 0.0 && Re > 4000.0) {
    f = friction_colebrook(Re, e_D);
  } else {
    f = friction_petukhov(std::max(Re, 3000.0));
  }

  // Apply f_multiplier before Nu computation so that correlations where
  // f appears in the Nu formula (Gnielinski, Petukhov) use the multiplied
  // friction factor consistently with the Jacobian FD stencils.
  f *= f_multiplier;

  // Nusselt number
  double Nu;
  if (Re < 2300.0) {
    Nu = heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
  } else if (correlation == "gnielinski") {
    Nu = nusselt_gnielinski(Re, Pr, f);
  } else if (correlation == "dittus_boelter") {
    if (Re < 10000.0) {
      throw std::invalid_argument(
          "channel_smooth: dittus_boelter requires Re > 10000 (got " +
          std::to_string(Re) + "). Use gnielinski for transition region.");
    }
    Nu = nusselt_dittus_boelter(Re, Pr, heating);
  } else if (correlation == "sieder_tate") {
    if (Re < 10000.0) {
      throw std::invalid_argument(
          "channel_smooth: sieder_tate requires Re > 10000 (got " +
          std::to_string(Re) + "). Use gnielinski for transition region.");
    }
    Nu = nusselt_sieder_tate(Re, Pr, mu_ratio);
  } else if (correlation == "petukhov") {
    if (Re < 10000.0) {
      throw std::invalid_argument(
          "channel_smooth: petukhov requires Re > 10000 (got " +
          std::to_string(Re) + "). Use gnielinski for transition region.");
    }
    Nu = nusselt_petukhov(Re, Pr, f);
  } else {
    throw std::invalid_argument("channel_smooth: unknown correlation '" +
                                correlation + "'");
  }

  // Apply Nu multiplier
  Nu *= Nu_multiplier;

  double h = htc_from_nusselt(Nu, k, diameter);
  double dP = (velocity > 0.0)
                  ? f * (length / diameter) * (rho * velocity * velocity / 2.0)
                  : 0.0;

  // Mach number and T_aw — always continuous
  double a = speed_of_sound(T, X);
  double M = (a > 0.0) ? velocity / a : 0.0;
  const bool turbulent_flow = true;
  double T_aw = T_adiabatic_wall(T, velocity, T, P, X, turbulent_flow);

  // Compute Jacobians
  // Note: For network solver, mdot is the independent variable, not velocity
  // velocity = mdot / (rho * A_cross), where A_cross = pi/4 * D^2
  // We compute derivatives w.r.t. mdot and T

  ChannelResult result = make_channel_result(h, Nu, Re, Pr, f, dP, M, T_aw, T_wall);

  // Only compute Jacobians if flow exists
  if (velocity > 0.0 && Re > 0.0) {
    double A_cross = M_PI / 4.0 * diameter * diameter;
    double mdot = rho * velocity * A_cross;

    // Thermal conductivity and Prandtl derivatives via central FD
    const double eps_T_props = 1e-3;
    double k_plus = thermal_conductivity(T + eps_T_props, P, X);
    double k_minus = thermal_conductivity(T - eps_T_props, P, X);
    double dk_dT = (k_plus - k_minus) / (2.0 * eps_T_props);

    double Pr_plus = prandtl(T + eps_T_props, P, X);
    double Pr_minus = prandtl(T - eps_T_props, P, X);
    double dPr_dT = (Pr_plus - Pr_minus) / (2.0 * eps_T_props);

    // Chain rule: dRe/dmdot and dRe/dT
    // Re = rho * v * D / mu = mdot * D / (A_cross * mu)
    // At constant mdot: dRe/dT = -Re * dmu/dT / mu (no rho dependence)
    double dRe_dmdot = diameter / (A_cross * mu);
    double dRe_dT = -Re * dmu_dT / mu;

    // Get Nu and f derivatives w.r.t. Re and Pr (these depend on correlation)
    double dNu_dRe = 0.0;
    double dNu_dPr = 0.0;
    double df_dRe = 0.0;

    if (Re < 2300.0) {
      // Laminar: f = 64/Re, so df/dRe = -64/Re^2 = -f/Re
      if (Re > 0.0) {
        df_dRe = -f / Re;
      }
      // Nu is constant in laminar regime, dNu/dRe = 0
    } else {
      // Friction factor derivative (needed first for Gnielinski/Petukhov)
      if (e_D > 0.0 && Re > 4000.0) {
        // Colebrook - use finite difference
        double eps = std::max(1.0, Re * 1e-6);
        double f_plus = friction_colebrook(Re + eps, e_D);
        double f_minus = friction_colebrook(Re - eps, e_D);
        df_dRe = (f_plus - f_minus) / (2.0 * eps) * f_multiplier;
      } else {
        // Petukhov: f_raw = (0.79*ln(Re) - 1.64)^(-2)
        // df_raw/dRe = -2 * f_raw^1.5 * 0.79 / Re
        // df/dRe = f_mult * df_raw/dRe  (f_raw = f / f_multiplier)
        df_dRe = -2.0 * std::pow(f / f_multiplier, 1.5) * 0.79 / Re * f_multiplier;
      }

      // Turbulent correlations have Re and Pr derivatives
      if (correlation == "gnielinski") {
        // Gnielinski: Nu = Nu(Re, Pr, f), so dNu/dRe_total = ∂Nu/∂Re + ∂Nu/∂f · df/dRe
        // Use finite difference for total derivative
        double eps = std::max(1.0, Re * 1e-6);
        double f_plus = (e_D > 0.0 && Re + eps > 4000.0) ? friction_colebrook(Re + eps, e_D) : friction_petukhov(Re + eps);
        double f_minus = (e_D > 0.0 && Re - eps > 4000.0) ? friction_colebrook(Re - eps, e_D) : friction_petukhov(Re - eps);
        f_plus *= f_multiplier;
        f_minus *= f_multiplier;
        double Nu_plus = nusselt_gnielinski(Re + eps, Pr, f_plus);
        double Nu_minus = nusselt_gnielinski(Re - eps, Pr, f_minus);
        dNu_dRe = (Nu_plus - Nu_minus) / (2.0 * eps) * Nu_multiplier;

        // dNu/dPr via central FD (f already includes f_multiplier)
        double eps_Pr = std::max(1e-6, Pr * 1e-6);
        double Nu_Pr_plus = nusselt_gnielinski(Re, Pr + eps_Pr, f);
        double Nu_Pr_minus = nusselt_gnielinski(Re, Pr - eps_Pr, f);
        dNu_dPr = (Nu_Pr_plus - Nu_Pr_minus) / (2.0 * eps_Pr) * Nu_multiplier;
      } else if (correlation == "dittus_boelter") {
        // Nu = 0.023 * Re^0.8 * Pr^n  =>  dNu/dRe = 0.8 * Nu / Re, dNu/dPr = n * Nu / Pr
        dNu_dRe = 0.8 * Nu / Re;
        double n = heating ? 0.4 : 0.3;
        dNu_dPr = n * Nu / Pr;
      } else if (correlation == "sieder_tate") {
        // Nu = 0.027 * Re^0.8 * Pr^(1/3) * mu_ratio^0.14  =>  dNu/dRe = 0.8 * Nu / Re, dNu/dPr = (1/3) * Nu / Pr
        dNu_dRe = 0.8 * Nu / Re;
        dNu_dPr = (1.0 / 3.0) * Nu / Pr;
      } else if (correlation == "petukhov") {
        // Petukhov: Nu = Nu(Re, Pr, f), so dNu/dRe_total = ∂Nu/∂Re + ∂Nu/∂f · df/dRe
        double eps = std::max(1.0, Re * 1e-6);
        double f_plus = friction_petukhov(Re + eps) * f_multiplier;
        double f_minus = friction_petukhov(Re - eps) * f_multiplier;
        double Nu_plus = nusselt_petukhov(Re + eps, Pr, f_plus);
        double Nu_minus = nusselt_petukhov(Re - eps, Pr, f_minus);
        dNu_dRe = (Nu_plus - Nu_minus) / (2.0 * eps) * Nu_multiplier;

        // dNu/dPr via central FD (f already includes f_multiplier)
        double eps_Pr = std::max(1e-6, Pr * 1e-6);
        double Nu_Pr_plus = nusselt_petukhov(Re, Pr + eps_Pr, f);
        double Nu_Pr_minus = nusselt_petukhov(Re, Pr - eps_Pr, f);
        dNu_dPr = (Nu_Pr_plus - Nu_Pr_minus) / (2.0 * eps_Pr) * Nu_multiplier;
      }
    }
    // Laminar: Nu and f are constant or 1/Re (derivatives handled separately if needed)

    // dh/dmdot = (dh/dNu) * (dNu/dRe) * (dRe/dmdot) = (k/D) * dNu/dRe * dRe/dmdot
    result.dh_dmdot = (k / diameter) * dNu_dRe * dRe_dmdot;

    // dh/dT: Full chain rule h = Nu(Re, Pr) * k(T) / D
    // dh/dT = (1/D) * [k * (dNu/dRe * dRe/dT + dNu/dPr * dPr/dT) + Nu * dk/dT]
    result.dh_dT = (1.0 / diameter) * (k * (dNu_dRe * dRe_dT + dNu_dPr * dPr_dT) + Nu * dk_dT);

    // ddP/dmdot: dP = f * (L/D) * (rho * v^2 / 2)
    // v = mdot / (rho * A), so v^2 = mdot^2 / (rho^2 * A^2)
    // dP = f * (L/D) * mdot^2 / (2 * rho * A^2)
    // d(dP)/dmdot = (df/dRe * dRe/dmdot) * (L/D) * mdot^2/(2*rho*A^2) + f*(L/D)*mdot/(rho*A^2)
    double dP_factor = (length / diameter) / (2.0 * rho * A_cross * A_cross);
    result.ddP_dmdot = df_dRe * dRe_dmdot * dP_factor * mdot * mdot
                     + f * (length / diameter) * mdot / (rho * A_cross * A_cross);

    // d(dP)/dT: includes df/dT and drho/dT terms
    // dP = f * (L/D) * mdot^2 / (2*rho*A^2), so d(dP)/dT has df/dT and d(1/rho)/dT = -drho/dT/rho^2
    result.ddP_dT = df_dRe * dRe_dT * dP_factor * mdot * mdot
                  - f * (length / diameter) * velocity * velocity / 2.0 * drho_dT;

    // dq/dmdot and dq/dT: q = h * (T_aw - T_wall)
    // Full product rule: dq/dx = dh/dx * (T_aw - T_wall) + h * dT_aw/dx
    // T_aw derivatives — always computed so Python relay can use them
    // T_aw = T + r * v^2 / (2*cp), where v = mdot/(rho*A), r = Pr^(1/3)
    double cp_mass_val = cp_mass(T, X);
    double r = std::cbrt(Pr);

    // dT_aw/dmdot = r * v / (cp * rho * A)
    result.dT_aw_dmdot = r * velocity / (cp_mass_val * rho * A_cross);

    // dT_aw/dT at constant mdot via central FD
    const double eps_T_aw = 0.5;
    double rho_plus = density(T + eps_T_aw, P, X);
    double rho_minus = density(T - eps_T_aw, P, X);
    double v_plus = mdot / (rho_plus * A_cross);
    double v_minus = mdot / (rho_minus * A_cross);
    double T_aw_plus = T_adiabatic_wall(T + eps_T_aw, v_plus, T + eps_T_aw, P, X, turbulent_flow);
    double T_aw_minus = T_adiabatic_wall(T - eps_T_aw, v_minus, T - eps_T_aw, P, X, turbulent_flow);
    result.dT_aw_dT = (T_aw_plus - T_aw_minus) / (2.0 * eps_T_aw);

    if (std::isfinite(T_wall)) {
      double dT_diff = T_aw - T_wall;
      result.dq_dmdot = result.dh_dmdot * dT_diff + h * result.dT_aw_dmdot;
      result.dq_dT = result.dh_dT * dT_diff + h * result.dT_aw_dT;
      result.dq_dT_wall = -h;
    }
  }

  return result;
}

// -------------------------------------------------------------
// channel_ribbed
// -------------------------------------------------------------

ChannelResult channel_ribbed(double T, double P, const std::vector<double> &X,
                             double velocity, double diameter, double length,
                             double e_D, double pitch_to_height,
                             double alpha_deg, double T_wall, bool heating,
                             double Nu_multiplier, double f_multiplier) {
  // Get smooth-pipe baseline (Gnielinski, no roughness)
  // Pass multipliers to baseline - they will be applied after rib factors
  ChannelResult base = channel_smooth(T, P, X, velocity, diameter, length,
                                      T_wall, "gnielinski", heating, 1.0, 0.0,
                                      Nu_multiplier, f_multiplier);

  // Apply rib multipliers to the baseline (before user multipliers)
  double enh = combaero::cooling::rib_enhancement_factor(e_D, pitch_to_height,
                                                         alpha_deg);
  double fmul =
      combaero::cooling::rib_friction_multiplier(e_D, pitch_to_height);

  // Rib factors are applied to the smooth baseline, then user multipliers
  double Nu_rib = base.Nu * enh;
  double f_rib = base.f * fmul;

  double rho = density(T, P, X);
  double k = thermal_conductivity(T, P, X);
  double h_rib = htc_from_nusselt(Nu_rib, k, diameter);
  double dP_rib = (velocity > 0.0) ? f_rib * (length / diameter) *
                                         (rho * velocity * velocity / 2.0)
                                   : 0.0;

  ChannelResult result = make_channel_result(h_rib, Nu_rib, base.Re, base.Pr, f_rib, dP_rib,
                             base.M, base.T_aw, T_wall);

  // Propagate Jacobians: scale by rib enhancement factors
  // dh_rib/dmdot = enh * dh_base/dmdot
  result.dh_dmdot = enh * base.dh_dmdot;
  result.dh_dT = enh * base.dh_dT;
  result.ddP_dmdot = fmul * base.ddP_dmdot;
  result.ddP_dT = fmul * base.ddP_dT;
  result.dT_aw_dmdot = base.dT_aw_dmdot;
  result.dT_aw_dT = base.dT_aw_dT;
  result.dq_dmdot = enh * base.dq_dmdot;
  result.dq_dT = enh * base.dq_dT;
  result.dq_dT_wall = base.dq_dT_wall * enh;  // = -h_rib

  return result;
}

// -------------------------------------------------------------
// channel_dimpled
// -------------------------------------------------------------

ChannelResult channel_dimpled(double T, double P, const std::vector<double> &X,
                              double velocity, double diameter, double length,
                              double d_Dh, double h_d, double S_d,
                              double T_wall, bool heating,
                              double Nu_multiplier, double f_multiplier) {
  // Get smooth-pipe baseline
  ChannelResult base = channel_smooth(T, P, X, velocity, diameter, length,
                                      T_wall, "gnielinski", heating, 1.0, 0.0,
                                      Nu_multiplier, f_multiplier);

  // Apply dimple multipliers
  double enh =
      combaero::cooling::dimple_nusselt_enhancement(base.Re, d_Dh, h_d, S_d);
  double fmul =
      combaero::cooling::dimple_friction_multiplier(base.Re, d_Dh, h_d);

  double Nu_dim = base.Nu * enh;
  double f_dim = base.f * fmul;

  double rho = density(T, P, X);
  double k = thermal_conductivity(T, P, X);
  double h_dim = htc_from_nusselt(Nu_dim, k, diameter);
  double dP_dim = (velocity > 0.0) ? f_dim * (length / diameter) *
                                         (rho * velocity * velocity / 2.0)
                                   : 0.0;

  ChannelResult result = make_channel_result(h_dim, Nu_dim, base.Re, base.Pr, f_dim, dP_dim,
                             base.M, base.T_aw, T_wall);

  // Propagate Jacobians: scale by dimple enhancement factors
  result.dh_dmdot = enh * base.dh_dmdot;
  result.dh_dT = enh * base.dh_dT;
  result.ddP_dmdot = fmul * base.ddP_dmdot;
  result.ddP_dT = fmul * base.ddP_dT;
  result.dT_aw_dmdot = base.dT_aw_dmdot;
  result.dT_aw_dT = base.dT_aw_dT;
  result.dq_dmdot = enh * base.dq_dmdot;
  result.dq_dT = enh * base.dq_dT;
  result.dq_dT_wall = base.dq_dT_wall * enh;

  return result;
}

// -------------------------------------------------------------
// channel_pin_fin
// -------------------------------------------------------------

ChannelResult channel_pin_fin(double T, double P, const std::vector<double> &X,
                              double velocity, double channel_height,
                              double pin_diameter, double S_D, double X_D,
                              int N_rows, double T_wall, bool is_staggered,
                              double Nu_multiplier, double f_multiplier) {
  if (velocity < 0.0) {
    throw std::invalid_argument(
        "channel_pin_fin: velocity must be non-negative");
  }
  if (pin_diameter <= 0.0) {
    throw std::invalid_argument(
        "channel_pin_fin: pin_diameter must be positive");
  }
  if (channel_height <= 0.0) {
    throw std::invalid_argument(
        "channel_pin_fin: channel_height must be positive");
  }
  if (N_rows < 1) {
    throw std::invalid_argument("channel_pin_fin: N_rows must be >= 1");
  }
  if (S_D <= 1.0) {
    throw std::invalid_argument(
        "channel_pin_fin: S_D must be > 1 (minimum cross-section constraint)");
  }

  // Fluid properties and their derivatives
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  auto [mu, dmu_dT, dmu_dP] = solver::viscosity_and_jacobians(T, P, X);
  double k = thermal_conductivity(T, P, X);
  double Pr = prandtl(T, P, X);

  // Re based on pin diameter and approach velocity
  double Re_d = (velocity > 0.0) ? rho * velocity * pin_diameter / mu : 0.0;

  // L_D = channel_height / pin_diameter
  double L_D = channel_height / pin_diameter;

  // Nu from Metzger correlation
  double Nu =
      combaero::cooling::pin_fin_nusselt(Re_d, Pr, L_D, S_D, X_D, is_staggered);
  Nu *= Nu_multiplier;
  double h = htc_from_nusselt(Nu, k, pin_diameter);

  // Pressure drop: dP = N_rows * f_pin * (rho * v_max^2 / 2)
  // v_max at minimum cross-section: v_max = v * S_D / (S_D - 1)
  double v_max = (S_D > 1.0) ? velocity * S_D / (S_D - 1.0) : velocity;
  double f_pin = (Re_d > 0.0)
                     ? combaero::cooling::pin_fin_friction(Re_d, is_staggered)
                     : 0.0;
  f_pin *= f_multiplier;
  double dP = static_cast<double>(N_rows) * f_pin * (rho * v_max * v_max / 2.0);

  // Mach and T_aw based on approach velocity
  double a = speed_of_sound(T, X);
  double M = (a > 0.0) ? velocity / a : 0.0;
  const bool turbulent_flow = true;
  double T_aw = T_adiabatic_wall(T, velocity, T, P, X, turbulent_flow);

  ChannelResult result = make_channel_result(h, Nu, Re_d, Pr, f_pin, dP, M, T_aw, T_wall);

  // Compute Jacobians (simplified - captures dominant Re dependence)
  if (velocity > 0.0 && Re_d > 0.0) {
    double A_cross = channel_height * pin_diameter * (S_D - 1.0) / S_D;  // Approx flow area
    double mdot = rho * velocity * A_cross;

    // Thermal conductivity and Prandtl derivatives via central FD
    const double eps_T_props = 1e-3;
    double k_plus = thermal_conductivity(T + eps_T_props, P, X);
    double k_minus = thermal_conductivity(T - eps_T_props, P, X);
    double dk_dT = (k_plus - k_minus) / (2.0 * eps_T_props);

    double Pr_plus = prandtl(T + eps_T_props, P, X);
    double Pr_minus = prandtl(T - eps_T_props, P, X);
    double dPr_dT = (Pr_plus - Pr_minus) / (2.0 * eps_T_props);

    // dRe/dmdot and dRe/dT
    // At constant mdot: dRe/dT = -Re * dmu/dT / mu (no rho dependence)
    double dRe_dmdot = pin_diameter / (A_cross * mu);
    double dRe_dT = -Re_d * dmu_dT / mu;

    // Simplified: assume Nu ~ Re^0.7 and f ~ Re^(-0.2) (typical for pin fins)
    double dNu_dRe = 0.7 * Nu / Re_d;
    double dNu_dPr = 0.4 * Nu / Pr;  // Typical Pr exponent for pin fins
    double df_dRe = -0.2 * f_pin / Re_d;

    result.dh_dmdot = (k / pin_diameter) * dNu_dRe * dRe_dmdot;
    // Full chain rule: dh/dT = (1/D) * [k * (dNu/dRe * dRe/dT + dNu/dPr * dPr/dT) + Nu * dk/dT]
    result.dh_dT = (1.0 / pin_diameter) * (k * (dNu_dRe * dRe_dT + dNu_dPr * dPr_dT) + Nu * dk_dT);

    // ddP/dmdot for pin fins
    double v_max_factor = S_D / (S_D - 1.0);
    result.ddP_dmdot = static_cast<double>(N_rows) *
                       (df_dRe * dRe_dmdot * rho * v_max * v_max / 2.0 +
                        f_pin * v_max_factor * v_max_factor * mdot / (rho * A_cross * A_cross));
    result.ddP_dT = static_cast<double>(N_rows) * df_dRe * dRe_dT * rho * v_max * v_max / 2.0;

    // T_aw derivatives — always computed so Python relay can use them
    double cp_mass_val = cp_mass(T, X);
    double r = std::cbrt(Pr);
    result.dT_aw_dmdot = r * velocity / (cp_mass_val * rho * A_cross);

    const double eps_T_aw = 0.5;
    double rho_plus = density(T + eps_T_aw, P, X);
    double rho_minus = density(T - eps_T_aw, P, X);
    double v_plus = mdot / (rho_plus * A_cross);
    double v_minus = mdot / (rho_minus * A_cross);
    const bool turbulent_flow = true;
    double T_aw_plus = T_adiabatic_wall(T + eps_T_aw, v_plus, T + eps_T_aw, P, X, turbulent_flow);
    double T_aw_minus = T_adiabatic_wall(T - eps_T_aw, v_minus, T - eps_T_aw, P, X, turbulent_flow);
    result.dT_aw_dT = (T_aw_plus - T_aw_minus) / (2.0 * eps_T_aw);

    if (std::isfinite(T_wall)) {
      double dT_diff = T_aw - T_wall;
      result.dq_dmdot = result.dh_dmdot * dT_diff + h * result.dT_aw_dmdot;
      result.dq_dT = result.dh_dT * dT_diff + h * result.dT_aw_dT;
      result.dq_dT_wall = -h;
    }
  }

  return result;
}

// -------------------------------------------------------------
// channel_impingement
// -------------------------------------------------------------

ChannelResult channel_impingement(double T, double P,
                                  const std::vector<double> &X, double mdot_jet,
                                  double d_jet, double z_D, double x_D,
                                  double y_D, double A_target, double T_wall,
                                  double Cd_jet,
                                  double Nu_multiplier, double f_multiplier) {
  if (mdot_jet < 0.0) {
    throw std::invalid_argument(
        "channel_impingement: mdot_jet must be non-negative");
  }
  if (d_jet <= 0.0) {
    throw std::invalid_argument("channel_impingement: d_jet must be positive");
  }
  if (A_target <= 0.0) {
    throw std::invalid_argument(
        "channel_impingement: A_target must be positive");
  }
  if (Cd_jet <= 0.0 || Cd_jet > 1.0) {
    throw std::invalid_argument(
        "channel_impingement: Cd_jet must be in (0, 1]");
  }

  // Fluid properties and their derivatives
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  auto [mu, dmu_dT, dmu_dP] = solver::viscosity_and_jacobians(T, P, X);
  double k = thermal_conductivity(T, P, X);
  double Pr = prandtl(T, P, X);

  // Jet area and velocity
  double A_jet = M_PI / 4.0 * d_jet * d_jet;
  double v_jet = (rho > 0.0 && A_jet > 0.0) ? mdot_jet / (rho * A_jet) : 0.0;
  double Re_jet = (mu > 0.0) ? rho * v_jet * d_jet / mu : 0.0;

  // Nu from Florschuetz / Martin correlation
  double Nu = combaero::cooling::impingement_nusselt(Re_jet, Pr, z_D, x_D, y_D);
  Nu *= Nu_multiplier;

  // h averaged over target area: h = Nu * k / d_jet
  // (Nu is already an area-averaged value from the correlation)
  double h = htc_from_nusselt(Nu, k, d_jet);

  // Pressure drop across jet plate: dP = (v_jet/Cd)^2 * rho/2
  // Equivalent loss coefficient f = 1/Cd^2
  double f = 1.0 / (Cd_jet * Cd_jet);
  f *= f_multiplier;
  double dP = f * rho * v_jet * v_jet / 2.0;

  // Mach and T_aw based on jet velocity
  double a = speed_of_sound(T, X);
  double M = (a > 0.0) ? v_jet / a : 0.0;
  const bool turbulent_flow = true;
  double T_aw = T_adiabatic_wall(T, v_jet, T, P, X, turbulent_flow);

  ChannelResult result = make_channel_result(h, Nu, Re_jet, Pr, f, dP, M, T_aw, T_wall);

  // Compute Jacobians (simplified - captures dominant Re dependence)
  if (mdot_jet > 0.0 && Re_jet > 0.0) {
    // Thermal conductivity and Prandtl derivatives via central FD
    const double eps_T_props = 1e-3;
    double k_plus = thermal_conductivity(T + eps_T_props, P, X);
    double k_minus = thermal_conductivity(T - eps_T_props, P, X);
    double dk_dT = (k_plus - k_minus) / (2.0 * eps_T_props);

    double Pr_plus = prandtl(T + eps_T_props, P, X);
    double Pr_minus = prandtl(T - eps_T_props, P, X);
    double dPr_dT = (Pr_plus - Pr_minus) / (2.0 * eps_T_props);

    // dRe/dmdot: Re = mdot * d / (A_jet * mu)
    // At constant mdot: dRe/dT = -Re * dmu/dT / mu (no rho dependence)
    double dRe_dmdot = d_jet / (A_jet * mu);
    double dRe_dT = -Re_jet * dmu_dT / mu;

    // Simplified: assume Nu ~ Re^0.7 (typical for impingement)
    double dNu_dRe = 0.7 * Nu / Re_jet;
    double dNu_dPr = 0.4 * Nu / Pr;  // Typical Pr exponent for impingement

    result.dh_dmdot = (k / d_jet) * dNu_dRe * dRe_dmdot;
    // Full chain rule: dh/dT = (1/D) * [k * (dNu/dRe * dRe/dT + dNu/dPr * dPr/dT) + Nu * dk/dT]
    result.dh_dT = (1.0 / d_jet) * (k * (dNu_dRe * dRe_dT + dNu_dPr * dPr_dT) + Nu * dk_dT);

    // ddP/dmdot: dP = f * rho * (mdot/(rho*A_jet))^2 / 2 = f * mdot^2 / (2*rho*A_jet^2)
    result.ddP_dmdot = f * mdot_jet / (rho * A_jet * A_jet);
    // Fix: remove spurious /rho in ddP/dT
    result.ddP_dT = -f * v_jet * v_jet / 2.0 * drho_dT;

    // T_aw derivatives — always computed so Python relay can use them
    double cp_mass_val = cp_mass(T, X);
    double r = std::cbrt(Pr);
    result.dT_aw_dmdot = r * v_jet / (cp_mass_val * rho * A_jet);

    const double eps_T_aw = 0.5;
    double rho_plus = density(T + eps_T_aw, P, X);
    double rho_minus = density(T - eps_T_aw, P, X);
    double v_plus = mdot_jet / (rho_plus * A_jet);
    double v_minus = mdot_jet / (rho_minus * A_jet);
    const bool turbulent_flow = true;
    double T_aw_plus = T_adiabatic_wall(T + eps_T_aw, v_plus, T + eps_T_aw, P, X, turbulent_flow);
    double T_aw_minus = T_adiabatic_wall(T - eps_T_aw, v_minus, T - eps_T_aw, P, X, turbulent_flow);
    result.dT_aw_dT = (T_aw_plus - T_aw_minus) / (2.0 * eps_T_aw);

    if (std::isfinite(T_wall)) {
      double dT_diff = T_aw - T_wall;
      result.dq_dmdot = result.dh_dmdot * dT_diff + h * result.dT_aw_dmdot;
      result.dq_dT = result.dh_dT * dT_diff + h * result.dT_aw_dT;
      result.dq_dT_wall = -h;
    }
  }

  return result;
}

// -------------------------------------------------------------
// wall_coupling_and_jacobian
// -------------------------------------------------------------

WallCouplingResult wall_coupling_and_jacobian(
    double h_a, double T_aw_a,
    double h_b, double T_aw_b,
    double t_over_k,
    double A) {
  WallCouplingResult result;

  // Overall HTC: U = 1 / (1/h_a + t/k + 1/h_b)
  double R_total = 1.0 / h_a + t_over_k + 1.0 / h_b;
  double U = 1.0 / R_total;

  // Heat transfer rate: Q = U * A * (T_aw_a - T_aw_b)
  double dT = T_aw_a - T_aw_b;
  result.Q = U * A * dT;

  // Wall temperature: hot-side surface temperature [K]
  if (T_aw_a >= T_aw_b) {
    result.T_wall = T_aw_a - result.Q / (h_a * A);
  } else {
    result.T_wall = T_aw_b - (-result.Q) / (h_b * A);
  }

  double U2 = U * U;
  result.dQ_dh_a = (U2 / (h_a * h_a)) * A * dT;
  result.dQ_dh_b = (U2 / (h_b * h_b)) * A * dT;
  result.dQ_dT_aw_a = U * A;
  result.dQ_dT_aw_b = -U * A;

  return result;
}

WallCouplingResult wall_coupling_and_jacobian(
    double h_a, double T_aw_a,
    double h_b, double T_aw_b,
    const std::vector<double> &t_over_k_layers,
    double A,
    double R_fouling) {
  WallCouplingResult result;

  // Total conductive resistance
  double R_wall = R_fouling;
  for (double tk : t_over_k_layers) {
    R_wall += tk;
  }

  // Overall HTC: U = 1 / (1/h_a + R_wall + 1/h_b)
  double R_total = 1.0 / h_a + R_wall + 1.0 / h_b;
  double U = 1.0 / R_total;

  // Heat transfer rate: Q = U * A * (T_aw_a - T_aw_b)
  double dT = T_aw_a - T_aw_b;
  result.Q = U * A * dT;

  // T_wall on side A (hot-side surface)
  result.T_wall = T_aw_a - (result.Q / (h_a * A));

  double U2 = U * U;
  result.dQ_dh_a = (U2 / (h_a * h_a)) * A * dT;
  result.dQ_dh_b = (U2 / (h_b * h_b)) * A * dT;
  result.dQ_dT_aw_a = U * A;
  result.dQ_dT_aw_b = -U * A;

  return result;
}

} // namespace combaero
