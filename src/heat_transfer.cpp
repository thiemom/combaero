#include "../include/heat_transfer.h"
#include "../include/correlation_status.h"
#include "../include/state.h"
#include "../include/friction.h"           // friction_petukhov, friction_colebrook
#include "../include/thermo.h"             // density, cp_mass, speed_of_sound
#include "../include/transport.h"          // viscosity, thermal_conductivity, prandtl
#include "../include/stagnation.h"         // T_adiabatic_wall
#include "../include/cooling_correlations.h" // rib_*, dimple_*, pin_fin_*, impingement_*
#include "../include/math_constants.h"     // M_PI cross-platform
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>

// -------------------------------------------------------------
// Internal Flow Correlations
// -------------------------------------------------------------

double nusselt_dittus_boelter(double Re, double Pr, bool heating,
                              combaero::CorrelationStatus* status) {
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
    double n = heating ? 0.4 : 0.3;
    return 0.023 * std::pow(Re, 0.8) * std::pow(Pr, n);
}

double nusselt_gnielinski(double Re, double Pr, double f,
                          combaero::CorrelationStatus* status) {
    if (f <= 0.0) {
        throw std::invalid_argument("nusselt_gnielinski: friction factor must be positive");
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
    double Pr_23 = std::pow(Pr, 2.0/3.0);
    double denom = 1.0 + 12.7 * sqrt_f8 * (Pr_23 - 1.0);

    if (Re >= 2300.0) {
        return f8 * (Re - 1000.0) * Pr / denom;
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
    constexpr double Re0 = 1000.0;   // lower anchor: laminar Nu
    constexpr double Re1 = 2300.0;   // upper anchor: Gnielinski
    const double Nu0 = NU_LAMINAR_CONST_T;  // 3.66
    const double dNu0 = 0.0;                // laminar: flat w.r.t. Re
    const double Nu1 = f8 * (Re1 - 1000.0) * Pr / denom;
    const double dNu1 = f8 * Pr / denom;   // d(Nu)/d(Re) at Re1

    // Normalised coordinate t in [0, 1]
    double h = Re1 - Re0;
    double t = (Re - Re0) / h;
    t = std::max(0.0, std::min(1.0, t));  // clamp outside [Re0, Re1]

    // Cubic Hermite basis
    double t2 = t * t;
    double t3 = t2 * t;
    double h00 =  2.0*t3 - 3.0*t2 + 1.0;
    double h10 =      t3 - 2.0*t2 + t;
    double h01 = -2.0*t3 + 3.0*t2;
    double h11 =      t3 -     t2;

    return h00*Nu0 + h10*h*dNu0 + h01*Nu1 + h11*h*dNu1;
}

double nusselt_gnielinski(double Re, double Pr,
                          combaero::CorrelationStatus* status) {
    double f = friction_petukhov(std::max(Re, 3000.0));  // avoid log(0) in Petukhov
    return nusselt_gnielinski(Re, Pr, f, status);
}

double nusselt_sieder_tate(double Re, double Pr, double mu_ratio,
                           combaero::CorrelationStatus* status) {
    if (mu_ratio <= 0.0) {
        throw std::invalid_argument("nusselt_sieder_tate: mu_ratio must be positive");
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
    return 0.027 * std::pow(Re, 0.8) * std::pow(Pr, 1.0/3.0) * std::pow(mu_ratio, 0.14);
}

double nusselt_petukhov(double Re, double Pr, double f,
                        combaero::CorrelationStatus* status) {
    if (f <= 0.0) {
        throw std::invalid_argument("nusselt_petukhov: friction factor must be positive");
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
    double Pr_23 = std::pow(Pr, 2.0/3.0);
    return f8 * Re * Pr / (1.07 + 12.7 * sqrt_f8 * (Pr_23 - 1.0));
}

double nusselt_petukhov(double Re, double Pr,
                        combaero::CorrelationStatus* status) {
    double f = friction_petukhov(std::max(Re, 3000.0));  // avoid log(0) in Petukhov
    return nusselt_petukhov(Re, Pr, f, status);
}

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

double htc_from_nusselt(double Nu, double k, double L) {
    if (L <= 0) {
        throw std::invalid_argument("htc_from_nusselt: characteristic length must be positive");
    }
    if (k <= 0) {
        throw std::invalid_argument("htc_from_nusselt: thermal conductivity must be positive");
    }
    return Nu * k / L;
}

// -------------------------------------------------------------
// Overall Heat Transfer Coefficient
// -------------------------------------------------------------

double overall_htc(const std::vector<double>& h_values,
                   const std::vector<double>& t_over_k) {
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
            throw std::invalid_argument("overall_htc: t/k values must be non-negative");
        }
        R_total += tk;
    }

    return 1.0 / R_total;
}

double overall_htc_wall(double h_inner, double h_outer,
                        double t_wall, double k_wall) {
    if (k_wall <= 0) {
        throw std::invalid_argument("overall_htc_wall: k_wall must be positive");
    }
    return overall_htc({h_inner, h_outer}, {t_wall / k_wall});
}

double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double>& t_over_k_layers) {
    return overall_htc({h_inner, h_outer}, t_over_k_layers);
}

double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double>& t_over_k_layers,
                        double R_fouling) {
    if (R_fouling < 0) {
        throw std::invalid_argument("overall_htc_wall: R_fouling must be non-negative");
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
        throw std::invalid_argument("lmtd: temperature differences must be positive");
    }

    // If dT1 ≈ dT2, use arithmetic mean to avoid 0/0
    double ratio = dT1 / dT2;
    if (std::abs(ratio - 1.0) < 1e-6) {
        return (dT1 + dT2) / 2.0;
    }

    return (dT1 - dT2) / std::log(ratio);
}

double lmtd_counterflow(double T_hot_in, double T_hot_out,
                        double T_cold_in, double T_cold_out) {
    double dT1 = T_hot_in - T_cold_out;
    double dT2 = T_hot_out - T_cold_in;
    return lmtd(dT1, dT2);
}

double lmtd_parallelflow(double T_hot_in, double T_hot_out,
                         double T_cold_in, double T_cold_out) {
    double dT1 = T_hot_in - T_cold_in;
    double dT2 = T_hot_out - T_cold_out;
    return lmtd(dT1, dT2);
}

// -------------------------------------------------------------
// Temperature Profile Through Wall
// -------------------------------------------------------------

std::vector<double> wall_temperature_profile(
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k,
    double& q) {

    if (h_hot <= 0 || h_cold <= 0) {
        throw std::invalid_argument("wall_temperature_profile: HTCs must be positive");
    }

    // Total thermal resistance per unit area [m²·K/W]
    double R_total = 1.0 / h_hot + 1.0 / h_cold;
    for (double r : t_over_k) {
        if (r < 0) {
            throw std::invalid_argument("wall_temperature_profile: t/k values must be non-negative");
        }
        R_total += r;
    }

    // Heat flux (constant through all layers)
    q = (T_hot - T_cold) / R_total;

    // Build temperature profile
    // Start from hot side, subtract temperature drops
    std::vector<double> temps;
    temps.reserve(t_over_k.size() + 2);

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

std::vector<double> wall_temperature_profile(
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k) {
    double q_unused;
    return wall_temperature_profile(T_hot, T_cold, h_hot, h_cold, t_over_k, q_unused);
}

// -------------------------------------------------------------
// Heat Flux from Measured Temperature
// -------------------------------------------------------------

double heat_flux_from_T_at_edge(
    double T_measured, std::size_t edge_idx,
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k) {

    std::size_t n_layers = t_over_k.size();
    if (edge_idx > n_layers) {
        throw std::invalid_argument(
            "heat_flux_from_T_at_edge: edge_idx out of range (max = " +
            std::to_string(n_layers) + ")");
    }

    if (h_hot <= 0 || h_cold <= 0) {
        throw std::invalid_argument("heat_flux_from_T_at_edge: HTCs must be positive");
    }

    // Thermal resistance from hot bulk to edge
    double R_hot_to_edge = 1.0 / h_hot;  // convective resistance
    for (std::size_t i = 0; i < edge_idx; ++i) {
        R_hot_to_edge += t_over_k[i];
    }

    // Thermal resistance from edge to cold bulk
    double R_edge_to_cold = 0.0;
    for (std::size_t i = edge_idx; i < n_layers; ++i) {
        R_edge_to_cold += t_over_k[i];
    }
    R_edge_to_cold += 1.0 / h_cold;  // convective resistance

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

double heat_flux_from_T_at_depth(
    double T_measured, double depth_from_hot,
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& thicknesses,
    const std::vector<double>& conductivities) {

    if (thicknesses.size() != conductivities.size()) {
        throw std::invalid_argument(
            "heat_flux_from_T_at_depth: thicknesses and conductivities must have same size");
    }
    if (depth_from_hot < 0) {
        throw std::invalid_argument("heat_flux_from_T_at_depth: depth must be non-negative");
    }
    if (h_hot <= 0 || h_cold <= 0) {
        throw std::invalid_argument("heat_flux_from_T_at_depth: HTCs must be positive");
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
    double R_hot_to_point = 1.0 / h_hot;  // convective resistance
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

double bulk_T_from_edge_T_and_q(
    double T_measured, std::size_t edge_idx, double q,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k,
    const std::string& solve_for) {

    std::size_t n_layers = t_over_k.size();
    if (edge_idx > n_layers) {
        throw std::invalid_argument(
            "bulk_T_from_edge_T_and_q: edge_idx out of range (max = " +
            std::to_string(n_layers) + ")");
    }

    if (h_hot <= 0 || h_cold <= 0) {
        throw std::invalid_argument("bulk_T_from_edge_T_and_q: HTCs must be positive");
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

std::pair<double, double> dT_edge_dT_bulk(
    std::size_t edge_idx,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k) {

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
    // T_edge = T_hot * (R_edge_to_cold / R_total) + T_cold * (R_hot_to_edge / R_total)
    double dT_dT_hot = R_edge_to_cold / R_total;
    double dT_dT_cold = R_hot_to_edge / R_total;

    return {dT_dT_hot, dT_dT_cold};
}

double dT_edge_dT_hot(
    std::size_t edge_idx,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k) {
    return dT_edge_dT_bulk(edge_idx, h_hot, h_cold, t_over_k).first;
}

double dT_edge_dT_cold(
    std::size_t edge_idx,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k) {
    return dT_edge_dT_bulk(edge_idx, h_hot, h_cold, t_over_k).second;
}

double dT_edge_dq(
    std::size_t edge_idx,
    double h_hot,
    const std::vector<double>& t_over_k) {

    std::size_t n_layers = t_over_k.size();
    if (edge_idx > n_layers) {
        throw std::invalid_argument(
            "dT_edge_dq: edge_idx out of range (max = " +
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
        throw std::invalid_argument("effectiveness_counterflow: NTU must be non-negative");
    }
    if (C_r < 0 || C_r > 1) {
        throw std::invalid_argument("effectiveness_counterflow: C_r must be in [0, 1]");
    }

    // Special case: C_r = 1 (balanced heat exchanger)
    if (std::abs(C_r - 1.0) < 1e-10) {
        return NTU / (1.0 + NTU);
    }

    // Special case: C_r = 0 (one fluid has infinite capacity, e.g., condenser/evaporator)
    if (C_r < 1e-10) {
        return 1.0 - std::exp(-NTU);
    }

    // General case
    double exp_term = std::exp(-NTU * (1.0 - C_r));
    return (1.0 - exp_term) / (1.0 - C_r * exp_term);
}

double effectiveness_parallelflow(double NTU, double C_r) {
    if (NTU < 0) {
        throw std::invalid_argument("effectiveness_parallelflow: NTU must be non-negative");
    }
    if (C_r < 0 || C_r > 1) {
        throw std::invalid_argument("effectiveness_parallelflow: C_r must be in [0, 1]");
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

double nusselt_pipe(const State& s, double velocity, double diameter,
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
    // Turbulent (Re > 10000): Gnielinski is still good, or could use Dittus-Boelter

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

double htc_pipe(const State& s, double velocity, double diameter,
                bool heating, double roughness) {
    double Nu = nusselt_pipe(s, velocity, diameter, heating, roughness);
    double k = s.k();
    return htc_from_nusselt(Nu, k, diameter);
}

// Composite function: compute HTC from thermodynamic state
// Returns: (h [W/(m²·K)], Nu [-], Re [-])
std::tuple<double, double, double> htc_pipe(
    double T, double P, const std::vector<double>& X,
    double velocity, double diameter,
    const std::string& correlation,
    bool heating,
    double mu_ratio,
    double roughness)
{
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
    double cp = cp_mass(T, X);
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
    }
    else if (correlation == "dittus_boelter") {
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
    }
    else if (correlation == "sieder_tate") {
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
    }
    else if (correlation == "petukhov") {
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
    }
    else {
        throw std::invalid_argument(
            "htc_pipe: unknown correlation '" + correlation + "'. " +
            "Valid options: 'gnielinski', 'dittus_boelter', 'sieder_tate', 'petukhov'");
    }

    // Heat transfer coefficient
    double h = htc_from_nusselt(Nu, k, diameter);

    return std::make_tuple(h, Nu, Re);
}

// -------------------------------------------------------------
// Internal helper: build ChannelResult from common fields
// -------------------------------------------------------------

static ChannelResult make_channel_result(
    double h, double Nu, double Re, double Pr, double f, double dP,
    double M, double T_aw, double T_wall)
{
    double q = std::numeric_limits<double>::quiet_NaN();
    if (std::isfinite(T_wall)) {
        q = h * (T_aw - T_wall);
    }
    return ChannelResult{h, Nu, Re, Pr, f, dP, M, T_aw, q};
}

// -------------------------------------------------------------
// channel_smooth
// -------------------------------------------------------------

ChannelResult channel_smooth(
    double T, double P, const std::vector<double>& X,
    double velocity, double diameter, double length,
    double T_wall,
    const std::string& correlation,
    bool heating,
    double mu_ratio,
    double roughness)
{
    if (velocity < 0.0) {
        throw std::invalid_argument("channel_smooth: velocity must be non-negative");
    }
    if (diameter <= 0.0) {
        throw std::invalid_argument("channel_smooth: diameter must be positive");
    }
    if (length < 0.0) {
        throw std::invalid_argument("channel_smooth: length must be non-negative");
    }

    // Fluid properties — evaluated once
    double rho = density(T, P, X);
    double mu  = viscosity(T, P, X);
    double k   = thermal_conductivity(T, P, X);
    double Pr  = prandtl(T, P, X);

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
        throw std::invalid_argument(
            "channel_smooth: unknown correlation '" + correlation + "'");
    }

    double h  = htc_from_nusselt(Nu, k, diameter);
    double dP = (velocity > 0.0) ? f * (length / diameter) * (rho * velocity * velocity / 2.0) : 0.0;

    // Mach number and T_aw — always continuous
    double a = speed_of_sound(T, X);
    double M = (a > 0.0) ? velocity / a : 0.0;
    const bool turbulent_flow = true;
    double T_aw = T_adiabatic_wall(T, velocity, T, P, X, turbulent_flow);

    return make_channel_result(h, Nu, Re, Pr, f, dP, M, T_aw, T_wall);
}

// -------------------------------------------------------------
// channel_ribbed
// -------------------------------------------------------------

ChannelResult channel_ribbed(
    double T, double P, const std::vector<double>& X,
    double velocity, double diameter, double length,
    double e_D, double P_e, double alpha,
    double T_wall,
    bool heating)
{
    // Get smooth-pipe baseline (Gnielinski, no roughness)
    ChannelResult base = channel_smooth(T, P, X, velocity, diameter, length,
                                        T_wall, "gnielinski", heating);

    // Apply rib multipliers
    double enh = combaero::cooling::rib_enhancement_factor(e_D, P_e, alpha);
    double fmul = combaero::cooling::rib_friction_multiplier(e_D, P_e);

    double Nu_rib = base.Nu * enh;
    double f_rib  = base.f  * fmul;

    double rho = density(T, P, X);
    double k   = thermal_conductivity(T, P, X);
    double h_rib = htc_from_nusselt(Nu_rib, k, diameter);
    double dP_rib = (velocity > 0.0) ? f_rib * (length / diameter) * (rho * velocity * velocity / 2.0) : 0.0;

    return make_channel_result(h_rib, Nu_rib, base.Re, base.Pr, f_rib, dP_rib,
                               base.M, base.T_aw, T_wall);
}

// -------------------------------------------------------------
// channel_dimpled
// -------------------------------------------------------------

ChannelResult channel_dimpled(
    double T, double P, const std::vector<double>& X,
    double velocity, double diameter, double length,
    double d_Dh, double h_d, double S_d,
    double T_wall,
    bool heating)
{
    // Get smooth-pipe baseline
    ChannelResult base = channel_smooth(T, P, X, velocity, diameter, length,
                                        T_wall, "gnielinski", heating);

    // Apply dimple multipliers
    double enh  = combaero::cooling::dimple_nusselt_enhancement(base.Re, d_Dh, h_d, S_d);
    double fmul = combaero::cooling::dimple_friction_multiplier(base.Re, d_Dh, h_d);

    double Nu_dim = base.Nu * enh;
    double f_dim  = base.f  * fmul;

    double rho = density(T, P, X);
    double k   = thermal_conductivity(T, P, X);
    double h_dim = htc_from_nusselt(Nu_dim, k, diameter);
    double dP_dim = (velocity > 0.0) ? f_dim * (length / diameter) * (rho * velocity * velocity / 2.0) : 0.0;

    return make_channel_result(h_dim, Nu_dim, base.Re, base.Pr, f_dim, dP_dim,
                               base.M, base.T_aw, T_wall);
}

// -------------------------------------------------------------
// channel_pin_fin
// -------------------------------------------------------------

ChannelResult channel_pin_fin(
    double T, double P, const std::vector<double>& X,
    double velocity, double channel_height, double pin_diameter,
    double S_D, double X_D, int N_rows,
    double T_wall,
    bool is_staggered)
{
    if (velocity < 0.0) {
        throw std::invalid_argument("channel_pin_fin: velocity must be non-negative");
    }
    if (pin_diameter <= 0.0) {
        throw std::invalid_argument("channel_pin_fin: pin_diameter must be positive");
    }
    if (channel_height <= 0.0) {
        throw std::invalid_argument("channel_pin_fin: channel_height must be positive");
    }
    if (N_rows < 1) {
        throw std::invalid_argument("channel_pin_fin: N_rows must be >= 1");
    }
    if (S_D <= 1.0) {
        throw std::invalid_argument("channel_pin_fin: S_D must be > 1 (minimum cross-section constraint)");
    }

    // Fluid properties
    double rho = density(T, P, X);
    double mu  = viscosity(T, P, X);
    double k   = thermal_conductivity(T, P, X);
    double Pr  = prandtl(T, P, X);

    // Re based on pin diameter and approach velocity
    double Re_d = (velocity > 0.0) ? rho * velocity * pin_diameter / mu : 0.0;

    // L_D = channel_height / pin_diameter
    double L_D = channel_height / pin_diameter;

    // Nu from Metzger correlation
    double Nu = combaero::cooling::pin_fin_nusselt(Re_d, Pr, L_D, S_D, X_D, is_staggered);
    double h  = htc_from_nusselt(Nu, k, pin_diameter);

    // Pressure drop: dP = N_rows * f_pin * (rho * v_max^2 / 2)
    // v_max at minimum cross-section: v_max = v * S_D / (S_D - 1)
    double v_max = (S_D > 1.0) ? velocity * S_D / (S_D - 1.0) : velocity;
    double f_pin = (Re_d > 0.0) ? combaero::cooling::pin_fin_friction(Re_d, is_staggered) : 0.0;
    double dP    = static_cast<double>(N_rows) * f_pin * (rho * v_max * v_max / 2.0);

    // Mach and T_aw based on approach velocity
    double a   = speed_of_sound(T, X);
    double M   = (a > 0.0) ? velocity / a : 0.0;
    const bool turbulent_flow = true;
    double T_aw = T_adiabatic_wall(T, velocity, T, P, X, turbulent_flow);

    return make_channel_result(h, Nu, Re_d, Pr, f_pin, dP, M, T_aw, T_wall);
}

// -------------------------------------------------------------
// channel_impingement
// -------------------------------------------------------------

ChannelResult channel_impingement(
    double T, double P, const std::vector<double>& X,
    double mdot_jet, double d_jet,
    double z_D, double x_D, double y_D,
    double A_target,
    double T_wall,
    double Cd_jet)
{
    if (mdot_jet < 0.0) {
        throw std::invalid_argument("channel_impingement: mdot_jet must be non-negative");
    }
    if (d_jet <= 0.0) {
        throw std::invalid_argument("channel_impingement: d_jet must be positive");
    }
    if (A_target <= 0.0) {
        throw std::invalid_argument("channel_impingement: A_target must be positive");
    }
    if (Cd_jet <= 0.0 || Cd_jet > 1.0) {
        throw std::invalid_argument("channel_impingement: Cd_jet must be in (0, 1]");
    }

    // Fluid properties
    double rho = density(T, P, X);
    double mu  = viscosity(T, P, X);
    double k   = thermal_conductivity(T, P, X);
    double Pr  = prandtl(T, P, X);

    // Jet area and velocity
    double A_jet  = M_PI / 4.0 * d_jet * d_jet;
    double v_jet  = (rho > 0.0 && A_jet > 0.0) ? mdot_jet / (rho * A_jet) : 0.0;
    double Re_jet = (mu > 0.0) ? rho * v_jet * d_jet / mu : 0.0;

    // Nu from Florschuetz / Martin correlation
    double Nu = combaero::cooling::impingement_nusselt(Re_jet, Pr, z_D, x_D, y_D);

    // h averaged over target area: h = Nu * k / d_jet
    // (Nu is already an area-averaged value from the correlation)
    double h = htc_from_nusselt(Nu, k, d_jet);

    // Pressure drop across jet plate: dP = (v_jet/Cd)^2 * rho/2
    // Equivalent loss coefficient f = 1/Cd^2
    double f  = 1.0 / (Cd_jet * Cd_jet);
    double dP = f * rho * v_jet * v_jet / 2.0;

    // Mach and T_aw based on jet velocity
    double a   = speed_of_sound(T, X);
    double M   = (a > 0.0) ? v_jet / a : 0.0;
    const bool turbulent_flow = true;
    double T_aw = T_adiabatic_wall(T, v_jet, T, P, X, turbulent_flow);

    return make_channel_result(h, Nu, Re_jet, Pr, f, dP, M, T_aw, T_wall);
}
