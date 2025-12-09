#include "../include/heat_transfer.h"
#include "../include/state.h"
#include "../include/friction.h"  // friction_petukhov, friction_colebrook
#include <cmath>
#include <stdexcept>
#include <algorithm>

// -------------------------------------------------------------
// Internal Flow Correlations
// -------------------------------------------------------------

double nusselt_dittus_boelter(double Re, double Pr, bool heating) {
    if (Re < 10000) {
        throw std::invalid_argument("nusselt_dittus_boelter: Re must be > 10000 for turbulent flow");
    }
    if (Pr < 0.6 || Pr > 160) {
        throw std::invalid_argument("nusselt_dittus_boelter: Pr must be in range [0.6, 160]");
    }
    
    double n = heating ? 0.4 : 0.3;
    return 0.023 * std::pow(Re, 0.8) * std::pow(Pr, n);
}

double nusselt_gnielinski(double Re, double Pr, double f) {
    if (Re < 2300 || Re > 5e6) {
        throw std::invalid_argument("nusselt_gnielinski: Re must be in range [2300, 5e6]");
    }
    if (Pr < 0.5 || Pr > 2000) {
        throw std::invalid_argument("nusselt_gnielinski: Pr must be in range [0.5, 2000]");
    }
    if (f <= 0) {
        throw std::invalid_argument("nusselt_gnielinski: friction factor must be positive");
    }
    
    double f8 = f / 8.0;
    double sqrt_f8 = std::sqrt(f8);
    double Pr_23 = std::pow(Pr, 2.0/3.0);
    
    return f8 * (Re - 1000.0) * Pr / (1.0 + 12.7 * sqrt_f8 * (Pr_23 - 1.0));
}

double nusselt_gnielinski(double Re, double Pr) {
    // Use Petukhov friction factor for smooth pipe
    double f = friction_petukhov(std::max(Re, 3000.0));
    return nusselt_gnielinski(Re, Pr, f);
}

double nusselt_sieder_tate(double Re, double Pr, double mu_ratio) {
    if (Re < 10000) {
        throw std::invalid_argument("nusselt_sieder_tate: Re must be > 10000 for turbulent flow");
    }
    if (Pr < 0.7 || Pr > 16700) {
        throw std::invalid_argument("nusselt_sieder_tate: Pr must be in range [0.7, 16700]");
    }
    if (mu_ratio <= 0) {
        throw std::invalid_argument("nusselt_sieder_tate: mu_ratio must be positive");
    }
    
    return 0.027 * std::pow(Re, 0.8) * std::pow(Pr, 1.0/3.0) * std::pow(mu_ratio, 0.14);
}

double nusselt_petukhov(double Re, double Pr, double f) {
    if (Re < 1e4 || Re > 5e6) {
        throw std::invalid_argument("nusselt_petukhov: Re must be in range [1e4, 5e6]");
    }
    if (Pr < 0.5 || Pr > 2000) {
        throw std::invalid_argument("nusselt_petukhov: Pr must be in range [0.5, 2000]");
    }
    if (f <= 0) {
        throw std::invalid_argument("nusselt_petukhov: friction factor must be positive");
    }
    
    double f8 = f / 8.0;
    double sqrt_f8 = std::sqrt(f8);
    double Pr_23 = std::pow(Pr, 2.0/3.0);
    
    return f8 * Re * Pr / (1.07 + 12.7 * sqrt_f8 * (Pr_23 - 1.0));
}

double nusselt_petukhov(double Re, double Pr) {
    double f = friction_petukhov(Re);
    return nusselt_petukhov(Re, Pr, f);
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
