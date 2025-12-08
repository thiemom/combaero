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
