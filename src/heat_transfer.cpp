#include "../include/heat_transfer.h"
#include "../include/state.h"
#include "../include/friction.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

// -------------------------------------------------------------
// Friction factor helpers
// -------------------------------------------------------------

double friction_petukhov(double Re) {
    if (Re < 3000) {
        throw std::invalid_argument("friction_petukhov: Re must be > 3000");
    }
    double x = 0.790 * std::log(Re) - 1.64;
    return 1.0 / (x * x);
}

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
