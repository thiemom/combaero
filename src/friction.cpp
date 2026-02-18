// Darcy friction factor correlations for turbulent pipe flow.
//
// References:
// - Colebrook (1939): Journal of the ICE, 11
// - Haaland (1983): ASME J. Fluids Eng., 105(1)
// - Serghides (1984): Chemical Engineering, Nov. 5

#include "friction.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

// Haaland correlation (1983)
// 1/√f = -1.8 * log10( (ε/D / 3.7)^1.11 + 6.9/Re )
// Accuracy: ~2-3% vs Colebrook
double friction_haaland(double Re, double e_D) {
    if (Re <= 0.0) {
        throw std::invalid_argument("friction_haaland: Re must be positive");
    }
    if (e_D < 0.0) {
        throw std::invalid_argument("friction_haaland: e_D must be non-negative");
    }

    double term = std::pow(e_D / 3.7, 1.11) + 6.9 / Re;
    double inv_sqrt_f = -1.8 * std::log10(term);

    return 1.0 / (inv_sqrt_f * inv_sqrt_f);
}

// Serghides correlation (1984) - Steffensen acceleration on Colebrook
// A = -2 * log10( ε/D/3.7 + 12/Re )
// B = -2 * log10( ε/D/3.7 + 2.51*A/Re )
// C = -2 * log10( ε/D/3.7 + 2.51*B/Re )
// 1/√f = A - (B-A)² / (C - 2B + A)
// Accuracy: <0.3% vs Colebrook
double friction_serghides(double Re, double e_D) {
    if (Re <= 0.0) {
        throw std::invalid_argument("friction_serghides: Re must be positive");
    }
    if (e_D < 0.0) {
        throw std::invalid_argument("friction_serghides: e_D must be non-negative");
    }

    // Note: In Serghides, A/B/C are successive approximations to 1/√f
    // The term 2.51*A/Re corresponds to 2.51/(Re*√f) when √f ≈ 1/A
    double A = -2.0 * std::log10(e_D / 3.7 + 12.0 / Re);
    double B = -2.0 * std::log10(e_D / 3.7 + 2.51 * A / Re);
    double C = -2.0 * std::log10(e_D / 3.7 + 2.51 * B / Re);

    double denom = C - 2.0 * B + A;

    double inv_sqrt_f;
    if (std::abs(denom) < 1e-30) {
        inv_sqrt_f = C;  // Use last iteration if denominator is zero
    } else {
        inv_sqrt_f = A - (B - A) * (B - A) / denom;
    }

    return 1.0 / (inv_sqrt_f * inv_sqrt_f);
}

// Colebrook-White equation (1939) - implicit, solved iteratively
// 1/√f = -2 * log10( ε/D/3.7 + 2.51/(Re*√f) )
// Uses Serghides as initial guess, Newton-Raphson iteration.
double friction_colebrook(double Re, double e_D, double tol, int max_iter) {
    if (Re <= 0.0) {
        throw std::invalid_argument("friction_colebrook: Re must be positive");
    }
    if (e_D < 0.0) {
        throw std::invalid_argument("friction_colebrook: e_D must be non-negative");
    }

    // Initial guess from Serghides (<0.3% error vs Colebrook, better than Haaland ~2-3%)
    double f = friction_serghides(Re, e_D);

    // Newton-Raphson iteration on: F(x) = x + 2*log10(ε/D/3.7 + 2.51*x/Re) = 0
    // where x = 1/√f
    // Note: 2.51/(Re*√f) = 2.51/(Re/x) = 2.51*x/Re
    const double ln10 = std::log(10.0);
    double x = 1.0 / std::sqrt(f);

    for (int iter = 0; iter < max_iter; ++iter) {
        double arg = e_D / 3.7 + 2.51 * x / Re;
        double F = x + 2.0 * std::log10(arg);

        // dF/dx = 1 + 2/(ln10 * arg) * d(arg)/dx
        // d(arg)/dx = 2.51/Re
        double darg_dx = 2.51 / Re;
        double dF = 1.0 + 2.0 * darg_dx / (ln10 * arg);

        if (std::abs(dF) < 1e-30) break;

        double dx = -F / dF;
        x += dx;

        if (x <= 0.0) x = 0.1;

        if (std::abs(dx) < tol * std::abs(x)) {
            break;
        }
    }

    return 1.0 / (x * x);
}

// -------------------------------------------------------------
// Pipe Roughness Database
// -------------------------------------------------------------

static const std::unordered_map<std::string, double> ROUGHNESS_DATA = {
    // Smooth surfaces
    {"smooth",              0.0},
    {"drawn_tubing",        1.5e-6},
    {"pvc",                 1.5e-6},
    {"plastic",             1.5e-6},
    // Steel pipes
    {"commercial_steel",    4.5e-5},
    {"new_steel",           4.5e-5},
    {"wrought_iron",        4.5e-5},
    {"galvanized_iron",     1.5e-4},
    {"galvanized_steel",    1.5e-4},
    {"rusted_steel",        2.5e-4},
    // Cast iron
    {"cast_iron",           2.6e-4},
    {"asphalted_cast_iron", 1.2e-4},
    // Concrete
    {"concrete",            3.0e-4},
    {"rough_concrete",      3.0e-3},
    // Other materials
    {"riveted_steel",       9.0e-4},
    {"wood_stave",          1.8e-4},
    {"corrugated_metal",    4.5e-2},
};

std::unordered_map<std::string, double> standard_pipe_roughness() {
    return ROUGHNESS_DATA;
}

double pipe_roughness(const std::string& material) {
    std::string key = material;
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    auto it = ROUGHNESS_DATA.find(key);
    if (it != ROUGHNESS_DATA.end()) {
        return it->second;
    }
    throw std::invalid_argument(
        "pipe_roughness: unknown material '" + material + "'. "
        "Use standard_pipe_roughness() to see available materials.");
}

// Petukhov correlation (1970) - smooth pipes only
// f = (0.790 * ln(Re) - 1.64)^(-2)
// Valid for 3000 < Re < 5x10^6
double friction_petukhov(double Re) {
    if (Re < 3000) {
        throw std::invalid_argument("friction_petukhov: Re must be > 3000");
    }
    double x = 0.790 * std::log(Re) - 1.64;
    return 1.0 / (x * x);
}
