// Darcy friction factor correlations for turbulent pipe flow.
//
// References:
// - Colebrook (1939): Journal of the ICE, 11
// - Haaland (1983): ASME J. Fluids Eng., 105(1)
// - Serghides (1984): Chemical Engineering, Nov. 5

#include "friction.h"
#include <cmath>
#include <stdexcept>

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
// Uses Haaland as initial guess, Newton-Raphson iteration.
double friction_colebrook(double Re, double e_D, double tol, int max_iter) {
    if (Re <= 0.0) {
        throw std::invalid_argument("friction_colebrook: Re must be positive");
    }
    if (e_D < 0.0) {
        throw std::invalid_argument("friction_colebrook: e_D must be non-negative");
    }
    
    // Initial guess from Haaland
    double f = friction_haaland(Re, e_D);
    
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
