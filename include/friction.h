#ifndef FRICTION_H
#define FRICTION_H

// -------------------------------------------------------------
// Darcy friction factor correlations for turbulent pipe flow
// -------------------------------------------------------------
// All functions compute the Darcy friction factor f for use in:
//   ΔP = f * (L/D) * (ρ * v²/2)
//
// Parameters:
//   Re  : Reynolds number [-] (must be > 0)
//   e_D : relative roughness ε/D [-] (must be >= 0)
//
// Valid for turbulent flow (Re > ~2300). For laminar flow, use f = 64/Re.
//
// References:
// - Colebrook (1939): Journal of the ICE, 11
// - Haaland (1983): ASME J. Fluids Eng., 105(1)
// - Serghides (1984): Chemical Engineering, Nov. 5

// Haaland correlation (1983) - explicit approximation
// 1/√f = -1.8 * log10( (ε/D / 3.7)^1.11 + 6.9/Re )
// Accuracy: ~2-3% vs Colebrook-White
double friction_haaland(double Re, double e_D);

// Serghides correlation (1984) - explicit approximation
// Uses Steffensen acceleration on Colebrook iteration.
// Accuracy: <0.3% vs Colebrook-White (very accurate)
double friction_serghides(double Re, double e_D);

// Colebrook-White equation (1939) - implicit, solved iteratively
// 1/√f = -2 * log10( ε/D/3.7 + 2.51/(Re*√f) )
// The reference standard for turbulent friction factor.
// Uses Haaland as initial guess, Newton-Raphson iteration.
double friction_colebrook(double Re, double e_D, double tol = 1e-10, int max_iter = 20);

// Petukhov correlation (1970) - smooth pipes only
// f = (0.790 * ln(Re) - 1.64)^(-2)
// Valid for 3000 < Re < 5×10^6
// Often used with Gnielinski/Petukhov heat transfer correlations.
// Reference: Petukhov (1970), Advances in Heat Transfer, 6, 503
double friction_petukhov(double Re);

#endif // FRICTION_H
