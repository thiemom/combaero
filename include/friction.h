#ifndef FRICTION_H
#define FRICTION_H

#include <string>
#include <unordered_map>

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
// Valid for 3000 < Re < 5x10^6
// Often used with Gnielinski/Petukhov heat transfer correlations.
// Reference: Petukhov (1970), Advances in Heat Transfer, 6, 503
double friction_petukhov(double Re);

// -------------------------------------------------------------
// Pipe Roughness Database
// -------------------------------------------------------------
// Absolute roughness ε [m] for common pipe materials.
// Used as e_D = ε/D in friction factor correlations.
//
// References:
// - Moody (1944): Transactions of the ASME, 66(8)
// - White (2011): Fluid Mechanics (7th ed.), Table 6.1
// - Munson et al. (2013): Fundamentals of Fluid Mechanics (7th ed.), Table 8.1
// - Crane Co. (2009): Technical Paper No. 410, Table A-24

// Get absolute roughness for a standard pipe material [m].
// Material names are case-insensitive.
// Throws std::invalid_argument if material not found.
// Use standard_pipe_roughness() to list available materials.
double pipe_roughness(const std::string& material);

// Get all standard pipe roughness values.
// Returns: map of material name -> roughness [m]
std::unordered_map<std::string, double> standard_pipe_roughness();

#endif // FRICTION_H
