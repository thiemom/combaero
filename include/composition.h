#pragma once

#include <vector>

namespace combaero {

// Normalize a vector of fractions to sum to 1.0.
// Returns all zeros with a warning if input contains all zeros.
std::vector<double> normalize_fractions(const std::vector<double>& fractions);

// Convert mole fractions X_k to mass fractions Y_k.
std::vector<double> mole_to_mass(const std::vector<double>& X);

// Convert mass fractions Y_k to mole fractions X_k.
std::vector<double> mass_to_mole(const std::vector<double>& Y);

// Convert mole fractions to dry fractions (remove water vapor and normalize).
std::vector<double> convert_to_dry_fractions(const std::vector<double>& mole_fractions);

// Calculate mixture molecular weight [g/mol].
double mwmix(const std::vector<double>& X);

// A gradient in the mole fractions, g_X[k] = dg/dX_k, as the gradient in the
// mass fractions at the mixture X (#717). The step keeps sum X = 1:
//   dX_i/dY_k = (delta_ik - X_i) MW / W_k,
// so dg/dY_k = (g_X[k] - sum_i X_i g_X[i]) MW / W_k. Any extension of g off
// the simplex gives the same result (its normal component cancels).
std::vector<double> mole_gradient_to_mass(const std::vector<double>& X,
                                          const std::vector<double>& g_X);

} // namespace combaero
