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

} // namespace combaero
