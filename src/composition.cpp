#include "../include/composition.h"
#include "../include/thermo_transport_data.h" // For molar_masses
#include "../include/thermo.h"               // For species_index_from_name
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>

namespace combaero {

std::vector<double> normalize_fractions(const std::vector<double>& fractions) {
    std::vector<double> clamped = fractions;
    for (double& value : clamped) {
        if (value < 1.0e-15) value = 1.0e-15;
    }

    double sum = std::accumulate(clamped.begin(), clamped.end(), 0.0);
    std::vector<double> normalized = clamped;

    if (std::abs(sum) < 1.0e-10) {
        std::cerr << "Warning: normalize_fractions received all zeros. Returning all zeros." << std::endl;
        return fractions;
    }

    for (double& value : normalized) {
        value /= sum;
    }

    return normalized;
}

std::vector<double> mole_to_mass(const std::vector<double>& X) {
    if (X.size() != molar_masses.size()) {
        throw std::invalid_argument("mole_to_mass: size mismatch");
    }

    double denom = 0.0;
    for (std::size_t k = 0; k < X.size(); ++k) {
        denom += X[k] * molar_masses[k];
    }
    if (denom <= 0.0) {
        throw std::runtime_error("mole_to_mass: non-positive denominator");
    }

    std::vector<double> Y(X.size());
    for (std::size_t k = 0; k < X.size(); ++k) {
        Y[k] = X[k] * molar_masses[k] / denom;
    }
    return Y;
}

std::vector<double> mass_to_mole(const std::vector<double>& Y) {
    if (Y.size() != molar_masses.size()) {
        throw std::invalid_argument("mass_to_mole: size mismatch");
    }

    double denom = 0.0;
    for (std::size_t k = 0; k < Y.size(); ++k) {
        double Wk = molar_masses[k];
        if (Wk <= 0.0) {
            throw std::runtime_error("mass_to_mole: non-positive molar mass");
        }
        denom += Y[k] / Wk;
    }
    if (denom <= 0.0) {
        throw std::runtime_error("mass_to_mole: non-positive denominator");
    }

    std::vector<double> X(Y.size());
    for (std::size_t k = 0; k < Y.size(); ++k) {
        X[k] = (Y[k] / molar_masses[k]) / denom;
    }
    return X;
}

std::vector<double> convert_to_dry_fractions(const std::vector<double>& mole_fractions) {
    std::size_t h2o_idx = species_index_from_name("H2O");
    std::vector<double> dry_fractions = mole_fractions;

    double sum = 0.0;
    for (std::size_t i = 0; i < dry_fractions.size(); ++i) {
        if (i != h2o_idx) {
            sum += dry_fractions[i];
        }
    }

    if (std::abs(sum) < 1.0e-10) {
        std::cerr << "Warning: convert_to_dry_fractions received only water vapor. Returning all zeros." << std::endl;
        std::fill(dry_fractions.begin(), dry_fractions.end(), 0.0);
        return dry_fractions;
    }

    dry_fractions[h2o_idx] = 0.0;

    for (double& value : dry_fractions) {
        value /= sum;
    }

    return dry_fractions;
}

double mwmix(const std::vector<double>& X) {
    if (X.size() != molar_masses.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double sum = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        sum += X[i] * molar_masses[i];
    }
    return sum;
}

} // namespace combaero
