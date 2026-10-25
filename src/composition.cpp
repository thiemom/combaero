#include "../include/composition.h"
#include "../include/thermo_transport_data.h" // For molar_masses
#include "../include/thermo.h"               // For species_index_from_name
#include "../include/correlation_status.h" // For warn()
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace combaero {

std::vector<double> normalize_fractions(const std::vector<double>& fractions) {
    // Negatives (a solver iterate's round-off) clamp to 0; an absent species
    // stays exactly 0. The former 1e-15 floor invented a trace of every
    // species in every mixed stream, which every mixture sum then evaluated
    // (dry air: 15 species instead of 4), for no protection: the logs of a
    // fraction are guarded where they are taken (s(), the equilibrium Q's).
    std::vector<double> clamped = fractions;
    for (double& value : clamped) {
        if (value < 0.0) value = 0.0;
    }

    double sum = std::accumulate(clamped.begin(), clamped.end(), 0.0);
    std::vector<double> normalized = clamped;

    if (std::abs(sum) < 1.0e-10) {
        warn("normalize_fractions received all zeros. Returning all zeros.");
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
        warn("convert_to_dry_fractions received only water vapor. Returning all zeros.");
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

std::vector<double> mole_gradient_to_mass(const std::vector<double>& X,
                                          const std::vector<double>& g_X) {
    if (X.size() != molar_masses.size() || g_X.size() != X.size()) {
        throw std::invalid_argument("mole_gradient_to_mass: vector size does not match number of species");
    }
    const double mw = mwmix(X);
    double x_dot = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        x_dot += X[i] * g_X[i];
    }
    std::vector<double> out(X.size(), 0.0);
    for (std::size_t k = 0; k < X.size(); ++k) {
        out[k] = (g_X[k] - x_dot) * mw / molar_masses[k];
    }
    return out;
}

} // namespace combaero
