#include "../include/combustion.h"
#include "../include/humidair.h"
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>


// Calculate oxygen required for complete combustion [mol O2/mol fuel]
double oxygen_required_per_mol_fuel(std::size_t fuel_index) {
    // Validate fuel index
    if (fuel_index >= molecular_structures.size()) {
        throw std::runtime_error("Invalid fuel index");
    }

    // Get molecular structure
    const Molecular_Structure& fuel = molecular_structures[fuel_index];

    // For complete combustion:
    // C_x H_y O_z N_w + (x + y/4 - z/2) O2 -> x CO2 + (y/2) H2O + (w/2) N2
    // Oxygen required = x + y/4 - z/2 [mol O2/mol fuel]

    double oxygen_required = static_cast<double>(fuel.C) + static_cast<double>(fuel.H) / 4.0 - static_cast<double>(fuel.O) / 2.0;

    // Ensure non-negative value (for cases where fuel already has excess oxygen)
    return std::max(0.0, oxygen_required);
}

// Calculate oxygen required for complete combustion [kg O2/kg fuel]
double oxygen_required_per_kg_fuel(std::size_t fuel_index) {
    // Validate fuel index
    if (fuel_index >= molar_masses.size()) {
        throw std::runtime_error("Invalid fuel index");
    }

    // Get molar masses
    double fuel_molar_mass = molar_masses[fuel_index];  // g/mol
    double oxygen_molar_mass = molar_masses[species_index.at("O2")];  // g/mol

    // Calculate molar oxygen requirement
    double molar_oxygen_required = oxygen_required_per_mol_fuel(fuel_index);  // mol O2/mol fuel

    // Convert to mass basis: (mol O2/mol fuel) * (g O2/mol O2) / (g fuel/mol fuel)
    return molar_oxygen_required * oxygen_molar_mass / fuel_molar_mass;  // kg O2/kg fuel
}

// Calculate oxygen required for complete combustion of a mixture [mol O2/mol mixture]
double oxygen_required_per_mol_mixture(const std::vector<double>& X) {
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::runtime_error("Mole fractions must sum to 1.0");
    }

    double total_oxygen_required = 0.0;

    // Calculate weighted sum of oxygen requirements for each fuel component
    for (size_t i = 0; i < X.size(); ++i) {
        // Skip if mole fraction is zero
        if (X[i] <= 0.0) continue;

        // Check if this is a fuel (contains carbon or hydrogen)
        if (molecular_structures[i].C > 0 || molecular_structures[i].H > 0) {
            total_oxygen_required += X[i] * oxygen_required_per_mol_fuel(i);
        }
    }

    return total_oxygen_required;  // mol O2/mol mixture
}

// Calculate oxygen required for complete combustion of a mixture [kg O2/kg mixture]
double oxygen_required_per_kg_mixture(const std::vector<double>& X) {
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::runtime_error("Mole fractions must sum to 1.0");
    }

    // Calculate mixture molecular weight
    double mixture_mw = mwmix(X);  // g/mol

    // Get oxygen molecular weight
    double oxygen_mw = molar_masses[species_index.at("O2")];  // g/mol

    // Calculate molar oxygen requirement for the mixture
    double molar_oxygen_required = oxygen_required_per_mol_mixture(X);  // mol O2/mol mixture

    // If no oxygen required (no fuel in mixture), return 0
    if (molar_oxygen_required <= 0.0) {
        return 0.0;
    }

    // Convert to mass basis: (mol O2/mol mixture) * (g O2/mol O2) / (g mixture/mol mixture)
    return molar_oxygen_required * oxygen_mw / mixture_mw;  // kg O2/kg mixture
}

// -------------------------------------------------------------
// Dry air requirements (using standard dry air composition)
// -------------------------------------------------------------

// Calculate dry air required for complete combustion [mol air/mol fuel]
double dryair_required_per_mol_fuel(std::size_t fuel_index) {
    double O2_required = oxygen_required_per_mol_fuel(fuel_index);
    double X_O2_air = dry_air_composition.at("O2");
    return O2_required / X_O2_air;
}

// Calculate dry air required for complete combustion [kg air/kg fuel]
double dryair_required_per_kg_fuel(std::size_t fuel_index) {
    if (fuel_index >= molar_masses.size()) {
        throw std::runtime_error("Invalid fuel index");
    }

    double fuel_mw = molar_masses[fuel_index];  // g/mol
    std::vector<double> X_air = standard_dry_air_composition();
    double air_mw = mwmix(X_air);  // g/mol

    double molar_air_required = dryair_required_per_mol_fuel(fuel_index);
    return molar_air_required * air_mw / fuel_mw;
}

// Calculate dry air required for complete combustion [mol air/mol mixture]
double dryair_required_per_mol_mixture(const std::vector<double>& X) {
    double O2_required = oxygen_required_per_mol_mixture(X);
    double X_O2_air = dry_air_composition.at("O2");
    return O2_required / X_O2_air;
}

// Calculate dry air required for complete combustion [kg air/kg mixture]
double dryair_required_per_kg_mixture(const std::vector<double>& X) {
    double mixture_mw = mwmix(X);  // g/mol
    std::vector<double> X_air = standard_dry_air_composition();
    double air_mw = mwmix(X_air);  // g/mol

    double molar_air_required = dryair_required_per_mol_mixture(X);
    return molar_air_required * air_mw / mixture_mw;
}

std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double>& X, double& fuel_burn_fraction) {
    const double tol_sum = 1.0e-5;
    const double tol     = 1.0e-12;

    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > tol_sum) {
        throw std::runtime_error("Mole fractions must sum to 1.0");
    }

    if (X.size() != molecular_structures.size()) {
        throw std::runtime_error("Mixture size does not match number of species");
    }

    // 1 mol of mixture as basis
    std::vector<double> n = X;

    const std::size_t idx_O2  = species_index.at("O2");
    const std::size_t idx_CO2 = species_index.at("CO2");
    const std::size_t idx_H2O = species_index.at("H2O");
    const std::size_t idx_N2  = species_index.at("N2");

    double n_O2_available = n[idx_O2];
    if (n_O2_available <= tol) {
        // No oxygen --> no combustion
        fuel_burn_fraction = 0.0;
        return X;
    }

    // Total O2 needed to fully burn all fuels in this mixture (1 mol basis)
    double n_O2_required_total = oxygen_required_per_mol_mixture(X);
    if (n_O2_required_total <= tol) {
        // No fuel --> nothing to burn
        fuel_burn_fraction = 0.0;
        return X;
    }

    // Fraction of fuel that can burn:
    // f = 1 if O2 >= stoich, else f = O2_available / O2_required_total
    double f = 1.0;
    if (n_O2_available + tol < n_O2_required_total) {
        f = n_O2_available / n_O2_required_total;  // 0 < f < 1
    }
    fuel_burn_fraction = f;

    double delta_CO2 = 0.0;
    double delta_H2O = 0.0;
    double delta_N2  = 0.0;

    // Loop over species and burn a fraction f of every O2-consuming fuel
    for (size_t i = 0; i < n.size(); ++i) {
        if (n[i] <= 0.0) continue;

        // Only treat species that actually require O2 (CO2/H2O have requirement=0)
        double nu_O2 = oxygen_required_per_mol_fuel(i);
        if (nu_O2 <= 0.0) continue;

        const Molecular_Structure& sp = molecular_structures[i];

        double n_fuel_initial   = n[i];
        double n_fuel_reacted   = f * n_fuel_initial;
        double n_fuel_remaining = n_fuel_initial - n_fuel_reacted;

        n[i] = n_fuel_remaining;

        // C_x H_y O_z N_w + ... -> x CO2 + (y/2) H2O + (w/2) N2
        delta_CO2 += static_cast<double>(sp.C) * n_fuel_reacted;
        delta_H2O += 0.5 * static_cast<double>(sp.H) * n_fuel_reacted;
        delta_N2  += 0.5 * static_cast<double>(sp.N) * n_fuel_reacted;
    }

    // Consume O2: fully if O2-limited, partially if fuel-limited
    if (f >= 1.0 - tol) {
        n[idx_O2] -= n_O2_required_total;
        if (n[idx_O2] < 0.0) n[idx_O2] = 0.0;  // numeric safety
    } else {
        n[idx_O2] = 0.0;  // all available O2 consumed
    }

    // Add products
    n[idx_CO2] += delta_CO2;
    n[idx_H2O] += delta_H2O;
    n[idx_N2]  += delta_N2;

    // Back to mole fractions
    return normalize_fractions(n);
}

std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double>& X) {
    double fuel_burn_fraction = 0.0;
    return complete_combustion_to_CO2_H2O(X, fuel_burn_fraction);
}

// -------------------------------------------------------------
// Equivalence ratio helpers (mole basis, multi-species fuel)
// -------------------------------------------------------------

// We assume species_names, molar_masses and oxygen_required_per_mol_mixture(X)
// are defined in thermo_transport_data / this translation unit.

// Stoichiometric fuel-to-oxidizer mole ratio (F/O)_st for a fuel mixture
// (mole fractions X_fuel) and an oxidizer mixture X_ox (mole fractions).
//
// Uses:
//   - oxygen_required_per_mol_mixture(X_fuel): moles O2 per mole fuel mixture
//   - O2 mole fraction in oxidizer stream X_ox[O2_index]
static double stoich_f_over_o_mole(
    const std::vector<double>& X_fuel,
    const std::vector<double>& X_ox)
{
    if (X_fuel.size() != species_names.size() ||
        X_ox.size()   != species_names.size()) {
        throw std::invalid_argument("stoich_f_over_o_mole: size mismatch");
    }

    const double X_O2_ox = X_ox[species_index.at("O2")];
    if (X_O2_ox <= 0.0) {
        throw std::runtime_error(
            "stoich_f_over_o_mole: oxidizer has zero O2 mole fraction");
    }

    // moles O2 required per mole of fuel mixture
    const double nu_O2_req = oxygen_required_per_mol_mixture(X_fuel);
    if (nu_O2_req <= 0.0) {
        throw std::runtime_error(
            "stoich_f_over_o_mole: oxygen_required_per_mol_mixture returned non-positive");
    }

    // Let:
    //   nu_O2_req = X_O2_ox * n_ox_st  ->  n_ox_st = nu_O2_req / X_O2_ox
    // then stoich F/O (molar) = n_F / n_Ox = 1 / n_ox_st = X_O2_ox / nu_O2_req
    const double r_st = X_O2_ox / nu_O2_req;
    return r_st;
}

// Given a mixture X_mix that is formed *only* by mixing:
//   - a fuel stream with composition X_fuel
//   - an oxidizer stream with composition X_ox
// we can write:
//   X_mix = (r * X_fuel + 1 * X_ox) / (r + 1)
// where r = n_F / n_Ox is the fuel-to-oxidizer mole ratio.
//
// For any species j with X_fuel[j] != X_ox[j]:
//   X_mix[j] (r + 1) = r X_fuel[j] + X_ox[j]
//   r (X_mix[j] - X_fuel[j]) = X_ox[j] - X_mix[j]
//   r = (X_ox[j] - X_mix[j]) / (X_mix[j] - X_fuel[j])
//
// We pick the first j with |X_fuel[j] - X_ox[j]| large enough.
static double actual_f_over_o_mole(
    const std::vector<double>& X_mix,
    const std::vector<double>& X_fuel,
    const std::vector<double>& X_ox)
{
    const std::size_t n = species_names.size();
    if (X_mix.size()  != n ||
        X_fuel.size() != n ||
        X_ox.size()   != n) {
        throw std::invalid_argument("actual_f_over_o_mole: size mismatch");
    }

    const double eps = 1e-12;
    bool found = false;
    double r = 0.0;

    for (std::size_t j = 0; j < n; ++j) {
        const double xf = X_fuel[j];
        const double xo = X_ox[j];
        const double xm = X_mix[j];

        double denom = xm - xf;
        if (std::abs(xf - xo) < eps) {
            continue; // no contrast between fuel and oxidizer for this species
        }
        if (std::abs(denom) < eps) {
            continue; // would be numerically unstable
        }

        double r_j = (xo - xm) / denom;

        // Sanity: r_j must be >= 0 for physical mixing
        if (r_j < -1e-8) {
            continue; // skip clearly unphysical projection
        }

        r = r_j;
        found = true;
        break;
    }

    if (!found) {
        throw std::runtime_error(
            "actual_f_over_o_mole: could not determine mixing ratio "
            "from X_mix, X_fuel, X_ox (mixture not on mixing line?)");
    }

    return r;
}

double equivalence_ratio_mole(
    const std::vector<double>& X_mix,
    const std::vector<double>& X_fuel,
    const std::vector<double>& X_ox)
{
    const double r_act = actual_f_over_o_mole(X_mix, X_fuel, X_ox);
    const double r_st  = stoich_f_over_o_mole(X_fuel, X_ox);

    return r_act / r_st; // φ = (F/O)_act / (F/O)_st
}

std::vector<double> set_equivalence_ratio_mole(
    double phi,
    const std::vector<double>& X_fuel,
    const std::vector<double>& X_ox)
{
    const std::size_t n = species_names.size();
    if (X_fuel.size() != n || X_ox.size() != n) {
        throw std::invalid_argument("set_equivalence_ratio_mole: size mismatch");
    }
    if (phi <= 0.0) {
        throw std::invalid_argument("set_equivalence_ratio_mole: phi must be > 0");
    }

    const double r_st = stoich_f_over_o_mole(X_fuel, X_ox);
    const double r    = phi * r_st; // actual F/O = φ * (F/O)_st

    const double n_F  = r;
    const double n_Ox = 1.0;
    const double n_tot = n_F + n_Ox;

    std::vector<double> X_mix(n, 0.0);

    // Fuel-stream contribution
    for (std::size_t k = 0; k < n; ++k) {
        if (X_fuel[k] != 0.0) {
            X_mix[k] += (n_F * X_fuel[k]) / n_tot;
        }
    }

    // Oxidizer-stream contribution
    for (std::size_t k = 0; k < n; ++k) {
        if (X_ox[k] != 0.0) {
            X_mix[k] += (n_Ox * X_ox[k]) / n_tot;
        }
    }

    // Renormalize to guard against small numerical drifts
    double sum = std::accumulate(X_mix.begin(), X_mix.end(), 0.0);
    if (sum <= 0.0) {
        throw std::runtime_error(
            "set_equivalence_ratio_mole: resulting mixture has non-positive sum");
    }
    for (double& x : X_mix) {
        x /= sum;
    }

    return X_mix;
}


// -------------------------------------------------------------
// Stoichiometric F/O (mass basis)
// -------------------------------------------------------------

// Stoichiometric fuel-to-oxidizer mass ratio (F/O)_st for a fuel mixture
// (mass fractions Y_fuel) and an oxidizer mixture Y_ox (mass fractions).
//
// Uses:
//   - oxygen_required_per_kg_mixture(X_fuel): kg O2 per kg fuel mixture
//   - O2 mass fraction in oxidizer stream Y_ox[O2_index]
static double stoich_f_over_o_mass(
    const std::vector<double>& Y_fuel,
    const std::vector<double>& Y_ox)
{
    if (Y_fuel.size() != species_names.size() ||
        Y_ox.size()   != species_names.size()) {
        throw std::invalid_argument("stoich_f_over_o_mass: size mismatch");
    }

    const double Y_O2_ox = Y_ox[species_index.at("O2")];
    if (Y_O2_ox <= 0.0) {
        throw std::runtime_error(
            "stoich_f_over_o_mass: oxidizer has zero O2 mass fraction");
    }

    // Convert fuel stream from mass fractions to mole fractions for use in
    // oxygen_required_per_kg_mixture(X_fuel).
    std::vector<double> X_fuel = mass_to_mole(Y_fuel);

    // kg O2 per kg fuel mixture
    const double m_O2_req = oxygen_required_per_kg_mixture(X_fuel);
    if (m_O2_req <= 0.0) {
        throw std::runtime_error(
            "stoich_f_over_o_mass: oxygen_required_per_kg_mixture returned non-positive");
    }

    // Stoichiometric oxidizer mass (kg) per 1 kg fuel:
    //   m_O2_req = Y_O2_ox * m_ox_st  ->  m_ox_st = m_O2_req / Y_O2_ox
    // hence (F/O)_st (mass) = m_fuel / m_ox_st = 1 / (m_O2_req / Y_O2_ox)
    //                        = Y_O2_ox / m_O2_req
    const double r_st = Y_O2_ox / m_O2_req;
    return r_st;
}

// -------------------------------------------------------------
// Actual F/O (mass basis) from a given mixture Y_mix
// -------------------------------------------------------------

// For mass mixing of fuel and oxidizer:
//   Y_mix = (m_F Y_fuel + m_Ox Y_ox) / (m_F + m_Ox),
// let r = m_F / m_Ox. For any species j with Y_fuel[j] != Y_ox[j]:
//
//   Y_mix[j] (r + 1) = r Y_fuel[j] + Y_ox[j]
//   r (Y_mix[j] - Y_fuel[j]) = Y_ox[j] - Y_mix[j]
//   r = (Y_ox[j] - Y_mix[j]) / (Y_mix[j] - Y_fuel[j])
//
// We pick the first j that gives a numerically stable r.
static double actual_f_over_o_mass(
    const std::vector<double>& Y_mix,
    const std::vector<double>& Y_fuel,
    const std::vector<double>& Y_ox)
{
    const std::size_t n = species_names.size();
    if (Y_mix.size()  != n ||
        Y_fuel.size() != n ||
        Y_ox.size()   != n) {
        throw std::invalid_argument("actual_f_over_o_mass: size mismatch");
    }

    const double eps = 1e-12;
    bool found = false;
    double r = 0.0;

    for (std::size_t j = 0; j < n; ++j) {
        const double yf = Y_fuel[j];
        const double yo = Y_ox[j];
        const double ym = Y_mix[j];

        if (std::abs(yf - yo) < eps) {
            continue; // no contrast between fuel and oxidizer for this species
        }

        const double denom = ym - yf;
        if (std::abs(denom) < eps) {
            continue; // numerically unstable
        }

        const double r_j = (yo - ym) / denom;

        // Require physically reasonable (non-negative) r_j
        if (r_j < -1e-8) {
            continue;
        }

        r = r_j;
        found = true;
        break;
    }

    if (!found) {
        throw std::runtime_error(
            "actual_f_over_o_mass: could not determine mixing ratio "
            "from Y_mix, Y_fuel, Y_ox (mixture not on mixing line?)");
    }

    return r;
}

// -------------------------------------------------------------
// Public mass-basis φ functions
// -------------------------------------------------------------

double equivalence_ratio_mass(
    const std::vector<double>& Y_mix,
    const std::vector<double>& Y_fuel,
    const std::vector<double>& Y_ox)
{
    const double r_act = actual_f_over_o_mass(Y_mix, Y_fuel, Y_ox);
    const double r_st  = stoich_f_over_o_mass(Y_fuel, Y_ox);

    return r_act / r_st; // φ = (F/O)_act / (F/O)_st (mass basis)
}

std::vector<double> set_equivalence_ratio_mass(
    double phi,
    const std::vector<double>& Y_fuel,
    const std::vector<double>& Y_ox)
{
    const std::size_t n = species_names.size();
    if (Y_fuel.size() != n || Y_ox.size() != n) {
        throw std::invalid_argument("set_equivalence_ratio_mass: size mismatch");
    }
    if (phi <= 0.0) {
        throw std::invalid_argument("set_equivalence_ratio_mass: phi must be > 0");
    }

    const double r_st = stoich_f_over_o_mass(Y_fuel, Y_ox);
    const double r    = phi * r_st; // (F/O)_act = φ * (F/O)_st

    const double m_F  = r;
    const double m_Ox = 1.0;
    const double m_tot = m_F + m_Ox;

    std::vector<double> Y_mix(n, 0.0);

    // Fuel-stream contribution
    for (std::size_t k = 0; k < n; ++k) {
        if (Y_fuel[k] != 0.0) {
            Y_mix[k] += (m_F * Y_fuel[k]) / m_tot;
        }
    }

    // Oxidizer-stream contribution
    for (std::size_t k = 0; k < n; ++k) {
        if (Y_ox[k] != 0.0) {
            Y_mix[k] += (m_Ox * Y_ox[k]) / m_tot;
        }
    }

    // Renormalize to guard against any small numerical drift
    double sum = std::accumulate(Y_mix.begin(), Y_mix.end(), 0.0);
    if (sum <= 0.0) {
        throw std::runtime_error(
            "set_equivalence_ratio_mass: resulting mixture has non-positive sum");
    }
    for (double& y : Y_mix) {
        y /= sum;
    }

    return Y_mix;
}

double bilger_stoich_mixture_fraction_mass(
    const std::vector<double>& Y_F,
    const std::vector<double>& Y_O)
{
    // r_st = (F/O)_st on a MASS basis
    // already implemented as part of equivalence_ratio_mass helpers:
    const double r_st = stoich_f_over_o_mass(Y_F, Y_O);

    if (r_st <= 0.0) {
        throw std::runtime_error(
            "bilger_stoich_mixture_fraction_mass: r_st <= 0");
    }

    // For 2-stream fuel/oxidizer mixing we always have:
    //
    //     Z = m_F / (m_F + m_O)
    //     r = m_F / m_O
    //
    //  =>  Z = r / (1 + r)
    //
    // At stoichiometry r = r_st, giving:
    //
    //     Z_st = r_st / (1 + r_st)
    //
    return r_st / (1.0 + r_st);
}

double bilger_Z_from_equivalence_ratio_mass(
    double phi,
    const std::vector<double>& Y_F,
    const std::vector<double>& Y_O)
{
    if (phi <= 0.0) {
        throw std::invalid_argument(
            "bilger_Z_from_equivalence_ratio_mass: phi must be > 0");
    }

    const double r_st = stoich_f_over_o_mass(Y_F, Y_O);
    if (r_st <= 0.0) {
        throw std::runtime_error(
            "bilger_Z_from_equivalence_ratio_mass: r_st <= 0");
    }

    // r = (F/O)_act (mass) = phi * r_st
    const double r = phi * r_st;

    // Two-stream mixing relation on a mass basis:
    //   Z = r / (1 + r)
    return r / (1.0 + r);
}

double equivalence_ratio_from_bilger_Z_mass(
    double Z,
    const std::vector<double>& Y_F,
    const std::vector<double>& Y_O)
{
    if (Z <= 0.0 || Z >= 1.0) {
        throw std::invalid_argument(
            "equivalence_ratio_from_bilger_Z_mass: Z must be in (0, 1)");
    }

    const double r_st = stoich_f_over_o_mass(Y_F, Y_O);
    if (r_st <= 0.0) {
        throw std::runtime_error(
            "equivalence_ratio_from_bilger_Z_mass: r_st <= 0");
    }

    // Invert Z = r / (1 + r) -> r = Z / (1 - Z)
    const double r = Z / (1.0 - Z);

    // phi = (F/O)_act / (F/O)_st = r / r_st
    return r / r_st;
}

// ============================
// Bilger mixture fraction
// ============================

// Bilger scalar beta(Y) based on element mass fractions of C, H, O.
// Uses:
//   molecular_structures[k].C/H/O  (atom counts per species)
//   molar_masses[k]                (species molar masses, g/mol)
//
// For each element e:
//   Z_e / W_e = sum_k Y_k * n_e,k / W_k
// so Bilger's beta can be written without explicit atomic weights:
//   beta = 2 * sum_k Y_k C_k / W_k
//        + 0.5 * sum_k Y_k H_k / W_k
//        -       sum_k Y_k O_k / W_k
double bilger_beta(const std::vector<double>& Y)
{
    if (Y.size() != molar_masses.size() ||
        Y.size() != molecular_structures.size()) {
        throw std::invalid_argument(
            "bilger_beta: mass-fraction vector size "
            "does not match species data");
    }

    double sum_C = 0.0;
    double sum_H = 0.0;
    double sum_O = 0.0;

    for (std::size_t k = 0; k < Y.size(); ++k) {
        double Yk = Y[k];
        if (Yk == 0.0) continue;

        const auto& ms = molecular_structures[k];
        double Wk = molar_masses[k]; // g/mol

        if (Wk <= 0.0) {
            throw std::runtime_error("bilger_beta: non-positive molar mass");
        }

        // Contributions to Z_e / W_e, see comment above
        sum_C += Yk * static_cast<double>(ms.C) / Wk;
        sum_H += Yk * static_cast<double>(ms.H) / Wk;
        sum_O += Yk * static_cast<double>(ms.O) / Wk;
    }

    // Bilger's scalar (proportional to 2 Z_C/W_C + 0.5 Z_H/W_H - Z_O/W_O)
    return 2.0 * sum_C + 0.5 * sum_H - sum_O;
}

// Bilger mixture fraction Z in [0,1] between two streams F and O,
// defined as a normalized Bilger scalar:
//   Z = (beta - beta_O) / (beta_F - beta_O)
double bilger_mixture_fraction(
    const std::vector<double>& Y,
    const std::vector<double>& Y_F,
    const std::vector<double>& Y_O)
{
    if (Y.size() != Y_F.size() || Y.size() != Y_O.size()) {
        throw std::invalid_argument(
            "bilger_mixture_fraction: all mass-fraction vectors "
            "must have the same size");
    }

    const double beta   = bilger_beta(Y);
    const double beta_F = bilger_beta(Y_F);
    const double beta_O = bilger_beta(Y_O);

    const double denom = beta_F - beta_O;
    if (std::abs(denom) < 1e-16) {
        throw std::runtime_error(
            "bilger_mixture_fraction: fuel and oxidizer have "
            "identical Bilger scalar (denominator ~ 0)");
    }

    double Z = (beta - beta_O) / denom;

    // Clamp to [0,1] to handle small numerical overshoots
    if (Z < 0.0)      Z = 0.0;
    else if (Z > 1.0) Z = 1.0;

    return Z;
}

// Convenience wrapper: inputs are mole fractions X, X_F, X_O.
// These are converted to mass fractions first and then passed to
// the mass-fraction-based Bilger implementation above.
double bilger_mixture_fraction_from_moles(
    const std::vector<double>& X,
    const std::vector<double>& X_F,
    const std::vector<double>& X_O)
{
    if (X.size() != X_F.size() || X.size() != X_O.size()) {
        throw std::invalid_argument(
            "bilger_mixture_fraction_from_moles: all mole-fraction "
            "vectors must have the same size");
    }

    const auto Y   = mole_to_mass(X);
    const auto Y_F = mole_to_mass(X_F);
    const auto Y_O = mole_to_mass(X_O);

    return bilger_mixture_fraction(Y, Y_F, Y_O);
}

// -------------------------------------------------------------
// Stream-based equivalence ratio helpers
// -------------------------------------------------------------

Stream set_fuel_stream_for_phi(double phi, const Stream& fuel, const Stream& oxidizer) {
    if (phi <= 0.0) {
        throw std::invalid_argument("set_fuel_stream_for_phi: phi must be positive");
    }
    if (oxidizer.mdot <= 0.0) {
        throw std::invalid_argument("set_fuel_stream_for_phi: oxidizer.mdot must be positive");
    }

    // Convert to mass fractions
    std::vector<double> Y_fuel = mole_to_mass(fuel.X());
    std::vector<double> Y_ox = mole_to_mass(oxidizer.X());

    // Get mixture mass fractions for target phi
    std::vector<double> Y_mix = set_equivalence_ratio_mass(phi, Y_fuel, Y_ox);

    // Calculate fuel mass flow from mixture composition
    // Y_mix = (mdot_fuel * Y_fuel + mdot_ox * Y_ox) / (mdot_fuel + mdot_ox)
    //
    // For a fuel-only species k (where Y_ox[k] = 0):
    //   Y_mix[k] = mdot_fuel * Y_fuel[k] / (mdot_fuel + mdot_ox)
    //   => mdot_fuel / (mdot_fuel + mdot_ox) = Y_mix[k] / Y_fuel[k]
    //
    // Find a fuel species (present in fuel but not in oxidizer)
    double alpha = -1.0;
    for (std::size_t i = 0; i < Y_fuel.size(); ++i) {
        if (Y_fuel[i] > 1e-10 && Y_ox[i] < 1e-10) {
            // This is a fuel-only species
            alpha = Y_mix[i] / Y_fuel[i];
            break;
        }
    }

    if (alpha < 0.0) {
        throw std::runtime_error(
            "set_fuel_stream_for_phi: no fuel-only species found "
            "(fuel must contain species not present in oxidizer)");
    }

    // mdot_fuel / (mdot_fuel + mdot_ox) = alpha
    // => mdot_fuel = alpha * mdot_ox / (1 - alpha)
    if (alpha >= 1.0 || alpha <= 0.0) {
        throw std::runtime_error("set_fuel_stream_for_phi: invalid mass fraction ratio");
    }

    double mdot_fuel = alpha * oxidizer.mdot / (1.0 - alpha);

    // Return fuel stream with computed mdot
    Stream result = fuel;
    result.mdot = mdot_fuel;
    return result;
}

// -------------------------------------------------------------
// Inverse solvers for fuel/oxidizer streams (complete combustion only)
// -------------------------------------------------------------

namespace {

// Bisection solver for monotonic functions
// Finds x such that f(x) = target, given that f is monotonic on [x_lo, x_hi]
// tol_f: tolerance on function value |f(x) - target|
// tol_x_rel: relative tolerance on x, convergence when (x_hi - x_lo) < tol_x_rel * x_mid
template<typename Func>
double bisection_solve(Func f, double target, double x_lo, double x_hi,
                       double tol_f, std::size_t max_iter, double tol_x_rel = 1e-10) {
    double f_lo = f(x_lo);
    double f_hi = f(x_hi);

    // Check if target is within range
    double f_min = std::min(f_lo, f_hi);
    double f_max = std::max(f_lo, f_hi);
    if (target < f_min - tol_f || target > f_max + tol_f) {
        throw std::runtime_error("bisection_solve: target outside achievable range");
    }

    // Determine monotonicity direction
    bool increasing = (f_hi > f_lo);

    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        double x_mid = 0.5 * (x_lo + x_hi);
        double f_mid = f(x_mid);

        // Check convergence on function value
        if (std::abs(f_mid - target) < tol_f) {
            return x_mid;
        }

        // Check convergence on x (relative)
        if ((x_hi - x_lo) < tol_x_rel * x_mid) {
            return x_mid;
        }

        bool mid_below_target = (f_mid < target);
        if (increasing == mid_below_target) {
            x_lo = x_mid;
        } else {
            x_hi = x_mid;
        }
    }

    return 0.5 * (x_lo + x_hi);
}

// Estimate stoichiometric mass flow ratio (mdot_fuel / mdot_ox at stoich)
double stoich_fuel_ox_ratio(const std::vector<double>& X_fuel, const std::vector<double>& X_ox) {
    double O2_req_per_kg_fuel = oxygen_required_per_kg_mixture(X_fuel);
    double X_O2_ox = X_ox[species_index.at("O2")];
    double MW_ox = mwmix(X_ox);
    double MW_O2 = molar_masses[species_index.at("O2")];
    double Y_O2_ox = X_O2_ox * MW_O2 / MW_ox;

    if (O2_req_per_kg_fuel > 0.0 && Y_O2_ox > 0.0) {
        return Y_O2_ox / O2_req_per_kg_fuel;
    }
    return 0.05;  // fallback
}

// Helper: compute Tad for given fuel mdot (oxidizer mdot fixed)
double compute_Tad_fuel(double mdot_fuel, const Stream& fuel, const Stream& oxidizer) {
    Stream fuel_trial = fuel;
    fuel_trial.mdot = mdot_fuel;
    Stream mixed = mix({fuel_trial, oxidizer});
    State burned = complete_combustion(mixed.state);
    return burned.T;
}

// Helper: compute Tad for given oxidizer mdot (fuel mdot fixed)
double compute_Tad_ox(double mdot_ox, const Stream& fuel, const Stream& oxidizer) {
    Stream ox_trial = oxidizer;
    ox_trial.mdot = mdot_ox;
    Stream mixed = mix({fuel, ox_trial});
    State burned = complete_combustion(mixed.state);
    return burned.T;
}

// Helper: compute species mole fraction in burned products
// vary_fuel: true = vary fuel mdot, false = vary oxidizer mdot
// dry_basis: true = remove H2O before computing mole fraction
double compute_species_burned(double mdot_var, const Stream& fuel, const Stream& oxidizer,
                               const std::string& species, bool vary_fuel, bool dry_basis) {
    Stream fuel_trial = fuel;
    Stream ox_trial = oxidizer;
    if (vary_fuel) {
        fuel_trial.mdot = mdot_var;
    } else {
        ox_trial.mdot = mdot_var;
    }
    Stream mixed = mix({fuel_trial, ox_trial});
    std::vector<double> X_burned = complete_combustion_to_CO2_H2O(mixed.state.X);
    if (dry_basis) {
        X_burned = convert_to_dry_fractions(X_burned);
    }
    return X_burned[species_index.at(species)];
}

// Convenience wrappers for backward compatibility and readability
double compute_O2_burned_fuel(double mdot, const Stream& f, const Stream& o) {
    return compute_species_burned(mdot, f, o, "O2", true, false);
}
double compute_O2_burned_ox(double mdot, const Stream& f, const Stream& o) {
    return compute_species_burned(mdot, f, o, "O2", false, false);
}
double compute_CO2_burned_fuel(double mdot, const Stream& f, const Stream& o) {
    return compute_species_burned(mdot, f, o, "CO2", true, false);
}
double compute_CO2_burned_ox(double mdot, const Stream& f, const Stream& o) {
    return compute_species_burned(mdot, f, o, "CO2", false, false);
}
double compute_O2_dry_burned_fuel(double mdot, const Stream& f, const Stream& o) {
    return compute_species_burned(mdot, f, o, "O2", true, true);
}
double compute_O2_dry_burned_ox(double mdot, const Stream& f, const Stream& o) {
    return compute_species_burned(mdot, f, o, "O2", false, true);
}
double compute_CO2_dry_burned_fuel(double mdot, const Stream& f, const Stream& o) {
    return compute_species_burned(mdot, f, o, "CO2", true, true);
}
double compute_CO2_dry_burned_ox(double mdot, const Stream& f, const Stream& o) {
    return compute_species_burned(mdot, f, o, "CO2", false, true);
}

}  // anonymous namespace

// -------------------------------------------------------------
// Find fuel stream for target property (oxidizer mdot fixed)
// -------------------------------------------------------------

Stream set_fuel_stream_for_Tad(double T_ad_target, const Stream& fuel, const Stream& oxidizer,
                                double tol, std::size_t max_iter, bool lean, double phi_max) {
    if (oxidizer.mdot <= 0.0) {
        throw std::invalid_argument("set_fuel_stream_for_Tad: oxidizer.mdot must be positive");
    }
    if (phi_max <= 1.0) {
        throw std::invalid_argument("set_fuel_stream_for_Tad: phi_max must be > 1.0");
    }

    // Tad must be above oxidizer temperature (minimum achievable with zero fuel)
    double T_min = oxidizer.T();
    if (T_ad_target <= T_min) {
        throw std::invalid_argument("set_fuel_stream_for_Tad: T_ad_target must be > oxidizer temperature");
    }
    if (T_ad_target > 5000.0) {
        throw std::invalid_argument("set_fuel_stream_for_Tad: T_ad_target must be <= 5000 K");
    }

    // Estimate stoichiometric fuel mass flow
    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = oxidizer.mdot * ratio;

    // Peak Tad is at stoichiometric
    double mdot_peak = mdot_stoich;
    double T_peak = compute_Tad_fuel(mdot_peak, fuel, oxidizer);

    if (T_ad_target > T_peak + tol) {
        throw std::runtime_error("set_fuel_stream_for_Tad: target T_ad exceeds maximum achievable (stoichiometric)");
    }

    double mdot_result = 0.0;
    auto f = [&](double mdot) { return compute_Tad_fuel(mdot, fuel, oxidizer); };

    if (lean) {
        // Search on lean side: low fuel to stoichiometric
        // Tad increases with fuel on lean side
        double mdot_lo = mdot_stoich * 0.01;
        double T_lo = compute_Tad_fuel(mdot_lo, fuel, oxidizer);

        if (T_ad_target < T_lo - tol) {
            throw std::runtime_error("set_fuel_stream_for_Tad: target T_ad below minimum achievable (lean)");
        }

        mdot_result = bisection_solve(f, T_ad_target, mdot_lo, mdot_peak, tol, max_iter);
    } else {
        // Search on rich side: stoichiometric to phi_max
        // Tad decreases with fuel on rich side
        double mdot_hi = mdot_stoich * phi_max;
        double T_hi = compute_Tad_fuel(mdot_hi, fuel, oxidizer);

        if (T_ad_target < T_hi - tol) {
            throw std::runtime_error("set_fuel_stream_for_Tad: target T_ad below minimum achievable (rich)");
        }

        mdot_result = bisection_solve(f, T_ad_target, mdot_peak, mdot_hi, tol, max_iter);
    }

    Stream result = fuel;
    result.mdot = mdot_result;
    return result;
}

Stream set_fuel_stream_for_O2(double X_O2_target, const Stream& fuel, const Stream& oxidizer,
                               double tol, std::size_t max_iter) {
    if (oxidizer.mdot <= 0.0) {
        throw std::invalid_argument("set_fuel_stream_for_O2: oxidizer.mdot must be positive");
    }

    double X_O2_ox = oxidizer.X()[species_index.at("O2")];
    if (X_O2_target <= 0.0 || X_O2_target >= X_O2_ox) {
        throw std::invalid_argument("set_fuel_stream_for_O2: X_O2_target must be in (0, X_O2_oxidizer)");
    }

    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = oxidizer.mdot * ratio;

    double mdot_lo = mdot_stoich * 0.001;
    double mdot_hi = mdot_stoich * 0.999;

    auto f = [&](double mdot) { return compute_O2_burned_fuel(mdot, fuel, oxidizer); };
    double mdot_result = bisection_solve(f, X_O2_target, mdot_lo, mdot_hi, tol, max_iter);

    Stream result = fuel;
    result.mdot = mdot_result;
    return result;
}

Stream set_fuel_stream_for_O2_dry(double X_O2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                   double tol, std::size_t max_iter) {
    if (oxidizer.mdot <= 0.0) {
        throw std::invalid_argument("set_fuel_stream_for_O2_dry: oxidizer.mdot must be positive");
    }

    // Dry O2 upper bound is from dry oxidizer composition
    std::vector<double> X_ox_dry = convert_to_dry_fractions(oxidizer.X());
    double X_O2_ox_dry = X_ox_dry[species_index.at("O2")];
    if (X_O2_dry_target <= 0.0 || X_O2_dry_target >= X_O2_ox_dry) {
        throw std::invalid_argument("set_fuel_stream_for_O2_dry: X_O2_dry_target must be in (0, X_O2_dry_oxidizer)");
    }

    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = oxidizer.mdot * ratio;

    double mdot_lo = mdot_stoich * 0.001;
    double mdot_hi = mdot_stoich * 0.999;

    auto f = [&](double mdot) { return compute_O2_dry_burned_fuel(mdot, fuel, oxidizer); };
    double mdot_result = bisection_solve(f, X_O2_dry_target, mdot_lo, mdot_hi, tol, max_iter);

    Stream result = fuel;
    result.mdot = mdot_result;
    return result;
}

Stream set_fuel_stream_for_CO2(double X_CO2_target, const Stream& fuel, const Stream& oxidizer,
                                double tol, std::size_t max_iter) {
    if (oxidizer.mdot <= 0.0) {
        throw std::invalid_argument("set_fuel_stream_for_CO2: oxidizer.mdot must be positive");
    }
    if (X_CO2_target <= 0.0) {
        throw std::invalid_argument("set_fuel_stream_for_CO2: X_CO2_target must be positive");
    }

    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = oxidizer.mdot * ratio;

    double mdot_lo = mdot_stoich * 0.001;
    double mdot_hi = mdot_stoich * 0.999;

    auto f = [&](double mdot) { return compute_CO2_burned_fuel(mdot, fuel, oxidizer); };
    double mdot_result = bisection_solve(f, X_CO2_target, mdot_lo, mdot_hi, tol, max_iter);

    Stream result = fuel;
    result.mdot = mdot_result;
    return result;
}

Stream set_fuel_stream_for_CO2_dry(double X_CO2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                    double tol, std::size_t max_iter) {
    if (oxidizer.mdot <= 0.0) {
        throw std::invalid_argument("set_fuel_stream_for_CO2_dry: oxidizer.mdot must be positive");
    }
    if (X_CO2_dry_target <= 0.0) {
        throw std::invalid_argument("set_fuel_stream_for_CO2_dry: X_CO2_dry_target must be positive");
    }

    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = oxidizer.mdot * ratio;

    double mdot_lo = mdot_stoich * 0.001;
    double mdot_hi = mdot_stoich * 0.999;

    auto f = [&](double mdot) { return compute_CO2_dry_burned_fuel(mdot, fuel, oxidizer); };
    double mdot_result = bisection_solve(f, X_CO2_dry_target, mdot_lo, mdot_hi, tol, max_iter);

    Stream result = fuel;
    result.mdot = mdot_result;
    return result;
}

// -------------------------------------------------------------
// Find oxidizer stream for target property (fuel mdot fixed)
// -------------------------------------------------------------

Stream set_oxidizer_stream_for_Tad(double T_ad_target, const Stream& fuel, const Stream& oxidizer,
                                    double tol, std::size_t max_iter, bool lean, double phi_max) {
    if (fuel.mdot <= 0.0) {
        throw std::invalid_argument("set_oxidizer_stream_for_Tad: fuel.mdot must be positive");
    }
    if (phi_max <= 1.0) {
        throw std::invalid_argument("set_oxidizer_stream_for_Tad: phi_max must be > 1.0");
    }

    // Tad must be above fuel temperature (minimum achievable with infinite air)
    double T_min = fuel.T();
    if (T_ad_target <= T_min) {
        throw std::invalid_argument("set_oxidizer_stream_for_Tad: T_ad_target must be > fuel temperature");
    }
    if (T_ad_target > 5000.0) {
        throw std::invalid_argument("set_oxidizer_stream_for_Tad: T_ad_target must be <= 5000 K");
    }

    // Estimate stoichiometric oxidizer mass flow
    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

    // Peak Tad is at stoichiometric
    double mdot_peak = mdot_stoich;
    double T_peak = compute_Tad_ox(mdot_peak, fuel, oxidizer);

    if (T_ad_target > T_peak + tol) {
        throw std::runtime_error("set_oxidizer_stream_for_Tad: target T_ad exceeds maximum achievable (stoichiometric)");
    }

    double mdot_result = 0.0;
    auto f = [&](double mdot) { return compute_Tad_ox(mdot, fuel, oxidizer); };

    if (lean) {
        // Search on lean side: stoichiometric to high oxidizer (1/phi_min factor)
        // Tad decreases with oxidizer on lean side
        double mdot_hi = mdot_stoich * phi_max;  // phi_max used as dilution factor
        double T_hi = compute_Tad_ox(mdot_hi, fuel, oxidizer);

        if (T_ad_target < T_hi - tol) {
            throw std::runtime_error("set_oxidizer_stream_for_Tad: target T_ad below minimum achievable (lean)");
        }

        mdot_result = bisection_solve(f, T_ad_target, mdot_peak, mdot_hi, tol, max_iter);
    } else {
        // Search on rich side: low oxidizer (1/phi_max) to stoichiometric
        // Tad increases with oxidizer on rich side
        double mdot_lo = mdot_stoich / phi_max;
        double T_lo = compute_Tad_ox(mdot_lo, fuel, oxidizer);

        if (T_ad_target < T_lo - tol) {
            throw std::runtime_error("set_oxidizer_stream_for_Tad: target T_ad below minimum achievable (rich)");
        }

        mdot_result = bisection_solve(f, T_ad_target, mdot_lo, mdot_peak, tol, max_iter);
    }

    Stream result = oxidizer;
    result.mdot = mdot_result;
    return result;
}

Stream set_oxidizer_stream_for_O2(double X_O2_target, const Stream& fuel, const Stream& oxidizer,
                                   double tol, std::size_t max_iter) {
    if (fuel.mdot <= 0.0) {
        throw std::invalid_argument("set_oxidizer_stream_for_O2: fuel.mdot must be positive");
    }

    double X_O2_ox = oxidizer.X()[species_index.at("O2")];
    if (X_O2_target <= 0.0 || X_O2_target >= X_O2_ox) {
        throw std::invalid_argument("set_oxidizer_stream_for_O2: X_O2_target must be in (0, X_O2_oxidizer)");
    }

    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

    double mdot_lo = mdot_stoich * 1.001;  // just above stoich -> low O2
    double mdot_hi = mdot_stoich * 100.0;  // very lean -> high O2

    // O2 increases with increasing oxidizer (on lean side)
    auto f = [&](double mdot) { return compute_O2_burned_ox(mdot, fuel, oxidizer); };
    double mdot_result = bisection_solve(f, X_O2_target, mdot_lo, mdot_hi, tol, max_iter);

    Stream result = oxidizer;
    result.mdot = mdot_result;
    return result;
}

Stream set_oxidizer_stream_for_O2_dry(double X_O2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                       double tol, std::size_t max_iter) {
    if (fuel.mdot <= 0.0) {
        throw std::invalid_argument("set_oxidizer_stream_for_O2_dry: fuel.mdot must be positive");
    }

    std::vector<double> X_ox_dry = convert_to_dry_fractions(oxidizer.X());
    double X_O2_ox_dry = X_ox_dry[species_index.at("O2")];
    if (X_O2_dry_target <= 0.0 || X_O2_dry_target >= X_O2_ox_dry) {
        throw std::invalid_argument("set_oxidizer_stream_for_O2_dry: X_O2_dry_target must be in (0, X_O2_dry_oxidizer)");
    }

    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

    double mdot_lo = mdot_stoich * 1.001;
    double mdot_hi = mdot_stoich * 100.0;

    auto f = [&](double mdot) { return compute_O2_dry_burned_ox(mdot, fuel, oxidizer); };
    double mdot_result = bisection_solve(f, X_O2_dry_target, mdot_lo, mdot_hi, tol, max_iter);

    Stream result = oxidizer;
    result.mdot = mdot_result;
    return result;
}

Stream set_oxidizer_stream_for_CO2(double X_CO2_target, const Stream& fuel, const Stream& oxidizer,
                                    double tol, std::size_t max_iter) {
    if (fuel.mdot <= 0.0) {
        throw std::invalid_argument("set_oxidizer_stream_for_CO2: fuel.mdot must be positive");
    }
    if (X_CO2_target <= 0.0) {
        throw std::invalid_argument("set_oxidizer_stream_for_CO2: X_CO2_target must be positive");
    }

    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

    double mdot_lo = mdot_stoich * 1.001;  // just above stoich -> high CO2
    double mdot_hi = mdot_stoich * 100.0;  // very lean -> low CO2

    // CO2 decreases with increasing oxidizer (dilution)
    auto f = [&](double mdot) { return compute_CO2_burned_ox(mdot, fuel, oxidizer); };
    double mdot_result = bisection_solve(f, X_CO2_target, mdot_lo, mdot_hi, tol, max_iter);

    Stream result = oxidizer;
    result.mdot = mdot_result;
    return result;
}

Stream set_oxidizer_stream_for_CO2_dry(double X_CO2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                        double tol, std::size_t max_iter) {
    if (fuel.mdot <= 0.0) {
        throw std::invalid_argument("set_oxidizer_stream_for_CO2_dry: fuel.mdot must be positive");
    }
    if (X_CO2_dry_target <= 0.0) {
        throw std::invalid_argument("set_oxidizer_stream_for_CO2_dry: X_CO2_dry_target must be positive");
    }

    double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
    double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

    double mdot_lo = mdot_stoich * 1.001;
    double mdot_hi = mdot_stoich * 100.0;

    auto f = [&](double mdot) { return compute_CO2_dry_burned_ox(mdot, fuel, oxidizer); };
    double mdot_result = bisection_solve(f, X_CO2_dry_target, mdot_lo, mdot_hi, tol, max_iter);

    Stream result = oxidizer;
    result.mdot = mdot_result;
    return result;
}

// -------------------------------------------------------------
// State-based combustion functions
// -------------------------------------------------------------

// Complete combustion (isothermal) - keeps input temperature
State complete_combustion_isothermal(const State& in) {
    State out;
    out.T = in.T;
    out.P = in.P;
    out.X = complete_combustion_to_CO2_H2O(in.X);
    return out;
}

// Complete combustion (adiabatic) - solves for flame temperature
State complete_combustion(const State& in) {
    // Get burned composition at constant T first
    std::vector<double> X_burned = complete_combustion_to_CO2_H2O(in.X);

    // Calculate inlet enthalpy
    double H_in = h(in.T, in.X);

    // Solve for adiabatic flame temperature: h(T_ad, X_burned) = H_in
    double T_ad = calc_T_from_h(H_in, X_burned, in.T);

    State out;
    out.T = T_ad;
    out.P = in.P;
    out.X = X_burned;
    return out;
}
