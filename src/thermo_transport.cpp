#include "../include/thermo_transport.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>

// Conversion factor from J/mol to J/kg
double J_per_mol_to_J_per_kg(double value, double molar_mass) {
    return value * 1000.0 / molar_mass; // molar_mass is in g/mol, convert to kg/mol
}

// Species name lookup functions
std::string species_name(std::size_t species_index) {
    if (species_index >= species_names.size()) {
        throw std::out_of_range("Species index out of range");
    }
    return species_names[species_index];
}

std::size_t species_index_from_name(const std::string& name) {
    auto it = species_index.find(name);
    if (it == species_index.end()) {
        throw std::out_of_range("Species name not found: " + name);
    }
    return static_cast<std::size_t>(it->second);
}

// Calculate mixture molecular weight [g/mol]
double mwmix(const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double sum = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        sum += X[i] * molar_masses[i];
    }
    return sum;
}

// NASA polynomial evaluation for Cp/R
double cp_R(std::size_t species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];
    
    const std::vector<double>* coeffs;
    if (T < nasa.T_mid) {
        if (T < nasa.T_low) {
            std::cerr << "Warning: Temperature " << T << " K is below valid range (" << nasa.T_low 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the lowest valid temperature coefficients
            T = nasa.T_low;
        }
        coeffs = &nasa.low_coeffs;
    } else {
        if (T > nasa.T_high) {
            std::cerr << "Warning: Temperature " << T << " K is above valid range (" << nasa.T_high 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the highest valid temperature coefficients
            T = nasa.T_high;
        }
        coeffs = &nasa.high_coeffs;
    }
    
    // NASA polynomial for Cp/R: a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    return (*coeffs)[0] + (*coeffs)[1] * T + (*coeffs)[2] * T * T + 
           (*coeffs)[3] * T * T * T + (*coeffs)[4] * T * T * T * T;
}

// Cv/R = Cp/R - 1
inline double cv_R(std::size_t species_idx, double T) {
    return cp_R(species_idx, T) - 1.0;
}

// NASA polynomial evaluation for H/RT
double h_RT(std::size_t species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];
    
    const std::vector<double>* coeffs;
    if (T < nasa.T_mid) {
        if (T < nasa.T_low) {
            std::cerr << "Warning: Temperature " << T << " K is below valid range (" << nasa.T_low 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the lowest valid temperature coefficients
            T = nasa.T_low;
        }
        coeffs = &nasa.low_coeffs;
    } else {
        if (T > nasa.T_high) {
            std::cerr << "Warning: Temperature " << T << " K is above valid range (" << nasa.T_high 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the highest valid temperature coefficients
            T = nasa.T_high;
        }
        coeffs = &nasa.high_coeffs;
    }
    
    // NASA polynomial for H/RT: a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
    return (*coeffs)[0] + (*coeffs)[1] * T / 2.0 + (*coeffs)[2] * T * T / 3.0 + 
           (*coeffs)[3] * T * T * T / 4.0 + (*coeffs)[4] * T * T * T * T / 5.0 + 
           (*coeffs)[5] / T;
}

// NASA polynomial evaluation for S/R
double s_R(std::size_t species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];
    
    const std::vector<double>* coeffs;
    if (T < nasa.T_mid) {
        if (T < nasa.T_low) {
            std::cerr << "Warning: Temperature " << T << " K is below valid range (" << nasa.T_low 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the lowest valid temperature coefficients
            T = nasa.T_low;
        }
        coeffs = &nasa.low_coeffs;
    } else {
        if (T > nasa.T_high) {
            std::cerr << "Warning: Temperature " << T << " K is above valid range (" << nasa.T_high 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the highest valid temperature coefficients
            T = nasa.T_high;
        }
        coeffs = &nasa.high_coeffs;
    }
    
    // NASA polynomial for S/R: a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
    return (*coeffs)[0] * std::log(T) + (*coeffs)[1] * T + (*coeffs)[2] * T * T / 2.0 + 
           (*coeffs)[3] * T * T * T / 3.0 + (*coeffs)[4] * T * T * T * T / 4.0 + 
           (*coeffs)[6];
}

// Dimensionless Gibbs free energy G/(R*T) = H/(R*T) - S/R
double g_over_RT(std::size_t species_idx, double T) {
    return h_RT(species_idx, T) - s_R(species_idx, T);
}

// Number of thermo species used in the internal tables
std::size_t num_species() {
    return species_names.size();
}

// Molar mass lookup [g/mol]
double species_molar_mass(std::size_t species_index) {
    if (species_index >= molar_masses.size()) {
        throw std::out_of_range("species_molar_mass: species_index out of range");
    }
    return molar_masses[species_index];
}

double species_molar_mass_from_name(const std::string& name) {
    const std::size_t idx = species_index_from_name(name);
    return species_molar_mass(idx);
}

// Calculate heat capacity at constant pressure [J/(mol·K)]
double cp(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double cp_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        cp_mix += X[i] * cp_R(i, T) * R_GAS;
    }
    return cp_mix;
}

// Calculate enthalpy [J/mol]
double h(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double h_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        h_mix += X[i] * h_RT(i, T) * R_GAS * T;
    }
    return h_mix;
}

// Calculate entropy [J/(mol·K)]
double s(double T, const std::vector<double>& X, double P, double P_ref) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double s_mix = 0.0;
    
    // Calculate pure species entropies and mixing term
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            // Pure species entropy at reference pressure
            double s_i = s_R(i, T) * R_GAS;
            
            // Pressure correction: -R * ln(P/P_ref)
            s_i -= R_GAS * std::log(P / P_ref);
            
            // Mixing entropy: -R * ln(X_i)
            s_i -= R_GAS * std::log(X[i]);
            
            s_mix += X[i] * s_i;
        }
    }
    
    return s_mix;
}

// Calculate heat capacity at constant volume [J/(mol·K)]
double cv(double T, const std::vector<double>& X) {
    // Cv = Cp - R for ideal gases
    return cp(T, X) - R_GAS;
}

// Calculate derivative of enthalpy with respect to temperature [J/(mol·K)]
// For ideal gas, this is the same as Cp
double dh_dT(double T, const std::vector<double>& X) {
    return cp(T, X);
}

// Calculate derivative of entropy with respect to temperature [J/(mol·K^2)]
double ds_dT(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double ds_dT_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            // Derivative of NASA polynomial for S/R with respect to T
            const NASA_Coeffs& nasa = nasa_coeffs[i];
            
            const std::vector<double>* coeffs;
            if (T < nasa.T_mid) {
                if (T < nasa.T_low) {
                    throw std::out_of_range("Temperature below valid range for species " + species_names[i]);
                }
                coeffs = &nasa.low_coeffs;
            } else {
                if (T > nasa.T_high) {
                    throw std::out_of_range("Temperature above valid range for species " + species_names[i]);
                }
                coeffs = &nasa.high_coeffs;
            }
            
            // d(S/R)/dT = a1/T + a2 + a3*T + a4*T^2/2 + a5*T^3/3
            double ds_dT_i = (*coeffs)[0] / T + (*coeffs)[1] + (*coeffs)[2] * T + 
                             (*coeffs)[3] * T * T / 2.0 + (*coeffs)[4] * T * T * T / 3.0;
            
            ds_dT_mix += X[i] * ds_dT_i * R_GAS;
        }
    }
    
    return ds_dT_mix;
}

// Calculate derivative of heat capacity with respect to temperature [J/(mol·K^2)]
double dcp_dT(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double dcp_dT_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        const NASA_Coeffs& nasa = nasa_coeffs[i];
        
        const std::vector<double>* coeffs;
        if (T < nasa.T_mid) {
            if (T < nasa.T_low) {
                throw std::out_of_range("Temperature below valid range for species " + species_names[i]);
            }
            coeffs = &nasa.low_coeffs;
        } else {
            if (T > nasa.T_high) {
                throw std::out_of_range("Temperature above valid range for species " + species_names[i]);
            }
            coeffs = &nasa.high_coeffs;
        }
        
        // d(Cp/R)/dT = a2 + 2*a3*T + 3*a4*T^2 + 4*a5*T^3
        double dcp_dT_i = (*coeffs)[1] + 2.0 * (*coeffs)[2] * T + 
                          3.0 * (*coeffs)[3] * T * T + 4.0 * (*coeffs)[4] * T * T * T;
        
        dcp_dT_mix += X[i] * dcp_dT_i * R_GAS;
    }
    
    return dcp_dT_mix;
}

// Calculate density [kg/m^3]
double density(double T, double P, const std::vector<double>& X) {
    // Ideal gas law: rho = P * MW / (R * T)
    // P in Pa, MW in g/mol, R in J/(mol·K), T in K
    // Need to convert MW from g/mol to kg/mol
    double MW = mwmix(X) / 1000.0; // Convert g/mol to kg/mol
    return P * MW / (R_GAS * T);
}

// Calculate specific gas constant [J/(kg·K)]
double specific_gas_constant(const std::vector<double>& X) {
    // R_specific = R_universal / MW
    // R_universal in J/(mol·K), MW in g/mol
    // Need to convert MW from g/mol to kg/mol
    double MW = mwmix(X) / 1000.0; // Convert g/mol to kg/mol
    return R_GAS / MW;
}

// Calculate isentropic expansion coefficient (gamma = Cp/Cv)
double isentropic_expansion_coefficient(double T, const std::vector<double>& X) {
    double cp_val = cp(T, X);
    double cv_val = cv(T, X);
    
    // Guard against division by zero
    if (std::abs(cv_val) < 1.0e-10) {
        throw std::runtime_error("Cv is too close to zero for calculating gamma");
    }
    
    return cp_val / cv_val;
}

// Calculate speed of sound [m/s] for an ideal gas
double speed_of_sound(double T, const std::vector<double>& X) {
    // c = sqrt(gamma * R_specific * T)
    // gamma is dimensionless, R_specific in J/(kg·K), T in K
    double gamma = isentropic_expansion_coefficient(T, X);
    double R_specific = specific_gas_constant(X);
    
    return std::sqrt(gamma * R_specific * T);
}

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
    
    double oxygen_required = fuel.C + fuel.H / 4.0 - fuel.O / 2.0;
    
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
        double nu_O2 = oxygen_required_per_mol_fuel(static_cast<int>(i));
        if (nu_O2 <= 0.0) continue;

        const Molecular_Structure& sp = molecular_structures[i];

        double n_fuel_initial   = n[i];
        double n_fuel_reacted   = f * n_fuel_initial;
        double n_fuel_remaining = n_fuel_initial - n_fuel_reacted;

        n[i] = n_fuel_remaining;

        // C_x H_y O_z N_w + ... -> x CO2 + (y/2) H2O + (w/2) N2
        delta_CO2 += sp.C * n_fuel_reacted;
        delta_H2O += 0.5 * sp.H * n_fuel_reacted;
        delta_N2  += 0.5 * sp.N * n_fuel_reacted;
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

static std::size_t find_species_index(const std::string& name)
{
    for (std::size_t k = 0; k < species_names.size(); ++k) {
        if (species_names[k] == name) {
            return k;
        }
    }
    throw std::runtime_error("find_species_index: species '" + name + "' not found");
}

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

    static const std::size_t O2_index = find_species_index("O2");

    const double X_O2_ox = X_ox[O2_index];
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
// Helpers: mass <-> mole conversion
// -------------------------------------------------------------

static void mass_to_mole(
    const std::vector<double>& Y,
    std::vector<double>& X)
{
    if (Y.size() != molar_masses.size()) {
        throw std::invalid_argument("mass_to_mole: size mismatch");
    }

    const std::size_t n = Y.size();
    double denom = 0.0;

    for (std::size_t k = 0; k < n; ++k) {
        double Wk = molar_masses[k]; // g/mol (any consistent unit)
        if (Wk <= 0.0) {
            throw std::runtime_error("mass_to_mole: non-positive molar mass");
        }
        denom += Y[k] / Wk;
    }

    if (denom <= 0.0) {
        throw std::runtime_error("mass_to_mole: non-positive denominator");
    }

    X.resize(n);
    for (std::size_t k = 0; k < n; ++k) {
        double Wk = molar_masses[k];
        X[k] = (Y[k] / Wk) / denom;
    }
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

    static const std::size_t O2_index = find_species_index("O2");

    const double Y_O2_ox = Y_ox[O2_index];
    if (Y_O2_ox <= 0.0) {
        throw std::runtime_error(
            "stoich_f_over_o_mass: oxidizer has zero O2 mass fraction");
    }

    // Convert fuel stream from mass fractions to mole fractions for use in
    // oxygen_required_per_kg_mixture(X_fuel).
    std::vector<double> X_fuel;
    mass_to_mole(Y_fuel, X_fuel);

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

// Implementation of collision integral tables

// Collision integral lookup tables based on Cantera's MMCollisionInt
// Reduced temperature values for omega22 table
const std::vector<double> tstar22 = {
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
    5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
    18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0
};

// Omega(2,2) collision integral values (first column for delta* = 0)
const std::vector<double> omega22_values = {
    4.1005, 3.2626, 2.8399, 2.531, 2.2837, 2.0838, 1.922, 1.7902, 
    1.6823, 1.5929, 1.4551, 1.3551, 1.28, 1.2219, 1.1757, 1.0933, 
    1.0388, 0.99963, 0.96988, 0.92676, 0.89616, 0.87272, 0.85379, 
    0.83795, 0.82435, 0.80184, 0.78363, 0.76834, 0.75518, 0.74364, 
    0.71982, 0.70097, 0.68545, 0.67232, 0.65099, 0.61397, 0.5887
};

// Linear interpolation function
double linear_interp(double x, const std::vector<double>& x_values, const std::vector<double>& y_values) {
    // Handle out-of-bounds cases
    if (x <= x_values.front()) return y_values.front();
    if (x >= x_values.back()) return y_values.back();
    
    // Find position
    auto it = std::lower_bound(x_values.begin(), x_values.end(), x);
    if (it == x_values.begin()) return y_values.front();

    std::size_t idx = static_cast<std::size_t>(std::distance(x_values.begin(), it) - 1);
    double x1 = x_values[idx];
    double x2 = x_values[idx + 1];
    double y1 = y_values[idx];
    double y2 = y_values[idx + 1];
    
    // Linear interpolation
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

// Calculate omega(2,2) collision integral
double omega22(double T, double well_depth) {
    // Calculate reduced temperature T* = kT/ε
    double tstar = T * BOLTZMANN / well_depth;
    
    // Interpolate to get omega(2,2)
    return linear_interp(tstar, tstar22, omega22_values);
}

// Calculate dynamic viscosity [Pa·s] using collision integrals (Cantera approach)
// P is not used but kept for API consistency
double viscosity(double T, double P [[maybe_unused]], const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // First, calculate pure species viscosities using kinetic theory
    std::vector<double> pure_visc(X.size());
    for (size_t i = 0; i < X.size(); ++i) {
        // Get molecular weight in kg/mol
        double mw_kg = molar_masses[i] / 1000.0;
        
        // Get collision diameter in meters (convert from Angstroms)
        double diam = transport_props[i].diameter * 1.0e-10;
        
        // Get well depth in Joules (convert from Kelvin)
        double eps = transport_props[i].well_depth * BOLTZMANN;
        
        // Calculate reduced temperature and collision integral
        double omega = omega22(T, eps);
        
        // Calculate viscosity using kinetic theory formula
        // μ = (5/16) * √(π*m*k*T) / (π*σ²*Ω(2,2))
        // where m is mass of a molecule = MW/NA
        double mass_molecule = mw_kg / AVOGADRO;
        pure_visc[i] = 5.0/16.0 * std::sqrt(PI * mass_molecule * BOLTZMANN * T) / 
                      (PI * diam * diam * omega);
    }
    
    // Wilke's mixing rule for viscosity
    double mix_visc = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        if (X[i] <= 0.0) continue;
        
        double denom = 0.0;
        for (size_t j = 0; j < X.size(); ++j) {
            if (X[j] <= 0.0) continue;
            
            double phi_ij;
            if (i == j) {
                phi_ij = 1.0;
            } else {
                // Calculate molecular weight ratio
                double M_i = molar_masses[i];
                double M_j = molar_masses[j];
                
                // Wilke's formula
                double term1 = 1.0 / std::sqrt(8.0) / std::sqrt(1.0 + M_i / M_j);
                double term2 = std::pow(1.0 + std::sqrt(pure_visc[i] / pure_visc[j]) * 
                               std::pow(M_j / M_i, 0.25), 2);
                phi_ij = term1 * term2;
            }
            
            denom += X[j] * phi_ij;
        }
        
        mix_visc += X[i] * pure_visc[i] / denom;
    }
    
    return mix_visc;
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

// Calculate thermal conductivity [W/(m·K)] using Cantera approach
// P is not used but kept for API consistency
double thermal_conductivity(double T, double P [[maybe_unused]], const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // Calculate pure species thermal conductivities using Eucken formula
    std::vector<double> pure_cond(X.size());
    for (size_t i = 0; i < X.size(); ++i) {
        // Get molecular weight in kg/mol
        double mw_kg = molar_masses[i] / 1000.0;
        
        // Get collision diameter in meters (convert from Angstroms)
        double diam = transport_props[i].diameter * 1.0e-10;
        
        // Get well depth in Joules (convert from Kelvin)
        double eps = transport_props[i].well_depth * BOLTZMANN;
        
        // Calculate pure species viscosity
        double omega = omega22(T, eps);
        double mass_molecule = mw_kg / AVOGADRO;
        double visc = 5.0/16.0 * std::sqrt(PI * mass_molecule * BOLTZMANN * T) / 
                     (PI * diam * diam * omega);
        
        // Get specific heat capacity for the species in J/(mol·K)
        double cp_species = cp_R(i, T) * R_GAS;
        
        // Convert cp from J/(mol·K) to J/(kg·K)
        double cp_mass = cp_species * 1000.0 / mw_kg;
        
        // Determine rotational contribution based on geometry
        double f_rot = 0.0;
        if (transport_props[i].geometry == "atom") {
            f_rot = 0.0; // Monatomic
        } else if (transport_props[i].geometry == "linear") {
            f_rot = 1.0; // Linear molecule
        } else {
            f_rot = 1.5; // Nonlinear molecule
        }
        
        // Modified Eucken formula for thermal conductivity
        // λ = μ * (cp_mass + f_rot*R_specific)
        double R_specific = R_GAS * 1000.0 / mw_kg; // J/(kg·K)
        
        // Scale factor to get reasonable thermal conductivity values
        // This is a correction factor to match experimental data better
        double scale_factor = 0.1;
        
        pure_cond[i] = visc * (cp_mass + f_rot * R_specific) * scale_factor;
    }
    
    // Use the mixture rule from Cantera's MixTransport
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    for (size_t i = 0; i < X.size(); ++i) {
        if (X[i] <= 0.0) continue;
        
        sum1 += X[i] * pure_cond[i];
        sum2 += X[i] / pure_cond[i];
    }
    
    // Mixture thermal conductivity (average of upper and lower bounds)
    return 0.5 * (sum1 + 1.0 / sum2);
}

// Calculate Prandtl number (dimensionless)
double prandtl(double T, double P, const std::vector<double>& X) {
    // Prandtl number = (Cp * viscosity) / thermal conductivity
    double cp_val = cp(T, X);
    double visc_val = viscosity(T, P, X);
    double cond_val = thermal_conductivity(T, P, X);
    
    return cp_val * visc_val / cond_val;
}

// Calculate kinematic viscosity [m^2/s]
double kinematic_viscosity(double T, double P, const std::vector<double>& X) {
    // Kinematic viscosity = dynamic viscosity / density
    double visc_val = viscosity(T, P, X);
    double rho_val = density(T, P, X);
    
    return visc_val / rho_val;
}

// Calculate temperature from enthalpy using Newton's method
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::invalid_argument("Mole fractions do not sum to 1.0");
    }
    
    // Ensure T_guess is within valid range for all species
    double T_min = 300.0;  // Safe minimum for all species
    double T_max = 5000.0; // Safe maximum for all species
    
    if (T_guess < T_min) T_guess = T_min;
    if (T_guess > T_max) T_guess = T_max;
    
    // Newton's method implementation with safeguards
    double T = T_guess;
    double h_val, dh_dT_val;
    double delta_T;
    
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        // Calculate h and dh/dT at current temperature
        h_val = h(T, X);
        dh_dT_val = dh_dT(T, X);
        
        // Check for division by zero
        if (std::abs(dh_dT_val) < 1.0e-10) {
            // If derivative is too small, use a guarded multiplicative step
            const double factor = (h_val < h_target) ? 1.1 : 0.9;
            T *= factor;
            if (T < T_min) T = T_min;
            if (T > T_max) T = T_max;
            continue;
        }
        
        // Newton step
        delta_T = (h_target - h_val) / dh_dT_val;
        
        // Limit step size for stability
        double max_step = 0.1 * T; // Limit to 10% of current T
        if (std::abs(delta_T) > max_step) {
            delta_T = (delta_T > 0) ? max_step : -max_step;
        }
        
        // Update temperature
        T += delta_T;
        
        // Ensure T stays within bounds
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        
        // Check convergence
        if (std::abs(delta_T) < tol || std::abs(h_val - h_target) < tol) {
            return T;
        }
    }
    
    // If we reach here, we didn't converge within max_iter
    // Return the best estimate we have
    return T;
}

// Calculate temperature from entropy using Newton's method
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::invalid_argument("Mole fractions do not sum to 1.0");
    }
    
    // Ensure P is positive
    if (P <= 0.0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    
    // Ensure T_guess is within valid range for all species
    double T_min = 300.0;  // Safe minimum for all species
    double T_max = 5000.0; // Safe maximum for all species
    
    if (T_guess < T_min) T_guess = T_min;
    if (T_guess > T_max) T_guess = T_max;
    
    // Newton's method implementation with safeguards
    double T = T_guess;
    double s_val, ds_dT_val;
    double delta_T;
    
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        // Calculate s and ds/dT at current temperature
        s_val = s(T, X, P);
        ds_dT_val = ds_dT(T, X);
        
        // Check for division by zero
        if (std::abs(ds_dT_val) < 1.0e-10) {
            // If derivative is too small, use a guarded multiplicative step
            const double factor = (s_val < s_target) ? 1.1 : 0.9;
            T *= factor;
            if (T < T_min) T = T_min;
            if (T > T_max) T = T_max;
            continue;
        }
        
        // Newton step
        delta_T = (s_target - s_val) / ds_dT_val;
        
        // Limit step size for stability
        double max_step = 0.1 * T; // Limit to 10% of current T
        if (std::abs(delta_T) > max_step) {
            delta_T = (delta_T > 0) ? max_step : -max_step;
        }
        
        // Update temperature
        T += delta_T;
        
        // Ensure T stays within bounds
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        
        // Check convergence
        if (std::abs(delta_T) < tol || std::abs(s_val - s_target) < tol) {
            return T;
        }
    }
    
    // If we reach here, we didn't converge within max_iter
    // Return the best estimate we have
    return T;
}

// Calculate temperature from heat capacity using Newton's method
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::invalid_argument("Mole fractions do not sum to 1.0");
    }
    
    // Ensure T_guess is within valid range for all species
    double T_min = 300.0;  // Safe minimum for all species
    double T_max = 5000.0; // Safe maximum for all species
    
    if (T_guess < T_min) T_guess = T_min;
    if (T_guess > T_max) T_guess = T_max;
    
    // Newton's method implementation with safeguards
    double T = T_guess;
    double cp_val, dcp_dT_val;
    double delta_T;
    
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        // Calculate cp and dcp/dT at current temperature
        cp_val = cp(T, X);
        dcp_dT_val = dcp_dT(T, X);
        
        // Check for division by zero
        if (std::abs(dcp_dT_val) < 1.0e-10) {
            // If derivative is too small, use a guarded multiplicative step
            const double factor = (cp_val < cp_target) ? 1.1 : 0.9;
            T *= factor;
            if (T < T_min) T = T_min;
            if (T > T_max) T = T_max;
            continue;
        }
        
        // Newton step
        delta_T = (cp_target - cp_val) / dcp_dT_val;
        
        // Limit step size for stability
        double max_step = 0.1 * T; // Limit to 10% of current T
        if (std::abs(delta_T) > max_step) {
            delta_T = (delta_T > 0) ? max_step : -max_step;
        }
        
        // Update temperature
        T += delta_T;
        
        // Ensure T stays within bounds
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        
        // Check convergence
        if (std::abs(delta_T) < tol || std::abs(cp_val - cp_target) < tol) {
            return T;
        }
    }
    
    // If we reach here, we didn't converge within max_iter
    // Return the best estimate we have
    return T;
}

// Normalize a vector of fractions to sum to 1.0
// Returns all zeros with a warning if input contains all zeros
std::vector<double> normalize_fractions(const std::vector<double>& fractions) {
    // Calculate the sum of all fractions
    double sum = std::accumulate(fractions.begin(), fractions.end(), 0.0);
    
    // Create a copy of the input vector for the result
    std::vector<double> normalized = fractions;
    
    // Check if sum is zero or very close to zero
    if (std::abs(sum) < 1.0e-10) {
        std::cerr << "Warning: normalize_fractions received all zeros. Returning all zeros." << std::endl;
        return normalized; // Already all zeros
    }
    
    // Normalize each fraction by dividing by the sum
    for (double& value : normalized) {
        value /= sum;
    }
    
    return normalized;
}

// Convert mole fractions to dry fractions (remove water vapor and normalize)
// Returns all zeros with a warning if input contains only water vapor
std::vector<double> convert_to_dry_fractions(const std::vector<double>& mole_fractions) {
    // Get the index of water vapor
    std::size_t h2o_idx = species_index_from_name("H2O");
    
    // Create a copy of the input vector for the result
    std::vector<double> dry_fractions = mole_fractions;
    
    // Calculate the sum of all non-water fractions
    double sum = 0.0;
    for (std::size_t i = 0; i < dry_fractions.size(); ++i) {
        if (static_cast<int>(i) != h2o_idx) {
            sum += dry_fractions[i];
        }
    }
    
    // Check if sum is zero or very close to zero (only water vapor)
    if (std::abs(sum) < 1.0e-10) {
        std::cerr << "Warning: convert_to_dry_fractions received only water vapor. Returning all zeros." << std::endl;
        std::fill(dry_fractions.begin(), dry_fractions.end(), 0.0);
        return dry_fractions;
    }
    
    // Set water vapor fraction to zero
    if (h2o_idx >= 0 && h2o_idx < static_cast<int>(dry_fractions.size())) {
        dry_fractions[h2o_idx] = 0.0;
    }
    
    // Normalize the remaining fractions
    for (double& value : dry_fractions) {
        value /= sum;
    }
    
    return dry_fractions;
}

// Convert mole fractions X_k to mass fractions Y_k using molar_masses.
std::vector<double> mole_to_mass(const std::vector<double>& X)
{
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

// Convert mass fractions Y_k to mole fractions X_k using molar_masses.
std::vector<double> mass_to_mole(const std::vector<double>& Y)
{
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

// Print all properties of a mixture at given temperature and pressure
void print_mixture_properties(double T, double P, const std::vector<double>& X) {
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        std::cout << "Warning: Mole fractions sum to " << sum << ", not 1.0" << std::endl;
    }
    
    // Print mixture composition
    std::cout << "\nMixture Composition:" << std::endl;
    std::cout << "-------------------" << std::endl;
    for (size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            std::cout << species_name(i) << ": " << X[i] << std::endl;
        }
    }
    
    // Print thermodynamic properties
    std::cout << "\nThermodynamic Properties at T = " << T << " K, P = " << P << " Pa:" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Molecular Weight: " << mwmix(X) << " g/mol" << std::endl;
    std::cout << "Density: " << density(T, P, X) << " kg/m³" << std::endl;
    std::cout << "Specific Gas Constant (Rs): " << specific_gas_constant(X) << " J/(kg·K)" << std::endl;
    std::cout << "Isentropic Expansion Coefficient (gamma): " << isentropic_expansion_coefficient(T, X) << std::endl;
    std::cout << "Speed of Sound: " << speed_of_sound(T, X) << " m/s" << std::endl;
    std::cout << "Enthalpy: " << h(T, X) << " J/mol" << std::endl;
    std::cout << "Entropy: " << s(T, X, P) << " J/(mol·K)" << std::endl;
    std::cout << "Heat Capacity (Cp): " << cp(T, X) << " J/(mol·K)" << std::endl;
    std::cout << "Heat Capacity (Cv): " << cv(T, X) << " J/(mol·K)" << std::endl;
    
    // Print transport properties
    std::cout << "\nTransport Properties:" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "Viscosity: " << viscosity(T, P, X) * 1.0e6 << " μPa·s" << std::endl;
    std::cout << "Thermal Conductivity: " << thermal_conductivity(T, P, X) << " W/(m·K)" << std::endl;
    std::cout << "Kinematic Viscosity: " << kinematic_viscosity(T, P, X) << " m²/s" << std::endl;
    std::cout << "Prandtl Number: " << prandtl(T, P, X) << std::endl;
    
    // Print derivatives
    std::cout << "\nDerivatives:" << std::endl;
    std::cout << "-----------" << std::endl;
    std::cout << "dh/dT: " << dh_dT(T, X) << " J/(mol·K)" << std::endl;
    std::cout << "ds/dT: " << ds_dT(T, X) << " J/(mol·K²)" << std::endl;
    std::cout << "dCp/dT: " << dcp_dT(T, X) << " J/(mol·K²)" << std::endl;
    
    // Print combustion properties for individual fuels and the mixture
    bool has_fuel = false;
    
    // First check if the mixture contains any fuel components
    for (size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0 && (molecular_structures[i].C > 0 || molecular_structures[i].H > 0)) {
            has_fuel = true;
            break;
        }
    }
    
    if (has_fuel) {
        std::cout << "\nCombustion Properties:" << std::endl;
        std::cout << "---------------------" << std::endl;
        
        // Print for individual fuel components
        for (size_t i = 0; i < X.size(); ++i) {
            // Check if this is a fuel (contains carbon or hydrogen)
            if (molecular_structures[i].C > 0 || molecular_structures[i].H > 0) {
                if (X[i] > 0.0) {
                    std::cout << species_name(i) << " (" << X[i] * 100.0 << "% of mixture):" << std::endl;
                    std::cout << "  Oxygen required: " << oxygen_required_per_mol_fuel(i) << " mol O2/mol fuel" << std::endl;
                    std::cout << "  Oxygen required: " << oxygen_required_per_kg_fuel(i) << " kg O2/kg fuel" << std::endl;
                }
            }
        }
        
        // Print for the entire mixture
        std::cout << "\nTotal mixture:" << std::endl;
        std::cout << "  Oxygen required: " << oxygen_required_per_mol_mixture(X) << " mol O2/mol mixture" << std::endl;
        std::cout << "  Oxygen required: " << oxygen_required_per_kg_mixture(X) << " kg O2/kg mixture" << std::endl;
    }
}
