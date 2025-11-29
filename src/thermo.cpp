#include "../include/thermo.h"
#include <cmath>
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

// Derivative of dimensionless Gibbs free energy with respect to temperature
double dg_over_RT_dT(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double dg_dT_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            // d(G/RT)/dT = d(H/RT)/dT - d(S/R)/dT
            // d(H/RT)/dT = -H/(R*T^2) + (1/RT)*dH/dT = -H/(R*T^2) + Cp/(R*T)
            double h_RT_val = h_RT(i, T);
            double cp_R_val = cp_R(i, T);
            double dh_RT_dT = -h_RT_val / T + cp_R_val / T;

            // d(S/R)/dT from NASA polynomial
            const NASA_Coeffs& nasa = nasa_coeffs[i];
            const std::vector<double>* coeffs;
            if (T < nasa.T_mid) {
                coeffs = &nasa.low_coeffs;
            } else {
                coeffs = &nasa.high_coeffs;
            }
            double ds_R_dT = (*coeffs)[0] / T + (*coeffs)[1] + (*coeffs)[2] * T +
                             (*coeffs)[3] * T * T / 2.0 + (*coeffs)[4] * T * T * T / 3.0;

            dg_dT_mix += X[i] * (dh_RT_dT - ds_R_dT);
        }
    }

    return dg_dT_mix;
}

// Calculate temperature from enthalpy using Newton's method
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::invalid_argument("Mole fractions do not sum to 1.0");
    }

    double T_min = 300.0;
    double T_max = 5000.0;

    if (T_guess < T_min) T_guess = T_min;
    if (T_guess > T_max) T_guess = T_max;

    double T = T_guess;

    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        double h_val = h(T, X);
        double dh_dT_val = dh_dT(T, X);

        if (std::abs(dh_dT_val) < 1.0e-10) {
            const double factor = (h_val < h_target) ? 1.1 : 0.9;
            T *= factor;
            if (T < T_min) T = T_min;
            if (T > T_max) T = T_max;
            continue;
        }

        double delta_T = (h_target - h_val) / dh_dT_val;
        double max_step = 0.1 * T;
        if (std::abs(delta_T) > max_step) {
            delta_T = (delta_T > 0) ? max_step : -max_step;
        }

        T += delta_T;
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;

        if (std::abs(delta_T) < tol || std::abs(h_val - h_target) < tol) {
            return T;
        }
    }

    return T;
}

// Calculate temperature from entropy using Newton's method
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::invalid_argument("Mole fractions do not sum to 1.0");
    }

    if (P <= 0.0) {
        throw std::invalid_argument("Pressure must be positive");
    }

    double T_min = 300.0;
    double T_max = 5000.0;

    if (T_guess < T_min) T_guess = T_min;
    if (T_guess > T_max) T_guess = T_max;

    double T = T_guess;

    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        double s_val = s(T, X, P);
        double ds_dT_val = ds_dT(T, X);

        if (std::abs(ds_dT_val) < 1.0e-10) {
            const double factor = (s_val < s_target) ? 1.1 : 0.9;
            T *= factor;
            if (T < T_min) T = T_min;
            if (T > T_max) T = T_max;
            continue;
        }

        double delta_T = (s_target - s_val) / ds_dT_val;
        double max_step = 0.1 * T;
        if (std::abs(delta_T) > max_step) {
            delta_T = (delta_T > 0) ? max_step : -max_step;
        }

        T += delta_T;
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;

        if (std::abs(delta_T) < tol || std::abs(s_val - s_target) < tol) {
            return T;
        }
    }

    return T;
}

// Calculate temperature from heat capacity using Newton's method
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::invalid_argument("Mole fractions do not sum to 1.0");
    }

    double T_min = 300.0;
    double T_max = 5000.0;

    if (T_guess < T_min) T_guess = T_min;
    if (T_guess > T_max) T_guess = T_max;

    double T = T_guess;

    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        double cp_val = cp(T, X);
        double dcp_dT_val = dcp_dT(T, X);

        if (std::abs(dcp_dT_val) < 1.0e-10) {
            const double factor = (cp_val < cp_target) ? 1.1 : 0.9;
            T *= factor;
            if (T < T_min) T = T_min;
            if (T > T_max) T = T_max;
            continue;
        }

        double delta_T = (cp_target - cp_val) / dcp_dT_val;
        double max_step = 0.1 * T;
        if (std::abs(delta_T) > max_step) {
            delta_T = (delta_T > 0) ? max_step : -max_step;
        }

        T += delta_T;
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;

        if (std::abs(delta_T) < tol || std::abs(cp_val - cp_target) < tol) {
            return T;
        }
    }

    return T;
}

// Normalize a vector of fractions to sum to 1.0
std::vector<double> normalize_fractions(const std::vector<double>& fractions) {
    double sum = std::accumulate(fractions.begin(), fractions.end(), 0.0);
    std::vector<double> normalized = fractions;

    if (std::abs(sum) < 1.0e-10) {
        std::cerr << "Warning: normalize_fractions received all zeros. Returning all zeros." << std::endl;
        return normalized;
    }

    for (double& value : normalized) {
        value /= sum;
    }

    return normalized;
}

// Convert mole fractions to dry fractions (remove water vapor and normalize)
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

// Print all properties of a mixture at given temperature and pressure
void print_mixture_properties(double T, double P, const std::vector<double>& X) {
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        std::cout << "Warning: Mole fractions sum to " << sum << ", not 1.0" << std::endl;
    }

    std::cout << "\nMixture Composition:" << std::endl;
    std::cout << "-------------------" << std::endl;
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            std::cout << species_name(i) << ": " << X[i] << std::endl;
        }
    }

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

    std::cout << "\nDerivatives:" << std::endl;
    std::cout << "-----------" << std::endl;
    std::cout << "dh/dT: " << dh_dT(T, X) << " J/(mol·K)" << std::endl;
    std::cout << "ds/dT: " << ds_dT(T, X) << " J/(mol·K²)" << std::endl;
    std::cout << "dCp/dT: " << dcp_dT(T, X) << " J/(mol·K²)" << std::endl;
}
