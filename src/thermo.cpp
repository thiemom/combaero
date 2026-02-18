#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"
#include "../include/humidair.h"
#include "../include/transport.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

using combaero::thermo::R_GAS;
using combaero::thermo::BOLTZMANN;
using combaero::thermo::AVOGADRO;

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

// Relative extrapolation tolerance beyond NASA-9 range boundaries.
// 5% of the range span is allowed silently (polynomial is smooth at boundaries).
// Beyond that, throws std::out_of_range.
static constexpr double NASA9_EXTRAP_TOL = 0.05;

// Helper: find the correct interval for NASA-9 coefficients
// Returns pointer to the 10-coefficient array for the given temperature.
// Clamps T to the boundary interval for small extrapolations; throws for large ones.
static const std::array<double, 10>* find_nasa9_interval(
    const NASA_Coeffs& nasa, double& T, std::size_t species_idx) {

    const auto& intervals = nasa.intervals;
    if (intervals.empty()) {
        throw std::runtime_error("No NASA-9 intervals for species " + species_names[species_idx]);
    }

    const double T_min    = intervals.front().T_min;
    const double T_max    = intervals.back().T_max;
    const double tol_low  = NASA9_EXTRAP_TOL * T_min;   // 5% of 200 K = 10 K
    const double tol_high = NASA9_EXTRAP_TOL * T_max;   // 5% of 6000 K = 300 K

    if (T < T_min - tol_low) {
        throw std::out_of_range(
            "Temperature " + std::to_string(T) + " K is below extrapolation limit (" +
            std::to_string(T_min - tol_low) + " K) for species " + species_names[species_idx]);
    }
    if (T > T_max + tol_high) {
        throw std::out_of_range(
            "Temperature " + std::to_string(T) + " K is above extrapolation limit (" +
            std::to_string(T_max + tol_high) + " K) for species " + species_names[species_idx]);
    }

    // Clamp to boundary for small extrapolations
    T = std::max(T, T_min);
    T = std::min(T, T_max);

    // Find the interval containing T
    for (const auto& interval : intervals) {
        if (T >= interval.T_min && T <= interval.T_max) {
            return &interval.coeffs;
        }
    }

    // Fallback to last interval (shouldn't happen after clamping)
    return &intervals.back().coeffs;
}

// NASA-9 polynomial for Cp/R: a1/T^2 + a2/T + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
// Uses Horner's algorithm for efficient evaluation
static double cp_R_nasa9(const std::array<double, 10>& a, double T) {
    double T2 = T * T;
    double T_inv = 1.0 / T;
    // Negative powers using Horner: (a1/T + a2)/T = a1/T^2 + a2/T
    double neg_terms = (a[0] * T_inv + a[1]) * T_inv;
    // Positive powers using Horner: a3 + T*(a4 + T*(a5 + T*(a6 + T*a7)))
    double pos_terms = a[2] + T * (a[3] + T * (a[4] + T * (a[5] + T * a[6])));
    return neg_terms + pos_terms;
}

// NASA-9 polynomial for H/RT: -a1/T^2 + a2*ln(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + b1/T
// Uses Horner's algorithm for efficient evaluation
static double h_RT_nasa9(const std::array<double, 10>& a, double T) {
    double T_inv = 1.0 / T;
    double ln_T = std::log(T);
    // Negative powers: -a1/T^2 + a2*ln(T)/T + b1/T
    double neg_terms = (-a[0] * T_inv + a[1] * ln_T + a[8]) * T_inv;
    // Positive powers using Horner: a3 + T*(a4/2 + T*(a5/3 + T*(a6/4 + T*a7/5)))
    double pos_terms = a[2] + T * (a[3] / 2.0 + T * (a[4] / 3.0 + T * (a[5] / 4.0 + T * a[6] / 5.0)));
    return neg_terms + pos_terms;
}

// NASA-9 polynomial for S/R: -a1/(2*T^2) - a2/T + a3*ln(T) + a4*T + a5*T^2/2 + a6*T^3/3 + a7*T^4/4 + b2
// Uses Horner's algorithm for efficient evaluation
static double s_R_nasa9(const std::array<double, 10>& a, double T) {
    double T_inv = 1.0 / T;
    double ln_T = std::log(T);
    // Negative powers: -a1/(2*T^2) - a2/T
    double neg_terms = (-a[0] / 2.0 * T_inv - a[1]) * T_inv;
    // Positive powers using Horner: a3*ln(T) + a4*T + T^2*(a5/2 + T*(a6/3 + T*a7/4))
    double pos_terms = a[2] * ln_T + a[3] * T + T * T * (a[4] / 2.0 + T * (a[5] / 3.0 + T * a[6] / 4.0));
    return neg_terms + pos_terms + a[9];
}

// NASA-9 derivative of Cp/R with respect to T:
// d(Cp/R)/dT = -2*a1/T^3 - a2/T^2 + a4 + 2*a5*T + 3*a6*T^2 + 4*a7*T^3
// Uses Horner's algorithm for efficient evaluation
static double dcp_R_dT_nasa9(const std::array<double, 10>& a, double T) {
    double T_inv = 1.0 / T;
    // Negative powers using Horner: (-2*a1/T - a2)/T^2
    double neg_terms = (-2.0 * a[0] * T_inv - a[1]) * T_inv * T_inv;
    // Positive powers using Horner: a4 + T*(2*a5 + T*(3*a6 + T*4*a7))
    double pos_terms = a[3] + T * (2.0 * a[4] + T * (3.0 * a[5] + T * 4.0 * a[6]));
    return neg_terms + pos_terms;
}

// NASA-9 derivative of S/R with respect to T:
// d(S/R)/dT = a1/T^3 + a2/T^2 + a3/T + a4 + a5*T + a6*T^2 + a7*T^3
// Uses Horner's algorithm for efficient evaluation
static double ds_R_dT_nasa9(const std::array<double, 10>& a, double T) {
    double T_inv = 1.0 / T;
    // Negative powers using Horner: (a1/T + a2)/T^2 + a3/T
    double neg_terms = ((a[0] * T_inv + a[1]) * T_inv + a[2]) * T_inv;
    // Positive powers using Horner: a4 + T*(a5 + T*(a6 + T*a7))
    double pos_terms = a[3] + T * (a[4] + T * (a[5] + T * a[6]));
    return neg_terms + pos_terms;
}

// NASA polynomial evaluation for Cp/R
double cp_R(std::size_t species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];
    const auto* coeffs = find_nasa9_interval(nasa, T, species_idx);
    return cp_R_nasa9(*coeffs, T);
}

// NASA polynomial evaluation for H/RT
double h_RT(std::size_t species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];
    const auto* coeffs = find_nasa9_interval(nasa, T, species_idx);
    return h_RT_nasa9(*coeffs, T);
}

// NASA polynomial evaluation for S/R
double s_R(std::size_t species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];
    const auto* coeffs = find_nasa9_interval(nasa, T, species_idx);
    return s_R_nasa9(*coeffs, T);
}

// Dimensionless Gibbs free energy G/(R*T) = H/(R*T) - S/R
double g_over_RT(std::size_t species_idx, double T) {
    return h_RT(species_idx, T) - s_R(species_idx, T);
}

// Per-species dimensional properties
double cp_species(std::size_t species_idx, double T) {
    return cp_R(species_idx, T) * R_GAS;
}

double h_species(std::size_t species_idx, double T) {
    return h_RT(species_idx, T) * R_GAS * T;
}

double s_species(std::size_t species_idx, double T) {
    return s_R(species_idx, T) * R_GAS;
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

// Calculate internal energy [J/mol]
// For ideal gas: u = h - R*T
double u(double T, const std::vector<double>& X) {
    return h(T, X) - R_GAS * T;
}

// -------------------------------------------------------------
// Mass-specific thermodynamic properties
// -------------------------------------------------------------
// These provide a symmetric interface to the molar-basis functions.
// Conversion: property_mass = property_molar / mwmix * 1000
// (mwmix is in g/mol, so multiply by 1000 to get kg/mol)

// Heat capacity at constant pressure [J/(kg·K)]
double cp_mass(double T, const std::vector<double>& X) {
    double cp_molar = cp(T, X);  // [J/(mol·K)]
    double mw = mwmix(X);         // [g/mol]
    return cp_molar / mw * 1000.0;  // [J/(kg·K)]
}

// Heat capacity at constant volume [J/(kg·K)]
double cv_mass(double T, const std::vector<double>& X) {
    double cv_molar = cv(T, X);  // [J/(mol·K)]
    double mw = mwmix(X);         // [g/mol]
    return cv_molar / mw * 1000.0;  // [J/(kg·K)]
}

// Enthalpy [J/kg]
double h_mass(double T, const std::vector<double>& X) {
    double h_molar = h(T, X);    // [J/mol]
    double mw = mwmix(X);         // [g/mol]
    return h_molar / mw * 1000.0;  // [J/kg]
}

// Entropy [J/(kg·K)]
double s_mass(double T, const std::vector<double>& X, double P, double P_ref) {
    double s_molar = s(T, X, P, P_ref);  // [J/(mol·K)]
    double mw = mwmix(X);                 // [g/mol]
    return s_molar / mw * 1000.0;          // [J/(kg·K)]
}

// Internal energy [J/kg]
double u_mass(double T, const std::vector<double>& X) {
    double u_molar = u(T, X);    // [J/mol]
    double mw = mwmix(X);         // [g/mol]
    return u_molar / mw * 1000.0;  // [J/kg]
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
            const NASA_Coeffs& nasa = nasa_coeffs[i];
            double T_local = T;
            const auto* coeffs = find_nasa9_interval(nasa, T_local, i);
            double ds_dT_i = ds_R_dT_nasa9(*coeffs, T_local);
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
        double T_local = T;
        const auto* coeffs = find_nasa9_interval(nasa, T_local, i);
        double dcp_dT_i = dcp_R_dT_nasa9(*coeffs, T_local);
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

// Ideal gas molar volume [m³/mol]
double molar_volume(double T, double P) {
    return R_GAS * T / P;
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
            // d(H/RT)/dT = -H/(R*T^2) + Cp/(R*T)
            double h_RT_val = h_RT(i, T);
            double cp_R_val = cp_R(i, T);
            double dh_RT_dT = -h_RT_val / T + cp_R_val / T;

            // d(S/R)/dT from NASA-9 polynomial
            const NASA_Coeffs& nasa = nasa_coeffs[i];
            double T_local = T;
            const auto* coeffs = find_nasa9_interval(nasa, T_local, i);
            double ds_R_dT = ds_R_dT_nasa9(*coeffs, T_local);

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

    double T_min = 200.0;
    double T_max = 6000.0;

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

    double T_min = 200.0;
    double T_max = 6000.0;

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

    double T_min = 200.0;
    double T_max = 6000.0;

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

// -------------------------------------------------------------
// State-based overloads
// -------------------------------------------------------------

double mwmix(const State& s) {
    return mwmix(s.X);
}

double cp(const State& s) {
    return cp(s.T, s.X);
}

double h(const State& s) {
    return h(s.T, s.X);
}

double s(const State& st) {
    return s(st.T, st.X, st.P);
}

double cv(const State& s) {
    return cv(s.T, s.X);
}

double u(const State& s) {
    return u(s.T, s.X);
}

double density(const State& s) {
    return density(s.T, s.P, s.X);
}

double specific_gas_constant(const State& s) {
    return specific_gas_constant(s.X);
}

double isentropic_expansion_coefficient(const State& s) {
    return isentropic_expansion_coefficient(s.T, s.X);
}

double speed_of_sound(const State& s) {
    return speed_of_sound(s.T, s.X);
}

// -------------------------------------------------------------
// Air Properties Bundle
// -------------------------------------------------------------

AirProperties air_properties(double T, double P, double humidity) {
    if (T <= 0) {
        throw std::invalid_argument("air_properties: temperature must be positive");
    }
    if (P <= 0) {
        throw std::invalid_argument("air_properties: pressure must be positive");
    }
    if (humidity < 0.0 || humidity > 1.0) {
        throw std::invalid_argument("air_properties: humidity must be in range [0, 1]");
    }

    // Get air composition (dry or humid)
    std::vector<double> X;
    if (humidity > 0.0) {
        X = humid_air_composition(T, P, humidity);
    } else {
        X = standard_dry_air_composition();
    }

    // Compute all properties
    AirProperties props;
    props.T = T;
    props.P = P;

    // Thermodynamic properties
    props.rho = density(T, P, X);
    props.cp = cp_mass(T, X);
    props.cv = cv_mass(T, X);
    props.gamma = props.cp / props.cv;
    props.a = speed_of_sound(T, X);

    // Transport properties
    props.mu = viscosity(T, P, X);
    props.k = thermal_conductivity(T, P, X);

    // Derived properties
    props.nu = props.mu / props.rho;           // Kinematic viscosity
    props.alpha = props.k / (props.rho * props.cp);  // Thermal diffusivity
    props.Pr = prandtl(T, P, X);

    return props;
}

// -------------------------------------------------------------
// Thermodynamic State Bundle
// -------------------------------------------------------------

ThermoState thermo_state(double T, double P, const std::vector<double>& X, double P_ref) {
    if (T <= 0) {
        throw std::invalid_argument("thermo_state: temperature must be positive");
    }
    if (P <= 0) {
        throw std::invalid_argument("thermo_state: pressure must be positive");
    }
    if (P_ref <= 0) {
        throw std::invalid_argument("thermo_state: reference pressure must be positive");
    }
    if (X.empty()) {
        throw std::invalid_argument("thermo_state: mole fractions vector cannot be empty");
    }

    ThermoState state;

    // Echo back inputs
    state.T = T;
    state.P = P;

    // Molecular weight
    state.mw = mwmix(X);

    // Molar-basis thermodynamic properties
    state.cp = cp(T, X);
    state.cv = cv(T, X);
    state.h = h(T, X);
    state.s = s(T, X, P, P_ref);
    state.u = u(T, X);
    state.gamma = isentropic_expansion_coefficient(T, X);
    state.a = speed_of_sound(T, X);

    // Mass-basis thermodynamic properties
    state.cp_mass = cp_mass(T, X);
    state.cv_mass = cv_mass(T, X);
    state.h_mass = h_mass(T, X);
    state.s_mass = s_mass(T, X, P, P_ref);
    state.u_mass = u_mass(T, X);

    // Density
    state.rho = density(T, P, X);

    return state;
}

// -------------------------------------------------------------
// Complete State Bundle (Thermo + Transport)
// -------------------------------------------------------------

CompleteState complete_state(double T, double P, const std::vector<double>& X, double P_ref) {
    CompleteState state;

    // Compute thermodynamic state
    state.thermo = thermo_state(T, P, X, P_ref);

    // Compute transport state
    state.transport = transport_state(T, P, X);

    return state;
}
