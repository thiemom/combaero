#include "../include/humidair.h"
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

// Constants for dry air composition (mole fractions)
const std::unordered_map<std::string, double> dry_air_composition = {
    {"N2", 0.7808},
    {"O2", 0.2095},
    {"AR", 0.0093},  // Match case from generated header
    {"CO2", 0.0004}
};

// Saturation vapor pressure constants
// Reference: Huang (2018), doi:10.1175/JAMC-D-17-0334.1
// Base equation: Ps = exp(a - b/(t + d1)) / (t + d2)^c, t in °C, Ps in Pa

// For water vapor (t > 0°C) - Equation (17)
static constexpr double WATER_A  = 34.494;
static constexpr double WATER_B  = 4924.99;
static constexpr double WATER_D1 = 237.1;
static constexpr double WATER_D2 = 105.0;
static constexpr double WATER_C  = 1.57;

// For ice (t ≤ 0°C) - Equation (18)
static constexpr double ICE_A  = 43.494;
static constexpr double ICE_B  = 6545.8;
static constexpr double ICE_D1 = 278.0;
static constexpr double ICE_D2 = 868.0;
static constexpr double ICE_C  = 2.0;

// Validated range [K]
static constexpr double SVP_T_MIN = 173.15;  // -100°C
static constexpr double SVP_T_MAX = 373.15;  // +100°C

// Maximum relative extrapolation beyond validated range before throwing
// 5% of the range span (200 K) = 10 K
static constexpr double SVP_EXTRAP_TOL = 0.05;

// Water vapor saturation pressure [Pa] using Huang (2018)
// Validated: -100°C to +100°C (173.15–373.15 K)
// Extrapolation up to 5% of range span (~10 K) is allowed silently.
// Beyond that, throws std::out_of_range.
double saturation_vapor_pressure(double T) {
    const double span = SVP_T_MAX - SVP_T_MIN;
    const double tol  = SVP_EXTRAP_TOL * span;
    if (T < SVP_T_MIN - tol || T > SVP_T_MAX + tol) {
        throw std::out_of_range(
            "saturation_vapor_pressure: T=" + std::to_string(T) +
            " K is outside extrapolation range [" +
            std::to_string(SVP_T_MIN - tol) + ", " +
            std::to_string(SVP_T_MAX + tol) + "] K");
    }

    // Convert temperature from K to °C
    double t = T - 273.15;

    // Calculate saturation vapor pressure
    double P_ws;

    if (t <= 0.0) {
        P_ws = std::exp(ICE_A - ICE_B / (t + ICE_D1)) / std::pow(t + ICE_D2, ICE_C);
    } else {
        P_ws = std::exp(WATER_A - WATER_B / (t + WATER_D1)) / std::pow(t + WATER_D2, WATER_C);
    }

    return P_ws;  // Return P_ws in Pa
}

// Calculate actual vapor pressure from relative humidity [Pa]
double vapor_pressure(double T, double RH) {
    // Validate inputs
    if (RH < 0.0 || RH > 1.0) {
        throw std::runtime_error("Relative humidity must be between 0 and 1");
    }

    return RH * saturation_vapor_pressure(T);
}

// Calculate humidity ratio (kg water vapor per kg dry air)
double humidity_ratio(double T, double P, double RH) {
    // Calculate vapor pressure
    double P_w = vapor_pressure(T, RH);

    // Check that vapor pressure is not greater than total pressure
    if (P_w >= P) {
        throw std::runtime_error("Vapor pressure cannot exceed total pressure");
    }

    // Calculate humidity ratio
    // W = 0.621945 * P_w / (P - P_w)
    // 0.621945 = molecular weight of water / molecular weight of dry air
    return 0.621945 * P_w / (P - P_w);
}

// Calculate mole fraction of water vapor in humid air
double water_vapor_mole_fraction(double T, double P, double RH) {
    // Calculate vapor pressure
    double P_w = vapor_pressure(T, RH);

    // Check that vapor pressure is not greater than total pressure
    if (P_w >= P) {
        throw std::runtime_error("Vapor pressure cannot exceed total pressure");
    }

    // Calculate mole fraction of water vapor (P_w / P)
    return P_w / P;
}

// Get standard dry air composition as a vector in the order defined by species_index in thermo_transport_data.h
std::vector<double> standard_dry_air_composition() {
    // Initialize mole fractions vector with zeros
    std::vector<double> X(species_names.size(), 0.0);

    // Set dry air component mole fractions
    for (const auto& component : dry_air_composition) {
        auto it = species_index.find(component.first);
        if (it != species_index.end()) {
            X[it->second] = component.second;
        }
    }

    return X;
}

// Calculate humid air composition (mole fractions)
// Returns a vector of mole fractions in the order defined by species_index in thermo_transport_data.h
std::vector<double> humid_air_composition(double T, double P, double RH) {
    // Get water vapor mole fraction
    double x_w = water_vapor_mole_fraction(T, P, RH);

    // Initialize composition vector with zeros (size = number of species)
    std::vector<double> X(species_names.size(), 0.0);

    // Set water vapor mole fraction if H2O is in the species list
    auto h2o_it = species_index.find("H2O");
    if (h2o_it != species_index.end()) {
        X[h2o_it->second] = x_w;
    }

    // Adjust dry air component mole fractions
    double dry_air_factor = 1.0 - x_w;

    for (const auto& component : dry_air_composition) {
        auto it = species_index.find(component.first);
        if (it != species_index.end()) {
            X[it->second] = component.second * dry_air_factor;
        }
    }

    return X;
}

// Calculate dewpoint temperature from T, P, RH [K]
// P is not used (dewpoint depends negligibly on pressure at atmospheric conditions)
double dewpoint(double T, double P, double RH) {
    (void)P;

    if (RH <= 0.0) {
        throw std::invalid_argument("dewpoint: RH must be > 0");
    }

    // Calculate vapor pressure
    double P_w = vapor_pressure(T, RH);

    // Use iterative method (Newton-Raphson) to solve for Tdp

    // Initial guess for dewpoint (slightly lower than ambient temperature)
    double Tdp = T - 10.0;

    // Ensure initial guess is within valid range
    Tdp = std::max(Tdp, 173.15); // Minimum -100°C
    Tdp = std::min(Tdp, T);      // Maximum ambient temperature

    // Iterative solution
    const std::size_t max_iterations = 100;
    const double tolerance = 1e-6;

    for (std::size_t i = 0; i < max_iterations; i++) {
        // Calculate saturation vapor pressure at current Tdp estimate
        double P_ws_Tdp = saturation_vapor_pressure(Tdp);

        // Calculate approximate derivative of saturation vapor pressure with respect to T
        // Use central difference approximation
        const double delta = 0.01; // Small temperature difference for derivative calculation
        double P_ws_Tdp_plus = saturation_vapor_pressure(Tdp + delta);
        double P_ws_Tdp_minus = saturation_vapor_pressure(Tdp - delta);
        double dP_ws_dT = (P_ws_Tdp_plus - P_ws_Tdp_minus) / (2.0 * delta);

        // Update Tdp using Newton-Raphson method
        double Tdp_new = Tdp - (P_ws_Tdp - P_w) / dP_ws_dT;

        // Ensure Tdp stays within valid range
        Tdp_new = std::max(Tdp_new, 173.15); // Minimum -100°C
        Tdp_new = std::min(Tdp_new, T);      // Maximum ambient temperature

        // Check for convergence
        if (std::abs(Tdp_new - Tdp) < tolerance) {
            return Tdp_new;
        }

        Tdp = Tdp_new;
    }

    throw std::runtime_error("Dewpoint calculation did not converge");
}

// Calculate relative humidity from dewpoint temperature [fraction]
// P is not used (kept for API symmetry with dewpoint())
double relative_humidity_from_dewpoint(double T, double Tdp, double P) {
    (void)P;

    if (Tdp > T) {
        throw std::invalid_argument("relative_humidity_from_dewpoint: Tdp cannot exceed T");
    }

    // Calculate saturation vapor pressure at dewpoint temperature
    double P_w = saturation_vapor_pressure(Tdp);

    // Calculate saturation vapor pressure at ambient temperature
    double P_ws = saturation_vapor_pressure(T);

    // Calculate relative humidity
    return P_w / P_ws;
}

// Derivative of saturation_vapor_pressure with respect to T [Pa/K]
// Used internally by wet_bulb_temperature Newton-Raphson solver
static double d_saturation_vapor_pressure_dT(double T) {
    const double delta = 0.01;
    return (saturation_vapor_pressure(T + delta) - saturation_vapor_pressure(T - delta)) / (2.0 * delta);
}

// Calculate wet-bulb temperature [K]
// Solves the psychrometric equation using Newton-Raphson:
//   f(Twb) = cp_air*(T - Twb) - hfg(Twb)*(Ws(Twb) - W) = 0
// where hfg(Twb) = 2501000 - 2381*(Twb - 273.15) [J/kg] (linear fit, ASHRAE)
double wet_bulb_temperature(double T, double P, double RH) {
    const double W = humidity_ratio(T, P, RH);
    static constexpr double cp_air = 1005.0;  // J/(kg·K)

    // hfg as a function of Twb: linear approximation valid 0-100°C
    auto hfg = [](double Tw) { return 2501000.0 - 2381.0 * (Tw - 273.15); };

    // Ws(Twb): saturation humidity ratio at wet-bulb temperature
    auto Ws = [&](double Tw) {
        double P_ws = saturation_vapor_pressure(Tw);
        return 0.621945 * P_ws / (P - P_ws);
    };

    // f(Twb) = cp_air*(T - Twb) - hfg(Twb)*(Ws(Twb) - W)
    auto f = [&](double Tw) {
        return cp_air * (T - Tw) - hfg(Tw) * (Ws(Tw) - W);
    };

    // df/dTwb = -cp_air - dhfg/dTwb*(Ws - W) - hfg*dWs/dTwb
    auto df = [&](double Tw) {
        const double dWs_dTw = 0.621945 * d_saturation_vapor_pressure_dT(Tw)
                               / std::pow(P - saturation_vapor_pressure(Tw), 2.0) * P;
        const double dhfg_dTw = -2381.0;
        return -cp_air - dhfg_dTw * (Ws(Tw) - W) - hfg(Tw) * dWs_dTw;
    };

    double Twb = T * (1.0 - 0.5 * RH);  // initial guess: between Tdp and T
    Twb = std::max(Twb, SVP_T_MIN);
    Twb = std::min(Twb, T);

    const std::size_t max_iter = 50;
    const double tol = 1e-6;

    for (std::size_t i = 0; i < max_iter; ++i) {
        double fval = f(Twb);
        double dfval = df(Twb);
        if (std::abs(dfval) < 1e-15) break;
        double step = fval / dfval;
        Twb -= step;
        Twb = std::max(Twb, SVP_T_MIN);
        Twb = std::min(Twb, T);
        if (std::abs(step) < tol) return Twb;
    }

    throw std::runtime_error("wet_bulb_temperature: Newton-Raphson did not converge");
}

// Calculate enthalpy of humid air [J/kg]
double humid_air_enthalpy(double T, double P, double RH) {
    // Calculate humidity ratio
    double W = humidity_ratio(T, P, RH);

    // Calculate enthalpy of dry air
    double h_da = 1005.0 * (T - 273.15);  // cp_air * (T - T_ref), cp_air = 1005 J/(kg·K)

    // Calculate enthalpy of water vapor
    // h_wv = h_fg + cp_wv * (T - T_ref)
    // h_fg = 2501000 J/kg at 0°C, cp_wv = 1860 J/(kg·K)
    double h_wv = 2501000.0 + 1860.0 * (T - 273.15);

    // Calculate total enthalpy of humid air
    return h_da + W * h_wv;  // J/kg of dry air
}

// Calculate enthalpy of humid air using NASA-9 thermodynamic data [J/kg dry air]
// Uses h(T, X) from the NASA-9 polynomial database, converted to per-kg-dry-air basis.
double humid_air_enthalpy_nasa9(double T, double P, double RH) {
    const std::vector<double> X = humid_air_composition(T, P, RH);
    const double W = humidity_ratio(T, P, RH);
    // h(T, X) returns J/mol of mixture; convert to J/kg dry air
    // MW_mix [g/mol], 1 kg dry air = (1+W) kg humid air
    // h_per_kg_mix = h(T,X) * 1000 / mwmix(X)  [J/kg humid air]
    // h_per_kg_dry = h_per_kg_mix * (1 + W)     [J/kg dry air]
    const double h_mix_J_per_mol = h(T, X);
    const double MW_mix = mwmix(X);  // g/mol
    return h_mix_J_per_mol * 1000.0 / MW_mix * (1.0 + W);
}

// Calculate density of humid air [kg/m³]
double humid_air_density(double T, double P, double RH) {
    // Calculate humidity ratio
    double W = humidity_ratio(T, P, RH);

    // Calculate gas constant for humid air
    // R_ha = R_da * (1 + W/0.621945) / (1 + W)
    // R_da = 287.058 J/(kg·K)
    double R_ha = 287.058 * (1.0 + W/0.621945) / (1.0 + W);

    // Calculate density using ideal gas law
    return P / (R_ha * T);
}
