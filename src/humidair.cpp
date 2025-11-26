#include "../include/humidair.h"
#include "../include/thermo_transport_data.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>

// Constants for dry air composition (mole fractions)
const std::unordered_map<std::string, double> dry_air_composition = {
    {"N2", 0.7808},
    {"O2", 0.2095},
    {"AR", 0.0093},  // Match case from generated header
    {"CO2", 0.0004}
};

// Constants for saturation vapor pressure calculation
// Reference: "A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice"
// Author: Jianhua Huang
// DOI: https://doi.org/10.1175/JAMC-D-17-0334.1
// Publication Date: 01 Jun 2018, Pages: 1265–1272

// Base equation (16): Ps = exp(a - b/(t + d1)) / (t + d2)^c
// where t is in °C and Ps is in Pa

// For water vapor (t > 0°C) - Equation (17)
const double WATER_A = 34.494;
const double WATER_B = 4924.99;
const double WATER_D1 = 237.1;  
const double WATER_D2 = 105.0;
const double WATER_C = 1.57;

// For ice (t ≤ 0°C) - Equation (18)
const double ICE_A = 43.494;
const double ICE_B = 6545.8;
const double ICE_D1 = 278.0;
const double ICE_D2 = 868.0;
const double ICE_C = 2.0;

// Water vapor saturation pressure [Pa]
// Validated for -100°C ≤ T ≤ 100°C
// Implementation of validated equations from Huang (2018)
// Note: For temperatures outside the validated range, results are extrapolated
// and may be less accurate, but represent the best available estimate
double saturation_vapor_pressure(double T) {
    // Check temperature range - issue warning for temperatures outside validated range
    if (T < 173.15 || T > 373.15) {
        std::cerr << "Warning: Temperature " << T << " K (" << T - 273.15 << " °C) is outside validated range for saturation vapor pressure calculation (-100°C to 100°C)" << std::endl;
        std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
    }
    
    // Convert temperature from K to °C
    double t = T - 273.15;
    
    // Calculate saturation vapor pressure
    double P_ws;
    
    if (t <= 0.0) {
        // For ice (t ≤ 0°C) - Equation (18) with denominator term
        // Validated against reference values from the paper
        // Gives 611.29 Pa at 0°C (0.023% error from reference 611.153 Pa)
        P_ws = exp(ICE_A - ICE_B / (t + ICE_D1)) / pow(t + ICE_D2, ICE_C);
    } else {
        // For water vapor (t > 0°C) - Equation (17)
        // Validated against reference values from the paper
        // Gives 611.25 Pa at 0°C (0.015% error from reference 611.153 Pa)
        P_ws = exp(WATER_A - WATER_B / (t + WATER_D1)) / pow(t + WATER_D2, WATER_C);
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
    double x_w = water_vapor_mole_fraction(T, RH, P);
    
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
double dewpoint(double T, double P, double RH) {
    // P is unused in this implementation but kept for API consistency
    (void)P; // Suppress unused parameter warning
    
    // Validate inputs
    if (RH <= 0.0) {
        throw std::runtime_error("Relative humidity must be greater than 0 for dewpoint calculation");
    }
    
    if (T < 173.15 || T > 373.15) {
        throw std::runtime_error("Temperature out of valid range for dewpoint calculation (-100°C to 100°C)");
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
double relative_humidity_from_dewpoint(double T, double Tdp, double P) {
    // P is unused in this implementation but kept for API consistency
    (void)P; // Suppress unused parameter warning
    
    // Validate inputs
    if (Tdp > T) {
        throw std::runtime_error("Dewpoint temperature cannot exceed ambient temperature");
    }
    
    // Validate temperature ranges
    if (T < 273.15 || T > 473.15) {
        throw std::runtime_error("Temperature out of valid range for Hyland-Wexler equation (0°C to 200°C)");
    }
    
    if (Tdp < 273.15 || Tdp > 473.15) {
        throw std::runtime_error("Dewpoint temperature out of valid range for Hyland-Wexler equation (0°C to 200°C)");
    }
    
    // Calculate saturation vapor pressure at dewpoint temperature
    double P_w = saturation_vapor_pressure(Tdp);
    
    // Calculate saturation vapor pressure at ambient temperature
    double P_ws = saturation_vapor_pressure(T);
    
    // Calculate relative humidity
    return P_w / P_ws;
}

// Calculate wet-bulb temperature [K]
double wet_bulb_temperature(double T, double P, double RH) {
    // Use iterative method to find wet-bulb temperature
    // Wet-bulb temperature is the temperature at which evaporation of water
    // would saturate air at the same enthalpy
    
    // Initial guess for wet-bulb temperature
    double Twb = T;
    
    // Calculate humidity ratio at current conditions
    double W = humidity_ratio(T, P, RH);
    
    // Iterative solution
    const std::size_t max_iterations = 100;
    const double tolerance = 1e-6;
    
    for (std::size_t i = 0; i < max_iterations; i++) {
        // Calculate saturation vapor pressure at wet-bulb temperature
        double P_ws_Twb = saturation_vapor_pressure(Twb);
        
        // Calculate saturation humidity ratio at wet-bulb temperature
        double Ws = 0.621945 * P_ws_Twb / (P - P_ws_Twb);
        
        // Calculate psychrometric constant (approximation)
        double cp_air = 1005.0;  // J/(kg·K)
        double hfg = 2501000.0;  // J/kg at 0°C
        
        // Update wet-bulb temperature
        double Twb_new = T - (hfg * (W - Ws)) / cp_air;
        
        // Check for convergence
        if (std::abs(Twb_new - Twb) < tolerance) {
            return Twb_new;
        }
        
        Twb = Twb_new;
    }
    
    throw std::runtime_error("Wet-bulb temperature calculation did not converge");
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
