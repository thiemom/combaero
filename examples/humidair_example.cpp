#include <iostream>
#include <iomanip>
#include <vector>
#include "../include/thermo_transport.h"
#include "../include/humidair.h"

/**
 * Humid Air Properties Example
 * 
 * This example demonstrates the calculation of various properties
 * of humid air using the Hyland and Wexler equations.
 */
int main() {
    // Set precision for output
    std::cout << std::fixed << std::setprecision(6);
    
    // Temperature, pressure, and relative humidity conditions to test
    std::vector<double> temperatures = {273.15, 293.15, 303.15, 313.15, 323.15}; // K (0°C, 20°C, 30°C, 40°C, 50°C)
    double pressure = 101325.0; // Pa (standard atmospheric pressure)
    std::vector<double> rel_humidities = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0}; // Relative humidity (0% to 100%)
    
    // Print header
    std::cout << "==========================================" << std::endl;
    std::cout << "Humid Air Properties Calculator" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Example 1: Dry air composition
    std::cout << "\n1. Dry Air Composition" << std::endl;
    std::cout << "---------------------" << std::endl;
    
    std::vector<double> dry_air = standard_dry_air_composition();
    
    std::cout << "Component | Mole Fraction" << std::endl;
    std::cout << "----------------------" << std::endl;
    
    for (size_t i = 0; i < species_names.size(); i++) {
        if (dry_air[i] > 0.0) {
            std::cout << std::setw(10) << species_names[i] << " | " 
                      << std::setw(12) << dry_air[i] << std::endl;
        }
    }
    
    // Example 2: Saturation vapor pressure
    std::cout << "\n2. Saturation Vapor Pressure (Hyland-Wexler)" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Temperature (°C) | Saturation Pressure (Pa)" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    
    for (double T : temperatures) {
        std::cout << std::setw(15) << T - 273.15 << " | " 
                  << std::setw(22) << saturation_vapor_pressure(T) << std::endl;
    }
    
    // Example 3: Humid air properties at different conditions
    std::cout << "\n3. Humid Air Properties at P = " << pressure << " Pa" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    
    for (double T : temperatures) {
        std::cout << "\nTemperature: " << T - 273.15 << " °C" << std::endl;
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << "RH | Vapor Press (Pa) | Humidity Ratio | Dewpoint (°C) | Density (kg/m³)" << std::endl;
        std::cout << "-------------------------------------------" << std::endl;
        
        for (double RH : rel_humidities) {
            double vp = vapor_pressure(T, RH);
            double hr = humidity_ratio(T, pressure, RH);
            double rho = humid_air_density(T, pressure, RH);
            
            std::cout << std::setw(3) << RH * 100 << "% | " 
                      << std::setw(15) << vp << " | "
                      << std::setw(14) << hr << " | ";
            
            // Handle dewpoint calculation (only valid for RH > 0)
            if (RH > 0.0) {
                try {
                    double dp = dewpoint(T, pressure, RH) - 273.15; // Convert to °C
                    std::cout << std::setw(13) << dp << " | ";
                } catch (const std::exception& e) {
                    std::cout << std::setw(13) << "N/A" << " | ";
                }
            } else {
                std::cout << std::setw(13) << "N/A" << " | "; // Dewpoint undefined at RH=0
            }
            
            std::cout << std::setw(14) << rho << std::endl;
        }
    }
    
    // Example 4: Humid air composition
    std::cout << "\n4. Humid Air Composition Example" << std::endl;
    std::cout << "----------------------------" << std::endl;
    double example_T = 303.15; // 30°C
    double example_RH = 0.6;   // 60% RH
    
    std::cout << "Conditions: T = " << example_T - 273.15 << " °C, P = " << pressure << " Pa, RH = " << example_RH * 100 << "%" << std::endl;
    
    // Calculate humid air composition
    std::vector<double> X = humid_air_composition(example_T, pressure, example_RH);
    
    std::cout << "\nMole Fractions:" << std::endl;
    std::cout << "Component | Mole Fraction" << std::endl;
    std::cout << "----------------------" << std::endl;
    for (size_t i = 0; i < X.size(); i++) {
        if (X[i] > 0.001) { // Only show components with significant fractions
            std::cout << std::setw(10) << species_names[i] << " | " 
                      << std::setw(12) << X[i] << std::endl;
        }
    }
    
    // Example 5: Thermodynamic properties of humid air
    std::cout << "\n5. Thermodynamic Properties of Humid Air" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    
    std::cout << "Temperature: " << example_T - 273.15 << " °C, Pressure: " << pressure << " Pa, RH: " << example_RH * 100 << "%" << std::endl;
    
    // Calculate properties
    double cp_value = cp(example_T, X);
    double h_value = h(example_T, X);
    double s_value = s(example_T, pressure, X);
    double rho_value = density(example_T, pressure, X);
    double mu_value = viscosity(example_T, pressure, X);
    double k_value = thermal_conductivity(example_T, pressure, X);
    double pr_value = prandtl(example_T, pressure, X);
    
    std::cout << "Specific Heat (Cp): " << cp_value << " J/(mol·K)" << std::endl;
    std::cout << "Enthalpy (H): " << h_value << " J/mol" << std::endl;
    std::cout << "Entropy (S): " << s_value << " J/(mol·K)" << std::endl;
    std::cout << "Density: " << rho_value << " kg/m³" << std::endl;
    std::cout << "Viscosity: " << mu_value * 1.0e6 << " μPa·s" << std::endl;
    std::cout << "Thermal Conductivity: " << k_value << " W/(m·K)" << std::endl;
    std::cout << "Prandtl Number: " << pr_value << std::endl;
    
    return 0;
}
