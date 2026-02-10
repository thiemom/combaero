#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"
#include "../include/transport.h"
#include "../include/humidair.h"
#include <iostream>
#include <vector>
#include <iomanip>

// CombAero Library Main Example
//
// This file provides a simple demonstration of the combaero library's capabilities.
// For more detailed examples, see the examples directory:
// - examples/thermo_example.cpp: Demonstrates thermodynamic and transport properties
// - examples/humidair_example.cpp: Demonstrates humid air calculations

int main() {
    // Set precision for output
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "===========================================" << std::endl;
    std::cout << "        CombAero Library Demo" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "This is a simple demonstration of the combaero library." << std::endl;
    std::cout << "For more detailed examples, see the examples directory." << std::endl;

    // Define temperature and pressure
    double T = 300.0;    // K
    double P = 101325.0; // Pa (1 atm)

    // Create a vector for mole fractions (initialize with zeros)
    std::vector<double> X(species_names.size(), 0.0);

    // Set mole fractions for air (O2 = 0.21, N2 = 0.79)
    std::size_t o2_idx = species_index_from_name("O2");
    std::size_t n2_idx = species_index_from_name("N2");

    X[o2_idx] = 0.21;
    X[n2_idx] = 0.79;

    // Print basic properties of air
    std::cout << "\n1. Basic Air Properties (T = 300K, P = 1 atm):" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Specific Heat (Cp): " << cp(T, X) << " J/(mol·K)" << std::endl;
    std::cout << "Enthalpy (H): " << h(T, X) << " J/mol" << std::endl;
    std::cout << "Entropy (S): " << s(T, X, P) << " J/(mol·K)" << std::endl;
    std::cout << "Density: " << density(T, P, X) << " kg/m³" << std::endl;
    std::cout << "Viscosity: " << viscosity(T, P, X) * 1.0e6 << " μPa·s" << std::endl;
    std::cout << "Thermal Conductivity: " << thermal_conductivity(T, P, X) << " W/(m·K)" << std::endl;
    std::cout << "Prandtl Number: " << prandtl(T, P, X) << std::endl;

    // Demonstrate humid air calculations
    std::cout << "\n2. Humid Air Properties (T = 25°C, P = 1 atm):" << std::endl;
    std::cout << "------------------------------------------" << std::endl;

    double T_humid = 298.15; // 25°C
    double RH = 0.5;         // 50% relative humidity

    // Calculate humid air properties
    double vp = vapor_pressure(T_humid, RH);
    double hr = humidity_ratio(T_humid, P, RH);
    double dp = dewpoint(T_humid, P, RH);
    double rho = humid_air_density(T_humid, P, RH);

    std::cout << "Temperature: " << T_humid - 273.15 << " °C" << std::endl;
    std::cout << "Relative Humidity: " << RH * 100.0 << "%" << std::endl;
    std::cout << "Vapor Pressure: " << vp << " Pa" << std::endl;
    std::cout << "Humidity Ratio: " << hr << " kg H₂O/kg dry air" << std::endl;
    std::cout << "Dewpoint: " << dp - 273.15 << " °C" << std::endl;
    std::cout << "Density: " << rho << " kg/m³" << std::endl;

    // Get humid air composition
    std::vector<double> humid_comp = humid_air_composition(T_humid, P, RH);

    std::cout << "\nHumid Air Composition (Mole Fractions):" << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    for (size_t i = 0; i < species_names.size(); i++) {
        if (humid_comp[i] > 0.001) { // Only show significant components
            std::cout << species_names[i] << ": " << humid_comp[i] << std::endl;
        }
    }

    std::cout << "\nFor more examples, run:" << std::endl;
    std::cout << "  ./combaero_example - For thermodynamic and transport properties" << std::endl;
    std::cout << "  ./humidair_example - For humid air calculations" << std::endl;

    return 0;
}
