#include "../include/thermo_transport.h"
#include <iostream>
#include <vector>
#include <iomanip>

/**
 * Thermodynamic and Transport Properties Example
 * 
 * This example demonstrates the calculation of various thermodynamic
 * and transport properties for gas mixtures using the thermo library.
 */
int main() {
    // Set precision for output
    std::cout << std::fixed << std::setprecision(6);
    
    // Define temperature and pressure
    double T = 300.0;    // K
    double P = 101325.0; // Pa (1 atm)
    
    // Create a vector for mole fractions (initialize with zeros)
    std::vector<double> X(species_names.size(), 0.0);
    
    // Get species indices
    int o2_idx = species_index_from_name("O2");
    int n2_idx = species_index_from_name("N2");
    int ch4_idx = species_index_from_name("CH4");
    int h2_idx = species_index_from_name("H2");
    
    std::cout << "=========================================" << std::endl;
    std::cout << "Thermodynamic and Transport Properties" << std::endl;
    std::cout << "=========================================" << std::endl;
    
    // Example 1: Air mixture (O2 = 0.21, N2 = 0.79)
    std::cout << "\n1. Standard Air Mixture" << std::endl;
    std::cout << "----------------------" << std::endl;
    
    X[o2_idx] = 0.21;
    X[n2_idx] = 0.79;
    
    // Print all properties of the air mixture
    print_mixture_properties(T, P, X);
    
    // Example 2: Air-fuel mixture
    std::cout << "\n\n2. Air-Fuel Mixture" << std::endl;
    std::cout << "-------------------" << std::endl;
    
    // Reset mixture
    std::fill(X.begin(), X.end(), 0.0);
    
    // Create a mixture with methane and air
    X[ch4_idx] = 0.10;  // 10% methane
    X[o2_idx] = 0.19;   // 19% oxygen
    X[n2_idx] = 0.71;   // 71% nitrogen
    
    // Print all properties of the mixture
    print_mixture_properties(T, P, X);
    
    // Example 3: Temperature dependence of properties
    std::cout << "\n\n3. Temperature Dependence of Properties" << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Temp (K)   Cp (J/mol·K)   Gamma   Sound Speed (m/s)   Viscosity (μPa·s)   Thermal Cond (W/m·K)" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
    
    // Use temperature range valid for all species
    for (double temp = 300.0; temp <= 2600.0; temp += 300.0) {
        std::cout << std::setw(8) << temp << "   " 
                  << std::setw(12) << cp(temp, X) << "   "
                  << std::setw(7) << isentropic_expansion_coefficient(temp, X) << "   "
                  << std::setw(16) << speed_of_sound(temp, X) << "   "
                  << std::setw(16) << viscosity(temp, P, X) * 1.0e6 << "   "
                  << std::setw(20) << thermal_conductivity(temp, P, X) << std::endl;
    }
    
    // Example 4: Fuel properties
    std::cout << "\n4. Fuel Properties" << std::endl;
    std::cout << "----------------" << std::endl;
    
    // For each fuel in the system
    for (const auto& fuel_name : {"CH4", "H2"}) {
        int fuel_idx = species_index_from_name(fuel_name);
        std::cout << "Fuel: " << fuel_name << std::endl;
        std::cout << "  Molecular structure: C=" << molecular_structures[fuel_idx].C 
                  << ", H=" << molecular_structures[fuel_idx].H 
                  << ", O=" << molecular_structures[fuel_idx].O 
                  << ", N=" << molecular_structures[fuel_idx].N << std::endl;
        std::cout << "  Oxygen required: " << oxygen_required_per_mol_fuel(fuel_idx) 
                  << " mol O2/mol fuel" << std::endl;
        std::cout << "  Oxygen required: " << oxygen_required_per_kg_fuel(fuel_idx) 
                  << " kg O2/kg fuel" << std::endl;
        std::cout << std::endl;
    }
    
    // Example 5: Fuel mixture properties
    std::cout << "\n5. Fuel Mixture Properties" << std::endl;
    std::cout << "------------------------" << std::endl;
    
    // Create a mixture with both fuels
    std::fill(X.begin(), X.end(), 0.0);
    X[ch4_idx] = 0.7;  // 70% methane
    X[h2_idx] = 0.3;   // 30% hydrogen
    
    std::cout << "Mixture: 70% CH4, 30% H2" << std::endl;
    std::cout << "  Oxygen required: " << oxygen_required_per_mol_mixture(X) 
              << " mol O2/mol mixture" << std::endl;
    std::cout << "  Oxygen required: " << oxygen_required_per_kg_mixture(X) 
              << " kg O2/kg mixture" << std::endl;
    
    // Example 6: Inverse calculations
    std::cout << "\n6. Inverse Calculations (Finding Temperature)" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    
    // Select a reference temperature to test against
    double T_ref = 1200.0;
    
    // Calculate properties at the reference temperature
    double h_ref = h(T_ref, X);
    double s_ref = s(T_ref, P, X);
    double cp_ref = cp(T_ref, X);
    
    // Solve for temperature given enthalpy
    double T_from_h = calc_T_from_h(h_ref, X);
    
    // Solve for temperature given entropy
    double T_from_s = calc_T_from_s(s_ref, P, X);
    
    // Solve for temperature given heat capacity
    double T_from_cp = calc_T_from_cp(cp_ref, X);
    
    std::cout << "Reference Temperature: " << T_ref << " K\n";
    std::cout << "Temperature calculated from enthalpy: " << T_from_h << " K\n";
    std::cout << "Temperature calculated from entropy: " << T_from_s << " K\n";
    std::cout << "Temperature calculated from heat capacity: " << T_from_cp << " K\n";
    std::cout << "Relative error (enthalpy method): " << std::abs(T_from_h - T_ref) / T_ref * 100.0 << "%\n";
    std::cout << "Relative error (entropy method): " << std::abs(T_from_s - T_ref) / T_ref * 100.0 << "%\n";
    std::cout << "Relative error (heat capacity method): " << std::abs(T_from_cp - T_ref) / T_ref * 100.0 << "%\n";
    
    return 0;
}
