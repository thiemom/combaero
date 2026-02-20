#include "../include/combustion.h"
#include "../include/humidair.h"
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"
#include "../include/transport.h"
#include "../include/utils.h"
#include <iomanip>
#include <iostream>
#include <vector>

// Thermodynamic and Transport Properties Example
//
// This example demonstrates the calculation of various thermodynamic
// and transport properties for gas mixtures using the thermo library.

int main() {
    // Set precision for output
    std::cout << std::fixed << std::setprecision(6);

    // Define temperature and pressure
    double T = 300.0;    // K
    double P = 101325.0; // Pa (1 atm)

    // Get species indices
    const std::size_t ch4_idx = species_index_from_name("CH4");
    const std::size_t h2_idx  = species_index_from_name("H2");

    // Standard dry air composition (N2, O2, Ar, CO2)
    const std::vector<double> X_air = standard_dry_air_composition();

    std::cout << "=========================================" << std::endl;
    std::cout << "Thermodynamic and Transport Properties" << std::endl;
    std::cout << "=========================================" << std::endl;

    // Example 1: Air mixture
    std::cout << "\n1. Standard Air Mixture" << std::endl;
    std::cout << "----------------------" << std::endl;

    print_mixture_properties(T, P, X_air);

    // Example 2: Air-fuel mixture (10% CH4 blended into air)
    std::cout << "\n\n2. Air-Fuel Mixture" << std::endl;
    std::cout << "-------------------" << std::endl;

    // Blend 10% CH4 into air (renormalize)
    std::vector<double> X_af = X_air;
    for (auto& x : X_af) x *= 0.90;
    X_af[ch4_idx] += 0.10;

    print_mixture_properties(T, P, X_af);

    // Example 3: Temperature dependence of properties (air)
    std::cout << "\n\n3. Temperature Dependence of Properties" << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Temp (K)   Cp (J/mol·K)   Gamma   Sound Speed (m/s)   Viscosity (μPa·s)   Thermal Cond (W/m·K)" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;

    for (double temp = 300.0; temp <= 2600.0; temp += 300.0) {
        std::cout << std::setw(8) << temp << "   "
                  << std::setw(12) << cp(temp, X_air) << "   "
                  << std::setw(7) << isentropic_expansion_coefficient(temp, X_air) << "   "
                  << std::setw(16) << speed_of_sound(temp, X_air) << "   "
                  << std::setw(16) << viscosity(temp, P, X_air) * 1.0e6 << "   "
                  << std::setw(20) << thermal_conductivity(temp, P, X_air) << std::endl;
    }

    // Example 4: Fuel properties
    std::cout << "\n4. Fuel Properties" << std::endl;
    std::cout << "----------------" << std::endl;

    // For each fuel in the system
    for (const auto& fuel_name : {"CH4", "H2"}) {
        const std::size_t fuel_idx = species_index_from_name(fuel_name);
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
    std::vector<double> X_fuel(species_names.size(), 0.0);
    X_fuel[ch4_idx] = 0.7;  // 70% methane
    X_fuel[h2_idx]  = 0.3;  // 30% hydrogen

    std::cout << "Mixture: 70% CH4, 30% H2" << std::endl;
    std::cout << "  Oxygen required: " << oxygen_required_per_mol_mixture(X_fuel)
              << " mol O2/mol mixture" << std::endl;
    std::cout << "  Oxygen required: " << oxygen_required_per_kg_mixture(X_fuel)
              << " kg O2/kg mixture" << std::endl;

    // Example 6: Inverse calculations (using air)
    std::cout << "\n6. Inverse Calculations (Finding Temperature)" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    double T_ref = 1200.0;

    double h_ref  = h(T_ref, X_air);
    double s_ref  = s(T_ref, X_air, P);
    double cp_ref = cp(T_ref, X_air);

    double T_from_h  = calc_T_from_h(h_ref, X_air);
    double T_from_s  = calc_T_from_s(s_ref, P, X_air);
    double T_from_cp = calc_T_from_cp(cp_ref, X_air);

    std::cout << "Reference Temperature: " << T_ref << " K\n";
    std::cout << "Temperature calculated from enthalpy: " << T_from_h << " K\n";
    std::cout << "Temperature calculated from entropy: " << T_from_s << " K\n";
    std::cout << "Temperature calculated from heat capacity: " << T_from_cp << " K\n";
    std::cout << "Relative error (enthalpy method): " << std::abs(T_from_h - T_ref) / T_ref * 100.0 << "%\n";
    std::cout << "Relative error (entropy method): " << std::abs(T_from_s - T_ref) / T_ref * 100.0 << "%\n";
    std::cout << "Relative error (heat capacity method): " << std::abs(T_from_cp - T_ref) / T_ref * 100.0 << "%\n";

    return 0;
}
