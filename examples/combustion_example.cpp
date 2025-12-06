#include "../include/state.h"
#include "../include/thermo.h"
#include "../include/combustion.h"
#include "../include/equilibrium.h"
#include "../include/humidair.h"
#include "../include/thermo_transport_data.h"
#include <iostream>
#include <iomanip>
#include <vector>

/**
 * Combustion Example: Fuel + Humid Air -> Mix -> Combust
 *
 * Demonstrates:
 * - Using set_fuel_stream_for_phi to configure fuel stream for target phi
 * - Stream mixing with enthalpy balance
 * - Complete combustion (CO2 + H2O products)
 * - WGS equilibrium (CO + H2O <-> CO2 + H2)
 * - Varying equivalence ratio from lean to stoichiometric
 */
int main()
{
    std::cout << std::fixed << std::setprecision(4);

    // =========================================================================
    // Setup: Define fuel and oxidizer streams
    // =========================================================================

    const std::size_t n_species = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");

    // Fuel stream: pure methane at 300 K (mdot will be set by set_fuel_stream_for_phi)
    Stream fuel;
    fuel.state.T = 300.0;
    fuel.state.P = 101325.0;
    fuel.state.X = std::vector<double>(n_species, 0.0);
    fuel.state.X[idx_CH4] = 1.0;

    // Oxidizer stream: humid air at 25°C, 60% RH, 10 kg/s
    Stream air;
    air.state.T = 298.15;
    air.state.P = 101325.0;
    air.state.X = humid_air_composition(air.state.T, air.state.P, 0.60);
    air.mdot = 10.0;  // kg/s (fixed)

    std::cout << "=========================================================================\n";
    std::cout << "Combustion Example: CH4 + Humid Air\n";
    std::cout << "=========================================================================\n\n";

    std::cout << "Fuel: Pure CH4 at " << fuel.state.T << " K\n";
    std::cout << "Air:  Humid air at " << air.state.T << " K, 60% RH, "
              << air.mdot << " kg/s\n\n";

    // =========================================================================
    // Sweep equivalence ratio from 0.45 to 1.0
    // =========================================================================

    std::cout << "Equivalence Ratio Sweep\n";
    std::cout << "-----------------------\n\n";

    std::cout << std::setw(6) << "phi"
              << std::setw(10) << "mdot_f"
              << std::setw(10) << "T_mix"
              << std::setw(10) << "T_ad"
              << std::setw(10) << "T_wgs"
              << std::setw(12) << "mu_ad"
              << std::setw(12) << "Pr_ad"
              << "\n";
    std::cout << std::setw(6) << "[-]"
              << std::setw(10) << "[kg/s]"
              << std::setw(10) << "[K]"
              << std::setw(10) << "[K]"
              << std::setw(10) << "[K]"
              << std::setw(12) << "[uPa.s]"
              << std::setw(12) << "[-]"
              << "\n";
    std::cout << std::string(70, '-') << "\n";

    for (double phi = 0.45; phi <= 1.01; phi += 0.05) {
        // Use set_fuel_stream_for_phi to get fuel stream with correct mdot
        Stream fuel_phi = set_fuel_stream_for_phi(phi, fuel, air);

        // Mix streams (enthalpy balance gives mixed temperature)
        Stream mixed = mix({fuel_phi, air});

        // Complete combustion (adiabatic, CO2 + H2O products)
        State burned_complete = complete_combustion(mixed.state);

        // WGS equilibrium on burned products (shifts CO + H2O <-> CO2 + H2)
        State burned_wgs = wgs_equilibrium_adiabatic(burned_complete);

        // Output results
        std::cout << std::setw(6) << phi
                  << std::setw(10) << fuel_phi.mdot
                  << std::setw(10) << mixed.state.T
                  << std::setw(10) << burned_complete.T
                  << std::setw(10) << burned_wgs.T
                  << std::setw(12) << burned_complete.mu() * 1.0e6  // convert to μPa·s
                  << std::setw(12) << burned_complete.Pr()
                  << "\n";
    }

    std::cout << "\n";

    // =========================================================================
    // Detailed output for stoichiometric case
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "Detailed Results at Stoichiometric (phi = 1.0)\n";
    std::cout << "=========================================================================\n\n";

    // Use set_fuel_stream_for_phi for stoichiometric mixture
    Stream fuel_stoich = set_fuel_stream_for_phi(1.0, fuel, air);
    Stream mixed_stoich = mix({fuel_stoich, air});
    State burned_stoich = complete_combustion(mixed_stoich.state);
    State wgs_stoich = wgs_equilibrium_adiabatic(burned_stoich);

    std::cout << "Mixed Stream (before combustion):\n";
    std::cout << "  T = " << mixed_stoich.state.T << " K\n";
    std::cout << "  P = " << mixed_stoich.state.P << " Pa\n";
    std::cout << "  mdot = " << mixed_stoich.mdot << " kg/s\n";
    std::cout << "  rho = " << mixed_stoich.state.rho() << " kg/m³\n";
    std::cout << "  cp = " << mixed_stoich.state.cp() << " J/(mol·K)\n\n";

    std::cout << "Complete Combustion (CO2 + H2O only):\n";
    std::cout << "  T_ad = " << burned_stoich.T << " K\n";
    std::cout << "  rho = " << burned_stoich.rho() << " kg/m³\n";
    std::cout << "  cp = " << burned_stoich.cp() << " J/(mol·K)\n";
    std::cout << "  mu = " << burned_stoich.mu() * 1.0e6 << " μPa·s\n";
    std::cout << "  k = " << burned_stoich.k() << " W/(m·K)\n";
    std::cout << "  Pr = " << burned_stoich.Pr() << "\n";
    std::cout << "  a = " << burned_stoich.a() << " m/s\n\n";

    std::cout << "WGS Equilibrium (CO + H2O <-> CO2 + H2):\n";
    std::cout << "  T_ad = " << wgs_stoich.T << " K\n";
    std::cout << "  rho = " << wgs_stoich.rho() << " kg/m³\n";
    std::cout << "  mu = " << wgs_stoich.mu() * 1.0e6 << " μPa·s\n";
    std::cout << "  Pr = " << wgs_stoich.Pr() << "\n\n";

    // Show product composition
    std::cout << "Product Composition (WGS equilibrium, mole fractions > 0.001):\n";
    for (std::size_t i = 0; i < n_species; ++i) {
        if (wgs_stoich.X[i] > 0.001) {
            std::cout << "  " << std::setw(6) << species_names[i]
                      << ": " << std::setprecision(4) << wgs_stoich.X[i] << "\n";
        }
    }

    return 0;
}
