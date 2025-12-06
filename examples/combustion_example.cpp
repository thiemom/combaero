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
 * Combustion Example: Fuel + Humid Air -> Mix -> Combust -> Equilibrium
 *
 * Demonstrates:
 * - Using set_fuel_stream_for_phi to configure fuel stream for target phi
 * - Stream mixing with enthalpy balance
 * - combustion_equilibrium: one-step combustion + reforming + WGS equilibrium
 * - Varying equivalence ratio from lean to rich
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
    // Sweep equivalence ratio from 0.5 to 1.2 (lean to slightly rich)
    // =========================================================================

    std::cout << "Equivalence Ratio Sweep\n";
    std::cout << "-----------------------\n\n";

    std::cout << std::setw(6) << "phi"
              << std::setw(10) << "mdot_f"
              << std::setw(10) << "T_mix"
              << std::setw(10) << "T_eq"
              << std::setw(12) << "mu_eq"
              << std::setw(12) << "Pr_eq"
              << "\n";
    std::cout << std::setw(6) << "[-]"
              << std::setw(10) << "[kg/s]"
              << std::setw(10) << "[K]"
              << std::setw(10) << "[K]"
              << std::setw(12) << "[uPa.s]"
              << std::setw(12) << "[-]"
              << "\n";
    std::cout << std::string(60, '-') << "\n";

    for (double phi = 0.5; phi <= 1.21; phi += 0.1) {
        // Use set_fuel_stream_for_phi to get fuel stream with correct mdot
        Stream fuel_phi = set_fuel_stream_for_phi(phi, fuel, air);

        // Mix streams (enthalpy balance gives mixed temperature)
        Stream mixed = mix({fuel_phi, air});

        // One-step: combustion + reforming + WGS equilibrium
        State eq = combustion_equilibrium(mixed.state);

        // Output results
        std::cout << std::setw(6) << phi
                  << std::setw(10) << fuel_phi.mdot
                  << std::setw(10) << mixed.state.T
                  << std::setw(10) << eq.T
                  << std::setw(12) << eq.mu() * 1.0e6  // convert to μPa·s
                  << std::setw(12) << eq.Pr()
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
    
    // One-step: combustion + reforming + WGS equilibrium
    State eq_stoich = combustion_equilibrium(mixed_stoich.state);
    
    // For comparison, also show complete combustion (before equilibrium)
    State burned_stoich = complete_combustion(mixed_stoich.state);

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

    std::cout << "Reforming + WGS Equilibrium:\n";
    std::cout << "  T_eq = " << eq_stoich.T << " K\n";
    std::cout << "  rho = " << eq_stoich.rho() << " kg/m³\n";
    std::cout << "  mu = " << eq_stoich.mu() * 1.0e6 << " μPa·s\n";
    std::cout << "  Pr = " << eq_stoich.Pr() << "\n\n";

    // Show product composition
    std::cout << "Product Composition (Reforming+WGS equilibrium, mole fractions > 0.001):\n";
    for (std::size_t i = 0; i < n_species; ++i) {
        if (eq_stoich.X[i] > 0.001) {
            std::cout << "  " << std::setw(6) << species_names[i]
                      << ": " << std::setprecision(4) << eq_stoich.X[i] << "\n";
        }
    }

    // =========================================================================
    // Rich case (phi = 1.2) - Reforming converts hydrocarbons to CO + H2
    // =========================================================================

    std::cout << "\n=========================================================================\n";
    std::cout << "Rich Case (phi = 1.2) - Reforming + WGS Equilibrium\n";
    std::cout << "=========================================================================\n\n";

    Stream fuel_rich = set_fuel_stream_for_phi(1.2, fuel, air);
    Stream mixed_rich = mix({fuel_rich, air});
    
    // One-step: combustion + reforming + WGS equilibrium
    State eq_rich = combustion_equilibrium(mixed_rich.state);
    
    // For comparison, also show complete combustion (before equilibrium)
    State burned_rich = complete_combustion(mixed_rich.state);

    std::cout << "Complete Combustion (before equilibrium):\n";
    std::cout << "  T_ad = " << burned_rich.T << " K\n";
    for (std::size_t i = 0; i < n_species; ++i) {
        if (burned_rich.X[i] > 0.001) {
            std::cout << "  " << std::setw(6) << species_names[i]
                      << ": " << std::setprecision(4) << burned_rich.X[i] << "\n";
        }
    }

    std::cout << "\nReforming + WGS Equilibrium (hydrocarbons reformed to CO + H2):\n";
    std::cout << "  T_eq = " << eq_rich.T << " K\n";
    for (std::size_t i = 0; i < n_species; ++i) {
        if (eq_rich.X[i] > 0.001) {
            std::cout << "  " << std::setw(6) << species_names[i]
                      << ": " << std::setprecision(4) << eq_rich.X[i] << "\n";
        }
    }

    return 0;
}
