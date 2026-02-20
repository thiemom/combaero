#include "../include/combustion.h"
#include "../include/humidair.h"
#include "../include/state.h"
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"

#include <iomanip>
#include <iostream>
#include <vector>

// Combustion State Example
//
// Demonstrates:
// 1. combustion_state: phi-based call (complete and equilibrium)
// 2. combustion_state_from_streams: stream-based call
// 3. PressureLossCorrelation hook (loop-free, Newton-solver safe)
// 4. Inverse solvers: set_fuel_stream_for_Tad, set_fuel_stream_for_O2
// 5. Fuel LHV and Bilger mixture fraction

int main()
{
    std::cout << std::fixed << std::setprecision(4);

    const std::size_t n = species_names.size();
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    const std::size_t idx_H2  = species_index_from_name("H2");

    // Fuel: 80% CH4 + 20% H2 (mole fractions)
    std::vector<double> X_fuel(n, 0.0);
    X_fuel[idx_CH4] = 0.8;
    X_fuel[idx_H2]  = 0.2;

    // Oxidizer: standard dry air
    const std::vector<double> X_air = standard_dry_air_composition();

    // =========================================================================
    // 1. combustion_state — phi sweep (complete combustion)
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "1. combustion_state — phi sweep (complete combustion)\n";
    std::cout << "=========================================================================\n\n";

    const double T_react = 700.0;   // K  reactant temperature
    const double P_comb  = 3.0e5;   // Pa combustor pressure

    std::cout << "Fuel: 80% CH4 + 20% H2,  T_react=" << T_react << " K,  P="
              << P_comb / 1e5 << " bar\n\n";

    std::cout << std::setw(6)  << "phi"
              << std::setw(10) << "T_ad [K]"
              << std::setw(10) << "P_out [bar]"
              << std::setw(10) << "Z_Bilger"
              << std::setw(10) << "gamma"
              << "\n";
    std::cout << std::string(46, '-') << "\n";

    for (double phi = 0.5; phi <= 1.51; phi += 0.25) {
        CombustionState cs = combustion_state(X_fuel, X_air, phi, T_react, P_comb,
                                              "CH4/H2", CombustionMethod::Complete);
        std::cout << std::setw(6)  << phi
                  << std::setw(10) << cs.products.thermo.T
                  << std::setw(10) << cs.products.thermo.P / 1e5
                  << std::setw(10) << cs.mixture_fraction
                  << std::setw(10) << cs.products.thermo.gamma
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 2. combustion_state_from_streams — stream-based call
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "2. combustion_state_from_streams\n";
    std::cout << "=========================================================================\n\n";

    // Air stream: 10 kg/s
    Stream air_stream;
    air_stream.state.T = T_react;
    air_stream.state.P = P_comb;
    air_stream.state.X = X_air;
    air_stream.mdot    = 10.0;

    // Fuel stream: mdot set to achieve phi = 0.8
    Stream fuel_base_s;
    fuel_base_s.state.T = T_react;
    fuel_base_s.state.P = P_comb;
    fuel_base_s.state.X = X_fuel;
    Stream fuel_stream = set_fuel_stream_for_phi(0.8, fuel_base_s, air_stream);

    CombustionState cs_streams = combustion_state_from_streams(fuel_stream, air_stream,
                                                                "CH4/H2",
                                                                CombustionMethod::Complete);

    std::cout << "Streams: mdot_fuel=" << fuel_stream.mdot << " kg/s"
              << ",  mdot_air=" << air_stream.mdot << " kg/s\n";
    std::cout << "  phi (computed)  = " << cs_streams.phi << "\n";
    std::cout << "  T_ad            = " << cs_streams.products.thermo.T << " K\n";
    std::cout << "  P_out           = " << cs_streams.products.thermo.P / 1e5 << " bar\n";
    std::cout << "  mixture_fraction= " << cs_streams.mixture_fraction << "\n\n";

    // =========================================================================
    // 3. PressureLossCorrelation hook
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "3. PressureLossCorrelation Hook (loop-free)\n";
    std::cout << "=========================================================================\n\n";

    // Hook: 3% base loss + 0.5% per unit dimensionless temperature rise theta
    // theta = T_ad/T_in - 1  (0 for cold flow, ~2-3 for typical combustion)
    // All fields in PressureLossContext are loop-free — safe for Newton solvers.
    PressureLossCorrelation loss_fn = [](const PressureLossContext& ctx) {
        return 0.03 + 0.005 * ctx.theta;
    };

    std::cout << "Hook: dP/P = 0.03 + 0.005 * theta\n\n";
    std::cout << std::setw(6)  << "phi"
              << std::setw(10) << "theta"
              << std::setw(12) << "dP/P [%]"
              << std::setw(12) << "P_out [bar]"
              << std::setw(10) << "T_ad [K]"
              << "\n";
    std::cout << std::string(50, '-') << "\n";

    for (double phi : {0.5, 0.7, 0.9, 1.0, 1.2}) {
        CombustionState cs_hook = combustion_state(X_fuel, X_air, phi, T_react, P_comb,
                                                   "CH4/H2", CombustionMethod::Complete,
                                                   loss_fn);
        const double theta   = cs_hook.products.thermo.T / T_react - 1.0;
        const double dp_frac = 0.03 + 0.005 * theta;
        std::cout << std::setw(6)  << phi
                  << std::setw(10) << theta
                  << std::setw(12) << dp_frac * 100.0
                  << std::setw(12) << cs_hook.products.thermo.P / 1e5
                  << std::setw(10) << cs_hook.products.thermo.T
                  << "\n";
    }
    std::cout << "\n";

    // Hook also works with streams — mdot_fuel and mdot_air are populated
    PressureLossCorrelation loss_with_mdot = [](const PressureLossContext& ctx) {
        // Could use ctx.mdot_fuel / ctx.mdot_air for fuel-split-dependent losses
        return 0.04 + 0.002 * ctx.theta;
    };

    CombustionState cs_hook_streams = combustion_state_from_streams(
        fuel_stream, air_stream, "CH4/H2", CombustionMethod::Complete, loss_with_mdot);

    std::cout << "Stream hook (dP/P = 0.04 + 0.002*theta):\n";
    std::cout << "  mdot_fuel in ctx = " << fuel_stream.mdot << " kg/s\n";
    std::cout << "  P_out = " << cs_hook_streams.products.thermo.P / 1e5 << " bar\n\n";

    // =========================================================================
    // 4. Inverse solvers
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "4. Inverse Solvers\n";
    std::cout << "=========================================================================\n\n";

    // set_fuel_stream_for_Tad: find fuel mdot to hit a target adiabatic flame temperature
    const double T_ad_target = 1800.0;  // K

    Stream fuel_base;
    fuel_base.state.T = T_react;
    fuel_base.state.P = P_comb;
    fuel_base.state.X = X_fuel;

    Stream fuel_for_Tad = set_fuel_stream_for_Tad(T_ad_target, fuel_base, air_stream);
    CombustionState cs_Tad = combustion_state_from_streams(fuel_for_Tad, air_stream);

    std::cout << "set_fuel_stream_for_Tad (target T_ad=" << T_ad_target << " K):\n";
    std::cout << "  mdot_fuel = " << fuel_for_Tad.mdot << " kg/s\n";
    std::cout << "  phi       = " << cs_Tad.phi << "\n";
    std::cout << "  T_ad      = " << cs_Tad.products.thermo.T << " K\n\n";

    // set_fuel_stream_for_O2: find fuel mdot to hit a target O2 mole fraction in products
    const double X_O2_target = 0.05;  // 5% O2 in products (lean)

    Stream fuel_for_O2 = set_fuel_stream_for_O2(X_O2_target, fuel_base, air_stream);
    CombustionState cs_O2 = combustion_state_from_streams(fuel_for_O2, air_stream);

    std::cout << "set_fuel_stream_for_O2 (target X_O2=" << X_O2_target << "):\n";
    std::cout << "  mdot_fuel = " << fuel_for_O2.mdot << " kg/s\n";
    std::cout << "  phi       = " << cs_O2.phi << "\n";
    std::cout << "  T_ad      = " << cs_O2.products.thermo.T << " K\n";

    // Verify O2 in products by checking phi (lean => O2 present)
    std::cout << "  phi (lean, O2 in products) = " << cs_O2.phi << "\n\n";

    // =========================================================================
    // 5. Fuel LHV and Bilger mixture fraction
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "5. Fuel LHV and Bilger Mixture Fraction\n";
    std::cout << "=========================================================================\n\n";

    // LHV of the 80/20 CH4/H2 blend
    const double lhv_molar = fuel_lhv_molar(X_fuel);
    const double lhv_mass  = fuel_lhv_mass(X_fuel);

    std::cout << "Fuel blend (80% CH4, 20% H2):\n";
    std::cout << "  LHV = " << lhv_molar / 1e6 << " MJ/mol\n";
    std::cout << "  LHV = " << lhv_mass  / 1e6 << " MJ/kg\n\n";

    // Pure species for reference
    std::vector<double> X_pure_CH4(n, 0.0);
    std::vector<double> X_pure_H2(n, 0.0);
    X_pure_CH4[idx_CH4] = 1.0;
    X_pure_H2[idx_H2]   = 1.0;

    std::cout << "  LHV CH4 = " << fuel_lhv_mass(X_pure_CH4) / 1e6 << " MJ/kg\n";
    std::cout << "  LHV H2  = " << fuel_lhv_mass(X_pure_H2)  / 1e6 << " MJ/kg\n\n";

    // Bilger mixture fraction at stoichiometric and at phi=0.8
    const double Z_st = bilger_stoich_mixture_fraction_mass(X_fuel, X_air);
    std::cout << "Bilger stoichiometric mixture fraction Z_st = " << Z_st << "\n";

    // Verify: phi from Z_st should be 1.0
    const double phi_from_Zst = equivalence_ratio_from_bilger_Z_mass(Z_st, X_fuel, X_air);
    std::cout << "phi from Z_st = " << phi_from_Zst << "  (should be 1.0)\n";

    // Z at phi=0.8
    const double Z_08 = bilger_Z_from_equivalence_ratio_mass(0.8, X_fuel, X_air);
    std::cout << "Z at phi=0.8  = " << Z_08 << "\n\n";

    // Bilger mixture fraction from mole fractions (convenience overload)
    const std::vector<double> X_mix_stoich = set_equivalence_ratio_mole(1.0, X_fuel, X_air);
    const double Z_bilger = bilger_mixture_fraction_from_moles(X_mix_stoich, X_fuel, X_air);
    std::cout << "Bilger Z (from mole fractions at phi=1.0) = " << Z_bilger
              << "  (should equal Z_st=" << Z_st << ")\n\n";

    return 0;
}
