#include "../include/combustion.h"
#include "../include/humidair.h"
#include "../include/stagnation.h"
#include "../include/state.h"
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"
#include "../include/transport.h"

#include <iomanip>
#include <iostream>
#include <vector>

// Stagnation Example
//
// Demonstrates:
// 1. Static <-> stagnation conversions (T0, P0) at varying Mach numbers
// 2. Round-trip accuracy check
// 3. Variable-cp vs ideal-gas (constant gamma) comparison
// 4. Adiabatic wall temperature: laminar vs turbulent recovery
// 5. Practical case: turbine inlet stagnation from combustor exit

int main()
{
    std::cout << std::fixed << std::setprecision(4);

    const std::vector<double> X_air = standard_dry_air_composition();

    // =========================================================================
    // 1. Static -> stagnation conversions over Mach sweep
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "1. Static -> Stagnation Conversions (Mach sweep)\n";
    std::cout << "=========================================================================\n\n";

    const double T_static = 800.0;   // K  turbine inlet static temperature
    const double P_static = 4.0e5;   // Pa static pressure

    std::cout << "T_static=" << T_static << " K,  P_static=" << P_static / 1e5 << " bar\n\n";
    std::cout << std::setw(8)  << "M [-]"
              << std::setw(12) << "T0 [K]"
              << std::setw(12) << "P0 [bar]"
              << std::setw(12) << "T0/T [-]"
              << std::setw(12) << "P0/P [-]"
              << "\n";
    std::cout << std::string(56, '-') << "\n";

    for (double M : {0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0}) {
        const double T0 = T0_from_static(T_static, M, X_air);
        const double P0 = P0_from_static(P_static, T_static, M, X_air);
        std::cout << std::setw(8)  << M
                  << std::setw(12) << T0
                  << std::setw(12) << P0 / 1e5
                  << std::setw(12) << T0 / T_static
                  << std::setw(12) << P0 / P_static
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 2. Round-trip accuracy
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "2. Round-Trip Accuracy (static -> stagnation -> static)\n";
    std::cout << "=========================================================================\n\n";

    std::cout << std::setw(8)  << "M [-]"
              << std::setw(14) << "T_err [K]"
              << std::setw(14) << "T_err [%]"
              << std::setw(14) << "P_err [Pa]"
              << std::setw(14) << "P_err [%]"
              << "\n";
    std::cout << std::string(64, '-') << "\n";

    for (double M : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        const double T0 = T0_from_static(T_static, M, X_air);
        const double P0 = P0_from_static(P_static, T_static, M, X_air);
        const double T_back = T_from_stagnation(T0, M, X_air);
        const double P_back = P_from_stagnation(P0, T0, M, X_air);
        std::cout << std::setw(8)  << M
                  << std::setw(14) << std::abs(T_back - T_static)
                  << std::setw(14) << std::abs(T_back - T_static) / T_static * 100.0
                  << std::setw(14) << std::abs(P_back - P_static)
                  << std::setw(14) << std::abs(P_back - P_static) / P_static * 100.0
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 3. Variable-cp vs constant-gamma comparison
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "3. Variable-cp vs Constant-gamma (T0 comparison)\n";
    std::cout << "=========================================================================\n\n";

    std::cout << "Constant-gamma uses gamma at T_static; variable-cp iterates on h(T).\n\n";
    std::cout << std::setw(8)  << "T [K]"
              << std::setw(8)  << "M [-]"
              << std::setw(14) << "T0_varcp [K]"
              << std::setw(14) << "T0_constg [K]"
              << std::setw(12) << "diff [K]"
              << "\n";
    std::cout << std::string(56, '-') << "\n";

    for (double T_s : {300.0, 600.0, 1000.0, 1500.0, 2000.0}) {
        for (double M_s : {0.3, 0.8}) {
            const double T0_var = T0_from_static(T_s, M_s, X_air);

            // Constant-gamma approximation: T0 = T * (1 + (gamma-1)/2 * M^2)
            const double gamma_s = isentropic_expansion_coefficient(T_s, X_air);
            const double T0_cg = T_s * (1.0 + 0.5 * (gamma_s - 1.0) * M_s * M_s);

            std::cout << std::setw(8)  << T_s
                      << std::setw(8)  << M_s
                      << std::setw(14) << T0_var
                      << std::setw(14) << T0_cg
                      << std::setw(12) << T0_var - T0_cg
                      << "\n";
        }
    }
    std::cout << "\n";

    // =========================================================================
    // 4. Adiabatic wall temperature: laminar vs turbulent
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "4. Adiabatic Wall Temperature (laminar vs turbulent)\n";
    std::cout << "=========================================================================\n\n";

    const double T_bulk = 1200.0;  // K  hot gas static temperature
    const double P_bulk = 8.0e5;   // Pa
    const double Pr_val = prandtl(T_bulk, P_bulk, X_air);

    std::cout << "T_static=" << T_bulk << " K,  P=" << P_bulk / 1e5 << " bar\n";
    std::cout << "Pr=" << Pr_val << "\n\n";

    std::cout << std::setw(8)  << "M [-]"
              << std::setw(12) << "T0 [K]"
              << std::setw(14) << "T_aw_turb [K]"
              << std::setw(14) << "T_aw_lam [K]"
              << std::setw(12) << "r_turb"
              << std::setw(12) << "r_lam"
              << "\n";
    std::cout << std::string(72, '-') << "\n";

    for (double M_aw : {0.1, 0.2, 0.3, 0.5, 0.7, 0.9}) {
        const double T0_aw   = T0_from_static(T_bulk, M_aw, X_air);
        const double T_aw_t  = T_adiabatic_wall_mach(T_bulk, M_aw, T_bulk, P_bulk, X_air, true);
        const double T_aw_l  = T_adiabatic_wall_mach(T_bulk, M_aw, T_bulk, P_bulk, X_air, false);
        const double r_turb  = recovery_factor(Pr_val, true);
        const double r_lam   = recovery_factor(Pr_val, false);
        std::cout << std::setw(8)  << M_aw
                  << std::setw(12) << T0_aw
                  << std::setw(14) << T_aw_t
                  << std::setw(14) << T_aw_l
                  << std::setw(12) << r_turb
                  << std::setw(12) << r_lam
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 5. Practical case: turbine inlet from combustor exit
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "5. Practical Case: Turbine Inlet from Combustor Exit\n";
    std::cout << "=========================================================================\n\n";

    // Combustor: CH4/air at phi=0.6 (lean, typical gas turbine)
    const std::size_t idx_CH4 = species_index_from_name("CH4");
    std::vector<double> X_CH4(species_names.size(), 0.0);
    X_CH4[idx_CH4] = 1.0;

    Stream fuel_t;
    fuel_t.state.T = 600.0;   // K  preheated fuel
    fuel_t.state.P = 20.0e5;  // Pa
    fuel_t.state.X = X_CH4;

    Stream air_t;
    air_t.state.T = 700.0;    // K  compressed air
    air_t.state.P = 20.0e5;
    air_t.state.X = X_air;
    air_t.mdot    = 50.0;     // kg/s

    Stream fuel_lean = set_fuel_stream_for_phi(0.6, fuel_t, air_t);
    Stream mixed_t   = mix({fuel_lean, air_t});
    State  burned_t  = complete_combustion(mixed_t.state);

    std::cout << "Combustor: CH4/air, phi=0.6, P=" << fuel_t.state.P / 1e5 << " bar\n";
    std::cout << "  T_ad (combustor exit) = " << burned_t.T << " K\n";
    std::cout << "  mdot_total            = " << (fuel_lean.mdot + air_t.mdot) << " kg/s\n\n";

    // Turbine nozzle guide vane: M=0.9 at exit
    const double M_nozzle = 0.9;
    const double T0_turbine = burned_t.T;  // stagnation T = T_ad (low velocity in combustor)
    const double P0_turbine = burned_t.P;  // stagnation P = combustor pressure

    const double T_nozzle_exit = T_from_stagnation(T0_turbine, M_nozzle, burned_t.X);
    const double P_nozzle_exit = P_from_stagnation(P0_turbine, T0_turbine, M_nozzle, burned_t.X);
    const double T_aw_nozzle   = T_adiabatic_wall_mach(T_nozzle_exit, M_nozzle,
                                                         T_nozzle_exit, P_nozzle_exit,
                                                         burned_t.X, true);

    std::cout << "Turbine nozzle guide vane (M=" << M_nozzle << "):\n";
    std::cout << "  T0_in   = " << T0_turbine << " K  (= T_ad, combustor exit)\n";
    std::cout << "  P0_in   = " << P0_turbine / 1e5 << " bar\n";
    std::cout << "  T_exit  = " << T_nozzle_exit << " K  (static)\n";
    std::cout << "  P_exit  = " << P_nozzle_exit / 1e5 << " bar  (static)\n";
    std::cout << "  T_aw    = " << T_aw_nozzle << " K  (use as hot-side T for cooling calc)\n";
    std::cout << "  dT_aw   = " << T_aw_nozzle - T_nozzle_exit
              << " K  (T_aw - T_static, recovery heating)\n\n";

    // Stagnation enthalpy and velocity
    const double h_static_mass = h(T_nozzle_exit, burned_t.X) / burned_t.mw() * 1000.0;
    const double a_exit = speed_of_sound(T_nozzle_exit, burned_t.X);
    const double v_exit = M_nozzle * a_exit;
    const double h0_val = h0_from_static(h_static_mass, v_exit);

    std::cout << "  v_exit  = " << v_exit << " m/s\n";
    std::cout << "  a_exit  = " << a_exit << " m/s\n";
    std::cout << "  h0      = " << h0_val / 1e3 << " kJ/kg\n\n";

    return 0;
}
