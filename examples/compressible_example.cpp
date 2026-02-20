#include "../include/combustion.h"
#include "../include/compressible.h"
#include "../include/humidair.h"
#include "../include/stagnation.h"
#include "../include/state.h"
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"

#include <iomanip>
#include <iostream>
#include <vector>

// Compressible Flow Example
//
// Demonstrates:
// 1. Isentropic nozzle flow: forward problem (given geometry, find mdot)
// 2. Inverse problems: solve for area, back pressure, or stagnation pressure
// 3. Converging-diverging nozzle with axial profile
// 4. Fanno pipe flow (adiabatic friction)
// 5. Rocket nozzle thrust calculation
// 6. Stagnation / static conversions and adiabatic wall temperature

int main()
{
    std::cout << std::fixed << std::setprecision(4);

    const std::vector<double> X_air = standard_dry_air_composition();

    // =========================================================================
    // 1. Isentropic nozzle flow — forward problem
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "1. Isentropic Nozzle Flow (forward problem)\n";
    std::cout << "=========================================================================\n\n";

    const double T0    = 800.0;     // K  stagnation temperature
    const double P0    = 5.0e5;     // Pa stagnation pressure
    const double A_eff = 1.5e-4;    // m² effective throat area (A * Cd)

    // Subsonic case: back pressure above critical
    const double P_back_sub = 4.0e5;
    CompressibleFlowSolution sol_sub = nozzle_flow(T0, P0, P_back_sub, A_eff, X_air);

    std::cout << "Subsonic (P_back = " << P_back_sub / 1e5 << " bar):\n";
    std::cout << "  mdot   = " << sol_sub.mdot << " kg/s\n";
    std::cout << "  M_exit = " << sol_sub.M << "\n";
    std::cout << "  T_exit = " << sol_sub.outlet.T << " K\n";
    std::cout << "  P_exit = " << sol_sub.outlet.P / 1e5 << " bar\n";
    std::cout << "  choked = " << (sol_sub.choked ? "yes" : "no") << "\n\n";

    // Choked case: back pressure at or below critical
    const double P_back_choked = 2.0e5;
    CompressibleFlowSolution sol_choked = nozzle_flow(T0, P0, P_back_choked, A_eff, X_air);

    std::cout << "Choked (P_back = " << P_back_choked / 1e5 << " bar):\n";
    std::cout << "  mdot   = " << sol_choked.mdot << " kg/s\n";
    std::cout << "  M_exit = " << sol_choked.M << "\n";
    std::cout << "  T_exit = " << sol_choked.outlet.T << " K\n";
    std::cout << "  P_exit = " << sol_choked.outlet.P / 1e5 << " bar\n";
    std::cout << "  choked = " << (sol_choked.choked ? "yes" : "no") << "\n\n";

    // Critical pressure ratio
    const double pr_crit = critical_pressure_ratio(T0, P0, X_air);
    std::cout << "Critical pressure ratio P*/P0 = " << pr_crit << "\n";
    std::cout << "Critical back pressure       = " << pr_crit * P0 / 1e5 << " bar\n\n";

    // =========================================================================
    // 2. Inverse problems
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "2. Inverse Problems\n";
    std::cout << "=========================================================================\n\n";

    const double mdot_target = sol_sub.mdot;

    // Find area for given mdot
    const double A_solved = solve_A_eff_from_mdot(T0, P0, P_back_sub, mdot_target, X_air);
    std::cout << "solve_A_eff_from_mdot: A_eff = " << A_solved * 1e4 << " cm²"
              << "  (input: " << A_eff * 1e4 << " cm²)\n";

    // Find back pressure for given mdot
    const double P_back_solved = solve_P_back_from_mdot(T0, P0, A_eff, mdot_target, X_air);
    std::cout << "solve_P_back_from_mdot: P_back = " << P_back_solved / 1e5 << " bar"
              << "  (input: " << P_back_sub / 1e5 << " bar)\n";

    // Find stagnation pressure for given mdot
    const double P0_solved = solve_P0_from_mdot(T0, P_back_sub, A_eff, mdot_target, X_air);
    std::cout << "solve_P0_from_mdot: P0 = " << P0_solved / 1e5 << " bar"
              << "  (input: " << P0 / 1e5 << " bar)\n\n";

    // =========================================================================
    // 3. Converging-diverging nozzle with axial profile
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "3. Converging-Diverging Nozzle (quasi-1D profile)\n";
    std::cout << "=========================================================================\n\n";

    // Geometry: inlet 4 cm², throat 1 cm², exit 2.5 cm²
    const double A_in     = 4.0e-4;   // m²
    const double A_throat = 1.0e-4;   // m²
    const double A_exit   = 2.5e-4;   // m²
    const double x_throat = 0.05;     // m
    const double x_exit   = 0.12;     // m

    // Subsonic solution (P_exit > P_critical)
    NozzleSolution noz = nozzle_cd(T0, P0, 3.5e5, A_in, A_throat, A_exit,
                                   x_throat, x_exit, X_air, 80);

    std::cout << "Converging-diverging nozzle (subsonic, P_exit=3.5 bar):\n";
    std::cout << "  mdot     = " << noz.mdot << " kg/s\n";
    std::cout << "  choked   = " << (noz.choked ? "yes" : "no") << "\n";
    std::cout << "  T0       = " << noz.T0 << " K\n";
    std::cout << "  P0       = " << noz.P0 / 1e5 << " bar\n";
    std::cout << "  T_exit   = " << noz.outlet.T << " K\n";
    std::cout << "  P_exit   = " << noz.outlet.P / 1e5 << " bar\n\n";

    // Print axial profile (every 10th station)
    std::cout << std::setw(8) << "x [m]"
              << std::setw(10) << "A [cm²]"
              << std::setw(8) << "M [-]"
              << std::setw(10) << "T [K]"
              << std::setw(10) << "P [bar]"
              << "\n";
    std::cout << std::string(46, '-') << "\n";
    for (std::size_t i = 0; i < noz.profile.size(); i += 10) {
        const auto& st = noz.profile[i];
        std::cout << std::setw(8) << st.x
                  << std::setw(10) << st.A * 1e4
                  << std::setw(8) << st.M
                  << std::setw(10) << st.T
                  << std::setw(10) << st.P / 1e5
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 4. Fanno pipe flow (adiabatic friction)
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "4. Fanno Pipe Flow (adiabatic, friction)\n";
    std::cout << "=========================================================================\n\n";

    const double T_in  = 600.0;   // K  inlet static temperature
    const double P_in  = 3.0e5;   // Pa inlet static pressure
    const double u_in  = 80.0;    // m/s inlet velocity
    const double L     = 2.0;     // m  pipe length
    const double D     = 0.02;    // m  pipe diameter
    const double rough = 5.0e-5;  // m  wall roughness (commercial steel)

    FannoSolution fanno = fanno_pipe_rough(T_in, P_in, u_in, L, D, rough, X_air,
                                           "haaland", 200, true);

    std::cout << "Fanno pipe (L=" << L << " m, D=" << D * 1e3 << " mm, rough="
              << rough * 1e6 << " μm):\n";
    std::cout << "  Re_in  = " << fanno.Re_in << "\n";
    std::cout << "  f_avg  = " << fanno.f_avg << "\n";
    std::cout << "  choked = " << (fanno.choked ? "yes" : "no") << "\n";
    std::cout << "  T_out  = " << fanno.outlet.T << " K\n";
    std::cout << "  P_out  = " << fanno.outlet.P / 1e5 << " bar\n";
    std::cout << "  dP     = " << (P_in - fanno.outlet.P) / 1e3 << " kPa\n\n";

    // Maximum pipe length before choking
    const double f_est = fanno.f_avg;
    const double L_max = fanno_max_length(T_in, P_in, u_in, D, f_est, X_air);
    std::cout << "  L_max (before choking) = " << L_max << " m\n\n";

    // Axial profile (every 20th station)
    std::cout << std::setw(8) << "x [m]"
              << std::setw(8) << "M [-]"
              << std::setw(10) << "T [K]"
              << std::setw(10) << "P [bar]"
              << std::setw(8) << "Re [-]"
              << "\n";
    std::cout << std::string(44, '-') << "\n";
    for (std::size_t i = 0; i < fanno.profile.size(); i += 20) {
        const auto& st = fanno.profile[i];
        std::cout << std::setw(8) << st.x
                  << std::setw(8) << st.M
                  << std::setw(10) << st.T
                  << std::setw(10) << st.P / 1e5
                  << std::setw(8) << static_cast<int>(st.Re)
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 5. Rocket nozzle thrust
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "5. Rocket Nozzle Thrust\n";
    std::cout << "=========================================================================\n\n";

    // Chamber gas: H2 + air combustion at phi=1.0, P=70 bar
    // This gives real product composition (H2O, N2, O2 traces) and T_ad as T0.
    const double P0_rocket = 70.0e5;   // Pa  chamber pressure
    const std::size_t idx_H2 = species_index_from_name("H2");

    std::vector<double> X_H2(species_names.size(), 0.0);
    X_H2[idx_H2] = 1.0;

    Stream fuel_r;
    fuel_r.state.T = 300.0;
    fuel_r.state.P = P0_rocket;
    fuel_r.state.X = X_H2;

    Stream air_r;
    air_r.state.T = 300.0;
    air_r.state.P = P0_rocket;
    air_r.state.X = standard_dry_air_composition();
    air_r.mdot    = 1.0;  // kg/s reference

    Stream fuel_stoich_r = set_fuel_stream_for_phi(1.0, fuel_r, air_r);

    // Mix reactants, combust to get product composition and T_ad
    Stream mixed_r         = mix({fuel_stoich_r, air_r});
    State  burned_r        = complete_combustion(mixed_r.state);
    const double T0_rocket         = burned_r.T;    // adiabatic flame T = chamber T0
    const std::vector<double>& X_products = burned_r.X;  // product mole fractions

    std::cout << "Chamber gas (H2/air, phi=1.0, P0=" << P0_rocket / 1e5 << " bar):\n";
    std::cout << "  T_ad (= T0) = " << T0_rocket << " K\n";
    std::cout << "  mdot_total  = " << (fuel_stoich_r.mdot + air_r.mdot) << " kg/s\n\n";

    const double P_amb_sl  = 101325.0;  // Pa  sea level
    const double P_amb_vac = 500.0;     // Pa  near-vacuum design point (P_amb=0 for thrust)

    const double A_in_r  = 1.0e-2;   // m²
    const double A_th_r  = 5.0e-4;   // m²
    const double A_ex_r  = 3.0e-3;   // m²
    const double x_th_r  = 0.15;     // m
    const double x_ex_r  = 0.40;     // m

    // nozzle_thrust(T0, P0, P_design, P_amb, ...)
    // P_design: nozzle exit pressure for Mach profile (must be > 0)
    // P_amb:    ambient for thrust equation (may be 0 for vacuum)
    ThrustResult thr_sl  = nozzle_thrust(T0_rocket, P0_rocket, P_amb_sl,  P_amb_sl,
                                          A_in_r, A_th_r, A_ex_r, x_th_r, x_ex_r, X_products);
    ThrustResult thr_vac = nozzle_thrust(T0_rocket, P0_rocket, P_amb_vac, 0.0,
                                          A_in_r, A_th_r, A_ex_r, x_th_r, x_ex_r, X_products);

    std::cout << "Rocket nozzle (H2/air, T0=" << T0_rocket << " K, P0=" << P0_rocket / 1e5 << " bar):\n\n";
    std::cout << std::setw(20) << "" << std::setw(14) << "Sea level" << std::setw(14) << "Vacuum\n";
    std::cout << std::string(48, '-') << "\n";
    std::cout << std::setw(20) << "Thrust [N]"
              << std::setw(14) << thr_sl.thrust
              << std::setw(14) << thr_vac.thrust << "\n";
    std::cout << std::setw(20) << "Isp [s]"
              << std::setw(14) << thr_sl.specific_impulse
              << std::setw(14) << thr_vac.specific_impulse << "\n";
    std::cout << std::setw(20) << "C_F [-]"
              << std::setw(14) << thr_sl.thrust_coefficient
              << std::setw(14) << thr_vac.thrust_coefficient << "\n";
    std::cout << std::setw(20) << "mdot [kg/s]"
              << std::setw(14) << thr_sl.mdot
              << std::setw(14) << thr_vac.mdot << "\n";
    std::cout << std::setw(20) << "u_exit [m/s]"
              << std::setw(14) << thr_sl.u_exit
              << std::setw(14) << thr_vac.u_exit << "\n\n";

    // =========================================================================
    // 6. Stagnation / static conversions and adiabatic wall temperature
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "6. Stagnation Conversions and Adiabatic Wall Temperature\n";
    std::cout << "=========================================================================\n\n";

    const double T_static = 700.0;   // K
    const double P_static = 2.0e5;   // Pa
    const double M_flow   = 0.5;     // Mach number

    const double T0_calc = T0_from_static(T_static, M_flow, X_air);
    const double P0_calc = P0_from_static(P_static, T_static, M_flow, X_air);

    std::cout << "Static -> Stagnation (M=" << M_flow << "):\n";
    std::cout << "  T_static = " << T_static << " K  ->  T0 = " << T0_calc << " K\n";
    std::cout << "  P_static = " << P_static / 1e5 << " bar  ->  P0 = " << P0_calc / 1e5 << " bar\n\n";

    // Round-trip check
    const double T_back = T_from_stagnation(T0_calc, M_flow, X_air);
    const double P_back2 = P_from_stagnation(P0_calc, T0_calc, M_flow, X_air);
    std::cout << "Round-trip check:\n";
    std::cout << "  T0 -> T_static = " << T_back << " K  (error: "
              << std::abs(T_back - T_static) / T_static * 100.0 << " %)\n";
    std::cout << "  P0 -> P_static = " << P_back2 / 1e5 << " bar  (error: "
              << std::abs(P_back2 - P_static) / P_static * 100.0 << " %)\n\n";

    // Adiabatic wall temperature (important for turbine cooling at high Mach)
    const double T_aw = T_adiabatic_wall(T_static, M_flow * speed_of_sound(T_static, X_air),
                                          T_static, P_static, X_air, true);
    std::cout << "Adiabatic wall temperature (turbulent, M=" << M_flow << "):\n";
    std::cout << "  T_static = " << T_static << " K\n";
    std::cout << "  T0       = " << T0_calc << " K\n";
    std::cout << "  T_aw     = " << T_aw << " K  (use as hot-side T in htc calculations)\n\n";

    return 0;
}
