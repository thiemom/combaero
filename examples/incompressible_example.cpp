#include "../include/friction.h"
#include "../include/humidair.h"
#include "../include/incompressible.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Incompressible Flow Example
//
// Demonstrates:
// 1. Scalar primitives: Bernoulli, orifice, pipe Darcy-Weisbach
// 2. Pressure loss coefficient (zeta / Cd conversions)
// 3. Thermo-aware orifice flow with fixed Cd and user-supplied Cd hook
// 4. Thermo-aware pipe flow with roughness correlation and K-loss hook
// 5. Friction factor correlations (Haaland, Serghides, Colebrook)

int main()
{
    std::cout << std::fixed << std::setprecision(4);

    const std::vector<double> X_air = standard_dry_air_composition();

    // =========================================================================
    // 1. Scalar primitives
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "1. Scalar Primitives (Bernoulli, Orifice, Pipe)\n";
    std::cout << "=========================================================================\n\n";

    const double rho = 1.2;     // kg/m³  (air at ~300 K, 1 atm)
    const double P1  = 101325.0;
    const double v1  = 5.0;     // m/s upstream velocity
    const double v2  = 20.0;    // m/s downstream velocity (area contraction)

    const double P2_bern = bernoulli_P2(P1, v1, v2, rho);
    std::cout << "Bernoulli: P1=" << P1 << " Pa, v1=" << v1 << " m/s, v2=" << v2 << " m/s\n";
    std::cout << "  P2 = " << P2_bern << " Pa  (dP = " << P1 - P2_bern << " Pa)\n\n";

    // Orifice primitives
    const double A_or  = 1.0e-4;  // m²  orifice area
    const double Cd_or = 0.62;    // sharp-edge discharge coefficient
    const double dP_or = 5000.0;  // Pa  differential pressure

    const double mdot_or = orifice_mdot(P1, P1 - dP_or, A_or, Cd_or, rho);
    const double Q_or    = orifice_Q(P1, P1 - dP_or, A_or, Cd_or, rho);
    const double dP_back = orifice_dP(mdot_or, A_or, Cd_or, rho);

    std::cout << "Orifice (A=" << A_or * 1e4 << " cm², Cd=" << Cd_or << ", dP=" << dP_or << " Pa):\n";
    std::cout << "  mdot = " << mdot_or << " kg/s\n";
    std::cout << "  Q    = " << Q_or * 1e3 << " L/s\n";
    std::cout << "  dP (back-calc) = " << dP_back << " Pa\n\n";

    // Pipe Darcy-Weisbach
    const double v_pipe = 3.0;    // m/s
    const double L_pipe = 5.0;    // m
    const double D_pipe = 0.025;  // m
    const double f_pipe = 0.02;   // Darcy friction factor

    const double dP_pipe = pipe_dP(v_pipe, L_pipe, D_pipe, f_pipe, rho);
    std::cout << "Pipe Darcy-Weisbach (v=" << v_pipe << " m/s, L=" << L_pipe
              << " m, D=" << D_pipe * 1e3 << " mm, f=" << f_pipe << "):\n";
    std::cout << "  dP = " << dP_pipe << " Pa\n\n";

    // =========================================================================
    // 2. Pressure loss coefficient conversions
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "2. Pressure Loss Coefficient (zeta / Cd)\n";
    std::cout << "=========================================================================\n\n";

    for (double Cd : {0.5, 0.62, 0.7, 0.8, 1.0}) {
        const double zeta = zeta_from_Cd(Cd);
        const double Cd_back = Cd_from_zeta(zeta);
        std::cout << "  Cd=" << Cd << "  ->  zeta=" << zeta
                  << "  ->  Cd_back=" << Cd_back << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 3. Thermo-aware orifice flow
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "3. Thermo-Aware Orifice Flow\n";
    std::cout << "=========================================================================\n\n";

    const double T_or   = 600.0;    // K  hot gas temperature
    const double P_or   = 3.0e5;    // Pa upstream pressure
    const double P_back = 2.8e5;    // Pa downstream pressure
    const double A_hole = 2.0e-4;   // m² orifice area

    // Fixed Cd
    IncompressibleFlowSolution sol_fixed = orifice_flow_thermo(T_or, P_or, X_air,
                                                                P_back, A_hole, 0.62);
    std::cout << "Fixed Cd=0.62:\n";
    std::cout << "  mdot = " << sol_fixed.mdot << " kg/s\n";
    std::cout << "  v    = " << sol_fixed.v << " m/s\n";
    std::cout << "  Re   = " << sol_fixed.Re << "\n";
    std::cout << "  rho  = " << sol_fixed.rho << " kg/m³\n\n";

    // User-supplied Cd hook: Re-dependent (simplified Idelchik-style)
    // Cd rises from ~0.6 at low Re to ~0.65 at high Re
    IncompressibleCdFn cd_re = [](double, double,
                                   const std::vector<double>&, double Re) {
        return 0.65 - 0.05 * std::exp(-Re / 5000.0);
    };

    IncompressibleFlowSolution sol_hook = orifice_flow_thermo(T_or, P_or, X_air,
                                                               P_back, A_hole, cd_re);
    std::cout << "Re-dependent Cd hook (Cd = 0.65 - 0.05*exp(-Re/5000)):\n";
    std::cout << "  mdot = " << sol_hook.mdot << " kg/s\n";
    std::cout << "  v    = " << sol_hook.v << " m/s\n";
    std::cout << "  Re   = " << sol_hook.Re << "\n";
    std::cout << "  Cd   = " << sol_hook.f << "  (stored in .f field)\n\n";

    // Sweep: compare fixed vs hook over a range of dP
    std::cout << "dP sweep (fixed Cd=0.62 vs Re-dependent hook):\n";
    std::cout << std::setw(10) << "dP [Pa]"
              << std::setw(12) << "mdot_fixed"
              << std::setw(12) << "mdot_hook"
              << std::setw(10) << "Cd_hook"
              << "\n";
    std::cout << std::string(44, '-') << "\n";
    for (double dP : {1000.0, 5000.0, 10000.0, 20000.0, 50000.0}) {
        auto s1 = orifice_flow_thermo(T_or, P_or, X_air, P_or - dP, A_hole, 0.62);
        auto s2 = orifice_flow_thermo(T_or, P_or, X_air, P_or - dP, A_hole, cd_re);
        std::cout << std::setw(10) << dP
                  << std::setw(12) << s1.mdot
                  << std::setw(12) << s2.mdot
                  << std::setw(10) << s2.f
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 4. Thermo-aware pipe flow with roughness and K-loss hook
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "4. Thermo-Aware Pipe Flow (roughness + K-loss hook)\n";
    std::cout << "=========================================================================\n\n";

    const double T_pipe  = 400.0;    // K
    const double P_pipe  = 2.0e5;    // Pa
    const double u_flow  = 15.0;     // m/s
    const double L_seg   = 3.0;      // m
    const double D_seg   = 0.03;     // m
    const double rough   = pipe_roughness("commercial_steel");

    std::cout << "Pipe roughness (commercial steel): " << rough * 1e6 << " μm\n\n";

    // Smooth pipe (roughness = 0)
    IncompressibleFlowSolution sol_smooth = pipe_flow_rough(T_pipe, P_pipe, X_air,
                                                             u_flow, L_seg, D_seg);
    // Rough pipe
    IncompressibleFlowSolution sol_rough = pipe_flow_rough(T_pipe, P_pipe, X_air,
                                                            u_flow, L_seg, D_seg, rough);

    std::cout << "Smooth pipe (roughness=0):\n";
    std::cout << "  dP = " << sol_smooth.dP << " Pa,  f = " << sol_smooth.f
              << ",  Re = " << sol_smooth.Re << "\n\n";
    std::cout << "Rough pipe (commercial steel):\n";
    std::cout << "  dP = " << sol_rough.dP << " Pa,  f = " << sol_rough.f
              << ",  Re = " << sol_rough.Re << "\n\n";

    // K-loss hook: 90° elbow (K ≈ 0.9) + sudden expansion (K ≈ 0.5)
    // Returns total additional K independent of Re for this simple model
    IncompressibleKLossFn k_fittings = [](double, double,
                                           const std::vector<double>&, double) {
        const double K_elbow     = 0.9;
        const double K_expansion = 0.5;
        return K_elbow + K_expansion;
    };

    IncompressibleFlowSolution sol_fittings = pipe_flow_rough(T_pipe, P_pipe, X_air,
                                                               u_flow, L_seg, D_seg,
                                                               rough, "haaland", k_fittings);
    std::cout << "Rough pipe + fittings (K_elbow=0.9, K_expansion=0.5):\n";
    std::cout << "  dP = " << sol_fittings.dP << " Pa  (extra dP from fittings: "
              << sol_fittings.dP - sol_rough.dP << " Pa)\n\n";

    // Correlation comparison sweep (velocity)
    std::cout << "Friction correlation comparison (D=" << D_seg * 1e3 << " mm, rough="
              << rough * 1e6 << " μm):\n";
    std::cout << std::setw(10) << "u [m/s]"
              << std::setw(10) << "Re"
              << std::setw(12) << "f_haaland"
              << std::setw(12) << "f_serghides"
              << std::setw(12) << "f_colebrook"
              << "\n";
    std::cout << std::string(56, '-') << "\n";
    for (double u : {5.0, 10.0, 20.0, 40.0, 80.0}) {
        auto sh = pipe_flow_rough(T_pipe, P_pipe, X_air, u, L_seg, D_seg, rough, "haaland");
        auto ss = pipe_flow_rough(T_pipe, P_pipe, X_air, u, L_seg, D_seg, rough, "serghides");
        auto sc = pipe_flow_rough(T_pipe, P_pipe, X_air, u, L_seg, D_seg, rough, "colebrook");
        std::cout << std::setw(10) << u
                  << std::setw(10) << static_cast<int>(sh.Re)
                  << std::setw(12) << sh.f
                  << std::setw(12) << ss.f
                  << std::setw(12) << sc.f
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 5. Friction factor correlations
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "5. Friction Factor Correlations (Moody chart)\n";
    std::cout << "=========================================================================\n\n";

    const double e_D = rough / D_seg;
    std::cout << "Relative roughness e/D = " << e_D << "\n\n";
    std::cout << std::setw(12) << "Re"
              << std::setw(12) << "Haaland"
              << std::setw(12) << "Serghides"
              << std::setw(12) << "Colebrook"
              << std::setw(12) << "Petukhov"
              << "\n";
    std::cout << std::string(60, '-') << "\n";
    for (double Re : {3000.0, 1e4, 5e4, 1e5, 5e5, 1e6}) {
        std::cout << std::setw(12) << static_cast<int>(Re)
                  << std::setw(12) << friction_haaland(Re, e_D)
                  << std::setw(12) << friction_serghides(Re, e_D)
                  << std::setw(12) << friction_colebrook(Re, e_D)
                  << std::setw(12) << friction_petukhov(Re)
                  << "\n";
    }
    std::cout << "\n";

    return 0;
}
