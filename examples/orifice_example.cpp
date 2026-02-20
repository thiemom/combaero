#include "../include/humidair.h"
#include "../include/orifice.h"
#include "../include/thermo.h"
#include "../include/transport.h"

#include <iomanip>
#include <iostream>
#include <vector>

// Orifice Example
//
// Demonstrates:
// 1. Cd correlations: Reader-Harris/Gallagher, Stolz, Miller (sharp thin-plate)
// 2. Thick-plate and rounded-entry corrections
// 3. Iterative Cd-Re solver (solve_orifice_mdot)
// 4. Full orifice_flow bundle (OrificeFlowResult)
// 5. Expansibility factor for compressible flow
// 6. Inverse: Cd from measurement

int main()
{
    std::cout << std::fixed << std::setprecision(5);

    const std::vector<double> X_air = standard_dry_air_composition();

    // =========================================================================
    // 1. Cd correlations — sharp thin-plate sweep over beta and Re
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "1. Sharp Thin-Plate Cd Correlations (beta and Re sweep)\n";
    std::cout << "=========================================================================\n\n";

    // Fixed pipe D = 100 mm (ISO 5167 minimum)
    const double D_pipe = 0.100;

    std::cout << "Pipe D=" << D_pipe * 1e3 << " mm\n\n";
    std::cout << std::setw(8)  << "beta"
              << std::setw(10) << "Re_D"
              << std::setw(12) << "RHG"
              << std::setw(10) << "Stolz"
              << std::setw(10) << "Miller"
              << "\n";
    std::cout << std::string(50, '-') << "\n";

    for (double beta : {0.3, 0.5, 0.65, 0.75}) {
        for (double Re_D : {1.0e4, 1.0e5, 1.0e6}) {
            OrificeGeometry geom;
            geom.d = beta * D_pipe;
            geom.D = D_pipe;

            OrificeState state;
            state.Re_D = Re_D;
            state.dP   = 1000.0;   // Pa (not used by Cd correlations directly)
            state.rho  = 1.2;
            state.mu   = 1.8e-5;

            const double cd_rhg   = orifice::Cd_ReaderHarrisGallagher(beta, Re_D, D_pipe);
            const double cd_stolz = orifice::Cd_Stolz(beta, Re_D);
            const double cd_miller = orifice::Cd_Miller(beta, Re_D);

            std::cout << std::setw(8)  << beta
                      << std::setw(10) << static_cast<int>(Re_D)
                      << std::setw(12) << cd_rhg
                      << std::setw(10) << cd_stolz
                      << std::setw(10) << cd_miller
                      << "\n";
        }
    }
    std::cout << "\n";

    // =========================================================================
    // 2. Thick-plate and rounded-entry corrections
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "2. Thick-Plate and Rounded-Entry Corrections\n";
    std::cout << "=========================================================================\n\n";

    const double beta_ref = 0.5;
    const double Re_ref   = 1.0e5;

    OrificeGeometry geom_base;
    geom_base.d = beta_ref * D_pipe;
    geom_base.D = D_pipe;

    OrificeState state_ref;
    state_ref.Re_D = Re_ref;
    state_ref.dP   = 5000.0;
    state_ref.rho  = 1.2;
    state_ref.mu   = 1.8e-5;

    const double cd_thin = orifice::Cd_ReaderHarrisGallagher(beta_ref, Re_ref, D_pipe);
    std::cout << "Reference (thin plate, beta=" << beta_ref << ", Re=" << Re_ref << "):\n";
    std::cout << "  Cd_thin = " << cd_thin << "\n\n";

    // Thickness correction: t/d from 0 to 2
    std::cout << "Thick-plate correction (t/d sweep):\n";
    std::cout << std::setw(10) << "t/d" << std::setw(12) << "correction" << std::setw(12) << "Cd_eff" << "\n";
    std::cout << std::string(34, '-') << "\n";
    for (double t_over_d : {0.0, 0.25, 0.5, 1.0, 1.5, 2.0}) {
        const double corr = orifice::thickness_correction(t_over_d, beta_ref,
                                                           state_ref.Re_d(beta_ref));
        std::cout << std::setw(10) << t_over_d
                  << std::setw(12) << corr
                  << std::setw(12) << cd_thin * corr
                  << "\n";
    }

    // Rounded-entry: r/d from 0 to 0.2
    std::cout << "\nRounded-entry correction (r/d sweep):\n";
    std::cout << std::setw(10) << "r/d" << std::setw(12) << "Cd" << "\n";
    std::cout << std::string(22, '-') << "\n";
    for (double r_over_d : {0.0, 0.02, 0.05, 0.10, 0.15, 0.20}) {
        const double cd_r = orifice::Cd_rounded(r_over_d, beta_ref, Re_ref);
        std::cout << std::setw(10) << r_over_d << std::setw(12) << cd_r << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 3. Iterative Cd-Re solver
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "3. Iterative Cd-Re Solver (solve_orifice_mdot)\n";
    std::cout << "=========================================================================\n\n";

    // Hot gas conditions: 600 K, 3 bar
    const double T_hot  = 600.0;
    const double P_hot  = 3.0e5;
    const double rho_hot = density(T_hot, P_hot, X_air);
    const double mu_hot  = viscosity(T_hot, P_hot, X_air);

    OrificeGeometry geom_iter;
    geom_iter.d = 0.020;   // 20 mm orifice
    geom_iter.D = 0.050;   // 50 mm pipe

    std::cout << "Hot gas: T=" << T_hot << " K, P=" << P_hot / 1e5 << " bar\n";
    std::cout << "  rho=" << rho_hot << " kg/m³,  mu=" << mu_hot * 1e6 << " μPa·s\n";
    std::cout << "Orifice: d=" << geom_iter.d * 1e3 << " mm, D=" << geom_iter.D * 1e3 << " mm"
              << ", beta=" << geom_iter.beta() << "\n\n";

    std::cout << "dP sweep (iterative Cd-Re vs fixed Cd=0.61):\n";
    std::cout << std::setw(10) << "dP [Pa]"
              << std::setw(12) << "mdot_iter"
              << std::setw(10) << "Cd_conv"
              << std::setw(12) << "mdot_fixed"
              << std::setw(10) << "diff [%]"
              << "\n";
    std::cout << std::string(54, '-') << "\n";

    for (double dP : {500.0, 2000.0, 5000.0, 10000.0, 20000.0}) {
        const double mdot_iter = solve_orifice_mdot(geom_iter, dP, rho_hot, mu_hot,
                                                     P_hot, 0.0,
                                                     CdCorrelation::ReaderHarrisGallagher);
        // Fixed Cd for comparison
        const double mdot_fixed = orifice_mdot(geom_iter, 0.61, dP, rho_hot);
        // Back-compute converged Cd
        const double cd_conv = orifice_Cd_from_measurement(geom_iter, mdot_iter, dP, rho_hot);
        const double diff_pct = (mdot_iter - mdot_fixed) / mdot_fixed * 100.0;

        std::cout << std::setw(10) << dP
                  << std::setw(12) << mdot_iter
                  << std::setw(10) << cd_conv
                  << std::setw(12) << mdot_fixed
                  << std::setw(10) << diff_pct
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 4. Full orifice_flow bundle
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "4. Full orifice_flow Bundle (OrificeFlowResult)\n";
    std::cout << "=========================================================================\n\n";

    const double dP_ref = 5000.0;
    OrificeFlowResult res = orifice_flow(geom_iter, dP_ref, T_hot, P_hot, mu_hot,
                                          1.0, X_air, 0.0,
                                          CdCorrelation::ReaderHarrisGallagher);

    std::cout << "orifice_flow (dP=" << dP_ref << " Pa):\n";
    std::cout << "  mdot         = " << res.mdot << " kg/s\n";
    std::cout << "  v            = " << res.v << " m/s\n";
    std::cout << "  Re_D         = " << res.Re_D << "\n";
    std::cout << "  Re_d         = " << res.Re_d << "\n";
    std::cout << "  Cd           = " << res.Cd << "\n";
    std::cout << "  epsilon      = " << res.epsilon << "  (1.0 = incompressible)\n";
    std::cout << "  rho_corrected= " << res.rho_corrected << " kg/m³\n\n";

    // =========================================================================
    // 5. Expansibility factor for compressible flow
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "5. Expansibility Factor (compressible correction)\n";
    std::cout << "=========================================================================\n\n";

    const double kappa = 1.38;  // cp/cv for hot air at 600 K
    std::cout << "beta=" << geom_iter.beta() << ", P_up=" << P_hot / 1e5
              << " bar, kappa=" << kappa << "\n\n";
    std::cout << std::setw(12) << "dP/P [%]"
              << std::setw(12) << "epsilon"
              << std::setw(14) << "mdot_incomp"
              << std::setw(14) << "mdot_comp"
              << "\n";
    std::cout << std::string(52, '-') << "\n";

    for (double dP_frac : {0.01, 0.05, 0.10, 0.15, 0.20, 0.25}) {
        const double dP_c = dP_frac * P_hot;
        const double eps  = expansibility_factor(geom_iter.beta(), dP_c, P_hot, kappa);
        const double mdot_inc = orifice_mdot(geom_iter, 0.62, dP_c, rho_hot);
        const double mdot_cmp = eps * mdot_inc;
        std::cout << std::setw(12) << dP_frac * 100.0
                  << std::setw(12) << eps
                  << std::setw(14) << mdot_inc
                  << std::setw(14) << mdot_cmp
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 6. Inverse: Cd from measurement
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "6. Inverse: Cd from Measurement\n";
    std::cout << "=========================================================================\n\n";

    // Simulate a calibration measurement: known mdot, measure dP
    const double mdot_measured = 0.050;   // kg/s
    const double dP_measured   = 8200.0;  // Pa

    const double Cd_measured = orifice_Cd_from_measurement(geom_iter, mdot_measured,
                                                             dP_measured, rho_hot);
    std::cout << "Calibration measurement:\n";
    std::cout << "  mdot = " << mdot_measured << " kg/s\n";
    std::cout << "  dP   = " << dP_measured << " Pa\n";
    std::cout << "  Cd   = " << Cd_measured << "\n\n";

    // Verify round-trip
    const double mdot_check = orifice_mdot(geom_iter, Cd_measured, dP_measured, rho_hot);
    std::cout << "Round-trip check: mdot_check = " << mdot_check
              << " (error: " << std::abs(mdot_check - mdot_measured) / mdot_measured * 100.0
              << " %)\n\n";

    // K-loss / zeta conversion
    const double K = orifice::K_from_Cd(Cd_measured, geom_iter.beta());
    const double Cd_back = orifice::Cd_from_K(K, geom_iter.beta());
    std::cout << "K-loss conversion:\n";
    std::cout << "  Cd=" << Cd_measured << "  ->  K=" << K
              << "  ->  Cd_back=" << Cd_back << "\n\n";

    return 0;
}
