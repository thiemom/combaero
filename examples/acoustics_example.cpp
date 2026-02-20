#include "../include/acoustics.h"
#include "../include/geometry.h"
#include "../include/humidair.h"
#include "../include/thermo.h"

#include <iomanip>
#include <iostream>
#include <vector>

// Acoustics Example
//
// Demonstrates:
// 1. Tube axial modes (closed-closed, open-closed)
// 2. Annulus axial + azimuthal modes (combustor can geometry)
// 3. Mean-flow frequency splitting
// 4. SDOF liner impedance and absorption sweep
// 5. 2DOF serial liner vs SDOF comparison
// 6. Can-annular coupling eigenfrequencies

int main()
{
    std::cout << std::fixed << std::setprecision(2);

    // Representative combustor conditions
    const std::vector<double> X_air = standard_dry_air_composition();
    const double T_hot = 1800.0;   // K  hot-gas temperature
    const double P_comb = 15.0e5;  // Pa combustor pressure
    const double c = speed_of_sound(T_hot, X_air);
    const double rho = density(T_hot, P_comb, X_air);

    std::cout << "Combustor conditions: T=" << T_hot << " K, P=" << P_comb / 1e5
              << " bar\n";
    std::cout << "  c = " << c << " m/s,  rho = " << rho << " kg/m³\n\n";

    using namespace combaero;

    // =========================================================================
    // 1. Tube axial modes
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "1. Tube Axial Modes\n";
    std::cout << "=========================================================================\n\n";

    // Combustor liner tube: L=0.4 m, D=0.12 m
    Tube liner{0.4, 0.12};

    // Closed-closed (both ends rigid — pressure antinodes at both ends)
    auto modes_cc = tube_axial_modes(liner, c,
                                     BoundaryCondition::Closed,
                                     BoundaryCondition::Closed, 4);

    std::cout << "Tube L=" << liner.L << " m, D=" << liner.D << " m, c=" << c << " m/s\n\n";
    std::cout << "Closed-Closed (rigid walls):\n";
    std::cout << std::setw(8) << "Mode" << std::setw(12) << "f [Hz]" << "\n";
    std::cout << std::string(20, '-') << "\n";
    for (const auto& m : modes_cc) {
        std::cout << std::setw(8) << m.label() << std::setw(12) << m.frequency << "\n";
    }

    // Open-closed (inlet open, outlet closed — quarter-wave resonator)
    auto modes_oc = tube_axial_modes(liner, c,
                                     BoundaryCondition::Open,
                                     BoundaryCondition::Closed, 4);
    std::cout << "\nOpen-Closed (inlet open, outlet rigid):\n";
    std::cout << std::setw(8) << "Mode" << std::setw(12) << "f [Hz]" << "\n";
    std::cout << std::string(20, '-') << "\n";
    for (const auto& m : modes_oc) {
        std::cout << std::setw(8) << m.label() << std::setw(12) << m.frequency << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 2. Annulus axial + azimuthal modes
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "2. Annulus Modes (combustor annulus)\n";
    std::cout << "=========================================================================\n\n";

    // Annular combustor: L=0.35 m, D_inner=0.5 m, D_outer=0.7 m
    Annulus annulus{0.35, 0.5, 0.7};

    auto modes_ann = annulus_modes(annulus, c,
                                   BoundaryCondition::Closed,
                                   BoundaryCondition::Closed, 3, 3);

    std::cout << "Annulus: L=" << annulus.L << " m, D_in=" << annulus.D_inner
              << " m, D_out=" << annulus.D_outer << " m\n\n";
    std::cout << std::setw(10) << "Mode"
              << std::setw(12) << "f [Hz]"
              << std::setw(8) << "n_ax"
              << std::setw(8) << "m_az"
              << "\n";
    std::cout << std::string(38, '-') << "\n";
    for (const auto& m : modes_ann) {
        std::cout << std::setw(10) << m.label()
                  << std::setw(12) << m.frequency
                  << std::setw(8) << m.n_axial
                  << std::setw(8) << m.n_azimuthal
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 3. Mean-flow frequency splitting
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "3. Mean-Flow Frequency Splitting\n";
    std::cout << "=========================================================================\n\n";

    const double M_flow = 0.08;  // typical combustor Mach number
    std::cout << "Mean flow Mach number M = " << M_flow << "\n\n";
    std::cout << std::setw(10) << "Mode"
              << std::setw(14) << "f0 [Hz]"
              << std::setw(14) << "f_up [Hz]"
              << std::setw(14) << "f_down [Hz]"
              << "\n";
    std::cout << std::string(52, '-') << "\n";
    for (const auto& m : modes_cc) {
        auto [f_up, f_dn] = axial_mode_split(m.frequency, M_flow);
        std::cout << std::setw(10) << m.label()
                  << std::setw(14) << m.frequency
                  << std::setw(14) << f_up
                  << std::setw(14) << f_dn
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 4. SDOF liner impedance and absorption sweep
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "4. SDOF Liner Absorption Sweep\n";
    std::cout << "=========================================================================\n\n";

    // Typical perforated liner panel
    LinerOrificeGeometry sdof_orifice{0.002, 0.001, 0.05, 0.76};  // d, l, sigma, Cd
    LinerCavity sdof_cavity{0.025};                                 // depth 25 mm
    LinerFlowState sdof_flow{2.0, 15.0};                            // u_bias, u_grazing [m/s]
    AcousticMedium medium{rho, c};

    // Frequency sweep 100–1000 Hz
    std::vector<double> freqs;
    for (double f = 100.0; f <= 1000.0; f += 50.0) {
        freqs.push_back(f);
    }
    auto alpha_sdof = sweep_liner_sdof_absorption(freqs, sdof_orifice, sdof_cavity,
                                                   sdof_flow, medium);

    std::cout << "SDOF liner: d=" << sdof_orifice.d_orifice * 1e3 << " mm"
              << ", l=" << sdof_orifice.l_orifice * 1e3 << " mm"
              << ", sigma=" << sdof_orifice.porosity
              << ", depth=" << sdof_cavity.depth * 1e3 << " mm\n";
    std::cout << "Flow: u_bias=" << sdof_flow.u_bias << " m/s"
              << ", u_grazing=" << sdof_flow.u_grazing << " m/s\n\n";
    std::cout << std::setw(10) << "f [Hz]" << std::setw(12) << "alpha [-]" << "\n";
    std::cout << std::string(22, '-') << "\n";
    double alpha_max = 0.0;
    double f_peak = 0.0;
    for (std::size_t i = 0; i < freqs.size(); ++i) {
        std::cout << std::setw(10) << freqs[i] << std::setw(12) << alpha_sdof[i] << "\n";
        if (alpha_sdof[i] > alpha_max) {
            alpha_max = alpha_sdof[i];
            f_peak = freqs[i];
        }
    }
    std::cout << "\nPeak absorption: alpha=" << alpha_max << " at f=" << f_peak << " Hz\n\n";

    // =========================================================================
    // 5. 2DOF serial liner vs SDOF
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "5. 2DOF Serial Liner vs SDOF\n";
    std::cout << "=========================================================================\n\n";

    // 2DOF: face sheet + septum, two cavities
    LinerOrificeGeometry face{0.002, 0.001, 0.05, 0.76};
    LinerOrificeGeometry septum{0.003, 0.001, 0.10, 0.76};
    LinerFlowState face_flow{2.0, 15.0};
    LinerFlowState septum_flow{0.5, 5.0};

    auto alpha_2dof = sweep_liner_2dof_serial_absorption(freqs, face, septum,
                                                          0.015, 0.025,
                                                          face_flow, septum_flow, medium);

    std::cout << "2DOF: depth1=15 mm, depth2=25 mm\n\n";
    std::cout << std::setw(10) << "f [Hz]"
              << std::setw(12) << "SDOF"
              << std::setw(12) << "2DOF"
              << "\n";
    std::cout << std::string(34, '-') << "\n";
    for (std::size_t i = 0; i < freqs.size(); ++i) {
        std::cout << std::setw(10) << freqs[i]
                  << std::setw(12) << alpha_sdof[i]
                  << std::setw(12) << alpha_2dof[i]
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // 6. Acoustic properties at peak frequency
    // =========================================================================

    std::cout << "=========================================================================\n";
    std::cout << "6. Acoustic Properties at Peak Frequency\n";
    std::cout << "=========================================================================\n\n";

    const double p_rms = 500.0;  // Pa  representative combustion noise level
    AcousticProperties props = acoustic_properties(f_peak, rho, c, p_rms);

    std::cout << "At f=" << f_peak << " Hz, p_rms=" << p_rms << " Pa:\n";
    std::cout << "  wavelength        = " << props.wavelength << " m\n";
    std::cout << "  impedance (rho*c) = " << props.impedance << " Pa·s/m\n";
    std::cout << "  particle velocity = " << props.particle_velocity << " m/s\n";
    std::cout << "  SPL               = " << props.spl << " dB\n\n";

    // Mode closest to peak
    const AcousticMode* nearest = closest_mode(modes_cc, f_peak);
    if (nearest) {
        std::cout << "Nearest tube mode to peak: " << nearest->label()
                  << " at " << nearest->frequency << " Hz"
                  << "  (separation: " << std::abs(nearest->frequency - f_peak) << " Hz)\n";
    }

    return 0;
}
