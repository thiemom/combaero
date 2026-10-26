#include "../include/acoustics.h"
#include "../include/acoustics_internal.h"
#include "../include/can_annular_solvers.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include "../include/math_constants.h"

namespace combaero {

// -------------------------------------------------------------
// AcousticMode methods
// -------------------------------------------------------------

std::string AcousticMode::label() const {
    std::string result;

    if (n_axial > 0) {
        result += std::to_string(n_axial) + "L";
    }
    if (n_azimuthal > 0) {
        result += std::to_string(n_azimuthal) + "T";
    }

    // Handle bulk mode (0L0T) - shouldn't normally occur
    if (result.empty()) {
        result = "0";
    }

    return result;
}

// -------------------------------------------------------------
// Tube Acoustic Modes
// -------------------------------------------------------------

std::vector<AcousticMode> tube_axial_modes(
    const Tube& tube, double c,
    BoundaryCondition upstream, BoundaryCondition downstream,
    int n_max) {

    if (tube.L <= 0 || tube.D <= 0) {
        throw std::invalid_argument("tube_axial_modes: tube dimensions must be positive");
    }
    if (c <= 0) {
        throw std::invalid_argument("tube_axial_modes: speed of sound must be positive");
    }
    if (n_max < 1) {
        throw std::invalid_argument("tube_axial_modes: n_max must be >= 1");
    }

    std::vector<AcousticMode> modes;
    modes.reserve(static_cast<std::size_t>(n_max));

    bool same_bc = (upstream == downstream);

    for (int n = 1; n <= n_max; ++n) {
        AcousticMode mode;
        mode.n_axial = n;
        mode.n_azimuthal = 0;

        if (same_bc) {
            // Closed-Closed or Open-Open: f_n = n * c / (2L)
            mode.frequency = n * c / (2.0 * tube.L);
        } else {
            // Open-Closed (quarter-wave): f_n = (2n-1) * c / (4L)
            mode.frequency = (2 * n - 1) * c / (4.0 * tube.L);
        }

        modes.push_back(mode);
    }

    return modes;
}

// -------------------------------------------------------------
// Annulus Acoustic Modes
// -------------------------------------------------------------

std::vector<AcousticMode> annulus_axial_modes(
    const Annulus& annulus, double c,
    BoundaryCondition upstream, BoundaryCondition downstream,
    int n_max) {

    if (annulus.L <= 0 || annulus.D_outer <= annulus.D_inner || annulus.D_inner < 0) {
        throw std::invalid_argument("annulus_axial_modes: invalid annulus dimensions");
    }
    if (c <= 0) {
        throw std::invalid_argument("annulus_axial_modes: speed of sound must be positive");
    }
    if (n_max < 1) {
        throw std::invalid_argument("annulus_axial_modes: n_max must be >= 1");
    }

    std::vector<AcousticMode> modes;
    modes.reserve(static_cast<std::size_t>(n_max));

    bool same_bc = (upstream == downstream);

    for (int n = 1; n <= n_max; ++n) {
        AcousticMode mode;
        mode.n_axial = n;
        mode.n_azimuthal = 0;

        if (same_bc) {
            mode.frequency = n * c / (2.0 * annulus.L);
        } else {
            mode.frequency = (2 * n - 1) * c / (4.0 * annulus.L);
        }

        modes.push_back(mode);
    }

    return modes;
}

std::vector<AcousticMode> annulus_azimuthal_modes(
    const Annulus& annulus, double c,
    int m_max) {

    if (annulus.D_outer <= annulus.D_inner || annulus.D_inner < 0) {
        throw std::invalid_argument("annulus_azimuthal_modes: invalid annulus dimensions");
    }
    if (c <= 0) {
        throw std::invalid_argument("annulus_azimuthal_modes: speed of sound must be positive");
    }
    if (m_max < 1) {
        throw std::invalid_argument("annulus_azimuthal_modes: m_max must be >= 1");
    }

    std::vector<AcousticMode> modes;
    modes.reserve(static_cast<std::size_t>(m_max));

    double D_mean = annulus.D_mean();

    for (int m = 1; m <= m_max; ++m) {
        AcousticMode mode;
        mode.n_axial = 0;
        mode.n_azimuthal = m;

        // f_m = m * c / (π * D_mean)
        mode.frequency = m * c / (M_PI * D_mean);

        modes.push_back(mode);
    }

    return modes;
}

std::vector<AcousticMode> annulus_modes(
    const Annulus& annulus, double c,
    BoundaryCondition upstream, BoundaryCondition downstream,
    int n_max, int m_max) {

    // Get pure axial and azimuthal modes
    auto axial = annulus_axial_modes(annulus, c, upstream, downstream, n_max);
    auto azimuthal = annulus_azimuthal_modes(annulus, c, m_max);

    std::vector<AcousticMode> modes;

    // Add pure axial modes (mT = 0)
    for (const auto& mode : axial) {
        modes.push_back(mode);
    }

    // Add pure azimuthal modes (nL = 0)
    for (const auto& mode : azimuthal) {
        modes.push_back(mode);
    }

    // Add combined modes
    for (const auto& ax : axial) {
        for (const auto& az : azimuthal) {
            AcousticMode combined;
            combined.n_axial = ax.n_axial;
            combined.n_azimuthal = az.n_azimuthal;
            // Combined frequency: f = sqrt(f_axial² + f_azimuthal²)
            combined.frequency = std::sqrt(ax.frequency * ax.frequency +
                                           az.frequency * az.frequency);
            modes.push_back(combined);
        }
    }

    // Sort by frequency
    std::sort(modes.begin(), modes.end(),
              [](const AcousticMode& a, const AcousticMode& b) {
                  return a.frequency < b.frequency;
              });

    return modes;
}

// -------------------------------------------------------------
// Utility Functions
// -------------------------------------------------------------

std::vector<AcousticMode> modes_in_range(
    const std::vector<AcousticMode>& modes,
    double f_min, double f_max) {

    std::vector<AcousticMode> result;

    for (const auto& mode : modes) {
        if (mode.frequency >= f_min && mode.frequency <= f_max) {
            result.push_back(mode);
        }
    }

    return result;
}

const AcousticMode* closest_mode(
    const std::vector<AcousticMode>& modes,
    double f_target) {

    if (modes.empty()) {
        return nullptr;
    }

    const AcousticMode* closest = &modes[0];
    double min_diff = std::abs(modes[0].frequency - f_target);

    for (const auto& mode : modes) {
        double diff = std::abs(mode.frequency - f_target);
        if (diff < min_diff) {
            min_diff = diff;
            closest = &mode;
        }
    }

    return closest;
}

double min_mode_separation(const std::vector<AcousticMode>& modes) {
    if (modes.size() < 2) {
        return 0.0;
    }

    // Make a sorted copy
    std::vector<double> freqs;
    freqs.reserve(modes.size());
    for (const auto& mode : modes) {
        freqs.push_back(mode.frequency);
    }
    std::sort(freqs.begin(), freqs.end());

    double min_sep = freqs[1] - freqs[0];
    for (std::size_t i = 2; i < freqs.size(); ++i) {
        double sep = freqs[i] - freqs[i - 1];
        if (sep < min_sep) {
            min_sep = sep;
        }
    }

    return min_sep;
}

// -------------------------------------------------------------
// Mean Flow Correction
// -------------------------------------------------------------

double axial_mode_upstream(double f0, double M) {
    if (f0 <= 0) {
        throw std::invalid_argument("axial_mode_upstream: f0 must be positive");
    }
    if (M < 0 || M >= 1) {
        throw std::invalid_argument("axial_mode_upstream: M must be in [0, 1)");
    }
    return f0 / (1.0 - M);
}

double axial_mode_downstream(double f0, double M) {
    if (f0 <= 0) {
        throw std::invalid_argument("axial_mode_downstream: f0 must be positive");
    }
    if (M < 0) {
        throw std::invalid_argument("axial_mode_downstream: M must be >= 0");
    }
    return f0 / (1.0 + M);
}

std::pair<double, double> axial_mode_split(double f0, double M) {
    return {axial_mode_upstream(f0, M), axial_mode_downstream(f0, M)};
}

// -------------------------------------------------------------
// Helmholtz Resonator
// -------------------------------------------------------------

double helmholtz_frequency(double V, double A_neck, double L_neck, double c,
                           double end_correction) {
    if (V <= 0 || A_neck <= 0 || L_neck <= 0 || c <= 0) {
        throw std::invalid_argument("helmholtz_frequency: all parameters must be positive");
    }
    if (end_correction < 0) {
        throw std::invalid_argument("helmholtz_frequency: end_correction must be >= 0");
    }

    // Effective neck length with end correction
    // d_neck = 2 * sqrt(A_neck / π)
    double d_neck = 2.0 * std::sqrt(A_neck / M_PI);
    double L_eff = L_neck + end_correction * d_neck;

    // f = (c / 2π) * sqrt(A / (V * L_eff))
    return (c / (2.0 * M_PI)) * std::sqrt(A_neck / (V * L_eff));
}

// -------------------------------------------------------------
// Strouhal Number
// -------------------------------------------------------------

double strouhal(double f, double L, double u) {
    if (f <= 0 || L <= 0 || u <= 0) {
        throw std::invalid_argument("strouhal: f, L, and u must be positive");
    }
    return f * L / u;
}

double frequency_from_strouhal(double St, double L, double u) {
    if (St <= 0 || L <= 0 || u <= 0) {
        throw std::invalid_argument("frequency_from_strouhal: St, L, and u must be positive");
    }
    return St * u / L;
}

// -------------------------------------------------------------
// Convenience Functions
// -------------------------------------------------------------

double quarter_wave_frequency(double L, double c) {
    if (L <= 0 || c <= 0) {
        throw std::invalid_argument("quarter_wave_frequency: L and c must be positive");
    }
    return c / (4.0 * L);
}

double half_wave_frequency(double L, double c) {
    if (L <= 0 || c <= 0) {
        throw std::invalid_argument("half_wave_frequency: L and c must be positive");
    }
    return c / (2.0 * L);
}

// -------------------------------------------------------------
// Viscothermal Boundary Layers
// -------------------------------------------------------------

double stokes_layer(double nu, double f) {
    if (nu <= 0 || f <= 0) {
        throw std::invalid_argument("stokes_layer: nu and f must be positive");
    }
    double omega = 2.0 * M_PI * f;
    return std::sqrt(2.0 * nu / omega);
}

double thermal_layer(double alpha, double f) {
    if (alpha <= 0 || f <= 0) {
        throw std::invalid_argument("thermal_layer: alpha and f must be positive");
    }
    double omega = 2.0 * M_PI * f;
    return std::sqrt(2.0 * alpha / omega);
}

double effective_viscothermal_layer(double delta_nu, double delta_kappa, double gamma) {
    if (delta_nu <= 0 || delta_kappa <= 0) {
        throw std::invalid_argument("effective_viscothermal_layer: deltas must be positive");
    }
    if (gamma <= 1.0) {
        throw std::invalid_argument("effective_viscothermal_layer: gamma must be > 1");
    }
    return delta_nu + (gamma - 1.0) * delta_kappa;
}

// -------------------------------------------------------------
// Quality Factor (Screening Estimates)
// -------------------------------------------------------------

double helmholtz_Q(double V, double A_neck, double L_neck,
                   double nu, double alpha, double gamma, double f) {
    if (V <= 0 || A_neck <= 0 || L_neck <= 0) {
        throw std::invalid_argument("helmholtz_Q: geometry parameters must be positive");
    }
    if (nu <= 0 || alpha <= 0 || f <= 0) {
        throw std::invalid_argument("helmholtz_Q: fluid/frequency parameters must be positive");
    }
    if (gamma <= 1.0) {
        throw std::invalid_argument("helmholtz_Q: gamma must be > 1");
    }

    // Boundary layer thicknesses
    double delta_nu = stokes_layer(nu, f);
    double delta_kappa = thermal_layer(alpha, f);
    double delta_eff = effective_viscothermal_layer(delta_nu, delta_kappa, gamma);

    // Neck diameter
    double d_neck = 2.0 * std::sqrt(A_neck / M_PI);

    // Effective neck length (with standard end correction)
    double L_eff = L_neck + kEndCorrectionFlangedBothFaces * d_neck;

    // Lumped element: Q = omega * M / R, with inertance M = rho * L_eff / A
    // and the Kirchhoff wall resistance of the neck, rho * omega * delta_eff
    // / (a * A) per unit length, over the physical neck length L_neck:
    //   Q = a * L_eff / (delta_eff * L_neck),  a = d_neck / 2.
    // V only sets the frequency (pass helmholtz_frequency as f); at a fixed f
    // the cavity volume does not enter. Cavity-wall thermal loss is neglected.
    double Q = 0.5 * d_neck * L_eff / (delta_eff * L_neck);

    return Q;
}

double tube_Q(double L, double D, double nu, double alpha, double gamma, double f) {
    if (L <= 0 || D <= 0) {
        throw std::invalid_argument("tube_Q: L and D must be positive");
    }
    if (nu <= 0 || alpha <= 0 || f <= 0) {
        throw std::invalid_argument("tube_Q: fluid/frequency parameters must be positive");
    }
    if (gamma <= 1.0) {
        throw std::invalid_argument("tube_Q: gamma must be > 1");
    }

    // Boundary layer thicknesses
    double delta_nu = stokes_layer(nu, f);
    double delta_kappa = thermal_layer(alpha, f);
    double delta_eff = effective_viscothermal_layer(delta_nu, delta_kappa, gamma);

    // For a tube, losses occur along the entire length
    // Kirchhoff wide-tube attenuation: alpha = k * delta_eff / D.
    // A line resonator has Q = k / (2 alpha), so Q = D / (2 * delta_eff).
    double Q = D / (2.0 * delta_eff);

    return Q;
}

double damping_ratio(double Q) {
    if (Q <= 0) {
        throw std::invalid_argument("damping_ratio: Q must be positive");
    }
    return 1.0 / (2.0 * Q);
}

double bandwidth(double f0, double Q) {
    if (f0 <= 0 || Q <= 0) {
        throw std::invalid_argument("bandwidth: f0 and Q must be positive");
    }
    return f0 / Q;
}

// -------------------------------------------------------------
// Acoustic properties bundle
// -------------------------------------------------------------

AcousticProperties acoustic_properties(
    double f,
    [[maybe_unused]] double rho,
    double c,
    double p_rms,
    double p_ref)
{
    // Input validation
    if (f <= 0.0) {
        throw std::invalid_argument("acoustic_properties: frequency must be positive");
    }
    if (rho <= 0.0) {
        throw std::invalid_argument("acoustic_properties: density must be positive");
    }
    if (c <= 0.0) {
        throw std::invalid_argument("acoustic_properties: speed of sound must be positive");
    }
    if (p_rms < 0.0) {
        throw std::invalid_argument("acoustic_properties: p_rms must be non-negative");
    }
    if (p_ref <= 0.0) {
        throw std::invalid_argument("acoustic_properties: p_ref must be positive");
    }

    AcousticProperties props;
    props.frequency = f;
    props.wavelength = c / f;
    props.impedance = rho * c;
    props.particle_velocity = p_rms / (rho * c);
    props.spl = 20.0 * std::log10(p_rms / p_ref);

    return props;
}

// -------------------------------------------------------------
// Transfer Matrix Method
// -------------------------------------------------------------

// Transfer matrix multiplication for cascading elements
TransferMatrix TransferMatrix::operator*(const TransferMatrix& other) const {
    return {
        T11 * other.T11 + T12 * other.T21,
        T11 * other.T12 + T12 * other.T22,
        T21 * other.T11 + T22 * other.T21,
        T21 * other.T12 + T22 * other.T22
    };
}

double absorption_from_impedance_norm(const std::complex<double>& z_norm) {
    const std::complex<double> one(1.0, 0.0);
    const std::complex<double> denom = z_norm + one;
    if (std::abs(denom) < 1e-12) {
        return 0.0;
    }

    const std::complex<double> r = (z_norm - one) / denom;
    const double alpha = 1.0 - std::norm(r);
    // Use manual clamp instead of std::clamp for maximum compatibility
    return std::max(0.0, std::min(1.0, alpha));
}

std::complex<double> liner_sdof_impedance_norm(
    double freq,
    const LinerOrificeGeometry& orifice,
    const LinerCavity& cavity,
    const LinerFlowState& flow,
    const AcousticMedium& medium
) {
    if (cavity.depth <= 0.0) {
        throw std::runtime_error("Cavity depth must be positive");
    }

    const std::complex<double> z_neck = orifice_impedance_with_flow(
        freq,
        flow.u_bias,
        flow.u_grazing,
        orifice.d_orifice,
        orifice.l_orifice,
        orifice.porosity,
        orifice.Cd,
        medium.rho,
        medium.c
    );

    const double k = (2.0 * M_PI * freq) / medium.c;
    const std::complex<double> i(0.0, 1.0);
    const std::complex<double> z_cavity = -i / std::tan(k * cavity.depth);
    return z_neck + z_cavity;
}

double liner_sdof_absorption(
    double freq,
    const LinerOrificeGeometry& orifice,
    const LinerCavity& cavity,
    const LinerFlowState& flow,
    const AcousticMedium& medium
) {
    const std::complex<double> z_norm =
        liner_sdof_impedance_norm(freq, orifice, cavity, flow, medium);
    return absorption_from_impedance_norm(z_norm);
}

std::vector<double> sweep_liner_sdof_absorption(
    const std::vector<double>& freqs,
    const LinerOrificeGeometry& orifice,
    const LinerCavity& cavity,
    const LinerFlowState& flow,
    const AcousticMedium& medium
) {
    std::vector<double> alpha;
    alpha.reserve(freqs.size());
    for (double freq : freqs) {
        alpha.push_back(liner_sdof_absorption(freq, orifice, cavity, flow, medium));
    }
    return alpha;
}

std::complex<double> liner_2dof_serial_impedance_norm(
    double freq,
    const LinerOrificeGeometry& face_orifice,
    const LinerOrificeGeometry& septum_orifice,
    double depth_1,
    double depth_2,
    const LinerFlowState& face_flow,
    const LinerFlowState& septum_flow,
    const AcousticMedium& medium
) {
    if (depth_1 <= 0.0 || depth_2 <= 0.0) {
        throw std::runtime_error("Serial cavity depths must be positive");
    }

    const std::complex<double> z_face = orifice_impedance_with_flow(
        freq,
        face_flow.u_bias,
        face_flow.u_grazing,
        face_orifice.d_orifice,
        face_orifice.l_orifice,
        face_orifice.porosity,
        face_orifice.Cd,
        medium.rho,
        medium.c
    );

    const std::complex<double> z_septum = orifice_impedance_with_flow(
        freq,
        septum_flow.u_bias,
        septum_flow.u_grazing,
        septum_orifice.d_orifice,
        septum_orifice.l_orifice,
        septum_orifice.porosity,
        septum_orifice.Cd,
        medium.rho,
        medium.c
    );

    const double k = (2.0 * M_PI * freq) / medium.c;
    const std::complex<double> i(0.0, 1.0);

    const std::complex<double> z_vol2 = -i / std::tan(k * depth_2);
    const std::complex<double> z_load = z_septum + z_vol2;

    // Transmission line transform through cavity 1 with normalized Z_char = 1
    const std::complex<double> tan_kL1 = std::tan(k * depth_1);
    std::complex<double> denom = 1.0 + i * z_load * tan_kL1;
    if (std::abs(denom) < 1e-12) {
        denom = std::complex<double>(1e12, 0.0);
    }
    const std::complex<double> z_in_vol1 = (z_load + i * tan_kL1) / denom;

    return z_face + z_in_vol1;
}

double liner_2dof_serial_absorption(
    double freq,
    const LinerOrificeGeometry& face_orifice,
    const LinerOrificeGeometry& septum_orifice,
    double depth_1,
    double depth_2,
    const LinerFlowState& face_flow,
    const LinerFlowState& septum_flow,
    const AcousticMedium& medium
) {
    const std::complex<double> z_norm = liner_2dof_serial_impedance_norm(
        freq,
        face_orifice,
        septum_orifice,
        depth_1,
        depth_2,
        face_flow,
        septum_flow,
        medium
    );
    return absorption_from_impedance_norm(z_norm);
}

std::vector<double> sweep_liner_2dof_serial_absorption(
    const std::vector<double>& freqs,
    const LinerOrificeGeometry& face_orifice,
    const LinerOrificeGeometry& septum_orifice,
    double depth_1,
    double depth_2,
    const LinerFlowState& face_flow,
    const LinerFlowState& septum_flow,
    const AcousticMedium& medium
) {
    std::vector<double> alpha;
    alpha.reserve(freqs.size());
    for (double freq : freqs) {
        alpha.push_back(liner_2dof_serial_absorption(
            freq,
            face_orifice,
            septum_orifice,
            depth_1,
            depth_2,
            face_flow,
            septum_flow,
            medium
        ));
    }
    return alpha;
}

// Orifice impedance with bias and grazing flow
std::complex<double> orifice_impedance_with_flow(
    double freq,
    double u_bias,
    double u_grazing,
    double d_orifice,
    double l_orifice,
    double porosity,
    double Cd,
    [[maybe_unused]] double rho,
    double c
) {
    // Validate parameters
    if (freq <= 0.0) {
        throw std::runtime_error("Frequency must be positive, got " + std::to_string(freq));
    }
    if (d_orifice <= 0.0 || l_orifice <= 0.0) {
        throw std::runtime_error("Orifice dimensions must be positive");
    }
    if (porosity <= 0.0 || porosity > 0.5) {
        throw std::runtime_error("Porosity must be in range (0, 0.5], got " + std::to_string(porosity));
    }
    if (Cd <= 0.0 || Cd > 1.0) {
        throw std::runtime_error("Discharge coefficient must be in range (0, 1], got " + std::to_string(Cd));
    }
    if (!(c > 0.0) || !std::isfinite(c)) {
        throw std::runtime_error("Speed of sound must be positive and finite, got " + std::to_string(c));
    }
    if (!(rho > 0.0) || !std::isfinite(rho)) {
        throw std::runtime_error("Density must be positive and finite, got " + std::to_string(rho));
    }

    // Check Mach numbers (should be << 1 for acoustic theory)
    double M_bias = std::abs(u_bias) / c;
    double M_grazing = std::abs(u_grazing) / c;
    if (M_bias > 0.3) {
        throw std::runtime_error("Bias Mach number too high: " + std::to_string(M_bias) + " > 0.3");
    }
    if (M_grazing > 0.3) {
        throw std::runtime_error("Grazing Mach number too high: " + std::to_string(M_grazing) + " > 0.3");
    }

    double omega = 2.0 * M_PI * freq;

    // 1. Resistance (Real Part)
    // Bias flow resistance (Rogers & Marble 1956)
    // R_bias = M_b / (σ * Cd²)
    double R_bias = M_bias / (porosity * Cd * Cd);

    // Grazing flow resistance (Howe 1979, Bourquard & Noiray 2019)
    // R_grazing = M_g / (2σ)
    double R_grazing = M_grazing / (2.0 * porosity);

    // Total resistance (linear superposition for small Mach numbers)
    double R = R_bias + R_grazing;

    // 2. Reactance (Imaginary Part)
    // Rayleigh end correction modified by flow: (8/(3 pi)) r = 0.85 r on EACH
    // face of the orifice, 0.85 d in total -- the same as helmholtz_frequency.
    double delta_0 = kEndCorrectionFlangedBothFaces * d_orifice;

    // Flow reduces effective end correction (frequency shift)
    // Empirical: δ_eff = δ_0 / (1 + M_g)
    double flow_modifier = 1.0 / (1.0 + M_grazing);
    double l_eff = l_orifice + delta_0 * flow_modifier;

    // Inertance (mass reactance): X = ω*l_eff / (c*σ)
    double X = (omega * l_eff) / (c * porosity);

    // Return normalized impedance Z/(ρc)
    return std::complex<double>(R, X);
}

// Quarter-wave resonator transfer matrix
TransferMatrix quarter_wave_resonator_tmm(
    double freq,
    double L_tube,
    double A_duct,
    double A_tube,
    double d_orifice,
    double l_orifice,
    double porosity,
    double Cd,
    double u_bias,
    double u_grazing,
    double rho,
    double c
) {
    // Validate parameters
    if (L_tube <= 0.0) {
        throw std::runtime_error("Tube length must be positive, got " + std::to_string(L_tube));
    }
    if (A_duct <= 0.0 || A_tube <= 0.0) {
        throw std::runtime_error("Areas must be positive");
    }

    double k = (2.0 * M_PI * freq) / c;

    // Check that we're below cutoff (kL < π for quarter-wave)
    if (k * L_tube > M_PI) {
        throw std::runtime_error("Frequency too high: kL = " + std::to_string(k * L_tube) +
                                 " > π (above quarter-wave cutoff)");
    }

    // 1. Neck impedance (with flow effects)
    std::complex<double> Z_neck = orifice_impedance_with_flow(
        freq, u_bias, u_grazing, d_orifice, l_orifice, porosity, Cd, rho, c
    );

    // 2. Tube backing impedance (closed-end quarter-wave)
    // For a closed-end tube: Z_tube = -i*cot(kL)
    // This is the normalized impedance (Z/(ρc))
    std::complex<double> i(0.0, 1.0);
    std::complex<double> Z_tube = -i / std::tan(k * L_tube);

    // 3. Series combination. The two parts are normalized on different
    // areas, so they are denormalized by their own reference area first.
    // Z_neck is normalized on the panel velocity (orifice area / porosity =
    // A_duct, see orifice_impedance_with_flow); Z_tube on the tube's own
    // characteristic impedance rho*c/A_tube. Acoustic impedance p/(u*S):
    //   Z_branch = rho*c * (Z_neck / A_duct + Z_tube / A_tube)
    std::complex<double> Z_branch_acoustic =
        rho * c * (Z_neck / A_duct + Z_tube / A_tube);
    std::complex<double> Y_branch = 1.0 / Z_branch_acoustic;

    // 5. Transfer matrix for side-branch
    // For a shunt element (parallel admittance):
    // [p_u]   [1    0  ] [p_d]
    // [u_u*S] = [Y_br  1  ] [u_d*S]
    //
    // where Y_br is the branch admittance
    return TransferMatrix{
        std::complex<double>(1.0, 0.0),  // T11
        std::complex<double>(0.0, 0.0),  // T12
        Y_branch,                         // T21
        std::complex<double>(1.0, 0.0)   // T22
    };
}

// Whistling risk assessment based on Strouhal number
bool is_whistling_risk(
    double freq,
    double u_bias,
    double d_orifice
) {
    // Validate inputs
    if (freq < 0.0 || d_orifice <= 0.0) {
        throw std::runtime_error("Invalid parameters for whistling risk assessment");
    }

    // No flow, no whistling
    if (std::abs(u_bias) < 0.01) {
        return false;
    }

    // Strouhal number: St = f * d / U
    double St = (freq * d_orifice) / std::abs(u_bias);

    // Critical Strouhal range for vortex lock-in: 0.2 < St < 0.5
    // This is when vortex shedding frequency matches acoustic frequency
    return (St >= 0.2 && St <= 0.5);
}

// -------------------------------------------------------------
// Can-Annular Combustor Acoustics (Bloch-Floquet Theory)
// -------------------------------------------------------------

// Adapter: Convert CanAnnularFlowGeometry to CanAnnularGeometry
CanAnnularGeometry to_acoustic_geometry(
    const CanAnnularFlowGeometry& flow_geom,
    int n_cans,
    double L_can,
    double D_can
) {
    // Validate inputs
    if (n_cans < 1) {
        throw std::invalid_argument("to_acoustic_geometry: n_cans must be positive");
    }
    if (L_can <= 0.0) {
        throw std::invalid_argument("to_acoustic_geometry: L_can must be positive");
    }
    if (D_can <= 0.0) {
        throw std::invalid_argument("to_acoustic_geometry: D_can must be positive");
    }

    CanAnnularGeometry geom;
    geom.n_cans = n_cans;
    geom.length_can = L_can;
    geom.area_can = M_PI * (D_can / 2.0) * (D_can / 2.0);
    geom.radius_plenum = flow_geom.D_mean() / 2.0;
    geom.area_plenum = flow_geom.area();

    return geom;
}

// BlochMode symmetry classification
std::string BlochMode::symmetry_type() const {
    if (m_azimuthal == 0) {
        return "Push-Push (all cans in phase)";
    } else if (n_cans % 2 == 0 && m_azimuthal == n_cans / 2) {
        // Only an even can count has an m = N/2 (alternating) Bloch order.
        return "Push-Pull (adjacent cans 180deg out of phase)";
    } else {
        return "Spinning/Standing wave (m=" + std::to_string(m_azimuthal) + ")";
    }
}

// Internal helper: dispersion relation of the can-annular Bloch problem.
// Overload for complex omega (needed for Argument Principle).
//
// Volume-flow balance at one can junction, outflow / p_j (same sign
// convention for every branch):
//   can (closed top):  Y_can = i Yc tan(k_c L)        Yc = A_can / (rho_c c_c)
//   can (open top):    Y_can = -i Yc cot(k_c L)
//   annulus sectors:   Y_ann = i 2Ya (cos(psi) - cos(k_p s)) / sin(k_p s)
//                      Ya = A_plenum / (rho_p c_p), s = 2 pi R / N,
//                      Bloch: p_{j+1} + p_{j-1} = 2 cos(psi) p_j, psi = 2 pi m / N
// Low-frequency check: m = 0 gives Y_ann -> i w (A s) / (rho c^2), the
// compliance of the sector volume, as a closed-top can gives i w (A L)/(rho c^2).
//
// Returned in POLE-FREE form, D = (Y_can + Y_ann) * den_can * sin(k_p s):
// the admittance sum has poles where the denominators vanish, and modes with
// a pressure node at every junction (sin(k_p s) = 0 with cos(psi) = cos(k_p s))
// are roots of this product but not of the bare sum. Entire in omega, so the
// Argument Principle counts zeros only. Units of admittance [m^3/(Pa s)].
std::complex<double> dispersion_relation_complex(
    std::complex<double> omega,
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    BoundaryCondition bc_top
) {
    const std::complex<double> i(0.0, 1.0);
    const std::complex<double> k_can = omega / c_can;
    const double Yc = geom.area_can / (rho_can * c_can);
    const std::complex<double> sin_kl_can = std::sin(k_can * geom.length_can);
    const std::complex<double> cos_kl_can = std::cos(k_can * geom.length_can);

    // Y_can = can_num / can_den
    std::complex<double> can_num;
    std::complex<double> can_den;
    if (bc_top == BoundaryCondition::Closed) {
        can_num = i * Yc * sin_kl_can;
        can_den = cos_kl_can;
    } else {
        can_num = -i * Yc * cos_kl_can;
        can_den = sin_kl_can;
    }

    // Single can: no annulus coupling, D = Y_can * can_den
    if (geom.n_cans <= 1) {
        return can_num;
    }

    const std::complex<double> k_plenum = omega / c_plenum;
    const double Ya = geom.area_plenum / (rho_plenum * c_plenum);
    const double sector_len = 2.0 * M_PI * geom.radius_plenum / geom.n_cans;
    const double psi = 2.0 * M_PI * m / geom.n_cans;
    const std::complex<double> ks = k_plenum * sector_len;

    // Y_ann = ann_num / sin(ks)
    const std::complex<double> ann_num = i * (2.0 * Ya) * (std::cos(psi) - std::cos(ks));
    return can_num * std::sin(ks) + ann_num * can_den;
}

// Real omega version (for backward compatibility)
std::complex<double> dispersion_relation(
    double omega,
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    BoundaryCondition bc_top
) {
    return dispersion_relation_complex(
        std::complex<double>(omega, 0.0), m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_top
    );
}

// Internal helper: Find approximate zero locations by scanning
std::vector<double> find_zero_guesses(
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    double f_min,
    double f_max,
    BoundaryCondition bc_top,
    int n_zeros_expected
) {
    std::vector<double> guesses;

    // Scan for local minima of |D(ω)|
    int n_scan = 500;
    std::vector<double> freqs;
    std::vector<double> mags;

    for (int i = 0; i < n_scan; ++i) {
        double f = f_min + (f_max - f_min) * static_cast<double>(i) / static_cast<double>(n_scan - 1);
        double omega = 2.0 * M_PI * f;
        auto D = dispersion_relation(omega, m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_top);

        freqs.push_back(f);
        mags.push_back(std::abs(D));
    }

    // Find local minima
    for (int i = 1; i < n_scan - 1; ++i) {
        if (mags[i] < mags[i-1] && mags[i] < mags[i+1]) {
            guesses.push_back(freqs[i]);
            if (static_cast<int>(guesses.size()) >= n_zeros_expected) break;
        }
    }

    // Fallback: if no minima found, use band center
    if (guesses.empty()) {
        guesses.push_back((f_min + f_max) / 2.0);
    }

    return guesses;
}

// Internal helper: Muller's method for complex root refinement
// More robust than Newton-Raphson for thermoacoustic problems
double refine_root_muller(
    double f_guess,
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    BoundaryCondition bc_top,
    int max_iter = 50
) {
    // Initialize three starting points around guess
    double h = 0.001 * f_guess;  // Step size
    std::complex<double> x0(2.0 * M_PI * (f_guess - h), 0.0);
    std::complex<double> x1(2.0 * M_PI * f_guess, 0.0);
    std::complex<double> x2(2.0 * M_PI * (f_guess + h), 0.0);

    auto f0 = dispersion_relation(x0.real(), m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_top);
    auto f1 = dispersion_relation(x1.real(), m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_top);
    auto f2 = dispersion_relation(x2.real(), m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_top);

    for (int iter = 0; iter < max_iter; ++iter) {
        // Muller's method: fit parabola through 3 points
        std::complex<double> q = (x2 - x1) / (x1 - x0);
        std::complex<double> A = q * f2 - q * (1.0 + q) * f1 + q * q * f0;
        std::complex<double> B = (2.0 * q + 1.0) * f2 - (1.0 + q) * (1.0 + q) * f1 + q * q * f0;
        std::complex<double> C = (1.0 + q) * f2;

        // Solve quadratic: choose root with larger denominator (more stable)
        std::complex<double> disc = std::sqrt(B * B - 4.0 * A * C);
        std::complex<double> denom1 = B + disc;
        std::complex<double> denom2 = B - disc;
        std::complex<double> denom = (std::abs(denom1) > std::abs(denom2)) ? denom1 : denom2;

        if (std::abs(denom) < 1e-14) break;  // Avoid division by zero

        std::complex<double> dx = -(x2 - x1) * 2.0 * C / denom;
        std::complex<double> x3 = x2 + dx;

        // Check convergence
        if (std::abs(dx) < 1e-6) {
            return x3.real() / (2.0 * M_PI);
        }

        // Update points
        x0 = x1;
        x1 = x2;
        x2 = x3;
        f0 = f1;
        f1 = f2;
        f2 = dispersion_relation(x2.real(), m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_top);
    }

    return x2.real() / (2.0 * M_PI);
}

// Main solver: Find all eigenmodes (delegates to solver infrastructure)
std::vector<BlochMode> can_annular_eigenmodes(
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    double f_max,
    BoundaryCondition bc_can_top
) {
    // Validate inputs
    if (geom.n_cans < 1) {
        throw std::runtime_error("Number of cans must be positive");
    }
    if (geom.length_can <= 0.0 || geom.area_can <= 0.0) {
        throw std::runtime_error("Can dimensions must be positive");
    }
    if (geom.radius_plenum <= 0.0 || geom.area_plenum <= 0.0) {
        throw std::runtime_error("Plenum dimensions must be positive");
    }
    if (c_can <= 0.0 || c_plenum <= 0.0) {
        throw std::runtime_error("Speed of sound must be positive");
    }
    if (rho_can <= 0.0 || rho_plenum <= 0.0) {
        throw std::runtime_error("Density must be positive");
    }
    if (f_max <= 0.0) {
        throw std::runtime_error("Maximum frequency must be positive");
    }

    // Delegate to the canonical magnitude-minimization solver.
    return solve_magnitude_minimization(geom, c_can, c_plenum, rho_can, rho_plenum, f_max, bc_can_top);
}

// -------------------------------------------------------------
// Utility functions (already in namespace)
// -------------------------------------------------------------

double wavelength(double f, double c) {
    if (f <= 0.0) {
        throw std::invalid_argument("wavelength: frequency must be positive");
    }
    if (c <= 0.0) {
        throw std::invalid_argument("wavelength: speed of sound must be positive");
    }
    return c / f;
}

double frequency_from_wavelength(double lambda, double c) {
    if (lambda <= 0.0) {
        throw std::invalid_argument("frequency_from_wavelength: wavelength must be positive");
    }
    if (c <= 0.0) {
        throw std::invalid_argument("frequency_from_wavelength: speed of sound must be positive");
    }
    return c / lambda;
}

double acoustic_impedance(double rho, double c) {
    if (rho <= 0.0) {
        throw std::invalid_argument("acoustic_impedance: density must be positive");
    }
    if (c <= 0.0) {
        throw std::invalid_argument("acoustic_impedance: speed of sound must be positive");
    }
    return rho * c;
}

double sound_pressure_level(double p_rms, double p_ref) {
    if (p_rms < 0.0) {
        throw std::invalid_argument("sound_pressure_level: p_rms must be non-negative");
    }
    if (p_ref <= 0.0) {
        throw std::invalid_argument("sound_pressure_level: p_ref must be positive");
    }
    if (p_rms == 0.0) {
        return -std::numeric_limits<double>::infinity();  // -∞ dB for zero pressure
    }
    return 20.0 * std::log10(p_rms / p_ref);
}

double particle_velocity(double p, double rho, double c) {
    if (p < 0.0) {
        throw std::invalid_argument("particle_velocity: pressure must be non-negative");
    }
    if (rho <= 0.0) {
        throw std::invalid_argument("particle_velocity: density must be positive");
    }
    if (c <= 0.0) {
        throw std::invalid_argument("particle_velocity: speed of sound must be positive");
    }
    return p / (rho * c);
}

} // namespace combaero
