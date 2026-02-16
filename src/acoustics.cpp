#include "../include/acoustics.h"
#include "../include/math_constants.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

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

    // Neck diameter and perimeter
    double d_neck = 2.0 * std::sqrt(A_neck / M_PI);
    double perimeter = M_PI * d_neck;

    // Effective neck length (with standard end correction)
    double L_eff = L_neck + 0.85 * d_neck;

    // Losses occur in the neck boundary layer
    // Q ≈ (neck volume) / (boundary layer volume in neck)
    // Q ≈ (A_neck * L_eff) / (perimeter * L_eff * delta_eff)
    // Q ≈ A_neck / (perimeter * delta_eff)
    // Q ≈ d_neck / (4 * delta_eff)  for circular neck
    //
    // But cavity also stores energy, so scale by V/(A_neck * L_eff)
    // This gives the classic result: Q ~ V / (A_neck * delta_eff * factor)

    // Simplified model: Q ≈ d_neck / (4 * delta_eff) * sqrt(V / (A_neck * L_eff))
    // This captures both neck losses and cavity energy storage
    double Q = (d_neck / (4.0 * delta_eff)) * std::sqrt(V / (A_neck * L_eff));

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
    // Q ≈ (tube cross-section) / (boundary layer area)
    // Q ≈ (π D²/4) / (π D * delta_eff)
    // Q ≈ D / (4 * delta_eff)
    double Q = D / (4.0 * delta_eff);

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
    double rho,
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

// Orifice impedance with bias and grazing flow
std::complex<double> orifice_impedance_with_flow(
    double freq,
    double u_bias,
    double u_grazing,
    double d_orifice,
    double l_orifice,
    double porosity,
    double Cd,
    double rho,
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
    // Rayleigh end correction modified by flow
    double delta_0 = 0.85 * (d_orifice / 2.0);  // ~0.425 * d
    
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
    
    // 3. Total branch impedance
    // Series combination: Z_branch = Z_neck + Z_tube
    std::complex<double> Z_branch_norm = Z_neck + Z_tube;
    
    // 4. Scale by area ratio and denormalize
    // The orifice area is: A_orifice = porosity * A_duct
    // Branch admittance in main duct: Y = A_orifice / (ρc * Z_branch)
    double A_orifice = porosity * A_duct;
    std::complex<double> Y_branch = A_orifice / (rho * c * Z_branch_norm);
    
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
