#include "../include/acoustics.h"
#include "../include/math_constants.h"
#include <algorithm>
#include <cmath>
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
