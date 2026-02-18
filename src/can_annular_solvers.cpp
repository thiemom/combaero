#include "../include/can_annular_solvers.h"
#include "../include/math_constants.h"
#include <algorithm>
#include <cmath>
#include <complex>

namespace combaero {

// Forward declarations of helpers from acoustics.cpp
extern std::complex<double> dispersion_relation_complex(
    std::complex<double> omega,
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    BoundaryCondition bc_top
);

extern std::complex<double> dispersion_relation(
    double omega,
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    BoundaryCondition bc_top
);

extern std::vector<double> find_zero_guesses(
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
);

extern double refine_root_muller(
    double f_guess,
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    BoundaryCondition bc_top,
    int max_iter
);


// -------------------------------------------------------------
// Main Solver with Method Selection
// -------------------------------------------------------------

std::vector<BlochMode> can_annular_eigenmodes_with_method(
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    double f_max,
    BoundaryCondition bc_can_top,
    CanAnnularSolverMethod method
) {
    if (method == CanAnnularSolverMethod::ArgumentPrinciple) {
        return solve_argument_principle(geom, c_can, c_plenum, rho_can, rho_plenum, f_max, bc_can_top);
    } else {
        return solve_magnitude_minimization(geom, c_can, c_plenum, rho_can, rho_plenum, f_max, bc_can_top);
    }
}

// -------------------------------------------------------------
// Method 1: Magnitude Minimization (CURRENT WORKING METHOD)
// -------------------------------------------------------------

std::vector<BlochMode> solve_magnitude_minimization(
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    double f_max,
    BoundaryCondition bc_can_top
) {
    std::vector<BlochMode> modes;
    int m_max = geom.n_cans / 2;

    // Continuous sliding window scan
    double scan_window_width = 50.0;  // Hz per window

    for (int m = 0; m <= m_max; ++m) {
        for (double f_start = 10.0; f_start < f_max; f_start += scan_window_width) {
            double f_end = std::min(f_start + scan_window_width, f_max);

            auto guesses = find_zero_guesses(
                m, geom, c_can, c_plenum, rho_can, rho_plenum,
                f_start, f_end, bc_can_top, 2
            );

            for (double f_guess : guesses) {
                double f_refined = refine_root_muller(
                    f_guess, m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_can_top, 50
                );

                if (f_refined > 1.0 && f_refined < f_max) {
                    double omega_refined = 2.0 * M_PI * f_refined;
                    auto D_refined = dispersion_relation(
                        omega_refined, m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_can_top
                    );

                    if (std::abs(D_refined) < 0.1) {
                        bool is_duplicate = std::any_of(modes.begin(), modes.end(),
                            [&](const BlochMode& e) {
                                return e.m_azimuthal == m &&
                                       std::abs(e.frequency - f_refined) < 1.0;
                            });

                        if (!is_duplicate) {
                            modes.push_back({m, f_refined, geom.n_cans});
                        }
                    }
                }
            }
        }
    }

    std::sort(modes.begin(), modes.end(),
              [](const BlochMode& a, const BlochMode& b) { return a.frequency < b.frequency; });

    return modes;
}

// -------------------------------------------------------------
// Method 2: Argument Principle (RESEARCH-GRADE)
// -------------------------------------------------------------

std::vector<double> identify_pole_frequencies(
    const CanAnnularGeometry& geom,
    double c_can,
    double f_min,
    double f_max,
    BoundaryCondition bc_top
) {
    std::vector<double> poles;

    // For closed boundary: Y = i*Y_char*tan(kL)
    // Poles when cos(kL) = 0, i.e., kL = (n + 1/2)*π
    // f_pole = (n + 1/2) * c / (2*L)

    if (bc_top == BoundaryCondition::Closed) {
        for (int n = 0; n < 50; ++n) {
            double f_pole = (n + 0.5) * c_can / (2.0 * geom.length_can);
            if (f_pole >= f_min && f_pole <= f_max) {
                poles.push_back(f_pole);
            }
        }
    } else {
        // Open boundary: Y = -i*Y_char*cot(kL)
        // Poles when sin(kL) = 0, i.e., kL = n*π
        // f_pole = n * c / (2*L)
        for (int n = 1; n < 50; ++n) {
            double f_pole = n * c_can / (2.0 * geom.length_can);
            if (f_pole >= f_min && f_pole <= f_max) {
                poles.push_back(f_pole);
            }
        }
    }

    return poles;
}

int count_zeros_argument_principle_improved(
    double f_min,
    double f_max,
    double imag_min,
    double imag_max,
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    BoundaryCondition bc_top
) {
    // Increase resolution for reliability
    const int n_steps_f = 40;
    const int n_steps_i = 20;

    std::vector<std::complex<double>> contour;
    contour.reserve(2U * (static_cast<std::size_t>(n_steps_f) + static_cast<std::size_t>(n_steps_i)));

    // Build rectangular contour (counter-clockwise)
    // Bottom edge (real axis approximation)
    for (int i = 0; i < n_steps_f; ++i) {
        double f = f_min + (f_max - f_min) * i / n_steps_f;
        contour.push_back(std::complex<double>(2.0 * M_PI * f, imag_min));
    }
    // Right edge
    for (int i = 0; i < n_steps_i; ++i) {
        double imag = imag_min + (imag_max - imag_min) * i / n_steps_i;
        contour.push_back(std::complex<double>(2.0 * M_PI * f_max, imag));
    }
    // Top edge
    for (int i = 0; i < n_steps_f; ++i) {
        double f = f_max - (f_max - f_min) * i / n_steps_f;
        contour.push_back(std::complex<double>(2.0 * M_PI * f, imag_max));
    }
    // Left edge
    for (int i = 0; i < n_steps_i; ++i) {
        double imag = imag_max - (imag_max - imag_min) * i / n_steps_i;
        contour.push_back(std::complex<double>(2.0 * M_PI * f_min, imag));
    }

    // Compute winding number
    double total_phase_change = 0.0;
    std::complex<double> prev_val = dispersion_relation_complex(
        contour[0], m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_top
    );

    for (size_t i = 1; i <= contour.size(); ++i) {
        std::complex<double> curr_z = contour[i % contour.size()];
        std::complex<double> curr_val = dispersion_relation_complex(
            curr_z, m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_top
        );

        // Robust phase accumulation: Im(log(z2/z1))
        // Handles branch cuts better than std::arg(z2) - std::arg(z1)
        if (std::abs(curr_val) < 1e-9) {
            return 1;  // Accidental hit on a root
        }

        std::complex<double> ratio = curr_val / prev_val;
        if (std::abs(ratio) < 1e-9 || std::isinf(ratio.real())) {
            return 0;  // Hit pole or singularity - abort this box
        }

        total_phase_change += std::arg(ratio);
        prev_val = curr_val;
    }

    // N = Δθ / 2π
    return static_cast<int>(std::round(total_phase_change / (2.0 * M_PI)));
}

std::vector<BlochMode> solve_argument_principle(
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    double f_max,
    BoundaryCondition bc_can_top
) {
    std::vector<BlochMode> modes;
    int m_max = geom.n_cans / 2;

    // Grid configuration
    double box_width = 50.0;   // Hz (width of search box)
    double imag_min = -10.0;   // Look slightly below real axis (stable modes)
    double imag_max = 200.0;   // Look deep into complex plane (damped modes)

    for (int m = 0; m <= m_max; ++m) {
        for (double f_start = 10.0; f_start < f_max; f_start += box_width) {
            double f_end = std::min(f_start + box_width, f_max);

            // Use Argument Principle to count zeros
            int n_zeros = count_zeros_argument_principle_improved(
                f_start, f_end, imag_min, imag_max,
                m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_can_top
            );

            if (n_zeros > 0) {
                // Use Muller's method starting from center of box
                double f_center = (f_start + f_end) / 2.0;

                // Can run Muller multiple times with perturbations if n_zeros > 1
                double f_refined = refine_root_muller(
                    f_center, m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_can_top, 50
                );

                if (f_refined > 0.0 && f_refined < f_max) {
                    bool is_duplicate = std::any_of(modes.begin(), modes.end(),
                        [&](const BlochMode& e) {
                            return e.m_azimuthal == m &&
                                   std::abs(e.frequency - f_refined) < 1.0;
                        });

                    if (!is_duplicate) {
                        modes.push_back({m, f_refined, geom.n_cans});
                    }
                }
            }
        }
    }

    std::sort(modes.begin(), modes.end(),
              [](const BlochMode& a, const BlochMode& b) { return a.frequency < b.frequency; });

    return modes;
}

// -------------------------------------------------------------
// Annular Duct Modes (Pure Annular Geometry)
// -------------------------------------------------------------

// Helper: Annular duct dispersion relation
// For a simple annular duct: D(ω) = Y_duct(ω, m)
static std::complex<double> annular_duct_dispersion(
    std::complex<double> omega,
    int m,
    const Annulus& geom,
    double c,
    double rho,
    BoundaryCondition bc_ends
) {
    std::complex<double> k = omega / c;
    std::complex<double> i(0.0, 1.0);

    // Thin-annulus approximation: split total wavenumber into azimuthal and axial parts
    const double r_mean = 0.5 * (geom.radius_inner() + geom.radius_outer());
    const double k_theta = (r_mean > 0.0) ? static_cast<double>(m) / r_mean : 0.0;
    const std::complex<double> k_axial = std::sqrt(k * k - k_theta * k_theta);

    if (std::abs(k) < 1e-12 || std::abs(k_axial) < 1e-12) {
        return std::complex<double>(1e12, 0.0);
    }

    // Duct admittance (similar to can admittance)
    std::complex<double> Y_char = (geom.area() / (rho * c)) * (k_axial / k);
    std::complex<double> sin_kl = std::sin(k_axial * geom.L);
    std::complex<double> cos_kl = std::cos(k_axial * geom.L);

    std::complex<double> Y_duct;
    if (bc_ends == BoundaryCondition::Closed) {
        // Y = i * Y_char * tan(kL)
        if (std::abs(cos_kl) < 1e-12) {
            return std::complex<double>(1e12, 0.0);  // Pole
        }
        Y_duct = i * Y_char * (sin_kl / cos_kl);
    } else {
        // Y = -i * Y_char * cot(kL)
        if (std::abs(sin_kl) < 1e-12) {
            return std::complex<double>(1e12, 0.0);  // Pole
        }
        Y_duct = -i * Y_char * (cos_kl / sin_kl);
    }

    // For annular duct, dispersion relation is the modal admittance
    // (zeros occur at resonance frequencies)
    return Y_duct;
}

// Helper: Count zeros for annular duct using Argument Principle
static int count_zeros_annular_duct(
    double f_min,
    double f_max,
    double imag_min,
    double imag_max,
    int m,
    const Annulus& geom,
    double c,
    double rho,
    BoundaryCondition bc_ends
) {
    const int n_steps_f = 40;
    const int n_steps_i = 20;

    std::vector<std::complex<double>> contour;
    contour.reserve(2U * (static_cast<std::size_t>(n_steps_f) + static_cast<std::size_t>(n_steps_i)));

    // Build rectangular contour
    for (int i = 0; i < n_steps_f; ++i) {
        double f = f_min + (f_max - f_min) * i / n_steps_f;
        contour.push_back(std::complex<double>(2.0 * M_PI * f, imag_min));
    }
    for (int i = 0; i < n_steps_i; ++i) {
        double imag = imag_min + (imag_max - imag_min) * i / n_steps_i;
        contour.push_back(std::complex<double>(2.0 * M_PI * f_max, imag));
    }
    for (int i = 0; i < n_steps_f; ++i) {
        double f = f_max - (f_max - f_min) * i / n_steps_f;
        contour.push_back(std::complex<double>(2.0 * M_PI * f, imag_max));
    }
    for (int i = 0; i < n_steps_i; ++i) {
        double imag = imag_max - (imag_max - imag_min) * i / n_steps_i;
        contour.push_back(std::complex<double>(2.0 * M_PI * f_min, imag));
    }

    // Compute winding number
    double total_phase_change = 0.0;
    std::complex<double> prev_val = annular_duct_dispersion(
        contour[0], m, geom, c, rho, bc_ends
    );

    for (size_t i = 1; i <= contour.size(); ++i) {
        std::complex<double> curr_z = contour[i % contour.size()];
        std::complex<double> curr_val = annular_duct_dispersion(
            curr_z, m, geom, c, rho, bc_ends
        );

        if (std::abs(curr_val) < 1e-9) {
            return 1;  // Hit a root
        }

        std::complex<double> ratio = curr_val / prev_val;
        if (std::abs(ratio) < 1e-9 || std::isinf(ratio.real())) {
            return 0;  // Hit pole - abort
        }

        total_phase_change += std::arg(ratio);
        prev_val = curr_val;
    }

    return static_cast<int>(std::round(total_phase_change / (2.0 * M_PI)));
}

std::vector<AnnularMode> annular_duct_eigenmodes(
    const Annulus& geom,
    double c,
    double rho,
    double f_max,
    BoundaryCondition bc_ends
) {
    // Validate inputs
    if (geom.L <= 0.0) {
        throw std::runtime_error("Duct length must be positive");
    }
    if (geom.radius_inner() < 0.0 || geom.radius_outer() <= geom.radius_inner()) {
        throw std::runtime_error("Invalid annular geometry");
    }
    if (c <= 0.0 || rho <= 0.0) {
        throw std::runtime_error("Speed of sound and density must be positive");
    }
    if (f_max <= 0.0) {
        throw std::runtime_error("Maximum frequency must be positive");
    }

    std::vector<AnnularMode> modes;

    // Grid search parameters
    double box_width = 50.0;   // Hz
    double imag_min = -10.0;
    double imag_max = 200.0;

    // Loop over azimuthal modes
    for (int m = 0; m <= geom.n_azimuthal_max; ++m) {
        int n_axial = 0;  // Track axial mode number

        for (double f_start = 10.0; f_start < f_max; f_start += box_width) {
            double f_end = std::min(f_start + box_width, f_max);

            // Use Argument Principle to count zeros
            int n_zeros = count_zeros_annular_duct(
                f_start, f_end, imag_min, imag_max,
                m, geom, c, rho, bc_ends
            );

            if (n_zeros > 0) {
                // Found a mode - estimate frequency
                // For closed-closed: f ≈ n * c / (2*L)
                // For closed-open: f ≈ (2n-1) * c / (4*L)
                double f_estimate = (f_start + f_end) / 2.0;

                // Simple refinement: scan for minimum |D|
                double f_best = f_estimate;
                double min_mag = 1e10;

                for (int i = 0; i < 100; ++i) {
                    double f_test = f_start + (f_end - f_start) * i / 99.0;
                    double omega_test = 2.0 * M_PI * f_test;
                    auto D_test = annular_duct_dispersion(
                        std::complex<double>(omega_test, 0.0),
                        m, geom, c, rho, bc_ends
                    );
                    double mag = std::abs(D_test);
                    if (mag < min_mag) {
                        min_mag = mag;
                        f_best = f_test;
                    }
                }

                // Check if it's a good zero
                if (min_mag < 0.1) {
                    n_axial++;

                    bool is_duplicate = std::any_of(modes.begin(), modes.end(),
                        [&](const AnnularMode& e) {
                            return e.m_azimuthal == m &&
                                   std::abs(e.frequency - f_best) < 1.0;
                        });

                    if (!is_duplicate) {
                        modes.push_back({m, n_axial, f_best});
                    }
                }
            }
        }
    }

    std::sort(modes.begin(), modes.end(),
              [](const AnnularMode& a, const AnnularMode& b) { return a.frequency < b.frequency; });

    return modes;
}

std::vector<AnnularMode> annular_duct_modes_analytical(
    const Annulus& geom,
    double c,
    double f_max,
    BoundaryCondition bc_ends
) {
    if (geom.L <= 0.0) {
        throw std::runtime_error("Duct length must be positive");
    }
    if (geom.radius_inner() < 0.0 || geom.radius_outer() <= geom.radius_inner()) {
        throw std::runtime_error("Invalid annular geometry");
    }
    if (c <= 0.0) {
        throw std::runtime_error("Speed of sound must be positive");
    }
    if (f_max <= 0.0) {
        throw std::runtime_error("Maximum frequency must be positive");
    }

    std::vector<AnnularMode> modes;
    const double r_mean = 0.5 * (geom.radius_inner() + geom.radius_outer());
    const double d_mean = 2.0 * r_mean;

    for (int m = 0; m <= geom.n_azimuthal_max; ++m) {
        const double f_az = (m == 0) ? 0.0 : (m * c) / (M_PI * d_mean);

        for (int n = 1; n < 200; ++n) {
            double f_axial = 0.0;
            if (bc_ends == BoundaryCondition::Closed) {
                f_axial = n * c / (2.0 * geom.L);
            } else {
                f_axial = (2.0 * n - 1.0) * c / (4.0 * geom.L);
            }

            const double frequency = std::sqrt(f_axial * f_axial + f_az * f_az);
            if (frequency > f_max) {
                break;
            }

            modes.push_back({m, n, frequency});
        }
    }

    std::sort(
        modes.begin(),
        modes.end(),
        [](const AnnularMode& a, const AnnularMode& b) { return a.frequency < b.frequency; }
    );

    return modes;
}

}  // namespace combaero
