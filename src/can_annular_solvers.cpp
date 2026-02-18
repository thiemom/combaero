#include "../include/acoustics_internal.h"
#include "../include/can_annular_solvers.h"
#include "../include/math_constants.h"
#include <algorithm>
#include <cmath>
#include <complex>

namespace combaero {

namespace {

// -------------------------------------------------------------
// Generic rectangular-contour winding-number counter
// -------------------------------------------------------------
// Counts zeros of f(z) inside the rectangle
//   [2π·f_min, 2π·f_max] × [imag_min, imag_max]
// using the Argument Principle: N = ΔArg[f] / 2π.
//
// Template parameter Func: callable std::complex<double>(std::complex<double>)
// Returns the zero count, or 0 if a pole/singularity is hit on the contour.
template<typename Func>
int count_zeros_in_contour(
    double f_min, double f_max,
    double imag_min, double imag_max,
    Func dispersion_fn
) {
    std::vector<std::complex<double>> contour;
    contour.reserve(static_cast<std::size_t>(
        2 * (kContourStepsFreq + kContourStepsImag)));

    // Build counter-clockwise rectangle, no duplicate corners.
    // Bottom edge: f_min → f_max at imag_min
    for (int i = 0; i < kContourStepsFreq; ++i) {
        double f = f_min + (f_max - f_min) * i / (kContourStepsFreq - 1);
        contour.emplace_back(2.0 * M_PI * f, imag_min);
    }
    // Right edge: imag_min → imag_max at f_max (skip first corner)
    for (int i = 1; i < kContourStepsImag; ++i) {
        double im = imag_min + (imag_max - imag_min) * i / (kContourStepsImag - 1);
        contour.emplace_back(2.0 * M_PI * f_max, im);
    }
    // Top edge: f_max → f_min at imag_max (skip first corner)
    for (int i = 1; i < kContourStepsFreq; ++i) {
        double f = f_max - (f_max - f_min) * i / (kContourStepsFreq - 1);
        contour.emplace_back(2.0 * M_PI * f, imag_max);
    }
    // Left edge: imag_max → imag_min at f_min (skip first and last corner)
    for (int i = 1; i < kContourStepsImag - 1; ++i) {
        double im = imag_max - (imag_max - imag_min) * i / (kContourStepsImag - 1);
        contour.emplace_back(2.0 * M_PI * f_min, im);
    }

    double total_phase = 0.0;
    std::complex<double> prev = dispersion_fn(contour[0]);

    for (std::size_t i = 1; i <= contour.size(); ++i) {
        std::complex<double> curr = dispersion_fn(contour[i % contour.size()]);

        if (std::abs(curr) < 1e-9) {
            return 1;  // Accidental hit on a root — report one zero
        }
        std::complex<double> ratio = curr / prev;
        if (std::abs(ratio) < 1e-9 || std::isnan(ratio.real()) || std::isinf(ratio.real())) {
            return 0;  // Hit pole or singularity — abort this box
        }

        total_phase += std::arg(ratio);
        prev = curr;
    }

    return static_cast<int>(std::round(total_phase / (2.0 * M_PI)));
}

// -------------------------------------------------------------
// Duplicate-mode guard with sliding window
// -------------------------------------------------------------
// Returns true if a mode at (m, f_candidate) is already represented
// in `modes` within kDuplicateGuardHz.
bool is_duplicate_mode(const std::vector<BlochMode>& modes, int m, double f_candidate) {
    return std::any_of(modes.begin(), modes.end(), [&](const BlochMode& e) {
        return e.m_azimuthal == m &&
               std::abs(e.frequency - f_candidate) < kDuplicateGuardHz;
    });
}

bool is_duplicate_annular(const std::vector<AnnularMode>& modes, int m, double f_candidate) {
    return std::any_of(modes.begin(), modes.end(), [&](const AnnularMode& e) {
        return e.m_azimuthal == m &&
               std::abs(e.frequency - f_candidate) < kDuplicateGuardHz;
    });
}

}  // anonymous namespace

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
    const int m_max = geom.n_cans / 2;

    for (int m = 0; m <= m_max; ++m) {
        for (double f_start = kScanStartHz; f_start < f_max; f_start += kScanWindowHz) {
            double f_end = std::min(f_start + kScanWindowHz, f_max);

            auto guesses = find_zero_guesses(
                m, geom, c_can, c_plenum, rho_can, rho_plenum,
                f_start, f_end, bc_can_top, 2
            );

            for (double f_guess : guesses) {
                double f_refined = refine_root_muller(
                    f_guess, m, geom, c_can, c_plenum, rho_can, rho_plenum,
                    bc_can_top, kMullerMaxIter
                );

                if (f_refined > 1.0 && f_refined < f_max) {
                    const double omega = 2.0 * M_PI * f_refined;
                    const auto D = dispersion_relation(
                        omega, m, geom, c_can, c_plenum, rho_can, rho_plenum, bc_can_top
                    );

                    if (std::abs(D) < kDispersionTolerance &&
                        !is_duplicate_mode(modes, m, f_refined)) {
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
    return count_zeros_in_contour(f_min, f_max, imag_min, imag_max,
        [&](std::complex<double> z) {
            return dispersion_relation_complex(z, m, geom, c_can, c_plenum,
                                              rho_can, rho_plenum, bc_top);
        });
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
    const int m_max = geom.n_cans / 2;

    for (int m = 0; m <= m_max; ++m) {
        for (double f_start = kScanStartHz; f_start < f_max; f_start += kScanWindowHz) {
            const double f_end = std::min(f_start + kScanWindowHz, f_max);

            const int n_zeros = count_zeros_in_contour(
                f_start, f_end, kContourImagMin, kContourImagMax,
                [&](std::complex<double> z) {
                    return dispersion_relation_complex(z, m, geom, c_can, c_plenum,
                                                      rho_can, rho_plenum, bc_can_top);
                });

            if (n_zeros <= 0) continue;

            // Run Muller from box center plus kMullerMultiZoneRestarts perturbed
            // starting points to recover all zeros when n_zeros > 1.
            const double box_width = f_end - f_start;
            const int n_starts = std::max(n_zeros, kMullerMultiZoneRestarts);
            for (int k = 0; k < n_starts; ++k) {
                double f_guess = f_start + box_width *
                    (0.5 + kMullerPerturbFraction * (k - n_starts / 2.0) / n_starts);
                f_guess = std::max(f_start + 1.0, std::min(f_end - 1.0, f_guess));

                double f_refined = refine_root_muller(
                    f_guess, m, geom, c_can, c_plenum, rho_can, rho_plenum,
                    bc_can_top, kMullerMaxIter
                );

                if (f_refined > kScanStartHz && f_refined < f_max &&
                    !is_duplicate_mode(modes, m, f_refined)) {
                    modes.push_back({m, f_refined, geom.n_cans});
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

// Muller's method for annular duct dispersion relation.
// Mirrors refine_root_muller from acoustics.cpp but uses annular_duct_dispersion.
static double refine_root_muller_annular(
    double f_guess,
    int m,
    const Annulus& geom,
    double c,
    double rho,
    BoundaryCondition bc_ends
) {
    const double h = 0.01 * f_guess;
    std::complex<double> x0(2.0 * M_PI * (f_guess - h), 0.0);
    std::complex<double> x1(2.0 * M_PI * f_guess,       0.0);
    std::complex<double> x2(2.0 * M_PI * (f_guess + h), 0.0);

    auto eval = [&](std::complex<double> z) {
        return annular_duct_dispersion(z, m, geom, c, rho, bc_ends);
    };

    auto f0 = eval(x0);
    auto f1 = eval(x1);
    auto f2 = eval(x2);

    for (int iter = 0; iter < kMullerMaxIter; ++iter) {
        const std::complex<double> q = (x2 - x1) / (x1 - x0);
        const std::complex<double> A = q * f2 - q * (1.0 + q) * f1 + q * q * f0;
        const std::complex<double> B = (2.0 * q + 1.0) * f2 - (1.0 + q) * (1.0 + q) * f1 + q * q * f0;
        const std::complex<double> C = (1.0 + q) * f2;

        const std::complex<double> disc   = std::sqrt(B * B - 4.0 * A * C);
        const std::complex<double> denom1 = B + disc;
        const std::complex<double> denom2 = B - disc;
        const std::complex<double> denom  = (std::abs(denom1) > std::abs(denom2)) ? denom1 : denom2;

        if (std::abs(denom) < 1e-14) break;

        const std::complex<double> dx = -(x2 - x1) * 2.0 * C / denom;
        const std::complex<double> x3 = x2 + dx;

        if (std::abs(dx) < 1e-6) {
            return x3.real() / (2.0 * M_PI);
        }

        x0 = x1; x1 = x2; x2 = x3;
        f0 = f1; f1 = f2; f2 = eval(x2);
    }

    return x2.real() / (2.0 * M_PI);
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

    for (int m = 0; m <= geom.n_azimuthal_max; ++m) {
        int n_axial = 0;

        for (double f_start = kScanStartHz; f_start < f_max; f_start += kScanWindowHz) {
            const double f_end = std::min(f_start + kScanWindowHz, f_max);

            const int n_zeros = count_zeros_in_contour(
                f_start, f_end, kContourImagMin, kContourImagMax,
                [&](std::complex<double> z) {
                    return annular_duct_dispersion(z, m, geom, c, rho, bc_ends);
                });

            if (n_zeros <= 0) continue;

            // Refine with Muller's method (replaces 100-point linear scan)
            const double f_guess = 0.5 * (f_start + f_end);
            const double f_best = refine_root_muller_annular(f_guess, m, geom, c, rho, bc_ends);

            if (f_best > kScanStartHz && f_best < f_max) {
                const double omega = 2.0 * M_PI * f_best;
                const double mag = std::abs(
                    annular_duct_dispersion(std::complex<double>(omega, 0.0), m, geom, c, rho, bc_ends)
                );

                if (mag < kDispersionTolerance && !is_duplicate_annular(modes, m, f_best)) {
                    ++n_axial;
                    modes.push_back({m, n_axial, f_best});
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
