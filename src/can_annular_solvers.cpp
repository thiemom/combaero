#include "../include/acoustics_internal.h"
#include "../include/can_annular_solvers.h"
#include "../include/math_constants.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace combaero {

namespace {

// -------------------------------------------------------------
// Generic rectangular-contour winding-number counter
// -------------------------------------------------------------
// Winding number of f(z) around the rectangle
//   [2pi*f_min, 2pi*f_max] x [imag_min, imag_max]
// by the Argument Principle: Delta Arg[f] / 2pi = Z - P, zeros minus poles
// inside. A caller wanting Z must add the poles back.
//
// Template parameter Func: callable std::complex<double>(std::complex<double>)
// Returns Z - P, or 0 if a pole/singularity is hit on the contour.
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
    // Bottom edge: f_min -> f_max at imag_min
    for (int i = 0; i < kContourStepsFreq; ++i) {
        double f = f_min + (f_max - f_min) * i / (kContourStepsFreq - 1);
        contour.emplace_back(2.0 * M_PI * f, imag_min);
    }
    // Right edge: imag_min -> imag_max at f_max (skip first corner)
    for (int i = 1; i < kContourStepsImag; ++i) {
        double im = imag_min + (imag_max - imag_min) * i / (kContourStepsImag - 1);
        contour.emplace_back(2.0 * M_PI * f_max, im);
    }
    // Top edge: f_max -> f_min at imag_max (skip first corner)
    for (int i = 1; i < kContourStepsFreq; ++i) {
        double f = f_max - (f_max - f_min) * i / (kContourStepsFreq - 1);
        contour.emplace_back(2.0 * M_PI * f, imag_max);
    }
    // Left edge: imag_max -> imag_min at f_min (skip first and last corner)
    for (int i = 1; i < kContourStepsImag - 1; ++i) {
        double im = imag_max - (imag_max - imag_min) * i / (kContourStepsImag - 1);
        contour.emplace_back(2.0 * M_PI * f_min, im);
    }

    // Phase increment along one segment. A zero or pole close to the contour
    // turns the argument by up to pi between two coarse samples, and arg()
    // then aliases it by 2 pi; bisect until each step turns less than
    // kContourMaxPhaseStep (or kContourMaxBisections is reached).
    bool hit_root = false;
    bool hit_pole = false;
    auto segment_phase = [&](auto&& self, std::complex<double> z0, std::complex<double> d0,
                             std::complex<double> z1, std::complex<double> d1,
                             int depth) -> double {
        const std::complex<double> ratio = d1 / d0;
        if (std::abs(ratio) < 1e-9 || std::isnan(ratio.real()) || std::isinf(ratio.real())) {
            hit_pole = true;
            return 0.0;
        }
        const double dphi = std::arg(ratio);
        if (std::abs(dphi) <= kContourMaxPhaseStep || depth >= kContourMaxBisections) {
            return dphi;
        }
        const std::complex<double> zm = 0.5 * (z0 + z1);
        const std::complex<double> dm = dispersion_fn(zm);
        if (std::abs(dm) < 1e-9) {
            hit_root = true;
            return 0.0;
        }
        return self(self, z0, d0, zm, dm, depth + 1) + self(self, zm, dm, z1, d1, depth + 1);
    };

    double total_phase = 0.0;
    std::complex<double> prev = dispersion_fn(contour[0]);

    for (std::size_t i = 1; i <= contour.size(); ++i) {
        std::complex<double> curr = dispersion_fn(contour[i % contour.size()]);

        if (std::abs(curr) < 1e-9) {
            return 1;  // Accidental hit on a root -- report one zero
        }
        total_phase += segment_phase(segment_phase, contour[i - 1], prev,
                                     contour[i % contour.size()], curr, 0);
        if (hit_root) return 1;
        if (hit_pole) return 0;  // Hit pole or singularity -- abort this box
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

// Characteristic admittance of the can/plenum junction [m^4 s/kg]: the scale
// that makes |D| dimensionless for the kDispersionTolerance acceptance test.
double junction_admittance_scale(
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum
) {
    double y = geom.area_can / (rho_can * c_can);
    if (geom.n_cans > 1) {
        y += 2.0 * geom.area_plenum / (rho_plenum * c_plenum);
    }
    return y;
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
    const double y_scale = junction_admittance_scale(geom, c_can, c_plenum, rho_can, rho_plenum);

    for (int m = 0; m <= m_max; ++m) {
        for (double f_start = kScanStartHz; f_start < f_max; f_start += kScanWindowHz) {
            double f_end = std::min(f_start + kScanWindowHz, f_max);

            auto guesses = find_zero_guesses(
                m, geom, c_can, c_plenum, rho_can, rho_plenum,
                f_start, f_end, bc_can_top, kMagnitudeScanMaxGuesses
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

                    if (std::abs(D) < kDispersionTolerance * y_scale &&
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
    // Poles when cos(kL) = 0, i.e., kL = (n + 1/2)*pi
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
        // Poles when sin(kL) = 0, i.e., kL = n*pi
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
    const double y_scale = junction_admittance_scale(geom, c_can, c_plenum, rho_can, rho_plenum);

    for (int m = 0; m <= m_max; ++m) {
        for (double f_start = kScanStartHz; f_start < f_max; f_start += kScanWindowHz) {
            const double f_end = std::min(f_start + kScanWindowHz, f_max);

            // dispersion_relation_complex is pole-free (entire in omega), so the
            // winding number counts zeros only.
            const double f_lo = f_start;
            const double f_hi = f_end;
            const int n_zeros = count_zeros_in_contour(
                f_lo, f_hi, kContourImagMin, kContourImagMax,
                [&](std::complex<double> z) {
                    return dispersion_relation_complex(z, m, geom, c_can, c_plenum,
                                                      rho_can, rho_plenum, bc_can_top);
                });

            if (n_zeros <= 0) continue;

            // Muller starts: the |D| minima of the box (a root near a box edge
            // is far from the centre), then the box centre plus
            // kMullerMultiZoneRestarts perturbed points.
            const double box_width = f_hi - f_lo;
            const int n_starts = std::max(n_zeros, kMullerMultiZoneRestarts);
            std::vector<double> starts = find_zero_guesses(
                m, geom, c_can, c_plenum, rho_can, rho_plenum,
                f_lo, f_hi, bc_can_top, n_zeros);
            for (int k = 0; k < n_starts; ++k) {
                double f_guess = f_lo + box_width *
                    (0.5 + kMullerPerturbFraction * (k - n_starts / 2.0) / n_starts);
                starts.push_back(std::max(f_lo + 1.0, std::min(f_hi - 1.0, f_guess)));
            }
            for (double f_guess : starts) {
                const double f_refined = refine_root_muller(
                    f_guess, m, geom, c_can, c_plenum, rho_can, rho_plenum,
                    bc_can_top, kMullerMaxIter
                );
                if (f_refined <= kScanStartHz || f_refined >= f_max) continue;

                // Muller returns its last iterate whether or not it converged:
                // accept only a point that is a root.
                const auto D = dispersion_relation(
                    2.0 * M_PI * f_refined, m, geom, c_can, c_plenum, rho_can, rho_plenum,
                    bc_can_top);
                if (std::abs(D) < kDispersionTolerance * y_scale &&
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
// For a simple annular duct: D(omega) = Y_duct(omega, m)
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

    if (std::abs(k) < 1e-12) {
        return std::complex<double>(1e12, 0.0);
    }

    // Near the cut-on (k_axial -> 0) evaluate k_x tan(k_x L) and k_x cot(k_x L)
    // by their series: both are even in k_x and finite there, and the
    // rigid-rigid n = 0 azimuthal mode is a root AT k_x = 0.
    const std::complex<double> x = k_axial * geom.L;
    if (std::abs(x) < kAnnularCutOnSeries) {
        const std::complex<double> x2 = x * x;
        const double y0 = geom.area() / (rho * c);
        if (bc_ends == BoundaryCondition::Closed) {
            // i Y0 (k_x / k) tan(k_x L) = i Y0 x^2 (1 + x^2/3) / (k L)
            return i * y0 * x2 * (1.0 + x2 / 3.0) / (k * geom.L);
        }
        // -i Y0 (k_x / k) cot(k_x L) = -i Y0 (1 - x^2/3) / (k L)
        return -i * y0 * (1.0 - x2 / 3.0) / (k * geom.L);
    }

    // Duct admittance (similar to can admittance)
    std::complex<double> Y_char = (geom.area() / (rho * c)) * (k_axial / k);
    std::complex<double> sin_kl = std::sin(x);
    std::complex<double> cos_kl = std::cos(x);

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
// h [Hz] is the half-width of the initial three-point stencil: it must stay
// inside the root's pole-free interval (the cut-on root and the next pole can
// be a few Hz apart), so callers pass their scan spacing.
static double refine_root_muller_annular(
    double f_guess,
    double h,
    int m,
    const Annulus& geom,
    double c,
    double rho,
    BoundaryCondition bc_ends
) {
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

// Number of real poles of annular_duct_dispersion with f_lo < f < f_hi [Hz].
// Y_duct ~ k_x tan(k_x L) (closed) or k_x cot(k_x L) (open), k_x the axial
// wavenumber at f^2 = f_az^2 + f_x^2:
//   closed : cos(k_x L) = 0  ->  f_x = (2n + 1) c / (4 L), n >= 0
//   open   : sin(k_x L) = 0  ->  f_x = n c / (2 L),        n >= 1
// (k_x = 0 is not a pole: the k_x prefactor cancels cot's 1/k_x.)
static int count_annular_duct_poles(
    double f_lo, double f_hi, double f_az, double L, double c, bool closed
) {
    int count = 0;
    for (int n = closed ? 0 : 1;; ++n) {
        const double f_x = closed ? (2.0 * n + 1.0) * c / (4.0 * L) : n * c / (2.0 * L);
        const double f = std::sqrt(f_x * f_x + f_az * f_az);
        if (f >= f_hi) break;
        if (f > f_lo) ++count;
    }
    return count;
}

std::vector<AnnularMode> annular_duct_eigenmodes(
    const Annulus& geom,
    double c,
    double rho,
    double f_max,
    BoundaryCondition bc_ends
) {
    // Validate inputs
    if (!(geom.L > 0.0)) {
        throw std::runtime_error("Duct length must be positive");
    }
    if (!(geom.radius_inner() >= 0.0 && geom.radius_outer() > geom.radius_inner())) {
        throw std::runtime_error("Invalid annular geometry");
    }
    if (!(c > 0.0 && rho > 0.0)) {
        throw std::runtime_error("Speed of sound and density must be positive");
    }
    if (!(f_max > 0.0) || !std::isfinite(f_max)) {
        throw std::runtime_error("Maximum frequency must be positive and finite");
    }

    std::vector<AnnularMode> modes;
    const bool closed = (bc_ends == BoundaryCondition::Closed);
    const double d_mean = geom.D_mean();
    const double y_scale = geom.area() / (rho * c);

    for (int m = 0; m <= geom.n_azimuthal_max; ++m) {
        const double f_az = (m == 0) ? 0.0 : (m * c) / (M_PI * d_mean);
        auto eval = [&](std::complex<double> z) {
            return annular_duct_dispersion(z, m, geom, c, rho, bc_ends);
        };

        for (double f_start = kScanStartHz; f_start < f_max; f_start += kScanWindowHz) {
            const double f_end = std::min(f_start + kScanWindowHz, f_max);

            // Same pole handling as solve_argument_principle: widen an edge
            // that sits on a pole, then add the poles back to Z - P.
            double f_lo = f_start;
            double f_hi = f_end;
            if (count_annular_duct_poles(f_lo - kContourPoleClearanceHz,
                                         f_lo + kContourPoleClearanceHz, f_az, geom.L, c,
                                         closed) > 0) {
                f_lo -= 2.0 * kContourPoleClearanceHz;
            }
            if (count_annular_duct_poles(f_hi - kContourPoleClearanceHz,
                                         f_hi + kContourPoleClearanceHz, f_az, geom.L, c,
                                         closed) > 0) {
                f_hi += 2.0 * kContourPoleClearanceHz;
            }
            const int n_zeros =
                count_zeros_in_contour(f_lo, f_hi, kContourImagMin, kContourImagMax, eval) +
                count_annular_duct_poles(f_lo, f_hi, f_az, geom.L, c, closed);

            if (n_zeros <= 0) continue;

            // Muller starts: |D| minima across the box, then the box centre.
            std::vector<double> starts;
            const double df = (f_hi - f_lo) / (kAnnularScanPoints - 1);
            double mag_prev = std::abs(eval(2.0 * M_PI * f_lo));
            double mag_curr = std::abs(eval(2.0 * M_PI * (f_lo + df)));
            for (int j = 2; j < kAnnularScanPoints; ++j) {
                const double mag_next = std::abs(eval(2.0 * M_PI * (f_lo + j * df)));
                if (mag_curr < mag_prev && mag_curr < mag_next) {
                    starts.push_back(f_lo + (j - 1) * df);
                }
                mag_prev = mag_curr;
                mag_curr = mag_next;
            }
            starts.push_back(0.5 * (f_lo + f_hi));

            for (double f_guess : starts) {
                const double f_best = refine_root_muller_annular(f_guess, df, m, geom, c, rho, bc_ends);
                if (f_best <= kScanStartHz || f_best >= f_max * (1.0 - kFmaxEdgeRelTol)) continue;

                const double mag = std::abs(eval(std::complex<double>(2.0 * M_PI * f_best, 0.0)));
                if (mag < kDispersionTolerance * y_scale &&
                    !is_duplicate_annular(modes, m, f_best)) {
                    // Axial order from the root itself, not from discovery order.
                    const double f_ax = std::sqrt(std::max(0.0, f_best * f_best - f_az * f_az));
                    const int n_axial = closed
                        ? static_cast<int>(std::lround(2.0 * geom.L * f_ax / c))
                        : static_cast<int>(std::lround(0.5 * (4.0 * geom.L * f_ax / c + 1.0)));
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
    if (!(geom.L > 0.0)) {
        throw std::runtime_error("Duct length must be positive");
    }
    if (!(geom.radius_inner() >= 0.0 && geom.radius_outer() > geom.radius_inner())) {
        throw std::runtime_error("Invalid annular geometry");
    }
    if (!(c > 0.0)) {
        throw std::runtime_error("Speed of sound must be positive");
    }
    if (!(f_max > 0.0) || !std::isfinite(f_max)) {
        throw std::runtime_error("Maximum frequency must be positive and finite");
    }

    std::vector<AnnularMode> modes;
    const double r_mean = 0.5 * (geom.radius_inner() + geom.radius_outer());
    const double d_mean = 2.0 * r_mean;

    const bool closed = (bc_ends == BoundaryCondition::Closed);
    for (int m = 0; m <= geom.n_azimuthal_max; ++m) {
        const double f_az = (m == 0) ? 0.0 : (m * c) / (M_PI * d_mean);

        // Rigid-rigid admits n = 0 for m >= 1: the axially uniform, purely
        // azimuthal mode at f = f_az. Rigid-open starts at the quarter wave.
        for (int n = (closed && m > 0) ? 0 : 1;; ++n) {
            const double f_axial = closed ? n * c / (2.0 * geom.L)
                                          : (2.0 * n - 1.0) * c / (4.0 * geom.L);
            const double frequency = std::sqrt(f_axial * f_axial + f_az * f_az);
            if (frequency >= f_max) {
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
