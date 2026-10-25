#pragma once

#include "acoustics.h"
#include "geometry.h"  // For Annulus
#include "math_constants.h"  // MSVC compatibility for M_PI
#include <vector>
#include <complex>

namespace combaero {

// -------------------------------------------------------------
// Solver tuning constants
// -------------------------------------------------------------
// All magic numbers are collected here with physical justification.

// Width of each frequency search window [Hz].
// Narrower windows reduce the chance of two modes sharing a box;
// 50 Hz is a safe default for combustor-scale geometries (L ~ 0.3-1 m).
inline constexpr double kScanWindowHz = 50.0;

// Minimum frequency to start scanning [Hz].
// Avoids DC / near-zero numerical artefacts in the dispersion relation.
inline constexpr double kScanStartHz = 10.0;

// Imaginary-axis bounds for the Argument Principle contour [rad/s].
// imag_min slightly below real axis to capture marginally stable modes.
// imag_max large enough to enclose all physically relevant damped modes.
inline constexpr double kContourImagMin = -10.0;
inline constexpr double kContourImagMax = 200.0;

// Minimum distance [Hz] between a real-axis pole of D and a vertical contour
// edge. A box edge closer than this to a pole is moved outward by twice this
// value, so the pole sits well inside the box and is counted.
inline constexpr double kContourPoleClearanceHz = 5.0;

// Points of the |D| scan that seeds Muller inside one annular-duct box
// (local minima become starting points). ~0.25 Hz spacing on a 50 Hz box.
inline constexpr int kAnnularScanPoints = 201;

// Magnitude minimization: every |D| minimum of a scan window seeds Muller
// (previously only the first two, which dropped the third root of a busy
// window).
inline constexpr int kMagnitudeScanMaxGuesses = 1000;

// |k_x L| below which the annular-duct admittance uses its cut-on series
// (two terms; relative truncation error ~ x^4 < 1e-16 here).
inline constexpr double kAnnularCutOnSeries = 1e-4;

// f_max is exclusive. A refined root within this relative distance below it is
// treated as AT f_max (Muller lands on 999.9999999 for a 1000 Hz root), so the
// numerical and analytical annular solvers return the same set.
inline constexpr double kFmaxEdgeRelTol = 1e-9;

// Number of contour steps along the frequency and imaginary axes.
// Higher values reduce the chance of missing a winding due to coarse sampling.
inline constexpr int kContourStepsFreq = 40;
inline constexpr int kContourStepsImag = 20;

// Adaptive contour refinement: a segment whose argument turns by more than
// kContourMaxPhaseStep [rad] is bisected, at most kContourMaxBisections deep.
// Without it a root within a few rad/s of a box edge aliases by 2 pi and the
// box is skipped.
inline constexpr double kContourMaxPhaseStep = M_PI / 4.0;
inline constexpr int kContourMaxBisections = 12;

// Acceptance threshold: |D(omega)| / Y_char must be below this for a refined
// candidate to be accepted as a genuine mode. Y_char is the sum of the branch
// characteristic admittances A / (rho c) at the junction, so the test is
// dimensionless. A converged Muller root sits at ~1e-12; a non-converged
// candidate next to a root or a pole sits at 1e-2 or above.
inline constexpr double kDispersionTolerance = 1e-6;

// Duplicate-mode guard: two modes with the same m and |Delta f| < this [Hz]
// are considered identical.  Set to half the minimum expected mode spacing.
inline constexpr double kDuplicateGuardHz = 1.0;

// Maximum iterations for Muller root refinement.
inline constexpr int kMullerMaxIter = 50;

// Number of Muller restarts with frequency perturbations when a box contains
// more than one zero (multi-zero case in Argument Principle solver).
inline constexpr int kMullerMultiZoneRestarts = 4;

// Perturbation fraction of box width used for multi-zero Muller restarts.
inline constexpr double kMullerPerturbFraction = 0.25;

// -------------------------------------------------------------
// Annular Duct Geometry (Pure Annular, No Cans)
// -------------------------------------------------------------

// [[deprecated("Use Annulus from geometry.h instead")]]
// AnnularDuctGeometry is now an alias for Annulus from geometry.h
// The Annulus struct has been extended with n_azimuthal_max and radius methods
// for full compatibility with acoustic solvers.
using AnnularDuctGeometry [[deprecated("Use Annulus from geometry.h instead")]] = Annulus;

// -------------------------------------------------------------
// Can-Annular Acoustic Solver Methods
// -------------------------------------------------------------

enum class CanAnnularSolverMethod {
    MagnitudeMinimization,  // Fast, reliable (default)
    ArgumentPrinciple       // Research-grade, finds ALL modes
};

// -------------------------------------------------------------
// Solver Interface
// -------------------------------------------------------------

// Main solver with method selection
std::vector<BlochMode> can_annular_eigenmodes_with_method(
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    double f_max,
    BoundaryCondition bc_can_top,
    CanAnnularSolverMethod method = CanAnnularSolverMethod::MagnitudeMinimization
);

// -------------------------------------------------------------
// Method 1: Magnitude Minimization (WORKING - DEFAULT)
// -------------------------------------------------------------
// Fast, reliable solver using sliding window + magnitude minimization
// Finds modes by scanning for local minima of |D(omega)|
std::vector<BlochMode> solve_magnitude_minimization(
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    double f_max,
    BoundaryCondition bc_can_top
);

// -------------------------------------------------------------
// Method 2: Argument Principle (RESEARCH-GRADE)
// -------------------------------------------------------------
// Robust solver using Argument Principle (Nyquist contour)
// Counts zeros per frequency box by winding number (D is in pole-free
// product form, so the count is zeros only) and refines each counted box with
// Muller from its |D| minima, accepting only points with |D|/Y_char below
// kDispersionTolerance.
//
// References:
//   - Silva et al. (2013): Argument Principle for thermoacoustics
//   - Nicoud et al. (2007): Complex frequency analysis
std::vector<BlochMode> solve_argument_principle(
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    double f_max,
    BoundaryCondition bc_can_top
);

// -------------------------------------------------------------
// Helper: Improved Argument Principle with Pole Handling
// -------------------------------------------------------------
// Winding number of D around the rectangular contour: zeros MINUS poles inside
// (D is pole-free, so this is the zero count). Returns 0 if the
// contour hits a pole, 1 if a contour point hits a root.
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
);

// -------------------------------------------------------------
// Helper: Identify Pole Frequencies
// -------------------------------------------------------------
// Returns the CAN-term poles of the dispersion relation in [f_min, f_max]
// (first 50 only). The plenum term adds poles at n c_plenum / (2 s) that this
// helper cannot see (it takes no c_plenum). Not used by the solvers:
// dispersion_relation_complex is in pole-free product form.
std::vector<double> identify_pole_frequencies(
    const CanAnnularGeometry& geom,
    double c_can,
    double f_min,
    double f_max,
    BoundaryCondition bc_top
);

// -------------------------------------------------------------
// Annular Duct Modes (Pure Annular Geometry, No Cans)
// -------------------------------------------------------------

// Annular duct mode (azimuthal + axial mode numbers)
struct AnnularMode {
    int m_azimuthal;    // Azimuthal mode number (0, 1, 2, ...)
    int n_axial;        // Axial mode number (1, 2, 3, ...)
    double frequency;   // Frequency [Hz]

    // Mode type description
    std::string mode_type() const {
        if (m_azimuthal == 0) {
            return "Axisymmetric (m=0)";
        } else {
            return "Azimuthal (m=" + std::to_string(m_azimuthal) + ")";
        }
    }
};

// Find all annular duct eigenmodes using Argument Principle
//
// Pure annular duct (no cans) with Bloch-Floquet periodic boundary conditions
// Finds modes by solving dispersion relation for annular waveguide
//
// Parameters:
//   geom    : Annulus geometry (from geometry.h)
//   c       : speed of sound [m/s]
//   rho     : density [kg/m^3]
//   f_max   : maximum frequency [Hz]
//   bc_ends : boundary condition at the far end x = L (default: Closed). The
//             end x = 0 is always rigid: Closed gives closed-closed modes
//             f_ax = n c / (2 L); Open gives closed-open (quarter-wave) modes
//             f_ax = (2n - 1) c / (4 L), NOT open-open.
//
// Returns: Vector of AnnularMode sorted by frequency
//
// Method: Argument Principle (Nyquist contour) for robust root finding
// Finds ALL modes in frequency range, no missed roots
// Assumes uniform annular geometry and homogeneous acoustic medium.
std::vector<AnnularMode> annular_duct_eigenmodes(
    const Annulus& geom,
    double c,
    double rho,
    double f_max,
    BoundaryCondition bc_ends = BoundaryCondition::Closed
);

// Analytical annular mode estimate used for validation and lightweight examples
// Uses thin-annulus approximation with combined axial and azimuthal wavenumbers.
// Assumes uniform annular geometry over full axial length.
std::vector<AnnularMode> annular_duct_modes_analytical(
    const Annulus& geom,
    double c,
    double f_max,
    BoundaryCondition bc_ends = BoundaryCondition::Closed
);

}  // namespace combaero
