#pragma once

#include "acoustics.h"
#include "math_constants.h"  // MSVC compatibility for M_PI
#include <vector>
#include <complex>

namespace combaero {

// -------------------------------------------------------------
// Annular Duct Geometry (Pure Annular, No Cans)
// -------------------------------------------------------------

struct AnnularDuctGeometry {
    double length;          // Axial length [m]
    double radius_inner;    // Inner radius [m]
    double radius_outer;    // Outer radius [m]
    int n_azimuthal_max;    // Maximum azimuthal mode number to search

    // Derived: cross-sectional area
    double area() const {
        return M_PI * (radius_outer * radius_outer - radius_inner * radius_inner);
    }
};

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
// Finds modes by scanning for local minima of |D(ω)|
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
// Guarantees finding ALL modes by counting zeros in complex plane
// Handles poles explicitly, works for complex frequencies
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
// Counts zeros inside rectangular contour with adaptive pole avoidance
// Returns: Number of zeros (or -1 if contour is unreliable)
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
// Returns frequencies where dispersion relation has poles
// Used to shift contours away from singularities
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
//   geom    : annular duct geometry
//   c       : speed of sound [m/s]
//   rho     : density [kg/m³]
//   f_max   : maximum frequency [Hz]
//   bc_ends : boundary condition at duct ends (default: Closed)
//
// Returns: Vector of AnnularMode sorted by frequency
//
// Method: Argument Principle (Nyquist contour) for robust root finding
// Finds ALL modes in frequency range, no missed roots
std::vector<AnnularMode> annular_duct_eigenmodes(
    const AnnularDuctGeometry& geom,
    double c,
    double rho,
    double f_max,
    BoundaryCondition bc_ends = BoundaryCondition::Closed
);

// Analytical annular mode estimate used for validation and lightweight examples
// Uses thin-annulus approximation with combined axial and azimuthal wavenumbers.
std::vector<AnnularMode> annular_duct_modes_analytical(
    const AnnularDuctGeometry& geom,
    double c,
    double f_max,
    BoundaryCondition bc_ends = BoundaryCondition::Closed
);

}  // namespace combaero
