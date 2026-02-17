#pragma once

#include "acoustics.h"
#include <vector>
#include <complex>

namespace combaero {

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

}  // namespace combaero
