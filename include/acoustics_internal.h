#pragma once

// Internal helpers shared between acoustics.cpp and can_annular_solvers.cpp.
// Not part of the public API — do not include from python/ or external code.

#include "acoustics.h"
#include <complex>
#include <vector>

namespace combaero {

// Dispersion relation D(ω, m) = Y_can + Y_annulus evaluated at complex ω.
// Used by Argument Principle contour integration.
std::complex<double> dispersion_relation_complex(
    std::complex<double> omega,
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    BoundaryCondition bc_top
);

// Real-ω overload (delegates to complex version with zero imaginary part).
std::complex<double> dispersion_relation(
    double omega,
    int m,
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    BoundaryCondition bc_top
);

// Scan [f_min, f_max] for local minima of |D(ω)| and return up to
// n_zeros_expected frequency guesses suitable for Muller refinement.
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
);

// Muller's method: refine a frequency guess to a root of D(ω, m) = 0.
// Returns refined frequency [Hz].
double refine_root_muller(
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

}  // namespace combaero
