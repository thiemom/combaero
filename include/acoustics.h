#ifndef ACOUSTICS_H
#define ACOUSTICS_H

#include "geometry.h"
#include <string>
#include <vector>

// -------------------------------------------------------------
// Acoustic Mode Analysis
// -------------------------------------------------------------
// Compute natural frequencies of acoustic cavities.
// All frequencies in Hz, dimensions in meters, speed of sound in m/s.
//
// References:
// - Culick, F.E.C. "Unsteady Motions in Combustion Chambers"
// - Dowling, A.P. & Ffowcs Williams, J.E. "Sound and Sources of Sound"
// - Poinsot, T. & Veynante, D. "Theoretical and Numerical Combustion"

// -------------------------------------------------------------
// Boundary Conditions
// -------------------------------------------------------------

enum class BoundaryCondition {
    Open,    // Pressure release (p' = 0), velocity antinode
    Closed   // Rigid wall (u' = 0), pressure antinode
};

// -------------------------------------------------------------
// Acoustic Mode
// -------------------------------------------------------------

struct AcousticMode {
    int n_axial;       // Longitudinal mode number (1, 2, 3, ...)
    int n_azimuthal;   // Tangential/azimuthal mode number (0, 1, 2, ...)
    double frequency;  // Natural frequency [Hz]
    
    // Human-readable mode label
    // Examples: "1L", "2L", "1T", "2T", "1L1T"
    std::string label() const;
};

// -------------------------------------------------------------
// Tube Acoustic Modes
// -------------------------------------------------------------
// Longitudinal (axial) modes in a cylindrical tube.
//
// Frequency formulas:
//   Closed-Closed or Open-Open: f_n = n * c / (2 * L)
//   Open-Closed:                f_n = (2n - 1) * c / (4 * L)
//
// Parameters:
//   tube      : tube geometry (L, D)
//   c         : speed of sound [m/s]
//   upstream  : boundary condition at x = 0
//   downstream: boundary condition at x = L
//   n_max     : maximum mode number to compute (default: 5)
//
// Returns: vector of acoustic modes sorted by frequency

std::vector<AcousticMode> tube_axial_modes(
    const Tube& tube, double c,
    BoundaryCondition upstream, BoundaryCondition downstream,
    int n_max = 5);

// Convenience overload: compute modes from State (uses speed_of_sound)
// Requires: #include "thermo.h"
// std::vector<AcousticMode> tube_axial_modes(
//     const Tube& tube, const State& state,
//     BoundaryCondition upstream, BoundaryCondition downstream,
//     int n_max = 5);

// -------------------------------------------------------------
// Annulus Acoustic Modes
// -------------------------------------------------------------
// Modes in an annular duct (thin annulus approximation).
//
// Longitudinal modes: same as tube
// Azimuthal (tangential) modes: f_m = m * c / (π * D_mean)
//
// Parameters:
//   annulus   : annular geometry (L, D_inner, D_outer)
//   c         : speed of sound [m/s]
//   upstream  : boundary condition at x = 0
//   downstream: boundary condition at x = L
//   n_max     : maximum longitudinal mode number (default: 5)
//   m_max     : maximum azimuthal mode number (default: 4)
//
// Returns: vector of acoustic modes sorted by frequency

std::vector<AcousticMode> annulus_axial_modes(
    const Annulus& annulus, double c,
    BoundaryCondition upstream, BoundaryCondition downstream,
    int n_max = 5);

std::vector<AcousticMode> annulus_azimuthal_modes(
    const Annulus& annulus, double c,
    int m_max = 4);

// Combined axial + azimuthal modes
// f_combined = sqrt(f_axial² + f_azimuthal²)
std::vector<AcousticMode> annulus_modes(
    const Annulus& annulus, double c,
    BoundaryCondition upstream, BoundaryCondition downstream,
    int n_max = 5, int m_max = 4);

// -------------------------------------------------------------
// Utility Functions
// -------------------------------------------------------------

// Find modes within a frequency band
// Useful for: "Are there any modes near my expected instability frequency?"
std::vector<AcousticMode> modes_in_range(
    const std::vector<AcousticMode>& modes,
    double f_min, double f_max);

// Find the mode closest to a target frequency
// Returns: pointer to closest mode, or nullptr if modes is empty
const AcousticMode* closest_mode(
    const std::vector<AcousticMode>& modes,
    double f_target);

// Minimum frequency separation between adjacent modes
// Useful for assessing mode clustering
double min_mode_separation(const std::vector<AcousticMode>& modes);

#endif // ACOUSTICS_H
