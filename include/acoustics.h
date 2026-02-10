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

// -------------------------------------------------------------
// Mean Flow Correction (Axial Modes)
// -------------------------------------------------------------
// Convective frequency shift for waves propagating with/against mean flow.
// f± = f₀ / (1 ∓ M)
//
// Valid for: M < 0.3 (linear acoustics)
// Note: For azimuthal modes with swirl, use specialized tools (not covered here)

// Upstream-propagating wave (against flow): f+ = f0 / (1 - M)
double axial_mode_upstream(double f0, double M);

// Downstream-propagating wave (with flow): f- = f0 / (1 + M)
double axial_mode_downstream(double f0, double M);

// Both directions at once: {f_upstream, f_downstream}
std::pair<double, double> axial_mode_split(double f0, double M);

// -------------------------------------------------------------
// Helmholtz Resonator
// -------------------------------------------------------------
// Classic "bottle resonance" frequency.
// f = (c / 2π) * √(A / (V * L_eff))
//
// where L_eff = L_neck + end_correction * d_neck
//
// Parameters:
//   V              : cavity volume [m³]
//   A_neck         : neck cross-sectional area [m²]
//   L_neck         : neck length [m]
//   c              : speed of sound [m/s]
//   end_correction : end correction factor [-] (default: 0.85 for flanged)
//                    Use ~0.6 for unflanged opening
//
// Applications: acoustic dampers, fuel system tuning, side-branch resonators
double helmholtz_frequency(double V, double A_neck, double L_neck, double c,
                           double end_correction = 0.85);

// -------------------------------------------------------------
// Strouhal Number
// -------------------------------------------------------------
// Dimensionless frequency relating oscillation to flow.
// St = f * L / u
//
// Applications:
// - Vortex shedding: St ≈ 0.2 for cylinders
// - Flame dynamics: relates heat release oscillation to convective time
// - Scaling between test rigs and full scale

double strouhal(double f, double L, double u);
double frequency_from_strouhal(double St, double L, double u);

// -------------------------------------------------------------
// Convenience Functions
// -------------------------------------------------------------
// Named helpers for common frequency calculations.
// Trivial but self-documenting.

// Quarter-wave resonator fundamental: f = c / (4L)
// (Open-closed tube, e.g., side-branch damper)
double quarter_wave_frequency(double L, double c);

// Half-wave resonator fundamental: f = c / (2L)
// (Open-open or closed-closed tube)
double half_wave_frequency(double L, double c);

// -------------------------------------------------------------
// Viscothermal Boundary Layers
// -------------------------------------------------------------
// Penetration depths for viscous and thermal effects at walls.
// These determine acoustic losses in resonators and ducts.
//
// Educational screening tools - expect ±50% accuracy.
// Real designs require experimental validation.

// Stokes (viscous) layer thickness: δ_ν = √(2ν/ω)
// nu : kinematic viscosity [m²/s]
// f  : frequency [Hz]
double stokes_layer(double nu, double f);

// Thermal penetration depth: δ_κ = √(2α/ω)
// alpha : thermal diffusivity [m²/s] (use thermal_diffusivity from transport.h)
// f     : frequency [Hz]
double thermal_layer(double alpha, double f);

// Effective viscothermal layer for acoustic losses
// δ_eff = δ_ν + (γ-1)·δ_κ
// Combines viscous and thermal boundary layer effects.
double effective_viscothermal_layer(double delta_nu, double delta_kappa, double gamma);

// -------------------------------------------------------------
// Quality Factor (Screening Estimates)
// -------------------------------------------------------------
// Idealized Q estimates for resonators with viscothermal losses.
// Assumes: smooth walls, no mean flow, linear acoustics, isothermal walls.
//
// These are ORDER OF MAGNITUDE screening tools.
// For design: add 50% margin, validate experimentally.
//
// References:
// - Lieuwen, Unsteady Combustor Physics (lumped element models)
// - Bellucci et al., JEGTP (Helmholtz dampers)
// - Bourquard & Noiray, JSV 2018 (HR/QW comparison)

// Helmholtz resonator quality factor (viscothermal losses only)
// Q ≈ V / (A_neck · δ_eff · correction_factor)
//
// V      : cavity volume [m³]
// A_neck : neck cross-sectional area [m²]
// L_neck : neck length [m]
// nu     : kinematic viscosity [m²/s]
// alpha  : thermal diffusivity [m²/s]
// gamma  : specific heat ratio [-]
// f      : frequency [Hz] (use helmholtz_frequency result)
double helmholtz_Q(double V, double A_neck, double L_neck,
                   double nu, double alpha, double gamma, double f);

// Quarter-wave tube quality factor (viscothermal losses only)
// Q ≈ D / (2 · δ_eff)
//
// L     : tube length [m]
// D     : tube diameter [m]
// nu    : kinematic viscosity [m²/s]
// alpha : thermal diffusivity [m²/s]
// gamma : specific heat ratio [-]
// f     : frequency [Hz]
double tube_Q(double L, double D, double nu, double alpha, double gamma, double f);

// Damping ratio from quality factor: ζ = 1/(2Q)
double damping_ratio(double Q);

// Half-power bandwidth: Δf = f₀/Q
double bandwidth(double f0, double Q);

#endif // ACOUSTICS_H
