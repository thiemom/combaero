#ifndef ACOUSTICS_H
#define ACOUSTICS_H

#include "geometry.h"
#include <complex>
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
// Acoustic properties bundle
// -------------------------------------------------------------

// Bundle of acoustic properties for convenient access
// All properties computed from (f, rho, c, p_rms) in a single call
struct AcousticProperties {
    double wavelength;         // Acoustic wavelength [m]
    double frequency;          // Frequency [Hz]
    double impedance;          // Characteristic acoustic impedance [Pa·s/m]
    double particle_velocity;  // Particle velocity amplitude [m/s]
    double spl;                // Sound pressure level [dB]
};

// Compute all acoustic properties at once
// Parameters:
//   f     : frequency [Hz]
//   rho   : density [kg/m³]
//   c     : speed of sound [m/s]
//   p_rms : RMS pressure amplitude [Pa] (default: 20e-6 Pa = 0 dB reference)
//   p_ref : reference pressure for SPL [Pa] (default: 20e-6 Pa in air)
// Returns: AcousticProperties struct with all properties
AcousticProperties acoustic_properties(
    double f,
    double rho,
    double c,
    double p_rms = 20e-6,
    double p_ref = 20e-6);

// -------------------------------------------------------------
// Transfer Matrix Method (TMM) for Acoustic Networks
// -------------------------------------------------------------

// 2x2 Transfer matrix relating [p; u*S] at upstream and downstream ports
// [p_u]   [T11 T12] [p_d]
// [u_u*S] = [T21 T22] [u_d*S]
//
// Used for cascading acoustic elements (orifices, tubes, resonators)
struct TransferMatrix {
    std::complex<double> T11, T12, T21, T22;

    // Matrix multiplication for cascading elements
    TransferMatrix operator*(const TransferMatrix& other) const;
};

// Convenience bundles for liner and damper calculations
// These reduce long argument lists and make units explicit.
struct AcousticMedium {
    double rho;  // Density [kg/m^3]
    double c;    // Speed of sound [m/s]
};

struct LinerOrificeGeometry {
    double d_orifice;  // Orifice diameter [m]
    double l_orifice;  // Orifice length [m]
    double porosity;   // Open area ratio [-]
    double Cd;         // Discharge coefficient [-]
};

struct LinerFlowState {
    double u_bias;     // Bias flow velocity [m/s]
    double u_grazing;  // Grazing flow velocity [m/s]
};

struct LinerCavity {
    double depth;  // Cavity depth [m]
};

// Reflection and absorption helpers using normalized impedance z = Z/(rho*c)
double absorption_from_impedance_norm(const std::complex<double>& z_norm);

// Single-cavity liner (SDOF) helpers
std::complex<double> liner_sdof_impedance_norm(
    double freq,
    const LinerOrificeGeometry& orifice,
    const LinerCavity& cavity,
    const LinerFlowState& flow,
    const AcousticMedium& medium
);

double liner_sdof_absorption(
    double freq,
    const LinerOrificeGeometry& orifice,
    const LinerCavity& cavity,
    const LinerFlowState& flow,
    const AcousticMedium& medium
);

std::vector<double> sweep_liner_sdof_absorption(
    const std::vector<double>& freqs,
    const LinerOrificeGeometry& orifice,
    const LinerCavity& cavity,
    const LinerFlowState& flow,
    const AcousticMedium& medium
);

// Two-cavity serial liner (2-DOF) helpers
std::complex<double> liner_2dof_serial_impedance_norm(
    double freq,
    const LinerOrificeGeometry& face_orifice,
    const LinerOrificeGeometry& septum_orifice,
    double depth_1,
    double depth_2,
    const LinerFlowState& face_flow,
    const LinerFlowState& septum_flow,
    const AcousticMedium& medium
);

double liner_2dof_serial_absorption(
    double freq,
    const LinerOrificeGeometry& face_orifice,
    const LinerOrificeGeometry& septum_orifice,
    double depth_1,
    double depth_2,
    const LinerFlowState& face_flow,
    const LinerFlowState& septum_flow,
    const AcousticMedium& medium
);

std::vector<double> sweep_liner_2dof_serial_absorption(
    const std::vector<double>& freqs,
    const LinerOrificeGeometry& face_orifice,
    const LinerOrificeGeometry& septum_orifice,
    double depth_1,
    double depth_2,
    const LinerFlowState& face_flow,
    const LinerFlowState& septum_flow,
    const AcousticMedium& medium
);

// Orifice/neck impedance with bias and grazing flow effects
// Combines resistance from bias flow and grazing flow
//
// Parameters:
//   freq       : frequency [Hz]
//   u_bias     : bias flow velocity through orifice [m/s]
//   u_grazing  : grazing flow velocity past orifice [m/s]
//   d_orifice  : orifice diameter [m]
//   l_orifice  : orifice/neck length [m]
//   porosity   : open area ratio A_orifice/A_duct [-]
//   Cd         : discharge coefficient [-]
//   rho        : fluid density [kg/m³]
//   c          : speed of sound [m/s]
//
// Returns: Complex impedance Z/ρc (normalized) [-]
//
// References:
//   - Bourquard & Noiray (2019): Resistance with flow
//   - Rogers & Marble (1956): Bias flow resistance
//   - Howe (1979): Grazing flow effects
//
// Valid: M_bias < 0.3, M_grazing < 0.3, porosity = 0.05-0.3
// Accuracy: +/-20-30% (flow-acoustic interaction is complex)
std::complex<double> orifice_impedance_with_flow(
    double freq,
    double u_bias,
    double u_grazing,
    double d_orifice,
    double l_orifice,
    double porosity,
    double Cd,
    double rho,
    double c
);

// Quarter-wave resonator transfer matrix with flow effects
// Distributed parameter model for side-branch tube resonator
//
// Three distinct areas:
//   A_duct     : main duct cross-sectional area [m²]
//   A_orifice  : neck/orifice area (= porosity * A_duct) [m²]
//   A_tube     : resonator tube cross-sectional area [m²]
//
// Parameters:
//   freq       : frequency [Hz]
//   L_tube     : resonator tube length [m]
//   A_duct     : main duct area [m²]
//   A_tube     : resonator tube area [m²]
//   d_orifice  : neck/orifice diameter [m]
//   l_orifice  : neck length [m]
//   porosity   : open area ratio A_orifice/A_duct [-]
//   Cd         : discharge coefficient [-]
//   u_bias     : bias flow through orifice [m/s]
//   u_grazing  : grazing flow in main duct [m/s]
//   rho        : fluid density [kg/m³]
//   c          : speed of sound [m/s]
//
// Returns: 2x2 TransferMatrix for the side-branch resonator
//
// References:
//   - Munjal (1987): Acoustics of Ducts and Mufflers
//   - Dowling & Ffowcs Williams (1983): Sound and Sources of Sound
//
// Valid: kL < π (below cutoff), M < 0.3
// Accuracy: +/-15-20%
TransferMatrix quarter_wave_resonator_tmm(
    double freq,
    double L_tube,
    double A_duct,
    double A_tube,
    double d_orifice,
    double l_orifice,
    double porosity,
    double Cd,
    double u_bias,
    double u_grazing,
    double rho,
    double c
);

// Strouhal-based whistling risk assessment for orifices
// Vortex shedding can lock-in with acoustic modes causing whistling
//
// Parameters:
//   freq      : acoustic frequency [Hz]
//   u_bias    : flow velocity through orifice [m/s]
//   d_orifice : orifice diameter [m]
//
// Returns: true if in critical Strouhal range (0.2 < St < 0.5)
//
// References:
//   - Rockwell & Naudascher (1979): Self-sustained oscillations
//   - Bruggeman et al. (1991): Flow-induced pulsations
//
// Critical Strouhal range: 0.2-0.5 (vortex lock-in regime)
bool is_whistling_risk(
    double freq,
    double u_bias,
    double d_orifice
);

// -------------------------------------------------------------
// Can-Annular Combustor Acoustics (Bloch-Floquet Theory)
// -------------------------------------------------------------

// Geometry of a can-annular combustor system
// N cans coupled by an annular plenum.
//
// Uniform-section assumption (current model):
// - Each can is represented by a single uniform acoustic segment
//   with constant area_can and length_can.
// - The annular plenum is represented by a uniform ring with
//   constant mean radius radius_plenum and area_plenum.
// - No axial variation/taper/segmentation is represented in this struct.
struct CanAnnularGeometry {
    int n_cans;            // Total number of cans (N)
    double length_can;     // Length of one can [m]
    double area_can;       // Cross-sectional area of one can [m²]
    double radius_plenum;  // Mean radius of the annular plenum [m]
    double area_plenum;    // Cross-sectional area of the plenum [m²]
};

// Adapter: Convert CanAnnularFlowGeometry (flow/geometry) to CanAnnularGeometry (acoustic solver)
//
// This function bridges the two geometry representations:
// - CanAnnularFlowGeometry (from geometry.h): Flow-focused with primary/transition/annular sections
// - CanAnnularGeometry (above): Acoustic-focused with cans + annular plenum
//
// Parameters:
//   flow_geom : CanAnnularFlowGeometry with flow geometry parameters
//   n_cans    : Number of cans (required for acoustic model)
//   L_can     : Can length [m] (acoustic model uses single section per can)
//   D_can     : Can diameter [m] (for circular can cross-section)
//
// Returns: CanAnnularGeometry for acoustic solvers
//
// Notes:
// - Plenum radius derived from flow_geom.D_mean() / 2
// - Plenum area from flow_geom.area()
// - Can area from π*(D_can/2)²
//
// Example:
//   auto flow_geom = CanAnnularFlowGeometry{...};
//   auto acoustic_geom = to_acoustic_geometry(flow_geom, 24, 0.5, 0.1);
//   auto modes = can_annular_eigenmodes(acoustic_geom, c_can, c_plenum, ...);
CanAnnularGeometry to_acoustic_geometry(
    const CanAnnularFlowGeometry& flow_geom,
    int n_cans,
    double L_can,
    double D_can
);

// Bloch mode (eigenmode of periodic can-annular system)
// Uses Bloch-Floquet theory to reduce N-can problem to single sector
struct BlochMode {
    int m_azimuthal;       // Azimuthal mode number (0 to N/2)
    double frequency;      // Eigenfrequency [Hz]
    int n_cans;            // Number of cans (needed for symmetry classification)

    // Mode symmetry classification
    // m=0: "Push-Push" (all cans in phase)
    // m=N/2: "Push-Pull" (adjacent cans 180° out of phase)
    // 0 < m < N/2: "Spinning/Standing wave"
    std::string symmetry_type() const;
};

// Solve for passive acoustic eigenmodes of can-annular system
// Uses Argument Principle Method (Nyquist contour) for robust root finding
//
// Finds ALL modes in frequency range by counting zeros of dispersion relation:
//   D(ω) = Y_can(ω) + Y_annulus(ω, m) = 0
//
// Modeling assumption: geometry is uniform within each region defined by
// CanAnnularGeometry (single can section + single annulus section).
//
// Parameters:
//   geom       : can-annular geometry
//   c_can      : speed of sound in can [m/s]
//   c_plenum   : speed of sound in plenum [m/s] (often different!)
//   rho_can    : density in can [kg/m³]
//   rho_plenum : density in plenum [kg/m³]
//   f_max      : maximum frequency to search [Hz] (default: 2000 Hz)
//   bc_can_top : boundary condition at can top (default: Closed)
//
// Returns: Vector of BlochMode sorted by frequency
//
// References:
//   - Noiray et al. (2011): "A unified framework for nonlinear combustion instability analysis"
//   - Evesque & Polifke (2005): "Low-order acoustic modelling for annular combustors"
//   - Silva et al. (2013): "Combining Helmholtz solver with flame describing function"
//
// Method: Argument Principle (counts zeros inside contour in complex plane)
// Accuracy: Finds all modes, no missed roots
// Valid: Passive acoustics (no flame response)
std::vector<BlochMode> can_annular_eigenmodes(
    const CanAnnularGeometry& geom,
    double c_can,
    double c_plenum,
    double rho_can,
    double rho_plenum,
    double f_max = 2000.0,
    BoundaryCondition bc_can_top = BoundaryCondition::Closed
);

// -------------------------------------------------------------
// Utility functions
// -------------------------------------------------------------

// Acoustic wavelength from frequency
// λ = c / f
double wavelength(double f, double c);

// Frequency from wavelength
// f = c / λ
double frequency_from_wavelength(double lambda, double c);

// Characteristic acoustic impedance
// Z = ρ · c
double acoustic_impedance(double rho, double c);

// Sound pressure level in dB
// SPL = 20 · log₁₀(p_rms / p_ref)
// Default p_ref = 20 μPa (standard reference in air)
double sound_pressure_level(double p_rms, double p_ref = 20e-6);

// Particle velocity from pressure amplitude
// u = p / (ρ · c)
double particle_velocity(double p, double rho, double c);

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
