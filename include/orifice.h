#ifndef ORIFICE_H
#define ORIFICE_H

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

// -------------------------------------------------------------
// Orifice Discharge Coefficient (Cd) Correlations
// -------------------------------------------------------------
//
// This module provides Cd correlations for various orifice geometries:
//   - Sharp thin-plate orifices (ISO 5167, Reader-Harris/Gallagher)
//   - Thick-plate orifices (t/d correction per Idelchik/Bohl)
//   - Rounded-entry orifices (r/d based, Idelchik)
//   - User-defined (tabulated or custom function)
//
// The discharge coefficient Cd relates actual to ideal flow:
//   mdot_actual = Cd * mdot_ideal
//   mdot_ideal  = A * sqrt(2 * rho * dP)
//
// Compressible flow note:
//   For sharp-edged orifices, Cd is weakly dependent on Mach number and
//   remains dominated by geometry and Reynolds number. These incompressible
//   correlations provide adequate Cd values into the choked-flow regime when
//   combined with compressible mass-flow relations (see compressible.h).
//   Pressure-ratio corrections are typically < 5% and may be added via
//   empirical factors (e.g., Spink) or user-supplied test data.
//
// References:
// - ISO 5167-2:2003 - Orifice plates
// - Reader-Harris & Gallagher (1998) - NEL/ASME correlation
// - Idelchik, I.E. - Handbook of Hydraulic Resistance (3rd ed.)
// - Bohl, W. - Technische Stroemungslehre
// - Spink, L.K. - Principles and Practice of Flow Meter Engineering

// -------------------------------------------------------------
// Orifice geometry
// -------------------------------------------------------------

enum class OrificeType {
    SharpThinPlate,   // ISO 5167-type sharp-edged thin plate
    ThickPlate,       // Finite thickness plate (t/d > 0)
    RoundedEntry,     // Rounded inlet edge (r/d > 0)
    Conical,          // Conical (beveled) inlet
    QuarterCircle,    // Quarter-circle (quadrant) edge
    UserDefined       // Custom correlation
};

struct OrificeGeometry {
    double d = 0.0;       // Orifice bore diameter [m]
    double D = 0.0;       // Pipe diameter [m]
    double t = 0.0;       // Plate thickness [m] (for thick plate)
    double r = 0.0;       // Inlet edge radius [m] (for rounded entry)
    double bevel = 0.0;   // Bevel angle [rad] (for conical)

    // Derived quantities
    double beta() const;          // Diameter ratio d/D [-]
    double area() const;          // Orifice area [m^2]
    double t_over_d() const;      // Thickness ratio t/d [-]
    double r_over_d() const;      // Radius ratio r/d [-]

    // Validation
    bool is_valid() const;
};

// -------------------------------------------------------------
// Flow state at orifice
// -------------------------------------------------------------

struct OrificeState {
    double Re_D = 0.0;    // Pipe Reynolds number (based on D) [-]
    double dP = 0.0;      // Differential pressure across orifice [Pa]
    double rho = 0.0;     // Fluid density [kg/m³]
    double mu = 0.0;      // Dynamic viscosity [Pa·s]

    // Derived quantities
    double Re_d(double beta) const;  // Orifice Reynolds number (based on d)
};

// -------------------------------------------------------------
// Correlation identifiers
// -------------------------------------------------------------

enum class CdCorrelation {
    // Sharp thin-plate correlations
    ReaderHarrisGallagher,  // ISO 5167-2 / ASME MFC-3M (most accurate)
    Stolz,                  // ISO 5167:1980 (older, simpler)
    Miller,                 // Miller (1996) - simplified

    // Thick-plate corrections
    IdelchikThick,          // Idelchik thick-plate correction
    BohlThick,              // Bohl thick-plate correction

    // Rounded-entry correlations
    IdelchikRounded,        // Idelchik rounded-entry
    BohlRounded,            // Bohl rounded-entry

    // Special
    Constant,               // Fixed Cd value (for testing/simple cases)
    UserFunction            // User-provided function
};

// -------------------------------------------------------------
// Cd computation - free functions (simple interface)
// -------------------------------------------------------------

// Sharp thin-plate orifice (ISO 5167-2, Reader-Harris/Gallagher)
// Valid for: 0.1 <= beta <= 0.75, Re_D >= 5000, D >= 50mm
double Cd_sharp_thin_plate(const OrificeGeometry& geom, const OrificeState& state);

// Thick-plate orifice (sharp edges, finite thickness)
// Applies thickness correction to thin-plate Cd
// Valid for: 0 < t/d < ~3
double Cd_thick_plate(const OrificeGeometry& geom, const OrificeState& state);

// Rounded-entry orifice
// Valid for: 0 < r/d <= 0.2 (typical)
double Cd_rounded_entry(const OrificeGeometry& geom, const OrificeState& state);

// Convenience: auto-select correlation based on geometry
double Cd(const OrificeGeometry& geom, const OrificeState& state);

// -------------------------------------------------------------
// Individual correlation implementations
// -------------------------------------------------------------

namespace orifice {

// Reader-Harris/Gallagher (1998) - ISO 5167-2
// The standard correlation for sharp-edged orifices
double Cd_ReaderHarrisGallagher(double beta, double Re_D, double D);

// Stolz (1978) - older ISO 5167 correlation
double Cd_Stolz(double beta, double Re_D);

// Miller (1996) - simplified correlation
double Cd_Miller(double beta, double Re_D);

// Thickness correction factor (multiplies thin-plate Cd)
// Based on Idelchik: accounts for flow reattachment in thick plates
double thickness_correction(double t_over_d, double beta);

// Rounded-entry Cd (Idelchik-based)
// For well-rounded entries, Cd approaches 1.0
double Cd_rounded(double r_over_d, double beta, double Re_D);

// Loss coefficient K from Cd: K = (1/Cd^2 - 1) * (1 - beta^4)
double K_from_Cd(double Cd, double beta);

// Cd from loss coefficient K
double Cd_from_K(double K, double beta);

} // namespace orifice

// -------------------------------------------------------------
// Orifice correlation class (for polymorphic use)
// -------------------------------------------------------------

class OrificeCorrelationBase {
public:
    virtual ~OrificeCorrelationBase() = default;
    virtual double Cd(const OrificeGeometry& geom, const OrificeState& state) const = 0;
    virtual std::string name() const = 0;
};

// Factory function to create correlation objects
// Returns nullptr for UserFunction (use make_user_correlation instead)
std::unique_ptr<OrificeCorrelationBase> make_correlation(CdCorrelation id);

// User-defined correlation from function
using CdFunction = std::function<double(const OrificeGeometry&, const OrificeState&)>;
std::unique_ptr<OrificeCorrelationBase> make_user_correlation(
    CdFunction fn,
    const std::string& name = "UserDefined");

// User-defined correlation from tabulated data (beta, Re_D, Cd)
// Interpolates linearly in beta and log(Re_D)
std::unique_ptr<OrificeCorrelationBase> make_tabulated_correlation(
    const std::vector<double>& beta_values,
    const std::vector<double>& Re_values,
    const std::vector<std::vector<double>>& Cd_table,
    const std::string& name = "Tabulated");

// -------------------------------------------------------------
// Orifice flow calculations (uses incompressible.h internally)
// -------------------------------------------------------------

// Mass flow through orifice given Cd
// mdot = Cd * A * sqrt(2 * rho * dP)
double orifice_mdot(const OrificeGeometry& geom, double Cd, double dP, double rho);

// Pressure drop for given mass flow
// dP = (mdot / (Cd * A))^2 / (2 * rho)
double orifice_dP(const OrificeGeometry& geom, double Cd, double mdot, double rho);

// Solve for Cd given measured mdot and dP
// Cd = mdot / (A * sqrt(2 * rho * dP))
double orifice_Cd_from_measurement(const OrificeGeometry& geom,
                                    double mdot, double dP, double rho);

// -------------------------------------------------------------
// Iterative solver for Cd-Re coupling
// -------------------------------------------------------------

// Solve for mass flow rate accounting for Cd-Re_D dependency.
//
// This function iteratively solves the coupled system:
//   mdot = ε · Cd(Re_D) · A · √(2 · ρ · ΔP)
//   Re_D = 4 · mdot / (π · D · μ)
//
// The iteration continues until mdot converges within the specified tolerance.
//
// Algorithm:
//   1. Start with initial Cd guess (typically 0.61 for sharp orifices)
//   2. Calculate mdot using current Cd and expansibility factor ε
//   3. Update Re_D based on new mdot
//   4. Recalculate Cd using updated Re_D
//   5. Repeat until |mdot_new - mdot_old| < tol * mdot_old
//
// Parameters:
//   geom        : Orifice geometry
//   dP          : Differential pressure [Pa]
//   rho         : Upstream density [kg/m³]
//   mu          : Dynamic viscosity [Pa·s]
//   P_upstream  : Absolute upstream pressure [Pa] (for compressibility)
//   kappa       : Isentropic exponent cp/cv [-] (0 = incompressible)
//   correlation : Cd correlation to use (default: ReaderHarrisGallagher)
//   tol         : Relative convergence tolerance (default: 1e-6)
//   max_iter    : Maximum iterations (default: 20)
//
// Returns: Converged mass flow rate [kg/s]
//
// Throws: std::runtime_error if iteration fails to converge
//
// Note: For incompressible flow, set kappa = 0 (expansibility ε = 1)
//       For typical applications, convergence occurs in 3-5 iterations
double solve_orifice_mdot(
    const OrificeGeometry& geom,
    double dP,
    double rho,
    double mu,
    double P_upstream = 101325.0,
    double kappa = 0.0,
    CdCorrelation correlation = CdCorrelation::ReaderHarrisGallagher,
    double tol = 1e-6,
    int max_iter = 20);

// -------------------------------------------------------------
// Compressible flow correction
// -------------------------------------------------------------

// Expansibility factor for compressible gas flow through orifices.
//
// The expansibility factor ε accounts for gas expansion as it accelerates
// through the orifice. For incompressible flow, ε = 1.0.
//
// Formula from ISO 5167-2:2003, Section 5.3.2.2:
//   ε = 1 - (0.351 + 0.256·β⁴ + 0.93·β⁸) · [1 - (1 - τ)^(1/κ)]
//
// where:
//   τ = ΔP / P_upstream  (pressure ratio)
//   β = d/D              (diameter ratio)
//   κ = cp/cv            (isentropic exponent)
//
// Usage in compressible mass flow:
//   mdot = ε · Cd · A · √(2 · ρ_upstream · ΔP)
//
// Valid for:
//   - 0.1 ≤ β ≤ 0.75
//   - τ ≤ 0.25 (ΔP/P ≤ 25%)
//   - Ideal gas behavior
//
// Parameters:
//   beta        : Diameter ratio d/D [-]
//   dP          : Differential pressure [Pa]
//   P_upstream  : Absolute upstream static pressure [Pa]
//   kappa       : Isentropic exponent cp/cv [-]
//
// Returns: Expansibility factor ε [-] (dimensionless, 0 < ε ≤ 1)
//
// Reference: ISO 5167-2:2003, Measurement of fluid flow by means of
//            pressure differential devices inserted in circular cross-section
//            conduits running full
double expansibility_factor(double beta, double dP, double P_upstream, double kappa);

#endif // ORIFICE_H
