#ifndef INCOMPRESSIBLE_H
#define INCOMPRESSIBLE_H

#include "geometry.h" // Re-export hydraulic_diameter functions
#include "thermo.h"   // Required for State struct in Solutions

#include <functional>
#include <string>
#include <tuple>
#include <vector>

// Incompressible flow equations
//
// This header provides two layers:
//
//   Layer 1 — scalar primitives (constant density):
//     Bernoulli, orifice, Darcy-Weisbach, zeta/Cd conversions.
//     Inputs: rho, v, P, A, Cd, f (user-supplied fluid properties).
//
//   Layer 2 — thermo-aware high-level API:
//     Takes (T, P, X) and evaluates rho/mu internally.
//     Returns IncompressibleFlowSolution with all relevant quantities.
//     Mirrors the compressible.h API structure.
//
// Valid for:
//   - Mach number < 0.3 (density change < 5%)
//   - Pressure ratio close to 1 (ΔP/P << 1)
//
// For higher Mach numbers or large pressure ratios, use compressible.h

// -------------------------------------------------------------
// Bernoulli equation
// -------------------------------------------------------------

namespace combaero {

// Bernoulli equation: P1 + ½ρv1² + ρgh1 = P2 + ½ρv2² + ρgh2
//
// Solve for P2 given P1, velocities, and elevation change.
//
// Inputs:
//   P1  : upstream pressure [Pa]
//   v1  : upstream velocity [m/s]
//   v2  : downstream velocity [m/s]
//   rho : fluid density [kg/m³]
//   dz  : elevation change z2 - z1 [m] (positive = upward), default 0
//   g   : gravitational acceleration [m/s²], default 9.80665
//
// Returns: downstream pressure P2 [Pa]
double bernoulli_P2(double P1, double v1, double v2, double rho,
                    double dz = 0.0, double g = 9.80665);

// Solve for v2 given pressures and v1.
// Returns: downstream velocity v2 [m/s]
double bernoulli_v2(double P1, double P2, double v1, double rho,
                    double dz = 0.0, double g = 9.80665);

// -------------------------------------------------------------
// Orifice / restriction flow
// -------------------------------------------------------------

// Orifice mass flow rate.
//
// ṁ = Cd · A · √(2 · ρ · ΔP)
//
// Inputs:
//   P1  : upstream pressure [Pa]
//   P2  : downstream pressure [Pa]
//   A   : orifice area [m²]
//   Cd  : discharge coefficient [-] (typically 0.6-0.65 for sharp-edge)
//   rho : fluid density [kg/m³]
//
// Returns: mass flow rate [kg/s]
double orifice_mdot(double P1, double P2, double A, double Cd, double rho);

// Orifice volumetric flow rate.
// Q = Cd · A · √(2 · ΔP / ρ)
// Returns: volumetric flow rate [m³/s]
double orifice_Q(double P1, double P2, double A, double Cd, double rho);

// Ideal orifice velocity (no losses).
// v = √(2 · ΔP / ρ)
// Returns: velocity [m/s]
double orifice_velocity(double P1, double P2, double rho);

// Solve for required orifice area given mass flow.
// A = ṁ / (Cd · √(2 · ρ · ΔP))
// Returns: orifice area [m²]
double orifice_area(double mdot, double P1, double P2, double Cd, double rho);

// Solve for pressure drop given mass flow and area.
// ΔP = (ṁ / (Cd · A))² / (2 · ρ)
// Returns: pressure drop P1 - P2 [Pa]
double orifice_dP(double mdot, double A, double Cd, double rho);

// -------------------------------------------------------------
// Pipe flow (Darcy-Weisbach)
// -------------------------------------------------------------

// Channel pressure drop using Darcy-Weisbach equation.
//
// ΔP = f · (L/D) · (ρ · v² / 2)
//
// Inputs:
//   v   : flow velocity [m/s]
//   L   : length [m]
//   D   : diameter [m]
//   f   : Darcy friction factor [-] (use friction_colebrook etc.)
//   rho : fluid density [kg/m³]
//
// Returns: pressure drop [Pa]
double channel_dP(double v, double L, double D, double f, double rho);


// Channel pressure drop from mass flow rate.
// v = ṁ / (ρ · A) where A = π·D²/4
// Returns: pressure drop [Pa]
double channel_dP_mdot(double mdot, double L, double D, double f, double rho);


// Channel velocity from mass flow rate.
// v = ṁ / (ρ · π · D² / 4)
// Returns: velocity [m/s]
double channel_velocity(double mdot, double D, double rho);


// Channel mass flow from velocity.
// ṁ = ρ · v · π · D² / 4
// Returns: mass flow rate [kg/s]
double channel_mdot(double v, double D, double rho);


// -------------------------------------------------------------
// Hydraulic utilities
// -------------------------------------------------------------

// Dynamic pressure (velocity head).
// q = ½ · ρ · v²
// Returns: dynamic pressure [Pa]
double dynamic_pressure(double v, double rho);

// Velocity from dynamic pressure.
// v = √(2 · q / ρ)
// Returns: velocity [m/s]
double velocity_from_q(double q, double rho);

// Hydraulic diameter functions are in geometry.h (included above)
// - hydraulic_diameter(A, P_wetted)
// - hydraulic_diameter_rect(a, b)
// - hydraulic_diameter_annulus(D_outer, D_inner)

// -------------------------------------------------------------
// Pressure loss coefficient (zeta / K)
// -------------------------------------------------------------
//
// The pressure loss coefficient ζ (also written K) relates pressure drop
// to dynamic pressure:
//
//   ΔP = ζ · ½ρv²
//
// where v is the reference velocity (at the throat / restriction area).
// This form is geometry-independent and applies to any local restriction:
// bends, valves, sudden expansions, orifices, etc.
//
// Relationship to discharge coefficient Cd (throat area as reference):
//
//   ζ = 1 / Cd²     Cd = 1 / √ζ
//
// Note: orifice.h provides K_from_Cd(Cd, beta) / Cd_from_K(K, beta) which
// include the area-ratio (β⁴) correction for pipe orifices where the
// reference velocity is the upstream pipe velocity, not the throat velocity.
// Use the functions below when the reference area IS the throat area.

// Pressure drop from velocity and loss coefficient.
// ΔP = ζ · ½ρv²
// Returns: pressure drop [Pa]
double pressure_loss(double v, double rho, double zeta);

// Velocity from pressure drop and loss coefficient.
// v = √(2·ΔP / (ζ·ρ))
// Returns: velocity [m/s]
double velocity_from_pressure_loss(double dP, double rho, double zeta);

// Loss coefficient from discharge coefficient (throat area as reference).
// ζ = 1 / Cd²
double zeta_from_Cd(double Cd);

// Discharge coefficient from loss coefficient (throat area as reference).
// Cd = 1 / √ζ
double Cd_from_zeta(double zeta);

// -------------------------------------------------------------
// High-level thermo-aware API
// -------------------------------------------------------------
// These functions take (T, P, X) and evaluate density and viscosity
// internally, returning a rich solution struct. They mirror the
// structure of compressible.h (nozzle_flow, fanno_pipe) so that
// the two regimes share a common calling convention.

// Result of incompressible internal flow profile at a single station
struct IncompressibleStation {
  double x = 0.0;   // Position along channel [m]
  double P = 0.0;   // Static pressure [Pa]
  double T = 0.0;   // Static temperature [K]
  double rho = 0.0; // Density [kg/m³]
  double v = 0.0;   // Velocity [m/s]
  double M = 0.0;   // Mach number [-]
  double h = 0.0;   // Specific enthalpy [J/kg]
};

// High-level wrapper for evaluating friction loss while fully respecting
// thermodynamic state (T, P, X). Calculates density internally.
// Function type for optional K-factor loss additions (bends, valves, etc.).
using IncompressibleKLossFn = std::function<double(
    double T, double P, const std::vector<double> &X, double Re)>;

// Result of a thermo-aware incompressible flow calculation.
struct IncompressibleFlowSolution {
  State inlet;       // Inlet thermodynamic state
  State outlet;      // Outlet thermodynamic state
  double mdot = 0.0; // Mass flow rate [kg/s]
  double v = 0.0;    // Bulk velocity [m/s]
  double dP = 0.0;   // Total pressure drop [Pa]
  double Re = 0.0;   // Reynolds number [-]
  double rho = 0.0;  // Inlet bulk density [kg/m³]
  double f = 0.0;    // Darcy friction factor or Cd [-]
  std::vector<IncompressibleStation>
      profile; // Symmetric interpolated spatial array
};

// User-supplied Cd function: f(T, P, X, Re) → Cd [-]
// Re is computed internally from the inlet conditions and velocity.
// Loop-free: uses only inlet state — safe for Newton solvers.
using IncompressibleCdFn = std::function<double(
    double T, double P, const std::vector<double> &X, double Re)>;

// Thermo-aware orifice flow.
//
// Evaluates rho from (T, P, X) internally, then applies the
// incompressible orifice equation: mdot = Cd * A * sqrt(2 * rho * dP).
//
// Inputs:
//   T      : temperature [K]
//   P      : upstream pressure [Pa]
//   X      : mole fractions [-]
//   P_back : downstream pressure [Pa]  (dP = P - P_back)
//   A      : orifice area [m²]
//   Cd     : discharge coefficient [-] (default: 1.0)
//
// Returns: IncompressibleFlowSolution
IncompressibleFlowSolution orifice_flow_thermo(double T, double P,
                                               const std::vector<double> &X,
                                               double P_back, double A,
                                               double Cd = 1.0);

// Overload with user-supplied Cd correlation.
// cd_fn(T, P, X, Re) is called with the inlet Re computed from the
// velocity implied by the fixed pressure difference.
IncompressibleFlowSolution orifice_flow_thermo(double T, double P,
                                               const std::vector<double> &X,
                                               double P_back, double A,
                                               const IncompressibleCdFn &cd_fn);

IncompressibleFlowSolution channel_flow(double T, double P,
                                        const std::vector<double> &X, double u,
                                        double L, double D, double f,
                                        std::size_t n_steps = 100,
                                        bool store_profile = false);


// 2) With roughness-based variable friction factor:
IncompressibleFlowSolution
channel_flow_rough(double T, double P, const std::vector<double> &X, double u,
                   double L, double D, double roughness = 0.0,
                   const std::string &correlation = "haaland",
                   std::size_t n_steps = 100, bool store_profile = false);


// 3) With roughness-based friction + optional K-factor loss:
IncompressibleFlowSolution
channel_flow_rough(double T, double P, const std::vector<double> &X, double u,
                   double L, double D, double roughness,
                   const std::string &correlation,
                   const IncompressibleKLossFn &k_loss_fn,
                   std::size_t n_steps = 10, bool store_profile = false);

inline IncompressibleFlowSolution
pipe_flow_rough(double T, double P, const std::vector<double> &X, double u,
                double L, double D, double roughness = 0.0,
                const std::string &correlation = "haaland",
                std::size_t n_steps = 100, bool store_profile = false) {
  return channel_flow_rough(T, P, X, u, L, D, roughness, correlation, n_steps, store_profile);
}

inline IncompressibleFlowSolution
pipe_flow_rough(double T, double P, const std::vector<double> &X, double u,
                double L, double D, double roughness,
                const std::string &correlation,
                const IncompressibleKLossFn &k_loss_fn,
                std::size_t n_steps = 10, bool store_profile = false) {
  return channel_flow_rough(T, P, X, u, L, D, roughness, correlation, k_loss_fn,
                            n_steps, store_profile);
}

// Composite channel pressure drop (legacy / convenience).
//
// Combines density, viscosity, Reynolds number, friction factor, and
// Darcy-Weisbach into a single call. Returns (dP, Re, f) tuple.
//
// Inputs:
//   T           : temperature [K]
//   P           : pressure [Pa]
//   X           : mole fractions [-]
//   v           : flow velocity [m/s]
//   D           : diameter [m]
//   L           : length [m]
//   roughness   : absolute roughness [m] (default: 0.0)
//   correlation : friction correlation (default: "haaland")
//
// Returns: tuple of (dP [Pa], Re [-], f [-])
std::tuple<double, double, double>
channel_pressure_drop(double T, double P, const std::vector<double> &X, double v,
                      double D, double L, double roughness = 0.0,
                      const std::string &correlation = "haaland");


} // namespace combaero

#endif // INCOMPRESSIBLE_H
