#ifndef INCOMPRESSIBLE_H
#define INCOMPRESSIBLE_H

#include "geometry.h"  // Re-export hydraulic_diameter functions

// Incompressible flow equations
//
// These functions assume constant density (incompressible flow).
// Valid for:
//   - Mach number < 0.3 (density change < 5%)
//   - Pressure ratio close to 1 (ΔP/P << 1)
//
// For higher Mach numbers or large pressure ratios, use compressible.h

// -------------------------------------------------------------
// Bernoulli equation
// -------------------------------------------------------------

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

// Pipe pressure drop using Darcy-Weisbach equation.
//
// ΔP = f · (L/D) · (ρ · v² / 2)
//
// Inputs:
//   v   : flow velocity [m/s]
//   L   : pipe length [m]
//   D   : pipe diameter [m]
//   f   : Darcy friction factor [-] (use friction_colebrook etc.)
//   rho : fluid density [kg/m³]
//
// Returns: pressure drop [Pa]
double pipe_dP(double v, double L, double D, double f, double rho);

// Pipe pressure drop from mass flow rate.
// v = ṁ / (ρ · A) where A = π·D²/4
// Returns: pressure drop [Pa]
double pipe_dP_mdot(double mdot, double L, double D, double f, double rho);

// Pipe velocity from mass flow rate.
// v = ṁ / (ρ · π · D² / 4)
// Returns: velocity [m/s]
double pipe_velocity(double mdot, double D, double rho);

// Pipe mass flow from velocity.
// ṁ = ρ · v · π · D² / 4
// Returns: mass flow rate [kg/s]
double pipe_mdot(double v, double D, double rho);

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

#endif // INCOMPRESSIBLE_H
