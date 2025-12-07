#ifndef COMPRESSIBLE_H
#define COMPRESSIBLE_H

#include "state.h"
#include <vector>
#include <cstddef>

// -------------------------------------------------------------
// Compressible flow for ideal gas with variable cp(T)
// -------------------------------------------------------------
// Uses ideal gas equation of state (PV = nRT) with temperature-dependent
// thermodynamic properties (cp, h, s) from NASA polynomial fits.
// This is sometimes called "thermally perfect" gas modeling.
//
// Includes:
// - Isentropic nozzle flow (solved numerically since gamma varies with T)
// - Fanno flow (adiabatic pipe flow with friction)

// -------------------------------------------------------------
// Compressible flow solution
// -------------------------------------------------------------

// Result of a compressible flow calculation
struct CompressibleFlowSolution {
    State stagnation;    // Stagnation state (T0, P0, h0, s0, ...)
    State outlet;        // Outlet static state (T, P, rho, h, s, a)
    double v = 0.0;      // Outlet velocity [m/s]
    double M = 0.0;      // Outlet Mach number [-]
    double mdot = 0.0;   // Mass flow rate [kg/s]
    bool choked = false; // True if flow is choked (M = 1 at throat)
};

// -------------------------------------------------------------
// Forward problem: given geometry and pressures, find mass flow
// -------------------------------------------------------------

// Isentropic nozzle flow: given stagnation conditions, back pressure, and
// effective flow area (A * Cd), compute mass flow and outlet state.
//
// T0      : stagnation temperature [K]
// P0      : stagnation pressure [Pa]
// P_back  : back pressure [Pa]
// A_eff   : effective flow area [m²] (= A * Cd)
// X       : mole fractions
// tol     : convergence tolerance (default: 1e-8)
// max_iter: maximum iterations (default: 50)
//
// Returns CompressibleFlowSolution with outlet state at throat/exit.
// If P_back <= P_critical, flow is choked and outlet is at sonic conditions.
CompressibleFlowSolution nozzle_flow(
    double T0, double P0, double P_back, double A_eff,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Inverse problems: solve for unknown given mass flow
// -------------------------------------------------------------

// Find effective area for given mass flow rate.
// Throws if mdot_target exceeds choked mass flow.
double solve_A_eff_from_mdot(
    double T0, double P0, double P_back, double mdot_target,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Find back pressure for given mass flow rate.
// Returns the subsonic solution (P_back > P_critical).
// Throws if mdot_target exceeds choked mass flow.
double solve_P_back_from_mdot(
    double T0, double P0, double A_eff, double mdot_target,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Find stagnation pressure for given mass flow rate.
// Throws if no solution exists.
double solve_P0_from_mdot(
    double T0, double P_back, double A_eff, double mdot_target,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Utility functions
// -------------------------------------------------------------

// Compute critical (sonic) pressure ratio P*/P0 for given stagnation state.
// Uses isentropic expansion to find where M = 1.
double critical_pressure_ratio(
    double T0, double P0, const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Compute Mach number from pressure ratio P/P0 (isentropic).
double mach_from_pressure_ratio(
    double T0, double P0, double P,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Compute mass flux G = rho * v [kg/(m²·s)] at given pressure for isentropic
// expansion from stagnation conditions.
double mass_flux_isentropic(
    double T0, double P0, double P,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Fanno flow: adiabatic pipe flow with friction
// -------------------------------------------------------------
// Solves compressible pipe flow with wall friction, conserving:
//   - Mass: ρ * u * A = const = mdot
//   - Energy: h + u²/2 = h0 (adiabatic, no work)
//   - Momentum: dp/dx = -f/(2D) * ρ * u² (friction loss)
//
// Integration via RK4 on p(x), solving for T from energy at each step.

// Result of Fanno flow calculation at a single station
struct FannoStation {
    double x = 0.0;      // Position along pipe [m]
    double P = 0.0;      // Static pressure [Pa]
    double T = 0.0;      // Static temperature [K]
    double rho = 0.0;    // Density [kg/m³]
    double u = 0.0;      // Velocity [m/s]
    double M = 0.0;      // Mach number [-]
    double h = 0.0;      // Specific enthalpy [J/kg]
    double s = 0.0;      // Specific entropy [J/(kg·K)]
};

// Result of Fanno flow pipe segment calculation
struct FannoSolution {
    State inlet;                        // Inlet thermodynamic state
    State outlet;                       // Outlet thermodynamic state
    double mdot = 0.0;                  // Mass flow rate [kg/s]
    double h0 = 0.0;                    // Stagnation enthalpy [J/kg]
    double L = 0.0;                     // Pipe length [m]
    double D = 0.0;                     // Pipe diameter [m]
    double f = 0.0;                     // Darcy friction factor [-]
    bool choked = false;                // True if flow reached M=1
    double L_choke = 0.0;               // Length to choking (if choked) [m]
    std::vector<FannoStation> profile;  // Axial profile (optional)
};

// Solve Fanno flow through a pipe segment.
//
// Inputs:
//   T_in    : Inlet static temperature [K]
//   P_in    : Inlet static pressure [Pa]
//   u_in    : Inlet velocity [m/s]
//   L       : Pipe length [m]
//   D       : Pipe diameter [m]
//   f       : Darcy friction factor [-] (constant along pipe)
//   X       : Mole fractions [-]
//   n_steps : Number of integration steps (default: 100)
//   store_profile : If true, store axial profile in solution
//
// Returns FannoSolution with outlet conditions.
// If flow chokes (M→1) before reaching L, integration stops and choked=true.
FannoSolution fanno_pipe(
    double T_in, double P_in, double u_in,
    double L, double D, double f,
    const std::vector<double>& X,
    std::size_t n_steps = 100,
    bool store_profile = false);

// Convenience: solve given inlet State and velocity
FannoSolution fanno_pipe(
    const State& inlet, double u_in,
    double L, double D, double f,
    std::size_t n_steps = 100,
    bool store_profile = false);

// Compute maximum pipe length before choking (L*) for given inlet conditions.
// Returns the length at which M would reach 1.0.
double fanno_max_length(
    double T_in, double P_in, double u_in,
    double D, double f,
    const std::vector<double>& X,
    double tol = 1e-6, std::size_t max_iter = 100);

#endif // COMPRESSIBLE_H
