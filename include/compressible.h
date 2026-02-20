#ifndef COMPRESSIBLE_H
#define COMPRESSIBLE_H

#include "state.h"
#include <vector>
#include <cstddef>
#include <functional>
#include <utility>

// -------------------------------------------------------------
// Compressible flow for ideal gas with variable cp(T)
// -------------------------------------------------------------
// Uses ideal gas equation of state (PV = nRT) with temperature-dependent
// thermodynamic properties (cp, h, s) from NASA polynomial fits.
// This is sometimes called "thermally perfect" gas modeling.
//
// Includes:
// - Isentropic nozzle flow (solved numerically since gamma varies with T)
// - Quasi-1D nozzle flow with variable area A(x)
// - Fanno flow (adiabatic pipe flow with friction)

// -------------------------------------------------------------
// Compressible flow solution
// -------------------------------------------------------------

// Result of a compressible flow calculation
struct CompressibleFlowSolution {
    // Stagnation state at (T0, P0). This is the one place in CombAero where
    // a State struct holds total rather than static conditions.
    // stagnation.h() returns h(T0) = stagnation enthalpy [J/mol].
    State stagnation;    // Stagnation state (T0, P0, h0, s0, ...)
    State outlet;        // Outlet static state (T, P, rho, h, s, a)
                         // For heat transfer at high Mach: use T_adiabatic_wall()
                         // from stagnation.h, not outlet.T, as the hot-side
                         // driving temperature in htc_pipe() / cooled_wall_heat_flux()
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
// Quasi-1D nozzle flow with variable area A(x)
// -------------------------------------------------------------
// Solves quasi-1D compressible flow through a nozzle with area variation:
//   - Mass: mdot = ρ * u * A(x) = const
//   - Energy: h + u²/2 = h0 (adiabatic)
//   - Momentum: dp/dx + ρ*u*du/dx = 0 (frictionless) or with loss model
//
// Marches along x, solving 2x2 nonlinear system (T, P) at each station
// using Newton iteration with analytic derivatives from thermo backend.

// Function type for area distribution A(x)
// x: axial position [m], returns area [m²]
using AreaFunction = std::function<double(double)>;

// Station data for quasi-1D nozzle flow
struct NozzleStation {
    double x = 0.0;      // Axial position [m]
    double A = 0.0;      // Area [m²]
    double P = 0.0;      // Static pressure [Pa]
    double T = 0.0;      // Static temperature [K]
    double rho = 0.0;    // Density [kg/m³]
    double u = 0.0;      // Velocity [m/s]
    double M = 0.0;      // Mach number [-]
    double h = 0.0;      // Specific enthalpy [J/kg]
};

// Result of quasi-1D nozzle calculation
struct NozzleSolution {
    State inlet;                         // Inlet thermodynamic state
    State outlet;                        // Outlet thermodynamic state
    double mdot = 0.0;                   // Mass flow rate [kg/s]
    double h0 = 0.0;                     // Stagnation enthalpy [J/kg]
    double T0 = 0.0;                     // Stagnation temperature [K]
    double P0 = 0.0;                     // Stagnation pressure [Pa]
    bool choked = false;                 // True if throat is sonic
    double x_throat = 0.0;               // Throat position [m]
    double A_throat = 0.0;               // Throat area [m²]
    std::vector<NozzleStation> profile;  // Axial profile
};

// Solve quasi-1D nozzle flow with given area distribution.
//
// Inputs:
//   T0, P0    : Stagnation conditions [K, Pa]
//   P_exit    : Exit static pressure [Pa] (determines if subsonic/supersonic)
//   area_func : Function A(x) returning area [m²] at position x [m]
//   x_start   : Start position [m]
//   x_end     : End position [m]
//   X         : Mole fractions [-]
//   n_stations: Number of axial stations (default: 100)
//
// Returns NozzleSolution with axial profile.
// For subsonic flow throughout, P decreases then increases.
// For choked flow, throat is sonic (M=1) at minimum area.
NozzleSolution nozzle_quasi1d(
    double T0, double P0, double P_exit,
    const AreaFunction& area_func,
    double x_start, double x_end,
    const std::vector<double>& X,
    std::size_t n_stations = 100,
    double tol = 1e-8, std::size_t max_iter = 50);

// Convenience: solve with area given as vector of (x, A) pairs
// Linearly interpolates between points.
NozzleSolution nozzle_quasi1d(
    double T0, double P0, double P_exit,
    const std::vector<std::pair<double, double>>& area_profile,  // (x, A) pairs
    const std::vector<double>& X,
    std::size_t n_stations = 100,
    double tol = 1e-8, std::size_t max_iter = 50);

// Convenience: solve with polynomial area distribution
// A(x) = a0 + a1*x + a2*x² + ...
NozzleSolution nozzle_quasi1d_poly(
    double T0, double P0, double P_exit,
    const std::vector<double>& area_coeffs,  // Polynomial coefficients [a0, a1, a2, ...]
    double x_start, double x_end,
    const std::vector<double>& X,
    std::size_t n_stations = 100,
    double tol = 1e-8, std::size_t max_iter = 50);

// Convenience: solve for converging-diverging nozzle with given geometry
// A(x) = A_inlet at x=0, A_throat at x=x_throat, A_exit at x=x_exit
// Uses smooth cosine interpolation between sections.
NozzleSolution nozzle_cd(
    double T0, double P0, double P_exit,
    double A_inlet, double A_throat, double A_exit,
    double x_throat, double x_exit,
    const std::vector<double>& X,
    std::size_t n_stations = 100,
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
    double x   = 0.0;    // Position along pipe [m]
    double P   = 0.0;    // Static pressure [Pa]
    double T   = 0.0;    // Static temperature [K]
    double rho = 0.0;    // Density [kg/m³]
    double u   = 0.0;    // Velocity [m/s]
    double M   = 0.0;    // Mach number [-]
    double h   = 0.0;    // Specific enthalpy [J/kg]
    double s   = 0.0;    // Specific entropy [J/(kg·K)]
    double f   = 0.0;    // Local Darcy friction factor [-]
    double Re  = 0.0;    // Local Reynolds number [-]
};

// Result of Fanno flow pipe segment calculation
struct FannoSolution {
    State inlet;                        // Inlet thermodynamic state
    State outlet;                       // Outlet thermodynamic state
    double mdot    = 0.0;               // Mass flow rate [kg/s]
    double h0      = 0.0;               // Stagnation enthalpy [J/kg]
    double L       = 0.0;               // Pipe length [m]
    double D       = 0.0;               // Pipe diameter [m]
    double f       = 0.0;               // Darcy friction factor (constant-f overload) [-]
    double f_avg   = 0.0;               // Average Darcy friction factor over pipe [-]
    double Re_in   = 0.0;               // Inlet Reynolds number [-]
    bool choked    = false;             // True if flow reached M=1
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

// Solve Fanno flow with roughness-based variable friction factor.
//
// At each RK4 stage, the local friction factor is recomputed from the
// local Reynolds number and wall roughness using the specified correlation.
// This is more physical than a constant f when Re varies significantly
// along the pipe (e.g., large temperature or density changes).
//
// Inputs:
//   T_in        : Inlet static temperature [K]
//   P_in        : Inlet static pressure [Pa]
//   u_in        : Inlet velocity [m/s]
//   L           : Pipe length [m]
//   D           : Pipe diameter [m]
//   roughness   : Absolute wall roughness [m]
//   X           : Mole fractions [-]
//   correlation : Friction correlation (default: "haaland")
//                 Options: "haaland", "serghides", "colebrook"
//   n_steps     : Number of integration steps (default: 100)
//   store_profile : If true, store axial profile in solution
//
// Returns FannoSolution with f_avg and Re_in populated.
FannoSolution fanno_pipe_rough(
    double T_in, double P_in, double u_in,
    double L, double D, double roughness,
    const std::vector<double>& X,
    const std::string& correlation = "haaland",
    std::size_t n_steps = 100,
    bool store_profile = false);

// Convenience: solve given inlet State and velocity (roughness-based)
FannoSolution fanno_pipe_rough(
    const State& inlet, double u_in,
    double L, double D, double roughness,
    const std::string& correlation = "haaland",
    std::size_t n_steps = 100,
    bool store_profile = false);

// Compute maximum pipe length before choking (L*) for given inlet conditions.
// Returns the length at which M would reach 1.0.
double fanno_max_length(
    double T_in, double P_in, double u_in,
    double D, double f,
    const std::vector<double>& X,
    double tol = 1e-6, std::size_t max_iter = 100);

// -------------------------------------------------------------
// Rocket nozzle thrust
// -------------------------------------------------------------

// Thrust calculation result
struct ThrustResult {
    double thrust = 0.0;              // Thrust force [N]
    double specific_impulse = 0.0;    // Specific impulse Isp [s]
    double thrust_coefficient = 0.0;  // Thrust coefficient C_F [-]
    double mdot = 0.0;                // Mass flow rate [kg/s]
    double u_exit = 0.0;              // Exit velocity [m/s]
    double P_exit = 0.0;              // Exit pressure [Pa]
};

// Calculate rocket nozzle thrust from nozzle solution and ambient pressure.
//
// F = mdot * u_e + (P_e - P_amb) * A_e
// Isp = F / (mdot * g0)
// C_F = F / (P0 * A_throat)
//
// Inputs:
//   sol   : NozzleSolution from nozzle_cd() or similar
//   P_amb : Ambient pressure [Pa] (e.g., 101325 for sea level, 0 for vacuum)
//
// Returns ThrustResult with thrust, Isp, and thrust coefficient.
ThrustResult nozzle_thrust(const NozzleSolution& sol, double P_amb);

// Convenience: calculate thrust directly from nozzle parameters.
//
// P_design is the nozzle exit design pressure used to compute the Mach profile
// (must be > 0; use the nozzle exit pressure at the operating point).
// P_amb is the ambient pressure used for the thrust equation (may be 0 for vacuum).
// For a design-point calculation pass P_design == P_amb.
ThrustResult nozzle_thrust(
    double T0, double P0, double P_design, double P_amb,
    double A_inlet, double A_throat, double A_exit,
    double x_throat, double x_exit,
    const std::vector<double>& X,
    std::size_t n_stations = 100,
    double tol = 1e-8, std::size_t max_iter = 50);

#endif // COMPRESSIBLE_H
