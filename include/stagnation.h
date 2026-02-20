#ifndef STAGNATION_H
#define STAGNATION_H

#include <cstddef>
#include <vector>

// -------------------------------------------------------------
// Isentropic stagnation / static conversion utilities
// -------------------------------------------------------------
// These functions convert between static and stagnation (total) conditions
// for a thermally perfect gas using the variable-cp thermo backend (thermo.h).
//
// Conventions used throughout CombAero:
//   - thermo.h functions h(T,X), cp(T,X), s(T,P,X) operate on STATIC conditions.
//   - Stagnation enthalpy:  h0 = h(T_static) + v^2/2  [J/kg]
//   - Stagnation temperature T0: defined by h(T0) = h0  (iterative for variable cp)
//   - Stagnation pressure P0:    isentropic relation s(T0,P0) = s(T,P)
//
// Adiabatic wall temperature (for convective heat transfer at high Mach):
//   q = h_conv * (T_aw - T_wall)          <- correct driving temperature
//   T_aw = T_static + r * v^2 / (2*cp)
//   r = Pr^(1/3)  turbulent  (default)
//   r = Pr^(1/2)  laminar
//
//   T_aw is always between T_static and T_total:
//     T_static < T_aw < T_total  (since r < 1 for air)
//   Using T_static underestimates heat load at high M.
//   Using T_total overestimates (r < 1, so T_aw < T0 always).
//
// References:
//   Anderson, J.D. (2003) Modern Compressible Flow, 3rd ed., Ch. 3
//   Incropera et al. (2007) Fundamentals of Heat and Mass Transfer, 6th ed.

// -------------------------------------------------------------
// Kinetic energy
// -------------------------------------------------------------

// Specific kinetic energy: v^2/2  [J/kg]
inline double kinetic_energy(double v) { return 0.5 * v * v; }

// -------------------------------------------------------------
// Stagnation enthalpy (exact, no iteration)
// -------------------------------------------------------------

// Stagnation enthalpy from static enthalpy and velocity.
// h0 = h_static + v^2/2  [J/kg]
// h_static must be in [J/kg] (use h_mass() from thermo.h).
double h0_from_static(double h_static_J_per_kg, double v);

// Velocity from stagnation and static enthalpy.
// v = sqrt(2*(h0 - h_static))  [m/s]
// Returns 0 if h0 <= h_static.
double v_from_h0(double h0, double h_static_J_per_kg);

// -------------------------------------------------------------
// Mach number
// -------------------------------------------------------------

// Mach number from velocity and static state.
// M = v / a(T, X)
double mach_number(double v, double T, const std::vector<double>& X);

// -------------------------------------------------------------
// Stagnation temperature (iterative — variable cp)
// -------------------------------------------------------------

// Stagnation temperature T0 from static T and Mach number.
// Solves h(T0) = h(T) + v^2/2 iteratively using Newton on h(T).
// v = M * a(T, X)  is computed internally.
//
// Parameters:
//   T   : static temperature [K]
//   M   : Mach number [-]
//   X   : mole fractions [-]
// Returns: stagnation temperature T0 [K]
double T0_from_static(double T, double M, const std::vector<double>& X,
                      double tol = 1e-8, std::size_t max_iter = 50);

// Convenience: T0 from static T and velocity v [m/s]
double T0_from_static_v(double T, double v, const std::vector<double>& X,
                        double tol = 1e-8, std::size_t max_iter = 50);

// Static temperature from stagnation T0 and Mach number.
// Solves h(T0) = h(T) + M^2*a(T)^2/2 iteratively.
// Returns: static temperature T [K]
double T_from_stagnation(double T0, double M, const std::vector<double>& X,
                         double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Stagnation pressure (isentropic — uses entropy conservation)
// -------------------------------------------------------------

// Stagnation pressure P0 from static (P, T) and Mach number.
// Uses isentropic relation: s(T0, P0) = s(T, P)
// Solves for P0 given T0 = T0_from_static(T, M, X).
//
// Parameters:
//   P   : static pressure [Pa]
//   T   : static temperature [K]
//   M   : Mach number [-]
//   X   : mole fractions [-]
// Returns: stagnation pressure P0 [Pa]
double P0_from_static(double P, double T, double M,
                      const std::vector<double>& X,
                      double tol = 1e-8, std::size_t max_iter = 50);

// Static pressure from stagnation (P0, T0) and Mach number.
// Uses isentropic relation: s(T, P) = s(T0, P0)
// Returns: static pressure P [Pa]
double P_from_stagnation(double P0, double T0, double M,
                         const std::vector<double>& X,
                         double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Adiabatic wall temperature
// -------------------------------------------------------------
// For convective heat transfer, the correct hot-side driving temperature is
// T_aw, not T_static or T_total.  Pass T_aw as T_hot to:
//   heat_transfer.h : htc_pipe(), nusselt_pipe()
//   cooling_correlations.h : cooled_wall_heat_flux()
//
// For M < 0.3: T_aw ~ T_static (< 1% error, correction negligible).
// For M = 0.3-0.8 (combustor liner, turbine cooling): correction is 2-15%.

// Adiabatic wall temperature from static temperature and velocity.
// T_aw = T_static + r * v^2 / (2 * cp(T, X))
// r = Pr^(1/3) for turbulent flow (default), Pr^(1/2) for laminar.
//
// Parameters:
//   T_static  : bulk static temperature [K]
//   v         : bulk flow velocity [m/s]
//   T         : temperature for property evaluation [K] (use T_static)
//   P         : pressure [Pa]
//   X         : mole fractions [-]
//   turbulent : true = turbulent (r = Pr^1/3), false = laminar (r = Pr^1/2)
// Returns: adiabatic wall temperature T_aw [K]
double T_adiabatic_wall(double T_static, double v,
                        double T, double P, const std::vector<double>& X,
                        bool turbulent = true);

// Convenience: T_aw from static temperature and Mach number.
// Computes v = M * a(T, X) internally.
double T_adiabatic_wall_mach(double T_static, double M,
                              double T, double P, const std::vector<double>& X,
                              bool turbulent = true);

// Recovery factor r [-]
// r = Pr^(1/3) for turbulent, Pr^(1/2) for laminar
double recovery_factor(double Pr, bool turbulent = true);

#endif // STAGNATION_H
