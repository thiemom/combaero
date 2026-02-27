#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#include <tuple>
#include <vector>
// -----------------------------------------------------------------------------
// Network Solver Fast-Path Native Interface
// -----------------------------------------------------------------------------
// Functions in this module execute highly optimized `(f, J)` combined
// evaluations required by the global Python `scipy.optimize.root` solver.
//
// Every function returns a `std::tuple<double, double>`, corresponding to:
//    (Value, Exact Analytical Jacobian / Derivative)
//
// This architecture natively evaluates chain-rule analytical gradients where
// possible, and tightly bounded finite-differences where strictly necessary,
// entirely shielding the Python runtime from grid discontinuity noise and
// locking the C++/Python return overhead to ~0.2 microseconds.

namespace combaero {
namespace solver {

// -----------------------------------------------------------------------------
// 1. Incompressible Flow Components
// -----------------------------------------------------------------------------

// Calculate Mass Flow and its derivative with respect to differential pressure.
//
// Method: Analytical Chain Rule
// Inputs:
//   dP   : Differential pressure drop across the orifice [Pa]
//   rho  : Fluid density [kg/m³]
//   Cd   : Discharge coefficient [-]
//   area : Orifice bore area [m²]
//
// Returns:
//   tuple(mdot, d_mdot_d_dP)
//   mdot        : Mass flow rate [kg/s]
//   d_mdot_d_dP : Jacobian gradient [kg/(s*Pa)]
std::tuple<double, double> orifice_mdot_and_jacobian(double dP, double rho,
                                                     double Cd, double area);

// Calculate Pressure Drop and its derivative with respect to velocity.
//
// Method: Analytical Chain Rule
// Inputs:
//   v   : Fluid velocity [m/s]
//   rho : Fluid density [kg/m³]
//   K   : Forward loss coefficient [-]
//
// Returns:
//   tuple(dP, d_dP_d_v)
//   d_dP_d_v : Jacobian gradient [Pa/(m/s)]
std::tuple<double, double> pressure_loss_and_jacobian(double v, double rho,
                                                      double K);

// Calculate Lossless Connection Ideal Pressure and Jacobian
//
// Method: Analytical Difference
// Inputs:
//   P_in  : Total pressure entering the element [Pa]
//   P_out : Total pressure exiting the element [Pa]
//
// Returns:
//   tuple(Residual, d_Residual_d_P_out)
//   Residual           : Difference between P_in and P_out [Pa]
//   d_Residual_d_P_out : Jacobian gradient [-]
std::tuple<double, double> lossless_pressure_and_jacobian(double P_in,
                                                          double P_out);

// -----------------------------------------------------------------------------
// 2. Heat Transfer Components
// -----------------------------------------------------------------------------

// Calculate Nusselt Number and its derivative with respect to Reynolds Number.
// Uses the Dittus-Boelter fundamental relation: Nu = 0.023 * Re^0.8 * Pr^n
//
// Method: Analytical Chain Rule
// Inputs:
//   Re      : Reynolds number [-]
//   Pr      : Prandtl number [-]
//   heating : True if fluid is being heated (n=0.4), False if cooled (n=0.3)
//
// Returns:
//   tuple(Nu, d_Nu_d_Re)
//   Nu        : Nusselt number [-]
//   d_Nu_d_Re : Jacobian gradient [-]
std::tuple<double, double>
nusselt_and_jacobian_dittus_boelter(double Re, double Pr, bool heating = true);

// Calculate Friction Factor and its derivative with respect to Reynolds Number.
// Uses the Haaland explicit pipe friction relation.
//
// Method: Analytical Chain Rule
// Inputs:
//   Re  : Reynolds number [-]
//   e_D : Relative roughness [-]
//
// Returns:
//   tuple(f, d_f_d_Re)
//   f        : Darcy friction factor [-]
//   d_f_d_Re : Jacobian gradient [-]
std::tuple<double, double> friction_and_jacobian_haaland(double Re, double e_D);

// -----------------------------------------------------------------------------
// 3. Thermodynamic & Transport Components
// -----------------------------------------------------------------------------

// Calculate Mass Density and its analytical derivatives w.r.t. T and P.
// Method: Ideal Gas Law Analytical Chain Rule
// Inputs:
//   T : Temperature [K]
//   P : Pressure [Pa]
//   X : Mole fractions [-]
//
// Returns:
//   tuple(rho, d_rho_d_T, d_rho_d_P)
//   rho       : Density [kg/m³]
//   d_rho_d_T : Jacobian gradient w.r.t Temperature [(kg/m³)/K]
//   d_rho_d_P : Jacobian gradient w.r.t Pressure [(kg/m³)/Pa]
std::tuple<double, double, double>
density_and_jacobians(double T, double P, const std::vector<double> &X);

// Calculate Mass Enthalpy and its analytical derivative w.r.t. T.
// Method: dh/dT is identically the mass-specific heat capacity (cp_mass).
// Inputs:
//   T : Temperature [K]
//   X : Mole fractions [-]
//
// Returns:
//   tuple(h, d_h_d_T)
//   h       : Mass Enthalpy [J/kg]
//   d_h_d_T : Jacobian gradient w.r.t Temperature [J/(kg*K)]
std::tuple<double, double> enthalpy_and_jacobian(double T,
                                                 const std::vector<double> &X);

// Calculate Dynamic Viscosity and its derivatives w.r.t. T and P.
// Method: Central Finite Difference (tightly bracketed)
// Inputs:
//   T : Temperature [K]
//   P : Pressure [Pa]
//   X : Mole fractions [-]
//
// Returns:
//   tuple(mu, d_mu_d_T, d_mu_d_P)
//   mu       : Dynamic viscosity [Pa*s]
//   d_mu_d_T : Jacobian gradient w.r.t Temperature [(Pa*s)/K]
//   d_mu_d_P : Jacobian gradient w.r.t Pressure [(Pa*s)/Pa]
std::tuple<double, double, double>
viscosity_and_jacobians(double T, double P, const std::vector<double> &X);

} // namespace solver
} // namespace combaero

#endif // SOLVER_INTERFACE_H
