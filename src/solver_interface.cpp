#include "solver_interface.h"

#include "friction.h"
#include "heat_transfer.h"
#include "math_constants.h"
#include "thermo.h"
#include "transport.h"
#include <cmath>
#include <stdexcept>

namespace combaero {
namespace solver {

// -----------------------------------------------------------------------------
// 1. Incompressible Flow Components
// -----------------------------------------------------------------------------

std::tuple<double, double> orifice_mdot_and_jacobian(double dP, double rho,
                                                     double Cd, double area) {
  // Basic orifice equation: mdot = Cd * A * sqrt(2 * rho * dP)
  // We strictly assume dP >= 0. For reverse flow, the solver wrapper handles
  // sign direction.
  if (dP < 0.0) {
    throw std::invalid_argument(
        "solver_interface::orifice_mdot requires dP >= 0");
  }

  // To cleanly avoid a divide-by-zero singularity in the gradient at dP = 0,
  // we return standard limits (mdot=0, d_mdot = infinity/large number) or bound
  // it. In practice, solvers rarely probe exactly at 0 if dP bounds are set.
  if (dP < 1e-12) {
    return {0.0, 1e12}; // Infinitely steep gradient at origin
  }

  double constant = Cd * area * std::sqrt(2.0 * rho);
  double mdot = constant * std::sqrt(dP);

  // Analytical derivative via chain rule:
  // d(mdot)/d(dP) = (Cd * A * sqrt(2*rho)) / (2 * sqrt(dP))
  //               = (Cd * A * sqrt(rho/2)) / sqrt(dP)
  // Alternatively: d(mdot)/d(dP) = mdot / (2 * dP)
  double jacobian = mdot / (2.0 * dP);

  return {mdot, jacobian};
}

std::tuple<double, double> pressure_loss_and_jacobian(double v, double rho,
                                                      double K) {
  // Standard K-factor loss: dP = K * 0.5 * rho * v^2
  double constant = K * 0.5 * rho;
  double dP = constant * v * v;

  // Analytical derivative via chain rule:
  // d(dP)/d(v) = K * rho * v
  double jacobian = K * rho * v;

  return {dP, jacobian};
}

// -----------------------------------------------------------------------------
// 2. Heat Transfer Components
// -----------------------------------------------------------------------------

std::tuple<double, double>
nusselt_and_jacobian_dittus_boelter(double Re, double Pr, bool heating) {
  // Dittus-Boelter fundamental relation: Nu = coeff_outer *
  // Re^coeff_reynolds_exp * Pr^n
  double n = heating ? dittus_boelter::coeff_prandtl_heating_exp
                     : dittus_boelter::coeff_prandtl_cooling_exp;

  if (Re <= 0.0) {
    return {0.0, 0.0};
  }

  double constant = dittus_boelter::coeff_outer * std::pow(Pr, n);
  double Nu = constant * std::pow(Re, dittus_boelter::coeff_reynolds_exp);

  // Analytical derivative via chain rule:
  // d(Nu)/d(Re) = coeff_outer * Pr^n * coeff_reynolds_exp *
  // Re^(coeff_reynolds_exp - 1)
  //             = (coeff_reynolds_exp * Nu) / Re
  double jacobian = (dittus_boelter::coeff_reynolds_exp * Nu) / Re;

  return {Nu, jacobian};
}

std::tuple<double, double> friction_and_jacobian_haaland(double Re,
                                                         double e_D) {
  // Explicit Haaland equation approximation of friction factor f (Darcy)
  // 1/sqrt(f) = coeff_outer * log10( (e_D/coeff_roughness)^coeff_exponent +
  // coeff_reynolds/Re )

  if (Re <= 0.0) {
    // Fallback or explicit failure;
    throw std::invalid_argument(
        "solver_interface::friction_haaland requires Re > 0");
  }

  // Isolate core logarithmic argument
  double a = std::pow(e_D / haaland::coeff_roughness, haaland::coeff_exponent);
  double b = haaland::coeff_reynolds / Re;
  double arg = a + b;

  // 1/sqrt(f)
  double inv_sqrt_f = haaland::coeff_outer * std::log10(arg);

  // f = [ coeff_outer * log10( (e_D/coeff_roughness)^coeff_exponent +
  // coeff_reynolds/Re ) ]^(-2)
  double f = 1.0 / (inv_sqrt_f * inv_sqrt_f);

  // Analytical derivative via chain rule:
  // d(f)/d(Re)
  // Let U = inv_sqrt_f = coeff_outer * log10(arg) = (coeff_outer / ln(10)) *
  // ln(arg) f = U^(-2) df/dU = -2 * U^(-3)
  //
  // dU/dRe = dU/d(arg) * d(arg)/dRe
  // dU/d(arg) = (coeff_outer / ln(10)) * (1 / arg)
  // d(arg)/dRe = -coeff_reynolds / Re^2
  //
  // Combining:
  // df/dRe = (-2 * U^(-3)) * (coeff_outer / (ln(10) * arg)) * (-coeff_reynolds
  // / Re^2)

  double dU_darg = haaland::coeff_outer / (M_LN10 * arg);
  double darg_dRe = -haaland::coeff_reynolds / (Re * Re);

  double dU_dRe = dU_darg * darg_dRe;
  double df_dU = -2.0 / (inv_sqrt_f * inv_sqrt_f * inv_sqrt_f);

  double jacobian = df_dU * dU_dRe;

  return {f, jacobian};
}

// -----------------------------------------------------------------------------
// 3. Thermodynamic & Transport Components
// -----------------------------------------------------------------------------

std::tuple<double, double, double>
density_and_jacobians(double T, double P, const std::vector<double> &X) {
  if (T <= 0.0 || P <= 0.0) {
    throw std::invalid_argument(
        "solver_interface::density requires T > 0 and P > 0");
  }

  // Calculate base density natively
  double rho = ::density(T, P, X);

  // Ideal gas law: rho = P * W / (R * T)
  // d(rho)/dT = - P * W / (R * T^2) = -rho / T
  double d_rho_d_T = -rho / T;

  // d(rho)/dP = W / (R * T) = rho / P
  double d_rho_d_P = rho / P;

  return {rho, d_rho_d_T, d_rho_d_P};
}

std::tuple<double, double> enthalpy_and_jacobian(double T,
                                                 const std::vector<double> &X) {
  if (T <= 0.0) {
    throw std::invalid_argument("solver_interface::enthalpy requires T > 0");
  }

  // h(T) = base mass enthalpy
  double h = ::h_mass(T, X);

  // Analytical derivative of enthalpy with respect to temperature is Cp
  double d_h_d_T = ::cp_mass(T, X);

  return {h, d_h_d_T};
}

std::tuple<double, double, double>
viscosity_and_jacobians(double T, double P, const std::vector<double> &X) {
  if (T <= 0.0 || P <= 0.0) {
    throw std::invalid_argument(
        "solver_interface::viscosity requires T > 0 and P > 0");
  }

  // Base viscosity
  double mu = ::viscosity(T, P, X);

  // We utilize a highly constrained central finite difference inside the C++
  // compiled kernel since the analytical collision integral derivations are
  // excessively complex. This masks the iteration noise from the global solver.

  // 1. Temperature Derivative
  double dT = 1e-3; // 1 mK perturbation
  double mu_T_plus = ::viscosity(T + dT, P, X);
  double mu_T_minus = ::viscosity(T - dT, P, X);
  double d_mu_d_T = (mu_T_plus - mu_T_minus) / (2.0 * dT);

  // 2. Pressure Derivative
  // Current models evaluate viscosity independently of pressure for ideal
  // gases, but we provide the generalized interface for completeness.
  double dP = 1.0; // 1 Pascal perturbation
  double mu_P_plus = ::viscosity(T, P + dP, X);
  double mu_P_minus = ::viscosity(T, P - dP, X);
  double d_mu_d_P = (mu_P_plus - mu_P_minus) / (2.0 * dP);

  return {mu, d_mu_d_T, d_mu_d_P};
}

} // namespace solver
} // namespace combaero
