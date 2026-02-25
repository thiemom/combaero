#include "solver_interface.h"

#include "math_constants.h"
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
  // Dittus-Boelter fundamental relation: Nu = 0.023 * Re^0.8 * Pr^n
  double n = heating ? 0.4 : 0.3;

  if (Re <= 0.0) {
    return {0.0, 0.0};
  }

  double constant = 0.023 * std::pow(Pr, n);
  double Nu = constant * std::pow(Re, 0.8);

  // Analytical derivative via chain rule:
  // d(Nu)/d(Re) = 0.023 * Pr^n * 0.8 * Re^(-0.2)
  //             = (0.8 * Nu) / Re
  double jacobian = (0.8 * Nu) / Re;

  return {Nu, jacobian};
}

std::tuple<double, double> friction_and_jacobian_haaland(double Re,
                                                         double e_D) {
  // Explicit Haaland equation approximation of friction factor f (Darcy)
  // 1/sqrt(f) = -1.8 * log10( (e_D/3.7)^1.11 + 6.9/Re )

  if (Re <= 0.0) {
    // Fallback or explicit failure;
    throw std::invalid_argument(
        "solver_interface::friction_haaland requires Re > 0");
  }

  // Isolate core logarithmic argument
  double a = std::pow(e_D / 3.7, 1.11);
  double b = 6.9 / Re;
  double arg = a + b;

  // 1/sqrt(f)
  double inv_sqrt_f = -1.8 * std::log10(arg);

  // f = [ -1.8 * log10( (e_D/3.7)^1.11 + 6.9/Re ) ]^(-2)
  double f = 1.0 / (inv_sqrt_f * inv_sqrt_f);

  // Analytical derivative via chain rule:
  // d(f)/d(Re)
  // Let U = inv_sqrt_f = -1.8 * log10(arg) = (-1.8 / ln(10)) * ln(arg)
  // f = U^(-2)
  // df/dU = -2 * U^(-3)
  //
  // dU/dRe = dU/d(arg) * d(arg)/dRe
  // dU/d(arg) = (-1.8 / ln(10)) * (1 / arg)
  // d(arg)/dRe = -6.9 / Re^2
  //
  // Combining:
  // df/dRe = (-2 * U^(-3)) * (-1.8 / (ln(10) * arg)) * (-6.9 / Re^2)
  //        = -3.14915... * (U^-3 / arg) * (-6.9 / Re^2) ...
  //        = -(12.42 * U^(-3)) / (ln(10) * arg * Re^2)

  double dU_darg = -1.8 / (M_LN10 * arg);
  double darg_dRe = -6.9 / (Re * Re);

  double dU_dRe = dU_darg * darg_dRe;
  double df_dU = -2.0 / (inv_sqrt_f * inv_sqrt_f * inv_sqrt_f);

  double jacobian = df_dU * dU_dRe;

  return {f, jacobian};
}

} // namespace solver
} // namespace combaero
