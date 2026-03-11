#include "solver_interface.h"

#include <algorithm>
#include <cmath>
#include "../include/composition.h"
using combaero::mwmix;
using combaero::mole_to_mass;
using combaero::mass_to_mole;
using combaero::normalize_fractions;
#include <numeric>
#include <stdexcept>

#include "combustion.h"
#include "cooling_correlations.h"
#include "equilibrium.h"
#include "friction.h"
#include "heat_transfer.h"
#include "math_constants.h" // MSVC compatibility for M_PI, M_LN10, etc.
#include "registry.h"
#include "stagnation.h"
#include "thermo.h"
#include "transport.h"

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

std::tuple<double, double> lossless_pressure_and_jacobian(double P_in,
                                                          double P_out) {
  // A lossless connection preserves total pressure natively.
  // Residual: P_in - P_out = 0.
  // Analytical derivative w.r.t P_out: -1.0
  return {P_in - P_out, -1.0};
}

// -----------------------------------------------------------------------------
// 2. Heat Transfer Components
// -----------------------------------------------------------------------------

CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_dittus_boelter(double Re, double Pr, bool heating) {
  if (Re <= 0.0) {
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "nusselt_dittus_boelter requires Re > 0"};
  }

  double n = heating ? dittus_boelter::coeff_prandtl_heating_exp
                     : dittus_boelter::coeff_prandtl_cooling_exp;

  double constant = dittus_boelter::coeff_outer * std::pow(Pr, n);
  double Nu = constant * std::pow(Re, dittus_boelter::coeff_reynolds_exp);

  // dNu_dRe = coeff_outer * coeff_reynolds_exp * Re^(coeff_reynolds_exp - 1) *
  // Pr^n
  //         = coeff_reynolds_exp * Nu / Re
  double dNu_dRe = dittus_boelter::coeff_reynolds_exp * Nu / Re;

  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_gnielinski(double Re, double Pr, double f) {
  if (Re <= 0.0) {
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "nusselt_gnielinski requires Re > 0"};
  }
  // Gnielinski uses Hermite blending and is complex; use tight analytical-like
  // FD
  const double eps =
      std::max(1e-6, Re * 1e-6); // tightly bracketed minimum perturbation
  double Nu_plus = nusselt_gnielinski(Re + eps, Pr, f);
  double Nu_minus = nusselt_gnielinski(Re - eps, Pr, f);
  double dNu_dRe = (Nu_plus - Nu_minus) / (2.0 * eps);
  double Nu = nusselt_gnielinski(Re, Pr, f); // center evaluation
  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_sieder_tate(double Re, double Pr, double mu_ratio) {
  if (Re <= 0.0) {
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "nusselt_sieder_tate requires Re > 0"};
  }
  double Nu = nusselt_sieder_tate(Re, Pr, mu_ratio);

  // Nu = 0.027 * Re^0.8 * Pr^(1/3) * mu_ratio^0.14
  // dNu/dRe = 0.8 * Nu / Re
  double dNu_dRe = sieder_tate::coeff_reynolds_exp * Nu / Re;
  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_petukhov(double Re, double Pr, double f) {
  if (Re <= 0.0) {
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "nusselt_petukhov requires Re > 0"};
  }
  double Nu = nusselt_petukhov(Re, Pr, f);

  // Nu = (f/8) * Re * Pr / D
  // Approximating dNu/dRe by assuming f is constant, standard in sequenced
  // networks where friction and heat transfer are computed sequentially but
  // decoupled in the local element Jacobian block
  double dNu_dRe = Nu / Re;
  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_haaland(double Re, double e_D) {
  // Explicit Haaland equation approximation of friction factor f (Darcy)
  // 1/sqrt(f) = coeff_outer * log10( (e_D/coeff_roughness)^coeff_exponent +
  // coeff_reynolds/Re )

  if (Re <= 0.0) {
    // Graceful invalidation for unphysical states
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "solver_interface::friction_haaland requires Re > 0"};
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

  return {{f, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_serghides(double Re, double e_D) {
  if (Re <= 0.0) {
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "solver_interface::friction_serghides requires Re > 0"};
  }

  // Serghides Steffensen form: 1/sqrt(f) = A - (B-A)^2 / (C - 2B + A)
  // A = -2 log10(e_D/3.7 + 12/Re)
  // B = -2 log10(e_D/3.7 + 2.51*A/Re)
  // C = -2 log10(e_D/3.7 + 2.51*B/Re)
  const double r = serghides::coeff_roughness;
  const double cA = serghides::coeff_reynolds_A;
  const double cBC = serghides::coeff_reynolds_BC;
  const double ln10 = M_LN10;

  double argA = e_D / r + cA / Re;
  double A = -2.0 * std::log10(argA);

  double argB = e_D / r + cBC * A / Re;
  double B = -2.0 * std::log10(argB);

  double argC = e_D / r + cBC * B / Re;
  double C = -2.0 * std::log10(argC);

  double denom = C - 2.0 * B + A;
  double inv_sqrt_f =
      (std::abs(denom) < 1e-30) ? C : A - (B - A) * (B - A) / denom;
  double f = 1.0 / (inv_sqrt_f * inv_sqrt_f);

  // Analytical derivative via chain rule through A, B, C, inv_sqrt_f, f.
  // dA/dRe = (2 / (ln10 * argA)) * (cA / Re^2)
  double dA_dRe = (2.0 / (ln10 * argA)) * (cA / (Re * Re));

  // dB/dRe = (2 / (ln10 * argB)) * (cBC/Re) * (A/Re - dA_dRe)
  //        = (2 / (ln10 * argB)) * cBC * (A/Re^2 - dA_dRe/Re) ... expand:
  // d(argB)/dRe = cBC * (dA_dRe * Re - A) / Re^2 = cBC * (dA_dRe/Re - A/Re^2)
  double dargB_dRe = cBC * (dA_dRe / Re - A / (Re * Re));
  double dB_dRe = (-2.0 / (ln10 * argB)) * dargB_dRe;

  double dargC_dRe = cBC * (dB_dRe / Re - B / (Re * Re));
  double dC_dRe = (-2.0 / (ln10 * argC)) * dargC_dRe;

  // d(inv_sqrt_f)/dRe from Steffensen formula
  double d_inv_sqrt_f_dRe;
  if (std::abs(denom) < 1e-30) {
    d_inv_sqrt_f_dRe = dC_dRe;
  } else {
    // inv_sqrt_f = A - (B-A)^2 / denom,  denom = C - 2B + A
    double BmA = B - A;
    double d_BmA_dRe = dB_dRe - dA_dRe;
    double d_denom_dRe = dC_dRe - 2.0 * dB_dRe + dA_dRe;
    // d((B-A)^2/denom)/dRe = (2*(B-A)*d(B-A) * denom - (B-A)^2 * d_denom) /
    // denom^2
    d_inv_sqrt_f_dRe =
        dA_dRe - (2.0 * BmA * d_BmA_dRe * denom - BmA * BmA * d_denom_dRe) /
                     (denom * denom);
  }

  // f = inv_sqrt_f^(-2)  => df/dRe = -2 * inv_sqrt_f^(-3) * d_inv_sqrt_f_dRe
  double jacobian =
      (-2.0 / (inv_sqrt_f * inv_sqrt_f * inv_sqrt_f)) * d_inv_sqrt_f_dRe;

  return {{f, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_colebrook(double Re, double e_D) {
  if (Re <= 0.0) {
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "solver_interface::friction_colebrook requires Re > 0"};
  }

  // Colebrook converges internally; use Serghides as starting guess.
  // After convergence, the implicit derivative is:
  //   F(x, Re) = x + 2*log10(e_D/3.7 + 2.51*x/Re) = 0,  x = 1/sqrt(f)
  //   dF/dx  = 1 + 2/(ln10 * arg) * 2.51/Re  = 1 + 2*2.51/(ln10*arg*Re)
  //   dF/dRe = 2/(ln10 * arg) * (-2.51*x/Re^2)
  //   dx/dRe = -(dF/dRe) / (dF/dx)
  //   df/dRe = -2 * x^(-3) * dx/dRe

  // Get converged x = 1/sqrt(f) from the colebrook scalar function
  double f_val = ::friction_colebrook(Re, e_D);
  double x = 1.0 / std::sqrt(f_val);

  const double ln10 = M_LN10;
  double arg =
      e_D / serghides::coeff_roughness + serghides::coeff_reynolds_BC * x / Re;

  double dF_dx = 1.0 + 2.0 * serghides::coeff_reynolds_BC / (ln10 * arg * Re);
  double dF_dRe =
      (-2.0 / (ln10 * arg)) * (serghides::coeff_reynolds_BC * x / (Re * Re));

  double dx_dRe = -dF_dRe / dF_dx;
  double jacobian = (-2.0 / (x * x * x)) * dx_dRe;

  return {{f_val, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_petukhov(double Re) {
  if (Re < 3000.0) {
    return {{0.0, 0.0},
            CorrelationValidity::EXTRAPOLATED,
            "solver_interface::friction_petukhov requires Re >= 3000"};
  }

  // f = (coeff_a * ln(Re) - coeff_b)^(-2)
  double x = petukhov::coeff_a * std::log(Re) - petukhov::coeff_b;
  double f = 1.0 / (x * x);

  // df/dRe = -2 * x^(-3) * dx/dRe
  // dx/dRe = coeff_a / Re
  double dx_dRe = petukhov::coeff_a / Re;
  double jacobian = (-2.0 / (x * x * x)) * dx_dRe;

  return {{f, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian(const std::string &tag, double Re, double e_D) {
  if (Registry::get().has_friction_model(tag)) {
    return Registry::get().get_friction_model(tag)(Re, e_D);
  }
  return {{0.0, 0.0},
          CorrelationValidity::INVALID,
          "friction_and_jacobian: unknown tag '" + tag +
              "'. Must be a registered friction model."};
}

// -----------------------------------------------------------------------------
// 3. Cooling Correlations
// -----------------------------------------------------------------------------

CorrelationResult<std::tuple<double, double>>
pin_fin_nusselt_and_jacobian(double Re_d, double Pr, double L_D, double S_D,
                             double X_D, bool is_staggered) {
  const double eps = std::max(1e-6, Re_d * 1e-6);
  double Nu_plus =
      cooling::pin_fin_nusselt(Re_d + eps, Pr, L_D, S_D, X_D, is_staggered);
  double Nu_minus =
      cooling::pin_fin_nusselt(Re_d - eps, Pr, L_D, S_D, X_D, is_staggered);
  double dNu_dRe = (Nu_plus - Nu_minus) / (2.0 * eps);
  double Nu = cooling::pin_fin_nusselt(Re_d, Pr, L_D, S_D, X_D, is_staggered);
  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
pin_fin_friction_and_jacobian(double Re_d, bool is_staggered) {
  const double eps = std::max(1e-6, Re_d * 1e-6);
  double f_plus = cooling::pin_fin_friction(Re_d + eps, is_staggered);
  double f_minus = cooling::pin_fin_friction(Re_d - eps, is_staggered);
  double df_dRe = (f_plus - f_minus) / (2.0 * eps);
  double f = cooling::pin_fin_friction(Re_d, is_staggered);
  return {{f, df_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
dimple_nusselt_enhancement_and_jacobian(double Re_Dh, double d_Dh, double h_d,
                                        double S_d) {
  const double eps = std::max(1e-6, Re_Dh * 1e-6);
  double Nu_plus =
      cooling::dimple_nusselt_enhancement(Re_Dh + eps, d_Dh, h_d, S_d);
  double Nu_minus =
      cooling::dimple_nusselt_enhancement(Re_Dh - eps, d_Dh, h_d, S_d);
  double dNu_dRe = (Nu_plus - Nu_minus) / (2.0 * eps);
  double Nu = cooling::dimple_nusselt_enhancement(Re_Dh, d_Dh, h_d, S_d);
  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
dimple_friction_multiplier_and_jacobian(double Re_Dh, double d_Dh, double h_d) {
  const double eps = std::max(1e-6, Re_Dh * 1e-6);
  double f_plus = cooling::dimple_friction_multiplier(Re_Dh + eps, d_Dh, h_d);
  double f_minus = cooling::dimple_friction_multiplier(Re_Dh - eps, d_Dh, h_d);
  double df_dRe = (f_plus - f_minus) / (2.0 * eps);
  double f = cooling::dimple_friction_multiplier(Re_Dh, d_Dh, h_d);
  return {{f, df_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
rib_enhancement_factor_high_re_and_jacobian(double e_D, double pitch_to_height,
                                            double alpha_deg, double Re) {
  const double eps = std::max(1e-6, Re * 1e-6);
  double E_plus = cooling::rib_enhancement_factor_high_re(e_D, pitch_to_height,
                                                          alpha_deg, Re + eps);
  double E_minus = cooling::rib_enhancement_factor_high_re(e_D, pitch_to_height,
                                                           alpha_deg, Re - eps);
  double dE_dRe = (E_plus - E_minus) / (2.0 * eps);
  double E = cooling::rib_enhancement_factor_high_re(e_D, pitch_to_height,
                                                     alpha_deg, Re);
  return {{E, dE_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
impingement_nusselt_and_jacobian(double Re_jet, double Pr, double z_D,
                                 double x_D, double y_D) {
  const double eps = std::max(1e-6, Re_jet * 1e-6);
  double Nu_plus =
      cooling::impingement_nusselt(Re_jet + eps, Pr, z_D, x_D, y_D);
  double Nu_minus =
      cooling::impingement_nusselt(Re_jet - eps, Pr, z_D, x_D, y_D);
  double dNu_dRe = (Nu_plus - Nu_minus) / (2.0 * eps);
  double Nu = cooling::impingement_nusselt(Re_jet, Pr, z_D, x_D, y_D);
  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
film_cooling_effectiveness_and_jacobian(double x_D, double M, double DR,
                                        double alpha_deg) {
  const double eps = std::max(1e-6, M * 1e-6); // wrt Blowing Ratio M
  double eta_plus =
      cooling::film_cooling_effectiveness(x_D, M + eps, DR, alpha_deg);
  double eta_minus =
      cooling::film_cooling_effectiveness(x_D, M - eps, DR, alpha_deg);
  double deta_dM = (eta_plus - eta_minus) / (2.0 * eps);
  double eta = cooling::film_cooling_effectiveness(x_D, M, DR, alpha_deg);
  return {{eta, deta_dM}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
effusion_effectiveness_and_jacobian(double x_D, double M, double DR,
                                    double porosity, double s_D,
                                    double alpha_deg) {
  const double eps = std::max(1e-6, M * 1e-6); // wrt Blowing Ratio M
  double eta_plus = cooling::effusion_effectiveness(x_D, M + eps, DR, porosity,
                                                    s_D, alpha_deg);
  double eta_minus = cooling::effusion_effectiveness(x_D, M - eps, DR, porosity,
                                                     s_D, alpha_deg);
  double deta_dM = (eta_plus - eta_minus) / (2.0 * eps);
  double eta =
      cooling::effusion_effectiveness(x_D, M, DR, porosity, s_D, alpha_deg);
  return {{eta, deta_dM}, CorrelationValidity::VALID, ""};
}

// -----------------------------------------------------------------------------
// 4. Thermodynamic & Transport Components
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

std::tuple<double, double>
mach_number_and_jacobian_v(double v, double T, const std::vector<double> &X) {
  double a = speed_of_sound(T, X);
  if (a <= 0.0) {
    throw std::invalid_argument(
        "mach_number_and_jacobian_v: speed of sound is non-positive");
  }
  double M = v / a;
  double dM_dv = 1.0 / a;
  return {M, dM_dv};
}

std::tuple<double, double>
T_adiabatic_wall_and_jacobian_v(double T_static, double v, double T, double P,
                                const std::vector<double> &X, bool turbulent) {
  if (T_static <= 0.0) {
    throw std::invalid_argument(
        "T_adiabatic_wall_and_jacobian_v: T_static must be positive");
  }
  if (v < 0.0) {
    throw std::invalid_argument(
        "T_adiabatic_wall_and_jacobian_v: v must be non-negative");
  }

  double Pr = prandtl(T, P, X);
  double r = turbulent ? std::cbrt(Pr) : std::sqrt(Pr);

  // mwmix(X) is g/mol. cp(T,X) is J/(mol*K).
  // cp_mass_at = cp(T, X) / mwmix(X) * 1000.0
  double mw_g = mwmix(X);
  double cp_val = cp(T, X) / mw_g * 1000.0;

  // T_aw = T_static + r * 0.5 * v^2 / cp_val
  // dT_aw_dv = r * v / cp_val
  double T_aw = T_static + r * 0.5 * v * v / cp_val;
  double dT_aw_dv = r * v / cp_val;

  return {T_aw, dT_aw_dv};
}

std::tuple<double, double>
T0_from_static_and_jacobian_M(double T, double M,
                              const std::vector<double> &X) {
  if (M < 0.0) {
    throw std::invalid_argument(
        "T0_from_static_and_jacobian_M: M must be non-negative");
  }
  const double eps = std::max(1e-6, M * 1e-6);
  double T0_plus = T0_from_static(T, M + eps, X);
  double T0_minus = T0_from_static(T, std::max(0.0, M - eps), X);

  double dT0_dM;
  if (M - eps < 0.0) {
    // Forward difference near zero
    dT0_dM = (T0_plus - T0_from_static(T, 0.0, X)) / eps;
  } else {
    dT0_dM = (T0_plus - T0_minus) / (2.0 * eps);
  }

  double T0 = T0_from_static(T, M, X);
  return {T0, dT0_dM};
}

std::tuple<double, double>
P0_from_static_and_jacobian_M(double P, double T, double M,
                              const std::vector<double> &X) {
  if (M < 0.0) {
    throw std::invalid_argument(
        "P0_from_static_and_jacobian_M: M must be non-negative");
  }
  const double eps = std::max(1e-6, M * 1e-6);
  double P0_plus = P0_from_static(P, T, M + eps, X);
  double P0_minus = P0_from_static(P, T, std::max(0.0, M - eps), X);

  double dP0_dM;
  if (M - eps < 0.0) {
    // Forward difference near zero
    dP0_dM = (P0_plus - P0_from_static(P, T, 0.0, X)) / eps;
  } else {
    dP0_dM = (P0_plus - P0_minus) / (2.0 * eps);
  }

  double P0 = P0_from_static(P, T, M, X);
  return {P0, dP0_dM};
}

// -----------------------------------------------------------------------------
// 6. Combustion Interfaces
// -----------------------------------------------------------------------------

std::tuple<double, double, std::vector<double>>
adiabatic_T_complete_and_jacobian_T(double T_in, double P,
                                    const std::vector<double> &X_in) {
  if (T_in <= 0.0 || P <= 0.0) {
    throw std::invalid_argument(
        "solver_interface::adiabatic_T_complete requires T_in > 0 and P > 0");
  }

  State in_state;
  in_state.set_TPX(T_in, P, X_in);
  State out_state = complete_combustion(in_state, true);
  double T_ad = out_state.T;
  std::vector<double> X_products = out_state.X;

  // Jacobian wrt T_in via heavily bracketed Central Finite Difference
  double dT_eps = std::max(1e-3, T_in * 1e-4);

  State in_plus, in_minus;
  in_plus.set_TPX(T_in + dT_eps, P, X_in);
  in_minus.set_TPX(std::max(0.1, T_in - dT_eps), P, X_in);

  double T_ad_plus = complete_combustion(in_plus, true).T;
  double T_ad_minus = complete_combustion(in_minus, true).T;

  double dT_ad_dT_in = (T_ad_plus - T_ad_minus) / (2.0 * dT_eps);

  return {T_ad, dT_ad_dT_in, X_products};
}

std::tuple<double, double, double, std::vector<double>>
adiabatic_T_equilibrium_and_jacobians(double T_in, double P,
                                      const std::vector<double> &X_in) {
  if (T_in <= 0.0 || P <= 0.0) {
    throw std::invalid_argument("solver_interface::adiabatic_T_equilibrium "
                                "requires T_in > 0 and P > 0");
  }

  State in_state;
  in_state.set_TPX(T_in, P, X_in);

  // Combustion Equilibrium
  EquilibriumResult eq_res = combustion_equilibrium(in_state, true);
  double T_ad = eq_res.state.T;
  std::vector<double> X_products = eq_res.state.X;

  // Jacobians via heavily bracketed Central Finite Difference
  double dT_eps = std::max(1e-3, T_in * 1e-4);
  double dP_eps = std::max(1.0, P * 1e-4);

  // wrt T_in
  State in_T_plus, in_T_minus;
  in_T_plus.set_TPX(T_in + dT_eps, P, X_in);
  in_T_minus.set_TPX(std::max(0.1, T_in - dT_eps), P, X_in);

  double T_ad_T_plus = combustion_equilibrium(in_T_plus, true).state.T;
  double T_ad_T_minus = combustion_equilibrium(in_T_minus, true).state.T;
  double dT_ad_dT_in = (T_ad_T_plus - T_ad_T_minus) / (2.0 * dT_eps);

  // wrt P
  State in_P_plus, in_P_minus;
  in_P_plus.set_TPX(T_in, P + dP_eps, X_in);
  in_P_minus.set_TPX(T_in, std::max(0.1, P - dP_eps), X_in);

  double T_ad_P_plus = combustion_equilibrium(in_P_plus, true).state.T;
  double T_ad_P_minus = combustion_equilibrium(in_P_minus, true).state.T;
  double dT_ad_dP = (T_ad_P_plus - T_ad_P_minus) / (2.0 * dP_eps);

  return {T_ad, dT_ad_dT_in, dT_ad_dP, X_products};
}

// -----------------------------------------------------------------------------
// 7. Stream-Based Network Solvers
// -----------------------------------------------------------------------------

MixerResult
mixer_from_streams_and_jacobians(const std::vector<Stream> &streams) {
  std::size_t n_streams = streams.size();
  std::size_t n_species = num_species();

  double mdot_tot = 0.0;
  std::vector<double> Y_mix(n_species, 0.0);
  double H_tot = 0.0;

  std::vector<double> h_stream(n_streams, 0.0);
  std::vector<double> cp_stream(n_streams, 0.0);
  std::vector<std::vector<double>> hk_mass(n_streams,
                                           std::vector<double>(n_species, 0.0));

  for (std::size_t i = 0; i < n_streams; ++i) {
    mdot_tot += streams[i].m_dot;

    for (std::size_t k = 0; k < n_species; ++k) {
      Y_mix[k] += streams[i].m_dot * streams[i].Y[k];
      hk_mass[i][k] = J_per_mol_to_J_per_kg(h_species(k, streams[i].T),
                                            species_molar_mass(k));
      h_stream[i] += streams[i].Y[k] * hk_mass[i][k];

      double cpk_mass = J_per_mol_to_J_per_kg(cp_species(k, streams[i].T),
                                              species_molar_mass(k));
      cp_stream[i] += streams[i].Y[k] * cpk_mass;
    }
    H_tot += streams[i].m_dot * h_stream[i];
  }

  if (mdot_tot > 0.0) {
    for (std::size_t k = 0; k < n_species; ++k) {
      Y_mix[k] /= mdot_tot;
    }
  } else {
    if (n_streams > 0) {
      Y_mix = streams[0].Y;
    }
  }

  double h_mix = (mdot_tot > 0.0) ? (H_tot / mdot_tot)
                                  : (n_streams > 0 ? h_stream[0] : 0.0);

  std::vector<double> X_mix = mass_to_mole(normalize_fractions(Y_mix));
  double T_guess = 300.0;
  if (n_streams > 0)
    T_guess = streams[0].T;
  double T_mix = calc_T_from_h_mass(h_mix, X_mix, T_guess);

  std::vector<std::vector<StreamJacobian>> dY_mix_d_stream(
      n_species, std::vector<StreamJacobian>(n_streams));
  for (std::size_t k = 0; k < n_species; ++k) {
    for (std::size_t i = 0; i < n_streams; ++i) {
      StreamJacobian &jac = dY_mix_d_stream[k][i];
      jac.d_Y.assign(n_species, 0.0);
      if (mdot_tot > 0.0) {
        jac.d_mdot = (streams[i].Y[k] - Y_mix[k]) / mdot_tot;
        jac.d_T = 0.0;
        jac.d_Y[k] = streams[i].m_dot / mdot_tot;
      } else {
        jac.d_mdot = 0.0;
        jac.d_T = 0.0;
      }
    }
  }

  double cp_mix = 0.0;
  std::vector<double> hk_mix(n_species, 0.0);
  for (std::size_t k = 0; k < n_species; ++k) {
    hk_mix[k] =
        J_per_mol_to_J_per_kg(h_species(k, T_mix), species_molar_mass(k));
    double cpk_mix =
        J_per_mol_to_J_per_kg(cp_species(k, T_mix), species_molar_mass(k));
    cp_mix += Y_mix[k] * cpk_mix;
  }

  std::vector<StreamJacobian> dT_mix_d_stream(n_streams);
  for (std::size_t i = 0; i < n_streams; ++i) {
    StreamJacobian &jac = dT_mix_d_stream[i];
    jac.d_Y.assign(n_species, 0.0);

    if (mdot_tot > 0.0) {
      double dh_mix_d_mdot = (h_stream[i] - h_mix) / mdot_tot;
      double dh_mix_d_Ti = (streams[i].m_dot * cp_stream[i]) / mdot_tot;

      jac.d_T = dh_mix_d_Ti / cp_mix;

      double sum_hk_dY = 0.0;
      for (std::size_t k = 0; k < n_species; ++k) {
        sum_hk_dY += hk_mix[k] * dY_mix_d_stream[k][i].d_mdot;
      }
      jac.d_mdot = (dh_mix_d_mdot - sum_hk_dY) / cp_mix;

      for (std::size_t k = 0; k < n_species; ++k) {
        double dh_mix_d_Yki = (streams[i].m_dot * hk_mass[i][k]) / mdot_tot;
        double dY_mix_k_d_Yki = dY_mix_d_stream[k][i].d_Y[k];
        jac.d_Y[k] = (dh_mix_d_Yki - hk_mix[k] * dY_mix_k_d_Yki) / cp_mix;
      }
    } else {
      jac.d_mdot = 0.0;
      jac.d_T = 0.0;
    }
  }

  return {T_mix, Y_mix, dT_mix_d_stream, dY_mix_d_stream};
}

MixerResult adiabatic_T_complete_and_jacobian_T_from_streams(
    const std::vector<Stream> &streams, double P) {
  MixerResult mix = mixer_from_streams_and_jacobians(streams);

  std::vector<double> X_mix = mass_to_mole(normalize_fractions(mix.Y_mix));
  auto [T_ad, dT_ad_dT_mix, X_b] =
      adiabatic_T_complete_and_jacobian_T(mix.T_mix, P, X_mix);
  std::vector<double> Y_b = mole_to_mass(normalize_fractions(X_b));

  std::size_t n_streams = streams.size();
  std::size_t n_species = num_species();

  std::vector<double> dT_ad_dY_mix(n_species, 0.0);
  std::vector<std::vector<double>> dYb_dY_mix(
      n_species, std::vector<double>(n_species, 0.0));
  const double eps_Y = 1e-6;

  for (std::size_t k = 0; k < n_species; ++k) {
    std::vector<double> Y_p = mix.Y_mix;
    Y_p[k] += eps_Y;
    std::vector<double> X_p = mass_to_mole(normalize_fractions(Y_p));
    auto [T_ad_p, dummy, X_b_p] =
        adiabatic_T_complete_and_jacobian_T(mix.T_mix, P, X_p);
    std::vector<double> Y_b_p = mole_to_mass(normalize_fractions(X_b_p));

    dT_ad_dY_mix[k] = (T_ad_p - T_ad) / eps_Y;
    for (std::size_t j = 0; j < n_species; ++j) {
      dYb_dY_mix[j][k] = (Y_b_p[j] - Y_b[j]) / eps_Y;
    }
  }

  std::vector<double> dYb_dT_mix(n_species, 0.0);

  MixerResult res;
  res.T_mix = T_ad;
  res.Y_mix = Y_b;
  res.dT_mix_d_stream.resize(n_streams);
  res.dY_mix_d_stream.assign(n_species, std::vector<StreamJacobian>(n_streams));

  for (std::size_t i = 0; i < n_streams; ++i) {
    StreamJacobian &dT_stream = res.dT_mix_d_stream[i];
    dT_stream.d_Y.assign(n_species, 0.0);

    auto compute_chain = [&](double d_prop_T,
                             const std::vector<double> &d_prop_Y) {
      double sum = dT_ad_dT_mix * d_prop_T;
      for (std::size_t k = 0; k < n_species; ++k) {
        sum += dT_ad_dY_mix[k] * d_prop_Y[k];
      }
      return sum;
    };

    std::vector<double> dY_mix_dmdot(n_species), dY_mix_dT(n_species);
    for (std::size_t k = 0; k < n_species; ++k) {
      dY_mix_dmdot[k] = mix.dY_mix_d_stream[k][i].d_mdot;
      dY_mix_dT[k] = mix.dY_mix_d_stream[k][i].d_T;
    }

    dT_stream.d_mdot =
        compute_chain(mix.dT_mix_d_stream[i].d_mdot, dY_mix_dmdot);
    dT_stream.d_T = compute_chain(mix.dT_mix_d_stream[i].d_T, dY_mix_dT);
    for (std::size_t j = 0; j < n_species; ++j) {
      std::vector<double> dY_j(n_species);
      for (std::size_t k = 0; k < n_species; ++k)
        dY_j[k] = mix.dY_mix_d_stream[k][i].d_Y[j];
      dT_stream.d_Y[j] = compute_chain(mix.dT_mix_d_stream[i].d_Y[j], dY_j);
    }

    for (std::size_t tgt = 0; tgt < n_species; ++tgt) {
      StreamJacobian &dY_stream = res.dY_mix_d_stream[tgt][i];
      dY_stream.d_Y.assign(n_species, 0.0);

      auto compute_Y_chain = [&](double d_prop_T,
                                 const std::vector<double> &d_prop_Y) {
        double sum = dYb_dT_mix[tgt] * d_prop_T;
        for (std::size_t k = 0; k < n_species; ++k) {
          sum += dYb_dY_mix[tgt][k] * d_prop_Y[k];
        }
        return sum;
      };

      dY_stream.d_mdot =
          compute_Y_chain(mix.dT_mix_d_stream[i].d_mdot, dY_mix_dmdot);
      dY_stream.d_T = compute_Y_chain(mix.dT_mix_d_stream[i].d_T, dY_mix_dT);
      for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> dY_j(n_species);
        for (std::size_t k = 0; k < n_species; ++k)
          dY_j[k] = mix.dY_mix_d_stream[k][i].d_Y[j];
        dY_stream.d_Y[j] = compute_Y_chain(mix.dT_mix_d_stream[i].d_Y[j], dY_j);
      }
    }
  }
  return res;
}

MixerResult adiabatic_T_equilibrium_and_jacobians_from_streams(
    const std::vector<Stream> &streams, double P) {
  MixerResult mix = mixer_from_streams_and_jacobians(streams);

  std::vector<double> X_mix = mass_to_mole(normalize_fractions(mix.Y_mix));
  auto [T_ad, dT_ad_dT_mix, dT_ad_dP_mix, X_b] =
      adiabatic_T_equilibrium_and_jacobians(mix.T_mix, P, X_mix);
  std::vector<double> Y_b = mole_to_mass(normalize_fractions(X_b));

  std::size_t n_streams = streams.size();
  std::size_t n_species = num_species();

  std::vector<double> dT_ad_dY_mix(n_species, 0.0);
  std::vector<std::vector<double>> dYb_dY_mix(
      n_species, std::vector<double>(n_species, 0.0));
  std::vector<double> dYb_dT_mix(n_species, 0.0);
  const double eps_Y = 1e-6;
  const double eps_T = 1.0;

  for (std::size_t k = 0; k < n_species; ++k) {
    std::vector<double> Y_p = mix.Y_mix;
    Y_p[k] += eps_Y;
    std::vector<double> X_p = mass_to_mole(normalize_fractions(Y_p));
    auto [T_ad_p, dummy, dummy2, X_b_p] =
        adiabatic_T_equilibrium_and_jacobians(mix.T_mix, P, X_p);
    std::vector<double> Y_b_p = mole_to_mass(normalize_fractions(X_b_p));

    dT_ad_dY_mix[k] = (T_ad_p - T_ad) / eps_Y;
    for (std::size_t j = 0; j < n_species; ++j) {
      dYb_dY_mix[j][k] = (Y_b_p[j] - Y_b[j]) / eps_Y;
    }
  }

  // Also compute dY_b/dT_mix for equilibrium! Equilibrium products shift with
  // inlet enthalpy/temperature.
  auto [T_ad_p_T, dummy_T, dummy2_T, X_b_p_T] =
      adiabatic_T_equilibrium_and_jacobians(mix.T_mix + eps_T, P, X_mix);
  std::vector<double> Y_b_p_T = mole_to_mass(normalize_fractions(X_b_p_T));
  for (std::size_t j = 0; j < n_species; ++j) {
    dYb_dT_mix[j] = (Y_b_p_T[j] - Y_b[j]) / eps_T;
  }

  MixerResult res;
  res.T_mix = T_ad;
  res.Y_mix = Y_b;
  res.dT_mix_d_stream.resize(n_streams);
  res.dY_mix_d_stream.assign(n_species, std::vector<StreamJacobian>(n_streams));

  for (std::size_t i = 0; i < n_streams; ++i) {
    StreamJacobian &dT_stream = res.dT_mix_d_stream[i];
    dT_stream.d_Y.assign(n_species, 0.0);

    auto compute_chain = [&](double d_prop_T,
                             const std::vector<double> &d_prop_Y) {
      double sum = dT_ad_dT_mix * d_prop_T;
      for (std::size_t k = 0; k < n_species; ++k) {
        sum += dT_ad_dY_mix[k] * d_prop_Y[k];
      }
      return sum;
    };

    std::vector<double> dY_mix_dmdot(n_species), dY_mix_dT(n_species);
    for (std::size_t k = 0; k < n_species; ++k) {
      dY_mix_dmdot[k] = mix.dY_mix_d_stream[k][i].d_mdot;
      dY_mix_dT[k] = mix.dY_mix_d_stream[k][i].d_T;
    }

    dT_stream.d_mdot =
        compute_chain(mix.dT_mix_d_stream[i].d_mdot, dY_mix_dmdot);
    dT_stream.d_T = compute_chain(mix.dT_mix_d_stream[i].d_T, dY_mix_dT);
    for (std::size_t j = 0; j < n_species; ++j) {
      std::vector<double> dY_j(n_species);
      for (std::size_t k = 0; k < n_species; ++k)
        dY_j[k] = mix.dY_mix_d_stream[k][i].d_Y[j];
      dT_stream.d_Y[j] = compute_chain(mix.dT_mix_d_stream[i].d_Y[j], dY_j);
    }

    for (std::size_t tgt = 0; tgt < n_species; ++tgt) {
      StreamJacobian &dY_stream = res.dY_mix_d_stream[tgt][i];
      dY_stream.d_Y.assign(n_species, 0.0);

      auto compute_Y_chain = [&](double d_prop_T,
                                 const std::vector<double> &d_prop_Y) {
        double sum = dYb_dT_mix[tgt] * d_prop_T;
        for (std::size_t k = 0; k < n_species; ++k) {
          sum += dYb_dY_mix[tgt][k] * d_prop_Y[k];
        }
        return sum;
      };

      dY_stream.d_mdot =
          compute_Y_chain(mix.dT_mix_d_stream[i].d_mdot, dY_mix_dmdot);
      dY_stream.d_T = compute_Y_chain(mix.dT_mix_d_stream[i].d_T, dY_mix_dT);
      for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> dY_j(n_species);
        for (std::size_t k = 0; k < n_species; ++k)
          dY_j[k] = mix.dY_mix_d_stream[k][i].d_Y[j];
        dY_stream.d_Y[j] = compute_Y_chain(mix.dT_mix_d_stream[i].d_Y[j], dY_j);
      }
    }
  }
  return res;
}

} // namespace solver
} // namespace combaero
