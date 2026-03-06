#include "solver_interface.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "cooling_correlations.h"
#include "friction.h"
#include "heat_transfer.h"
#include "math_constants.h" // MSVC compatibility for M_PI, M_LN10, etc.
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
  if (tag == "haaland")
    return friction_and_jacobian_haaland(Re, e_D);
  if (tag == "serghides")
    return friction_and_jacobian_serghides(Re, e_D);
  if (tag == "colebrook")
    return friction_and_jacobian_colebrook(Re, e_D);
  if (tag == "petukhov")
    return friction_and_jacobian_petukhov(Re);
  return {{0.0, 0.0},
          CorrelationValidity::INVALID,
          "friction_and_jacobian: unknown tag '" + tag +
              "'. Must be haaland|serghides|colebrook|petukhov."};
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

} // namespace solver
} // namespace combaero
