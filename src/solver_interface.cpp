#include "../include/solver_interface.h"
#include "../include/area_change.h"
#include "../include/tee_junction.h"
#include "../include/composition.h"
#include "../include/transport.h"
#include "../include/combustion.h"
#include "../include/compressible.h"
#include "../include/cooling_correlations.h"
#include "../include/equilibrium.h"
#include "../include/friction.h"
#include "../include/heat_transfer.h"
#include "../include/registry.h"
#include "../include/stagnation.h"
#include "../include/thermo.h"
#include "../include/math_constants.h"
#include "../include/thermo_transport_data.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace combaero {
namespace solver {

namespace {
// Get global minimum temperature across all species from NASA9 thermo data
// Used to derive smooth floor limits for finite difference calculations
double get_global_thermo_T_min() {
  double T_min = std::numeric_limits<double>::max();
  for (const auto& nasa : nasa_coeffs) {
    if (!nasa.intervals.empty()) {
      T_min = std::min(T_min, nasa.intervals.front().T_min);
    }
  }
  return T_min;
}

// Smooth floor limit for temperature finite differences
// Set to 2/3 of global TMIN to provide large cushion (e.g., 200K for TMIN=300K)
const double THERMO_T_MIN = get_global_thermo_T_min();
const double T_SMOOTH_FLOOR_LIMIT = THERMO_T_MIN * 0.67;
} // anonymous namespace

// -----------------------------------------------------------------------------
// 1. Incompressible Flow Components
// -----------------------------------------------------------------------------

std::tuple<double, double, double> orifice_mdot_and_jacobian(double dP, double rho,
                                                              double Cd, double area,
                                                              double beta) {
  // Regularized orifice equation: mdot = Cd * E * A * f_reg(dP, rho)
  // where E = 1/sqrt(1 - beta^4) is the velocity-of-approach factor.
  // We use a smooth sqrt to avoid the infinite gradient at dP = 0.
  double E = 1.0;
  if (beta > 0.0 && beta < 1.0) {
    double b2 = beta * beta;
    E = 1.0 / std::sqrt(1.0 - b2 * b2);
  }

  /// Smooth clamp: rho_eff ~ max(rho, L)
  constexpr double L = 1e-9;
  constexpr double delta = 1e-4;
  double diff = rho - L;
  double sq_dist = std::sqrt(diff * diff + delta);
  double rho_eff = 0.5 * (rho + L + sq_dist);
  double d_rho_eff_d_rho = 0.5 * (1.0 + diff / sq_dist);

  // Pressure regularization
  constexpr double eps2 = 1.0; // Pa^2 (regularisation scale eps = 1 Pa)
  double dP2 = dP * dP;
  double U = dP2 + eps2;
  double U_pow_025 = std::sqrt(std::sqrt(U));
  double U_pow_125 = U * U_pow_025;

  // Mass Flow Calculation
  // mdot = (Cd * E * A * sqrt(2 * rho_eff)) * (dP / U^0.25)
  double geom_const = Cd * E * area;
  double sqrt_2rho = std::sqrt(2.0 * rho_eff);
  double mdot = geom_const * sqrt_2rho * (dP / U_pow_025);

  // Jacobian w.r.t dP
  // d(mdot)/d(dP) = constant * (0.5 * dP^2 + eps^2) / (dP^2 + eps^2)^1.25
  double d_mdot_dp = geom_const * sqrt_2rho * (0.5 * dP2 + eps2) / U_pow_125;

  // Jacobian w.r.t rho
  // d(mdot)/d(rho) = (mdot / (2 * rho_eff)) * d(rho_eff)/d(rho)
  double d_mdot_rho = (mdot / (2.0 * rho_eff)) * d_rho_eff_d_rho;

  return {mdot, d_mdot_dp, d_mdot_rho};
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

  // Central finite difference in Re (holding f fixed), consistent with
  // Gnielinski Jacobian treatment.
  const double eps = std::max(1.0, std::abs(Re) * 1e-6);
  double Nu_plus = nusselt_petukhov(Re + eps, Pr, f);
  double Nu_minus = nusselt_petukhov(Re - eps, Pr, f);
  double dNu_dRe = (Nu_plus - Nu_minus) / (2.0 * eps);
  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_haaland(double Re, double e_D) {
  // Explicit Haaland equation approximation of friction factor f (Darcy)
  // 1/sqrt(f) = coeff_outer * log10( (e_D/coeff_roughness)^coeff_exponent +
  // coeff_reynolds/Re )

  // Regularized Re to avoid singularity at Re=0 and sync with base friction evaluation
  const double re_eps = 64.0;
  double Re_eff = std::sqrt(Re * Re + re_eps * re_eps);

  // Isolate core logarithmic argument
  double a = std::pow(e_D / haaland::coeff_roughness, haaland::coeff_exponent);
  double b = haaland::coeff_reynolds / Re_eff;
  double arg = a + b;

  // 1/sqrt(f) for turbulent branch
  double inv_sqrt_f = haaland::coeff_outer * std::log10(arg);

  // Analytical derivative via chain rule for turbulent branch:
  // Let U = inv_sqrt_f = coeff_outer * log10(arg) = (coeff_outer / ln(10)) * ln(arg)
  double dU_darg = haaland::coeff_outer / (M_LN10 * arg);
  // d(arg)/dRe_eff = -coeff_reynolds / Re_eff^2
  double darg_dRe_eff = -haaland::coeff_reynolds / (Re_eff * Re_eff);
  double dRe_eff_dRe = Re / Re_eff;
  double darg_dRe = darg_dRe_eff * dRe_eff_dRe;
  double dU_dRe = dU_darg * darg_dRe;

  // Laminar branch
  double f_lam = 64.0 / Re_eff;
  double df_lam_dRe = (-64.0 / (Re_eff * Re_eff)) * dRe_eff_dRe;

  // Turbulent branch (Haaland)
  double f_turb = 1.0 / (inv_sqrt_f * inv_sqrt_f);
  double df_turb_dU = -2.0 / (inv_sqrt_f * inv_sqrt_f * inv_sqrt_f);
  double df_turb_dRe = df_turb_dU * dU_dRe;

  // Blending function
  double w_arg = (Re_eff - 2300.0) / 400.0;
  double w = 0.5 * (1.0 + std::tanh(w_arg));
  double dw_dRe_eff = 0.5 * (1.0 - std::tanh(w_arg) * std::tanh(w_arg)) / 400.0;
  double dw_dRe = dw_dRe_eff * dRe_eff_dRe;

  // Blended friction
  double f = (1.0 - w) * f_lam + w * f_turb;

  // Blended Jacobian (Product Rule)
  double jacobian = (1.0 - w) * df_lam_dRe + w * df_turb_dRe + dw_dRe * (f_turb - f_lam);

  return {{f, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_serghides(double Re, double e_D) {
  const double re_eps = 64.0;
  double Re_eff = std::sqrt(Re * Re + re_eps * re_eps);

  // Serghides Steffensen form: 1/sqrt(f) = A - (B-A)^2 / (C - 2B + A)
  // A = -2 log10(e_D/3.7 + 12/Re)
  // B = -2 log10(e_D/3.7 + 2.51*A/Re)
  // C = -2 log10(e_D/3.7 + 2.51*B/Re)
  const double r = serghides::coeff_roughness;
  const double cA = serghides::coeff_reynolds_A;
  const double cBC = serghides::coeff_reynolds_BC;
  const double ln10 = M_LN10;

  double argA = e_D / r + cA / Re_eff;
  double A = -2.0 * std::log10(argA);

  double argB = e_D / r + cBC * A / Re_eff;
  double B = -2.0 * std::log10(argB);

  double argC = e_D / r + cBC * B / Re_eff;
  double C = -2.0 * std::log10(argC);

  double denom = C - 2.0 * B + A;
  double inv_sqrt_f =
      (std::abs(denom) < 1e-30) ? C : A - (B - A) * (B - A) / denom;

  // Analytical derivative via chain rule through A, B, C, inv_sqrt_f, f.
  // dA/dRe = (2 / (ln10 * argA)) * (cA / Re^2)
  // dA/dRe = dA/dRe_eff * dRe_eff/dRe
  double dRe_eff_dRe = Re / Re_eff;
  double dA_dRe = (2.0 / (ln10 * argA)) * (cA / (Re_eff * Re_eff)) * dRe_eff_dRe;

  // dB/dRe = dB/d(argB) * d(argB)/dRe
  // d(argB)/dRe = cBC * d(A/Re_eff)/dRe = cBC * (dA/dRe * Re_eff - A * dRe_eff/dRe) / Re_eff^2
  double dargB_dRe = cBC * (dA_dRe * Re_eff - A * dRe_eff_dRe) / (Re_eff * Re_eff);
  double dB_dRe = (-2.0 / (ln10 * argB)) * dargB_dRe;

  double dargC_dRe = cBC * (dB_dRe * Re_eff - B * dRe_eff_dRe) / (Re_eff * Re_eff);
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
  // Laminar branch
  double f_lam = 64.0 / Re_eff;
  double df_lam_dRe = (-64.0 / (Re_eff * Re_eff)) * dRe_eff_dRe;

  // Turbulent branch (Serghides)
  double f_turb = 1.0 / (inv_sqrt_f * inv_sqrt_f);
  double df_turb_dRe = (-2.0 / (inv_sqrt_f * inv_sqrt_f * inv_sqrt_f)) * d_inv_sqrt_f_dRe;

  // Blending function
  double w_arg = (Re_eff - 2300.0) / 400.0;
  double w = 0.5 * (1.0 + std::tanh(w_arg));
  double dw_dRe_eff = 0.5 * (1.0 - std::tanh(w_arg) * std::tanh(w_arg)) / 400.0;
  double dw_dRe = dw_dRe_eff * dRe_eff_dRe;

  // Blended friction
  double f = (1.0 - w) * f_lam + w * f_turb;

  // Blended Jacobian (Product Rule)
  double jacobian = (1.0 - w) * df_lam_dRe + w * df_turb_dRe + dw_dRe * (f_turb - f_lam);

  return {{f, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_colebrook(double Re, double e_D) {
  const double re_eps = 64.0;
  double Re_eff = std::sqrt(Re * Re + re_eps * re_eps);

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
      e_D / serghides::coeff_roughness + serghides::coeff_reynolds_BC * x / Re_eff;

  double dF_dx = 1.0 + 2.0 * serghides::coeff_reynolds_BC / (ln10 * arg * Re_eff);
  // dx/dRe = -(dF/dRe) / (dF/dx)
  // dF/dRe = dF/dRe_eff * dRe_eff/dRe
  // dF/dRe_eff = 2/(ln10 * arg) * (-t_BC * x / Re_eff^2)
  double dF_dRe_eff = (-2.0 / (ln10 * arg)) * (serghides::coeff_reynolds_BC * x / (Re_eff * Re_eff));
  double dx_dRe = -(dF_dRe_eff * (Re / Re_eff)) / dF_dx;
  double jacobian = (-2.0 / (x * x * x)) * dx_dRe;

  return {{f_val, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_petukhov(double Re) {
  // f = (coeff_a * ln(Re) - coeff_b)^(-2)
  double x = petukhov::coeff_a * std::log(Re) - petukhov::coeff_b;
  double f = 1.0 / (x * x);

  // df/dRe = -2 * x^(-3) * dx/dRe
  // dx/dRe = coeff_a / Re
  double dx_dRe = petukhov::coeff_a / Re;
  double jacobian = (-2.0 / (x * x * x)) * dx_dRe;

  CorrelationValidity validity =
      (Re >= 3000.0) ? CorrelationValidity::VALID : CorrelationValidity::EXTRAPOLATED;
  return {{f, jacobian}, validity, ""};
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
  // Robust clamping with leaky gradient to avoid solver divergence into negative space
  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));
  double P_eff = (P > 1.0) ? P : (1.0 + 0.1 * (P - 1.0));
  double dTeff_dT = (T > 10.0) ? 1.0 : 0.1;
  double dPeff_dP = (P > 1.0) ? 1.0 : 0.1;

  // Calculate base density natively
  double rho = combaero::density(T_eff, P_eff, X);

  // Ideal gas law: rho = P * W / (R * T)
  // d_rho_dT = d(rho)/d(T_eff) * d(T_eff)/dT
  double d_rho_d_T = (-rho / T_eff) * dTeff_dT;

  // d_rho_dP = d(rho)/d(P_eff) * d(P_eff)/dP
  double d_rho_d_P = (rho / P_eff) * dPeff_dP;

  return {rho, d_rho_d_T, d_rho_d_P};
}

std::tuple<double, double> enthalpy_and_jacobian(double T,
                                                  const std::vector<double> &X) {
  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));
  double dTeff_dT = (T > 10.0) ? 1.0 : 0.1;

  // h(T) = base mass enthalpy
  double h_val = combaero::h_mass(T_eff, X);

  // Analytical derivative of enthalpy with respect to temperature is Cp
  double d_h_d_T = combaero::cp_mass(T_eff, X) * dTeff_dT;

  return {h_val, d_h_d_T};
}

std::tuple<double, double, double>
viscosity_and_jacobians(double T, double P, const std::vector<double> &X) {

  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));
  double P_eff = (P > 1.0) ? P : (1.0 + 0.1 * (P - 1.0));
  double dTeff_dT = (T > 10.0) ? 1.0 : 0.1;
  double dPeff_dP = (P > 1.0) ? 1.0 : 0.1;

  // Base viscosity
  double mu = combaero::viscosity(T_eff, P_eff, X);

  // 1. Temperature Derivative
  double dT_p = 1e-3; // 1 mK perturbation
  double mu_plus = combaero::viscosity(T_eff + dT_p, P_eff, X);
  double mu_minus = combaero::viscosity(T_eff - dT_p, P_eff, X);
  double d_mu_d_T = ((mu_plus - mu_minus) / (2.0 * dT_p)) * dTeff_dT;

  // 2. Pressure Derivative
  double dP_p = 1.0; // 1 Pascal perturbation
  double mu_P_plus = combaero::viscosity(T_eff, P_eff + dP_p, X);
  double mu_P_minus = combaero::viscosity(T_eff, P_eff - dP_p, X);
  double d_mu_d_P = ((mu_P_plus - mu_P_minus) / (2.0 * dP_p)) * dPeff_dP;

  return {mu, d_mu_d_T, d_mu_d_P};
}

std::tuple<double, double>
mach_number_and_jacobian_v(double v, double T, const std::vector<double> &X) {
  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));

  double a = combaero::speed_of_sound(T_eff, X);
  if (a <= 1.0) a = 1.0; // Minimal fallback

  double M = v / a;
  // dM/dv = 1/a
  double dM_dv = 1.0 / a;

  return {M, dM_dv};
}

std::tuple<double, double>
T_adiabatic_wall_and_jacobian_v(double T_static, double v, double T, double P,
                                 const std::vector<double> &X, bool turbulent) {
  double T_s_eff = (T_static > 10.0) ? T_static : (10.0 + 0.1 * (T_static - 10.0));
  double v_eff = (v >= 0.0) ? v : (0.1 * v);
  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));
  double P_eff = (P > 1.0) ? P : (1.0 + 0.1 * (P - 1.0));

  double Pr = combaero::prandtl(T_eff, P_eff, X);
  double r = turbulent ? std::cbrt(Pr) : std::sqrt(Pr);

  double mw_g = combaero::mwmix(X);
  double cp_val = combaero::cp(T_eff, X) / mw_g * 1000.0;

  double T_aw = T_s_eff + r * 0.5 * v_eff * v_eff / cp_val;
  // Use v_eff for derivative to keep it leaky
  double dT_aw_dv = r * v_eff / cp_val;

  return {T_aw, dT_aw_dv};
}

std::tuple<double, double>
T0_from_static_and_jacobian_M(double T, double M,
                               const std::vector<double> &X) {
  // Enforce non-negative M (T0_from_static requires M >= 0).
  // Use a smooth floor so the derivative remains informative near M = 0.
  // For solver purposes the useful range is M >= 0; negative values are
  // mapped to zero with a leaky 0.01 gradient.
  double M_eff = (M >= 0.0) ? M : 0.01 * M;
  double dMeff_dM = (M >= 0.0) ? 1.0 : 0.01;

  double T0 = T0_from_static(T, M_eff, X);

  // Analytical derivative via the definition h(T0, X) = h(T, X) + 0.5 * v^2
  // where v = M_eff * a(T, X).  Differentiating w.r.t. M_eff:
  //   cp_mass(T0) * dT0/dM_eff = M_eff * a(T)^2
  // => dT0/dM_eff = M_eff * a^2 / cp_mass(T0)
  double a = combaero::speed_of_sound(T, X);
  double cp_T0 = combaero::cp_mass(T0, X);
  double dT0_dMeff = (cp_T0 > 1.0) ? (M_eff * a * a / cp_T0) : 0.0;

  double dT0_dM = dT0_dMeff * dMeff_dM;

  return {T0, dT0_dM};
}

std::tuple<double, double>
P0_from_static_and_jacobian_M(double P, double T, double M,
                               const std::vector<double> &X) {
  // Enforce non-negative M with a leaky floor for differentiability.
  double M_eff = (M >= 0.0) ? M : 0.01 * M;
  double dMeff_dM = (M >= 0.0) ? 1.0 : 0.01;

  double P0 = P0_from_static(P, T, M_eff, X);

  // Chain rule: P0 = P * exp((s(T0) - s(T,P)) / R_specific)
  // where s(T0) = s(T0, X, P_REF) and T0 = T0_from_static(T, M_eff, X).
  //
  // dP0/dM_eff = P0 * d/dM_eff [ln(P0 / P)]
  //            = P0 * d/dT0 [(s(T0) - const) / R_specific] * dT0/dM_eff
  //            = P0 * [cp_mass(T0) / T0 / R_specific] * (dT0/dM_eff)
  //
  // The positive sign comes from: s increases with T (at const P_ref) for ideal
  // gas, so raising T0 increases the entropy difference, increasing P0/P.
  //
  // dT0/dM_eff (analytical, from T0 definition): M_eff * a(T)^2 / cp_mass(T0)
  double a = combaero::speed_of_sound(T, X);
  double T0_val = T0_from_static(T, M_eff, X);
  double cp_T0 = combaero::cp_mass(T0_val, X);
  double R_spec = combaero::specific_gas_constant(X); // J/(kg*K)

  double dT0_dMeff = (cp_T0 > 1.0) ? (M_eff * a * a / cp_T0) : 0.0;

  // dP0/dT0 = P0 * (+cp_mass(T0) / T0) / R_specific
  double dP0_dT0 = (T0_val > 1.0 && R_spec > 1.0)
                       ? P0 * (cp_T0 / T0_val) / R_spec
                       : 0.0;

  double dP0_dM = dP0_dT0 * dT0_dMeff * dMeff_dM;

  return {P0, dP0_dM};
}

// -----------------------------------------------------------------------------
// 6. Combustion Interfaces
// -----------------------------------------------------------------------------

std::tuple<double, double, std::vector<double>>
adiabatic_T_complete_and_jacobian_T(double T_in, double P,
                                     const std::vector<double> &X_in) {
  // Clamp to safe physical floor instead of throwing: the solver may
  // probe unphysical (T<=0, P<=0) intermediates during Newton steps.
  // Clamping keeps residuals finite with correct gradient direction so
  // the trust-region can step back toward physical territory.
  T_in = std::max(T_in, T_SMOOTH_FLOOR_LIMIT);
  P = std::max(P, 100.0);

  State in_state;
  in_state.set_TPX(T_in, P, X_in);
  State out_state = combaero::complete_combustion(in_state, true, true);
  double T_ad = out_state.T;
  std::vector<double> X_products = out_state.X;

  // Jacobian wrt T_in via Central Finite Difference with smooth floor
  double dT_eps = std::max(1e-3, T_in * 1e-4);

  // Smooth floor to avoid thermo data limits
  const double T_limit = T_SMOOTH_FLOOR_LIMIT;
  constexpr double delta_T = 1e-3;
  double T_target_minus = T_in - dT_eps;
  double diff_T = T_target_minus - T_limit;
  double sq_dist_T = std::sqrt(diff_T * diff_T + delta_T);
  double T_minus = 0.5 * (T_target_minus + T_limit + sq_dist_T);
  double T_plus = T_in + dT_eps;
  double actual_dT = T_plus - T_minus;

  State in_plus, in_minus;
  in_plus.set_TPX(T_plus, P, X_in);
  in_minus.set_TPX(T_minus, P, X_in);

  double T_ad_plus = combaero::complete_combustion(in_plus, true, true).T;
  double T_ad_minus = combaero::complete_combustion(in_minus, true, true).T;

  double dT_ad_dT_in = (T_ad_plus - T_ad_minus) / actual_dT;

  return {T_ad, dT_ad_dT_in, X_products};
}

std::tuple<double, double, double, std::vector<double>>
adiabatic_T_equilibrium_and_jacobians(double T_in, double P,
                                       const std::vector<double> &X_in) {
  T_in = std::max(T_in, T_SMOOTH_FLOOR_LIMIT);
  P = std::max(P, 100.0);

  State in_state;
  in_state.set_TPX(T_in, P, X_in);
  // Combustion Equilibrium
  EquilibriumResult eq_res = combaero::combustion_equilibrium(in_state);
  double T_ad = eq_res.state.T;
  std::vector<double> X_products = eq_res.state.X;

  // Jacobians via Central Finite Difference with smooth floor clamping
  double dT_eps = std::max(1e-3, T_in * 1e-4);
  double dP_eps = std::max(1.0, P * 1e-4);

  // wrt T_in: smooth floor to avoid thermo data limits
  const double T_limit = T_SMOOTH_FLOOR_LIMIT;
  constexpr double delta_T = 1e-3;
  double T_target_minus = T_in - dT_eps;
  double diff_T = T_target_minus - T_limit;
  double sq_dist_T = std::sqrt(diff_T * diff_T + delta_T);
  double T_minus = 0.5 * (T_target_minus + T_limit + sq_dist_T);
  double T_plus = T_in + dT_eps;
  double actual_dT = T_plus - T_minus;

  State in_T_plus, in_T_minus;
  in_T_plus.set_TPX(T_plus, P, X_in);
  in_T_minus.set_TPX(T_minus, P, X_in);

  double T_ad_T_plus = combaero::combustion_equilibrium(in_T_plus).state.T;
  double T_ad_T_minus = combaero::combustion_equilibrium(in_T_minus).state.T;
  double dT_ad_dT_in = (T_ad_T_plus - T_ad_T_minus) / actual_dT;

  // wrt P: smooth floor to avoid negative pressure
  constexpr double P_limit = 1.0;
  constexpr double delta_P = 1e-3;
  double P_target_minus = P - dP_eps;
  double diff_P = P_target_minus - P_limit;
  double sq_dist_P = std::sqrt(diff_P * diff_P + delta_P);
  double P_minus = 0.5 * (P_target_minus + P_limit + sq_dist_P);
  double P_plus = P + dP_eps;
  double actual_dP = P_plus - P_minus;

  State in_P_plus, in_P_minus;
  in_P_plus.set_TPX(T_in, P_plus, X_in);
  in_P_minus.set_TPX(T_in, P_minus, X_in);

  double T_ad_P_plus = combaero::combustion_equilibrium(in_P_plus).state.T;
  double T_ad_P_minus = combaero::combustion_equilibrium(in_P_minus).state.T;
  double dT_ad_dP = (T_ad_P_plus - T_ad_P_minus) / actual_dP;

  return {T_ad, dT_ad_dT_in, dT_ad_dP, X_products};
}

// -----------------------------------------------------------------------------
// 7. Stream-Based Network Solvers
// -----------------------------------------------------------------------------

MixerResult
mixer_from_streams_and_jacobians(const std::vector<Stream> &streams,
                                double Q, double fraction) {
  std::size_t n_streams = streams.size();
  std::size_t n_species = combaero::num_species();

  double mdot_tot = 0.0;
  std::vector<double> Y_mix(n_species, 0.0);
  double H_tot = 0.0;
  double P_total_tot = 0.0;

  std::vector<double> h_stream(n_streams, 0.0);
  std::vector<double> cp_stream(n_streams, 0.0);
  std::vector<std::vector<double>> hk_mass(n_streams,
                                            std::vector<double>(n_species, 0.0));

  for (std::size_t i = 0; i < n_streams; ++i) {
    mdot_tot += streams[i].m_dot;

    for (std::size_t k = 0; k < n_species; ++k) {
      Y_mix[k] += streams[i].m_dot * streams[i].Y[k];
      hk_mass[i][k] = J_per_mol_to_J_per_kg(combaero::h_species(k, streams[i].T),
                                             combaero::species_molar_mass(k));
      h_stream[i] += streams[i].Y[k] * hk_mass[i][k];

      double cpk_mass = J_per_mol_to_J_per_kg(combaero::cp_species(k, streams[i].T),
                                               combaero::species_molar_mass(k));
      cp_stream[i] += streams[i].Y[k] * cpk_mass;
    }
    H_tot += streams[i].m_dot * h_stream[i];
    P_total_tot += streams[i].m_dot * streams[i].P_total;
  }

  const double mdot_eps = 1e-12;
  if (std::abs(mdot_tot) > mdot_eps) {
    for (std::size_t k = 0; k < n_species; ++k) {
      Y_mix[k] /= mdot_tot;
    }
  } else {
    if (n_streams > 0) {
      Y_mix = streams[0].Y;
    }
  }

  // Compute adiabatic mass-averaged enthalpy
  double h_mix_base = (mdot_tot > 0.0) ? (H_tot / mdot_tot)
                                        : (n_streams > 0 ? h_stream[0] : 0.0);

  // Use effective mass flow for heat addition to avoid singularity when mdot_tot -> 0
  // mdot_eff is strictly positive and smooth through 0. Minimum bound is 1e-3 kg/s
  double mdot_eff = std::sqrt(mdot_tot * mdot_tot + 1e-6);

  // Apply energy transfer: delta_h = Q/mdot_eff + fraction * h_mix_base
  double delta_h = Q / mdot_eff + fraction * h_mix_base;

  double h_mix = h_mix_base + delta_h;
  double P_total_mix = (mdot_tot > 0.0)
                            ? (P_total_tot / mdot_tot)
                            : (n_streams > 0 ? streams[0].P_total : 0.0);

  std::vector<double> normalized_Y_mix = combaero::normalize_fractions(Y_mix);
  std::vector<double> X_mix = combaero::mass_to_mole(normalized_Y_mix);
  double T_guess = 300.0;
  if (n_streams > 0)
    T_guess = streams[0].T;
  double T_mix = combaero::calc_T_from_h_mass(h_mix, X_mix, T_guess);

  MixerResult res;
  res.T_mix = T_mix;
  res.P_total_mix = P_total_mix;
  res.Y_mix = normalized_Y_mix;
  res.dT_mix_d_stream.resize(n_streams);
  res.dP_total_mix_d_stream.resize(n_streams);
  res.dY_mix_d_stream.assign(n_species, std::vector<StreamJacobian>(n_streams));

  double cp_mix = 0.0;
  std::vector<double> hk_mix(n_species, 0.0);
  for (std::size_t k = 0; k < n_species; ++k) {
    hk_mix[k] =
        J_per_mol_to_J_per_kg(combaero::h_species(k, T_mix), combaero::species_molar_mass(k));
    double cpk_mix =
        J_per_mol_to_J_per_kg(combaero::cp_species(k, T_mix), combaero::species_molar_mass(k));
    cp_mix += Y_mix[k] * cpk_mix;
  }

  // Store dT/d(delta_h) for reference (not directly used in Jacobians)
  res.dT_mix_d_delta_h = (cp_mix > 0.0) ? (1.0 / cp_mix) : 0.0;

  for (std::size_t i = 0; i < n_streams; ++i) {
    StreamJacobian &dT_jac = res.dT_mix_d_stream[i];
    StreamJacobian &dP_jac = res.dP_total_mix_d_stream[i];
    dT_jac.d_Y.assign(n_species, 0.0);
    dP_jac.d_Y.assign(n_species, 0.0);

    if (std::abs(mdot_tot) > mdot_eps) {
      // P_total sensitivities
      dT_jac.d_P_total = 0.0;
      dP_jac.d_P_total = streams[i].m_dot / mdot_tot;
      dP_jac.d_mdot = (streams[i].P_total - P_total_mix) / mdot_tot;
      dP_jac.d_T = 0.0;

      // d(T_mix)/d(T_i)
      // Base contribution
      double dT_base_dT = (streams[i].m_dot * cp_stream[i]) / (mdot_tot * cp_mix);

      // Fraction contribution: d(fraction * h_mix_base)/d(T_i)
      // = fraction * (mdot_i / mdot_tot) * cp_i
      double dT_frac_dT = (fraction * streams[i].m_dot * cp_stream[i]) / (mdot_tot * cp_mix);

      dT_jac.d_T = dT_base_dT + dT_frac_dT;

      // d(T_mix)/d(m_dot_i)
      // Base contribution: d(h_mix_base)/d(mdot_i)
      double h_diff = (h_stream[i] - h_mix_base);
      double y_sum = 0.0;
      for (std::size_t k = 0; k < n_species; ++k) {
        y_sum += hk_mix[k] * (streams[i].Y[k] - Y_mix[k]);
      }
      double dh_base_dmdot = (h_diff - y_sum) / mdot_tot;

      // Q contribution: d(Q/mdot_eff)/d(mdot_i) = -Q * mdot_tot / (mdot_eff^3)
      double mdot_eff_3 = mdot_eff * mdot_eff * mdot_eff;
      double dh_Q_dmdot = -Q * mdot_tot / mdot_eff_3;

      // Fraction contribution: d(fraction * h_mix_base)/d(mdot_i)
      // = fraction * d(h_mix_base)/d(mdot_i)
      double dh_frac_dmdot = fraction * dh_base_dmdot;

      // Total: dT/dmdot = (1/cp_mix) * (dh_base + dh_Q + dh_frac) / dmdot
      dT_jac.d_mdot = (dh_base_dmdot + dh_Q_dmdot + dh_frac_dmdot) / cp_mix;

      // d(T_mix)/d(Y_i,k)
      for (std::size_t k = 0; k < n_species; ++k) {
        // Base contribution
        double dT_base_dY = streams[i].m_dot * (hk_mass[i][k] - hk_mix[k] + h_mix_base) / (mdot_tot * cp_mix);

        // Fraction contribution: d(fraction * h_mix_base)/d(Y_i,k)
        // = fraction * (mdot_i / mdot_tot) * (hk_i - hk_mix + h_mix_base)
        double dT_frac_dY = fraction * streams[i].m_dot * (hk_mass[i][k] - hk_mix[k] + h_mix_base) / (mdot_tot * cp_mix);

        dT_jac.d_Y[k] = dT_base_dY + dT_frac_dY;
      }

      // Y_mix sensitivities
      for (std::size_t k = 0; k < n_species; ++k) {
        StreamJacobian &dY_jac_k = res.dY_mix_d_stream[k][i];
        dY_jac_k.d_Y.resize(n_species);
        dY_jac_k.d_mdot = (streams[i].Y[k] - Y_mix[k]) / mdot_tot;
        dY_jac_k.d_T = 0.0;
        dY_jac_k.d_P_total = 0.0;
        for (std::size_t j = 0; j < n_species; ++j) {
            double delta = (k == j) ? 1.0 : 0.0;
            dY_jac_k.d_Y[j] = (streams[i].m_dot / mdot_tot) * (delta - Y_mix[k]);
        }
      }
    } else {
      dT_jac.d_mdot = 0.0;
      dT_jac.d_T = 0.0;
      dT_jac.d_P_total = 0.0;
      dP_jac.d_mdot = 0.0;
      dP_jac.d_T = 0.0;
      dP_jac.d_P_total = 0.0;
      for (std::size_t k = 0; k < n_species; ++k) {
        StreamJacobian &dY_jac_k = res.dY_mix_d_stream[k][i];
        dY_jac_k.d_mdot = 0.0;
        dY_jac_k.d_T = 0.0;
        dY_jac_k.d_P_total = 0.0;
        dY_jac_k.d_Y.assign(n_species, 0.0);
      }
    }
  }

  return res;
}

MixerResult adiabatic_T_complete_and_jacobian_T_from_streams(
    const std::vector<Stream> &streams, double P, double Q, double fraction) {
  // First mix adiabatically to get base state
  MixerResult mix = mixer_from_streams_and_jacobians(streams);

  // Compute total enthalpy flow for fraction calculation
  double mdot_tot = 0.0;
  double H_tot = 0.0;
  std::size_t n_species = combaero::num_species();
  for (const auto &s : streams) {
    mdot_tot += s.m_dot;
    double h_stream = 0.0;
    for (std::size_t k = 0; k < n_species; ++k) {
      double hk = J_per_mol_to_J_per_kg(combaero::h_species(k, s.T),
                                        combaero::species_molar_mass(k));
      h_stream += s.Y[k] * hk;
    }
    H_tot += s.m_dot * h_stream;
  }
  double h_mix_base = (mdot_tot > 1e-12) ? (H_tot / mdot_tot) : 0.0;

  // Use effective mass flow for heat addition to avoid singularity when mdot_tot -> 0
  double mdot_eff = std::sqrt(mdot_tot * mdot_tot + 1e-6);

  // Compute delta_h from Q and fraction
  double delta_h = Q / mdot_eff + fraction * h_mix_base;

  std::vector<double> X_mix = combaero::mass_to_mole(combaero::normalize_fractions(mix.Y_mix));
  auto complete_res =
      adiabatic_T_complete_and_jacobian_T(mix.T_mix, P, X_mix);
  double T_ad = std::get<0>(complete_res);
  double dT_ad_dT_mix = std::get<1>(complete_res);
  std::vector<double> X_b = std::get<2>(complete_res);
  std::vector<double> Y_b = combaero::mole_to_mass(combaero::normalize_fractions(X_b));

  // Post-combustion enthalpy shift
  double dT_out_dT_ad = 1.0;
  double cp_out = combaero::cp_mass(T_ad, X_b);
  double T_out = T_ad;
  if (std::abs(delta_h) > 0.0) {
    double h_ad = combaero::h_mass(T_ad, X_b);
    T_out = combaero::calc_T_from_h_mass(h_ad + delta_h, X_b, T_ad);
    cp_out = combaero::cp_mass(T_out, X_b);
    double cp_ad = combaero::cp_mass(T_ad, X_b);
    dT_out_dT_ad = cp_ad / cp_out;
  }
  T_ad = T_out;

  std::size_t n_streams = streams.size();

  std::vector<double> dT_ad_dY_mix(n_species, 0.0);
  std::vector<std::vector<double>> dYb_dY_mix(
      n_species, std::vector<double>(n_species, 0.0));
  const double eps_Y = 1e-6;

  for (std::size_t k = 0; k < n_species; ++k) {
    std::vector<double> Y_p = mix.Y_mix;
    Y_p[k] += eps_Y;
    std::vector<double> X_p = combaero::mass_to_mole(combaero::normalize_fractions(Y_p));
    auto cp_res =
        adiabatic_T_complete_and_jacobian_T(mix.T_mix, P, X_p);
    double T_ad_p = std::get<0>(cp_res);
    std::vector<double> X_b_p = std::get<2>(cp_res);
    std::vector<double> Y_b_p = combaero::mole_to_mass(combaero::normalize_fractions(X_b_p));

    dT_ad_dY_mix[k] = (T_ad_p - T_ad) / eps_Y;
    for (std::size_t j = 0; j < n_species; ++j) {
      dYb_dY_mix[j][k] = (Y_b_p[j] - Y_b[j]) / eps_Y;
    }
  }

  std::vector<double> dYb_dT_mix(n_species, 0.0);

  // Scale adiabatic derivatives by dT_out/dT_ad for post-combustion shift
  dT_ad_dT_mix *= dT_out_dT_ad;
  for (std::size_t k = 0; k < n_species; ++k) {
    dT_ad_dY_mix[k] *= dT_out_dT_ad;
  }

  MixerResult res;
  res.T_mix = T_ad;
  res.P_total_mix = mix.P_total_mix;
  res.Y_mix = Y_b;
  res.dT_mix_d_delta_h = (cp_out > 0.0) ? (1.0 / cp_out) : 0.0;
  res.dT_mix_d_stream.resize(n_streams);
  res.dP_total_mix_d_stream = mix.dP_total_mix_d_stream;
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

    std::vector<double> dY_mix_dmdot(n_species), dY_mix_dT(n_species),
        dY_mix_dP_total(n_species);
    for (std::size_t k = 0; k < n_species; ++k) {
      dY_mix_dmdot[k] = mix.dY_mix_d_stream[k][i].d_mdot;
      dY_mix_dT[k] = mix.dY_mix_d_stream[k][i].d_T;
      dY_mix_dP_total[k] = mix.dY_mix_d_stream[k][i].d_P_total;
    }

    dT_stream.d_mdot =
        compute_chain(mix.dT_mix_d_stream[i].d_mdot, dY_mix_dmdot);
    dT_stream.d_T = compute_chain(mix.dT_mix_d_stream[i].d_T, dY_mix_dT);
    dT_stream.d_P_total =
        compute_chain(mix.dT_mix_d_stream[i].d_P_total, dY_mix_dP_total);

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
      dY_stream.d_P_total =
          compute_Y_chain(mix.dT_mix_d_stream[i].d_P_total, dY_mix_dP_total);

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
    const std::vector<Stream> &streams, double P, double Q, double fraction) {
  // First mix adiabatically to get base state
  MixerResult mix = mixer_from_streams_and_jacobians(streams);

  // Compute total enthalpy flow for fraction calculation
  double mdot_tot = 0.0;
  double H_tot = 0.0;
  std::size_t n_species = combaero::num_species();
  for (const auto &s : streams) {
    mdot_tot += s.m_dot;
    double h_stream = 0.0;
    for (std::size_t k = 0; k < n_species; ++k) {
      double hk = J_per_mol_to_J_per_kg(combaero::h_species(k, s.T),
                                        combaero::species_molar_mass(k));
      h_stream += s.Y[k] * hk;
    }
    H_tot += s.m_dot * h_stream;
  }
  double h_mix_base = (mdot_tot > 1e-12) ? (H_tot / mdot_tot) : 0.0;

  // Use effective mass flow for heat addition to avoid singularity when mdot_tot -> 0
  double mdot_eff = std::sqrt(mdot_tot * mdot_tot + 1e-6);

  // Compute delta_h from Q and fraction
  double delta_h = Q / mdot_eff + fraction * h_mix_base;

  std::vector<double> X_mix = combaero::mass_to_mole(combaero::normalize_fractions(mix.Y_mix));
  auto eq_base_res =
      adiabatic_T_equilibrium_and_jacobians(mix.T_mix, P, X_mix);
  double T_ad = std::get<0>(eq_base_res);
  double dT_ad_dT_mix = std::get<1>(eq_base_res);
  std::vector<double> X_b = std::get<3>(eq_base_res);
  std::vector<double> Y_b = combaero::mole_to_mass(combaero::normalize_fractions(X_b));

  // Post-combustion enthalpy shift
  double dT_out_dT_ad = 1.0;
  double cp_out = combaero::cp_mass(T_ad, X_b);
  double T_out = T_ad;
  if (std::abs(delta_h) > 0.0) {
    double h_ad = combaero::h_mass(T_ad, X_b);
    T_out = combaero::calc_T_from_h_mass(h_ad + delta_h, X_b, T_ad);
    cp_out = combaero::cp_mass(T_out, X_b);
    double cp_ad = combaero::cp_mass(T_ad, X_b);
    dT_out_dT_ad = cp_ad / cp_out;
  }
  T_ad = T_out;

  std::size_t n_streams = streams.size();

  std::vector<double> dT_ad_dY_mix(n_species, 0.0);
  std::vector<std::vector<double>> dYb_dY_mix(
      n_species, std::vector<double>(n_species, 0.0));
  std::vector<double> dYb_dT_mix(n_species, 0.0);
  const double eps_Y = 1e-4;
  const double eps_T = 2.0;

  for (std::size_t k = 0; k < n_species; ++k) {
    std::vector<double> Y_p = mix.Y_mix;
    Y_p[k] += eps_Y;
    std::vector<double> X_p = combaero::mass_to_mole(combaero::normalize_fractions(Y_p));
    auto eq_p_res =
        adiabatic_T_equilibrium_and_jacobians(mix.T_mix, P, X_p);
    double T_ad_p = std::get<0>(eq_p_res);
    std::vector<double> X_b_p = std::get<3>(eq_p_res);
    std::vector<double> Y_b_p = combaero::mole_to_mass(combaero::normalize_fractions(X_b_p));

    dT_ad_dY_mix[k] = (T_ad_p - T_ad) / eps_Y;
    for (std::size_t j = 0; j < n_species; ++j) {
      dYb_dY_mix[j][k] = (Y_b_p[j] - Y_b[j]) / eps_Y;
    }
  }

  // dY_b/dT_mix: equilibrium products shift with inlet enthalpy/temperature
  auto eq_T_res =
      adiabatic_T_equilibrium_and_jacobians(mix.T_mix + eps_T, P, X_mix);
  std::vector<double> X_b_p_T = std::get<3>(eq_T_res);
  std::vector<double> Y_b_p_T = combaero::mole_to_mass(combaero::normalize_fractions(X_b_p_T));
  for (std::size_t j = 0; j < n_species; ++j) {
    dYb_dT_mix[j] = (Y_b_p_T[j] - Y_b[j]) / eps_T;
  }

  // Scale adiabatic derivatives by dT_out/dT_ad for post-combustion shift
  dT_ad_dT_mix *= dT_out_dT_ad;
  for (std::size_t k = 0; k < n_species; ++k) {
    dT_ad_dY_mix[k] *= dT_out_dT_ad;
  }

  MixerResult res;
  res.T_mix = T_ad;
  res.P_total_mix = mix.P_total_mix;
  res.Y_mix = Y_b;
  res.dT_mix_d_delta_h = (cp_out > 0.0) ? (1.0 / cp_out) : 0.0;
  res.dT_mix_d_stream.resize(n_streams);
  res.dP_total_mix_d_stream = mix.dP_total_mix_d_stream;
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

    std::vector<double> dY_mix_dmdot(n_species), dY_mix_dT(n_species),
        dY_mix_dP_total(n_species);
    for (std::size_t k = 0; k < n_species; ++k) {
      dY_mix_dmdot[k] = mix.dY_mix_d_stream[k][i].d_mdot;
      dY_mix_dT[k] = mix.dY_mix_d_stream[k][i].d_T;
      dY_mix_dP_total[k] = mix.dY_mix_d_stream[k][i].d_P_total;
    }

    dT_stream.d_mdot =
        compute_chain(mix.dT_mix_d_stream[i].d_mdot, dY_mix_dmdot);
    dT_stream.d_T = compute_chain(mix.dT_mix_d_stream[i].d_T, dY_mix_dT);
    dT_stream.d_P_total =
        compute_chain(mix.dT_mix_d_stream[i].d_P_total, dY_mix_dP_total);

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
      dY_stream.d_P_total =
          compute_Y_chain(mix.dT_mix_d_stream[i].d_P_total, dY_mix_dP_total);

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

MixerResult combustor_residuals_and_jacobians(
    const std::vector<Stream> &streams, double P, double Q, double fraction,
    const PressureLossCorrelation &pressure_loss, bool use_equilibrium) {

  // 1. Get unburned mixture state (T_in, P_in)
  MixerResult mix_in = mixer_from_streams_and_jacobians(streams);
  double T_in = mix_in.T_mix;

  // 2. Get burned state (T_ad, Y_products) including all stream jacobians
  MixerResult res = use_equilibrium
      ? adiabatic_T_equilibrium_and_jacobians_from_streams(streams, P, Q, fraction)
      : adiabatic_T_complete_and_jacobian_T_from_streams(streams, P, Q, fraction);

  double T_ad = res.T_mix;
  double theta = (T_in > 0.0) ? (T_ad / T_in - 1.0) : 0.0;

  // 3. Form Context and evaluate Hook
  State reactant_state;
  std::vector<double> X_in = mass_to_mole(mix_in.Y_mix);
  reactant_state.set_TPX(T_in, P, X_in);

  // Compute X_products from Y_products for the context
  std::vector<double> X_products = mass_to_mole(res.Y_mix);

  double mdot_total = 0.0;
  for (const auto &s : streams)
    mdot_total += s.m_dot;

  PressureLossContext ctx{
      reactant_state,
      0.0, // phi
      T_ad,
      X_products,
      res.Y_mix, // Y_products
      theta,
      0.0, 0.0, // mdot fuel/air
      mdot_total,
  };

  auto loss_res = pressure_loss(ctx);
  double zeta = std::get<0>(loss_res);
  double dzeta_dtheta = std::get<1>(loss_res);

  // 4. Update Result with P_out
  double P_out = P * (1.0 - zeta);
  res.P_total_mix = P_out;

  // 5. Analytical Chain Rule for dP_out/d_stream
  // dP_out/dx = dP/dx * (1 - zeta) - P * dzeta/dx
  // dzeta/dx = dzeta/dtheta * dtheta/dx
  // dtheta/dx = d(T_ad/T_in - 1)/dx = (1/T_in) * dT_ad/dx - (T_ad/T_in^2) * dT_in/dx

  std::size_t n_streams = streams.size();
  std::size_t n_species = num_species();

  for (std::size_t i = 0; i < n_streams; ++i) {
      auto update_jac = [&](StreamJacobian& dP_out_jac,
                          const StreamJacobian& dT_ad_jac,
                          const StreamJacobian& dT_in_jac,
                          const StreamJacobian& dP_in_jac) {

          auto compute_dtheta = [&](double d_Tad, double d_Tin) {
              return (T_in > 0.0) ? (d_Tad / T_in - (T_ad / (T_in * T_in)) * d_Tin) : 0.0;
          };

          double dtheta_dmdot = compute_dtheta(dT_ad_jac.d_mdot, dT_in_jac.d_mdot);
          double dtheta_dT = compute_dtheta(dT_ad_jac.d_T, dT_in_jac.d_T);
          double dtheta_dP = compute_dtheta(dT_ad_jac.d_P_total, dT_in_jac.d_P_total);

          dP_out_jac.d_mdot = dP_in_jac.d_mdot * (1.0 - zeta) - P * dzeta_dtheta * dtheta_dmdot;
          dP_out_jac.d_T = dP_in_jac.d_T * (1.0 - zeta) - P * dzeta_dtheta * dtheta_dT;
          dP_out_jac.d_P_total = dP_in_jac.d_P_total * (1.0 - zeta) - P * dzeta_dtheta * dtheta_dP;

          for (std::size_t k = 0; k < n_species; ++k) {
              double dtheta_dY = compute_dtheta(dT_ad_jac.d_Y[k], dT_in_jac.d_Y[k]);
              dP_out_jac.d_Y[k] = dP_in_jac.d_Y[k] * (1.0 - zeta) - P * dzeta_dtheta * dtheta_dY;
          }
      };

      update_jac(res.dP_total_mix_d_stream[i],
                res.dT_mix_d_stream[i], // dT_ad
                mix_in.dT_mix_d_stream[i], // dT_in
                mix_in.dP_total_mix_d_stream[i]); // dP_in
  }

  return res;
}

OrificeResult orifice_residuals_and_jacobian(double m_dot, double P_total_up,
                                             double P_static_up, double T_up,
                                             const std::vector<double> &Y_up,
                                             double P_static_down, double Cd,
                                             double area, double beta) {
  (void)m_dot;
  double dP = P_total_up - P_static_down;
  const std::vector<double> X_up = combaero::mass_to_mole(combaero::normalize_fractions(Y_up));
  auto [rho, drho_dT, drho_dP] = density_and_jacobians(T_up, P_static_up, X_up);

  // Get mass flow and analytical Jacobians from low-level function
  auto [mdot_calc, dmdot_ddP, dmdot_drho] = orifice_mdot_and_jacobian(dP, rho, Cd, area, beta);

  OrificeResult res;
  res.m_dot_calc = mdot_calc;
  res.d_mdot_dP_total_up = dmdot_ddP;
  res.d_mdot_dP_static_down = -dmdot_ddP;

  // Chain rule: d(mdot)/d(P_static_up) = d(mdot)/d(rho) * d(rho)/d(P)
  res.d_mdot_dP_static_up = dmdot_drho * drho_dP;
  res.d_mdot_dT_up = dmdot_drho * drho_dT;

  // Analytical Y derivatives via chain rule: d(mdot)/d(Y[i]) = d(mdot)/d(rho) * d(rho)/d(Y[i])
  // Compute d(rho)/d(Y[i]) using finite differences (no analytical drho_dY available)
  res.d_mdot_dY_up.resize(Y_up.size(), 0.0);
  const double eps_Y = 1e-6;
  for (std::size_t i = 0; i < Y_up.size(); ++i) {
    std::vector<double> Y_plus = Y_up;
    Y_plus[i] += eps_Y;
    auto X_plus = combaero::mass_to_mole(combaero::normalize_fractions(Y_plus));
    double rho_plus = std::get<0>(density_and_jacobians(T_up, P_static_up, X_plus));
    double drho_dY = (rho_plus - rho) / eps_Y;
    res.d_mdot_dY_up[i] = dmdot_drho * drho_dY;
  }

  return res;
}

ChannelResult channel_residuals_and_jacobian(double m_dot, double P_total_up,
                                       double P_static_up, double T_up,
                                       const std::vector<double> &Y_up,
                                       double P_static_down, double L, double D,
                                       double roughness,
                                       const std::string &friction_model,
                                       double f_multiplier) {
  (void)P_total_up;
  (void)P_static_down;
  const std::vector<double> X_up = combaero::mass_to_mole(combaero::normalize_fractions(Y_up));
  auto [rho, drho_dT, drho_dP] = density_and_jacobians(T_up, P_static_up, X_up);
  auto [mu, dmu_dT, dmu_dP] = viscosity_and_jacobians(T_up, P_static_up, X_up);

  double area = 0.25 * M_PI * D * D;
  double v = m_dot / (rho * area);
  double abs_v = std::abs(v);
  double Re = rho * abs_v * D / mu;

  auto [f, df_dRe] = friction_and_jacobian(friction_model, Re, roughness / D).result;
  f *= f_multiplier;
  df_dRe *= f_multiplier;

  double dP_calc = f * (L / D) * (0.5 * rho * v * abs_v);

  ChannelResult res;
  res.dP_calc = dP_calc;

  double dRe_dmdot = D / (mu * area);
  res.d_dP_d_mdot = (L / D) * (1.0 / (rho * area * area)) * m_dot * (f + 0.5 * m_dot * df_dRe * dRe_dmdot);

  const double eps_P = 1.0;
  auto [rho_p, drho_dT_p, drho_dP_p] = density_and_jacobians(T_up, P_static_up + eps_P, X_up);
  auto [mu_p, dmu_dT_p, dmu_dP_p] = viscosity_and_jacobians(T_up, P_static_up + eps_P, X_up);
  auto [rho_m, drho_dT_m, drho_dP_m] = density_and_jacobians(T_up, P_static_up - eps_P, X_up);
  auto [mu_m, dmu_dT_m, dmu_dP_m] = viscosity_and_jacobians(T_up, P_static_up - eps_P, X_up);
  double v_p = m_dot / (rho_p * area);
  double v_m = m_dot / (rho_m * area);
  double Re_p = rho_p * std::abs(v_p) * D / mu_p;
  double Re_m = rho_m * std::abs(v_m) * D / mu_m;
  auto [f_p, df_dRe_p] = friction_and_jacobian(friction_model, Re_p, roughness / D).result;
  auto [f_m, df_dRe_m] = friction_and_jacobian(friction_model, Re_m, roughness / D).result;
  f_p *= f_multiplier;
  f_m *= f_multiplier;
  double dP_p = f_p * (L / D) * (0.5 * rho_p * v_p * std::abs(v_p));
  double dP_m = f_m * (L / D) * (0.5 * rho_m * v_m * std::abs(v_m));
  res.d_dP_dP_static_up = (dP_p - dP_m) / (2.0 * eps_P);

  const double eps_T = 1e-3;
  auto [rho_T_p, drho_dT_T_p, drho_dP_T_p] = density_and_jacobians(T_up + eps_T, P_static_up, X_up);
  auto [mu_T_p, dmu_dT_T_p, dmu_dP_T_p] = viscosity_and_jacobians(T_up + eps_T, P_static_up, X_up);
  auto [rho_T_m, drho_dT_T_m, drho_dP_T_m] = density_and_jacobians(T_up - eps_T, P_static_up, X_up);
  auto [mu_T_m, dmu_dT_T_m, dmu_dP_T_m] = viscosity_and_jacobians(T_up - eps_T, P_static_up, X_up);
  double v_T_p = m_dot / (rho_T_p * area);
  double v_T_m = m_dot / (rho_T_m * area);
  double Re_T_p = rho_T_p * std::abs(v_T_p) * D / mu_T_p;
  double Re_T_m = rho_T_m * std::abs(v_T_m) * D / mu_T_m;
  auto [f_T_p, df_dRe_T_p] = friction_and_jacobian(friction_model, Re_T_p, roughness / D).result;
  auto [f_T_m, df_dRe_T_m] = friction_and_jacobian(friction_model, Re_T_m, roughness / D).result;
  f_T_p *= f_multiplier;
  f_T_m *= f_multiplier;
  double dP_T_p = f_T_p * (L / D) * (0.5 * rho_T_p * v_T_p * std::abs(v_T_p));
  double dP_T_m = f_T_m * (L / D) * (0.5 * rho_T_m * v_T_m * std::abs(v_T_m));
  res.d_dP_dT_up = (dP_T_p - dP_T_m) / (2.0 * eps_T);

  res.d_dP_dY_up.resize(Y_up.size(), 0.0);
  const double eps_Y = 1e-6;
  for (std::size_t i = 0; i < Y_up.size(); ++i) {
    std::vector<double> Y_plus = Y_up;
    Y_plus[i] += eps_Y;
    auto X_plus = combaero::mass_to_mole(combaero::normalize_fractions(Y_plus));
    double rho_plus = std::get<0>(density_and_jacobians(T_up, P_static_up, X_plus));
    double mu_plus = std::get<0>(viscosity_and_jacobians(T_up, P_static_up, X_plus));
    double v_plus = m_dot / (rho_plus * area);
    double Re_plus = rho_plus * std::abs(v_plus) * D / mu_plus;
    auto [f_plus, dummy] = friction_and_jacobian(friction_model, Re_plus, roughness / D).result;
    f_plus *= f_multiplier;
    double dP_plus = f_plus * (L / D) * (0.5 * rho_plus * v_plus * std::abs(v_plus));
    res.d_dP_dY_up[i] = (dP_plus - dP_calc) / eps_Y;
  }

  return res;
}

MomentumChamberResult momentum_chamber_residual_and_jacobian(
    double P, double P_total, double m_dot, double T,
    const std::vector<double> &Y, double area) {
  // Momentum chamber: P_total = P_static + 0.5 * rho * v^2
  // where v = m_dot / (rho * A)
  // Residual: P_total - P_static - 0.5 * rho * v^2 = 0

  // Soft normalization to handle edge cases (all-zero Y)
  std::vector<double> Y_safe = Y;
  double sum_Y = 0.0;
  for (double y : Y) {
    sum_Y += y;
  }
  // Use softmax-like approach: if sum is too small, use standard air
  const double eps = 1e-10;
  if (sum_Y < eps) {
    // Use standard dry air composition as fallback (N2=0.767, O2=0.233 mass fractions)
    Y_safe = std::vector<double>(Y.size(), 0.0);
    std::size_t n2_idx = combaero::species_index_from_name("N2");
    std::size_t o2_idx = combaero::species_index_from_name("O2");
    Y_safe[n2_idx] = 0.767;
    Y_safe[o2_idx] = 0.233;
  }

  const std::vector<double> X = combaero::mass_to_mole(combaero::normalize_fractions(Y_safe));
  auto [rho, drho_dT, drho_dP] = density_and_jacobians(T, P, X);

  // Velocity from mass flow
  double v = m_dot / (rho * area);

  // Dynamic pressure
  double q_dynamic = 0.5 * rho * v * v;

  // Residual: P_total - (P + q_dynamic) = 0
  MomentumChamberResult res;
  res.residual = P_total - P - q_dynamic;

  // Analytical Jacobians
  // d(res)/d(P_total) = 1.0
  res.d_res_dP_total = 1.0;

  // d(res)/d(P) = -1.0 - d(q_dynamic)/d(P)
  // q_dynamic = 0.5 * m_dot^2 / (rho * A^2)
  // d(q_dynamic)/d(P) = d(q_dynamic)/d(rho) * d(rho)/d(P)
  //                   = -0.5 * m_dot^2 / (rho^2 * A^2) * drho_dP
  double dq_dP = -0.5 * m_dot * m_dot / (rho * rho * area * area) * drho_dP;
  res.d_res_dP = -1.0 - dq_dP;

  // d(res)/d(m_dot) = -d(q_dynamic)/d(m_dot)
  // q_dynamic = 0.5 * m_dot^2 / (rho * A^2)
  // d(q_dynamic)/d(m_dot) = m_dot / (rho * A^2)
  double dq_dmdot = m_dot / (rho * area * area);
  res.d_res_dmdot = -dq_dmdot;

  return res;
}

// -----------------------------------------------------------------------------
// Compressible Flow Components
// -----------------------------------------------------------------------------

std::tuple<double, double, double, double> orifice_compressible_mdot_and_jacobian(
    double T0, double P0, double P_back,
    const std::vector<double>& X,
    double Cd, double area, double beta) {

  // Compute effective area with velocity-of-approach factor
  double E = 1.0;
  if (beta > 0.0 && beta < 1.0) {
    double b2 = beta * beta;
    E = 1.0 / std::sqrt(1.0 - b2 * b2);
  }
  double A_eff = Cd * E * area;

  // Handle reverse flow (P_back > P0)
  bool reverse_flow = (P_back > P0);
  double T0_fwd = T0;  // Temperature same in both directions
  double P0_fwd = reverse_flow ? P_back : P0;
  double P_back_fwd = reverse_flow ? P0 : P_back;

  // Call nozzle_flow for forward direction
  auto sol = combaero::nozzle_flow(T0_fwd, P0_fwd, P_back_fwd, A_eff, X);
  double mdot = reverse_flow ? -sol.mdot : sol.mdot;

  // Compute Jacobians analytically using isentropic relations
  // G = rho * v = (P/RT) * sqrt(2*(h0 - h))
  // For isentropic flow: dT/dP = 1/(rho*cp), d_rho/dP = 1/a^2
  // dG/dP = 1/v * (M^2 - 1)

  double T_outlet = sol.outlet.T;
  double v_outlet = sol.v;
  double rho_outlet = sol.outlet.rho();
  double M_outlet = sol.M;
  double cp_outlet = combaero::cp_mass(T_outlet, X);
  double cp0 = combaero::cp_mass(T0, X);

  // Unchoked derivatives (subsonic branch)
  double dmdot_dP_back_sub = 0.0;
  if (v_outlet > 1e-6) {
    dmdot_dP_back_sub = A_eff * (1.0 / v_outlet) * (M_outlet * M_outlet - 1.0);
  }

  // Choked derivatives (constant mass flow)
  // mdot_choked = A_eff * G_crit
  // G_crit ~ P0 / sqrt(T0)
  double dmdot_dP_back_choked = 0.0;

  // Apply smooth transition (C1 continuous)
  double PR = P_back_fwd / P0_fwd;
  double PR_crit = combaero::critical_pressure_ratio(T0_fwd, P0_fwd, X);
  const double delta_PR = 0.01;
  double t = (PR - PR_crit) / delta_PR;
  double choke_factor = 1.0;
  if (t >= 1.0) choke_factor = 1.0;
  else if (t <= -1.0) choke_factor = 0.0;
  else choke_factor = 0.5 * (1.0 + t * (3.0 - t * t));

  double d_mdot_dP_back_fwd = dmdot_dP_back_sub * choke_factor + dmdot_dP_back_choked * (1.0 - choke_factor);

  // Upstream pressure derivative: P0 only affects PR and scales G linearly if choked
  // dmdot/dP0 = mdot/P0 - (Pb/P0) * dmdot/dPb
  double d_mdot_dP0_fwd = (P0_fwd > 1e-6) ? (sol.mdot / P0_fwd - PR * d_mdot_dP_back_fwd) : 0.0;

  // Upstream temperature derivative: mdot ~ 1/sqrt(T0) roughly
  // More precisely: dG/dT0 = rho * cp0 / v * [(1 - T/T0) - v^2/(cp*T0)]
  double d_mdot_dT0_fwd = 0.0;
  if (v_outlet > 1e-6 && T0 > 1e-6) {
    d_mdot_dT0_fwd = A_eff * rho_outlet * cp0 / v_outlet * ((1.0 - T_outlet / T0) - v_outlet * v_outlet / (cp_outlet * T0));
  } else if (T0 > 1e-6) {
    // Choked limit: dG/dT0 = -G / (2*T0)
    d_mdot_dT0_fwd = -sol.mdot / (2.0 * T0);
  }

  // Assign results
  double d_mdot_dP0 = d_mdot_dP0_fwd;
  double d_mdot_dP_back = d_mdot_dP_back_fwd;
  double d_mdot_dT0 = d_mdot_dT0_fwd;

  // For reverse flow, swap Jacobian signs appropriately.
  // mdot = -mdot_fwd(P_back, P0), so all three derivatives are negated,
  // and P0/P_back roles are swapped.
  if (reverse_flow) {
    std::swap(d_mdot_dP0, d_mdot_dP_back);
    d_mdot_dP0 = -d_mdot_dP0;
    d_mdot_dP_back = -d_mdot_dP_back;
    d_mdot_dT0 = -d_mdot_dT0;
  }

  return {mdot, d_mdot_dP0, d_mdot_dP_back, d_mdot_dT0};
}

OrificeResult orifice_compressible_residuals_and_jacobian(
    double m_dot, double P_total_up, double T_up,
    const std::vector<double>& Y_up,
    double P_static_down, double Cd, double area, double beta) {

  (void)m_dot;  // Not used in compressible orifice (computed from pressures)
  const std::vector<double> X_up = combaero::mass_to_mole(combaero::normalize_fractions(Y_up));

  // Get mass flow and Jacobians from low-level function
  auto [mdot_calc, dmdot_dP0, dmdot_dP_back, dmdot_dT0] =
      orifice_compressible_mdot_and_jacobian(T_up, P_total_up, P_static_down, X_up, Cd, area, beta);

  OrificeResult res;
  res.m_dot_calc = mdot_calc;
  res.d_mdot_dP_total_up = dmdot_dP0;
  res.d_mdot_dP_static_down = dmdot_dP_back;
  res.d_mdot_dP_static_up = 0.0;  // Compressible model uses total pressure, not static
  res.d_mdot_dT_up = dmdot_dT0;

  // Compute d(mdot)/d(Y[i]) using finite differences
  res.d_mdot_dY_up.resize(Y_up.size(), 0.0);
  const double eps_Y = 1e-6;
  for (std::size_t i = 0; i < Y_up.size(); ++i) {
    std::vector<double> Y_plus = Y_up;
    Y_plus[i] += eps_Y;
    auto X_plus = combaero::mass_to_mole(combaero::normalize_fractions(Y_plus));
    auto [mdot_plus, dummy1, dummy2, dummy3] =
        orifice_compressible_mdot_and_jacobian(T_up, P_total_up, P_static_down, X_plus, Cd, area, beta);
    res.d_mdot_dY_up[i] = (mdot_plus - mdot_calc) / eps_Y;
  }

  return res;
}

std::tuple<double, double, double, double> channel_compressible_mdot_and_jacobian(
    double T_in, double P_in, double u_in,
    const std::vector<double>& X,
    double L, double D, double roughness,
    const std::string& friction_model,
    double f_multiplier,
    bool compute_jacobians) {

  // Handle reverse flow by swapping direction
  bool reverse_flow = (u_in < 0.0);
  double u_fwd = std::abs(u_in);

  // Call fanno_channel_rough for forward direction
  auto sol = combaero::fanno_channel_rough(T_in, P_in, u_fwd, L, D, roughness, X, friction_model, f_multiplier);
  double dP = P_in - sol.outlet.P;

  // For reverse flow, negate dP to maintain sign convention
  if (reverse_flow) {
    dP = -dP;
  }

  double d_dP_dP_in = 0.0;
  double d_dP_dT_in = 0.0;
  double d_dP_du_in = 0.0;

  if (compute_jacobians) {
    // Compute Jacobians via finite differences
    const double eps_P = std::max(1e-6, std::abs(P_in) * 1e-6);
    const double eps_T = std::max(1e-6, std::abs(T_in) * 1e-6);
    const double eps_u = std::max(1e-6, std::abs(u_in) * 1e-6);

    // Jacobian w.r.t. P_in
    auto sol_P_plus = combaero::fanno_channel_rough(T_in, P_in + eps_P, u_fwd, L, D, roughness, X, friction_model, f_multiplier);
    auto sol_P_minus = combaero::fanno_channel_rough(T_in, P_in - eps_P, u_fwd, L, D, roughness, X, friction_model, f_multiplier);
    double dP_P_plus = (P_in + eps_P) - sol_P_plus.outlet.P;
    double dP_P_minus = (P_in - eps_P) - sol_P_minus.outlet.P;
    if (reverse_flow) dP_P_plus = -dP_P_plus;
    if (reverse_flow) dP_P_minus = -dP_P_minus;
    d_dP_dP_in = (dP_P_plus - dP_P_minus) / (2.0 * eps_P);

    // Jacobian w.r.t. T_in
    auto sol_T_plus = combaero::fanno_channel_rough(T_in + eps_T, P_in, u_fwd, L, D, roughness, X, friction_model, f_multiplier);
    auto sol_T_minus = combaero::fanno_channel_rough(T_in - eps_T, P_in, u_fwd, L, D, roughness, X, friction_model, f_multiplier);
    double dP_T_plus = P_in - sol_T_plus.outlet.P;
    double dP_T_minus = P_in - sol_T_minus.outlet.P;
    if (reverse_flow) dP_T_plus = -dP_T_plus;
    if (reverse_flow) dP_T_minus = -dP_T_minus;
    d_dP_dT_in = (dP_T_plus - dP_T_minus) / (2.0 * eps_T);

    // Jacobian w.r.t. u_in (note: u_in is signed, but we use u_fwd for Fanno)
    auto sol_u_plus = combaero::fanno_channel_rough(T_in, P_in, u_fwd + eps_u, L, D, roughness, X, friction_model, f_multiplier);
    auto sol_u_minus = combaero::fanno_channel_rough(T_in, P_in, std::max(1e-9, u_fwd - eps_u), L, D, roughness, X, friction_model, f_multiplier);
    double dP_u_plus = P_in - sol_u_plus.outlet.P;
    double dP_u_minus = P_in - sol_u_minus.outlet.P;
    if (reverse_flow) dP_u_plus = -dP_u_plus;
    if (reverse_flow) dP_u_minus = -dP_u_minus;
    double d_dP_du_fwd = (dP_u_plus - dP_u_minus) / (2.0 * eps_u);
    d_dP_du_in = reverse_flow ? -d_dP_du_fwd : d_dP_du_fwd;

    // Apply smooth transition near choked conditions
    if (sol.choked || sol.L_choke < L * 1.2) {
      // Near choking - Jacobians may need smoothing
      // For now, keep as-is since Fanno solver handles this internally
    }
  }

  return {dP, d_dP_dP_in, d_dP_dT_in, d_dP_du_in};
}

ChannelResult channel_compressible_residuals_and_jacobian(
    double m_dot, double P_total_up, double T_up,
    const std::vector<double>& Y_up,
    double P_static_down, double L, double D, double roughness,
    const std::string& friction_model,
    double f_multiplier) {

  (void)P_static_down;  // Computed from Fanno flow, not an input

  const std::vector<double> X_up = combaero::mass_to_mole(combaero::normalize_fractions(Y_up));

  // Compute density at (T_up, P_total_up) and derive velocity from continuity.
  // Uses P_total as proxy for P_static (low-Mach approximation).
  // The Jacobians below include the full chain rule through rho(T,P) -> u -> dP
  // so that Newton converges correctly even at moderate Mach numbers.
  auto [rho, drho_dT, drho_dP] = density_and_jacobians(T_up, P_total_up, X_up);
  double area = 0.25 * M_PI * D * D;
  double u_in = m_dot / (rho * area);

  // Get pressure drop and partial Jacobians w.r.t. (P_in, T_in, u_in)
  auto [dP_calc, d_dP_dP_in, d_dP_dT_in, d_dP_du_in] =
      channel_compressible_mdot_and_jacobian(T_up, P_total_up, u_in, X_up, L, D, roughness, friction_model, f_multiplier, true);

  ChannelResult res;
  res.dP_calc = dP_calc;

  // --- Full chain rule through u = m_dot / (rho(T,P) * area) ---
  // d(u)/d(m_dot) = 1 / (rho * area)
  // d(u)/d(P)     = -m_dot / (rho^2 * area) * drho/dP
  // d(u)/d(T)     = -m_dot / (rho^2 * area) * drho/dT
  double du_dmdot = 1.0 / (rho * area);
  double du_dP = -m_dot / (rho * rho * area) * drho_dP;
  double du_dT = -m_dot / (rho * rho * area) * drho_dT;

  // d(dP)/d(m_dot) = d(dP)/d(u) * d(u)/d(m_dot)
  res.d_dP_d_mdot = d_dP_du_in * du_dmdot;

  // d(dP)/d(P_total) = d(dP)/d(P_in) + d(dP)/d(u) * d(u)/d(P)
  res.d_dP_dP_static_up = d_dP_dP_in + d_dP_du_in * du_dP;

  // d(dP)/d(T) = d(dP)/d(T_in) + d(dP)/d(u) * d(u)/d(T)
  res.d_dP_dT_up = d_dP_dT_in + d_dP_du_in * du_dT;

  // Compute d(dP)/d(Y[i]) using finite differences
  res.d_dP_dY_up.resize(Y_up.size(), 0.0);
  const double eps_Y = 1e-6;
  for (std::size_t i = 0; i < Y_up.size(); ++i) {
    std::vector<double> Y_plus = Y_up;
    Y_plus[i] += eps_Y;
    auto X_plus = combaero::mass_to_mole(combaero::normalize_fractions(Y_plus));
    auto [rho_plus, dummy1, dummy2] = density_and_jacobians(T_up, P_total_up, X_plus);
    double u_plus = m_dot / (rho_plus * area);

    auto [dP_plus, dummy3, dummy4, dummy5] =
        channel_compressible_mdot_and_jacobian(T_up, P_total_up, u_plus, X_plus, L, D, roughness, friction_model, f_multiplier, false);
    res.d_dP_dY_up[i] = (dP_plus - dP_calc) / eps_Y;
  }

  return res;
}

// -----------------------------------------------------------------------------
// Area Change Elements
// -----------------------------------------------------------------------------

AreaChangeElementResult area_change_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double T_up,
    const std::vector<double>& Y_up, double P_static_down,
    double F0, double F1, double m_scale, double D_h) {

  (void)P_total_up;
  (void)P_static_down;
  const std::vector<double> X_up =
      combaero::mass_to_mole(combaero::normalize_fractions(Y_up));
  auto [rho, drho_dT, drho_dP] =
      density_and_jacobians(T_up, P_static_up, X_up);
  auto [mu, dmu_dT, dmu_dP] =
      viscosity_and_jacobians(T_up, P_static_up, X_up);

  // Core physics evaluation
  auto core = combaero::sharp_area_change(m_dot, rho, mu, F0, F1, 0.0, m_scale, D_h);

  AreaChangeElementResult res;
  res.dP_calc = core.dP;
  res.d_dP_d_mdot = core.dS_dm;
  res.mach_clamped = core.mach_clamped;

  // Chain rule: d(dP)/dX = d(dP)/drho * drho/dX + d(dP)/dmu * dmu/dX
  res.d_dP_dP_static_up = core.dS_drho * drho_dP + core.dS_dmu * dmu_dP;
  res.d_dP_dT_up = core.dS_drho * drho_dT + core.dS_dmu * dmu_dT;

  // d(dP)/dY[i] via finite differences on rho(Y)
  res.d_dP_dY_up.resize(Y_up.size(), 0.0);
  const double eps_Y = 1e-6;
  for (std::size_t i = 0; i < Y_up.size(); ++i) {
    std::vector<double> Y_plus = Y_up;
    Y_plus[i] += eps_Y;
    auto X_plus =
        combaero::mass_to_mole(combaero::normalize_fractions(Y_plus));
    double rho_plus =
        std::get<0>(density_and_jacobians(T_up, P_static_up, X_plus));
    double mu_plus =
        std::get<0>(viscosity_and_jacobians(T_up, P_static_up, X_plus));
    auto core_plus =
        combaero::sharp_area_change(m_dot, rho_plus, mu_plus, F0, F1, 0.0, m_scale, D_h);
    res.d_dP_dY_up[i] = (core_plus.dP - core.dP) / eps_Y;
  }

  return res;
}

AreaChangeElementResult conical_area_change_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double T_up,
    const std::vector<double>& Y_up, double P_static_down,
    double F0, double F1, double length, double m_scale) {

  (void)P_total_up;
  (void)P_static_down;
  const std::vector<double> X_up =
      combaero::mass_to_mole(combaero::normalize_fractions(Y_up));
  auto [rho, drho_dT, drho_dP] =
      density_and_jacobians(T_up, P_static_up, X_up);
  auto [mu, dmu_dT, dmu_dP] =
      viscosity_and_jacobians(T_up, P_static_up, X_up);

  // Core physics evaluation
  auto core = combaero::conical_area_change(m_dot, rho, mu, F0, F1, length,
                                            0.0, m_scale);

  AreaChangeElementResult res;
  res.dP_calc = core.dP;
  res.d_dP_d_mdot = core.dS_dm;
  res.mach_clamped = core.mach_clamped;

  // Chain rule: d(dP)/dX = d(dP)/drho * drho/dX + d(dP)/dmu * dmu/dX
  // (dS_dmu is 0 for conical, but we include the terms for generality)
  res.d_dP_dP_static_up = core.dS_drho * drho_dP + core.dS_dmu * dmu_dP;
  res.d_dP_dT_up = core.dS_drho * drho_dT + core.dS_dmu * dmu_dT;

  // d(dP)/dY[i] via finite differences on rho(Y) and mu(Y)
  res.d_dP_dY_up.resize(Y_up.size(), 0.0);
  const double eps_Y = 1e-6;
  for (std::size_t i = 0; i < Y_up.size(); ++i) {
    std::vector<double> Y_plus = Y_up;
    Y_plus[i] += eps_Y;
    auto X_plus =
        combaero::mass_to_mole(combaero::normalize_fractions(Y_plus));
    double rho_plus =
        std::get<0>(density_and_jacobians(T_up, P_static_up, X_plus));
    double mu_plus =
        std::get<0>(viscosity_and_jacobians(T_up, P_static_up, X_plus));
    auto core_plus = combaero::conical_area_change(m_dot, rho_plus, mu_plus,
                                                   F0, F1, length, 0.0,
                                                   m_scale);
    res.d_dP_dY_up[i] = (core_plus.dP - core.dP) / eps_Y;
  }

  return res;
}

// -----------------------------------------------------------------------------
// 5. Tee Junction Elements (Bassett 2001)
// -----------------------------------------------------------------------------

namespace {
// Shared core computation for merging/branching tee wrappers.
// Returns residuals + analytical m_dot Jacobians given pre-computed K values.
struct TeeCore {
    double R_straight;
    double R_branch;
    double dRs_d_mcom;
    double dRs_d_mbranch;
    double dRb_d_mcom;
    double dRb_d_mbranch;
    double q_dyn;
    double q;
};

TeeCore tee_core(double m_dot_com, double m_dot_branch,
                 double dP0_straight, double dP0_branch,
                 double rho, double F_C,
                 double Ks, double dKs_dq,
                 double Kb, double dKb_dq) {
    // Regularised flow ratio: q = m_dot_branch/m_dot_com, safe at m_dot_com=0.
    // Uses q = m_b * m_c / (m_c^2 + eps^2) which gives 0 at m_c=0 and
    // correct sign for both positive and negative m_dot_com.
    constexpr double eps_m = 1e-10; // kg/s
    const double m_com_sq = m_dot_com * m_dot_com;
    const double m_com_sq_safe = m_com_sq + eps_m * eps_m;
    const double q = m_dot_branch * m_dot_com / m_com_sq_safe;

    // Regularised u_com^2 so residual is finite even at zero flow.
    constexpr double eps_v_sq = 1e-12; // (m/s)^2
    const double u_com_sq_safe =
        m_com_sq / (rho * rho * F_C * F_C) + eps_v_sq;
    const double q_dyn = 0.5 * rho * u_com_sq_safe;

    // Residuals
    const double R_straight = dP0_straight - Ks * q_dyn;
    const double R_branch = dP0_branch - Kb * q_dyn;

    // Analytical Jacobians for m_dot (psi/theta fixed; only q varies with flow)
    // dq/dm_dot_com = m_dot_branch * (eps^2 - m_com_sq) / m_com_sq_safe^2
    // dq/dm_dot_branch = m_dot_com / m_com_sq_safe
    const double dq_dmcom =
        m_dot_branch * (eps_m * eps_m - m_com_sq) / (m_com_sq_safe * m_com_sq_safe);
    const double dq_dmbranch = m_dot_com / m_com_sq_safe;

    // dq_dyn/dm_dot_com = m_dot_com / (rho * F_C^2)
    const double dqdyn_dmcom = m_dot_com / (rho * F_C * F_C);

    const double dRs_dmcom =
        -dKs_dq * dq_dmcom * q_dyn - Ks * dqdyn_dmcom;
    const double dRs_dmbranch = -dKs_dq * dq_dmbranch * q_dyn;
    const double dRb_dmcom =
        -dKb_dq * dq_dmcom * q_dyn - Kb * dqdyn_dmcom;
    const double dRb_dmbranch = -dKb_dq * dq_dmbranch * q_dyn;

    return {R_straight, R_branch, dRs_dmcom, dRs_dmbranch,
            dRb_dmcom, dRb_dmbranch, q_dyn, q};
}

// Finite-difference Jacobians for density-dependent variables (T, P_static, Y).
// Evaluates func(rho_perturbed) to get the FD derivative.
void tee_fd_density_jacobians(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static, double T,
    const std::vector<double>& X,
    double theta, double psi, double F_C, double blend_k,
    bool is_merging,
    double& dRs_dP, double& dRb_dP,
    double& dRs_dT, double& dRb_dT,
    std::vector<double>& dRs_dY, std::vector<double>& dRb_dY,
    const std::vector<double>& Y) {

    // K depends only on q/psi/theta -- compute once, reuse across all rho evals.
    constexpr double eps_m = 1e-10;
    const double q = m_dot_branch * m_dot_com
                     / (m_dot_com * m_dot_com + eps_m * eps_m);
    double Ks, dKs, Kb, dKb;
    if (is_merging) {
        Ks  = merging_tee_K_straight(q, psi, theta, blend_k);
        dKs = merging_tee_dK_straight_dq(q, psi, theta, blend_k);
        Kb  = merging_tee_K_branch(q, psi, theta, blend_k);
        dKb = merging_tee_dK_branch_dq(q, psi, theta, blend_k);
    } else {
        Ks  = branching_tee_K_straight(q, psi, theta, blend_k);
        dKs = branching_tee_dK_straight_dq(q, psi, theta, blend_k);
        Kb  = branching_tee_K_branch(q, psi, theta, blend_k);
        dKb = branching_tee_dK_branch_dq(q, psi, theta, blend_k);
    }

    auto eval = [&](double rho) {
        auto c = tee_core(m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
                          rho, F_C, Ks, dKs, Kb, dKb);
        return std::make_pair(c.R_straight, c.R_branch);
    };

    // Central FD w.r.t. P_static — scale with pressure so the step is
    // never negligibly small at high-altitude or high-pressure conditions.
    const double eps_P = std::max(1.0, std::abs(P_static) * 1e-4);
    auto [Rs_pp, Rb_pp] = eval(std::get<0>(density_and_jacobians(T, P_static + eps_P, X)));
    auto [Rs_pm, Rb_pm] = eval(std::get<0>(density_and_jacobians(T, P_static - eps_P, X)));
    dRs_dP = (Rs_pp - Rs_pm) / (2.0 * eps_P);
    dRb_dP = (Rb_pp - Rb_pm) / (2.0 * eps_P);

    // Central FD w.r.t. T
    const double eps_T = 1e-3;
    auto [Rs_tp, Rb_tp] = eval(std::get<0>(density_and_jacobians(T + eps_T, P_static, X)));
    auto [Rs_tm, Rb_tm] = eval(std::get<0>(density_and_jacobians(T - eps_T, P_static, X)));
    dRs_dT = (Rs_tp - Rs_tm) / (2.0 * eps_T);
    dRb_dT = (Rb_tp - Rb_tm) / (2.0 * eps_T);

    // Central FD w.r.t. Y[i] (composition)
    const double eps_Y = 1e-6;
    dRs_dY.assign(Y.size(), 0.0);
    dRb_dY.assign(Y.size(), 0.0);
    for (std::size_t i = 0; i < Y.size(); ++i) {
        std::vector<double> Y_p = Y, Y_m = Y;
        Y_p[i] += eps_Y;
        Y_m[i] -= eps_Y;
        auto X_p = combaero::mass_to_mole(combaero::normalize_fractions(Y_p));
        auto X_m = combaero::mass_to_mole(combaero::normalize_fractions(Y_m));
        auto [Rs_yp, Rb_yp] = eval(std::get<0>(density_and_jacobians(T, P_static, X_p)));
        auto [Rs_ym, Rb_ym] = eval(std::get<0>(density_and_jacobians(T, P_static, X_m)));
        dRs_dY[i] = (Rs_yp - Rs_ym) / (2.0 * eps_Y);
        dRb_dY[i] = (Rb_yp - Rb_ym) / (2.0 * eps_Y);
    }
}
} // anonymous namespace

TeeJunctionResult merging_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com,
    const std::vector<double>& Y_com,
    double theta, double psi, double F_C,
    double blend_k) {

    const auto X_com = combaero::mass_to_mole(
        combaero::normalize_fractions(Y_com));
    const auto [rho, drho_dT, drho_dP] =
        density_and_jacobians(T_com, P_static_com, X_com);

    // Flow ratio (regularised) and blended K values
    constexpr double eps_m = 1e-10;
    const double m_com_sq_safe =
        m_dot_com * m_dot_com + eps_m * eps_m;
    const double q = m_dot_branch * m_dot_com / m_com_sq_safe;

    const double Ks = merging_tee_K_straight(q, psi, theta, blend_k);
    const double dKs = merging_tee_dK_straight_dq(q, psi, theta, blend_k);
    const double Kb = merging_tee_K_branch(q, psi, theta, blend_k);
    const double dKb = merging_tee_dK_branch_dq(q, psi, theta, blend_k);

    auto core = tee_core(m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
                         rho, F_C, Ks, dKs, Kb, dKb);

    TeeJunctionResult res;
    res.R_straight = core.R_straight;
    res.R_branch = core.R_branch;
    res.dR_straight_d_mdot_com = core.dRs_d_mcom;
    res.dR_straight_d_mdot_branch = core.dRs_d_mbranch;
    res.dR_branch_d_mdot_com = core.dRb_d_mcom;
    res.dR_branch_d_mdot_branch = core.dRb_d_mbranch;

    tee_fd_density_jacobians(
        m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
        P_static_com, T_com, X_com, theta, psi, F_C, blend_k, true,
        res.dR_straight_dP_static_com, res.dR_branch_dP_static_com,
        res.dR_straight_dT_com, res.dR_branch_dT_com,
        res.dR_straight_dY_com, res.dR_branch_dY_com, Y_com);

    res.K_straight = Ks;
    res.K_branch = Kb;
    res.q = q;
    res.blend_w = blend_weight(q, blend_k);
    res.topology_valid = tee_topology_valid(q);
    auto inp = tee_check_inputs(q, psi, theta);
    res.status = inp.valid() ? CorrelationValidity::VALID
                             : CorrelationValidity::EXTRAPOLATED;

    return res;
}

TeeJunctionResult branching_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com,
    const std::vector<double>& Y_com,
    double theta, double psi, double F_C,
    double blend_k) {

    const auto X_com = combaero::mass_to_mole(
        combaero::normalize_fractions(Y_com));
    const auto [rho, drho_dT, drho_dP] =
        density_and_jacobians(T_com, P_static_com, X_com);

    constexpr double eps_m = 1e-10;
    const double m_com_sq_safe =
        m_dot_com * m_dot_com + eps_m * eps_m;
    const double q = m_dot_branch * m_dot_com / m_com_sq_safe;

    const double Ks = branching_tee_K_straight(q, psi, theta, blend_k);
    const double dKs = branching_tee_dK_straight_dq(q, psi, theta, blend_k);
    const double Kb = branching_tee_K_branch(q, psi, theta, blend_k);
    const double dKb = branching_tee_dK_branch_dq(q, psi, theta, blend_k);

    auto core = tee_core(m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
                         rho, F_C, Ks, dKs, Kb, dKb);

    TeeJunctionResult res;
    res.R_straight = core.R_straight;
    res.R_branch = core.R_branch;
    res.dR_straight_d_mdot_com = core.dRs_d_mcom;
    res.dR_straight_d_mdot_branch = core.dRs_d_mbranch;
    res.dR_branch_d_mdot_com = core.dRb_d_mcom;
    res.dR_branch_d_mdot_branch = core.dRb_d_mbranch;

    tee_fd_density_jacobians(
        m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
        P_static_com, T_com, X_com, theta, psi, F_C, blend_k, false,
        res.dR_straight_dP_static_com, res.dR_branch_dP_static_com,
        res.dR_straight_dT_com, res.dR_branch_dT_com,
        res.dR_straight_dY_com, res.dR_branch_dY_com, Y_com);

    res.K_straight = Ks;
    res.K_branch = Kb;
    res.q = q;
    res.blend_w = blend_weight(q, blend_k);
    res.topology_valid = tee_topology_valid(q);
    auto inp = tee_check_inputs(q, psi, theta);
    res.status = inp.valid() ? CorrelationValidity::VALID
                             : CorrelationValidity::EXTRAPOLATED;

    return res;
}

// =============================================================================
// Unified0D compressible tee junction (Section 5-7, docs/junction/)
// =============================================================================

CompressibleTeeResult compressible_branching_tee_rj(
    const BranchInput& com,
    const BranchInput& str,
    const BranchInput& bra)
{
    constexpr double eps_m = 1e-10;

    // Pseudodatum = common branch (single supplier for branching)
    const double m_com = com.m_dot;
    const double m_bra = bra.m_dot;
    const double m_str = str.m_dot;  // = m_com - m_bra, passed by caller
    const double Pc    = com.P_static;
    const double Tc    = com.T;
    const double Rc    = com.R_gas;
    const double Ac    = com.A;

    const double m_com_sq_safe = m_com * m_com + eps_m * eps_m;
    const double u_dat  = m_com * Rc * Tc / (Pc * Ac);
    const double q_ref  = 0.5 * (Pc / (Rc * Tc)) * u_dat * u_dat;
    // A_dat = Ac (derived: m_com / (rho_dat * u_dat) = m_com * Rc * Tc / (Pc * u_dat) = Ac)

    // Flow ratios (regularised, same pattern as tee_core)
    const double x_str = m_str * m_com / m_com_sq_safe * (Ac / str.A);
    const double x_bra = m_bra * m_com / m_com_sq_safe * (Ac / bra.A);

    // Effective inflow angles: phi_j = 0.75 * |theta_j - theta_dat|.
    // (The LaTeX formula 0.75*(pi-(theta_dat-theta_j)) applies the convention that
    // theta_dat is measured as the branch's outward direction; Python passes the
    // "same-axis = 0" convention so the equivalent form is the absolute difference.)
    const double phi_str = 0.75 * std::abs(str.theta - com.theta);
    const double phi_bra = 0.75 * std::abs(bra.theta - com.theta);

    // Both collector legs use one consistent blended turning-loss closure (Finding 5):
    //   R = Pt_com - Pt_leg - K_turn(theta)*q_leg - BETA_EXTRACT*x_leg*K_bc_leg*q_ref
    // where q_leg is the leg-local dynamic head (turning loss, vanishes at theta=0) and
    // K_bc_leg = 1 + x^2 - 2*x*cos(phi) is the incompressible Borda-Carnot kernel
    // (faithful to the FD prototype; phi is constant w.r.t. the unknowns, so its
    // x-derivative is the trivial 2x - 2cos(phi)). The straight leg has theta_str = 0,
    // so K_turn(0) = 0 and its closure collapses to the pure extraction term.
    const double BE = combaero::BETA_EXTRACT;

    // Straight collector (theta_str = 0 -> no turning loss).
    const double cos_phi_str = std::cos(phi_str);
    const double K_bc_str    = 1.0 + x_str * x_str - 2.0 * x_str * cos_phi_str;
    const double dKbc_str_dx = 2.0 * x_str - 2.0 * cos_phi_str;

    const double R_0 = com.Pt - str.Pt - BE * x_str * K_bc_str * q_ref;

    // Branch collector (theta_bra != 0 -> turning loss present).
    const double cos_phi_bra = std::cos(phi_bra);
    const double K_bc        = 1.0 + x_bra * x_bra - 2.0 * x_bra * cos_phi_bra;
    const double dKbc_dx     = 2.0 * x_bra - 2.0 * cos_phi_bra;
    const double K_turn      = combaero::K_turn_div(bra.theta);

    const double Pb    = bra.P_static;
    const double Tb    = bra.T;
    const double Rb    = bra.R_gas;
    const double Ab    = bra.A;
    const double q_bra = 0.5 * Rb * Tb * m_bra * m_bra / (Pb * Ab * Ab);

    const double R_1 = com.Pt - bra.Pt - K_turn * q_bra - BE * x_bra * K_bc * q_ref;

    // Pseudodatum derivatives (q_ref = 0.5 * m_com^2 * Rc * Tc / (Pc * Ac^2)).
    const double dq_dPc    = -q_ref / Pc;
    const double dq_dTc    = q_ref / Tc;
    const double dq_dm_com = m_com * Rc * Tc / (Pc * Ac * Ac);

    // Flow-ratio derivatives (regularised).
    // x_str = m_str * m_com / m_com_sq_safe * psi_str
    // In Python: m_str = m_com - m_branch, so dm_str/dm_com = 1, dm_str/dm_branch = -1.
    // Total dx_str/dm_com includes both the partial w.r.t. m_com and the chain through m_str.
    const double dx_str_dm_com_partial = m_str * (eps_m*eps_m - m_com*m_com) / (m_com_sq_safe*m_com_sq_safe) * (Ac / str.A);
    const double dx_str_dm_str = m_com / m_com_sq_safe * (Ac / str.A);  // dm_str/dm_com = 1
    const double dx_str_dm_com = dx_str_dm_com_partial + dx_str_dm_str;
    const double dx_str_dm_branch = -dx_str_dm_str;  // dm_str/dm_branch = -1
    const double dx_bra_dm_com    = m_bra * (eps_m*eps_m - m_com*m_com) / (m_com_sq_safe*m_com_sq_safe) * (Ac / bra.A);
    const double dx_bra_dm_branch = m_com / m_com_sq_safe * (Ac / bra.A);   // dm_bra/dm_branch = +1

    CompressibleTeeResult res{};
    res.R_0 = R_0;
    res.R_1 = R_1;

    // R_0 Jacobians: R_0 = Pt_com - Pt_str - BETA_EXTRACT*x_str*K_bc_str*q_ref.
    // Pure extraction term (straight leg, no turning loss); phi_str is constant w.r.t.
    // the unknowns. d(x_str*K_bc_str)/dx_str = K_bc_str + x_str*dKbc_str_dx.
    const double d_xKbc_str_dx = K_bc_str + x_str * dKbc_str_dx;
    res.dR0_dPt_com    = 1.0;
    res.dR0_dPt_str    = -1.0;
    res.dR0_dP_com     = -BE * x_str * K_bc_str * dq_dPc;
    res.dR0_dT_com     = -BE * x_str * K_bc_str * dq_dTc;
    res.dR0_dmdot_com  = -BE * (q_ref * d_xKbc_str_dx * dx_str_dm_com + x_str * K_bc_str * dq_dm_com);
    res.dR0_dmdot_branch = -BE * q_ref * d_xKbc_str_dx * dx_str_dm_branch;
    // R_0 has no bra node dependence for branching.
    res.dR0_dP_str = 0.0;
    res.dR0_dT_str = 0.0;
    res.dR0_dP_bra = 0.0;

    // R_1 Jacobians: R_1 = Pt_com - Pt_bra - K_turn*q_bra - BETA_EXTRACT*x_bra*K_bc*q_ref.
    // K_turn and phi_bra (hence K_bc's cos term) are constant w.r.t. the unknowns.
    // Extraction term E = BETA_EXTRACT*x_bra*K_bc*q_ref; d(x_bra*K_bc)/dx_bra = K_bc + x_bra*dKbc_dx.
    const double d_xKbc_dx      = K_bc + x_bra * dKbc_dx;
    const double dq_bra_dm_bra  = Rb * Tb * m_bra / (Pb * Ab * Ab);  // d(q_bra)/d(m_branch)

    res.dR1_dPt_com    = 1.0;
    res.dR1_dPt_bra    = -1.0;
    // Branch-local q_bra (q_bra ~ 1/P_bra, ~ T_bra); turning-loss term only.
    res.dR1_dP_bra     = K_turn * q_bra / Pb;
    res.dR1_dT_bra     = -K_turn * q_bra / Tb;
    // Common datum enters only through q_ref in the extraction term (K_bc is incompressible).
    res.dR1_dP_com     = -BE * x_bra * K_bc * dq_dPc;
    res.dR1_dT_com     = -BE * x_bra * K_bc * dq_dTc;
    // Flows: q_bra depends on m_branch only; x_bra and q_ref depend on m_com; x_bra on m_branch.
    res.dR1_dmdot_com  = -BE * (q_ref * d_xKbc_dx * dx_bra_dm_com + x_bra * K_bc * dq_dm_com);
    res.dR1_dmdot_branch = -K_turn * dq_bra_dm_bra
                           - BE * q_ref * d_xKbc_dx * dx_bra_dm_branch;
    // R_1 has no straight-node dependence for branching.
    res.dR1_dP_str = 0.0;
    res.dR1_dT_str = 0.0;
    res.dR1_dPt_str = 0.0;

    return res;
}

CompressibleTeeResult compressible_merging_tee_rj(
    const BranchInput& com,
    const BranchInput& str,
    const BranchInput& bra)
{
    using combaero::K_dat_j_closed;
    using combaero::cstep_deriv;

    const double eps_m = 1e-12;

    // Supplier branches: str and bra. Collector: com.
    const double m_str = str.m_dot;  // positive (supplier)
    const double m_bra = bra.m_dot;  // positive (supplier)
    const double m_dat = m_str + m_bra;
    const double m_dat_safe = (m_dat < eps_m) ? eps_m : m_dat;

    const double Ps = str.P_static, Ts = str.T, Rs = str.R_gas, gs = str.gamma_eff, As = str.A;
    const double Pb = bra.P_static, Tb = bra.T, Rb = bra.R_gas, gb = bra.gamma_eff, Ab = bra.A;

    const double cp_s = gs * Rs / (gs - 1.0);
    const double cp_b = gb * Rb / (gb - 1.0);
    const double cp_dat = cp_s;  // reference = straight supplier

    // Per-branch velocities and stagnation enthalpies
    const double u_str = m_str * Rs * Ts / (Ps * As);
    const double u_bra = m_bra * Rb * Tb / (Pb * Ab);
    const double h0_str = cp_s * Ts + 0.5 * u_str * u_str;
    const double h0_bra = cp_b * Tb + 0.5 * u_bra * u_bra;

    // Pseudodatum (mass-weighted average of suppliers)
    const double u_dat   = (m_str * u_str + m_bra * u_bra) / m_dat_safe;
    const double h0_dat  = (m_str * h0_str + m_bra * h0_bra) / m_dat_safe;
    const double T_dat   = (h0_dat - 0.5 * u_dat * u_dat) / cp_dat;
    const double p_dat   = Ps;  // reference = straight supplier static pressure
    const double g_dat   = gs;
    const double rho_dat = p_dat / (Rs * T_dat);
    const double M2_dat  = (u_dat > eps_m) ? u_dat * u_dat / (g_dat * Rs * T_dat) : 0.0;
    const double q_ref   = 0.5 * rho_dat * u_dat * u_dat;
    const double A_dat   = (u_dat > eps_m) ? m_dat_safe / (rho_dat * u_dat) : 1.0;

    // Pseudodatum angle (mass-weighted)
    const double S_S = m_str * std::sin(str.theta) + m_bra * std::sin(bra.theta);
    const double C_S = m_str * std::cos(str.theta) + m_bra * std::cos(bra.theta);
    const double SC_sq = S_S * S_S + C_S * C_S;
    const double theta_dat = std::atan2(S_S, C_S);

    // K_com: loss coefficient from pseudodatum to collector (common) branch.
    // x_com = A_dat / A_com  (m_com = m_dat so flow ratio = 1, area ratio gives x).
    const double x_com_v   = (u_dat > eps_m) ? A_dat / com.A : 0.0;
    const double phi_com_v = 0.75 * std::abs(com.theta - theta_dat);
    // dphi_com / d(theta_dat): sign depends on relative direction
    const double dphi_com_dtheta = (com.theta >= theta_dat) ? -0.75 : 0.75;
    const double K_com_v   = K_dat_j_closed(x_com_v, phi_com_v, M2_dat);

    // Stagnation pressure of pseudodatum: p0_dat = p_dat * Phi_dat^(gamma/(gamma-1))
    const double Phi_dat    = 1.0 + 0.5 * (g_dat - 1.0) * M2_dat;
    const double p0_dat_val = (u_dat > eps_m)
        ? p_dat * std::pow(Phi_dat, g_dat / (g_dat - 1.0)) : p_dat;

    // Residuals (LaTeX Sec 6, merging 2-supplier / 1-collector):
    //   R_sp      = p_str - p_bra              (equal static pressure at merge point)
    //   R_K,com   = p0_dat - p0_com - K * q    (stagnation loss from dat to collector)
    const double R_0 = Ps - Pb;
    const double R_1 = p0_dat_val - com.Pt - K_com_v * q_ref;

    // Complex-step K derivatives for K_com
    auto K_fn = [](auto x, auto phi, auto M2) { return K_dat_j_closed(x, phi, M2); };
    const double dKc_dx   = cstep_deriv(K_fn, 0, x_com_v, phi_com_v, M2_dat);
    const double dKc_dphi = cstep_deriv(K_fn, 1, x_com_v, phi_com_v, M2_dat);
    const double dKc_dM2  = cstep_deriv(K_fn, 2, x_com_v, phi_com_v, M2_dat);

    // Helper: compute d(-K_com * q_ref)/dv from pseudodatum chain.
    // m_j = m_dat, A_j = com.A, dm_j_dv = dm_d (since m_j = m_dat).
    auto dNeg_Kq = [&](double du_d, double dh0_d,
                        double dm_d, double dp_d, double dtheta_d) -> double
    {
        const double dT_d = (dh0_d - u_dat * du_d) / cp_dat;
        const double drho_d = rho_dat * (dp_d / p_dat - dT_d / T_dat);
        const double dM2_d = (u_dat > eps_m)
            ? M2_dat * (2.0 * du_d / u_dat - dT_d / T_dat) : 0.0;
        const double dq_d = (u_dat > eps_m)
            ? q_ref * (drho_d / rho_dat + 2.0 * du_d / u_dat) : 0.0;
        const double dA_d = (u_dat > eps_m)
            ? A_dat * (dm_d / m_dat_safe - drho_d / rho_dat - du_d / u_dat) : 0.0;
        // x_com = A_dat / A_com, so dx_com = dA_d / A_com
        const double dx_j  = dA_d / com.A;
        const double dphi_j = dphi_com_dtheta * dtheta_d;
        const double dK_j  = dKc_dx * dx_j + dKc_dphi * dphi_j + dKc_dM2 * dM2_d;
        return -(dK_j * q_ref + K_com_v * dq_d);
    };

    // Helper: d(p0_dat)/dv via pseudodatum chain.
    auto dp0dat_dv = [&](double du_d, double dh0_d, double dp_d) -> double {
        if (u_dat < eps_m) return dp_d;
        const double dT_d  = (dh0_d - u_dat * du_d) / cp_dat;
        const double dM2_d = M2_dat * (2.0 * du_d / u_dat - dT_d / T_dat);
        return (p0_dat_val / p_dat) * dp_d
             + p0_dat_val * g_dat / 2.0 / Phi_dat * dM2_d;
    };

    CompressibleTeeResult res{};
    res.R_0 = R_0;
    res.R_1 = R_1;

    // R_0 = Ps - Pb: trivial Jacobian (only static pressures of suppliers)
    res.dR0_dPt_str = 0.0; res.dR0_dPt_com = 0.0; res.dR0_dPt_bra = 0.0;
    res.dR0_dP_str  = 1.0; res.dR0_dP_com  = 0.0; res.dR0_dP_bra  = -1.0;
    res.dR0_dT_str  = 0.0; res.dR0_dT_com  = 0.0; res.dR0_dT_bra  = 0.0;
    res.dR0_dmdot_com = 0.0; res.dR0_dmdot_branch = 0.0;

    // R_1 = p0_dat - Pt_com - K_com*q_ref
    // Direct term: dR1/dPt_com = -1
    res.dR1_dPt_com = -1.0;
    res.dR1_dPt_str = 0.0; res.dR1_dPt_bra = 0.0;
    res.dR1_dP_com  = 0.0; res.dR1_dT_com  = 0.0;

    // dR1/dP_str: u_str depends on Ps; p_dat = Ps so dp_d = 1
    {
        const double du_s = -u_str / Ps;
        const double dh0_s = u_str * du_s;
        const double du_d = m_str / m_dat_safe * du_s;
        const double dh0_d = m_str / m_dat_safe * dh0_s;
        res.dR1_dP_str = dp0dat_dv(du_d, dh0_d, 1.0)
                       + dNeg_Kq(du_d, dh0_d, 0.0, 1.0, 0.0);
    }
    // dR1/dT_str
    {
        const double du_s = u_str / Ts;
        const double dh0_s = cp_s + u_str * u_str / Ts;
        const double du_d = m_str / m_dat_safe * du_s;
        const double dh0_d = m_str / m_dat_safe * dh0_s;
        res.dR1_dT_str = dp0dat_dv(du_d, dh0_d, 0.0)
                       + dNeg_Kq(du_d, dh0_d, 0.0, 0.0, 0.0);
    }
    // dR1/dP_bra: p_dat = Ps so dp_d = 0 for bra
    {
        const double du_b = -u_bra / Pb;
        const double dh0_b = u_bra * du_b;
        const double du_d = m_bra / m_dat_safe * du_b;
        const double dh0_d = m_bra / m_dat_safe * dh0_b;
        res.dR1_dP_bra = dp0dat_dv(du_d, dh0_d, 0.0)
                       + dNeg_Kq(du_d, dh0_d, 0.0, 0.0, 0.0);
    }
    // dR1/dT_bra
    {
        const double du_b = u_bra / Tb;
        const double dh0_b = cp_b + u_bra * u_bra / Tb;
        const double du_d = m_bra / m_dat_safe * du_b;
        const double dh0_d = m_bra / m_dat_safe * dh0_b;
        res.dR1_dT_bra = dp0dat_dv(du_d, dh0_d, 0.0)
                       + dNeg_Kq(du_d, dh0_d, 0.0, 0.0, 0.0);
    }
    // dR1/dm_dot_com: m_str = m_com - m_bra, dm_str/dm_com = 1, dm_dat/dm_com = 1
    {
        const double du_d = (2.0 * u_str - u_dat) / m_dat_safe;
        const double dh0_d = (u_str * u_str + h0_str - h0_dat) / m_dat_safe;
        const double dtheta_d = (SC_sq > 1e-30)
            ? (C_S * std::sin(str.theta) - S_S * std::cos(str.theta)) / SC_sq : 0.0;
        res.dR1_dmdot_com = dp0dat_dv(du_d, dh0_d, 0.0)
                          + dNeg_Kq(du_d, dh0_d, 1.0, 0.0, dtheta_d);
    }
    // dR1/dm_dot_branch: dm_str/dm_branch = -1, dm_bra/dm_branch = +1, dm_dat = 0
    {
        const double du_d = (-2.0 * u_str + 2.0 * u_bra) / m_dat_safe;
        const double dh0_d = (-h0_str - u_str*u_str + h0_bra + u_bra*u_bra) / m_dat_safe;
        const double dt_s = (SC_sq > 1e-30)
            ? (C_S * std::sin(str.theta) - S_S * std::cos(str.theta)) / SC_sq : 0.0;
        const double dt_b = (SC_sq > 1e-30)
            ? (C_S * std::sin(bra.theta) - S_S * std::cos(bra.theta)) / SC_sq : 0.0;
        const double dtheta_d = -dt_s + dt_b;
        res.dR1_dmdot_branch = dp0dat_dv(du_d, dh0_d, 0.0)
                             + dNeg_Kq(du_d, dh0_d, 0.0, 0.0, dtheta_d);
    }

    return res;
}

} // namespace solver
} // namespace combaero
