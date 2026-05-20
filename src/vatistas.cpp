// Vatistas n-vortex model.
// Reference: Vatistas, Kozel, Mih (1991), Exp. Fluids 11, 73-76.
//            DOI: 10.1007/BF00198434

#include "vatistas.h"
#include "math_constants.h"
#include "quadrature.h"
#include <cmath>
#include <stdexcept>
#include <tuple>

namespace combaero::solver {

// ---------------------------------------------------------------------------
// Normalised shape functions
// ---------------------------------------------------------------------------

double vatistas_v0_bar(double r_bar, double n) {
  return r_bar / std::pow(1.0 + std::pow(r_bar, 2.0 * n), 1.0 / n);
}

double vatistas_dv0_bar_drbar(double r_bar, double n) {
  double r_2n = std::pow(r_bar, 2.0 * n);
  double denom = std::pow(1.0 + r_2n, 1.0 + 1.0 / n);
  return (1.0 - r_2n) / denom;
}

double vatistas_vr_bar(double r_bar, double n) {
  // Eq. (3): Vr_bar = -2*(1+n)*r_bar^{2n-1} / (1 + r_bar^{2n})
  // Negative by definition (inward flow toward the vortex axis).
  double r_2n = std::pow(r_bar, 2.0 * n);
  double r_2n_1 = std::pow(r_bar, 2.0 * n - 1.0);
  return -2.0 * (1.0 + n) * r_2n_1 / (1.0 + r_2n);
}

double vatistas_dvr_bar_drbar(double r_bar, double n) {
  // d/dr_bar[ -2(1+n)*r^{2n-1} / (1+r^{2n}) ]
  // = -2(1+n) * [ (2n-1)*r^{2n-2}*(1+r^{2n}) - r^{2n-1}*2n*r^{2n-1} ]
  //             / (1+r^{2n})^2
  // = -2(1+n) * [ (2n-1)*r^{2n-2} - r^{4n-2} ] / (1+r^{2n})^2
  double r_2n = std::pow(r_bar, 2.0 * n);
  double denom = (1.0 + r_2n) * (1.0 + r_2n);
  double num = (2.0 * n - 1.0) * std::pow(r_bar, 2.0 * n - 2.0) -
               std::pow(r_bar, 4.0 * n - 2.0);
  return -2.0 * (1.0 + n) * num / denom;
}

// The integrand of the pressure integral: V0_bar(r')^2 / r'
// = r' / (1 + r'^{2n})^{2/n}
static double pressure_integrand(double r_bar, double n) {
  return r_bar / std::pow(1.0 + std::pow(r_bar, 2.0 * n), 2.0 / n);
}

double vatistas_pressure_integral(double r_bar, double n) {
  // I(r_bar; n) = integral_0^{r_bar} V0_bar(r')^2 / r' dr'
  // Dimensional delta_P = (Gamma/(2*pi*r_c))^2 * rho * I(r/r_c; n)
  if (r_bar <= 0.0) return 0.0;

  if (n == 1.0) {
    return r_bar * r_bar / (2.0 * (1.0 + r_bar * r_bar));
  }
  if (n == 2.0) {
    return std::atan(r_bar * r_bar) / 2.0;
  }

  // General case: numerical integration, split domain at r_bar = 1 so that
  // the velocity peak is resolved by both sub-intervals.
  auto integrand = [n](double rb) { return pressure_integrand(rb, n); };
  if (r_bar <= 1.0) {
    return composite_simpson(integrand, 0.0, r_bar);
  }
  return composite_simpson(integrand, 0.0, 1.0) +
         composite_simpson(integrand, 1.0, r_bar);
}

double vatistas_d_pressure_integral_drbar(double r_bar, double n) {
  // d(I)/d(r_bar) = V0_bar(r_bar)^2 / r_bar  [FTC, analytical for all n]
  if (r_bar <= 0.0) return 0.0;
  return pressure_integrand(r_bar, n);
}

// ---------------------------------------------------------------------------
// Dimensional functions
// ---------------------------------------------------------------------------

static void validate_vatistas(double Gamma, double r_c, double n) {
  if (n < vatistas::n_min) {
    throw std::invalid_argument(
        "vatistas: n must be >= 1 (singularities in V_r and V_z for n < 1)");
  }
  if (Gamma <= 0.0) {
    throw std::invalid_argument("vatistas: Gamma must be > 0");
  }
  if (r_c <= 0.0) {
    throw std::invalid_argument("vatistas: r_c must be > 0");
  }
}

double vatistas_v_theta(double r, double Gamma, double r_c, double n) {
  validate_vatistas(Gamma, r_c, n);
  double r_bar = r / r_c;
  double scale = Gamma / (2.0 * M_PI * r_c);
  return scale * vatistas_v0_bar(r_bar, n);
}

std::tuple<double, double, double, double>
vatistas_v_theta_and_jacobians(double r, double Gamma, double r_c, double n) {
  validate_vatistas(Gamma, r_c, n);
  double r_bar = r / r_c;
  double scale = Gamma / (2.0 * M_PI * r_c);
  double V0 = vatistas_v0_bar(r_bar, n);
  double dV0_drbar = vatistas_dv0_bar_drbar(r_bar, n);

  double V_theta = scale * V0;

  // dV_theta/dr: chain rule through r_bar = r/r_c
  //   = scale * dV0/dr_bar * (1/r_c)
  double dV_dr = scale / r_c * dV0_drbar;

  // dV_theta/dGamma: V_theta is linear in Gamma
  //   = V_theta / Gamma
  double dV_dGamma = V_theta / Gamma;

  // dV_theta/dr_c:
  //   d(scale)/dr_c = -scale/r_c
  //   d(V0_bar)/dr_c = dV0/dr_bar * d(r_bar)/dr_c = dV0/dr_bar * (-r/r_c^2)
  //                  = dV0/dr_bar * (-r_bar/r_c)
  //   dV/dr_c = d(scale)/dr_c * V0 + scale * d(V0)/dr_c
  //           = (-scale/r_c)*V0 + scale*dV0/dr_bar*(-r_bar/r_c)
  //           = -(scale/r_c)*(V0 + r_bar*dV0_drbar)
  double dV_drc = -(scale / r_c) * (V0 + r_bar * dV0_drbar);

  return {V_theta, dV_dr, dV_dGamma, dV_drc};
}

double vatistas_delta_p(double r, double rho, double Gamma, double r_c,
                        double n) {
  validate_vatistas(Gamma, r_c, n);
  double r_bar = r / r_c;
  double scale = Gamma / (2.0 * M_PI * r_c);
  double V_scale = scale * scale * rho;
  return V_scale * vatistas_pressure_integral(r_bar, n);
}

std::tuple<double, double, double, double>
vatistas_delta_p_and_jacobians(double r, double rho, double Gamma, double r_c,
                               double n) {
  validate_vatistas(Gamma, r_c, n);
  double r_bar = r / r_c;
  double scale = Gamma / (2.0 * M_PI * r_c);
  double V_scale = scale * scale * rho;
  double I = vatistas_pressure_integral(r_bar, n);
  double dI_drbar = vatistas_d_pressure_integral_drbar(r_bar, n);

  double dP = V_scale * I;

  // d(delta_P)/dr = V_scale * dI/dr_bar * d(r_bar)/dr = V_scale * dI/dr_bar / r_c
  // = rho * V_theta^2 / r  (centrifugal balance, always analytical)
  double d_dP_dr = V_scale * dI_drbar / r_c;

  // d(delta_P)/dGamma: delta_P is quadratic in Gamma => derivative = 2*delta_P/Gamma
  double d_dP_dGamma = 2.0 * dP / Gamma;

  // d(delta_P)/dr_c:
  //   V_scale = (Gamma/(2*pi))^2 * rho / r_c^2
  //   d(V_scale)/dr_c = -2*V_scale/r_c
  //   r_bar = r/r_c  => d(I)/dr_c = dI/dr_bar * d(r_bar)/dr_c = dI_drbar * (-r/r_c^2)
  //                   = -dI_drbar * r_bar / r_c
  //   d(delta_P)/dr_c = d(V_scale)/dr_c * I + V_scale * d(I)/dr_c
  //                   = -2*V_scale/r_c * I + V_scale * (-dI_drbar * r_bar / r_c)
  //                   = -(V_scale/r_c) * (2*I + dI_drbar * r_bar)
  double d_dP_drc = -(V_scale / r_c) * (2.0 * I + dI_drbar * r_bar);

  return {dP, d_dP_dr, d_dP_dGamma, d_dP_drc};
}

} // namespace combaero::solver
