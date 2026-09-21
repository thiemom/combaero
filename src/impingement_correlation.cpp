#include "impingement_correlation.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "math_constants.h"

namespace combaero {
namespace cooling {

namespace {

// Smooth floor so a power-law argument stays positive and finite even when a
// solver momentarily probes a zero, negative or reversed-flow input --
// mirrors EPLUS_FLOOR/ep_safe in rib_correlation.cpp. Both floors below are
// tiny relative to the correlations' real operating scales (Re_j ~ 1e3-1e5,
// z/d*Gc/Gj ~ 0-2.4), so they only matter in states that are already far
// outside any physical operating point.
constexpr double RE_FLOOR = 1.0;
constexpr double CROSSFLOW_ARG_FLOOR = 1e-3;

double smooth_magnitude(double x, double floor_value) {
  return std::sqrt(x * x + floor_value * floor_value);
}

// value^exponent, guarded the same way rib_correlation.cpp's power_term
// guards its ratio: a non-positive or non-finite base contributes a neutral
// 1.0 rather than NaN.
double guarded_power(double value, double exponent) {
  if (exponent == 0.0) {
    return 1.0;
  }
  if (!(value > 0.0) || !std::isfinite(value)) {
    return 1.0;
  }
  return std::pow(value, exponent);
}

double geometric_power_fit(const JetArrayGeometricFit &fit, double xn_d,
                           double yn_d, double z_d) {
  return fit.C * guarded_power(xn_d, fit.nx) * guarded_power(yn_d, fit.ny) *
         guarded_power(z_d, fit.nz);
}

bool outside(const ImpingementRange &r, double v) {
  if (!r.bounded()) {
    return r.lo > 0.0 && v < r.lo;
  }
  return v < r.lo || v > r.hi;
}

// beta = C_D * sqrt(2) * (pi/4) / [(yn/d)(z/d)], Florschuetz's own flow
// distribution parameter, introduced following their Eq. (7).
double crossflow_beta(double yn_d, double z_d, double C_D) {
  const double denom = yn_d * z_d;
  if (!(denom > 0.0) || !std::isfinite(denom) || !(C_D > 0.0)) {
    return 0.0;
  }
  return C_D * M_SQRT2 * (M_PI / 4.0) / denom;
}

void validate_fit(const JetArrayGeometricFit &fit, const char *name) {
  if (!std::isfinite(fit.C) || !std::isfinite(fit.nx) ||
      !std::isfinite(fit.ny) || !std::isfinite(fit.nz)) {
    throw std::invalid_argument(std::string("validate_jet_array_set: ") +
                                name + " has a non-finite field");
  }
}

}  // namespace

SingleJetImpingementSet goldstein_1986_single_jet() {
  SingleJetImpingementSet s;
  s.name = "goldstein_1986_single_jet";
  s.source = "Goldstein, Behbahani and Heppelmann (1986), IJHMT 29(8), 1227, "
             "via Han, Dutta & Ekkad (2012) 2nd ed. Eq. 4.1";
  s.A = 24.0;
  s.B = 533.0;
  s.C = 44.0;
  s.Re_exponent = 0.76;
  s.n_const_heat_flux = 1.285;
  s.n_const_wall_temp = 1.394;
  return s;
}

double single_jet_impingement_nu(const SingleJetImpingementSet &set,
                                 ImpingementThermalBC bc, double Re,
                                 double L_D, double R_D) {
  const double n = (bc == ImpingementThermalBC::ConstantHeatFlux)
                       ? set.n_const_heat_flux
                       : set.n_const_wall_temp;
  const double Re_safe = smooth_magnitude(Re, RE_FLOOR);
  const double numerator = set.A - std::abs(L_D - 7.75);
  const double denominator = set.B + set.C * guarded_power(R_D, n);
  return std::pow(Re_safe, set.Re_exponent) * numerator / denominator;
}

JetArrayCorrelationSet florschuetz_1981_inline() {
  JetArrayCorrelationSet s;
  s.name = "florschuetz_1981_inline";
  s.source = "Florschuetz, Truman and Metzger (1981), ASME J. Heat Transfer "
             "103, 337-342, Eq. (10a)/(10b)/Table 2, inline pattern, via "
             "Han, Dutta & Ekkad (2012) 2nd ed. Eq. 4.9/Table 4.1";
  s.validity_source =
      "Florschuetz, Truman and Metzger (1981), p. 337 -- NOT Han's "
      "reprinted Table 4.1 validity box, which is a misprint (see "
      "han_impingement.md item 19/I1)";
  s.pattern = JetHolePattern::Inline;

  s.A_fit = {1.18, -0.944, -0.642, 0.169};
  s.m_fit = {0.612, 0.059, 0.032, -0.022};
  s.B_fit = {0.437, -0.095, -0.219, 0.275};
  s.n_fit = {0.092, -0.005, 0.599, 1.04};

  s.valid_Re_j = {2500.0, 70000.0};
  s.valid_Gc_Gj = {0.0, 0.8};
  s.valid_xn_d = {5.0, 15.0};
  s.valid_yn_d = {4.0, 8.0};
  s.valid_z_d = {1.0, 3.0};
  s.valid_aspect_ratio = {0.625, 3.75};
  s.standard_error = 0.056;
  return s;
}

JetArrayCorrelationSet florschuetz_1981_staggered() {
  JetArrayCorrelationSet s;
  s.name = "florschuetz_1981_staggered";
  s.source = "Florschuetz, Truman and Metzger (1981), ASME J. Heat Transfer "
             "103, 337-342, Eq. (10a)/(10b)/Table 2, staggered pattern, via "
             "Han, Dutta & Ekkad (2012) 2nd ed. Eq. 4.9/Table 4.1";
  s.validity_source =
      "Florschuetz, Truman and Metzger (1981), p. 337 -- NOT Han's "
      "reprinted Table 4.1 validity box, which is a misprint (see "
      "han_impingement.md item 19/I1)";
  s.pattern = JetHolePattern::Staggered;

  s.A_fit = {1.87, -0.771, -0.999, -0.257};
  s.m_fit = {0.571, 0.028, 0.092, 0.039};
  s.B_fit = {1.03, -0.243, -0.307, 0.059};
  s.n_fit = {0.442, 0.098, -0.003, 0.304};

  s.valid_Re_j = {2500.0, 70000.0};
  s.valid_Gc_Gj = {0.0, 0.8};
  // Staggered's own tighter xn/d bound (5-10, vs inline's 5-15) -- item 20.
  s.valid_xn_d = {5.0, 10.0};
  s.valid_yn_d = {4.0, 8.0};
  s.valid_z_d = {1.0, 3.0};
  s.valid_aspect_ratio = {0.625, 3.75};
  s.standard_error = 0.061;
  return s;
}

JetArrayImpingementResult jet_array_impingement_nu(
    const JetArrayCorrelationSet &set, double Re_j, double Gc_Gj, double Pr,
    double xn_d, double yn_d, double z_d) {
  JetArrayImpingementResult out;

  const double A = geometric_power_fit(set.A_fit, xn_d, yn_d, z_d);
  const double m = geometric_power_fit(set.m_fit, xn_d, yn_d, z_d);
  const double B = geometric_power_fit(set.B_fit, xn_d, yn_d, z_d);
  const double n = geometric_power_fit(set.n_fit, xn_d, yn_d, z_d);

  const double Re_j_safe = smooth_magnitude(Re_j, RE_FLOOR);
  const double Re_term = std::pow(Re_j_safe, m);

  // z/d*Gc/Gj is a magnitude by construction (z/d > 0, Gc/Gj a ratio of two
  // mass fluxes) but floored the same way in case a solver probes a
  // momentarily negative Gc/Gj under reversed or transient flow -- see
  // RE_FLOOR's comment above.
  const double crossflow_arg = z_d * Gc_Gj;
  const double crossflow_safe = smooth_magnitude(crossflow_arg, CROSSFLOW_ARG_FLOOR);
  const double bracket = 1.0 - B * std::pow(crossflow_safe, n);

  const double Pr_term = (Pr > 0.0 && std::isfinite(Pr)) ? std::cbrt(Pr) : 1.0;

  out.Nu = A * Re_term * bracket * Pr_term;
  out.extrapolated = outside(set.valid_Re_j, std::abs(Re_j)) ||
                     outside(set.valid_Gc_Gj, Gc_Gj) ||
                     outside(set.valid_xn_d, xn_d) ||
                     outside(set.valid_yn_d, yn_d) ||
                     outside(set.valid_z_d, z_d) ||
                     outside(set.valid_aspect_ratio, xn_d / yn_d);
  return out;
}

double crossflow_to_jet_ratio_at_x(double yn_d, double z_d, double C_D,
                                   double x_over_xn) {
  const double beta = crossflow_beta(yn_d, z_d, C_D);
  if (!(C_D > 0.0) || !std::isfinite(C_D)) {
    return 0.0;
  }
  const double denom = M_SQRT2 * C_D * std::cosh(beta * x_over_xn);
  if (!(denom > 0.0) || !std::isfinite(denom)) {
    return 0.0;
  }
  return std::sinh(beta * (x_over_xn - 0.5)) / denom;
}

double crossflow_to_jet_ratio_at_row(double yn_d, double z_d, double C_D,
                                     int row) {
  return crossflow_to_jet_ratio_at_x(yn_d, z_d, C_D,
                                     static_cast<double>(row) - 0.5);
}

void validate_single_jet_set(const SingleJetImpingementSet &set) {
  if (!(set.A > 0.0) || !std::isfinite(set.A)) {
    throw std::invalid_argument("validate_single_jet_set: A must be positive");
  }
  if (!(set.B > 0.0) || !std::isfinite(set.B)) {
    throw std::invalid_argument("validate_single_jet_set: B must be positive");
  }
  if (!(set.C >= 0.0) || !std::isfinite(set.C)) {
    throw std::invalid_argument(
        "validate_single_jet_set: C must be non-negative");
  }
  if (!std::isfinite(set.Re_exponent)) {
    throw std::invalid_argument(
        "validate_single_jet_set: Re_exponent must be finite");
  }
  if (!std::isfinite(set.n_const_heat_flux)) {
    throw std::invalid_argument(
        "validate_single_jet_set: n_const_heat_flux must be finite");
  }
  if (!std::isfinite(set.n_const_wall_temp)) {
    throw std::invalid_argument(
        "validate_single_jet_set: n_const_wall_temp must be finite");
  }
}

void validate_jet_array_set(const JetArrayCorrelationSet &set) {
  validate_fit(set.A_fit, "A_fit");
  validate_fit(set.m_fit, "m_fit");
  validate_fit(set.B_fit, "B_fit");
  validate_fit(set.n_fit, "n_fit");
}

}  // namespace cooling
}  // namespace combaero
