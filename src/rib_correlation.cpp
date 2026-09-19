#include "rib_correlation.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace combaero {
namespace cooling {

RibCorrelationSet han_1988_orthogonal() {
  RibCorrelationSet s;
  s.name = "han_1988_orthogonal";
  s.source = "Han, J.C. (1988). ASME J. Heat Transfer 110, 321, via Han, "
             "Dutta & Ekkad (2012) 2nd ed. Eq. 4.15/4.16, Fig. 4.46";
  s.validity_source = s.source;
  s.provenance = RibProvenance::Extracted;
  // 90 deg orthogonal ribs: reversing the flow leaves the geometry unchanged.
  s.symmetric = true;

  // R = 3.2 (p/e / 10)^0.35, independent of e+.
  s.C_R = 3.2;
  s.R_pe = {0.35, 10.0};

  // G = 3.7 (e+)^0.28 at Pr ~ 0.7. No geometry terms: Han reports G as a
  // function of e+ alone for this configuration.
  s.C_G = 3.7;
  s.G_eplus_exponent = 0.28;

  s.valid_Re = {10000.0, 60000.0};
  s.valid_eD = {0.047, 0.078};
  s.valid_pe = {10.0, 20.0};
  s.valid_WH = {1.0, 4.0};
  s.valid_alpha = {90.0, 90.0};
  // The correlation is stated for e+ >= 50; no upper bound is given.
  s.valid_eplus = {50.0, 0.0};
  s.valid_Pr = 0.7;

  s.accuracy_R = 0.06;  // 95% of data within 6%
  s.accuracy_G = 0.08;  // 95% of data within 8%
  return s;
}

namespace {

void require_positive_reference(const RibTerm &t, const std::string &field) {
  if (!(t.reference > 0.0) || !std::isfinite(t.reference)) {
    throw std::invalid_argument(
        "rib correlation set: " + field +
        ".reference must be finite and positive, got " +
        std::to_string(t.reference) +
        ". A missing reference is indistinguishable from a deliberate 1.0, so "
        "it is rejected rather than defaulted.");
  }
  if (!std::isfinite(t.exponent)) {
    throw std::invalid_argument("rib correlation set: " + field +
                                ".exponent must be finite");
  }
}

}  // namespace

namespace {

// Smooth floor on e+ so the heat-transfer power law is total and its
// derivative continuous through zero. A hard abs() or max() would put a kink
// exactly where Newton iterates; sqrt(x^2 + eps^2) does not. At e+ = 1000 the
// distortion is 1.4e-07.
//
// The floor is only significant below e+ ~ 5, which for a typical rib geometry
// is Re below roughly 700 -- an order of magnitude under the bottom of Han's
// validity (Re = 10,000, e+ >= 50) and a regime where the velocity, the
// pressure drop and the Nusselt number are all negligible anyway. The guard
// exists to keep the solver's arithmetic well behaved while it passes through
// such states, not to predict them.
constexpr double EPLUS_FLOOR = 1.0;

// The bracket (2/f)^(1/2) must stay positive for f to be meaningful. It can
// reach zero for coefficient and geometry combinations a user can supply --
// with C_R = 1.0 it happens at e/D = 0.274, inside plausible input -- so it is
// floored rather than allowed to produce an infinite or negative f.
constexpr double BRACKET_FLOOR = 1e-3;

// The Stanton denominator 1 + (G - R)(f/2)^(1/2) thins as e/D and C_R grow and
// would give a NEGATIVE Stanton number if it crossed zero, which is worse than
// a crash because it looks like a number.
constexpr double ST_DENOM_FLOOR = 1e-3;

double power_term(const RibTerm &t, double value) {
  if (t.exponent == 0.0) {
    return 1.0;
  }
  const double x = value / t.reference;
  if (!(x > 0.0) || !std::isfinite(x)) {
    return 1.0;
  }
  return std::pow(x, t.exponent);
}

double softmin_floor(double x, double floor_value) {
  // Smooth one-sided floor: equals x well above floor_value, approaches
  // floor_value below it, differentiable throughout.
  return 0.5 * (x + std::sqrt(x * x + floor_value * floor_value));
}

bool outside(const RibRange &r, double v) {
  if (!r.bounded()) {
    // A half-open range: lo set, hi unset.
    return r.lo > 0.0 && v < r.lo;
  }
  return v < r.lo || v > r.hi;
}

}  // namespace

RibResult evaluate_rib(const RibCorrelationSet &set, const RibGeometry &geom,
                       double Re) {
  RibResult out;

  const double e_D = geom.e_D;
  const double W_H = geom.W_H;

  out.R = set.C_R * power_term(set.R_eD, e_D) * power_term(set.R_pe, geom.p_e) *
          power_term(set.R_WH, W_H) *
          power_term(set.R_alpha, geom.alpha_deg / 90.0 * set.R_alpha.reference);

  // Invert the wall law for f. The geometry group is
  // (2 e/D) * (2 W/(W + H)); with W_H = W/H it is (2 e/D) * (2 W_H/(W_H + 1)).
  const double geom_group =
      (2.0 * e_D) * (2.0 * W_H / (W_H + 1.0));
  double bracket = 0.0;
  if (geom_group > 0.0 && std::isfinite(geom_group)) {
    bracket = out.R - 2.5 * std::log(geom_group) - 2.5;
  }
  bracket = softmin_floor(bracket, BRACKET_FLOOR);
  out.f = 2.0 / (bracket * bracket);

  // e+ keeps the sign of the flow so a caller can tell direction; the power
  // law below uses the smoothed magnitude.
  out.e_plus = e_D * Re * std::sqrt(out.f / 2.0);
  const double ep_safe =
      std::sqrt(out.e_plus * out.e_plus + EPLUS_FLOOR * EPLUS_FLOOR);

  out.G = set.C_G * power_term(set.G_eD, e_D) * power_term(set.G_pe, geom.p_e) *
          power_term(set.G_WH, W_H) *
          power_term(set.G_alpha, geom.alpha_deg / 90.0 * set.G_alpha.reference) *
          std::pow(ep_safe, set.G_eplus_exponent);

  const double root_f2 = std::sqrt(out.f / 2.0);
  const double denom =
      softmin_floor(1.0 + (out.G - out.R) * root_f2, ST_DENOM_FLOOR);
  out.St_r = out.f / (2.0 * denom);

  // dSt/dRe by the chain rule. f is independent of Re for this family, so only
  // G moves: dG/dRe = G * n * (d ep_safe/dRe) / ep_safe, and
  // d ep_safe/dRe = e_plus * (e_D * sqrt(f/2)) / ep_safe.
  const double dep_dRe = e_D * root_f2;
  const double dep_safe_dRe = out.e_plus * dep_dRe / ep_safe;
  const double dG_dRe = set.G_eplus_exponent * out.G / ep_safe * dep_safe_dRe;
  out.dSt_dRe = -out.f * root_f2 * dG_dRe / (2.0 * denom * denom);

  out.extrapolated =
      outside(set.valid_Re, std::abs(Re)) || outside(set.valid_eD, e_D) ||
      outside(set.valid_pe, geom.p_e) || outside(set.valid_WH, W_H) ||
      outside(set.valid_alpha, geom.alpha_deg) ||
      outside(set.valid_eplus, std::abs(out.e_plus));

  return out;
}

void validate_rib_set(const RibCorrelationSet &set) {
  if (set.name.empty()) {
    throw std::invalid_argument("rib correlation set: name must not be empty");
  }
  if (set.source.empty()) {
    throw std::invalid_argument(
        "rib correlation set '" + set.name +
        "': source must not be empty. A coefficient whose origin is not "
        "recorded is what this rebuild exists to remove.");
  }
  if (!(set.C_R > 0.0) || !std::isfinite(set.C_R)) {
    throw std::invalid_argument("rib correlation set '" + set.name +
                                "': C_R must be finite and positive, got " +
                                std::to_string(set.C_R));
  }
  if (!(set.C_G > 0.0) || !std::isfinite(set.C_G)) {
    throw std::invalid_argument("rib correlation set '" + set.name +
                                "': C_G must be finite and positive, got " +
                                std::to_string(set.C_G));
  }
  if (!std::isfinite(set.G_eplus_exponent)) {
    throw std::invalid_argument("rib correlation set '" + set.name +
                                "': G_eplus_exponent must be finite");
  }
  require_positive_reference(set.R_eD, "R_eD");
  require_positive_reference(set.R_pe, "R_pe");
  require_positive_reference(set.R_WH, "R_WH");
  require_positive_reference(set.R_alpha, "R_alpha");
  require_positive_reference(set.G_eD, "G_eD");
  require_positive_reference(set.G_pe, "G_pe");
  require_positive_reference(set.G_WH, "G_WH");
  require_positive_reference(set.G_alpha, "G_alpha");
}

}  // namespace cooling
}  // namespace combaero
