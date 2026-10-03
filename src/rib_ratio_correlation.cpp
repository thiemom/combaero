#include "rib_ratio_correlation.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "heat_transfer.h"

namespace combaero {
namespace cooling {

namespace {

// Smooth one-sided floor for Pr: equals Pr well above PR_EPS, tends to it
// below, differentiable throughout. Pr is a fluid property and positive in
// practice; the guard only keeps Pr^n finite for a state no fluid reaches.
constexpr double PR_EPS = 1e-4;

double smooth_pr(double Pr) {
  return 0.5 * (Pr + std::sqrt(Pr * Pr + PR_EPS * PR_EPS));
}

// (value / reference)^exponent, 1 for a zero exponent; geometry terms only,
// whose values are fixed per evaluation and need no derivative.
double geometry_term(const RibTerm &t, double value) {
  if (t.exponent == 0.0) {
    return 1.0;
  }
  const double x = value / t.reference;
  if (!(x > 0.0) || !std::isfinite(x)) {
    return 1.0;
  }
  return std::pow(x, t.exponent);
}

struct ValueSlope {
  double v = 0.0;
  double d = 0.0;  // derivative with respect to the smoothed |Re|
};

// C * geometry * (x / Re_ref)^a, and its x-derivative. x > 0 always.
ValueSlope ratio_law(double C, double geometry, const RibTerm &re_term,
                     double x) {
  const double v =
      C * geometry *
      (re_term.exponent == 0.0
           ? 1.0
           : std::pow(x / re_term.reference, re_term.exponent));
  return {v, re_term.exponent * v / x};
}

ValueSlope baseline_law(const RatioBaseline &b, double x, double pr) {
  const double v =
      b.coeff * std::pow(x, b.re_exponent) * std::pow(pr, b.pr_exponent);
  return {v, b.re_exponent * v / x};
}

bool outside(const RibRange &r, double v) {
  if (!r.bounded()) {
    return r.lo > 0.0 && v < r.lo;
  }
  return v < r.lo || v > r.hi;
}

// C1 smoothstep weight in ln x over [x_lo, x_hi]: 0 at and below x_lo, 1 at
// and above x_hi, zero slope at both ends.
ValueSlope blend_weight(double x, double x_lo, double x_hi) {
  if (x >= x_hi) {
    return {1.0, 0.0};
  }
  if (x <= x_lo) {
    return {0.0, 0.0};
  }
  const double span = std::log(x_hi / x_lo);
  const double s = std::log(x / x_lo) / span;
  return {s * s * (3.0 - 2.0 * s), 6.0 * s * (1.0 - s) / (span * x)};
}

void require_term(const std::string &set_name, const RibTerm &t,
                  const char *field) {
  if (!(t.reference > 0.0) || !std::isfinite(t.reference) ||
      !std::isfinite(t.exponent)) {
    throw std::invalid_argument("rib ratio set '" + set_name + "': " +
                                std::string(field) +
                                " needs a finite positive reference and a "
                                "finite exponent");
  }
}

void require_baseline(const std::string &set_name, const RatioBaseline &b,
                      const char *field) {
  if (!(b.coeff > 0.0) || !std::isfinite(b.coeff) ||
      !std::isfinite(b.re_exponent) || !std::isfinite(b.pr_exponent)) {
    throw std::invalid_argument(
        "rib ratio set '" + set_name + "': " + std::string(field) +
        " must have a finite positive coeff and finite exponents. A ratio "
        "without the baseline it was fitted to is not a number.");
  }
}

void check_accuracy(const std::string &set_name, const char *field,
                    const StatedAccuracy &a) {
  const bool has_value = std::isfinite(a.value);
  const bool claims_value = a.provenance != AccuracyProvenance::Unstated;
  if (has_value != claims_value || (claims_value && !(a.value >= 0.0))) {
    throw std::invalid_argument(
        "rib ratio set '" + set_name + "': " + std::string(field) +
        " must carry a finite, non-negative value if and only if its "
        "provenance is not Unstated");
  }
}

}  // namespace

RibRatioResult evaluate_rib_ratio(const RibRatioSet &set,
                                  const RibGeometry &geom, double Re,
                                  double Pr, const RibRatioOptions &options) {
  RibRatioResult out;

  // Smoothed |Re| and its derivative: even values, odd derivatives.
  const double x = std::sqrt(Re * Re + RATIO_RE_EPS * RATIO_RE_EPS);
  const double dx_dRe = Re / x;
  const double pr = smooth_pr(Pr);

  const double gN = geometry_term(set.Nu_eD, geom.e_D) *
                    geometry_term(set.Nu_pe, geom.p_e) *
                    geometry_term(set.Nu_WH, geom.W_H) *
                    geometry_term(set.Nu_alpha,
                                  geom.alpha_deg / 90.0 * set.Nu_alpha.reference);
  const double gF = geometry_term(set.f_eD, geom.e_D) *
                    geometry_term(set.f_pe, geom.p_e) *
                    geometry_term(set.f_WH, geom.W_H) *
                    geometry_term(set.f_alpha,
                                  geom.alpha_deg / 90.0 * set.f_alpha.reference);

  const double x_hi = set.Re_floor;
  const double x_lo = set.Re_floor / RATIO_BLEND_FACTOR;
  const ValueSlope w = blend_weight(x, x_lo, x_hi);

  // ---- Nu: r * Nu0_source above the floor, k * Nu0_ext below ----
  const ValueSlope rN = ratio_law(set.C_Nu, gN, set.Nu_Re, x);
  const ValueSlope Ns = baseline_law(set.Nu0_source, x, pr);
  const ValueSlope A = {rN.v * Ns.v, rN.d * Ns.v + rN.v * Ns.d};
  out.ratio_Nu = rN.v;

  double nu = A.v;
  double dnu = A.d;
  if (w.v < 1.0) {
    auto ext = [&](double at) -> ValueSlope {
      switch (options.below_floor) {
        case RatioBelowFloor::Gnielinski: {
          const auto g = nusselt_gnielinski_smooth_with_derivative(at, pr);
          return {g.Nu, g.dNu_dRe};
        }
        case RatioBelowFloor::User:
          return baseline_law(options.user_Nu0, at, pr);
        case RatioBelowFloor::SourceBaseline:
        default:
          return baseline_law(set.Nu0_source, at, pr);
      }
    };
    ValueSlope E = ext(x);
    ValueSlope E_floor = ext(x_hi);
    if (!(E_floor.v > 0.0) || !std::isfinite(E_floor.v)) {
      // Unvalidated options (validate_rib_ratio_options rejects them): keep
      // the evaluator total by falling back to the source baseline.
      E = baseline_law(set.Nu0_source, x, pr);
      E_floor = baseline_law(set.Nu0_source, x_hi, pr);
    }
    const ValueSlope rN_floor = ratio_law(set.C_Nu, gN, set.Nu_Re, x_hi);
    const ValueSlope Ns_floor = baseline_law(set.Nu0_source, x_hi, pr);
    // k matches the source form's VALUE at the floor; the blend supplies C1.
    const double k = rN_floor.v * Ns_floor.v / E_floor.v;
    nu = w.v * A.v + (1.0 - w.v) * k * E.v;
    dnu = w.d * (A.v - k * E.v) + w.v * A.d + (1.0 - w.v) * k * E.d;
  }
  out.Nu = nu;
  out.dNu_dRe = dnu * dx_dRe;

  // ---- f: ratio held at its floor value below the range ----
  const ValueSlope rF = ratio_law(set.C_f, gF, set.f_Re, x);
  const ValueSlope rF_floor = ratio_law(set.C_f, gF, set.f_Re, x_hi);
  const double rf = w.v * rF.v + (1.0 - w.v) * rF_floor.v;
  const double drf = w.d * (rF.v - rF_floor.v) + w.v * rF.d;
  const ValueSlope F0 = baseline_law(set.f0_source, x, pr);
  out.ratio_f = rf;
  out.f = rf * F0.v;
  out.df_dRe = (drf * F0.v + rf * F0.d) * dx_dRe;

  out.below_floor = x < set.Re_floor;
  out.extrapolated = out.below_floor || outside(set.valid_Re, std::abs(Re)) ||
                     outside(set.valid_eD, geom.e_D) ||
                     outside(set.valid_pe, geom.p_e) ||
                     outside(set.valid_WH, geom.W_H) ||
                     outside(set.valid_alpha, geom.alpha_deg);
  return out;
}

void validate_rib_ratio_set(const RibRatioSet &set) {
  if (set.name.empty()) {
    throw std::invalid_argument("rib ratio set: name must not be empty");
  }
  if (set.source.empty()) {
    throw std::invalid_argument("rib ratio set '" + set.name +
                                "': source must not be empty");
  }
  if (!(set.C_Nu > 0.0) || !std::isfinite(set.C_Nu) || !(set.C_f > 0.0) ||
      !std::isfinite(set.C_f)) {
    throw std::invalid_argument(
        "rib ratio set '" + set.name +
        "': C_Nu and C_f must be finite and positive -- a non-positive "
        "multiplier makes Nu or f negative");
  }
  if (!(set.Re_floor > 0.0) || !std::isfinite(set.Re_floor)) {
    throw std::invalid_argument("rib ratio set '" + set.name +
                                "': Re_floor must be finite and positive");
  }
  require_baseline(set.name, set.Nu0_source, "Nu0_source");
  require_baseline(set.name, set.f0_source, "f0_source");
  require_term(set.name, set.Nu_Re, "Nu_Re");
  require_term(set.name, set.Nu_eD, "Nu_eD");
  require_term(set.name, set.Nu_pe, "Nu_pe");
  require_term(set.name, set.Nu_WH, "Nu_WH");
  require_term(set.name, set.Nu_alpha, "Nu_alpha");
  require_term(set.name, set.f_Re, "f_Re");
  require_term(set.name, set.f_eD, "f_eD");
  require_term(set.name, set.f_pe, "f_pe");
  require_term(set.name, set.f_WH, "f_WH");
  require_term(set.name, set.f_alpha, "f_alpha");
  check_accuracy(set.name, "accuracy_Nu", set.accuracy_Nu);
  check_accuracy(set.name, "accuracy_f", set.accuracy_f);
}

void validate_rib_ratio_options(const RibRatioOptions &options) {
  if (options.below_floor == RatioBelowFloor::User) {
    require_baseline("options", options.user_Nu0, "user_Nu0");
  }
}

}  // namespace cooling
}  // namespace combaero
