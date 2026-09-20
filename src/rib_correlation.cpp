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

RibCorrelationSet rallabandi_2009_high_re() {
  RibCorrelationSet s;
  s.name = "rallabandi_2009_high_re";
  s.source = "Rallabandi, A.P., Yang, H. and Han, J.-C. (2009). ASME J. Heat "
             "Transfer 131(7), 071703. Eq. (17)/(18)";
  s.validity_source = s.source;
  s.provenance = RibProvenance::Extracted;
  // 45 deg square/sharp-edged ribs, same family as Han and Park's angled
  // ribs: reversing the flow gives the mirror image, same magnitude in a
  // symmetric duct. Not orthogonal (unlike han_1988_orthogonal), but still
  // symmetric -- see the case A/B/C reverse-flow discussion on #334.
  s.symmetric = true;

  // R = 1.13 (e/D)^-0.17 (p/e)^0.38. Raw e/D and p/e, no /10 normaliser --
  // unlike han_1988_orthogonal's R_pe, which divides by 10 because that is
  // how Han's Eq. 4.16 is printed. Two different papers, two different
  // normalisations; carrying han_1988_orthogonal's convention here would be
  // silently wrong by (p/e)^0.38 evaluated at the wrong reference.
  s.C_R = 1.13;
  s.R_eD = {-0.17, 1.0};
  s.R_pe = {0.38, 1.0};

  // G = 1.24 (e/D)^0.014 (p/e)^-0.02 (e+)^0.42. The exponents on e/D and p/e
  // are the resolved reading (extraction item 31): the textbook reprint at
  // Fig. 4.193c prints 0.14 for the e/D exponent, a lost decimal place
  // against the paper's own 0.014, worth 22% in G. Confirmed against the
  // paper directly, not the reprint -- see han_ribbed_high_re.md, D4.
  s.C_G = 1.24;
  s.G_eD = {0.014, 1.0};
  s.G_pe = {-0.02, 1.0};
  s.G_eplus_exponent = 0.42;

  s.valid_Re = {30000.0, 400000.0};
  s.valid_eD = {0.1, 0.18};
  s.valid_pe = {5.0, 10.0};
  // Square channel only (W/H = 1) and 45 deg ribs only: both are how the
  // source's own experiments were run, not modelling choices, so they are
  // recorded as point ranges the same way han_1988_orthogonal pins alpha.
  s.valid_WH = {1.0, 1.0};
  s.valid_alpha = {45.0, 45.0};
  // Derived, not separately stated: the e+ reached at the corners of the
  // Re x e/D x p/e validity box above, through this set's own R -> f -> e+
  // chain. The source's item-8 anchor (e+ = 18,000 at Re = 400K, e/D = 0.18)
  // falls inside this box; it does not fix p/e, which this box does span.
  s.valid_eplus = {542.0, 25340.0};
  s.valid_Pr = 0.7;  // unchanged from Han; the source does not restate it

  // accuracy_G is NOT an author-stated band, unlike han_1988_orthogonal's
  // (Han's own "95% within X%" claims). Rallabandi et al. state no percentage
  // accuracy for Eq. (18) in the extracted text. This is the RMS this
  // project measured, of Eq. (18) itself against 38 digitised points from
  // Fig. 4.193c (validation/cooling/data/han2012/fig4.193c_G_scatter.csv):
  // mean pred/data 0.994, RMS 6.9% -- see han_ribbed_high_re.md, items 27-31.
  // accuracy_R is left at 0 (unstated): no R data was digitised for this
  // range, so there is nothing to measure it against.
  s.accuracy_G = 0.069;
  return s;
}

RibCorrelationSet han_park_1988_angled() {
  RibCorrelationSet s;
  s.name = "han_park_1988_angled";
  s.source = "Han, J.C. and Park, J.S. (1988). IJHMT 31(1), 183. Eq. "
             "4.17/4.18, via Han, Dutta & Ekkad (2012) 2nd ed. Fig. 4.47";
  s.validity_source = s.source;
  s.provenance = RibProvenance::Extracted;
  // Parallel angled ribs: reversing the flow gives the mirror image, same
  // magnitude in a symmetric duct, with the secondary-flow direction
  // flipped -- not the reversal case that changes the configuration (a
  // V-shaped rib reversed becomes an inverted V, a different measured
  // geometry). See the case A/B/C reverse-flow discussion on #334.
  s.symmetric = true;

  // R/[(P/e/10)^0.35 (W/H)^m] = 12.31 - 27.07(alpha/90) + 17.86(alpha/90)^2
  // m = 0 at alpha = 90 deg (any W/H collapses to 1), m = 0.35 otherwise.
  // The leading constant is 12.31 (Fig. 4.47), not p. 376's printed 12.3 --
  // decision D2 in han_ribbed.md, accepted in favour of the figure.
  s.R_alpha_shape = RibCorrelationSet::RAlphaShape::QuadraticAlpha;
  s.R_quad_c0 = 12.31;
  s.R_quad_c1 = -27.07;
  s.R_quad_c2 = 17.86;
  s.R_pe = {0.35, 10.0};
  s.R_quad_WH_exponent_at_90 = 0.0;
  s.R_quad_WH_exponent_off_90 = 0.35;
  s.R_quad_WH_cap = 2.0;  // Eq. 4.17: "if W/H > 2, set W/H = 2"

  // G = 2.24 (W/H)^0.1 (alpha/90)^m (p/e/10)^n (e+)^0.35. (W/H)^0.1 sits
  // OUTSIDE the square/rectangular switch and applies always -- carried on
  // G_WH, which the switch does not touch. m, n switch on channel shape:
  // square (W/H=1) gets m=0.35, n=0.1; rectangular gets m=n=0, the stated
  // consequence being that alpha and P/e stop mattering for G once the
  // channel is not square.
  s.C_G = 2.24;
  s.G_WH = {0.1, 1.0};
  s.G_shape_model = RibCorrelationSet::GShapeModel::SquareVsRectangular;
  s.G_shape_alpha_exponent_square = 0.35;
  s.G_shape_alpha_exponent_rect = 0.0;
  s.G_shape_pe_exponent_square = 0.1;
  s.G_shape_pe_exponent_rect = 0.0;
  s.G_eplus_exponent = 0.35;

  s.valid_Re = {10000.0, 60000.0};
  s.valid_eD = {0.047, 0.078};
  s.valid_pe = {10.0, 20.0};
  s.valid_WH = {1.0, 4.0};
  s.valid_alpha = {30.0, 90.0};
  s.valid_eplus = {50.0, 0.0};  // e+ >= 50, same floor as han_1988_orthogonal
  s.valid_Pr = 0.7;

  // Neither accuracy is author-stated for this pair the way Han (1988)'s
  // own R and G are (6%/8%, "95% of data within"). These are the RMS this
  // project measured scoring THROUGH evaluate_rib against figure 4.47's
  // own digitised data (39 points for R, 115 for G) with a representative
  // geometry (the cloud carries no legend) -- see
  // validation/cooling/extractions/han_ribbed.md and
  // validation/cooling/data/han2012/fig4.47_*.csv.
  s.accuracy_R = 0.105;
  s.accuracy_G = 0.088;
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

// True for an angle within floating-point noise of exactly 90 degrees.
// alpha_deg is a fixed geometry parameter in every current caller, not a
// Newton-iterated state, so an exact-valued switch does not sit in a
// solver's search path today; see evaluate_rib's header comment for what
// changes if that ever stops being true.
bool is_alpha_90(double alpha_deg) {
  return std::abs(alpha_deg - 90.0) < 1e-9;
}

// True for a square channel, within floating-point noise of W/H = 1.
bool is_square_channel(double W_H) { return std::abs(W_H - 1.0) < 1e-9; }

double han_park_R(const RibCorrelationSet &set, double e_D, double W_H,
                  double p_e, double alpha_deg) {
  (void)e_D;  // Eq. 4.17 carries no e/D term.
  const double u = alpha_deg / 90.0;
  const double alpha_poly =
      set.R_quad_c0 + set.R_quad_c1 * u + set.R_quad_c2 * u * u;
  double W_H_capped = W_H;
  if (set.R_quad_WH_cap > 0.0 && W_H_capped > set.R_quad_WH_cap) {
    W_H_capped = set.R_quad_WH_cap;
  }
  const double m = is_alpha_90(alpha_deg) ? set.R_quad_WH_exponent_at_90
                                          : set.R_quad_WH_exponent_off_90;
  const double WH_term =
      (W_H_capped > 0.0 && std::isfinite(W_H_capped) && m != 0.0)
          ? std::pow(W_H_capped, m)
          : 1.0;
  return alpha_poly * power_term(set.R_pe, p_e) * WH_term;
}

}  // namespace

RibResult evaluate_rib(const RibCorrelationSet &set, const RibGeometry &geom,
                       double Re) {
  RibResult out;

  const double e_D = geom.e_D;
  const double W_H = geom.W_H;

  if (set.R_alpha_shape == RibCorrelationSet::RAlphaShape::QuadraticAlpha) {
    out.R = han_park_R(set, e_D, W_H, geom.p_e, geom.alpha_deg);
  } else {
    out.R =
        set.C_R * power_term(set.R_eD, e_D) * power_term(set.R_pe, geom.p_e) *
        power_term(set.R_WH, W_H) *
        power_term(set.R_alpha, geom.alpha_deg / 90.0 * set.R_alpha.reference);
  }

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

  double G_alpha_term, G_pe_term;
  if (set.G_shape_model == RibCorrelationSet::GShapeModel::SquareVsRectangular) {
    const bool square = is_square_channel(W_H);
    const double m = square ? set.G_shape_alpha_exponent_square
                            : set.G_shape_alpha_exponent_rect;
    const double n = square ? set.G_shape_pe_exponent_square
                            : set.G_shape_pe_exponent_rect;
    G_alpha_term = (m != 0.0) ? std::pow(geom.alpha_deg / 90.0, m) : 1.0;
    G_pe_term = (n != 0.0 && geom.p_e / 10.0 > 0.0)
                    ? std::pow(geom.p_e / 10.0, n)
                    : 1.0;
  } else {
    G_alpha_term =
        power_term(set.G_alpha, geom.alpha_deg / 90.0 * set.G_alpha.reference);
    G_pe_term = power_term(set.G_pe, geom.p_e);
  }

  out.G = set.C_G * power_term(set.G_eD, e_D) * G_pe_term *
          power_term(set.G_WH, W_H) * G_alpha_term *
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
  if (set.R_alpha_shape == RibCorrelationSet::RAlphaShape::QuadraticAlpha) {
    // C_R plays no part in this shape -- the quadratic fields replace it
    // entirely -- so it is exempt from the positivity check below rather
    // than requiring a caller to set a dummy value nothing reads.
    if (!std::isfinite(set.R_quad_c0) || !std::isfinite(set.R_quad_c1) ||
        !std::isfinite(set.R_quad_c2)) {
      throw std::invalid_argument("rib correlation set '" + set.name +
                                  "': R_quad_c0/c1/c2 must be finite");
    }
    if (!std::isfinite(set.R_quad_WH_exponent_at_90) ||
        !std::isfinite(set.R_quad_WH_exponent_off_90)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': R_quad_WH_exponent_at_90/off_90 must be finite");
    }
    if (!(set.R_quad_WH_cap >= 0.0) || !std::isfinite(set.R_quad_WH_cap)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': R_quad_WH_cap must be finite and non-negative (0 = uncapped)");
    }
  } else if (!(set.C_R > 0.0) || !std::isfinite(set.C_R)) {
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
  if (set.G_shape_model == RibCorrelationSet::GShapeModel::SquareVsRectangular) {
    if (!std::isfinite(set.G_shape_alpha_exponent_square) ||
        !std::isfinite(set.G_shape_alpha_exponent_rect) ||
        !std::isfinite(set.G_shape_pe_exponent_square) ||
        !std::isfinite(set.G_shape_pe_exponent_rect)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': G_shape_* exponents must be finite");
    }
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
