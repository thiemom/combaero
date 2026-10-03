#pragma once

// Rib correlations in RATIO form: Nu/Nu0 and f/f0 as multipliers on a declared
// smooth-duct baseline, beside the R/G law-of-the-wall sets in
// rib_correlation.h (#444).
//
// WHY A SECOND FORM. R and G rest on local wall-normal shear similarity, which
// weakens once angled, V or crossed ribs drive counter-rotating secondary
// cells. Most independent sources publish Nu/Nu0 and f/f0 directly, and
// converting them to R/G to score them would impose the very assumption they
// could test.
//
// THE BASELINE IS PART OF THE DATA. A ratio means something only against the
// reference it was fitted to (Dittus-Boelter for Han and Taslim; 0.0176 Re^0.8
// for NASA CR-4396). Inside the fitted Re range the set evaluates
//
//     Nu = r(Re, geometry) * Nu0_source(Re, Pr)
//
// which reproduces the paper whatever baseline a caller prefers: re-referencing
// exactly to another baseline returns the same Nu, so the choice cancels.
//
// BELOW THE FITTED RANGE the source baseline is never evaluated alone. Over
// one factor of RATIO_BLEND_FACTOR below Re_floor the evaluator blends, C1 in
// ln Re, to an EXTRAPOLATION baseline scaled to match at the floor:
//
//     Nu = k * Nu0_ext(Re),   k = r(Re_floor) Nu0_source(Re_floor)
//                                 / Nu0_ext(Re_floor)
//
// Nu0_ext defaults to smooth-pipe Gnielinski, so Nu -> k * 3.66 (its laminar
// limit) as Re -> 0 instead of Dittus-Boelter's collapse to 0. This is a
// transfer, not a fit. Friction holds its ratio at the floor value instead.
//
// ROBUSTNESS, for shipped and user-supplied sets alike:
//   * Nu0 is never divided by: Nu is always r * Nu0 or k * Nu0_ext.
//   * Re enters as the smoothed magnitude sqrt(Re^2 + eps^2): Nu and f are
//     even in Re, their derivatives odd and continuous, all finite at Re = 0.
//   * Parameters that can only be mistakes are rejected once by
//     validate_rib_ratio_set; every operating point returns finite,
//     non-negative values -- the evaluator never throws.
//   * dNu/dRe and df/dRe are analytic (solver (f, J) rule).

#include <string>

#include "rib_correlation.h"

namespace combaero {
namespace cooling {

// Below-floor blend width: the source and extrapolation forms are blended over
// [Re_floor / RATIO_BLEND_FACTOR, Re_floor], C1 in ln Re.
constexpr double RATIO_BLEND_FACTOR = 2.0;

// Smoothing scale for |Re| -- negligible above Re ~ 100, keeps every power of
// Re finite and every derivative continuous through Re = 0.
constexpr double RATIO_RE_EPS = 1.0;

// A smooth-duct baseline of power-law form: value = coeff Re^re_exponent
// Pr^pr_exponent. Dittus-Boelter heating is {0.023, 0.8, 0.4}; CR-4396's
// Kays-Perkins reference {0.0176, 0.8, 0}; Blasius-form friction
// {0.046, -0.2, 0}.
struct RatioBaseline {
  double coeff = 0.0;
  double re_exponent = 0.0;
  double pr_exponent = 0.0;
};

// What Nu hands over to below the fitted range.
enum class RatioBelowFloor {
  // Smooth-pipe Gnielinski (Petukhov friction): laminar-limited at Re -> 0.
  Gnielinski,
  // Hold the ratio at its floor value on the source baseline. With a
  // Dittus-Boelter source this tends to Nu = 0 as Re -> 0.
  SourceBaseline,
  // The caller's own power-law Nu0 (RibRatioOptions::user_Nu0).
  User,
};

struct RibRatioOptions {
  RatioBelowFloor below_floor = RatioBelowFloor::Gnielinski;
  RatioBaseline user_Nu0;  // used only with RatioBelowFloor::User
};

struct RibRatioSet {
  std::string name;
  std::string source;
  std::string validity_source;
  RibProvenance provenance = RibProvenance::User;
  RibShape shape = RibShape::Unspecified;
  bool symmetric = true;

  // Nu/Nu0 = C_Nu * (Re/ref)^a * (e/D/ref)^.. * (p/e/ref)^.. * (W/H/ref)^..
  //          * (alpha/90)^..    -- RibTerm carries each exponent and its
  //          normaliser, exactly as in RibCorrelationSet.
  double C_Nu = 0.0;
  RibTerm Nu_Re, Nu_eD, Nu_pe, Nu_WH, Nu_alpha;
  // f/f0, same form. f has the definition of f0_source (Fanning or Darcy).
  double C_f = 0.0;
  RibTerm f_Re, f_eD, f_pe, f_WH, f_alpha;

  // The baselines the source normalised by. Required: a ratio without its
  // reference is not a number.
  RatioBaseline Nu0_source;
  RatioBaseline f0_source;

  // Bottom of the fitted Re range. Below it the handover applies.
  double Re_floor = 0.0;

  RibRange valid_Re, valid_eD, valid_pe, valid_WH, valid_alpha;
  double valid_Pr = 0.0;  // 0 means unstated
  StatedAccuracy accuracy_Nu, accuracy_f;
};

struct RibRatioResult {
  double Nu = 0.0;        // ribbed-surface Nusselt number
  double dNu_dRe = 0.0;   // analytic, odd in Re
  double f = 0.0;         // friction factor, f0_source's definition
  double df_dRe = 0.0;    // analytic, odd in Re
  double ratio_Nu = 0.0;  // the set's Nu multiplier at |Re| (before handover)
  double ratio_f = 0.0;   // the friction multiplier actually applied
  bool below_floor = false;   // |Re| < Re_floor: the handover is active
  bool extrapolated = false;  // outside advisory validity, or below floor
};

// Evaluate at a Reynolds number and Prandtl number. Never throws: Re may be
// negative or zero, Pr is guarded smoothly.
RibRatioResult evaluate_rib_ratio(const RibRatioSet &set,
                                  const RibGeometry &geom, double Re,
                                  double Pr,
                                  const RibRatioOptions &options = {});

// Reject a set that cannot be evaluated. Hard errors, named by field.
void validate_rib_ratio_set(const RibRatioSet &set);
void validate_rib_ratio_options(const RibRatioOptions &options);

}  // namespace cooling
}  // namespace combaero
