#pragma once

// Parametrised rib-roughness correlations.
//
// Rib correlations in the literature share one algebraic shape: a constant
// times power-law terms in the geometry, with the heat-transfer function
// carrying an extra term in the roughness Reynolds number.
//
//   R = C_R * (e/D / nD)^a * (p/e / nP)^b * (W/H / nW)^c * (alpha/90)^d
//   G = C_G * (same geometry terms) * (e+)^n
//
// Expressing them as data rather than code means a new correlation is a
// parameter set, and a user can supply their own -- which is the expected path
// for real hardware, since no published correlation is precise enough for a
// specific rig.
//
// THE NORMALISER IS PART OF THE DATA, NOT A CONVENTION. (p/e/10)^b and
// (p/e)^b are the same function with different constants: Han's
// R = 3.2 (p/e/10)^0.35 is identically 1.4294 (p/e)^0.35. Applying one
// source's constant under another's convention is wrong by a constant factor
// -- 10^0.35 = 2.24x here -- which never looks like a trend and so survives
// review. Every term therefore carries its own reference value.
//
// R CARRIES NO e+ TERM by construction, so the friction factor this family
// produces is independent of Reynolds number. That is true of both
// correlations in scope (Han 1988, Rallabandi et al. 2009) and it makes
// df/d(mdot) identically zero. A correlation whose R varies with e+ cannot be
// expressed here and would need its own implicit solve.

#include <string>

namespace combaero {
namespace cooling {

// How much a parameter set's numbers can be trusted, and on whose authority.
enum class RibProvenance {
  // An equation printed in a source, transcribed and checked against it.
  Extracted,
  // A form WE chose, fitted to data WE digitised from a published figure.
  // The exponents are ours; the source only supplied points.
  Fitted,
  // Supplied by the caller. Carries no claim at all.
  User,
};

// One power-law term. `reference` is what the variable is divided by before
// the exponent is applied; 1.0 means the raw value.
struct RibTerm {
  double exponent = 0.0;
  double reference = 1.0;
};

// An advisory range. Outside it the evaluator warns; it never refuses, because
// a user supplying their own coefficients knows their own hardware and a band
// belongs to the source's rig rather than to theirs.
struct RibRange {
  double lo = 0.0;
  double hi = 0.0;
  bool bounded() const { return hi > lo; }
};

struct RibCorrelationSet {
  std::string name;
  // Where the coefficients came from.
  std::string source;
  // Where the VALIDITY came from. Usually the same, but not always: a set may
  // borrow a range from a later study that showed the correlation holds
  // further than its author claimed.
  std::string validity_source;
  RibProvenance provenance = RibProvenance::User;

  // False when reversing the flow changes the configuration rather than only
  // its direction -- a forward V-rib becomes an inverted V, which the
  // literature measures as a different geometry. True for 90 deg orthogonal
  // and parallel angled ribs, where reversal is a mirror operation.
  bool symmetric = true;

  // Friction roughness function. Two shapes exist in the sources:
  //
  //   PowerLaw       R = C_R * (e/D)^.. * (p/e)^.. * (W/H)^.. * (alpha/90)^d
  //                  (Han 1988, Rallabandi 2009 -- both have d = 0)
  //
  //   QuadraticAlpha R / [(p/e/10)^0.35 * (W/H)^m] = R_quad_c0 +
  //                      R_quad_c1*(alpha/90) + R_quad_c2*(alpha/90)^2
  //                  where m switches on the rib angle itself:
  //                      m = R_quad_WH_exponent_at_90   if alpha == 90 deg
  //                      m = R_quad_WH_exponent_off_90  otherwise
  //                  (Han and Park 1988, Eq. 4.17). R_eD, R_pe, R_WH and
  //                  R_alpha are IGNORED in this shape -- the quadratic
  //                  fields replace them entirely, they do not compose.
  //
  // The switch on m is a genuine discontinuity the source states, not a
  // numerical artefact: at W/H = 2 it is a 27.5% jump in R exactly at
  // alpha = 90, at W/H = 4 it is 62.5%. No smooth interpolation is given by
  // the source, so none is invented here -- see evaluate_rib's comment for
  // what that means for a solver that traverses this exact angle.
  enum class RAlphaShape { PowerLaw, QuadraticAlpha };
  RAlphaShape R_alpha_shape = RAlphaShape::PowerLaw;
  double C_R = 0.0;
  RibTerm R_eD, R_pe, R_WH, R_alpha;
  double R_quad_c0 = 0.0, R_quad_c1 = 0.0, R_quad_c2 = 0.0;
  double R_quad_WH_exponent_at_90 = 0.0;
  double R_quad_WH_exponent_off_90 = 0.0;
  // Eq. 4.17's own cap: "if W/H > 2, set W/H = 2". 0 means uncapped.
  double R_quad_WH_cap = 0.0;

  // Heat-transfer roughness function. G_alpha and G_pe are either fixed
  // constants (Fixed, the default -- Han 1988, Rallabandi 2009) or switch
  // on whether the channel is square (Han and Park 1988, Eq. 4.18):
  //
  //   m = G_shape_alpha_exponent_square, n = G_shape_pe_exponent_square
  //       if W/H == 1 (square)
  //   m = G_shape_alpha_exponent_rect,   n = G_shape_pe_exponent_rect
  //       otherwise (rectangular)
  //
  // G_alpha and G_pe are IGNORED when G_shape_model is SquareVsRectangular.
  // Same caveat as R's switch: a genuine, unsmoothed discontinuity at
  // W/H = 1, up to 27% at alpha = 30, p/e = 20.
  enum class GShapeModel { Fixed, SquareVsRectangular };
  GShapeModel G_shape_model = GShapeModel::Fixed;
  double C_G = 0.0;
  RibTerm G_eD, G_pe, G_WH, G_alpha;
  double G_eplus_exponent = 0.0;
  double G_shape_alpha_exponent_square = 0.0, G_shape_alpha_exponent_rect = 0.0;
  double G_shape_pe_exponent_square = 0.0, G_shape_pe_exponent_rect = 0.0;

  // Advisory validity.
  RibRange valid_Re, valid_eD, valid_pe, valid_WH, valid_alpha, valid_eplus;
  double valid_Pr = 0.0;  // 0 means unstated

  // Stated accuracy, as a fraction (0.06 for "within 6%"). 0 means unstated.
  double accuracy_R = 0.0;
  double accuracy_G = 0.0;
};

// Han, J.C. (1988), ASME J. Heat Transfer 110, 321, for 90 deg orthogonal ribs
// in two-opposite-wall rectangular channels. Extracted and confirmed; see
// validation/cooling/extractions/han_ribbed.md.
RibCorrelationSet han_1988_orthogonal();

// Rallabandi, A.P., Yang, H. and Han, J.-C. (2009), ASME J. Heat Transfer
// 131(7), 071703, for 45 deg square/sharp-edged ribs in a square channel at
// Reynolds numbers an order of magnitude above han_1988_orthogonal's range --
// "typical of land-based turbines" per the source. A separate regime, not a
// revision: different experiments, and per the source the correlations do
// not agree with Han's in the extended range. Selected explicitly, never
// blended -- see validation/cooling/extractions/han_ribbed_high_re.md.
//
// Sharp-edged ribs only. The source reports round-edged ribs at the same
// conditions instead following Han's correlation, but "coincidentally" --
// its own word, meaning reported as fortuitous rather than physically
// grounded, so it is not encoded as a switch here. A caller with round-edged
// ribs above this Reynolds range chooses the correlation directly, same as
// any other set choice.
RibCorrelationSet rallabandi_2009_high_re();

// Han, J.C. and Park, J.S. (1988), IJHMT 31(1), 183, Eq. 4.17/4.18, for
// broad-aspect-ratio rectangular ducts with angled ribs -- alpha 30-90 deg,
// W/H 1-4. A different paper and configuration from han_1988_orthogonal
// (which is 90 deg only): the two R correlations happen to agree closely
// at alpha=90 (3.44%), but the two G correlations disagree by up to 22%
// there, because Fig. 4.46 and this set describe genuinely different rib
// configurations that are not obliged to agree. See
// validation/cooling/extractions/han_ribbed.md, item 23.
RibCorrelationSet han_park_1988_angled();

// Geometry of the ribbed channel, as the correlation sees it.
struct RibGeometry {
  double e_D = 0.0;        // rib height / hydraulic diameter
  double p_e = 0.0;        // rib pitch / rib height
  double W_H = 1.0;        // channel width / height
  double alpha_deg = 90.0; // rib angle to the flow
};

// The chain, evaluated. Every field is finite for every real input: the
// evaluator never throws from inside a residual, because the solver
// legitimately probes states that are not physical and a throw there kills the
// solve. Bad PARAMETERS are rejected once, by validate_rib_set; bad OPERATING
// POINTS are guarded smoothly.
struct RibResult {
  double R = 0.0;        // friction roughness function
  double f = 0.0;        // four-sided channel friction factor
  double e_plus = 0.0;   // roughness Reynolds number, signed with the flow
  double G = 0.0;        // heat-transfer roughness function
  double St_r = 0.0;     // ribbed-side Stanton number
  // d(St_r)/d(Re). df/dRe is identically zero for this family -- R carries no
  // e+ term -- so it is not reported.
  double dSt_dRe = 0.0;
  bool extrapolated = false;  // outside the set's advisory validity
};

// Evaluate at a Reynolds number. `Re` may be negative (reverse flow) or zero;
// the guards are smooth through both, so the derivative is continuous where
// Newton iterates. `geom.alpha_deg` is NOT guarded the same way: a set's R/G
// terms are polynomials in alpha with no built-in symmetry, so an angle
// outside `valid_alpha` is evaluated raw rather than folded or clamped. See
// validation/cooling/extractions/han_ribbed.md, "Eq. 4.17 outside its stated
// 30-90 degree range" for the numerical characterisation this rests on: the
// raw polynomial stays finite and positive well past the fitted range, but
// its vertex sits at 68 deg (not 90), so it is unbounded and physically
// unmotivated beyond 30-90, not a safe extrapolation. No guard exists yet;
// if one is added, reflection-folding (alpha_eff = 180 - alpha) is the
// bounded, `symmetric`-consistent choice documented there, not the raw value.
RibResult evaluate_rib(const RibCorrelationSet &set, const RibGeometry &geom,
                       double Re);

// Reject a set that cannot be evaluated. Throws with a message naming the
// field. Unlike the validity ranges, these are hard errors: they are mistakes,
// not operating points.
void validate_rib_set(const RibCorrelationSet &set);

}  // namespace cooling
}  // namespace combaero
