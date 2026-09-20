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

  // Friction roughness function.
  double C_R = 0.0;
  RibTerm R_eD, R_pe, R_WH, R_alpha;

  // Heat-transfer roughness function.
  double C_G = 0.0;
  RibTerm G_eD, G_pe, G_WH, G_alpha;
  double G_eplus_exponent = 0.0;

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
