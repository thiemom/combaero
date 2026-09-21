#pragma once

// Jet impingement heat transfer correlations.
//
// Two independent regimes, matching Han, Dutta and Ekkad (2012) Section 4.1's
// own split -- see validation/cooling/extractions/han_impingement.md for the
// full extraction and review this is built from.
//
//   Single jet: Goldstein, Behbahani and Heppelmann (1986). One free round
//   jet impinging on a flat plate. Nu_bar is a function of jet Reynolds
//   number, jet-to-plate spacing (L/D) and radial distance from the jet
//   centerline (R/D). No crossflow, no array.
//
//   Jet array with crossflow: Florschuetz, Truman and Metzger (1981). A
//   staggered or inline array of jets fed from a common plenum, where every
//   row's heat transfer is degraded by the crossflow (spent air) from every
//   row upstream of it. The defining input is the crossflow-to-jet mass flux
//   ratio Gc/Gj at the row of interest, which this header also derives in
//   closed form -- the paper's own Eq. (8) -- so a caller supplies geometry
//   and a discharge coefficient rather than a pre-computed ratio.
//
// SCOPE. This covers exactly what the extraction confirmed as ready: the two
// Nu correlations and the Gc/Gj closed form. It deliberately does NOT cover
// converting an array's mean/total mass flow into a per-row local Re_j (the
// paper's Eq. (7), the jet velocity distribution) -- that is a network
// element's job (how a total coolant flow splits across N holes and Nc rows
// is a design/wiring question, not a correlation), and is not yet decided.
// Leading-edge/curved-surface impingement (Section 4.1.4) and the simpler,
// looser forms (Eq. 4.6 Kercher-Tabakoff -- graphical, not closed-form
// anyway; Eq. 4.7/4.8 -- Florschuetz's own less-tight alternate) are also
// out of scope, deferred per the extraction's I3.

#include <string>

namespace combaero {
namespace cooling {

// An advisory range. Outside it the evaluator warns; it never refuses,
// mirroring RibRange in rib_correlation.h -- a duplicate of that tiny struct
// rather than a shared dependency between the two correlation families,
// which stay independent by design.
struct ImpingementRange {
  double lo = 0.0;
  double hi = 0.0;
  bool bounded() const { return hi > lo; }
};

// ---------------------------------------------------------------
// Single jet (Goldstein, Behbahani and Heppelmann, 1986)
// ---------------------------------------------------------------

// Which surface boundary condition the correlation's exponent was fitted
// under. The two branches share A, B, C and differ only in this exponent --
// see han_impingement.md item 1.
enum class ImpingementThermalBC { ConstantHeatFlux, ConstantWallTemperature };

struct SingleJetImpingementSet {
  std::string name;
  std::string source;

  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  // Re exponent is FIXED at 0.76 for both branches -- it is only the (R/D)
  // exponent below that switches with the boundary condition. Getting this
  // backwards (applying n to Re instead of R/D) misses the closed-form
  // check point by two orders of magnitude; see the test that pins it.
  double Re_exponent = 0.0;
  double n_const_heat_flux = 0.0;
  double n_const_wall_temp = 0.0;
};

// Goldstein, Behbahani and Heppelmann (1986), IJHMT 29(8), 1227-1235, as
// reprinted in Han, Dutta and Ekkad (2012) Eq. 4.1. Extracted and confirmed
// against two independent channels (page image and an independent OCR
// extraction); see han_impingement.md items 1-7.
SingleJetImpingementSet goldstein_1986_single_jet();

// Nu_bar for a single round jet impinging on a flat plate:
//   Nu_bar = Re^0.76 * (A - |L/D - 7.75|) / (B + C*(R/D)^n)
// with n = n_const_heat_flux or n_const_wall_temp depending on `bc`.
//
// Re: jet Reynolds number.
// L_D: jet-to-target-plate spacing / jet diameter.
// R_D: radial distance from the jet centerline / jet diameter.
//
// No stated validity range exists in the source beyond the closed-form check
// point (Re=25000, R/D=5, L/D=7.75 -> Nu=60 constant-heat-flux, 56
// constant-wall-temperature) -- see han_impingement.md item 4. The
// |L/D - 7.75| term is not floored: it goes negative for L/D far from the
// optimum, exactly as the source's own formula does, because Han's book
// states no floor and none is invented here (same policy as Eq. 4.17's
// unguarded angle range in rib_correlation.h).
double single_jet_impingement_nu(const SingleJetImpingementSet &set,
                                  ImpingementThermalBC bc, double Re,
                                  double L_D, double R_D);

// ---------------------------------------------------------------
// Jet array with crossflow (Florschuetz, Truman and Metzger, 1981)
// ---------------------------------------------------------------

enum class JetHolePattern { Inline, Staggered };

// One of Table 4.1's four geometry-dependent quantities (A, m, B or n):
// value = C * (xn/d)^nx * (yn/d)^ny * (z/d)^nz
struct JetArrayGeometricFit {
  double C = 0.0;
  double nx = 0.0;
  double ny = 0.0;
  double nz = 0.0;
};

struct JetArrayCorrelationSet {
  std::string name;
  std::string source;
  std::string validity_source;
  JetHolePattern pattern = JetHolePattern::Inline;

  JetArrayGeometricFit A_fit, m_fit, B_fit, n_fit;

  // Advisory validity, from the primary paper's own stated overall ranges
  // (han_impingement.md item 20) -- NOT the misprinted rib-validity box
  // Han's book carries under Table 4.1 (item 19/I1).
  ImpingementRange valid_Re_j, valid_Gc_Gj, valid_xn_d, valid_yn_d, valid_z_d;
  ImpingementRange valid_aspect_ratio;  // xn/yn

  // Standard error of the fit (0.056 for inline, 0.061 for staggered), NOT
  // the experimental measurement uncertainty (+-5% at 95% confidence) --
  // see han_impingement.md items 12 and 23, which are two different numbers.
  double standard_error = 0.0;
};

// Florschuetz, Truman and Metzger (1981), ASME J. Heat Transfer 103,
// 337-342, Eq. (10a)/(10b)/Table 2 (Han's Eq. 4.9/Table 4.1). Extracted and
// confirmed against the primary paper directly, page-image by page-image;
// see han_impingement.md items 14-20, 23.
JetArrayCorrelationSet florschuetz_1981_inline();
JetArrayCorrelationSet florschuetz_1981_staggered();

// The primary paper's own recommended default jet-plate discharge
// coefficient, "for jet plates similar to those utilized here", absent a
// measured value (han_impingement.md item 24). Measured values in the
// paper's own Table 1 range 0.73-0.85 depending on configuration. combaero
// has no correlation of its own for this yet -- see #375.
constexpr double FLORSCHUETZ_1981_DEFAULT_CD = 0.79;

struct JetArrayImpingementResult {
  double Nu = 0.0;
  bool extrapolated = false;  // outside the set's advisory validity
};

// Nu at one spanwise row, from Eq. (9)/Table 2 (Han's Eq. 4.9/Table 4.1).
//
// Re_j, Gc_Gj: the ROW's OWN local values, not an array mean -- the paper
// correlates "individual spanwise row jet Reynolds number" and the crossflow
// ratio immediately upstream of that row (see crossflow_to_jet_ratio_at_row
// below). Pr is the coolant's Prandtl number.
JetArrayImpingementResult jet_array_impingement_nu(
    const JetArrayCorrelationSet &set, double Re_j, double Gc_Gj, double Pr,
    double xn_d, double yn_d, double z_d);

// ---------------------------------------------------------------
// Crossflow-to-jet mass flux ratio (Florschuetz Eq. (8))
// ---------------------------------------------------------------
//
// A one-dimensional continuous-injection model (the paper's Eqs. 1-8,
// verified there against directly measured pressure-traverse data), reduced
// to its row-resolved closed form. Depends on (yn/d)(z/d) only -- NOT xn/d --
// which the paper states explicitly: "the flow distribution is independent
// of the streamwise hole spacing and hole pattern". See
// han_impingement.md's "The Gc/Gj closed form" section for the full
// derivation and why xn/d drops out.

// Gc/Gj at continuous streamwise position x/xn (x measured from the
// upstream end of the channel, in units of the streamwise hole spacing).
// C_D is the jet-plate discharge coefficient (FLORSCHUETZ_1981_DEFAULT_CD
// absent a measured value).
double crossflow_to_jet_ratio_at_x(double yn_d, double z_d, double C_D,
                                   double x_over_xn);

// Convenience form at discrete spanwise row `row`, 1-indexed counting from
// upstream, matching x = xn*(row - 1/2) -- how the paper itself applies
// Eq. (9)/Table 2, "Nusselt numbers resolved to one streamwise hole
// spacing". Row 1 always returns exactly 0 (Gc/Gj = 0 at the first row,
// Nu1's own definition), independent of geometry or C_D.
double crossflow_to_jet_ratio_at_row(double yn_d, double z_d, double C_D,
                                     int row);

// Reject a set that cannot be evaluated. Throws with a message naming the
// field. Unlike the validity ranges, these are hard errors: they are
// mistakes, not operating points.
void validate_single_jet_set(const SingleJetImpingementSet &set);
void validate_jet_array_set(const JetArrayCorrelationSet &set);

}  // namespace cooling
}  // namespace combaero
