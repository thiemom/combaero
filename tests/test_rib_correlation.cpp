#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "rib_correlation.h"

using combaero::cooling::evaluate_rib;
using combaero::cooling::han_1988_orthogonal;
using combaero::cooling::han_park_1988_angled;
using combaero::cooling::rallabandi_2009_high_re;
using combaero::cooling::RibCorrelationSet;
using combaero::cooling::RibGeometry;
using combaero::cooling::validate_rib_set;

namespace {

RibGeometry ref_geometry() {
  RibGeometry g;
  g.e_D = 0.047;
  g.p_e = 10.0;
  g.W_H = 1.0;
  g.alpha_deg = 90.0;
  return g;
}

RibGeometry rallabandi_ref_geometry() {
  RibGeometry g;
  g.e_D = 0.1;
  g.p_e = 10.0;
  g.W_H = 1.0;
  g.alpha_deg = 45.0;
  return g;
}

}  // namespace

// The numbers below are the worked values recorded in the CONFIRMED extraction
// at validation/cooling/extractions/han_ribbed.md. They come from the document,
// not from this implementation, so this test can fail.
TEST(RibCorrelationTest, ReproducesTheConfirmedExtraction) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();

  const auto lo = evaluate_rib(set, g, 10000.0);
  EXPECT_NEAR(lo.R, 3.2000, 1e-4);
  EXPECT_NEAR(lo.f, 0.04576, 1e-5);
  EXPECT_NEAR(lo.e_plus, 71.1, 0.1);
  EXPECT_NEAR(lo.G, 12.21, 0.01);
  EXPECT_NEAR(lo.St_r, 0.00968, 1e-5);

  const auto hi = evaluate_rib(set, g, 60000.0);
  EXPECT_NEAR(hi.e_plus, 426.6, 0.1);
  EXPECT_NEAR(hi.G, 20.16, 0.01);
  EXPECT_NEAR(hi.St_r, 0.00642, 1e-5);
}

// R carries no e+ term, so f cannot depend on Reynolds number. This is a
// structural property of the family, and it is what makes df/d(mdot) zero.
TEST(RibCorrelationTest, FrictionFactorIsIndependentOfReynoldsNumber) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();
  const double f_ref = evaluate_rib(set, g, 10000.0).f;
  for (double Re : {1.0, 1e3, 1e5, 1e7}) {
    EXPECT_NEAR(evaluate_rib(set, g, Re).f, f_ref, 1e-12);
  }
}

// The solver probes states that are not physical. Every one of these must
// return finite numbers rather than NaN, an infinity, or a complex result.
TEST(RibCorrelationTest, GuardsHoldForStatesTheSolverActuallyProbes) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();
  for (double Re : {-1e6, -1e4, -1.0, 0.0, 1e-12, 1.0, 1e8}) {
    const auto r = evaluate_rib(set, g, Re);
    EXPECT_TRUE(std::isfinite(r.f)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.G)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.St_r)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.dSt_dRe)) << "Re = " << Re;
    EXPECT_GT(r.f, 0.0) << "Re = " << Re;
    EXPECT_GT(r.St_r, 0.0) << "Re = " << Re;
  }
}

// Reversing the flow must not change the magnitude of a heat-transfer
// quantity: h does not care which way the fluid goes. Only the sign of e+
// records direction.
TEST(RibCorrelationTest, ReverseFlowIsSymmetricInMagnitude) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();
  const auto fwd = evaluate_rib(set, g, 30000.0);
  const auto rev = evaluate_rib(set, g, -30000.0);

  EXPECT_NEAR(fwd.G, rev.G, 1e-12);
  EXPECT_NEAR(fwd.St_r, rev.St_r, 1e-12);
  EXPECT_NEAR(fwd.e_plus, -rev.e_plus, 1e-9);
  // The derivative of an even function flips sign across zero.
  EXPECT_NEAR(fwd.dSt_dRe, -rev.dSt_dRe, 1e-12);
}

// The analytic derivative must match a central difference THROUGH the guard,
// not only in the far field where the guard is inert. A pin-fin Jacobian that
// was right everywhere except where it mattered is exactly what #328 fixed.
TEST(RibCorrelationTest, StantonDerivativeMatchesCentralDifferences) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();
  // The last two sit where the e+ floor is active.
  for (double Re : {50000.0, 10000.0, 100.0, 10.0}) {
    const double h = std::max(1e-6, std::abs(Re) * 1e-6);
    const double fd = (evaluate_rib(set, g, Re + h).St_r -
                       evaluate_rib(set, g, Re - h).St_r) /
                      (2.0 * h);
    const double analytic = evaluate_rib(set, g, Re).dSt_dRe;
    const double scale = std::max({std::abs(fd), std::abs(analytic), 1e-14});
    EXPECT_LT(std::abs(analytic - fd) / scale, 1e-5) << "Re = " << Re;
  }
}

// Outside the source's ranges the evaluator reports, it does not refuse: the
// band belongs to the source's rig, not to the caller's hardware.
TEST(RibCorrelationTest, ValidityIsAdvisoryNotEnforced) {
  const auto set = han_1988_orthogonal();
  auto g = ref_geometry();
  EXPECT_FALSE(evaluate_rib(set, g, 30000.0).extrapolated);

  g.e_D = 0.20;  // well past the 0.047-0.078 band
  const auto r = evaluate_rib(set, g, 30000.0);
  EXPECT_TRUE(r.extrapolated);
  EXPECT_TRUE(std::isfinite(r.St_r));
  EXPECT_GT(r.St_r, 0.0);
}

// Bad PARAMETERS are a mistake and are rejected once, loudly -- the opposite
// treatment from a bad operating point.
TEST(RibCorrelationTest, MalformedSetsAreRejected) {
  EXPECT_NO_THROW(validate_rib_set(han_1988_orthogonal()));

  auto no_source = han_1988_orthogonal();
  no_source.source.clear();
  EXPECT_THROW(validate_rib_set(no_source), std::invalid_argument);

  auto zero_reference = han_1988_orthogonal();
  zero_reference.R_pe.reference = 0.0;
  EXPECT_THROW(validate_rib_set(zero_reference), std::invalid_argument);

  auto nan_constant = han_1988_orthogonal();
  nan_constant.C_G = std::nan("");
  EXPECT_THROW(validate_rib_set(nan_constant), std::invalid_argument);

  auto negative = han_1988_orthogonal();
  negative.C_R = -1.0;
  EXPECT_THROW(validate_rib_set(negative), std::invalid_argument);
}

// The numbers below are hand-computed from Eq. (17)/(18) as printed in
// Rallabandi, Yang and Han (2009), independently of this implementation --
// see validation/cooling/extractions/han_ribbed_high_re.md. R carries no e+
// term here either, so it repeats across Re exactly as for han_1988_orthogonal.
TEST(RibCorrelationTest, ReproducesTheConfirmedRallabandiExtraction) {
  const auto set = rallabandi_2009_high_re();
  const auto g = rallabandi_ref_geometry();

  const auto lo = evaluate_rib(set, g, 30000.0);
  EXPECT_NEAR(lo.R, 4.00939, 1e-4);
  EXPECT_NEAR(lo.f, 0.06533, 1e-4);
  EXPECT_NEAR(lo.e_plus, 542.20, 0.1);
  EXPECT_NEAR(lo.G, 16.1350, 0.01);
  EXPECT_NEAR(lo.St_r, 0.010235, 1e-5);

  const auto hi = evaluate_rib(set, g, 400000.0);
  EXPECT_NEAR(hi.e_plus, 7229.37, 0.1);
  EXPECT_NEAR(hi.G, 47.8899, 0.01);
  EXPECT_NEAR(hi.St_r, 0.003658, 1e-5);
}

// han_1988_orthogonal (90 deg) and rallabandi_2009_high_re (45 deg, sharp
// ribs) are different papers a caller must choose between explicitly. Their
// Reynolds-number bands actually OVERLAP (10e3-60e3 vs 30e3-400e3, sharing
// 30e3-60e3) -- Re alone does not distinguish them, which is exactly why
// there is no auto-switching on Re. What is genuinely disjoint is alpha and
// e/D, confirmed here so a future edit cannot silently narrow either band
// into the other's territory without this test noticing.
TEST(RibCorrelationTest, TheTwoSetsDoNotOverlapInAlphaOrEd) {
  const auto han = han_1988_orthogonal();
  const auto rallabandi = rallabandi_2009_high_re();

  EXPECT_EQ(han.valid_alpha.lo, 90.0);
  EXPECT_EQ(han.valid_alpha.hi, 90.0);
  EXPECT_EQ(rallabandi.valid_alpha.lo, 45.0);
  EXPECT_EQ(rallabandi.valid_alpha.hi, 45.0);

  // Rallabandi's study was a square channel only.
  EXPECT_EQ(rallabandi.valid_WH.lo, 1.0);
  EXPECT_EQ(rallabandi.valid_WH.hi, 1.0);

  // e/D bands do not overlap: Han's is small ribs, Rallabandi's is large.
  EXPECT_LT(han.valid_eD.hi, rallabandi.valid_eD.lo);

  // Their Re bands DO overlap -- recorded as a fact, not a defect. See the
  // comment above: this is why the two sets are named and chosen, not
  // switched on Reynolds number.
  EXPECT_GT(han.valid_Re.hi, rallabandi.valid_Re.lo);
}

// Reversing the flow must not change the magnitude, same reasoning as the
// 90 deg set: a 45 deg parallel rib reversed is its mirror image.
TEST(RibCorrelationTest, RallabandiReverseFlowIsSymmetricInMagnitude) {
  const auto set = rallabandi_2009_high_re();
  const auto g = rallabandi_ref_geometry();
  const auto fwd = evaluate_rib(set, g, 100000.0);
  const auto rev = evaluate_rib(set, g, -100000.0);

  EXPECT_NEAR(fwd.G, rev.G, 1e-9);
  EXPECT_NEAR(fwd.St_r, rev.St_r, 1e-9);
  EXPECT_NEAR(fwd.e_plus, -rev.e_plus, 1e-6);
}

// Same guard discipline as the 90 deg set: every state a solver can probe
// must return finite, positive numbers, never NaN or a complex result.
TEST(RibCorrelationTest, RallabandiGuardsHoldForStatesTheSolverActuallyProbes) {
  const auto set = rallabandi_2009_high_re();
  const auto g = rallabandi_ref_geometry();
  for (double Re : {-1e6, -1e4, -1.0, 0.0, 1e-12, 1.0, 1e8}) {
    const auto r = evaluate_rib(set, g, Re);
    EXPECT_TRUE(std::isfinite(r.f)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.G)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.St_r)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.dSt_dRe)) << "Re = " << Re;
    EXPECT_GT(r.f, 0.0) << "Re = " << Re;
    EXPECT_GT(r.St_r, 0.0) << "Re = " << Re;
  }
}

// The analytic Stanton derivative must match a central difference through
// this set's own guards too -- not assumed to inherit correctness from the
// shared evaluate_rib code path just because han_1988_orthogonal's does.
TEST(RibCorrelationTest, RallabandiStantonDerivativeMatchesCentralDifferences) {
  const auto set = rallabandi_2009_high_re();
  const auto g = rallabandi_ref_geometry();
  for (double Re : {200000.0, 50000.0, 100.0, 10.0}) {
    const double h = std::max(1e-6, std::abs(Re) * 1e-6);
    const double fd = (evaluate_rib(set, g, Re + h).St_r -
                       evaluate_rib(set, g, Re - h).St_r) /
                      (2.0 * h);
    const double analytic = evaluate_rib(set, g, Re).dSt_dRe;
    const double scale = std::max({std::abs(fd), std::abs(analytic), 1e-14});
    EXPECT_LT(std::abs(analytic - fd) / scale, 1e-5) << "Re = " << Re;
  }
}

TEST(RibCorrelationTest, RallabandiSetPassesValidation) {
  EXPECT_NO_THROW(validate_rib_set(rallabandi_2009_high_re()));
}


// Hand-computed independently (Python, from Eq. 4.17/4.18 as printed), not
// from this implementation -- see
// validation/cooling/extractions/han_ribbed.md. e+ agrees only after
// applying the documented smooth floor (sqrt(e+^2 + 1)); a raw-e+
// computation differs at the 5th decimal in G by construction, not error.
TEST(RibCorrelationTest, ReproducesTheConfirmedHanParkExtraction) {
  const auto set = han_park_1988_angled();

  RibGeometry square90;
  square90.e_D = 0.0625;
  square90.p_e = 15.0;
  square90.W_H = 1.0;
  square90.alpha_deg = 90.0;
  const auto r1 = evaluate_rib(set, square90, 30000.0);
  EXPECT_NEAR(r1.R, 3.5726760161, 1e-6);
  EXPECT_NEAR(r1.f, 0.0508531225, 1e-6);
  EXPECT_NEAR(r1.e_plus, 298.9820302819, 1e-3);
  EXPECT_NEAR(r1.G, 17.1527507, 1e-4);  // matches only with the e+ floor
  EXPECT_NEAR(r1.St_r, 0.0080325541, 1e-6);

  RibGeometry angled;
  angled.e_D = 0.0625;
  angled.p_e = 15.0;
  angled.W_H = 2.0;
  angled.alpha_deg = 45.0;
  const auto r2 = evaluate_rib(set, angled, 30000.0);
  EXPECT_NEAR(r2.R, 4.7592382829, 1e-6);
  EXPECT_NEAR(r2.f, 0.0440439106, 1e-6);

  // W/H = 3, above Eq. 4.17's own cap: the R quadratic must use W/H = 2, not
  // 3, in its (W/H)^m term. If the cap were dropped this would not match.
  RibGeometry capped;
  capped.e_D = 0.06;
  capped.p_e = 12.0;
  capped.W_H = 3.0;
  capped.alpha_deg = 60.0;
  const auto r3 = evaluate_rib(set, capped, 40000.0);
  EXPECT_NEAR(r3.R, 2.9903078484, 1e-6);
  EXPECT_NEAR(r3.f, 0.0876323069, 1e-6);
}

// Eq. 4.17's own cap ("if W/H > 2, set W/H = 2") must actually engage: R at
// W/H = 3 and W/H = 5 should be identical to R at the capped value 2, not
// merely close, since above the cap the (W/H)^m term stops depending on
// W/H at all.
TEST(RibCorrelationTest, HanParkWHCapActuallyCaps) {
  const auto set = han_park_1988_angled();
  RibGeometry g;
  g.e_D = 0.06;
  g.p_e = 12.0;
  g.alpha_deg = 60.0;  // off 90 deg, so the (W/H)^m term is active (m=0.35)

  g.W_H = 2.0;
  const double R_at_cap = evaluate_rib(set, g, 40000.0).R;
  g.W_H = 3.0;
  const double R_above = evaluate_rib(set, g, 40000.0).R;
  g.W_H = 5.0;
  const double R_further_above = evaluate_rib(set, g, 40000.0).R;

  EXPECT_NEAR(R_above, R_at_cap, 1e-9);
  EXPECT_NEAR(R_further_above, R_at_cap, 1e-9);

  // Below the cap, R must actually vary with W/H -- otherwise the cap logic
  // could be silently always-on and this test would not distinguish it.
  g.W_H = 1.0;
  const double R_below = evaluate_rib(set, g, 40000.0).R;
  EXPECT_GT(std::abs(R_below - R_at_cap), 1e-6);
}

// The R and G switches are genuine discontinuities the source states, not
// numerical artefacts to be smoothed away -- pin the jump sizes so a future
// "helpful" smoothing change cannot land silently.
TEST(RibCorrelationTest, HanParkRAndGSwitchesAreGenuineDiscontinuities) {
  const auto set = han_park_1988_angled();
  RibGeometry g;
  g.e_D = 0.06;
  g.p_e = 12.0;
  g.W_H = 2.0;

  g.alpha_deg = 90.0;
  const double R_at_90 = evaluate_rib(set, g, 40000.0).R;
  g.alpha_deg = 90.0 - 1e-6;
  const double R_just_below_90 = evaluate_rib(set, g, 40000.0).R;
  // A ~28% jump at W/H=2 (2^0.35 - 1), not a rounding-sized difference.
  EXPECT_GT(std::abs(R_just_below_90 / R_at_90 - 1.0), 0.25);

  RibGeometry g2;
  g2.e_D = 0.06;
  g2.p_e = 20.0;
  g2.alpha_deg = 30.0;

  g2.W_H = 1.0;
  const double G_at_square = evaluate_rib(set, g2, 40000.0).G;
  g2.W_H = 1.0 - 1e-6;
  const double G_just_off_square = evaluate_rib(set, g2, 40000.0).G;
  EXPECT_GT(std::abs(G_just_off_square / G_at_square - 1.0), 0.20);
}

// Same reverse-flow symmetry as the other two sets: a parallel angled rib
// reversed is its mirror image.
TEST(RibCorrelationTest, HanParkReverseFlowIsSymmetricInMagnitude) {
  const auto set = han_park_1988_angled();
  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 15.0;
  g.W_H = 2.0;
  g.alpha_deg = 45.0;
  const auto fwd = evaluate_rib(set, g, 30000.0);
  const auto rev = evaluate_rib(set, g, -30000.0);

  EXPECT_NEAR(fwd.G, rev.G, 1e-9);
  EXPECT_NEAR(fwd.St_r, rev.St_r, 1e-9);
  EXPECT_NEAR(fwd.e_plus, -rev.e_plus, 1e-6);
}

// Same guard discipline: every state a solver can probe must return finite,
// positive numbers, across the full alpha and W/H range including the two
// switch boundaries themselves.
TEST(RibCorrelationTest, HanParkGuardsHoldForStatesTheSolverActuallyProbes) {
  const auto set = han_park_1988_angled();
  for (double alpha : {30.0, 45.0, 60.0, 89.999, 90.0, 90.001}) {
    for (double W_H : {0.9999, 1.0, 1.0001, 2.0, 4.0, 10.0}) {
      RibGeometry g;
      g.e_D = 0.0625;
      g.p_e = 15.0;
      g.W_H = W_H;
      g.alpha_deg = alpha;
      for (double Re : {-1e6, -1.0, 0.0, 1.0, 1e8}) {
        const auto r = evaluate_rib(set, g, Re);
        EXPECT_TRUE(std::isfinite(r.f))
            << "alpha=" << alpha << " W_H=" << W_H << " Re=" << Re;
        EXPECT_TRUE(std::isfinite(r.G))
            << "alpha=" << alpha << " W_H=" << W_H << " Re=" << Re;
        EXPECT_TRUE(std::isfinite(r.St_r))
            << "alpha=" << alpha << " W_H=" << W_H << " Re=" << Re;
        EXPECT_GT(r.f, 0.0) << "alpha=" << alpha << " W_H=" << W_H;
        EXPECT_GT(r.St_r, 0.0) << "alpha=" << alpha << " W_H=" << W_H;
      }
    }
  }
}

// The analytic Stanton derivative must match a central difference through
// this set's guards too. Perturbing alpha by h risks crossing the alpha=90
// switch itself if alpha is near it, which would corrupt a naive central
// difference through no fault of the Jacobian -- so this sweeps Re (which
// the switch does not depend on) at a fixed, safely-off-switch geometry.
TEST(RibCorrelationTest, HanParkStantonDerivativeMatchesCentralDifferences) {
  const auto set = han_park_1988_angled();
  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 15.0;
  g.W_H = 2.0;
  g.alpha_deg = 45.0;
  for (double Re : {50000.0, 20000.0, 100.0, 10.0}) {
    const double h = std::max(1e-6, std::abs(Re) * 1e-6);
    const double fd = (evaluate_rib(set, g, Re + h).St_r -
                       evaluate_rib(set, g, Re - h).St_r) /
                      (2.0 * h);
    const double analytic = evaluate_rib(set, g, Re).dSt_dRe;
    const double scale = std::max({std::abs(fd), std::abs(analytic), 1e-14});
    EXPECT_LT(std::abs(analytic - fd) / scale, 1e-5) << "Re = " << Re;
  }
}

TEST(RibCorrelationTest, HanParkSetPassesValidation) {
  EXPECT_NO_THROW(validate_rib_set(han_park_1988_angled()));

  // C_R plays no role in the quadratic-alpha shape, so it must NOT be
  // rejected for being unset (0.0) the way the plain power-law shape is.
  auto s = han_park_1988_angled();
  EXPECT_EQ(s.C_R, 0.0);
  EXPECT_NO_THROW(validate_rib_set(s));

  auto bad_cap = han_park_1988_angled();
  bad_cap.R_quad_WH_cap = -1.0;
  EXPECT_THROW(validate_rib_set(bad_cap), std::invalid_argument);
}

// Han (1988) (90 deg only) and Han & Park (1988) (angled, this set) both
// claim alpha = 90 deg, and their R correlations were confirmed to agree
// closely there (han_ribbed.md item 23: 3.44%) while their G correlations
// disagree by up to 22% -- different papers, not obliged to agree, and the
// asymmetry between the two functions is itself part of the confirmed
// record. Pinned so a future edit cannot silently "fix" G to match R's
// agreement, which would misrepresent what the sources actually say.
TEST(RibCorrelationTest, HanParkAgreesWithHan1988OnRButNotG) {
  const auto han = han_1988_orthogonal();
  const auto park = han_park_1988_angled();

  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 10.0;
  g.W_H = 1.0;
  g.alpha_deg = 90.0;

  const auto r_han = evaluate_rib(han, g, 30000.0);
  const auto r_park = evaluate_rib(park, g, 30000.0);
  const double R_ratio = r_park.R / r_han.R;
  EXPECT_NEAR(R_ratio, 1.0, 0.06);  // within ~6%, close per item 23

  const double G_ratio = r_park.G / r_han.G;
  EXPECT_GT(std::abs(G_ratio - 1.0), 0.05);  // genuinely disagree
}

// The normaliser is data, not a convention: the same correlation written
// against a raw p/e needs a different constant, and mixing them is wrong by a
// constant factor that never looks like a trend.
TEST(RibCorrelationTest, NormaliserAndConstantAreEquivalentWhenConsistent) {
  auto normalised = han_1988_orthogonal();  // 3.2 * (p/e / 10)^0.35

  auto raw = han_1988_orthogonal();
  raw.C_R = 3.2 * std::pow(10.0, -0.35);  // 1.4294...
  raw.R_pe.reference = 1.0;

  auto g = ref_geometry();
  for (double p_e : {5.0, 10.0, 20.0}) {
    g.p_e = p_e;
    EXPECT_NEAR(evaluate_rib(normalised, g, 30000.0).R,
                evaluate_rib(raw, g, 30000.0).R, 1e-12)
        << "p/e = " << p_e;
  }

  // And the mistake this prevents: the normalised constant used raw.
  auto mixed = han_1988_orthogonal();
  mixed.R_pe.reference = 1.0;  // constant left at 3.2
  g.p_e = 10.0;
  EXPECT_NEAR(evaluate_rib(mixed, g, 30000.0).R /
                  evaluate_rib(normalised, g, 30000.0).R,
              std::pow(10.0, 0.35), 1e-9);
}
