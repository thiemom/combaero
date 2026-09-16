#include <gtest/gtest.h>

#include <cmath>

#include "rib_correlation.h"

using combaero::cooling::evaluate_rib;
using combaero::cooling::han_1988_orthogonal;
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
