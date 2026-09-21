#include <gtest/gtest.h>

#include <cmath>

#include "impingement_correlation.h"

using combaero::cooling::crossflow_to_jet_ratio_at_row;
using combaero::cooling::crossflow_to_jet_ratio_at_x;
using combaero::cooling::florschuetz_1981_inline;
using combaero::cooling::florschuetz_1981_staggered;
using combaero::cooling::FLORSCHUETZ_1981_DEFAULT_CD;
using combaero::cooling::goldstein_1986_single_jet;
using combaero::cooling::ImpingementThermalBC;
using combaero::cooling::jet_array_impingement_nu;
using combaero::cooling::JetArrayCorrelationSet;
using combaero::cooling::single_jet_impingement_nu;
using combaero::cooling::validate_jet_array_set;
using combaero::cooling::validate_single_jet_set;

// -----------------------------------------------------------------
// Single jet (Goldstein, Behbahani and Heppelmann, 1986)
// -----------------------------------------------------------------

// The numbers below are the closed-form check point recorded in the
// CONFIRMED extraction at validation/cooling/extractions/han_impingement.md,
// item 4. They come from the document, not from this implementation, so
// this test can fail.
TEST(SingleJetImpingementTest, ReproducesTheConfirmedCheckPoint) {
  const auto set = goldstein_1986_single_jet();

  // Han's book states these as rounded integers (60, 56); the tolerance
  // matches that rounding precision, not the implementation.
  const double Nu_const_q = single_jet_impingement_nu(
      set, ImpingementThermalBC::ConstantHeatFlux, 25000.0, 7.75, 5.0);
  EXPECT_NEAR(Nu_const_q, 60.0, 0.5);

  const double Nu_const_Tw = single_jet_impingement_nu(
      set, ImpingementThermalBC::ConstantWallTemperature, 25000.0, 7.75, 5.0);
  EXPECT_NEAR(Nu_const_Tw, 56.0, 0.5);
}

// Falsification: swapping which exponent (0.76 vs n) applies to Re vs R/D
// misses the check point by two orders of magnitude -- this was a real bug
// caught while writing this test, not a hypothetical one. Pin the direction
// explicitly so a regression is loud rather than silently wrong.
TEST(SingleJetImpingementTest, ReExponentIsFixedNotTheBoundaryConditionOne) {
  auto set = goldstein_1986_single_jet();
  const double Re_exponent = set.Re_exponent;
  const double n = set.n_const_heat_flux;
  ASSERT_NE(Re_exponent, n);

  set.Re_exponent = n;  // deliberately wrong
  set.n_const_heat_flux = Re_exponent;
  const double wrong = single_jet_impingement_nu(
      set, ImpingementThermalBC::ConstantHeatFlux, 25000.0, 7.75, 5.0);
  EXPECT_GT(std::abs(wrong - 60.0), 100.0);
}

// Eq. 4.1's optimum is structural, not a separately fitted fact: A itself IS
// the numerator's value at L/D=7.75, for any R/D or Re.
TEST(SingleJetImpingementTest, OptimumSpacingIsExactlyWhereTheAbsTermIsZero) {
  const auto set = goldstein_1986_single_jet();
  const double at_optimum = single_jet_impingement_nu(
      set, ImpingementThermalBC::ConstantHeatFlux, 10000.0, 7.75, 3.0);
  const double away_from_optimum = single_jet_impingement_nu(
      set, ImpingementThermalBC::ConstantHeatFlux, 10000.0, 4.0, 3.0);
  EXPECT_GT(at_optimum, away_from_optimum);
}

TEST(SingleJetImpingementTest, ValidateAcceptsThePublishedSet) {
  EXPECT_NO_THROW(validate_single_jet_set(goldstein_1986_single_jet()));
}

TEST(SingleJetImpingementTest, ValidateRejectsNonPositiveA) {
  auto set = goldstein_1986_single_jet();
  set.A = 0.0;
  EXPECT_THROW(validate_single_jet_set(set), std::invalid_argument);
}

TEST(SingleJetImpingementTest, StaysFiniteForReversedFlow) {
  const auto set = goldstein_1986_single_jet();
  const double Nu = single_jet_impingement_nu(
      set, ImpingementThermalBC::ConstantHeatFlux, -25000.0, 7.75, 5.0);
  EXPECT_TRUE(std::isfinite(Nu));
}

// -----------------------------------------------------------------
// Jet array with crossflow (Florschuetz, Truman and Metzger, 1981)
// -----------------------------------------------------------------

TEST(JetArrayCorrelationTest, TableCoefficientsMatchTheConfirmedExtraction) {
  const auto inl = florschuetz_1981_inline();
  EXPECT_DOUBLE_EQ(inl.A_fit.C, 1.18);
  EXPECT_DOUBLE_EQ(inl.A_fit.nx, -0.944);
  EXPECT_DOUBLE_EQ(inl.n_fit.nz, 1.04);

  const auto stg = florschuetz_1981_staggered();
  EXPECT_DOUBLE_EQ(stg.A_fit.C, 1.87);
  EXPECT_DOUBLE_EQ(stg.B_fit.nz, 0.059);
}

// Table 1's own worked geometry: (xn/d, yn/d, z/d) = (5, 4, 1), a real tested
// configuration from the primary paper. This pins the four power-law fits
// evaluating to a finite, positive Nu at a real operating point rather than
// merely "compiles".
TEST(JetArrayCorrelationTest, EvaluatesAtARealTestedGeometry) {
  const auto set = florschuetz_1981_inline();
  const auto result =
      jet_array_impingement_nu(set, 10000.0, 0.3, 0.7, 5.0, 4.0, 1.0);
  EXPECT_GT(result.Nu, 0.0);
  EXPECT_TRUE(std::isfinite(result.Nu));
  EXPECT_FALSE(result.extrapolated);
}

// Increasing crossflow at fixed Re_j and geometry must decrease Nu -- this is
// the whole physical point of the bracket term, stated throughout the source
// (jets degrade as spent air from upstream rows accumulates).
TEST(JetArrayCorrelationTest, MoreCrossflowDegradesNu) {
  const auto set = florschuetz_1981_inline();
  const auto low_crossflow =
      jet_array_impingement_nu(set, 10000.0, 0.1, 0.7, 5.0, 8.0, 1.0);
  const auto high_crossflow =
      jet_array_impingement_nu(set, 10000.0, 0.6, 0.7, 5.0, 8.0, 1.0);
  EXPECT_LT(high_crossflow.Nu, low_crossflow.Nu);
}

TEST(JetArrayCorrelationTest, FlagsExtrapolationOutsideValidity) {
  const auto set = florschuetz_1981_inline();
  const auto inside =
      jet_array_impingement_nu(set, 10000.0, 0.3, 0.7, 10.0, 6.0, 2.0);
  EXPECT_FALSE(inside.extrapolated);

  const auto outside_re =
      jet_array_impingement_nu(set, 1000.0, 0.3, 0.7, 10.0, 6.0, 2.0);
  EXPECT_TRUE(outside_re.extrapolated);
}

// Staggered's xn/d bound is genuinely tighter than inline's (5-10 vs 5-15) --
// this is a stated fact from the primary paper, not a smoothed-over detail.
TEST(JetArrayCorrelationTest, StaggeredHasATighterXnDBoundThanInline) {
  const auto inl = florschuetz_1981_inline();
  const auto stg = florschuetz_1981_staggered();
  EXPECT_DOUBLE_EQ(inl.valid_xn_d.hi, 15.0);
  EXPECT_DOUBLE_EQ(stg.valid_xn_d.hi, 10.0);
}

TEST(JetArrayCorrelationTest, ValidateAcceptsBothPublishedSets) {
  EXPECT_NO_THROW(validate_jet_array_set(florschuetz_1981_inline()));
  EXPECT_NO_THROW(validate_jet_array_set(florschuetz_1981_staggered()));
}

TEST(JetArrayCorrelationTest, ValidateRejectsNonFiniteFit) {
  auto set = florschuetz_1981_inline();
  set.m_fit.nx = std::nan("");
  EXPECT_THROW(validate_jet_array_set(set), std::invalid_argument);
}

TEST(JetArrayCorrelationTest, StaysFiniteForReversedFlowAndOutOfRangeGcGj) {
  const auto set = florschuetz_1981_inline();
  const auto reversed =
      jet_array_impingement_nu(set, -10000.0, 0.3, 0.7, 10.0, 6.0, 2.0);
  EXPECT_TRUE(std::isfinite(reversed.Nu));

  const auto negative_crossflow =
      jet_array_impingement_nu(set, 10000.0, -0.3, 0.7, 10.0, 6.0, 2.0);
  EXPECT_TRUE(std::isfinite(negative_crossflow.Nu));
}

// -----------------------------------------------------------------
// Crossflow-to-jet ratio (Florschuetz Eq. 8)
// -----------------------------------------------------------------

// Nu1's own definition (Florschuetz nomenclature, item 21): Gc/Gj = 0 at the
// FIRST spanwise row, for any geometry or discharge coefficient.
TEST(CrossflowRatioTest, IsExactlyZeroAtTheFirstRow) {
  EXPECT_DOUBLE_EQ(crossflow_to_jet_ratio_at_row(8.0, 2.0, 0.79, 1), 0.0);
  EXPECT_DOUBLE_EQ(crossflow_to_jet_ratio_at_row(4.0, 1.0, 0.65, 1), 0.0);
}

// The paper states the crossflow ratio increases monotonically upstream to
// downstream for the tested geometries -- verify row 10 exceeds row 1 (which
// is exactly zero) for a representative array.
TEST(CrossflowRatioTest, IncreasesFromUpstreamToDownstreamRow) {
  const double row1 =
      crossflow_to_jet_ratio_at_row(8.0, 2.0, FLORSCHUETZ_1981_DEFAULT_CD, 1);
  const double row10 = crossflow_to_jet_ratio_at_row(
      8.0, 2.0, FLORSCHUETZ_1981_DEFAULT_CD, 10);
  EXPECT_GT(row10, row1);
}

// Stated explicitly by the paper: "the flow distribution is independent of
// the streamwise hole spacing and hole pattern, depending... only on the
// geometric parameter (yn/d)(z/d)". crossflow_to_jet_ratio takes no xn_d
// argument at all -- this test pins that the SAME (yn/d)(z/d) product gives
// the same ratio even though xn/d is not passed, i.e. the physical
// invariant the source states is what the function's signature encodes.
TEST(CrossflowRatioTest, DependsOnlyOnYnDTimesZD) {
  // (yn/d, z/d) = (8, 2) has the same product as (4, 4) and (16, 1).
  const double a = crossflow_to_jet_ratio_at_row(8.0, 2.0, 0.79, 6);
  const double b = crossflow_to_jet_ratio_at_row(4.0, 4.0, 0.79, 6);
  const double c = crossflow_to_jet_ratio_at_row(16.0, 1.0, 0.79, 6);
  EXPECT_NEAR(a, b, 1e-12);
  EXPECT_NEAR(a, c, 1e-12);
}

TEST(CrossflowRatioTest, RowConvenienceMatchesTheContinuousForm) {
  const double from_row =
      crossflow_to_jet_ratio_at_row(8.0, 2.0, 0.79, 4);
  const double from_x = crossflow_to_jet_ratio_at_x(8.0, 2.0, 0.79, 3.5);
  EXPECT_DOUBLE_EQ(from_row, from_x);
}

TEST(CrossflowRatioTest, StaysFiniteForDegenerateGeometry) {
  EXPECT_TRUE(std::isfinite(crossflow_to_jet_ratio_at_row(0.0, 0.0, 0.79, 5)));
  EXPECT_TRUE(std::isfinite(crossflow_to_jet_ratio_at_row(8.0, 2.0, 0.0, 5)));
}
