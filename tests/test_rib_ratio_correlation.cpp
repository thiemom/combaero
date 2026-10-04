#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "heat_transfer.h"
#include "rib_ratio_correlation.h"

using combaero::cooling::evaluate_rib_ratio;
using combaero::cooling::RatioBelowFloor;
using combaero::cooling::RibGeometry;
using combaero::cooling::RibRatioOptions;
using combaero::cooling::RibRatioSet;
using combaero::cooling::validate_rib_ratio_options;
using combaero::cooling::validate_rib_ratio_set;

namespace {

// A user-style set: ratio Nu/Nu_DB = 2.5 (Re/1e4)^-0.2, f/f0 = 8 (Re/1e4)^0.1,
// Dittus-Boelter and Blasius-form baselines, fitted down to Re 1e4. The
// NEGATIVE Re exponent is deliberate: held naively it would blow up at Re->0.
RibRatioSet user_set() {
  RibRatioSet s;
  s.name = "test_ratio";
  s.source = "unit test";
  s.C_Nu = 2.5;
  s.Nu_Re = {-0.2, 1.0e4};
  s.C_f = 8.0;
  s.f_Re = {0.1, 1.0e4};
  s.Nu0_source = {0.023, 0.8, 0.4};
  s.f0_source = {0.046, -0.2, 0.0};
  s.Re_floor = 1.0e4;
  s.valid_Re = {1.0e4, 9.0e4};
  return s;
}

RibGeometry geom() {
  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 10.0;
  g.W_H = 1.0;
  g.alpha_deg = 90.0;
  return g;
}

}  // namespace

// Inside the fitted range the set reproduces the paper -- r * Nu0_source --
// whatever extrapolation baseline the caller picks.
TEST(RibRatioTest, InRangeIsTheSourceProductForEveryBelowFloorChoice) {
  const auto s = user_set();
  for (auto choice : {RatioBelowFloor::Gnielinski, RatioBelowFloor::SourceBaseline,
                      RatioBelowFloor::User}) {
    RibRatioOptions o;
    o.below_floor = choice;
    o.user_Nu0 = {0.02, 0.8, 0.33};
    for (double Re : {1.0e4, 3.0e4, 9.0e4}) {
      const auto r = evaluate_rib_ratio(s, geom(), Re, 0.7, o);
      const double x = std::sqrt(Re * Re + 1.0);
      const double expect = 2.5 * std::pow(x / 1e4, -0.2) * 0.023 *
                            std::pow(x, 0.8) * std::pow(0.7, 0.4);
      EXPECT_NEAR(r.Nu, expect, 1e-7 * expect) << "Re=" << Re;
      EXPECT_FALSE(r.below_floor);
    }
  }
}

// Far below the floor, Nu is the extrapolation baseline scaled to match at
// the floor, and with Gnielinski it inherits the laminar limit at Re -> 0
// instead of Dittus-Boelter's collapse to zero.
TEST(RibRatioTest, BelowTheFloorHandsOverToGnielinski) {
  const auto s = user_set();
  const auto at_floor = evaluate_rib_ratio(s, geom(), 1.0e4, 0.7);
  const double k = at_floor.Nu / combaero::nusselt_gnielinski_smooth_with_derivative(1.0e4, 0.7).Nu;
  for (double Re : {3000.0, 500.0, 1.0, 0.0}) {
    const auto r = evaluate_rib_ratio(s, geom(), Re, 0.7);
    const double x = std::sqrt(Re * Re + 1.0);
    EXPECT_NEAR(r.Nu, k * combaero::nusselt_gnielinski_smooth_with_derivative(x, 0.7).Nu, 1e-8 * r.Nu)
        << "Re=" << Re;
    EXPECT_TRUE(r.below_floor && r.extrapolated);
  }
  const double laminar = evaluate_rib_ratio(s, geom(), 0.0, 0.7).Nu;
  EXPECT_NEAR(laminar, k * combaero::NU_LAMINAR_CONST_T, 1e-8 * laminar);

  // Holding the ratio on a Dittus-Boelter source instead collapses toward 0:
  // at Re = 0 (smoothed |Re| = 1) it is under 1% of the laminar handover.
  RibRatioOptions src;
  src.below_floor = RatioBelowFloor::SourceBaseline;
  EXPECT_LT(evaluate_rib_ratio(s, geom(), 0.0, 0.7, src).Nu, 0.01 * laminar);
}

// Every state the solver can reach: finite, non-negative, even in Re with an
// odd derivative -- including the negative-exponent ratio at Re = 0.
TEST(RibRatioTest, GuardsHoldAcrossZeroAndBothSigns) {
  const auto s = user_set();
  for (double Re : {-1e7, -2e4, -1e4, -7e3, -1.0, 0.0, 1.0, 7e3, 1e4, 2e4, 1e7}) {
    for (double Pr : {-1.0, 0.0, 0.7, 50.0}) {
      const auto p = evaluate_rib_ratio(s, geom(), Re, Pr);
      const auto m = evaluate_rib_ratio(s, geom(), -Re, Pr);
      EXPECT_TRUE(std::isfinite(p.Nu) && std::isfinite(p.f) &&
                  std::isfinite(p.dNu_dRe) && std::isfinite(p.df_dRe))
          << "Re=" << Re << " Pr=" << Pr;
      EXPECT_GE(p.Nu, 0.0);
      EXPECT_GE(p.f, 0.0);
      EXPECT_DOUBLE_EQ(p.Nu, m.Nu);
      EXPECT_DOUBLE_EQ(p.f, m.f);
      EXPECT_DOUBLE_EQ(p.dNu_dRe, -m.dNu_dRe);
      EXPECT_DOUBLE_EQ(p.df_dRe, -m.df_dRe);
    }
  }
}

// Analytic derivatives against central differences in every regime: in
// range, inside the blend band, far below it, near and at Re = 0 -- for all
// three below-floor choices -- including Re 3000, where smooth-pipe
// Gnielinski's friction clamp had a slope jump until #446 made it C1.
TEST(RibRatioTest, DerivativesMatchCentralDifferencesEverywhere) {
  const auto s = user_set();
  for (auto choice : {RatioBelowFloor::Gnielinski, RatioBelowFloor::SourceBaseline,
                      RatioBelowFloor::User}) {
    RibRatioOptions o;
    o.below_floor = choice;
    o.user_Nu0 = {0.02, 0.8, 0.33};
    for (double Re : {5.0e4, 1.2e4, 8.0e3, 6.0e3, 5.2e3, 3.5e3, 3.0e3, 2000.0, 0.5, -3.0e3}) {
      const double h = std::max(1e-4, std::abs(Re) * 1e-6);
      const auto c = evaluate_rib_ratio(s, geom(), Re, 0.7, o);
      const auto p = evaluate_rib_ratio(s, geom(), Re + h, 0.7, o);
      const auto m = evaluate_rib_ratio(s, geom(), Re - h, 0.7, o);
      const double fdN = (p.Nu - m.Nu) / (2.0 * h);
      const double fdF = (p.f - m.f) / (2.0 * h);
      EXPECT_LT(std::abs(c.dNu_dRe - fdN) / std::max({std::abs(fdN), 1e-9}), 1e-4)
          << "Nu, Re=" << Re;
      EXPECT_LT(std::abs(c.df_dRe - fdF) / std::max({std::abs(fdF), 1e-12}), 1e-4)
          << "f, Re=" << Re;
    }
  }
}

// C1 across the blend band ends: the value and slope just inside and just
// outside each end agree.
TEST(RibRatioTest, HandoverIsC1AtBothEndsOfTheBlend) {
  const auto s = user_set();
  for (double edge : {1.0e4, 5.0e3}) {
    const auto lo = evaluate_rib_ratio(s, geom(), edge * (1.0 - 1e-7), 0.7);
    const auto hi = evaluate_rib_ratio(s, geom(), edge * (1.0 + 1e-7), 0.7);
    EXPECT_NEAR(lo.Nu, hi.Nu, 1e-6 * hi.Nu) << edge;
    EXPECT_NEAR(lo.dNu_dRe, hi.dNu_dRe, 1e-4 * std::abs(hi.dNu_dRe)) << edge;
    EXPECT_NEAR(lo.f, hi.f, 1e-6 * hi.f) << edge;
    EXPECT_NEAR(lo.df_dRe, hi.df_dRe, 1e-4 * std::abs(hi.df_dRe) + 1e-15) << edge;
  }
}

TEST(RibRatioTest, MistakesAreRejectedOnce) {
  EXPECT_NO_THROW(validate_rib_ratio_set(user_set()));
  auto a = user_set();
  a.C_Nu = 0.0;
  EXPECT_THROW(validate_rib_ratio_set(a), std::invalid_argument);
  auto b = user_set();
  b.Nu0_source = {};
  EXPECT_THROW(validate_rib_ratio_set(b), std::invalid_argument);
  auto c = user_set();
  c.Re_floor = -1.0;
  EXPECT_THROW(validate_rib_ratio_set(c), std::invalid_argument);
  auto d = user_set();
  d.f_Re = {0.1, 0.0};
  EXPECT_THROW(validate_rib_ratio_set(d), std::invalid_argument);

  RibRatioOptions o;
  o.below_floor = RatioBelowFloor::User;
  EXPECT_THROW(validate_rib_ratio_options(o), std::invalid_argument);
  // ...and an evaluation with those unvalidated options still stays finite.
  const auto r = evaluate_rib_ratio(user_set(), geom(), 100.0, 0.7, o);
  EXPECT_TRUE(std::isfinite(r.Nu) && r.Nu >= 0.0);
}

// Taslim & Spring (1987): the shipped ratio sets reproduce their own inputs
// -- C at Re_ref, the Re-independent f across the friction range -- and
// refuse configurations the paper did not test.
TEST(RibRatioTest, TaslimSpring1987SetsCarryThePapersNumbers) {
  using combaero::cooling::taslim_spring_1987;
  struct Cfg {
    double ar, eD, C, f;
  };
  const Cfg cfgs[] = {{0.5, 0.125, 3.34316, 0.10146}, {0.5, 0.250, 3.90148, 0.51645},
                      {1.0, 0.083, 3.42519, 0.04680}, {1.0, 0.167, 3.99269, 0.11610},
                      {3.5, 0.053, 3.03740, 0.01498}, {3.5, 0.107, 3.12189, 0.02090},
                      {3.5, 0.161, 3.33653, 0.03162}};
  for (const auto &c : cfgs) {
    const auto s = taslim_spring_1987(c.ar, c.eD);
    EXPECT_NO_THROW(validate_rib_ratio_set(s)) << s.name;
    EXPECT_NEAR(s.valid_WH.lo, 1.0 / c.ar, 1e-12) << "Han W/H = 1 / Taslim AR";
    RibGeometry g;
    g.e_D = c.eD;
    g.p_e = 10.0;
    g.W_H = 1.0 / c.ar;
    g.alpha_deg = 90.0;
    // At Re = 1e4 * k with k chosen inside the fitted range, Nu/Nu_DB is
    // C (Re/1e4)^-0.2 exactly; f is the configuration's mean.
    const double Re = std::max(s.Re_floor, s.Re_floor_f) * 1.5;
    const auto r = evaluate_rib_ratio(s, g, Re, 0.70);
    const double x = std::sqrt(Re * Re + 1.0);
    const double nu_db = 0.023 * std::pow(x, 0.8) * std::pow(0.70, 0.4);
    EXPECT_NEAR(r.Nu, c.C * std::pow(x / 1e4, -0.2) * nu_db, 1e-6 * r.Nu) << s.name;
    EXPECT_NEAR(r.f, c.f, 1e-6 * c.f) << s.name;
    EXPECT_FALSE(r.extrapolated) << s.name;
  }
  EXPECT_THROW(taslim_spring_1987(1.0, 0.250), std::invalid_argument)
      << "AR 1.0 e/D 0.25 has friction but no Nu in the paper";
  EXPECT_THROW(taslim_spring_1987(2.0, 0.125), std::invalid_argument);
}
