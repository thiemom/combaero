#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

#include "math_constants.h"  // MSVC compatibility for M_PI
#include "pin_fin_correlation.h"

using namespace combaero::cooling;

namespace {

PinFinGeometry staggered(double S, double X, double H, int N) {
  PinFinGeometry g;
  g.S_D = S;
  g.X_D = X;
  g.H_D = H;
  g.N_rows = N;
  g.arrangement = PinArrangement::Staggered;
  return g;
}

double rel(double a, double b) { return std::abs(a / b - 1.0); }

}  // namespace

// Each shipped set reproduces its printed equation in its own basis.
TEST(PinFinTest, MetzgerReproducesArmstrongWinstanleyEq2) {
  const auto set = metzger_1986_staggered_nu();
  validate_pin_fin_nu_set(set);
  const auto g = staggered(2.5, 2.0, 1.0, 10);
  for (double Re : {3.0e3, 4.0e4}) {
    const auto r = evaluate_pin_fin_nu(set, g, Re, 0.71);
    const double x = std::sqrt(Re * Re + 1.0);
    EXPECT_LT(rel(r.Nu, 0.135 * std::pow(x, 0.69) * std::pow(2.0, -0.34)),
              1e-12);
    EXPECT_FALSE(r.extrapolated);
  }
}

TEST(PinFinTest, MetzgerFrictionIsExactOutsideTheBlendAndC1Across) {
  const auto set = metzger_1982_staggered_friction();
  validate_pin_fin_friction_set(set);
  const auto g = staggered(2.5, 2.5, 1.0, 10);
  const double lo = 1.0e4 / PIN_FIN_F_BLEND;
  const double hi = 1.0e4 * PIN_FIN_F_BLEND;
  for (double Re : {2.0e3, lo * 0.999}) {
    const double x = std::sqrt(Re * Re + 1.0);
    EXPECT_LT(rel(evaluate_pin_fin_friction(set, g, Re).f,
                  0.317 * std::pow(x, -0.132)),
              1e-12);
  }
  for (double Re : {hi * 1.001, 8.0e4}) {
    const double x = std::sqrt(Re * Re + 1.0);
    EXPECT_LT(rel(evaluate_pin_fin_friction(set, g, Re).f,
                  1.76 * std::pow(x, -0.318)),
              1e-12);
  }
  // Scan the band: no step and no slope jump.
  double prev = evaluate_pin_fin_friction(set, g, lo * 0.98).f;
  double prev_slope = evaluate_pin_fin_friction(set, g, lo * 0.98).df_dRe;
  for (double Re = lo * 0.98; Re <= hi * 1.02; Re *= 1.005) {
    const auto r = evaluate_pin_fin_friction(set, g, Re);
    EXPECT_LT(rel(r.f, prev), 0.002) << Re;
    // C1: the slope moves from -0.132 f/Re to -0.318 f/Re across the band
    // without a jump; each 0.5% step changes it by well under 10%.
    EXPECT_LT(std::abs(r.df_dRe - prev_slope),
              0.10 * std::abs(prev_slope)) << Re;
    prev = r.f;
    prev_slope = r.df_dRe;
  }
}

TEST(PinFinTest, VanFossenRoundTripsThroughDprime) {
  const auto set = vanfossen_1982_staggered_nu();
  validate_pin_fin_nu_set(set);
  // VanFossen's H/D 2 array: S/D 4, X/D = 4 sqrt(3)/2.
  const auto g = staggered(4.0, 2.0 * std::sqrt(3.0), 2.0, 4);
  EXPECT_NEAR(pin_fin_dprime_over_D(g), 3.225, 0.001);
  const double dp = pin_fin_dprime_over_D(g);
  const double ap = pin_fin_aprime_over_amin(g);
  const double Re_D = 1.0e4;
  const auto r = evaluate_pin_fin_nu(set, g, Re_D, 0.71);
  const double Re_p = std::sqrt(Re_D * Re_D + 1.0) * dp / ap;
  EXPECT_NEAR(r.Re_native, Re_D * dp / ap, 1e-9 * Re_D);
  EXPECT_LT(rel(r.Nu * dp, 0.153 * std::pow(Re_p, 0.685)), 1e-12);
  EXPECT_FALSE(r.extrapolated);
}

TEST(PinFinTest, DamerowConvertsFromNMinusOneRows) {
  const auto set = damerow_1972_staggered_friction();
  validate_pin_fin_friction_set(set);
  const auto g = staggered(4.24, 2.12, 2.0, 10);
  const double Re = 1.0e4;
  const double x = std::sqrt(Re * Re + 1.0);
  const double native = 2.06 * std::pow(4.24, -1.1) * std::pow(x, -0.16);
  EXPECT_LT(rel(evaluate_pin_fin_friction(set, g, Re).f, native * 9.0 / 10.0),
            1e-12);
  EXPECT_FALSE(evaluate_pin_fin_friction(set, g, Re).extrapolated);
}

TEST(PinFinTest, ChyuSetsAndTheirRatioModifierAgreeExactly) {
  const auto stag = chyu_1998_nu(PinArrangement::Staggered, PinNuSurface::Total);
  const auto inl = chyu_1998_nu(PinArrangement::Inline, PinNuSurface::Total);
  const auto mod = chyu_1998_inline_over_staggered();
  validate_pin_fin_nu_set(stag);
  validate_pin_fin_nu_set(inl);
  validate_pin_fin_modifier(mod);
  auto gs = staggered(2.5, 2.5, 1.0, 7);
  auto gi = gs;
  gi.arrangement = PinArrangement::Inline;
  for (double Re : {7.0e3, 2.0e4}) {
    const double ns = evaluate_pin_fin_nu(stag, gs, Re, 0.71).Nu;
    const double ni = evaluate_pin_fin_nu(inl, gi, Re, 0.71).Nu;
    const auto m = evaluate_pin_fin_modifier(mod, gi, Re);
    EXPECT_LT(rel(ns * m.ratio_Nu, ni), 1e-12);
    EXPECT_FALSE(m.has_f);
    EXPECT_EQ(m.ratio_f, 1.0);
    EXPECT_FALSE(m.extrapolated);
  }
  // Table 4.7's printed form, Nu / Pr^0.4 = a Re^b.
  const double Re = 1.0e4;
  EXPECT_LT(rel(evaluate_pin_fin_nu(stag, gs, Re, 0.7).Nu,
                0.320 * std::pow(std::sqrt(Re * Re + 1.0), 0.583) *
                    std::pow(0.7, 0.4)),
            1e-6);
  EXPECT_THROW(chyu_1998_nu(static_cast<PinArrangement>(7), PinNuSurface::Pin),
               std::invalid_argument);
}

// Derivatives are analytic and match central differences, through Re = 0.
TEST(PinFinTest, DerivativesMatchCentralDifferences) {
  const auto nu_sets = {metzger_1986_staggered_nu(),
                        vanfossen_1982_staggered_nu()};
  const auto f_sets = {metzger_1982_staggered_friction(),
                       damerow_1972_staggered_friction()};
  const auto g = staggered(3.0, 2.5, 1.5, 10);
  for (double Re : {-2.0e4, -5.0, 0.0, 3.0, 4.0e3, 1.0e4, 3.0e4}) {
    const double h = std::max(1e-6, 1e-6 * std::abs(Re));
    for (const auto &s : nu_sets) {
      const double fd = (evaluate_pin_fin_nu(s, g, Re + h, 0.7).Nu -
                         evaluate_pin_fin_nu(s, g, Re - h, 0.7).Nu) /
                        (2 * h);
      const double an = evaluate_pin_fin_nu(s, g, Re, 0.7).dNu_dRe;
      EXPECT_NEAR(an, fd, 1e-5 * std::max(1.0, std::abs(fd))) << s.name << Re;
    }
    for (const auto &s : f_sets) {
      const double fd = (evaluate_pin_fin_friction(s, g, Re + h).f -
                         evaluate_pin_fin_friction(s, g, Re - h).f) /
                        (2 * h);
      const double an = evaluate_pin_fin_friction(s, g, Re).df_dRe;
      EXPECT_NEAR(an, fd, 1e-5 * std::max(1e-3, std::abs(fd))) << s.name << Re;
      EXPECT_TRUE(std::isfinite(evaluate_pin_fin_friction(s, g, Re).f));
    }
  }
}

TEST(PinFinTest, MinimumAreaSwitchesToTheDiagonalForDenseRows) {
  // Transverse-limited at the sources' geometries.
  EXPECT_DOUBLE_EQ(pin_fin_min_gap_D(staggered(2.5, 1.5, 1.0, 10)), 1.5);
  // S 4, X 1: diagonal 2 (sqrt(4 + 1) - 1) = 2.472 < 3.
  const auto dense = staggered(4.0, 1.0, 1.0, 10);
  EXPECT_NEAR(pin_fin_min_gap_D(dense), 2.0 * (std::sqrt(5.0) - 1.0), 1e-12);
  auto inl = dense;
  inl.arrangement = PinArrangement::Inline;
  EXPECT_DOUBLE_EQ(pin_fin_min_gap_D(inl), 3.0);
}

TEST(PinFinTest, AreaFractionsCloseTheUnitCell) {
  const auto g = staggered(2.5, 2.5, 1.0, 10);
  const auto a = pin_fin_area_fractions(g);
  EXPECT_NEAR(a.endwall_exposed, (6.25 - M_PI / 4) / 6.25, 1e-14);
  EXPECT_NEAR(a.pin, (M_PI / 2) / 6.25, 1e-14);
  EXPECT_NEAR(a.pin_over_total, a.pin / (a.pin + a.endwall_exposed), 1e-14);
}

TEST(PinFinTest, FinEfficiencyLimitsAndDerivative) {
  const double D = 1e-3, H = 1e-3, Af = 0.3;
  // Infinite conductivity and zero h: a perfect fin.
  EXPECT_NEAR(pin_fin_array_efficiency(500.0, 1e12, D, H, Af).eta_t, 1.0,
              1e-12);
  EXPECT_NEAR(pin_fin_array_efficiency(0.0, 20.0, D, H, Af).eta_t, 1.0,
              1e-15);
  // h -> 0 slope limit: -H^2 / (3 k D).
  EXPECT_NEAR(pin_fin_array_efficiency(0.0, 20.0, D, H, Af).deta_fin_dh,
              -H * H / (3.0 * 20.0 * D), 1e-12);
  // Against central differences across the series/closed-form switch.
  // Spans both sides of the series switch (z = 0.02 at h = 8 here).
  for (double h : {1e-3, 0.5, 7.9, 8.1, 2.0e3, 2.0e4, 1.0e5}) {
    const double dh = 1e-3 * h;
    const double fd =
        (pin_fin_array_efficiency(h + dh, 20.0, D, H, Af).eta_t -
         pin_fin_array_efficiency(h - dh, 20.0, D, H, Af).eta_t) /
        (2 * dh);
    EXPECT_NEAR(pin_fin_array_efficiency(h, 20.0, D, H, Af).deta_t_dh, fd,
                1e-6 * std::abs(fd))
        << h;
  }
  // A real case: Inconel pins in a trailing edge are noticeably below 1.
  const auto e = pin_fin_array_efficiency(2.0e3, 20.0, 1e-3, 2e-3, Af);
  EXPECT_LT(e.eta_fin, 0.9);
  EXPECT_GT(e.eta_fin, 0.5);
}

TEST(PinFinTest, OutOfBoxIsFlaggedNotRefused) {
  const auto set = metzger_1986_staggered_nu();
  EXPECT_FALSE(evaluate_pin_fin_nu(set, staggered(2.5, 2.5, 1.0, 10), 1e4, 0.7)
                   .extrapolated);
  EXPECT_TRUE(evaluate_pin_fin_nu(set, staggered(2.5, 2.5, 4.0, 10), 1e4, 0.7)
                  .extrapolated);  // H/D > 3
  EXPECT_TRUE(evaluate_pin_fin_nu(set, staggered(2.5, 2.5, 1.0, 4), 1e4, 0.7)
                  .extrapolated);  // row count
  EXPECT_TRUE(evaluate_pin_fin_nu(set, staggered(2.5, 2.5, 1.0, 10), 500, 0.7)
                  .extrapolated);  // Re
  auto inl = staggered(2.5, 2.5, 1.0, 10);
  inl.arrangement = PinArrangement::Inline;
  EXPECT_TRUE(evaluate_pin_fin_nu(set, inl, 1e4, 0.7).extrapolated);
  // Chyu's single geometry: S/D 3 is outside it.
  const auto chyu = chyu_1998_nu(PinArrangement::Staggered, PinNuSurface::Total);
  EXPECT_TRUE(
      evaluate_pin_fin_nu(chyu, staggered(3.0, 2.5, 1.0, 7), 1e4, 0.7)
          .extrapolated);
}

TEST(PinFinTest, ValidatorsRejectMistakes) {
  auto s = metzger_1986_staggered_nu();
  s.C = 0.0;
  EXPECT_THROW(validate_pin_fin_nu_set(s), std::invalid_argument);
  auto f = metzger_1982_staggered_friction();
  f.C2 = -1.0;
  EXPECT_THROW(validate_pin_fin_friction_set(f), std::invalid_argument);
  f = metzger_1982_staggered_friction();
  f.accuracy_f = StatedAccuracy::unstated();
  f.accuracy_f.provenance = AccuracyProvenance::Stated;  // claim without value
  EXPECT_THROW(validate_pin_fin_friction_set(f), std::invalid_argument);
  EXPECT_THROW(validate_pin_fin_geometry(staggered(1.0, 2.0, 1.0, 10)),
               std::invalid_argument);
  // S 1.5, X 0.5: sqrt(0.75^2 + 0.5^2) < 1, adjacent rows' pins overlap.
  EXPECT_THROW(validate_pin_fin_geometry(staggered(1.5, 0.5, 1.0, 10)),
               std::invalid_argument);
  EXPECT_NO_THROW(validate_pin_fin_geometry(staggered(2.5, 2.5, 1.0, 10)));
}

// Chyu (1990): Table 2 reproduced, Fig. 6 fits converted by 1/4, and the
// fillet modifier is the exact ratio of the two Table 2 sets.
TEST(PinFinTest, Chyu1990TableTwoAndFigureSixFits) {
  const auto g_stag = staggered(2.5, 2.5, 1.0, 7);
  auto g_inl = g_stag;
  g_inl.arrangement = PinArrangement::Inline;
  const double Re = 1.5e4;
  const double x = std::sqrt(Re * Re + 1.0);
  struct Row {
    PinArrangement a;
    bool fillet;
    double A, B;
  };
  for (const Row &r : {Row{PinArrangement::Inline, false, 0.463, 0.537},
                       Row{PinArrangement::Inline, true, 0.403, 0.550},
                       Row{PinArrangement::Staggered, false, 0.690, 0.511},
                       Row{PinArrangement::Staggered, true, 0.234, 0.608}}) {
    const auto set = chyu_1990_nu(r.a, r.fillet);
    validate_pin_fin_nu_set(set);
    const auto &g = r.a == PinArrangement::Inline ? g_inl : g_stag;
    const auto nu = evaluate_pin_fin_nu(set, g, Re, 0.7);
    EXPECT_LT(rel(nu.Nu, r.A * std::pow(x, r.B) * std::pow(0.7, 0.4)), 1e-6)
        << set.name;
    EXPECT_FALSE(nu.extrapolated) << set.name;
    EXPECT_EQ(set.surface, PinNuSurface::Pin);
  }
  // Inline straight friction is Re-independent: 0.1693 in Chyu's basis.
  const auto fi = chyu_1990_friction(PinArrangement::Inline, false);
  validate_pin_fin_friction_set(fi);
  EXPECT_EQ(fi.provenance, RibProvenance::Fitted);
  EXPECT_NEAR(evaluate_pin_fin_friction(fi, g_inl, Re).f, 0.1693 / 4.0, 1e-12);
  EXPECT_FALSE(evaluate_pin_fin_friction(fi, g_inl, Re).extrapolated);
  const auto fs = chyu_1990_friction(PinArrangement::Staggered, false);
  EXPECT_LT(rel(evaluate_pin_fin_friction(fs, g_stag, Re).f,
                1.6163 * std::pow(x, -0.1867) / 4.0),
            1e-12);
  // The paper's own reading of its fillet effect: -25 / -17 / -8 %.
  const auto mod = chyu_1990_fillet_over_straight(PinArrangement::Staggered);
  validate_pin_fin_modifier(mod);
  EXPECT_NEAR(evaluate_pin_fin_modifier(mod, g_stag, 5.0e3).ratio_Nu, 0.775,
              0.005);
  EXPECT_NEAR(evaluate_pin_fin_modifier(mod, g_stag, 3.0e4).ratio_Nu, 0.92,
              0.005);
  for (double R : {6.0e3, 2.5e4}) {
    const double straight =
        evaluate_pin_fin_nu(chyu_1990_nu(PinArrangement::Staggered, false),
                            g_stag, R, 0.7)
            .Nu;
    const double fillet =
        evaluate_pin_fin_nu(chyu_1990_nu(PinArrangement::Staggered, true),
                            g_stag, R, 0.7)
            .Nu;
    EXPECT_LT(rel(straight * evaluate_pin_fin_modifier(mod, g_stag, R).ratio_Nu,
                  fillet),
              1e-12);
  }
}
