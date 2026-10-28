#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>

#include "crossflow_cd.h"

using combaero::CrossflowCd;
using combaero::florschuetz_metzger_su_1982_cd;
using combaero::florschuetz_metzger_su_1982_has_plate;
using combaero::rogers_hersh_1976_cd;

namespace {

constexpr double kGc = 40.0;  // crossflow mass velocity [kg/(m^2 s)]

// Central differences of a law in (G_ideal, G_c).
template <class F>
void expect_partials(F law, double Gi, double Gc, double rel) {
  const CrossflowCd r = law(Gi, Gc);
  const double hi = Gi * 1e-6;
  const double hc = Gc * 1e-6;
  const double fd_i = (law(Gi + hi, Gc).Cd - law(Gi - hi, Gc).Cd) / (2.0 * hi);
  const double fd_c = (law(Gi, Gc + hc).Cd - law(Gi, Gc - hc).Cd) / (2.0 * hc);
  EXPECT_NEAR(r.dCd_dG_ideal, fd_i, rel * std::abs(fd_i) + 1e-14);
  EXPECT_NEAR(r.dCd_dG_c, fd_c, rel * std::abs(fd_c) + 1e-14);
}

double zeta(double x) { return 1.06 / std::pow(x + 0.806, 0.602); }

}  // namespace

// ---- Rogers & Hersh (1976) --------------------------------------------------

TEST(RogersHersh1976, ResolvesThePapersFitInTheActualJetVelocity) {
  // The paper: Cd = a (Vi/V_inf)^b with Vi the ACTUAL jet velocity, Gj = Cd G_ideal.
  for (const bool suction : {false, true}) {
    const double a = suction ? 0.80 : 0.97;
    const double b = suction ? 0.47 : 0.534;
    for (const double Gi : {0.05, 0.5, 5.0, 20.0}) {
      const CrossflowCd r = rogers_hersh_1976_cd(0.85, Gi, kGc, suction);
      ASSERT_FALSE(r.capped);
      ASSERT_LT(r.Cd, 0.9 * 0.85);  // below the cap's blend band
      const double ratio = r.Cd * Gi / kGc;  // Gj / Gc
      EXPECT_NEAR(r.Cd, a * std::pow(ratio, b), 1e-12);
    }
  }
}

TEST(RogersHersh1976, HoldsTheHolesOwnCdUnderAWeakCrossflow) {
  EXPECT_TRUE(rogers_hersh_1976_cd(0.79, 200.0, kGc, false).capped);
  EXPECT_DOUBLE_EQ(rogers_hersh_1976_cd(0.79, 200.0, kGc, false).Cd, 0.79);
  EXPECT_EQ(rogers_hersh_1976_cd(0.79, 200.0, kGc, false).dCd_dG_c, 0.0);
  // No crossflow at all: the hole's own Cd, whatever the drive.
  EXPECT_TRUE(rogers_hersh_1976_cd(0.79, 1e-6, 0.0, false).capped);
  // Injection reaches the 0.79 cap at Gc/Gj = 1.47: Florschuetz's arrays
  // (Gc/Gj <= 0.77) never leave it.
  const double q = std::pow(0.79 / std::pow(0.97, 1.0 / 0.466), 0.466 / 0.534);
  EXPECT_NEAR(1.0 / (0.79 * q), 1.469, 2e-3);
}

TEST(RogersHersh1976, JoinsTheCapC1) {
  // Across the blend band the Cd and its slope are continuous, and it never
  // exceeds Cd0 or the fit.
  double prev = 0.0;
  double prev_slope = 0.0;
  for (double Gi = 25.0; Gi <= 60.0; Gi += 0.01) {
    const CrossflowCd r = rogers_hersh_1976_cd(0.85, Gi, kGc, false);
    const double fit = std::pow(0.97, 1.0 / 0.466) * std::pow(Gi / kGc, 0.534 / 0.466);
    EXPECT_LE(r.Cd, 0.85 + 1e-15);
    EXPECT_LE(r.Cd, fit + 1e-15);
    if (prev > 0.0) {
      // The secant over the step is the mean of the end slopes (C1, no jump).
      const double secant = (r.Cd - prev) / 0.01;
      EXPECT_NEAR(secant, 0.5 * (r.dCd_dG_ideal + prev_slope), 1e-5) << Gi;
    }
    prev = r.Cd;
    prev_slope = r.dCd_dG_ideal;
  }
}

TEST(RogersHersh1976, NoDriveNoFlow) {
  EXPECT_EQ(rogers_hersh_1976_cd(0.85, 0.0, kGc, false).Cd, 0.0);
  EXPECT_EQ(rogers_hersh_1976_cd(0.85, -1.0, kGc, true).Cd, 0.0);
}

TEST(RogersHersh1976, AtTheStagnantRowsOfCr3630) {
  // Gc/Gj ~ 125 (B(5,4,1) mc/mj 0.99, rows 2-4): Cd ~ 0.07-0.08 either way.
  for (const bool suction : {false, true}) {
    double Gi = 1.0;
    for (int it = 0; it < 60; ++it) {  // find Gi with Gj/Gc = 1/125
      const double cd = rogers_hersh_1976_cd(0.85, Gi, kGc, suction).Cd;
      Gi *= std::sqrt((kGc / 125.0) / (cd * Gi));
    }
    const double cd = rogers_hersh_1976_cd(0.85, Gi, kGc, suction).Cd;
    EXPECT_NEAR(cd, suction ? 0.8 * std::pow(0.008, 0.47) : 0.97 * std::pow(0.008, 0.534), 1e-6);
    EXPECT_LT(cd, 0.09);
  }
}

TEST(RogersHersh1976, AnalyticPartials) {
  for (const bool suction : {false, true}) {
    auto law = [suction](double Gi, double Gc) {
      return rogers_hersh_1976_cd(0.85, Gi, Gc, suction);
    };
    for (const double Gi : {0.02, 0.4, 8.0, 30.0, 38.0, 41.0}) {
      expect_partials(law, Gi, kGc, 1e-6);
    }
  }
}

TEST(RogersHersh1976, RefusesANonPositiveCd0) {
  EXPECT_THROW(rogers_hersh_1976_cd(0.0, 1.0, kGc, false), std::invalid_argument);
}

// ---- Florschuetz, Metzger & Su (1982), CR-3630 Table 4.1 --------------------

TEST(FlorschuetzMetzgerSu1982, CoversTheTestedPlatesOnly) {
  EXPECT_TRUE(florschuetz_metzger_su_1982_has_plate(5, 4, 1));
  EXPECT_TRUE(florschuetz_metzger_su_1982_has_plate(10, 8, 1));
  EXPECT_FALSE(florschuetz_metzger_su_1982_has_plate(5, 8, 2));
  EXPECT_FALSE(florschuetz_metzger_su_1982_has_plate(15, 4, 1));
  EXPECT_THROW(florschuetz_metzger_su_1982_cd(5, 8, 2, 1.0, kGc), std::invalid_argument);
}

// The law in x = Gc/Gj as implemented: the printed pieces joined
// continuously (see crossflow_cd.h).
double table41(double xn, double yn, double z, double x) {
  if (xn == 5 && yn == 4) {
    const double mult = z == 1 ? 1.0 : (z == 2 ? 1.08 : 1.11);
    return std::min(0.85, mult * zeta(x));
  }
  const bool b581 = xn == 5 && yn == 8;
  const double A = b581 ? -0.169 : -0.128;
  const double end = b581 ? 1.5 : 1.8;
  const double cd0 = b581 ? 0.80 : 0.76;
  const double B = zeta(end) - A * end;
  return x <= end ? std::min(cd0, A * x + B) : zeta(x);
}

TEST(FlorschuetzMetzgerSu1982, IsTheTableInGcOverGj) {
  struct Case {
    double xn, yn, z;
  };
  for (const Case c : {Case{5, 4, 1}, Case{5, 4, 2}, Case{5, 4, 3}, Case{5, 8, 1}, Case{10, 4, 1}}) {
    for (const double Gi : {0.05, 0.3, 2.0, 8.0, 30.0, 45.0, 52.0, 56.0, 60.0, 80.0, 400.0}) {
      const CrossflowCd r = florschuetz_metzger_su_1982_cd(c.xn, c.yn, c.z, Gi, kGc);
      const double x = kGc / (r.Cd * Gi);  // Gc / Gj at the solved Cd
      const double cd0 = c.xn == 10 ? 0.76 : (c.yn == 8 ? 0.80 : 0.85);
      if (r.capped) {
        EXPECT_DOUBLE_EQ(r.Cd, cd0);
      } else if (r.Cd <= 0.9 * cd0) {
        // Below the join band: the table's falling piece exactly.
        EXPECT_NEAR(r.Cd, table41(c.xn, c.yn, c.z, x), 1e-12)
            << c.xn << "," << c.yn << "," << c.z << " Gi " << Gi;
      } else {
        // In the band: between the piece and cd0, within 2.5% of either.
        EXPECT_LE(r.Cd, cd0);
        EXPECT_NEAR(r.Cd, table41(c.xn, c.yn, c.z, x), 0.025 * cd0);
      }
    }
  }
}

TEST(FlorschuetzMetzgerSu1982, StaysWithinThePrintedDigitsOfThePrintedPieces) {
  // The printed Table 4.1, piece by piece, against the continuous law.
  for (double x = 0.0; x <= 6.0; x += 0.01) {
    const double p541 = x <= 0.63 ? 0.85 : zeta(x);
    const double p581 = x <= 0.54 ? 0.80 : (x <= 1.5 ? -0.169 * x + 0.893 : zeta(x));
    const double p1041 = x <= 0.54 ? 0.76 : (x <= 1.8 ? -0.128 * x + 0.825 : zeta(x));
    EXPECT_NEAR(table41(5, 4, 1, x), p541, 3e-3) << x;
    EXPECT_NEAR(table41(5, 8, 1, x), p581, 3e-3) << x;
    EXPECT_NEAR(table41(10, 4, 1, x), p1041, 4e-3) << x;
  }
}

TEST(FlorschuetzMetzgerSu1982, AnalyticPartials) {
  auto law41 = [](double Gi, double Gc) { return florschuetz_metzger_su_1982_cd(5, 4, 1, Gi, Gc); };
  auto law81 = [](double Gi, double Gc) { return florschuetz_metzger_su_1982_cd(5, 8, 1, Gi, Gc); };
  for (const double Gi : {0.05, 1.0, 10.0, 40.0, 52.0, 56.0}) {
    expect_partials(law41, Gi, kGc, 1e-6);
  }
  // On B(5,8,1)'s linear piece (x 0.54-1.5) and on zeta beyond it.
  for (const double Gi : {35.0, 60.0, 5.0}) {
    expect_partials(law81, Gi, kGc, 1e-6);
  }
}
