// The C2 re-fit of the NASA-9 thermo (#549): continuity at every former join,
// equivalence with NASA-9 (identical outside the bridge windows, bounded
// inside), and h' = cp, T s' = cp inside the bridges. The continuity test
// runs against BOTH models: it must fail NASA-9 (whose dcp/dT steps at the
// joins) and pass the re-fit -- a check that cannot tell them apart would be
// vacuous.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "thermo.h"
#include "thermo_c2_refit_data.h"

using combaero::C2Bridge;
using combaero::c2_bridges;
using combaero::ThermoModel;

namespace {

// Restores the default model however a test exits.
struct ModelGuard {
  explicit ModelGuard(ThermoModel m) { combaero::set_thermo_model(m); }
  ~ModelGuard() { combaero::set_thermo_model(ThermoModel::NASA9); }
};

double H_R(std::size_t i, double T) { return combaero::h_RT(i, T) * T; }

// Largest relative dcp/dT jump across any join, and largest relative jump of
// d2cp/dT2 across any bridge edge, for the model currently selected.
struct Jumps {
  double c1 = 0.0;
  double c2 = 0.0;
};

Jumps measure_jumps() {
  Jumps out;
  for (std::size_t i = 0; i < c2_bridges.size(); ++i) {
    for (const C2Bridge& b : c2_bridges[i]) {
      const double e = 1e-7;
      const double lo = combaero::dcp_R_dT(i, b.T_join - e);
      const double hi = combaero::dcp_R_dT(i, b.T_join + e);
      const double scale = std::max({std::abs(lo), std::abs(hi), 1e-6});
      out.c1 = std::max(out.c1, std::abs(hi - lo) / scale);
      // d2cp/dT2 one-sided on each side of each window edge. Measured as its
      // effect on cp across the window, |jump| w^2 / cp: a relative jump
      // blows up where d2cp/dT2 is exactly 0 (AR's constant low-T cp).
      for (double edge : {b.T_join - b.half_width, b.T_join + b.half_width}) {
        const double h = 1e-4;
        const double inner =
            (combaero::dcp_R_dT(i, edge - h) - combaero::dcp_R_dT(i, edge - 2 * h)) / h;
        const double outer =
            (combaero::dcp_R_dT(i, edge + 2 * h) - combaero::dcp_R_dT(i, edge + h)) / h;
        out.c2 = std::max(out.c2, std::abs(outer - inner) * b.half_width * b.half_width /
                                      combaero::cp_R(i, edge));
      }
    }
  }
  return out;
}

}  // namespace

TEST(ThermoC2Refit, DefaultIsTheOriginalNasa9) {
  EXPECT_EQ(combaero::thermo_model(), ThermoModel::NASA9);
}

TEST(ThermoC2Refit, DcpDtIsContinuousAtEveryJoinOnlyInTheRefit) {
  {
    ModelGuard g(ThermoModel::NASA9);
    // The published polynomials step: 0.15 % (CH4) to 100 % (AR) at 1000 K,
    // up to 48 % (CO) at 6000 K. The check must see it.
    EXPECT_GT(measure_jumps().c1, 1e-2);
  }
  {
    ModelGuard g(ThermoModel::NASA9_C2);
    const Jumps j = measure_jumps();
    EXPECT_LT(j.c1, 1e-6);
    // d2cp/dT2 matched at the window edges by construction; the bound is the
    // one-sided difference's own truncation.
    EXPECT_LT(j.c2, 1e-6) << "c2 measure " << j.c2;
  }
}

TEST(ThermoC2Refit, IdenticalToNasa9OutsideTheWindows) {
  const double Ts[] = {250.0, 298.15, 500.0, 940.0, 1060.0, 3000.0, 5600.0, 6400.0, 15000.0};
  for (std::size_t i = 0; i < c2_bridges.size(); ++i) {
    for (double T : Ts) {
      bool inside = false;
      for (const C2Bridge& b : c2_bridges[i]) {
        inside = inside || std::abs(T - b.T_join) < b.half_width;
      }
      if (inside) {
        continue;
      }
      double cp0, h0, s0;
      {
        ModelGuard g(ThermoModel::NASA9);
        cp0 = combaero::cp_R(i, T);
        h0 = combaero::h_RT(i, T);
        s0 = combaero::s_R(i, T);
      }
      ModelGuard g(ThermoModel::NASA9_C2);
      EXPECT_EQ(combaero::cp_R(i, T), cp0) << i << " " << T;
      EXPECT_EQ(combaero::h_RT(i, T), h0) << i << " " << T;
      EXPECT_EQ(combaero::s_R(i, T), s0) << i << " " << T;
    }
  }
}

TEST(ThermoC2Refit, InsideTheWindowsTheDeviationIsBounded) {
  // Measured worst cases: cp 4.8e-4 (AR at 6000 K), h 1.1e-5 cp T, s 3.1e-5 R.
  for (std::size_t i = 0; i < c2_bridges.size(); ++i) {
    for (const C2Bridge& b : c2_bridges[i]) {
      for (int k = -40; k <= 40; ++k) {
        const double T = b.T_join + b.half_width * k / 41.0;
        double cp0, H0, s0;
        {
          ModelGuard g(ThermoModel::NASA9);
          cp0 = combaero::cp_R(i, T);
          H0 = H_R(i, T);
          s0 = combaero::s_R(i, T);
        }
        ModelGuard g(ThermoModel::NASA9_C2);
        EXPECT_LT(std::abs(combaero::cp_R(i, T) - cp0) / cp0, 1e-3) << i << " " << T;
        EXPECT_LT(std::abs(H_R(i, T) - H0) / (cp0 * T), 2e-5) << i << " " << T;
        EXPECT_LT(std::abs(combaero::s_R(i, T) - s0), 5e-5) << i << " " << T;
      }
    }
  }
}

TEST(ThermoC2Refit, EnthalpyAndEntropyAreTheIntegralsOfCpInsideTheBridges) {
  ModelGuard g(ThermoModel::NASA9_C2);
  for (std::size_t i = 0; i < c2_bridges.size(); ++i) {
    for (const C2Bridge& b : c2_bridges[i]) {
      for (int k = -9; k <= 9; ++k) {
        const double T = b.T_join + 0.1 * b.half_width * k;
        const double e = 1e-3;
        const double cp = combaero::cp_R(i, T);
        const double dH = (H_R(i, T + e) - H_R(i, T - e)) / (2 * e);
        const double dS = (combaero::s_R(i, T + e) - combaero::s_R(i, T - e)) / (2 * e);
        EXPECT_NEAR(dH, cp, 1e-7 * cp) << i << " " << T;
        EXPECT_NEAR(dS * T, cp, 1e-7 * cp) << i << " " << T;
      }
    }
  }
}
