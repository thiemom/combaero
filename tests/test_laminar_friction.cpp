// The network friction kernels are laminar below transition: f Re = 64
// (Hagen-Poiseuille) down to essentially zero flow, for every model, and the
// turbulent fits are untouched above it. The laminar branch used to share the
// turbulent fits' Re floor of 64, so f saturated at 1 and a channel's drop
// went as m^2 at low flow; Colebrook had no laminar branch at all.

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <tuple>

#include "friction.h"
#include "solver_interface.h"

using combaero::solver::friction_and_jacobian;

namespace {

const char* const kModels[] = {"haaland", "serghides", "colebrook"};

double f_of(const std::string& m, double Re, double eD) {
  return std::get<0>(friction_and_jacobian(m, Re, eD).result);
}

double df_of(const std::string& m, double Re, double eD) {
  return std::get<1>(friction_and_jacobian(m, Re, eD).result);
}

}  // namespace

TEST(LaminarFriction, IsHagenPoiseuilleBelowTransition) {
  for (const char* m : kModels) {
    for (double Re : {0.01, 0.1, 1.0, 10.0, 64.0, 200.0, 500.0, 1000.0}) {
      for (double eD : {0.0, 1e-3}) {
        // The blend's turbulent share at Re 1000 is 1e-3: the bound covers it.
        EXPECT_NEAR(f_of(m, Re, eD) * Re, 64.0, 64.0 * 2e-3) << m << " Re " << Re;
      }
    }
  }
}

TEST(LaminarFriction, LeavesTheTurbulentFitsUntouched) {
  for (double Re : {1e4, 1e5, 1e6}) {
    for (double eD : {0.0, 1e-4, 1e-2}) {
      EXPECT_NEAR(f_of("haaland", Re, eD), friction_haaland(Re, eD),
                  1e-12 * friction_haaland(Re, eD));
      EXPECT_NEAR(f_of("serghides", Re, eD), friction_serghides(Re, eD),
                  1e-12 * friction_serghides(Re, eD));
      EXPECT_NEAR(f_of("colebrook", Re, eD), friction_colebrook(Re, eD),
                  1e-9 * friction_colebrook(Re, eD));
    }
  }
}

// A channel's drop, f Re |m| ~ m: finite and linear as the flow stops.
TEST(LaminarFriction, KeepsTheDropLinearAsTheFlowStops) {
  for (const char* m : kModels) {
    const double ratio = (f_of(m, 1e-2, 0.0) * 1e-2) / (f_of(m, 1.0, 0.0) * 1.0);
    EXPECT_NEAR(ratio, 1.0, 1e-6) << m;
    EXPECT_TRUE(std::isfinite(f_of(m, 0.0, 0.0))) << m;
  }
}

TEST(LaminarFriction, DerivativesMatchCentralDifferences) {
  for (const char* m : kModels) {
    for (double Re : {0.002, 1.0, 50.0, 800.0, 1900.0, 2300.0, 2700.0, 3500.0, 2e4}) {
      for (double eD : {0.0, 1e-3}) {
        const double h = 1e-6 * Re;
        const double fd = (f_of(m, Re + h, eD) - f_of(m, Re - h, eD)) / (2.0 * h);
        EXPECT_NEAR(df_of(m, Re, eD), fd, 1e-6 * std::abs(fd) + 1e-12)
            << m << " Re " << Re << " eD " << eD;
      }
    }
  }
}
