#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "friction.h"
#include "heat_transfer.h"

using combaero::friction_channel_and_derivative;
using combaero::friction_turbulent_and_derivative;
using combaero::nusselt_channel_gnielinski_and_derivative;
using combaero::nusselt_gnielinski_with_derivative;
using combaero::NU_LAMINAR_CONST_T;

namespace {
constexpr double kRough[] = {0.0, 1e-4, 1e-3, 1e-2};
}

// #448: each regime is EXACT outside its band -- only the bands move.
TEST(ChannelTransitions, RegimesAreExactOutsideTheirBands) {
  for (double eD : kRough) {
    for (double Re : {100.0, 1500.0, 2300.0}) {
      EXPECT_DOUBLE_EQ(friction_channel_and_derivative(Re, eD).f, 64.0 / Re);
    }
    // Smooth Petukhov from 3000; a rough wall's band starts there too, so
    // only a smooth wall stays on Petukhov above it.
    for (double Re : {3000.0, 3500.0}) {
      if (eD > 0.0 && Re > 3000.0) {
        continue;
      }
      EXPECT_DOUBLE_EQ(friction_channel_and_derivative(Re, eD).f,
                       friction_petukhov_clamped(Re))
          << "smooth Petukhov, e/D " << eD << " Re " << Re;
    }
    for (double Re : {4000.0, 1.0e5}) {
      const double expect =
          eD > 0.0 ? friction_colebrook(Re, eD) : friction_petukhov_clamped(Re);
      EXPECT_DOUBLE_EQ(friction_channel_and_derivative(Re, eD).f, expect);
    }
  }
  EXPECT_EQ(friction_channel_and_derivative(0.0, 0.0).f, 0.0);
  EXPECT_EQ(friction_channel_and_derivative(-5.0, 0.0).f, 0.0);

  for (double Re : {500.0, 2300.0}) {
    EXPECT_DOUBLE_EQ(
        nusselt_channel_gnielinski_and_derivative(Re, 0.7, 0.04, 0.0, NU_LAMINAR_CONST_T).Nu,
        NU_LAMINAR_CONST_T);
  }
  for (double Re : {3000.0, 1.0e5}) {
    const auto ft = friction_turbulent_and_derivative(Re, 0.0);
    EXPECT_DOUBLE_EQ(
        nusselt_channel_gnielinski_and_derivative(Re, 0.7, ft.f, ft.df_dRe, NU_LAMINAR_CONST_T).Nu,
        nusselt_gnielinski_with_derivative(Re, 0.7, ft.f, ft.df_dRe).Nu);
  }
}

// The jumps the issue reported are gone: value AND slope continuous at every
// band edge, smooth and rough.
TEST(ChannelTransitions, NoJumpAtAnyBandEdge) {
  for (double eD : kRough) {
    for (double edge : {2300.0, 3000.0, 4000.0}) {
      const double h = 1e-7 * edge;
      const auto lo = friction_channel_and_derivative(edge - h, eD);
      const auto hi = friction_channel_and_derivative(edge + h, eD);
      EXPECT_NEAR(lo.f, hi.f, 1e-6 * hi.f) << "f, e/D " << eD << " Re " << edge;
      EXPECT_NEAR(lo.df_dRe, hi.df_dRe, 1e-4 * std::abs(hi.df_dRe) + 1e-11)
          << "df/dRe, e/D " << eD << " Re " << edge;

      const auto tlo = friction_turbulent_and_derivative(edge - h, eD);
      const auto thi = friction_turbulent_and_derivative(edge + h, eD);
      const auto nlo = nusselt_channel_gnielinski_and_derivative(
          edge - h, 0.7, tlo.f, tlo.df_dRe, NU_LAMINAR_CONST_T);
      const auto nhi = nusselt_channel_gnielinski_and_derivative(
          edge + h, 0.7, thi.f, thi.df_dRe, NU_LAMINAR_CONST_T);
      EXPECT_NEAR(nlo.Nu, nhi.Nu, 1e-6 * nhi.Nu) << "Nu, e/D " << eD << " Re " << edge;
      // Absolute floor: at a band's lower edge the slope is the smoothstep's,
      // O(h) and tending to 0 -- typical dNu/dRe here is ~1e-3.
      EXPECT_NEAR(nlo.dNu_dRe, nhi.dNu_dRe, 1e-4 * std::abs(nhi.dNu_dRe) + 1e-7)
          << "dNu/dRe, e/D " << eD << " Re " << edge;
    }
  }
}

TEST(ChannelTransitions, DerivativesMatchCentralDifferences) {
  for (double eD : kRough) {
    for (double Re : {1500.0, 2400.0, 2700.0, 2999.0, 3001.0, 3500.0, 3999.0,
                      4001.0, 2.0e4, 1.0e6}) {
      const double h = Re * 1e-6;
      const double fd_f = (friction_channel_and_derivative(Re + h, eD).f -
                           friction_channel_and_derivative(Re - h, eD).f) / (2.0 * h);
      const double an_f = friction_channel_and_derivative(Re, eD).df_dRe;
      EXPECT_LT(std::abs(an_f - fd_f) / std::max(std::abs(fd_f), 1e-14), 1e-5)
          << "f, e/D " << eD << " Re " << Re;

      auto nu = [&](double r) {
        const auto ft = friction_turbulent_and_derivative(r, eD);
        return nusselt_channel_gnielinski_and_derivative(r, 0.7, ft.f, ft.df_dRe,
                                                         NU_LAMINAR_CONST_T);
      };
      const double fd_n = (nu(Re + h).Nu - nu(Re - h).Nu) / (2.0 * h);
      const double an_n = nu(Re).dNu_dRe;
      EXPECT_LT(std::abs(an_n - fd_n) / std::max(std::abs(fd_n), 1e-12), 1e-5)
          << "Nu, e/D " << eD << " Re " << Re;
    }
  }
}

TEST(ChannelTransitions, ColebrookDerivativeIsExact) {
  for (double eD : {0.0, 1e-4, 1e-2}) {
    for (double Re : {5000.0, 1.0e5, 1.0e7}) {
      const double h = Re * 1e-6;
      const double fd = (friction_colebrook(Re + h, eD) - friction_colebrook(Re - h, eD)) /
                        (2.0 * h);
      EXPECT_LT(std::abs(friction_colebrook_dRe(Re, eD) - fd) / std::abs(fd), 1e-5)
          << "e/D " << eD << " Re " << Re;
    }
  }
}

// No jump ANYWHERE in the bands, not only at their edges: a fine scan from
// below 2300 to above 4000 bounds every step-to-step change. The old hard
// switches moved Nu by +78% and f by +64% in one step.
TEST(ChannelTransitions, ContinuousAcrossTheWholeTransitionRange) {
  for (double eD : kRough) {
    double f_prev = friction_channel_and_derivative(2200.0, eD).f;
    auto nu = [&](double r) {
      const auto ft = friction_turbulent_and_derivative(r, eD);
      return nusselt_channel_gnielinski_and_derivative(r, 0.7, ft.f, ft.df_dRe,
                                                       NU_LAMINAR_CONST_T)
          .Nu;
    };
    double nu_prev = nu(2200.0);
    for (double Re = 2200.5; Re <= 4200.0; Re += 0.5) {
      const double f = friction_channel_and_derivative(Re, eD).f;
      const double n = nu(Re);
      EXPECT_LT(std::abs(f / f_prev - 1.0), 2e-3) << "f, e/D " << eD << " Re " << Re;
      EXPECT_LT(std::abs(n / nu_prev - 1.0), 2e-3) << "Nu, e/D " << eD << " Re " << Re;
      f_prev = f;
      nu_prev = n;
    }
  }
}
