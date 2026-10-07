#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

#include "side_stream_momentum.h"

using combaero::side_stream_momentum_drop;

// Half of each adjacent station's drop: stations (m_arr -> m) and (m -> m_out).
TEST(SideStreamMomentum, IsHalfOfBothAdjacentStations) {
  const double m_arr = 0.015, m = 0.02, m_out = 0.026, rho = 1.2, A = 5e-4;
  const auto r = side_stream_momentum_drop(m_arr, m_out, rho, A);
  const double up = (m * m - m_arr * m_arr) / (rho * A * A);
  const double dn = (m_out * m_out - m * m) / (rho * A * A);
  EXPECT_NEAR(r.dP, 0.5 * (up + dn), 1e-9 * r.dP);
}

TEST(SideStreamMomentum, DerivativesMatchCentralDifferences) {
  const double rho = 1.1, A = 3e-4, e = 1e-8;
  // Off zero: m|m| is C1 there with zero slope, which the analytic form
  // returns exactly, but a central difference leaves an O(step) residue.
  for (double m_arr : {-0.004, 0.003, 0.01}) {
    for (double m_out : {-0.002, 0.005, 0.02}) {
      const auto r = side_stream_momentum_drop(m_arr, m_out, rho, A);
      auto f = [&](double a, double b, double p) {
        return side_stream_momentum_drop(a, b, p, A).dP;
      };
      const double fd_out = (f(m_arr, m_out + e, rho) - f(m_arr, m_out - e, rho)) / (2 * e);
      const double fd_arr = (f(m_arr + e, m_out, rho) - f(m_arr - e, m_out, rho)) / (2 * e);
      EXPECT_NEAR(r.d_dm_out, fd_out, 1e-6 * std::abs(fd_out));
      EXPECT_NEAR(r.d_dm_arr, fd_arr, 1e-6 * std::abs(fd_arr));
      EXPECT_NEAR(r.d_drho, (f(m_arr, m_out, rho + 1e-6) - f(m_arr, m_out, rho - 1e-6)) / 2e-6,
                  1e-6 + 1e-6 * std::abs(r.d_drho));
    }
  }
}

TEST(SideStreamMomentum, ZeroFlowHasZeroSlope) {
  const auto r = side_stream_momentum_drop(0.0, 0.0, 1.2, 1e-3);
  EXPECT_EQ(r.d_dm_out, 0.0);
  EXPECT_EQ(r.d_dm_arr, 0.0);
}

TEST(SideStreamMomentum, NoSideStreamNoChange) {
  // Arriving equals leaving with nothing joining: the stations cancel.
  // Near, not equal: FMA contraction leaves a 1e-16 residue in a - a.
  EXPECT_NEAR(side_stream_momentum_drop(0.01, 0.01, 1.2, 1e-3).dP, 0.0, 1e-12);
  EXPECT_THROW(side_stream_momentum_drop(0.01, 0.02, 0.0, 1e-3), std::invalid_argument);
  EXPECT_THROW(side_stream_momentum_drop(0.01, 0.02, 1.2, 0.0), std::invalid_argument);
}
