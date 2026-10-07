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

// Merge chamber (#471): the main-inlet face relative to the chamber state.
TEST(ChamberMerge, NoSideStreamsIsTheSingleInletChamber) {
  const auto r = combaero::chamber_merge_face_offset(0.3, 0.3, 0.0, 1.1, 0.01);
  EXPECT_NEAR(r.dP_face, 0.0, 1e-12);
  EXPECT_NEAR(r.dPt_face, 0.0, 1e-12);
}

TEST(ChamberMerge, NormalInjectionIsTheSideStreamStation) {
  // At 90 deg the static drop over one station is the full side-stream drop:
  // twice the centred segment's half, with m_arr = m_main, m_out = m_out.
  const double mm = 0.3, mo = 0.36, rho = 1.1, A = 0.01;
  const auto r = combaero::chamber_merge_face_offset(mm, mo, 0.0, rho, A);
  const auto half = combaero::side_stream_momentum_drop(mm, mo, rho, A);
  EXPECT_NEAR(r.dP_face, 2.0 * half.dP, 1e-9 * r.dP_face);
  // Mixing loss: Pt drops by half the static drop.
  EXPECT_NEAR(r.dPt_face, 0.5 * r.dP_face, 1e-9 * r.dP_face);
}

TEST(ChamberMerge, AxialSideMomentumPushesTheChamber) {
  const double mm = 0.3, mo = 0.36, rho = 1.1, A = 0.01;
  const double S = 0.06 * 40.0;  // 0.06 kg/s at 40 m/s, axial
  const auto r0 = combaero::chamber_merge_face_offset(mm, mo, 0.0, rho, A);
  const auto r = combaero::chamber_merge_face_offset(mm, mo, S, rho, A);
  EXPECT_NEAR(r0.dP_face - r.dP_face, S / A, 1e-9 * S / A);
}

TEST(ChamberMerge, DerivativesMatchCentralDifferences) {
  const double mm = 0.3, mo = 0.36, S = 1.5, rho = 1.1, A = 0.01, e = 1e-7;
  const auto r = combaero::chamber_merge_face_offset(mm, mo, S, rho, A);
  auto f = [&](double a, double b, double c, double d) {
    return combaero::chamber_merge_face_offset(a, b, c, d, A);
  };
  auto fdP = [&](double da, double db, double dc, double dd) {
    return (f(mm + da, mo + db, S + dc, rho + dd).dP_face -
            f(mm - da, mo - db, S - dc, rho - dd).dP_face) / (2 * e);
  };
  auto fdPt = [&](double da, double db, double dc, double dd) {
    return (f(mm + da, mo + db, S + dc, rho + dd).dPt_face -
            f(mm - da, mo - db, S - dc, rho - dd).dPt_face) / (2 * e);
  };
  EXPECT_NEAR(r.dP_dm_main, fdP(e, 0, 0, 0), 1e-6 * std::abs(r.dP_dm_main));
  EXPECT_NEAR(r.dP_dm_out, fdP(0, e, 0, 0), 1e-6 * std::abs(r.dP_dm_out));
  EXPECT_NEAR(r.dP_dS, fdP(0, 0, e, 0), 1e-6 * std::abs(r.dP_dS));
  EXPECT_NEAR(r.dP_drho, fdP(0, 0, 0, e), 1e-5 * std::abs(r.dP_drho));
  EXPECT_NEAR(r.dPt_dm_main, fdPt(e, 0, 0, 0), 1e-6 * std::abs(r.dPt_dm_main));
  EXPECT_NEAR(r.dPt_dm_out, fdPt(0, e, 0, 0), 1e-6 * std::abs(r.dPt_dm_out));
  EXPECT_NEAR(r.dPt_dS, fdPt(0, 0, e, 0), 1e-6 * std::abs(r.dPt_dS));
  EXPECT_NEAR(r.dPt_drho, fdPt(0, 0, 0, e), 1e-5 * std::abs(r.dPt_drho));
}
