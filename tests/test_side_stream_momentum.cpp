#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "compressible.h"
#include "side_stream_momentum.h"
#include "stagnation.h"
#include "tee_junction.h"
#include "thermo.h"

using combaero::chamber_merge_face_state;
using combaero::jet_impulse;

namespace {
std::vector<double> air() {
  std::vector<double> X(combaero::num_species(), 0.0);
  X[combaero::species_index_from_name("N2")] = 0.7808;
  X[combaero::species_index_from_name("O2")] = 0.2095;
  X[combaero::species_index_from_name("AR")] = 0.0097;
  return X;
}
}  // namespace

// ---- merge chamber face ----------------------------------------------------

TEST(ChamberMergeFace, NoSideStreamsIsTheChamberItself) {
  const auto X = air();
  const double m = 0.5, P = 2.0e5, T = 900.0, A = 0.005;
  const auto r = chamber_merge_face_state(m, T, X, m, P, T, X, 0.0, A);
  EXPECT_NEAR(r.P_face, P, 1e-9 * P);
  const double rho = combaero::density(T, P, X);
  const double M = m / (rho * A) / combaero::speed_of_sound(T, X);
  EXPECT_NEAR(r.Pt_face, combaero::P0_from_static(P, T, M, X), 1e-6 * P);
  EXPECT_FALSE(r.choked);
}

TEST(ChamberMergeFace, TheImpulseBalanceHoldsAtTheFaceDensity) {
  // Mach ~0.4 at the outlet: compressible, and exact by construction.
  const auto X = air();
  const double mm = 1.0, mo = 1.15, J = 25.0, P = 1.5e5, T = 1100.0, Tm = 1200.0;
  const double A = 0.01;
  const auto r = chamber_merge_face_state(mm, Tm, X, mo, P, T, X, J, A);
  const double rho_f = combaero::density(Tm, r.P_face, X);
  const double rho = combaero::density(T, P, X);
  const double lhs = r.P_face * A + mm * mm / (rho_f * A) + J;
  const double rhs = P * A + mo * mo / (rho * A);
  EXPECT_NEAR(lhs, rhs, 1e-9 * rhs);
  EXPECT_GT(r.M_face, 0.2);
}

TEST(ChamberMergeFace, LowMachIsTheIncompressibleMerge) {
  const auto X = air();
  const double mm = 0.03, mo = 0.036, J = 0.05, P = 1.0e5, T = 300.0, A = 0.01;
  const auto r = chamber_merge_face_state(mm, T, X, mo, P, T, X, J, A);
  const double rho = combaero::density(T, P, X);
  const double incompressible = (mo * mo - mm * mm) / (rho * A * A) - J / A;
  EXPECT_NEAR(r.P_face - P, incompressible, 1e-3 * std::abs(incompressible));
}

TEST(ChamberMergeFace, DerivativesMatchCentralDifferences) {
  const auto X = air();
  const double mm = 0.8, mo = 0.95, J = 12.0, P = 1.6e5, T = 1000.0, Tm = 1150.0;
  const double A = 0.01;
  const auto r = chamber_merge_face_state(mm, Tm, X, mo, P, T, X, J, A);
  auto f = [&](double a, double b, double c, double d, double e2, double g) {
    return chamber_merge_face_state(a, b, X, c, d, e2, X, g, A);
  };
  struct Probe {
    double dPf, dPtf, h;
    int which;
  };
  const Probe probes[] = {
      {r.dPf_dm_main, r.dPtf_dm_main, 1e-6, 0}, {r.dPf_dT_main, r.dPtf_dT_main, 1e-3, 1},
      {r.dPf_dm_out, r.dPtf_dm_out, 1e-6, 2},   {r.dPf_dP, r.dPtf_dP, 1.0, 3},
      {r.dPf_dT, r.dPtf_dT, 1e-3, 4},           {r.dPf_dJ, r.dPtf_dJ, 1e-4, 5}};
  for (const auto& p : probes) {
    double up[6] = {mm, Tm, mo, P, T, J}, dn[6] = {mm, Tm, mo, P, T, J};
    up[p.which] += p.h;
    dn[p.which] -= p.h;
    const auto a = f(up[0], up[1], up[2], up[3], up[4], up[5]);
    const auto b = f(dn[0], dn[1], dn[2], dn[3], dn[4], dn[5]);
    const double fdP = (a.P_face - b.P_face) / (2 * p.h);
    const double fdPt = (a.Pt_face - b.Pt_face) / (2 * p.h);
    EXPECT_NEAR(p.dPf, fdP, 1e-4 * std::abs(fdP) + 1e-6) << p.which;
    EXPECT_NEAR(p.dPtf, fdPt, 1e-4 * std::abs(fdPt) + 1e-6) << p.which;
  }
}

TEST(ChamberMergeFace, AnImpulseThatCannotBeCarriedSubsonicallyIsFlagged) {
  const auto X = air();
  // Side jets pushing harder than the main stream's impulse can absorb: the
  // quadratic P_f^2 - Pi P_f + c has no real root (Pi = 0.84e5 against
  // 2 sqrt(c) = 1.17e5 here), so the face would have to choke.
  const auto r = chamber_merge_face_state(1.0, 1200.0, X, 1.0, 1.0e5, 1200.0, X, 500.0, 0.01);
  EXPECT_TRUE(r.choked);
  EXPECT_TRUE(std::isfinite(r.P_face));
  EXPECT_TRUE(std::isfinite(r.Pt_face));
}

// ---- jet impulse -------------------------------------------------------------

TEST(JetImpulse, LowPressureRatioIsBernoulli) {
  const auto X = air();
  const double Pt = 1.001e5, Tt = 600.0, P = 1.0e5, m = 0.01;
  const auto r = jet_impulse(m, Pt, Tt, P, X);
  const double rho = combaero::density(Tt, P, X);
  EXPECT_NEAR(r.dJ_dm, std::sqrt(2.0 * (Pt - P) / rho), 1e-3 * r.dJ_dm);
  EXPECT_NEAR(r.J, m * r.dJ_dm, 1e-12 * r.J);
  EXPECT_FALSE(r.choked);
}

TEST(JetImpulse, PressureDerivativeIsEuler) {
  // Isentropic expansion: w dw = -dP / rho_static.
  const auto X = air();
  const double Pt = 1.4e5, Tt = 700.0, P = 1.0e5, m = 0.02;
  const auto r = jet_impulse(m, Pt, Tt, P, X);
  const double w = r.dJ_dm;
  const double M = combaero::mach_from_pressure_ratio(Tt, Pt, P, X);
  const double Ts = combaero::T_from_stagnation(Tt, M, X);
  const double rho_s = combaero::density(Ts, P, X);
  EXPECT_NEAR(r.dJ_dP / m, -1.0 / (rho_s * w), 1e-4 / (rho_s * w));
}

TEST(JetImpulse, IsContinuousAndSmoothAcrossChoking) {
  const auto X = air();
  const double Pt = 3.0e5, Tt = 600.0, m = 0.02;
  const double Pstar = combaero::critical_pressure_ratio(Tt, Pt, X) * Pt;
  const auto a = jet_impulse(m, Pt, Tt, Pstar * (1.0 + 1e-6), X);
  const auto b = jet_impulse(m, Pt, Tt, Pstar * (1.0 - 1e-6), X);
  EXPECT_FALSE(a.choked);
  EXPECT_TRUE(b.choked);
  EXPECT_NEAR(a.J, b.J, 1e-4 * a.J);
  EXPECT_NEAR(a.dJ_dP, b.dJ_dP, 1e-2 * std::abs(a.dJ_dP));
  // Below choking the jet still gains thrust as the back pressure falls.
  EXPECT_GT(jet_impulse(m, Pt, Tt, 0.5 * Pstar, X).J, b.J);
}

TEST(JetImpulse, NoOutflowNoMomentum) {
  const auto X = air();
  const auto r = jet_impulse(0.01, 1.0e5, 600.0, 1.0e5, X);
  EXPECT_EQ(r.J, 0.0);
  EXPECT_EQ(r.dJ_dm, 0.0);
}

// ---- one station, half at a time (#471) ------------------------------------

TEST(StationHalfDrop, BleedAtBassettKappaIsBassettK5ForEverySplit) {
  // Constant density: both halves at one state make the whole station.
  const auto X = air();
  const double P = 1.2e5, T = 500.0, A = 4e-3, m_a = 0.4;
  const double rho = combaero::density(T, P, X);
  for (double q : {0.0, 0.2, 0.5, 0.8, 0.95, 1.0}) {
    const double m_b = q * m_a;
    const auto h = combaero::station_half_drop(m_a, m_b, P, T, X, A,
                                               combaero::STATION_KAPPA_BLEED_BASSETT);
    const double drop = 2.0 * h.dP;  // P_a - P_b
    const double ua = m_a / (rho * A), ub = m_b / (rho * A);
    const double K = (drop + 0.5 * rho * (ua * ua - ub * ub)) / (0.5 * rho * ua * ua);
    EXPECT_NEAR(K, combaero::K5(q), 1e-12) << q;
  }
}

TEST(StationHalfDrop, NormalMergeIsTheImpingementTerm) {
  const auto X = air();
  const double P = 1.0e5, T = 300.0, A = 5e-4;
  const auto h = combaero::station_half_drop(0.015, 0.021, P, T, X, A,
                                             combaero::STATION_KAPPA_MERGE_NORMAL);
  const double rho = combaero::density(T, P, X);
  EXPECT_NEAR(h.dP, 0.5 * (0.021 * 0.021 - 0.015 * 0.015) / (rho * A * A), 1e-12 * h.dP);
}

TEST(StationHalfDrop, DerivativesMatchCentralDifferences) {
  const auto X = air();
  const double P = 1.1e5, T = 420.0, A = 2e-3, e = 1e-8;
  for (double kappa : {0.0, 0.75, 1.0}) {
    for (double ma : {-0.05, 0.08}) {
      for (double mb : {-0.03, 0.06}) {
        const auto r = combaero::station_half_drop(ma, mb, P, T, X, A, kappa);
        auto f = [&](double a, double b, double p, double t) {
          return combaero::station_half_drop(a, b, p, t, X, A, kappa).dP;
        };
        const double fa = (f(ma + e, mb, P, T) - f(ma - e, mb, P, T)) / (2 * e);
        const double fb = (f(ma, mb + e, P, T) - f(ma, mb - e, P, T)) / (2 * e);
        const double fP = (f(ma, mb, P + 1.0, T) - f(ma, mb, P - 1.0, T)) / 2.0;
        const double fT = (f(ma, mb, P, T + 1e-3) - f(ma, mb, P, T - 1e-3)) / 2e-3;
        // Absolute floor: round-off of a ~1e2 Pa value over a 1e-8 step,
        // where the analytic slope can be exactly 0 (2|m_b| = kappa m_a).
        EXPECT_NEAR(r.d_dm_a, fa, 1e-6 * std::abs(fa) + 1e-4);
        EXPECT_NEAR(r.d_dm_b, fb, 1e-6 * std::abs(fb) + 1e-4);
        EXPECT_NEAR(r.d_dP, fP, 1e-5 * std::abs(fP) + 1e-12);
        EXPECT_NEAR(r.d_dT, fT, 1e-5 * std::abs(fT) + 1e-12);
      }
    }
  }
}

TEST(ChannelEntryDrop, IsTheDynamicHeadPlusTheEntryLoss) {
  const auto X = air();
  const double P = 1.0e5, T = 300.0, A = 1e-3, m = 0.05;
  const double rho = combaero::density(T, P, X);
  const double u = m / (rho * A);
  const auto r = combaero::channel_entry_drop(m, P, T, X, A, 0.5);
  EXPECT_NEAR(r.dP, 1.5 * 0.5 * rho * u * u, 1e-9 * r.dP);
  const double e = 1e-8;
  const double fd = (combaero::channel_entry_drop(m + e, P, T, X, A, 0.5).dP -
                     combaero::channel_entry_drop(m - e, P, T, X, A, 0.5).dP) / (2 * e);
  EXPECT_NEAR(r.d_dm_a, fd, 1e-6 * fd);
}
