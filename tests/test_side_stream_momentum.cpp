#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "compressible.h"
#include "side_stream_momentum.h"
#include "stagnation.h"
#include "thermo.h"

using combaero::chamber_merge_face_state;
using combaero::jet_impulse;
using combaero::side_stream_momentum_drop;

namespace {
std::vector<double> air() {
  std::vector<double> X(combaero::num_species(), 0.0);
  X[combaero::species_index_from_name("N2")] = 0.7808;
  X[combaero::species_index_from_name("O2")] = 0.2095;
  X[combaero::species_index_from_name("AR")] = 0.0097;
  return X;
}
}  // namespace

// ---- centred crossflow segment -------------------------------------------

TEST(SideStreamMomentum, IsHalfOfBothAdjacentStationsWithTheirOwnDensities) {
  const auto X = air();
  const double m_arr = 0.015, m = 0.02, m_out = 0.026, A = 5e-4;
  const double Pa = 1.04e5, Ta = 320.0, Pm = 1.02e5, Tm = 330.0, Po = 1.0e5, To = 340.0;
  const auto r = side_stream_momentum_drop(m_arr, m_out, Pa, Ta, X, Po, To, X, A);
  const double ra = combaero::density(Ta, Pa, X), rm = combaero::density(Tm, Pm, X);
  const double ro = combaero::density(To, Po, X);
  const double up = (m * m / rm - m_arr * m_arr / ra) / (A * A);
  const double dn = (m_out * m_out / ro - m * m / rm) / (A * A);
  EXPECT_NEAR(r.dP, 0.5 * (up + dn), 1e-9 * r.dP);
}

TEST(SideStreamMomentum, DerivativesMatchCentralDifferences) {
  const auto X = air();
  const double A = 3e-4, e = 1e-8;
  const double Pa = 1.05e5, Ta = 310.0, Po = 1.0e5, To = 360.0;
  for (double m_arr : {-0.004, 0.003, 0.01}) {
    for (double m_out : {-0.002, 0.005, 0.02}) {
      const auto r = side_stream_momentum_drop(m_arr, m_out, Pa, Ta, X, Po, To, X, A);
      auto f = [&](double a, double b, double pa, double ta, double po, double to) {
        return side_stream_momentum_drop(a, b, pa, ta, X, po, to, X, A).dP;
      };
      const double hP = 1.0, hT = 1e-3;
      const double fo = (f(m_arr, m_out + e, Pa, Ta, Po, To) - f(m_arr, m_out - e, Pa, Ta, Po, To)) / (2 * e);
      const double fa = (f(m_arr + e, m_out, Pa, Ta, Po, To) - f(m_arr - e, m_out, Pa, Ta, Po, To)) / (2 * e);
      const double fPa = (f(m_arr, m_out, Pa + hP, Ta, Po, To) - f(m_arr, m_out, Pa - hP, Ta, Po, To)) / (2 * hP);
      const double fTa = (f(m_arr, m_out, Pa, Ta + hT, Po, To) - f(m_arr, m_out, Pa, Ta - hT, Po, To)) / (2 * hT);
      const double fPo = (f(m_arr, m_out, Pa, Ta, Po + hP, To) - f(m_arr, m_out, Pa, Ta, Po - hP, To)) / (2 * hP);
      const double fTo = (f(m_arr, m_out, Pa, Ta, Po, To + hT) - f(m_arr, m_out, Pa, Ta, Po, To - hT)) / (2 * hT);
      EXPECT_NEAR(r.d_dm_out, fo, 1e-6 * std::abs(fo));
      EXPECT_NEAR(r.d_dm_arr, fa, 1e-6 * std::abs(fa));
      EXPECT_NEAR(r.d_dP_arr, fPa, 1e-5 * std::abs(fPa) + 1e-12);
      EXPECT_NEAR(r.d_dT_arr, fTa, 1e-5 * std::abs(fTa) + 1e-12);
      EXPECT_NEAR(r.d_dP_out, fPo, 1e-5 * std::abs(fPo) + 1e-12);
      EXPECT_NEAR(r.d_dT_out, fTo, 1e-5 * std::abs(fTo) + 1e-12);
    }
  }
}

TEST(SideStreamMomentum, ZeroFlowHasZeroSlopeAndBadInputsThrow) {
  const auto X = air();
  const auto r = side_stream_momentum_drop(0.0, 0.0, 1e5, 300.0, X, 1e5, 300.0, X, 1e-3);
  EXPECT_EQ(r.d_dm_out, 0.0);
  EXPECT_EQ(r.d_dm_arr, 0.0);
  // Near, not equal: FMA contraction leaves a 1e-16 residue in a - a.
  EXPECT_NEAR(side_stream_momentum_drop(0.01, 0.01, 1e5, 300.0, X, 1e5, 300.0, X, 1e-3).dP,
              0.0, 1e-12);
  EXPECT_THROW(side_stream_momentum_drop(0.01, 0.02, 0.0, 300.0, X, 1e5, 300.0, X, 1e-3),
               std::invalid_argument);
  EXPECT_THROW(side_stream_momentum_drop(0.01, 0.02, 1e5, 300.0, X, 1e5, 300.0, X, 0.0),
               std::invalid_argument);
}

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
