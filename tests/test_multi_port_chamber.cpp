#include "multi_port_chamber.h"
#include "math_constants.h"
#include "solver_interface.h"
#include "composition.h"
#include "thermo.h"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

using combaero::border_carnot_L;
using combaero::dborder_carnot_L_ddelta;
using combaero::BC_LOSS_PREFACTOR;
using combaero::HAGER_FRACTION;
using combaero::solver::border_carnot_loss_residual_and_jacobian;
using combaero::solver::multi_port_chamber_residuals_and_jacobian;

namespace {

// Build a dry-air mass-fraction vector for tests. Composition layout matches
// combaero::species: indices look up by name to stay robust under reordering.
std::vector<double> dry_air_Y() {
  std::vector<double> Y(combaero::num_species(), 0.0);
  Y[combaero::species_index_from_name("N2")] = 0.767;
  Y[combaero::species_index_from_name("O2")] = 0.233;
  return Y;
}

// 1D central-difference helper for scalar functions of double.
template <typename Fn>
double fd_central(Fn f, double x, double h) {
  return (f(x + h) - f(x - h)) / (2.0 * h);
}

}  // namespace

// =============================================================================
// border_carnot_L kernel
// =============================================================================

TEST(BorderCarnotKernel, ZeroAngleGivesZeroLoss) {
  EXPECT_DOUBLE_EQ(border_carnot_L(0.0), 0.0);
}

TEST(BorderCarnotKernel, NinetyDegreeMatchesHagerSharpEdge) {
  // delta_geom = pi/2 -> theta_eff = (3/4) * pi/2 = 3*pi/8 (67.5 deg)
  // L = 4 * (1 - cos(67.5 deg))^2.
  const double delta = M_PI / 2.0;
  const double theta_eff = HAGER_FRACTION * delta;
  const double expected = BC_LOSS_PREFACTOR
                          * (1.0 - std::cos(theta_eff))
                          * (1.0 - std::cos(theta_eff));
  EXPECT_NEAR(border_carnot_L(delta), expected, 1e-14);
  // Sanity bound: this should be O(1.4) for the sharp-edged 90 deg lateral.
  EXPECT_NEAR(border_carnot_L(delta), 1.524, 1e-3);
}

TEST(BorderCarnotKernel, AnalyticDerivativeMatchesFD) {
  // Scan delta_geom across the lateral-branch range.
  const double h = 1e-7;
  for (double delta = 0.05; delta <= M_PI / 2.0; delta += 0.1) {
    const double analytic = dborder_carnot_L_ddelta(delta);
    const double fd = fd_central(
        [](double x) { return border_carnot_L(x); }, delta, h);
    EXPECT_NEAR(analytic, fd, 1e-7)
        << "delta_geom=" << delta;
  }
}

// =============================================================================
// multi_port_chamber: residual identities at the canonical sign convention
// =============================================================================

TEST(MultiPortChamber, MassResidualEqualsSignedSum) {
  // 3 ports, mass-balanced: +0.1 in, two symmetric outflows of -0.05 each.
  // Use a lateral port (theta != 0) so the impulse term is non-trivial under
  // the sin^2(theta) projection.
  const auto Y = dry_air_Y();
  const std::vector<double> P{1.0e5, 1.0e5, 1.0e5};
  const std::vector<double> mdot{+0.1, -0.05, -0.05};
  const std::vector<double> T{300.0, 300.0, 300.0};
  const std::vector<std::vector<double>> Yall{Y, Y, Y};
  const std::vector<double> A{1e-3, 1e-3, 1e-3};

  // All lateral (theta=pi/2 -> sin^2=1) for this test of the impulse scaling.
  const std::vector<double> theta(3, M_PI / 2.0);
  auto res = multi_port_chamber_residuals_and_jacobian(
      1.0e5, P, mdot, T, Yall, A, theta);
  EXPECT_NEAR(res.mass_residual, 0.0, 1e-14);

  // Symmetric outflow ports (same |mdot|, T, P, A) -> identical impulse residuals.
  // This is the sign-free property: u^2 doesn't care about flow direction.
  EXPECT_NEAR(res.impulse_residuals[1], res.impulse_residuals[2], 1e-12)
      << "symmetric ports (|mdot|=0.05) should give identical impulse residuals";
  // Larger |mdot| -> larger dynamic head -> larger residual at uniform P_jct=P_i.
  EXPECT_GT(res.impulse_residuals[0], res.impulse_residuals[1])
      << "larger |mdot| should give larger dynamic head";
}

TEST(MultiPortChamber, SignFreeInMdot) {
  // Same |mdot|, opposite sign at otherwise-identical ports gives same impulse residual.
  const auto Y = dry_air_Y();
  const std::vector<double> theta(2, 0.0);
  auto r_plus = multi_port_chamber_residuals_and_jacobian(
      1.0e5, {1.0e5, 1.0e5}, {+0.1, -0.1}, {300.0, 300.0}, {Y, Y}, {1e-3, 1e-3}, theta);
  auto r_minus = multi_port_chamber_residuals_and_jacobian(
      1.0e5, {1.0e5, 1.0e5}, {-0.1, +0.1}, {300.0, 300.0}, {Y, Y}, {1e-3, 1e-3}, theta);
  EXPECT_NEAR(r_plus.impulse_residuals[0], r_minus.impulse_residuals[0], 1e-12);
  EXPECT_NEAR(r_plus.impulse_residuals[1], r_minus.impulse_residuals[1], 1e-12);
  // Mass residual flips sign with the flow reversal.
  EXPECT_NEAR(r_plus.mass_residual, -r_minus.mass_residual, 1e-14);
}

TEST(MultiPortChamber, RejectsSizeMismatch) {
  const auto Y = dry_air_Y();
  EXPECT_THROW(multi_port_chamber_residuals_and_jacobian(
                   1.0e5, {1.0e5, 1.0e5}, {0.1, -0.05, -0.05},
                   {300.0, 300.0}, {Y, Y}, {1e-3, 1e-3}, {0.0, 0.0}),
               std::invalid_argument);
}

TEST(MultiPortChamber, RejectsTooFewPorts) {
  const auto Y = dry_air_Y();
  EXPECT_THROW(multi_port_chamber_residuals_and_jacobian(
                   1.0e5, {1.0e5}, {0.0}, {300.0}, {Y}, {1e-3}, {0.0}),
               std::invalid_argument);
}

TEST(MultiPortChamber, FourPortManifoldScalesWithN) {
  // N=4: 1 inflow, 3 equal outflows. Mass balance closes.
  const auto Y = dry_air_Y();
  const std::vector<double> P(4, 1.0e5);
  const std::vector<double> mdot{+0.15, -0.05, -0.05, -0.05};
  const std::vector<double> T(4, 300.0);
  const std::vector<std::vector<double>> Yall(4, Y);
  const std::vector<double> A(4, 1e-3);

  const std::vector<double> theta(4, 0.0);
  auto res = multi_port_chamber_residuals_and_jacobian(
      1.0e5, P, mdot, T, Yall, A, theta);
  EXPECT_EQ(res.impulse_residuals.size(), 4u);
  EXPECT_EQ(res.port_jac.size(), 4u);
  EXPECT_NEAR(res.mass_residual, 0.0, 1e-14);

  // Outflow ports (1, 2, 3) have identical mdot/T/P -> identical impulse residuals.
  EXPECT_NEAR(res.impulse_residuals[1], res.impulse_residuals[2], 1e-12);
  EXPECT_NEAR(res.impulse_residuals[2], res.impulse_residuals[3], 1e-12);
}

// =============================================================================
// multi_port_chamber: analytic vs FD Jacobians (per port)
// =============================================================================

TEST(MultiPortChamber, AnalyticJacobianMatchesFD) {
  const auto Y = dry_air_Y();
  const std::vector<double> P{1.0e5, 1.005e5, 0.995e5};
  const std::vector<double> mdot{+0.1, -0.05, -0.05};
  const std::vector<double> T{300.0, 310.0, 305.0};
  const std::vector<std::vector<double>> Yall{Y, Y, Y};
  const std::vector<double> A{1e-3, 1.2e-3, 0.8e-3};
  const double P_jct = 1.0e5;

  const std::vector<double> theta(3, 0.0);
  auto base = multi_port_chamber_residuals_and_jacobian(
      P_jct, P, mdot, T, Yall, A, theta);

  for (std::size_t k = 0; k < 3; ++k) {
    // dR_k/dP_k
    {
      auto P_p = P;
      P_p[k] += 1.0;
      auto P_m = P;
      P_m[k] -= 1.0;
      auto rp = multi_port_chamber_residuals_and_jacobian(P_jct, P_p, mdot, T, Yall, A, theta);
      auto rm = multi_port_chamber_residuals_and_jacobian(P_jct, P_m, mdot, T, Yall, A, theta);
      const double fd = (rp.impulse_residuals[k] - rm.impulse_residuals[k]) / 2.0;
      EXPECT_NEAR(fd, base.port_jac[k].dR_dP, 1e-6)
          << "port " << k << " dR/dP";
    }
    // dR_k/dmdot_k
    {
      auto md_p = mdot;
      md_p[k] += 1e-5;
      auto md_m = mdot;
      md_m[k] -= 1e-5;
      auto rp = multi_port_chamber_residuals_and_jacobian(P_jct, P, md_p, T, Yall, A, theta);
      auto rm = multi_port_chamber_residuals_and_jacobian(P_jct, P, md_m, T, Yall, A, theta);
      const double fd = (rp.impulse_residuals[k] - rm.impulse_residuals[k]) / 2e-5;
      const double rel = std::abs(fd - base.port_jac[k].dR_dmdot)
                         / std::max(std::abs(base.port_jac[k].dR_dmdot), 1e-12);
      EXPECT_LT(rel, 1e-6) << "port " << k << " dR/dmdot rel-error";
    }
    // dR_k/dT_k
    {
      auto T_p = T;
      T_p[k] += 0.01;
      auto T_m = T;
      T_m[k] -= 0.01;
      auto rp = multi_port_chamber_residuals_and_jacobian(P_jct, P, mdot, T_p, Yall, A, theta);
      auto rm = multi_port_chamber_residuals_and_jacobian(P_jct, P, mdot, T_m, Yall, A, theta);
      const double fd = (rp.impulse_residuals[k] - rm.impulse_residuals[k]) / 0.02;
      const double rel = std::abs(fd - base.port_jac[k].dR_dT)
                         / std::max(std::abs(base.port_jac[k].dR_dT), 1e-12);
      EXPECT_LT(rel, 1e-6) << "port " << k << " dR/dT rel-error";
    }
  }
}

// -----------------------------------------------------------------------------
// MPCE Jacobian: cross-port terms (axial-reference cross-coupling) + dR/dP_jct
// + mass row. theta != 0 so the cross-coupling formula is actually exercised.
// Tight tolerance (1e-6) so any future MPCE-v1 regression trips immediately.
// -----------------------------------------------------------------------------
TEST(MultiPortChamber, JacobianFullCoverageAtNonzeroTheta) {
  const auto Y = dry_air_Y();
  const std::vector<double> P{1.0e5, 1.0050e5, 0.9950e5};
  const std::vector<double> mdot{+0.10, -0.06, -0.04};
  const std::vector<double> T{300.0, 308.0, 305.0};
  const std::vector<std::vector<double>> Yall{Y, Y, Y};
  const std::vector<double> A{1.0e-3, 1.2e-3, 0.8e-3};
  const double P_jct = 1.0e5;
  // theta_0 is the axial reference (always 0); lateral port 2 at 45 deg
  // exercises the cross-coupling terms that the diagonal-only test misses.
  const std::vector<double> theta{0.0, 0.0, M_PI / 4.0};

  auto base = multi_port_chamber_residuals_and_jacobian(
      P_jct, P, mdot, T, Yall, A, theta);

  // dR_k/dP_jct: should be -1 by construction.
  for (std::size_t k = 0; k < 3; ++k) {
    auto Pj_p = multi_port_chamber_residuals_and_jacobian(
        P_jct + 1.0, P, mdot, T, Yall, A, theta);
    auto Pj_m = multi_port_chamber_residuals_and_jacobian(
        P_jct - 1.0, P, mdot, T, Yall, A, theta);
    const double fd =
        (Pj_p.impulse_residuals[k] - Pj_m.impulse_residuals[k]) / 2.0;
    EXPECT_NEAR(fd, -1.0, 1e-9) << "port " << k << " dR/dP_jct";
  }

  // Cross-port dR_k/dP_axial (where axial = port 0). For port 0 itself this
  // collapses to the diagonal; for ports 1, 2 it exercises the cross-coupling
  // contribution via cross_dR_dP_axial[k].
  for (std::size_t k = 1; k < 3; ++k) {
    auto P_p = P;
    P_p[0] += 1.0;
    auto P_m = P;
    P_m[0] -= 1.0;
    auto rp = multi_port_chamber_residuals_and_jacobian(P_jct, P_p, mdot, T, Yall, A, theta);
    auto rm = multi_port_chamber_residuals_and_jacobian(P_jct, P_m, mdot, T, Yall, A, theta);
    const double fd = (rp.impulse_residuals[k] - rm.impulse_residuals[k]) / 2.0;
    const double expected = base.cross_dR_dP_axial[k];
    EXPECT_NEAR(fd, expected, 1e-6)
        << "cross dR_" << k << "/dP_0 disagrees with FD";
  }

  // Cross-port dR_k/dmdot_0.
  for (std::size_t k = 1; k < 3; ++k) {
    auto md_p = mdot;
    md_p[0] += 1e-5;
    auto md_m = mdot;
    md_m[0] -= 1e-5;
    auto rp = multi_port_chamber_residuals_and_jacobian(P_jct, P, md_p, T, Yall, A, theta);
    auto rm = multi_port_chamber_residuals_and_jacobian(P_jct, P, md_m, T, Yall, A, theta);
    const double fd = (rp.impulse_residuals[k] - rm.impulse_residuals[k]) / 2e-5;
    const double expected = base.cross_dR_dmdot_axial[k];
    const double rel =
        std::abs(fd - expected) / std::max(std::abs(expected), 1e-12);
    EXPECT_LT(rel, 1e-5)
        << "cross dR_" << k << "/dmdot_0 disagrees with FD (rel=" << rel
        << ", fd=" << fd << ", analytic=" << expected << ")";
  }

  // Cross-port dR_k/dT_0.
  for (std::size_t k = 1; k < 3; ++k) {
    auto T_p = T;
    T_p[0] += 0.01;
    auto T_m = T;
    T_m[0] -= 0.01;
    auto rp = multi_port_chamber_residuals_and_jacobian(P_jct, P, mdot, T_p, Yall, A, theta);
    auto rm = multi_port_chamber_residuals_and_jacobian(P_jct, P, mdot, T_m, Yall, A, theta);
    const double fd = (rp.impulse_residuals[k] - rm.impulse_residuals[k]) / 0.02;
    const double expected = base.cross_dR_dT_axial[k];
    if (std::abs(expected) < 1e-12 && std::abs(fd) < 1e-6) continue;
    const double rel =
        std::abs(fd - expected) / std::max(std::abs(expected), 1e-12);
    EXPECT_LT(rel, 1e-5)
        << "cross dR_" << k << "/dT_0 disagrees with FD (rel=" << rel
        << ", fd=" << fd << ", analytic=" << expected << ")";
  }

  // Mass row: R_mass = sum_i mdot_i, so dR_mass/dmdot_k = 1 exactly.
  for (std::size_t k = 0; k < 3; ++k) {
    auto md_p = mdot;
    md_p[k] += 1e-5;
    auto md_m = mdot;
    md_m[k] -= 1e-5;
    auto rp = multi_port_chamber_residuals_and_jacobian(P_jct, P, md_p, T, Yall, A, theta);
    auto rm = multi_port_chamber_residuals_and_jacobian(P_jct, P, md_m, T, Yall, A, theta);
    const double fd = (rp.mass_residual - rm.mass_residual) / 2e-5;
    EXPECT_NEAR(fd, 1.0, 1e-9) << "dR_mass/dmdot_" << k;
  }
}

// =============================================================================
// border_carnot_loss_residual_and_jacobian: analytic vs FD
// =============================================================================

TEST(BorderCarnotLossElement, AnalyticJacobianMatchesFD) {
  const auto Y = dry_air_Y();
  const double mdot = 0.1;
  const double Pt_in = 1.05e5;
  const double Pt_out = 1.04e5;
  const double P_in = 1.0e5;
  const double T_in = 300.0;
  const double area = 1e-3;
  const double delta = M_PI / 2.0;  // sharp-edged 90 deg lateral

  auto base = border_carnot_loss_residual_and_jacobian(
      mdot, Pt_in, Pt_out, P_in, T_in, Y, area, delta);

  // Sanity: residual should be negative (loss > applied stagnation drop).
  EXPECT_LT(base.residual, 0.0);

  // d/dPt_in == +1, d/dPt_out == -1 (exact).
  EXPECT_DOUBLE_EQ(base.d_res_dPt_in, 1.0);
  EXPECT_DOUBLE_EQ(base.d_res_dPt_out, -1.0);

  // d/dmdot: FD
  {
    const double h = 1e-5;
    auto rp = border_carnot_loss_residual_and_jacobian(
        mdot + h, Pt_in, Pt_out, P_in, T_in, Y, area, delta);
    auto rm = border_carnot_loss_residual_and_jacobian(
        mdot - h, Pt_in, Pt_out, P_in, T_in, Y, area, delta);
    const double fd = (rp.residual - rm.residual) / (2.0 * h);
    const double rel = std::abs(fd - base.d_res_dmdot)
                       / std::max(std::abs(base.d_res_dmdot), 1e-12);
    EXPECT_LT(rel, 1e-6) << "d/dmdot rel-error";
  }
  // d/dP_in: FD
  {
    const double h = 10.0;
    auto rp = border_carnot_loss_residual_and_jacobian(
        mdot, Pt_in, Pt_out, P_in + h, T_in, Y, area, delta);
    auto rm = border_carnot_loss_residual_and_jacobian(
        mdot, Pt_in, Pt_out, P_in - h, T_in, Y, area, delta);
    const double fd = (rp.residual - rm.residual) / (2.0 * h);
    const double rel = std::abs(fd - base.d_res_dP_in)
                       / std::max(std::abs(base.d_res_dP_in), 1e-12);
    EXPECT_LT(rel, 1e-6) << "d/dP_in rel-error";
  }
  // d/dT_in: FD
  {
    const double h = 0.01;
    auto rp = border_carnot_loss_residual_and_jacobian(
        mdot, Pt_in, Pt_out, P_in, T_in + h, Y, area, delta);
    auto rm = border_carnot_loss_residual_and_jacobian(
        mdot, Pt_in, Pt_out, P_in, T_in - h, Y, area, delta);
    const double fd = (rp.residual - rm.residual) / (2.0 * h);
    const double rel = std::abs(fd - base.d_res_dT_in)
                       / std::max(std::abs(base.d_res_dT_in), 1e-12);
    EXPECT_LT(rel, 1e-6) << "d/dT_in rel-error";
  }
}

TEST(BorderCarnotLossElement, SignFreeInMdot) {
  // Reversing mdot sign should leave the residual unchanged (u^2 is sign-free).
  const auto Y = dry_air_Y();
  const double Pt_in = 1.05e5;
  const double Pt_out = 1.04e5;
  const double P_in = 1.0e5;
  const double T_in = 300.0;
  const double area = 1e-3;
  const double delta = M_PI / 2.0;
  auto rp = border_carnot_loss_residual_and_jacobian(
      +0.1, Pt_in, Pt_out, P_in, T_in, Y, area, delta);
  auto rm = border_carnot_loss_residual_and_jacobian(
      -0.1, Pt_in, Pt_out, P_in, T_in, Y, area, delta);
  EXPECT_DOUBLE_EQ(rp.residual, rm.residual);
  // Jacobian d/dmdot DOES flip sign.
  EXPECT_NEAR(rp.d_res_dmdot, -rm.d_res_dmdot, 1e-9);
}
