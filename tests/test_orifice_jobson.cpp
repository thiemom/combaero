#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <vector>

#include "composition.h"
#include "math_constants.h"
#include "orifice.h"
#include "solver_interface.h"
#include "thermo.h"

// Jobson (1955) compressible orifice and the ISO 5167-2 Eq. (5) kernel (#854).

namespace {

std::vector<double> air_Y() {
  std::vector<double> X(combaero::num_species(), 0.0);
  X[combaero::species_index_from_name("N2")] = 0.7808;
  X[combaero::species_index_from_name("O2")] = 0.2095;
  X[combaero::species_index_from_name("AR")] = 0.0097;
  return combaero::mole_to_mass(X);
}

// Jobson's Table 2: the slit, Ci = pi/(pi+2), n = 1.4, as printed.
struct Row {
  double r;
  double C;
};
const Row kTable2[] = {{0.932, 0.622}, {0.865, 0.635}, {0.805, 0.648}, {0.745, 0.663},
                       {0.690, 0.679}, {0.640, 0.696}, {0.590, 0.714}, {0.545, 0.734},
                       {0.528, 0.741}};

}  // namespace

TEST(Jobson1955, ReproducesTable2WithinItsPrintedDigits) {
  const double Ci = M_PI / (M_PI + 2.0);
  EXPECT_NEAR(jobson_contraction_coefficient(1.0, Ci, 1.4), 0.611, 5e-4);
  for (const Row& row : kTable2) {
    EXPECT_NEAR(jobson_contraction_coefficient(row.r, Ci, 1.4), row.C, 1e-3) << row.r;
  }
}

TEST(Jobson1955, ForceDefectIsEq5AtBetaZeroAndBordaIsZero) {
  for (double Ci : {0.5, 0.611, 0.7, 1.0}) {
    EXPECT_DOUBLE_EQ(jobson_force_defect(Ci), 1.0 / Ci - 1.0 / (2.0 * Ci * Ci));
  }
  EXPECT_DOUBLE_EQ(jobson_force_defect(0.5), 0.0);
}

TEST(Jobson1955, KernelFlowIsC1ThroughTheCriticalRatio) {
  // Both branches meet at r_c in value and in d(mdot)/dPb (dG/dp = 0 and
  // dB/dp = -1/Gc on both sides), so the solver sees no kink at choke.
  const auto Y = air_Y();
  const double Pt = 3.0e5, Tt = 400.0, Cd = 0.61, A = 1e-4;
  const auto X = combaero::mass_to_mole(Y);
  const double n = combaero::isentropic_expansion_coefficient(Tt, X);
  const double rc = std::pow(2.0 / (n + 1.0), n / (n - 1.0));
  for (double beta : {0.0, 0.5}) {
    const double P1 = beta > 0.0 ? 0.995 * Pt : Pt;
    const double h = 1e-7 * Pt;
    const auto sub = combaero::solver::orifice_jobson_residuals_and_jacobian(
        0.0, Pt, P1, Tt, Y, rc * Pt + h, Cd, A, beta, 0.98);
    const auto sup = combaero::solver::orifice_jobson_residuals_and_jacobian(
        0.0, Pt, P1, Tt, Y, rc * Pt - h, Cd, A, beta, 0.98);
    EXPECT_FALSE(sub.choked);
    EXPECT_TRUE(sup.choked);
    // No jump: the 2h step accounts for the whole difference (to O(h^2)).
    EXPECT_NEAR(sub.m_dot_calc - sup.m_dot_calc, 2.0 * h * sup.d_mdot_dP_static_down,
                1e-11 * sub.m_dot_calc);
    EXPECT_NEAR(sub.d_mdot_dP_static_down, sup.d_mdot_dP_static_down,
                1e-4 * std::abs(sup.d_mdot_dP_static_down));
    // Choked, the jet keeps widening: the flow still rises as Pb falls.
    EXPECT_LT(sup.d_mdot_dP_static_down, 0.0);
  }
}

TEST(Jobson1955, KernelJacobianMatchesCentralDifferences) {
  const auto Y = air_Y();
  const double A = 1e-4;
  using Fn = std::function<combaero::solver::OrificeModelResult(double, double, double, double,
                                                                double, double)>;
  const Fn jobson = [&](double Pt, double P1, double Tt, double Pb, double Cd, double beta) {
    return combaero::solver::orifice_jobson_residuals_and_jacobian(0.0, Pt, P1, Tt, Y, Pb, Cd,
                                                                   A, beta, 0.98);
  };
  const Fn iso = [&](double Pt, double P1, double Tt, double Pb, double Cd, double beta) {
    return combaero::solver::orifice_iso5167_residuals_and_jacobian(0.0, Pt, P1, Tt, Y, Pb, Cd,
                                                                    A, beta);
  };
  for (const Fn* fn : {&jobson, &iso}) {
    for (double beta : {0.0, 0.6}) {
      // Subcritical, near choke on either side, deep choke.
      for (double Pb : {1.8e5, 1.07e5, 1.04e5, 0.6e5}) {
        double x[5] = {2.0e5, beta > 0.0 ? 1.98e5 : 2.0e5, 600.0, Pb, 0.61};
        const auto r = (*fn)(x[0], x[1], x[2], x[3], x[4], beta);
        const double an[5] = {r.d_mdot_dP_total_up, r.d_mdot_dP_static_up, r.d_mdot_dTt_up,
                              r.d_mdot_dP_static_down, r.d_mdot_dCd};
        for (int k = 0; k < 5; ++k) {
          if (beta == 0.0 && k == 1) {
            EXPECT_EQ(an[k], 0.0);  // a reservoir does not read the static
            continue;
          }
          double xp[5], xm[5];
          for (int j = 0; j < 5; ++j) xp[j] = xm[j] = x[j];
          const double h = 1e-6 * x[k];
          xp[k] += h;
          xm[k] -= h;
          const double fd = ((*fn)(xp[0], xp[1], xp[2], xp[3], xp[4], beta).m_dot_calc -
                             (*fn)(xm[0], xm[1], xm[2], xm[3], xm[4], beta).m_dot_calc) /
                            (2.0 * h);
          const double scale = std::max(std::abs(fd), 1e-6 * r.m_dot_calc / x[k]);
          EXPECT_NEAR(an[k], fd, 1e-6 * scale) << "beta " << beta << " Pb " << Pb << " k " << k;
        }
      }
    }
  }
}

TEST(Jobson1955, IncompressibleLimitIsTheElementsOwnFlow) {
  // dp -> 0: mdot -> Cd E A sqrt(2 rho1 dp), whatever beta and Cv (Cv only
  // splits Cd into Cv Ci), so the model changes no incompressible answer.
  const auto Y = air_Y();
  for (double beta : {0.0, 0.5, 0.75}) {
    for (double Cv : {1.0, 0.982}) {
      // The approach static consistent with the flow it carries:
      // Pt - P1 = rho V1^2 / 2, V1 = mdot beta^2 / (rho A).
      const double Pt = 2.0e5, Tt = 300.0, A = 1e-4;
      const double rho = combaero::density(Tt, Pt, combaero::mass_to_mole(Y));
      const double b4 = beta * beta * beta * beta;
      double P1 = Pt;
      combaero::solver::OrificeModelResult r;
      for (int it = 0; it < 60; ++it) {
        r = combaero::solver::orifice_jobson_residuals_and_jacobian(0.0, Pt, P1, Tt, Y,
                                                                    P1 - 2.0, 0.6, A, beta, Cv);
        P1 = Pt - r.m_dot_calc * r.m_dot_calc * b4 / (2.0 * rho * A * A);
      }
      EXPECT_NEAR(r.expansibility, 1.0, 1e-4) << beta << " " << Cv;
    }
  }
}

TEST(Iso5167Kernel, IsEq1WithEq5AndFlagsBelowItsRange) {
  const auto Y = air_Y();
  const auto X = combaero::mass_to_mole(Y);
  const double Pt = 2.0e5, Tt = 300.0, beta = 0.5, A = 1e-3, Cd = 0.605;
  const double kappa = combaero::isentropic_expansion_coefficient(Tt, X);
  for (double r : {0.97, 0.8, 0.75, 0.6}) {
    const double p2 = r * Pt;
    // A plenum-fed plate here (beta 0.5 with P1 = Pt): rho1 at the total state.
    const auto res = combaero::solver::orifice_iso5167_residuals_and_jacobian(
        0.0, Pt, Pt, Tt, Y, p2, Cd, A, beta);
    const double eps = expansibility_factor(beta, Pt - p2, Pt, kappa);
    const double rho1 = combaero::density(Tt, Pt, X);
    const double E = velocity_of_approach_factor(beta);
    EXPECT_NEAR(res.expansibility, eps, 1e-12);
    EXPECT_NEAR(res.m_dot_calc, Cd * E * eps * A * std::sqrt(2.0 * rho1 * (Pt - p2)),
                1e-5 * res.m_dot_calc);  // the kernel's 1e-4 density clamp
    EXPECT_EQ(res.in_range, r >= orifice::iso5167_expansibility::min_pressure_ratio);
  }
}
