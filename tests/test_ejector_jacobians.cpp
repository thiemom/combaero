// Jacobian-accuracy check for the ejector C++ port (include/ejector.h):
// compares the analytic d/dp_g, d/dt_g, d/dp_e, d/dt_e against the
// pre-baked central-difference targets (h = 1e-6*|x|) in
// validation/ejector/data/huang1999_reference_data.h, for all 31 rows.
// Mirrors tests/test_solver_jacobians.cpp's OrificeCompressibleDerivatives
// EXPECT_NEAR idiom.

#include "ejector.h"
#include "validation/ejector/data/huang1999_reference_data.h"
#include <cmath>
#include <gtest/gtest.h>

using namespace combaero::solver;
using combaero::validation::ejector::kHuangCases;

namespace {

// Relative tolerance for analytic-vs-CD comparison. The reference file's
// own CD step (h = 1e-6*|x|) has truncation/round-off noise around
// 1e-8-1e-7 relative on these expressions (confirmed in
// scripts/derive_ejector_jacobian.py's numeric cross-check); 1e-5 leaves
// comfortable margin while still catching a real derivation error.
constexpr double kRelTol = 1.0e-5;
constexpr double kAbsFloor = 1.0e-9;

void expect_jacobian_close(double analytic, double target, const char* what,
                           const char* case_id) {
  double tol = kRelTol * std::abs(target) + kAbsFloor;
  EXPECT_NEAR(analytic, target, tol) << what << " mismatch for case " << case_id
                                     << " (analytic=" << analytic << ", target=" << target << ")";
}

} // namespace

TEST(EjectorJacobians, EntrainmentRatioMatchesCentralDifferenceForAllCases) {
  for (const auto& c : kHuangCases) {
    EjectorGeometry geom{c.area_ratio_nozzle, c.area_ratio_mix};
    EjectorEntrainmentJacobian jac =
        ejector_entrainment_ratio_and_jacobian(c.p_g_pa, c.t_g_k, c.p_e_pa, c.t_e_k, geom, c.gamma);

    expect_jacobian_close(jac.domega_dp_g, c.domega_dp_g, "domega_dp_g", c.id);
    expect_jacobian_close(jac.domega_dt_g, c.domega_dt_g, "domega_dt_g", c.id);
    expect_jacobian_close(jac.domega_dp_e, c.domega_dp_e, "domega_dp_e", c.id);
    expect_jacobian_close(jac.domega_dt_e, c.domega_dt_e, "domega_dt_e", c.id);
  }
}

TEST(EjectorJacobians, CriticalBackPressureMatchesCentralDifferenceForAllCases) {
  for (const auto& c : kHuangCases) {
    EjectorGeometry geom{c.area_ratio_nozzle, c.area_ratio_mix};
    EjectorCriticalPressureJacobian jac = ejector_critical_back_pressure_and_jacobian(
        c.p_g_pa, c.t_g_k, c.p_e_pa, c.t_e_k, geom, c.gamma, c.r_gas, c.recovery_efficiency);

    expect_jacobian_close(jac.dpc_dp_g, c.dpc_dp_g, "dpc_dp_g", c.id);
    expect_jacobian_close(jac.dpc_dt_g, c.dpc_dt_g, "dpc_dt_g", c.id);
    expect_jacobian_close(jac.dpc_dp_e, c.dpc_dp_e, "dpc_dp_e", c.id);
    expect_jacobian_close(jac.dpc_dt_e, c.dpc_dt_e, "dpc_dt_e", c.id);
  }
}

TEST(EjectorJacobians, ChokedMassFlowMatchesCentralDifference) {
  // huang1999_reference_data.h has no column for choked_mass_flow (the
  // reference model never needed the absolute mass flow, only the ratio
  // form) -- CD spot-check computed directly here instead.
  double gamma = 1.3;
  double r_gas = 287.0;
  double area_throat = 5.0e-5;
  double eta = 0.95;

  for (auto [p0, t0] : {std::pair{4.0e5, 400.0}, std::pair{1.2e5, 320.0}}) {
    auto [mdot, d_dp0, d_dt0] =
        ejector_choked_mass_flow_and_jacobian(p0, t0, area_throat, gamma, r_gas, eta);

    double eps_p = 1.0e-6 * p0;
    double fd_dp0 = (ejector_choked_mass_flow(p0 + eps_p, t0, area_throat, gamma, r_gas, eta) -
                     ejector_choked_mass_flow(p0 - eps_p, t0, area_throat, gamma, r_gas, eta)) /
                    (2.0 * eps_p);
    double eps_t = 1.0e-6 * t0;
    double fd_dt0 = (ejector_choked_mass_flow(p0, t0 + eps_t, area_throat, gamma, r_gas, eta) -
                     ejector_choked_mass_flow(p0, t0 - eps_t, area_throat, gamma, r_gas, eta)) /
                    (2.0 * eps_t);

    EXPECT_NEAR(d_dp0, fd_dp0, std::abs(fd_dp0) * 1e-8 + 1e-12);
    EXPECT_NEAR(d_dt0, fd_dt0, std::abs(fd_dt0) * 1e-8 + 1e-12);
    EXPECT_GT(mdot, 0.0);
  }
}

TEST(EjectorJacobians, CDNozzleAnalyticMatchesCentralDifference) {
  // C-D nozzle R0 residual: analytic (p0, t0, p_py) partials vs central
  // difference of the value function. Covers unchoked, near-threshold, and
  // choked (flat) points.
  const double A_t = 3.14e-5, A_e = 1.0e-4, g = 1.40, R = 287.0, eta = 0.95, eps = 1e-3;
  struct Pt {
    double p0, t0, p_py;
  };
  const Pt pts[] = {
      {101325.0, 300.0, 100207.0}, // unchoked (reported root)
      {101325.0, 300.0, 98700.0},  // near the choke threshold
      {101325.0, 300.0, 60000.0},  // choked (flat region)
      {604000.0, 368.0, 300000.0}, // higher motive, unchoked
  };
  auto val = [&](double p0, double t0, double p_py) {
    return ejector_cd_nozzle_mass_flow(p0, t0, p_py, A_t, A_e, g, R, eta, eps);
  };
  for (const auto& p : pts) {
    EjectorCDNozzleJacobian j =
        ejector_cd_nozzle_mass_flow_and_jacobian(p.p0, p.t0, p.p_py, A_t, A_e, g, R, eta, eps);
    double hp0 = 1e-6 * p.p0, ht0 = 1e-6 * p.t0, hpy = 1e-6 * p.p_py;
    double cd_p0 = (val(p.p0 + hp0, p.t0, p.p_py) - val(p.p0 - hp0, p.t0, p.p_py)) / (2 * hp0);
    double cd_t0 = (val(p.p0, p.t0 + ht0, p.p_py) - val(p.p0, p.t0 - ht0, p.p_py)) / (2 * ht0);
    double cd_py = (val(p.p0, p.t0, p.p_py + hpy) - val(p.p0, p.t0, p.p_py - hpy)) / (2 * hpy);
    EXPECT_NEAR(j.dmdot_dp0, cd_p0, 1e-6 * std::abs(cd_p0) + 1e-12) << "dp0 at p_py=" << p.p_py;
    EXPECT_NEAR(j.dmdot_dt0, cd_t0, 1e-6 * std::abs(cd_t0) + 1e-12) << "dt0 at p_py=" << p.p_py;
    EXPECT_NEAR(j.dmdot_dp_py, cd_py, 1e-5 * std::abs(cd_py) + 1e-12) << "dp_py at p_py=" << p.p_py;
  }
}

TEST(EjectorJacobians, JetPumpDischargeValueAndJacobian) {
  // R3 building block: value parity vs the Python mixing closure, and analytic
  // partials (p_g, t_g, p_e, t_e, p_py, omega) vs central difference.
  struct Case {
    double p_g, t_g, p_e, t_e, p_py, omega, gamma, rec, expected;
  };
  const Case cs[] = {
      {101325.0, 300.0, 101325.0, 300.0, 100207.0, 7.0, 1.4, 1.0, 101325.00000000015},
      {102435.0, 300.0, 101325.0, 288.15, 100000.0, 3.5, 1.4, 1.0, 101544.85636197214},
      {110000.0, 320.0, 100000.0, 300.0, 95000.0, 1.5, 1.4, 0.9, 92957.65455779647},
  };
  auto val = [](double pg, double tg, double pe, double te, double ppy, double om, double gm,
                double rec) {
    return ejector_jetpump_discharge_and_jacobian(pg, tg, pe, te, ppy, om, gm, rec).p03;
  };
  for (const auto& c : cs) {
    auto j = ejector_jetpump_discharge_and_jacobian(c.p_g, c.t_g, c.p_e, c.t_e, c.p_py, c.omega,
                                                    c.gamma, c.rec);
    EXPECT_NEAR(j.p03, c.expected, 1e-7 * c.expected) << "value at p_py=" << c.p_py;
    double h[6] = {1e-6 * c.p_g, 1e-6 * c.t_g, 1e-6 * c.p_e,
                   1e-6 * c.t_e, 1e-6 * c.p_py, 1e-6 * c.omega};
    double an[6] = {j.dp03_dp_g, j.dp03_dt_g, j.dp03_dp_e, j.dp03_dt_e, j.dp03_dp_py, j.dp03_domega};
    for (int k = 0; k < 6; ++k) {
      double a[8] = {c.p_g, c.t_g, c.p_e, c.t_e, c.p_py, c.omega, c.gamma, c.rec};
      double b[8];
      for (int m = 0; m < 8; ++m) b[m] = a[m];
      a[k] += h[k];
      b[k] -= h[k];
      double cd = (val(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]) -
                   val(b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7])) /
                  (2 * h[k]);
      // Abs floor 1e-5 tolerates central-difference round-off on partials that
      // are genuinely ~0 at the symmetric degenerate point (p_g=p_e, t_g=t_e);
      // real partials here are O(0.01)-O(1000), caught by the relative term.
      EXPECT_NEAR(an[k], cd, 1e-5 * std::abs(cd) + 1e-5) << "seed " << k << " at p_py=" << c.p_py;
    }
  }
}
