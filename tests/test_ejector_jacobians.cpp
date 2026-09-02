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

TEST(EjectorJacobians, CDNozzleExitStaticImplicitDerivatives) {
  // Inverse of the C-D nozzle: analytic implicit-function partials
  // (m_dot, p0, t0) vs central difference of the bisection value. Also checks
  // the round-trip cd_nozzle(exit_static(m_dot)) == m_dot.
  const double A_t = 3.14e-5, A_e = 1.0e-4, g = 1.40, R = 287.0, eta = 0.95, eps = 1e-3;
  struct Pt {
    double m_dot, p0, t0;
  };
  const Pt pts[] = {
      {0.0051, 101325.0, 300.0}, // reported
      {0.0030, 101325.0, 300.0}, // deeper unchoked
      {0.0065, 101325.0, 300.0}, // near choke
      {0.0300, 604000.0, 368.0}, // higher motive
  };
  auto val = [&](double md, double p0, double t0) {
    return ejector_cd_nozzle_exit_static_and_jacobian(md, p0, t0, A_t, A_e, g, R, eta, eps).p_py;
  };
  for (const auto& p : pts) {
    auto j = ejector_cd_nozzle_exit_static_and_jacobian(p.m_dot, p.p0, p.t0, A_t, A_e, g, R, eta,
                                                        eps);
    // round-trip
    double back = ejector_cd_nozzle_mass_flow(p.p0, p.t0, j.p_py, A_t, A_e, g, R, eta, eps);
    EXPECT_NEAR(back, p.m_dot, 1e-9 * p.m_dot) << "round-trip at m_dot=" << p.m_dot;
    // implicit derivatives vs FD
    double hm = 1e-7 * p.m_dot, hp = 1e-6 * p.p0, ht = 1e-6 * p.t0;
    double cd_m = (val(p.m_dot + hm, p.p0, p.t0) - val(p.m_dot - hm, p.p0, p.t0)) / (2 * hm);
    double cd_p = (val(p.m_dot, p.p0 + hp, p.t0) - val(p.m_dot, p.p0 - hp, p.t0)) / (2 * hp);
    double cd_t = (val(p.m_dot, p.p0, p.t0 + ht) - val(p.m_dot, p.p0, p.t0 - ht)) / (2 * ht);
    EXPECT_NEAR(j.dp_py_dm_dot, cd_m, 1e-4 * std::abs(cd_m) + 1e-6) << "dm at " << p.m_dot;
    EXPECT_NEAR(j.dp_py_dp0, cd_p, 1e-4 * std::abs(cd_p) + 1e-9) << "dp0 at " << p.m_dot;
    EXPECT_NEAR(j.dp_py_dt0, cd_t, 1e-4 * std::abs(cd_t) + 1e-9) << "dt0 at " << p.m_dot;
  }
}

TEST(EjectorJacobians, ElementResidualsMatchCentralDifference) {
  // Whole-element (f, J): every entry of the 4x9 analytic Jacobian vs a
  // central difference of the residual value, at operating points inside each
  // regime. This directly guards the DualN<9> assembly in
  // ejector_element_residuals_and_jacobian (the C++ home of the element
  // assembly; 1:1 with EjectorElement.residuals()). Points are chosen
  // ASYMMETRIC (p_g != p_e or t_g != t_e) and away from the s_choke/s_sub
  // seams so no partial is genuinely ~0 -- the same reason the Python FD test
  // avoids the symmetric degenerate root.
  const double A_t = 3.14e-5, A_e = 1.0e-4, A_mix = 8.0e-4, A_s = A_mix - A_e;
  const EjectorGeometry geom{A_e / A_t, A_mix / A_t};
  const double R = 287.0, ep = 0.95, es = 0.85, rec = 1.0, eps = 1e-3;
  const double sc_lo = 0.90, sc_hi = 0.999, ss_lo = 0.98, ss_hi = 1.02;
  // gamma is a frozen coefficient here (matches the element's frozen-property
  // Jacobian); use a fixed representative air value so the FD holds it fixed.
  const double g = 1.4;

  struct Case {
    const char* id;
    double u[9]; // mp, ms, mdot_out, p_g, t_g, p_e, t_e, p_out, p_py
  };
  const Case cases[] = {
      // Unchoked jet-pump regime (s_choke -> 0), asymmetric streams.
      {"jetpump", {0.0051, 0.030, 0.0351, 105000.0, 305.0, 100000.0, 295.0, 100500.0, 99000.0}},
      // Deeper unchoked, different values.
      {"jetpump2", {0.0040, 0.022, 0.0260, 130000.0, 310.0, 101325.0, 292.0, 101000.0, 98500.0}},
      // Double-choked critical regime (s_choke -> 1, outlet < P_c* -> s_sub 0).
      {"critical", {0.44, 0.13, 0.57, 700000.0, 440.0, 22850.0, 280.0, 40000.0, 100000.0}},
  };

  auto resid = [&](const double u[9], int row) {
    return ejector_element_residuals_and_jacobian(u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7],
                                                  u[8], geom, A_t, A_e, A_s, g, R, ep, es, rec, eps,
                                                  sc_lo, sc_hi, ss_lo, ss_hi)
        .residuals[row];
  };

  for (const auto& c : cases) {
    auto rj = ejector_element_residuals_and_jacobian(c.u[0], c.u[1], c.u[2], c.u[3], c.u[4], c.u[5],
                                                     c.u[6], c.u[7], c.u[8], geom, A_t, A_e, A_s, g,
                                                     R, ep, es, rec, eps, sc_lo, sc_hi, ss_lo, ss_hi);
    for (int row = 0; row < 4; ++row) {
      for (int k = 0; k < 9; ++k) {
        double a[9], b[9];
        for (int m = 0; m < 9; ++m) a[m] = b[m] = c.u[m];
        double h = 1e-6 * std::abs(c.u[k]);
        if (h == 0.0) h = 1e-6;
        a[k] += h;
        b[k] -= h;
        double cd = (resid(a, row) - resid(b, row)) / (2 * h);
        double an = rj.jacobian[row][k];
        // Residuals span mass rows (O(0.01-1)) and pressure rows (O(1e4-1e5)),
        // so the abs floor is scaled to the residual magnitude to tolerate CD
        // round-off on the pressure rows without masking mass-row errors.
        double floor = 1e-6 * std::abs(rj.residuals[row]) + 1e-7;
        EXPECT_NEAR(an, cd, 1e-5 * std::abs(cd) + floor)
            << "case " << c.id << " row " << row << " seed " << k;
      }
    }
  }
}
