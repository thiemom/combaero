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
