// Value-parity check: the C++ port (include/ejector.h) must reproduce the
// validated Python reference model's outputs for all 31 digitized Huang
// (1999) Table 3/4 rows. This is a re-implementation check, not a physics
// validation -- accuracy against Huang's published theory column was
// already established in the Python reference model (see
// validation/ejector/README.md).

#include "ejector.h"
#include "validation/ejector/data/huang1999_reference_data.h"
#include <cmath>
#include <gtest/gtest.h>

using namespace combaero::solver;
using combaero::validation::ejector::kHuangCases;

namespace {

void expect_close(double actual, double expected, double rel_tol, const char* what,
                  const char* case_id) {
  double tol = rel_tol * std::abs(expected) + 1e-12;
  EXPECT_NEAR(actual, expected, tol) << what << " mismatch for case " << case_id;
}

} // namespace

TEST(EjectorPhysics, EntrainmentRatioMatchesReferenceForAllCases) {
  for (const auto& c : kHuangCases) {
    EjectorGeometry geom{c.area_ratio_nozzle, c.area_ratio_mix};
    EjectorEntrainmentResult res =
        ejector_entrainment_ratio(c.p_g_pa, c.t_g_k, c.p_e_pa, c.t_e_k, geom, c.gamma);

    expect_close(res.omega, c.omega, 1e-9, "omega", c.id);
    expect_close(res.mach_nozzle_exit, c.mach_nozzle_exit, 1e-9, "mach_nozzle_exit", c.id);
    expect_close(res.mach_hypothetical_throat, c.mach_hypothetical_throat, 1e-9,
                "mach_hypothetical_throat", c.id);
    expect_close(res.area_ratio_primary_core, c.area_ratio_primary_core, 1e-9,
                "area_ratio_primary_core", c.id);
    expect_close(res.area_ratio_entrained, c.area_ratio_entrained, 1e-9, "area_ratio_entrained",
                c.id);
  }
}

TEST(EjectorPhysics, CriticalBackPressureMatchesReferenceForAllCases) {
  for (const auto& c : kHuangCases) {
    EjectorGeometry geom{c.area_ratio_nozzle, c.area_ratio_mix};
    EjectorCriticalPressureResult res = ejector_critical_back_pressure(
        c.p_g_pa, c.t_g_k, c.p_e_pa, c.t_e_k, geom, c.gamma, c.r_gas, c.recovery_efficiency);

    expect_close(res.p_c_pa, c.p_c_pa, 1e-9, "p_c_pa", c.id);
    expect_close(res.p_mixed_stagnation_pa, c.p_mixed_stagnation_pa, 1e-9,
                "p_mixed_stagnation_pa", c.id);
    expect_close(res.temp_mixed_stagnation_k, c.temp_mixed_stagnation_k, 1e-9,
                "temp_mixed_stagnation_k", c.id);
    expect_close(res.lambda_mixed, c.lambda_mixed, 1e-9, "lambda_mixed", c.id);
    expect_close(res.mach_mixed, c.mach_mixed, 1e-9, "mach_mixed", c.id);
  }
}

TEST(EjectorPhysics, ChokedMassFlowScalesLinearlyWithPressureAndArea) {
  // Sanity check independent of the reference table (choked_mass_flow has
  // no dedicated reference column): mdot = C * p0 * area / sqrt(t0).
  double gamma = 1.3;
  double r_gas = 287.0;
  double eta = 0.95;
  double mdot1 = ejector_choked_mass_flow(1.0e5, 300.0, 1.0e-4, gamma, r_gas, eta);
  double mdot2 = ejector_choked_mass_flow(2.0e5, 300.0, 1.0e-4, gamma, r_gas, eta);
  double mdot3 = ejector_choked_mass_flow(1.0e5, 300.0, 2.0e-4, gamma, r_gas, eta);

  EXPECT_NEAR(mdot2, 2.0 * mdot1, 1e-9 * mdot2);
  EXPECT_NEAR(mdot3, 2.0 * mdot1, 1e-9 * mdot3);
  EXPECT_GT(mdot1, 0.0);
}

TEST(EjectorPhysics, AreaOverAstarIsUnityAtMachOne) {
  EXPECT_NEAR(ejector_area_over_astar(1.0, 1.4), 1.0, 1e-12);
}

TEST(EjectorPhysics, MachFromAreaRatioSupersonicInvertsAreaOverAstar) {
  for (double area_ratio : {1.5, 2.0, 3.271, 8.0}) {
    double gamma = 1.3;
    double mach = ejector_mach_from_area_ratio_supersonic(area_ratio, gamma);
    EXPECT_GT(mach, 1.0);
    EXPECT_NEAR(ejector_area_over_astar(mach, gamma), area_ratio, 1e-8);
  }
}

TEST(EjectorPhysics, CDNozzleMassFlowMatchesPythonReference) {
  // Value parity vs cd_nozzle_mass_flow in _ejector_huang1999.py (the R0
  // primary residual). Cases: reported unchoked root, deeply choked, healthy
  // motive (choked), and a near-threshold point.
  struct Case {
    double p0, t0, p_py, at, ae, gamma, r_gas, eta, eps, expected;
  };
  const Case cases[] = {
      {101325.0, 300.0, 100207.0, 3.14e-05, 0.0001, 1.4, 287.0, 0.95, 0.001, 0.004970173368263851},
      {101325.0, 300.0, 40000.0, 3.14e-05, 0.0001, 1.4, 287.0, 0.95, 0.001, 0.007236469139827049},
      {604000.0, 368.0, 52828.0, 3.14e-05, 0.0001, 1.4, 287.0, 0.95, 0.001, 0.03894786052281535},
      {101325.0, 300.0, 101000.0, 3.14e-05, 0.0001, 1.4, 287.0, 0.95, 0.001, 0.0026910837365242053},
  };
  for (const auto& c : cases) {
    double mdot = ejector_cd_nozzle_mass_flow(c.p0, c.t0, c.p_py, c.at, c.ae, c.gamma,
                                              c.r_gas, c.eta, c.eps);
    EXPECT_NEAR(mdot, c.expected, 1e-9 * c.expected)
        << "cd_nozzle mismatch at p_py=" << c.p_py;
  }
}
