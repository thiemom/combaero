#include "solver_interface.h"
#include "thermo.h"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace combaero;
using namespace combaero::solver;

// Global flag to enable detailed difference reporting
bool g_print_jacobian_differences = false;

// File output settings (empty string means no file output configured)
std::string g_jacobian_output_file =
    ""; // Will be set by JACOBIAN_OUTPUT_FILE env var
bool g_write_jacobian_to_file = false;

// Helper to collect and report differences
struct JacobianDifference {
  std::string test_name;
  std::string derivative;
  double analytical;
  double finite_diff;
  double absolute_diff;
  double relative_diff;
};

std::vector<JacobianDifference> g_jacobian_differences;

void report_jacobian_difference(const std::string &test_name,
                                const std::string &derivative,
                                double analytical, double finite_diff) {
  double abs_diff = std::abs(analytical - finite_diff);
  double rel_diff =
      (std::abs(finite_diff) > 1e-15) ? abs_diff / std::abs(finite_diff) : 0.0;

  // Console output
  if (g_print_jacobian_differences) {
    std::cout << std::fixed << std::setprecision(6)
              << "JACOBIAN_DIFF: " << test_name << " | " << derivative
              << " | analytical=" << analytical << " | fd=" << finite_diff
              << " | abs_diff=" << abs_diff << " | rel_diff=" << rel_diff
              << std::endl;
  }

  // File output
  if (g_write_jacobian_to_file && !g_jacobian_output_file.empty()) {
    std::ofstream file(g_jacobian_output_file, std::ios::app);
    if (file.is_open()) {
      file << std::fixed << std::setprecision(12) << test_name << ","
           << derivative << "," << analytical << "," << finite_diff << ","
           << abs_diff << "," << rel_diff << std::endl;
      file.close();
    }
  }

  g_jacobian_differences.push_back(
      {test_name, derivative, analytical, finite_diff, abs_diff, rel_diff});
}

// Print summary of all collected differences
void print_jacobian_difference_summary() {
  if (g_jacobian_differences.empty())
    return;

  std::cout << "\n=== JACOBIAN DIFFERENCE SUMMARY ===" << std::endl;
  std::cout << std::setw(25) << "Test" << std::setw(20) << "Derivative"
            << std::setw(15) << "Abs Diff" << std::setw(15) << "Rel Diff"
            << std::endl;
  std::cout << std::string(75, '-') << std::endl;

  for (const auto &diff : g_jacobian_differences) {
    std::cout << std::fixed << std::setprecision(3) << std::setw(25)
              << diff.test_name << std::setw(20) << diff.derivative
              << std::setw(15) << diff.absolute_diff << std::setw(15)
              << diff.relative_diff << std::endl;
  }

  // Statistics
  double max_abs = 0, max_rel = 0;
  for (const auto &diff : g_jacobian_differences) {
    max_abs = std::max(max_abs, diff.absolute_diff);
    max_rel = std::max(max_rel, diff.relative_diff);
  }

  std::cout << std::string(75, '-') << std::endl;
  std::cout << "Maximum absolute difference: " << std::scientific << max_abs
            << std::endl;
  std::cout << "Maximum relative difference: " << std::scientific << max_rel
            << std::endl;
  std::cout << "Total comparisons: " << g_jacobian_differences.size()
            << std::endl;
  std::cout << "=== END SUMMARY ===" << std::endl;
}

// Check environment variables once at startup
bool check_jacobian_diff_flag() {
  const char *env = std::getenv("PRINT_JACOBIAN_DIFFERENCES");
  return env && (std::string(env) == "1" || std::string(env) == "true");
}

void initialize_jacobian_file_output() {
  const char *file_env = std::getenv("JACOBIAN_OUTPUT_FILE");
  if (file_env) {
    g_jacobian_output_file = std::string(file_env);
    g_write_jacobian_to_file = true;

    // Write CSV header
    std::ofstream file(g_jacobian_output_file, std::ios::app);
    if (file.is_open()) {
      // Check if file is empty, write header if needed
      file.seekp(0, std::ios::end);
      if (file.tellp() == 0) {
        file << "test_name,derivative,analytical,finite_diff,abs_diff,rel_diff"
             << std::endl;
      }
      file.close();
    }
  }
}

// Initialize flags
static bool jacobian_diff_flag_initialized = []() {
  g_print_jacobian_differences = check_jacobian_diff_flag();
  initialize_jacobian_file_output();
  return true;
}();

TEST(SolverJacobianTest, OrificeDerivatives) {
  // Boilerplate for species
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  // Assume air-like for testing
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      Y[i] = 0.77;
    else if (species_name(i) == "O2")
      Y[i] = 0.23;
  }

  double m_dot = 0.1;
  double P_total_up = 200000.0;
  double P_static_up = 190000.0;
  double T_up = 300.0;
  double P_static_down = 180000.0;
  double Cd = 0.62;
  double area = 0.001;

  double beta = 0.5;
  OrificeResult res = orifice_residuals_and_jacobian(
      m_dot, P_total_up, P_static_up, T_up, Y, P_static_down, Cd, area, beta);

  // FD checks
  double eps = 1.0;

  // d_mdot_dP_total_up
  auto res_p =
      orifice_residuals_and_jacobian(m_dot, P_total_up + eps, P_static_up, T_up,
                                     Y, P_static_down, Cd, area, beta);
  double fd_dP_total_up = (res_p.m_dot_calc - res.m_dot_calc) / eps;
  report_jacobian_difference("OrificeDerivatives", "d_mdot_dP_total_up",
                             res.d_mdot_dP_total_up, fd_dP_total_up);
  EXPECT_NEAR(res.d_mdot_dP_total_up, fd_dP_total_up, 1e-10);

  // d_mdot_dP_static_down
  auto res_down_p =
      orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y,
                                     P_static_down + eps, Cd, area, beta);
  double fd_dP_static_down = (res_down_p.m_dot_calc - res.m_dot_calc) / eps;
  report_jacobian_difference("OrificeDerivatives", "d_mdot_dP_static_down",
                             res.d_mdot_dP_static_down, fd_dP_static_down);
  EXPECT_NEAR(res.d_mdot_dP_static_down, fd_dP_static_down, 1e-10);

  // d_mdot_dP_static_up (via density)
  auto res_up_p =
      orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up + eps, T_up,
                                     Y, P_static_down, Cd, area, beta);
  double fd_dP_static_up = (res_up_p.m_dot_calc - res.m_dot_calc) / eps;
  report_jacobian_difference("OrificeDerivatives", "d_mdot_dP_static_up",
                             res.d_mdot_dP_static_up, fd_dP_static_up);
  EXPECT_NEAR(res.d_mdot_dP_static_up, fd_dP_static_up, 1e-12);

  // d_mdot_dT_up
  double eps_T = 1e-4;
  auto res_T_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up,
                                                T_up + eps_T, Y, P_static_down,
                                                Cd, area, beta);
  auto res_T_m = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up,
                                                T_up - eps_T, Y, P_static_down,
                                                Cd, area, beta);
  double fd_dT_up = (res_T_p.m_dot_calc - res_T_m.m_dot_calc) / (2.0 * eps_T);
  report_jacobian_difference("OrificeDerivatives", "d_mdot_dT_up",
                             res.d_mdot_dT_up, fd_dT_up);
  EXPECT_NEAR(res.d_mdot_dT_up, fd_dT_up, std::abs(fd_dT_up) * 1e-6);

  // Y derivatives
  double eps_Y = 1e-6;
  for (size_t i = 0; i < ns; ++i) {
    std::vector<double> Y_p = Y;
    Y_p[i] += eps_Y;
    auto res_Y_p =
        orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up,
                                       Y_p, P_static_down, Cd, area, beta);
    double fd_dYi = (res_Y_p.m_dot_calc - res.m_dot_calc) / eps_Y;
    EXPECT_NEAR(res.d_mdot_dY_up[i], fd_dYi, std::abs(fd_dYi) * 1e-3 + 1e-8);
  }
}

TEST(SolverJacobianTest, PipeDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      Y[i] = 0.77;
    else if (species_name(i) == "O2")
      Y[i] = 0.23;
  }

  double m_dot = 0.5;
  double P_total_up = 200000.0;
  double P_static_up = 190000.0;
  double T_up = 300.0;
  double P_static_down = 180000.0;
  double L = 1.0;
  double D = 0.05;
  double roughness = 1e-5;
  std::string model = "haaland";

  PipeResult res =
      channel_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y,
                                     P_static_down, L, D, roughness, model);

  double eps_m = 1e-4;
  auto res_m_p = channel_residuals_and_jacobian(
      m_dot + eps_m, P_total_up, P_static_up, T_up, Y, P_static_down, L, D,
      roughness, model);
  auto res_m_m = channel_residuals_and_jacobian(
      m_dot - eps_m, P_total_up, P_static_up, T_up, Y, P_static_down, L, D,
      roughness, model);
  double fd_dmdot = (res_m_p.dP_calc - res_m_m.dP_calc) / (2.0 * eps_m);
  report_jacobian_difference("PipeDerivatives", "d_dP_d_mdot", res.d_dP_d_mdot,
                             fd_dmdot);
  EXPECT_NEAR(res.d_dP_d_mdot, fd_dmdot, std::abs(fd_dmdot) * 1e-8);

  double eps_P = 1.0;
  auto res_P_p = channel_residuals_and_jacobian(
      m_dot, P_total_up, P_static_up + eps_P, T_up, Y, P_static_down, L, D,
      roughness, model);
  auto res_P_m = channel_residuals_and_jacobian(
      m_dot, P_total_up, P_static_up - eps_P, T_up, Y, P_static_down, L, D,
      roughness, model);
  double fd_dP_static = (res_P_p.dP_calc - res_P_m.dP_calc) / (2.0 * eps_P);
  report_jacobian_difference("PipeDerivatives", "d_dP_dP_static_up",
                             res.d_dP_dP_static_up, fd_dP_static);
  EXPECT_NEAR(res.d_dP_dP_static_up, fd_dP_static,
              std::abs(fd_dP_static) * 1e-8);

  double eps_T = 0.1;
  auto res_T_p = channel_residuals_and_jacobian(m_dot, P_total_up, P_static_up,
                                                T_up + eps_T, Y, P_static_down,
                                                L, D, roughness, model);
  auto res_T_m = channel_residuals_and_jacobian(m_dot, P_total_up, P_static_up,
                                                T_up - eps_T, Y, P_static_down,
                                                L, D, roughness, model);
  double fd_dT = (res_T_p.dP_calc - res_T_m.dP_calc) / (2.0 * eps_T);
  report_jacobian_difference("PipeDerivatives", "d_dP_dT_up", res.d_dP_dT_up,
                             fd_dT);
  EXPECT_NEAR(res.d_dP_dT_up, fd_dT, std::abs(fd_dT) * 1e-8);

  double eps_Y = 1e-6;
  for (size_t i = 0; i < ns; ++i) {
    std::vector<double> Y_p = Y;
    Y_p[i] += eps_Y;
    auto res_Y_p = channel_residuals_and_jacobian(
        m_dot, P_total_up, P_static_up, T_up, Y_p, P_static_down, L, D,
        roughness, model);
    double fd_dYi = (res_Y_p.dP_calc - res.dP_calc) / eps_Y;
    EXPECT_NEAR(res.d_dP_dY_up[i], fd_dYi, std::abs(fd_dYi) * 1e-2 + 1e-8);
  }
}

TEST(SolverJacobianTest, OrificeCompressibleDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      Y[i] = 0.77;
    else if (species_name(i) == "O2")
      Y[i] = 0.23;
  }

  double m_dot = 0.30;
  double P_total_up = 200000.0;
  double T_up = 300.0;
  double P_static_down = 150000.0;
  double Cd = 0.65;
  double area = 1e-4;
  double beta = 0.5;

  OrificeResult res = orifice_compressible_residuals_and_jacobian(
      m_dot, P_total_up, T_up, Y, P_static_down, Cd, area, beta);

  double eps_P = 1.0;
  auto res_P0_p = orifice_compressible_residuals_and_jacobian(
      m_dot, P_total_up + eps_P, T_up, Y, P_static_down, Cd, area, beta);
  auto res_P0_m = orifice_compressible_residuals_and_jacobian(
      m_dot, P_total_up - eps_P, T_up, Y, P_static_down, Cd, area, beta);
  double fd_dP0 = (res_P0_p.m_dot_calc - res_P0_m.m_dot_calc) / (2.0 * eps_P);
  report_jacobian_difference("OrificeCompressibleDerivatives",
                             "d_mdot_dP_total_up", res.d_mdot_dP_total_up,
                             fd_dP0);
  EXPECT_NEAR(res.d_mdot_dP_total_up, fd_dP0, std::abs(fd_dP0) * 1e-8 + 1e-12);

  auto res_Pb_p = orifice_compressible_residuals_and_jacobian(
      m_dot, P_total_up, T_up, Y, P_static_down + eps_P, Cd, area, beta);
  auto res_Pb_m = orifice_compressible_residuals_and_jacobian(
      m_dot, P_total_up, T_up, Y, P_static_down - eps_P, Cd, area, beta);
  double fd_dPb = (res_Pb_p.m_dot_calc - res_Pb_m.m_dot_calc) / (2.0 * eps_P);
  report_jacobian_difference("OrificeCompressibleDerivatives",
                             "d_mdot_dP_static_down", res.d_mdot_dP_static_down,
                             fd_dPb);
  EXPECT_NEAR(res.d_mdot_dP_static_down, fd_dPb,
              std::abs(fd_dPb) * 1e-7 + 1e-12);

  double eps_T = 1e-3;
  auto res_T_p = orifice_compressible_residuals_and_jacobian(
      m_dot, P_total_up, T_up + eps_T, Y, P_static_down, Cd, area, beta);
  auto res_T_m = orifice_compressible_residuals_and_jacobian(
      m_dot, P_total_up, T_up - eps_T, Y, P_static_down, Cd, area, beta);
  double fd_dT = (res_T_p.m_dot_calc - res_T_m.m_dot_calc) / (2.0 * eps_T);
  report_jacobian_difference("OrificeCompressibleDerivatives", "d_mdot_dT_up",
                             res.d_mdot_dT_up, fd_dT);
  EXPECT_NEAR(res.d_mdot_dT_up, fd_dT, std::abs(fd_dT) * 1e-8 + 1e-12);
}

TEST(SolverJacobianTest, PipeCompressibleDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      Y[i] = 0.77;
    else if (species_name(i) == "O2")
      Y[i] = 0.23;
  }

  double m_dot = 0.4;
  double P_total_up = 2.0e5;
  double T_up = 400.0;
  double P_static_down = 1.8e5;
  double L = 2.0;
  double D = 0.05;
  double roughness = 1e-4;
  std::string model = "haaland";

  PipeResult res = channel_compressible_residuals_and_jacobian(
      m_dot, P_total_up, T_up, Y, P_static_down, L, D, roughness, model);

  double eps_m = 1e-5;
  auto res_m_p = channel_compressible_residuals_and_jacobian(
      m_dot + eps_m, P_total_up, T_up, Y, P_static_down, L, D, roughness,
      model);
  auto res_m_m = channel_compressible_residuals_and_jacobian(
      m_dot - eps_m, P_total_up, T_up, Y, P_static_down, L, D, roughness,
      model);
  double fd_dmdot = (res_m_p.dP_calc - res_m_m.dP_calc) / (2.0 * eps_m);
  report_jacobian_difference("PipeCompressibleDerivatives", "d_dP_d_mdot",
                             res.d_dP_d_mdot, fd_dmdot);
  EXPECT_NEAR(res.d_dP_d_mdot, fd_dmdot, std::abs(fd_dmdot) * 1e-8 + 1e-12);

  double eps_P = 1.0;
  auto res_P_p = channel_compressible_residuals_and_jacobian(
      m_dot, P_total_up + eps_P, T_up, Y, P_static_down, L, D, roughness,
      model);
  auto res_P_m = channel_compressible_residuals_and_jacobian(
      m_dot, P_total_up - eps_P, T_up, Y, P_static_down, L, D, roughness,
      model);
  double fd_dP = (res_P_p.dP_calc - res_P_m.dP_calc) / (2.0 * eps_P);
  report_jacobian_difference("PipeCompressibleDerivatives", "d_dP_dP_static_up",
                             res.d_dP_dP_static_up, fd_dP);
  EXPECT_NEAR(res.d_dP_dP_static_up, fd_dP, std::abs(fd_dP) * 1e-8 + 1e-12);

  double eps_T = 1e-3;
  auto res_T_p = channel_compressible_residuals_and_jacobian(
      m_dot, P_total_up, T_up + eps_T, Y, P_static_down, L, D, roughness,
      model);
  auto res_T_m = channel_compressible_residuals_and_jacobian(
      m_dot, P_total_up, T_up - eps_T, Y, P_static_down, L, D, roughness,
      model);
  double fd_dT = (res_T_p.dP_calc - res_T_m.dP_calc) / (2.0 * eps_T);
  report_jacobian_difference("PipeCompressibleDerivatives", "d_dP_dT_up",
                             res.d_dP_dT_up, fd_dT);
  EXPECT_NEAR(res.d_dP_dT_up, fd_dT, std::abs(fd_dT) * 1e-8 + 1e-12);
}

TEST(SolverJacobianTest, MomentumChamberDerivatives) {
  // Test momentum chamber residual: P_total = P + 0.5 * rho * v^2
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    Y[i] = (i == species_index_from_name("N2"))   ? 0.767
           : (i == species_index_from_name("O2")) ? 0.233
                                                  : 0.0;
  }

  double P = 101325.0;
  double P_total = 105000.0;
  double m_dot = 0.5;
  double T = 600.0;
  double area = 0.01; // 0.01 m^2

  MomentumChamberResult res =
      momentum_chamber_residual_and_jacobian(P, P_total, m_dot, T, Y, area);

  // Finite difference checks
  double eps_P = 1.0;
  auto res_P_p = momentum_chamber_residual_and_jacobian(P + eps_P, P_total,
                                                        m_dot, T, Y, area);
  double fd_dP = (res_P_p.residual - res.residual) / eps_P;
  report_jacobian_difference("MomentumChamberDerivatives", "d_res_dP",
                             res.d_res_dP, fd_dP);
  EXPECT_NEAR(res.d_res_dP, fd_dP, std::abs(fd_dP) * 1e-6 + 1e-10);

  auto res_Pt_p = momentum_chamber_residual_and_jacobian(P, P_total + eps_P,
                                                         m_dot, T, Y, area);
  double fd_dP_total = (res_Pt_p.residual - res.residual) / eps_P;
  report_jacobian_difference("MomentumChamberDerivatives", "d_res_dP_total",
                             res.d_res_dP_total, fd_dP_total);
  EXPECT_NEAR(res.d_res_dP_total, fd_dP_total,
              std::abs(fd_dP_total) * 1e-8 + 1e-12);

  double eps_mdot = 1e-4;
  auto res_mdot_p = momentum_chamber_residual_and_jacobian(
      P, P_total, m_dot + eps_mdot, T, Y, area);
  auto res_mdot_m = momentum_chamber_residual_and_jacobian(
      P, P_total, m_dot - eps_mdot, T, Y, area);
  double fd_dmdot =
      (res_mdot_p.residual - res_mdot_m.residual) / (2.0 * eps_mdot);
  report_jacobian_difference("MomentumChamberDerivatives", "d_res_dmdot",
                             res.d_res_dmdot, fd_dmdot);
  EXPECT_NEAR(res.d_res_dmdot, fd_dmdot, std::abs(fd_dmdot) * 1e-8 + 1e-12);

  // Verify residual is correct: P_total - P - 0.5*rho*v^2 = 0
  std::vector<double> X = mass_to_mole(normalize_fractions(Y));
  double rho = density(T, P, X);
  double v = m_dot / (rho * area);
  double q_dynamic = 0.5 * rho * v * v;
  double expected_residual = P_total - P - q_dynamic;
  EXPECT_NEAR(res.residual, expected_residual, 1e-6);
}

TEST(SolverJacobianTest, AdiabaticTCompleteDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);

  // Fuel-rich mixture (CH4 + air, phi > 1)
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "CH4")
      X[i] = 0.12;
    else if (species_name(i) == "O2")
      X[i] = 0.18;
    else if (species_name(i) == "N2")
      X[i] = 0.70;
  }

  double T_in = 300.0;
  double P = 101325.0;

  auto [T_ad, dT_ad_dT_in, X_products] =
      adiabatic_T_complete_and_jacobian_T(T_in, P, X);

  // Central FD for dT_ad/dT_in
  double eps_T = 1e-3;
  auto [T_ad_plus, dummy1, dummy2] =
      adiabatic_T_complete_and_jacobian_T(T_in + eps_T, P, X);
  auto [T_ad_minus, dummy3, dummy4] =
      adiabatic_T_complete_and_jacobian_T(T_in - eps_T, P, X);
  double fd_dT_ad_dT_in = (T_ad_plus - T_ad_minus) / (2.0 * eps_T);

  report_jacobian_difference("AdiabaticTCompleteDerivatives", "dT_ad_dT_in",
                             dT_ad_dT_in, fd_dT_ad_dT_in);
  EXPECT_NEAR(dT_ad_dT_in, fd_dT_ad_dT_in,
              std::abs(fd_dT_ad_dT_in) * 1e-8 + 1e-12);
  EXPECT_GT(T_ad, T_in); // Sanity: adiabatic T should be elevated
}

TEST(SolverJacobianTest, AdiabaticTEquilibriumDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);

  // Stoichiometric mixture (CH4 + air)
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "CH4")
      X[i] = 0.095;
    else if (species_name(i) == "O2")
      X[i] = 0.19;
    else if (species_name(i) == "N2")
      X[i] = 0.715;
  }

  double T_in = 300.0;
  double P = 101325.0;

  auto [T_ad, dT_ad_dT_in, dT_ad_dP, X_products] =
      adiabatic_T_equilibrium_and_jacobians(T_in, P, X);

  // Central FD for dT_ad/dT_in
  double eps_T = 1e-3;
  auto [T_ad_T_plus, dummy1, dummy2, dummy3] =
      adiabatic_T_equilibrium_and_jacobians(T_in + eps_T, P, X);
  auto [T_ad_T_minus, dummy4, dummy5, dummy6] =
      adiabatic_T_equilibrium_and_jacobians(T_in - eps_T, P, X);
  double fd_dT_ad_dT_in = (T_ad_T_plus - T_ad_T_minus) / (2.0 * eps_T);

  report_jacobian_difference("AdiabaticTEquilibriumDerivatives", "dT_ad_dT_in",
                             dT_ad_dT_in, fd_dT_ad_dT_in);
  EXPECT_NEAR(dT_ad_dT_in, fd_dT_ad_dT_in,
              std::abs(fd_dT_ad_dT_in) * 1e-8 + 1e-12);

  // Central FD for dT_ad/dP
  double eps_P = 10.0;
  auto [T_ad_P_plus, dummy7, dummy8, dummy9] =
      adiabatic_T_equilibrium_and_jacobians(T_in, P + eps_P, X);
  auto [T_ad_P_minus, dummy10, dummy11, dummy12] =
      adiabatic_T_equilibrium_and_jacobians(T_in, P - eps_P, X);
  double fd_dT_ad_dP = (T_ad_P_plus - T_ad_P_minus) / (2.0 * eps_P);

  report_jacobian_difference("AdiabaticTEquilibriumDerivatives", "dT_ad_dP",
                             dT_ad_dP, fd_dT_ad_dP);
  EXPECT_NEAR(dT_ad_dP, fd_dT_ad_dP, std::abs(fd_dT_ad_dP) * 1e-8 + 1e-12);
  EXPECT_GT(T_ad, T_in); // Sanity: adiabatic T should be elevated
}

TEST(SolverJacobianTest, T0FromStaticDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double T = 300.0;
  double M = 0.5;

  auto [T0, dT0_dM] = T0_from_static_and_jacobian_M(T, M, X);

  // Central FD for dT0/dM
  double eps_M = 1e-5;
  auto [T0_plus, dummy1] = T0_from_static_and_jacobian_M(T, M + eps_M, X);
  auto [T0_minus, dummy2] = T0_from_static_and_jacobian_M(T, M - eps_M, X);
  double fd_dT0_dM = (T0_plus - T0_minus) / (2.0 * eps_M);

  report_jacobian_difference("T0FromStaticDerivatives", "dT0_dM", dT0_dM,
                             fd_dT0_dM);
  EXPECT_NEAR(dT0_dM, fd_dT0_dM, std::abs(fd_dT0_dM) * 1e-6 + 1e-10);
  EXPECT_GT(T0, T); // Sanity: T0 > T for M > 0

  // Test at higher Mach
  M = 0.8;
  auto [T0_high, dT0_dM_high] = T0_from_static_and_jacobian_M(T, M, X);
  auto [T0_high_plus, dummy3] = T0_from_static_and_jacobian_M(T, M + eps_M, X);
  auto [T0_high_minus, dummy4] = T0_from_static_and_jacobian_M(T, M - eps_M, X);
  double fd_dT0_dM_high = (T0_high_plus - T0_high_minus) / (2.0 * eps_M);

  report_jacobian_difference("T0FromStaticDerivatives", "dT0_dM_high",
                             dT0_dM_high, fd_dT0_dM_high);

  EXPECT_NEAR(dT0_dM_high, fd_dT0_dM_high,
              std::abs(fd_dT0_dM_high) * 1e-6 + 1e-10);
}

TEST(SolverJacobianTest, P0FromStaticDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double P = 101325.0;
  double T = 300.0;
  double M = 0.5;

  auto [P0, dP0_dM] = P0_from_static_and_jacobian_M(P, T, M, X);

  // Central FD for dP0/dM
  double eps_M = 1e-5;
  auto [P0_plus, dummy1] = P0_from_static_and_jacobian_M(P, T, M + eps_M, X);
  auto [P0_minus, dummy2] = P0_from_static_and_jacobian_M(P, T, M - eps_M, X);
  double fd_dP0_dM = (P0_plus - P0_minus) / (2.0 * eps_M);

  report_jacobian_difference("P0FromStaticDerivatives", "dP0_dM", dP0_dM,
                             fd_dP0_dM);
  EXPECT_NEAR(dP0_dM, fd_dP0_dM, std::abs(fd_dP0_dM) * 1e-6 + 1e-8);
  EXPECT_GT(P0, P); // Sanity: P0 > P for M > 0

  // Test at higher Mach
  M = 0.8;
  auto [P0_high, dP0_dM_high] = P0_from_static_and_jacobian_M(P, T, M, X);
  auto [P0_high_plus, dummy3] =
      P0_from_static_and_jacobian_M(P, T, M + eps_M, X);
  auto [P0_high_minus, dummy4] =
      P0_from_static_and_jacobian_M(P, T, M - eps_M, X);
  double fd_dP0_dM_high = (P0_high_plus - P0_high_minus) / (2.0 * eps_M);

  report_jacobian_difference("P0FromStaticDerivatives", "dP0_dM_high",
                             dP0_dM_high, fd_dP0_dM_high);
  EXPECT_NEAR(dP0_dM_high, fd_dP0_dM_high,
              std::abs(fd_dP0_dM_high) * 1e-6 + 1e-8);
}

TEST(SolverJacobianTest, OrificeCompressibleMdotDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double T0 = 350.0;
  double P0 = 200000.0;
  double P_back = 150000.0;
  double Cd = 0.65;
  double area = 1e-4;
  double beta = 0.5;

  auto [mdot, dmdot_dP0, dmdot_dP_back, dmdot_dT0] =
      orifice_compressible_mdot_and_jacobian(T0, P0, P_back, X, Cd, area, beta);

  // Central FD for dmdot/dT0
  double eps_T = 1e-3;
  auto [mdot_T_plus, dummy1, dummy2, dummy3] =
      orifice_compressible_mdot_and_jacobian(T0 + eps_T, P0, P_back, X, Cd,
                                             area, beta);
  auto [mdot_T_minus, dummy4, dummy5, dummy6] =
      orifice_compressible_mdot_and_jacobian(T0 - eps_T, P0, P_back, X, Cd,
                                             area, beta);
  double fd_dT0 = (mdot_T_plus - mdot_T_minus) / (2.0 * eps_T);
  EXPECT_NEAR(dmdot_dT0, fd_dT0, std::abs(fd_dT0) * 1e-6 + 1e-12);

  // Central FD for dmdot/dP0
  double eps_P = 1.0;
  auto [mdot_P_plus, dummy7, dummy8, dummy9] =
      orifice_compressible_mdot_and_jacobian(T0, P0 + eps_P, P_back, X, Cd,
                                             area, beta);
  auto [mdot_P_minus, dummy10, dummy11, dummy12] =
      orifice_compressible_mdot_and_jacobian(T0, P0 - eps_P, P_back, X, Cd,
                                             area, beta);
  double fd_dP0 = (mdot_P_plus - mdot_P_minus) / (2.0 * eps_P);
  EXPECT_NEAR(dmdot_dP0, fd_dP0, std::abs(fd_dP0) * 1e-6 + 1e-12);

  // Central FD for dmdot/dP_back
  auto [mdot_Pb_plus, dummy13, dummy14, dummy15] =
      orifice_compressible_mdot_and_jacobian(T0, P0, P_back + eps_P, X, Cd,
                                             area, beta);
  auto [mdot_Pb_minus, dummy16, dummy17, dummy18] =
      orifice_compressible_mdot_and_jacobian(T0, P0, P_back - eps_P, X, Cd,
                                             area, beta);
  double fd_dP_back = (mdot_Pb_plus - mdot_Pb_minus) / (2.0 * eps_P);
  EXPECT_NEAR(dmdot_dP_back, fd_dP_back, std::abs(fd_dP_back) * 1e-6 + 1e-12);
}

TEST(SolverJacobianTest, PipeCompressibleMdotDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double T_in = 400.0;
  double P_in = 200000.0;
  double u_in = 100.0;
  double L = 2.0;
  double D = 0.05;
  double roughness = 1e-4;
  std::string model = "haaland";

  auto [dP, ddP_dP_in, ddP_dT_in, ddP_du_in] =
      channel_compressible_mdot_and_jacobian(T_in, P_in, u_in, X, L, D,
                                             roughness, model);

  // Central FD for ddP/dT_in
  double eps_T = 1e-3;
  auto [dP_T_plus, dummy1, dummy2, dummy3] =
      channel_compressible_mdot_and_jacobian(T_in + eps_T, P_in, u_in, X, L, D,
                                             roughness, model);
  auto [dP_T_minus, dummy4, dummy5, dummy6] =
      channel_compressible_mdot_and_jacobian(T_in - eps_T, P_in, u_in, X, L, D,
                                             roughness, model);
  double fd_dT_in = (dP_T_plus - dP_T_minus) / (2.0 * eps_T);
  EXPECT_NEAR(ddP_dT_in, fd_dT_in, std::abs(fd_dT_in) * 1e-6 + 1e-8);

  // Central FD for ddP/dP_in
  double eps_P = 1.0;
  auto [dP_P_plus, dummy19, dummy20, dummy21] =
      channel_compressible_mdot_and_jacobian(T_in, P_in + eps_P, u_in, X, L, D,
                                             roughness, model);
  auto [dP_P_minus, dummy22, dummy23, dummy24] =
      channel_compressible_mdot_and_jacobian(T_in, P_in - eps_P, u_in, X, L, D,
                                             roughness, model);
  double fd_dP_in = (dP_P_plus - dP_P_minus) / (2.0 * eps_P);
  EXPECT_NEAR(ddP_dP_in, fd_dP_in, std::abs(fd_dP_in) * 1e-6 + 1e-10);

  // Central FD for ddP/du_in
  double eps_u = 1e-2;
  auto [dP_u_plus, dummy25, dummy26, dummy27] =
      channel_compressible_mdot_and_jacobian(T_in, P_in, u_in + eps_u, X, L, D,
                                             roughness, model);
  auto [dP_u_minus, dummy28, dummy29, dummy30] =
      channel_compressible_mdot_and_jacobian(T_in, P_in, u_in - eps_u, X, L, D,
                                             roughness, model);
  double fd_du_in = (dP_u_plus - dP_u_minus) / (2.0 * eps_u);
  EXPECT_NEAR(ddP_du_in, fd_du_in, std::abs(fd_du_in) * 1e-6 + 1e-7);
}

TEST(SolverJacobianTest, MachNumberDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double v = 200.0;
  double T = 300.0;

  auto [M, dM_dv] = mach_number_and_jacobian_v(v, T, X);

  // Central FD for dM/dv
  double eps_v = 1e-3;
  auto [M_plus, dummy1] = mach_number_and_jacobian_v(v + eps_v, T, X);
  auto [M_minus, dummy2] = mach_number_and_jacobian_v(v - eps_v, T, X);
  double fd_dM_dv = (M_plus - M_minus) / (2.0 * eps_v);
  EXPECT_NEAR(dM_dv, fd_dM_dv, std::abs(fd_dM_dv) * 1e-6 + 1e-14);
}

TEST(SolverJacobianTest, AdiabaticWallDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double T_static = 350.0;
  double v = 150.0;
  double T = 300.0;
  double P = 101325.0;

  auto [T_aw, dT_aw_dv] = T_adiabatic_wall_and_jacobian_v(T_static, v, T, P, X);

  // Central FD for dT_aw/dv
  double eps_v = 1e-3;
  auto [T_aw_plus, dummy1] =
      T_adiabatic_wall_and_jacobian_v(T_static, v + eps_v, T, P, X);
  auto [T_aw_minus, dummy2] =
      T_adiabatic_wall_and_jacobian_v(T_static, v - eps_v, T, P, X);
  double fd_dT_aw_dv = (T_aw_plus - T_aw_minus) / (2.0 * eps_v);
  report_jacobian_difference("AdiabaticWallDerivatives", "dT_aw_dv", dT_aw_dv,
                             fd_dT_aw_dv);
  EXPECT_NEAR(dT_aw_dv, fd_dT_aw_dv, std::abs(fd_dT_aw_dv) * 1e-6 + 1e-11);
}

// Test fixture to print summary at the end
class JacobianTestEnvironment : public ::testing::Environment {
public:
  void TearDown() override { print_jacobian_difference_summary(); }
};

// Register the environment
::testing::Environment *const jacobian_env =
    ::testing::AddGlobalTestEnvironment(new JacobianTestEnvironment());
