#include <gtest/gtest.h>
#include "humidair.h"
#include "solver_interface.h"
#include "state.h"
#include "thermo.h"
#include <cmath>
#include <vector>

using namespace combaero;
namespace sv = combaero::solver;

// Helper: build air-like mass fraction vector
static std::vector<double> air_Y() {
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (std::size_t k = 0; k < ns; ++k) {
    if (species_name(k) == "N2")
      Y[k] = 0.767;
    else if (species_name(k) == "O2")
      Y[k] = 0.233;
  }
  return Y;
}

// Helper: build CH4-like mass fraction vector
static std::vector<double> ch4_Y() {
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (std::size_t k = 0; k < ns; ++k) {
    if (species_name(k) == "CH4")
      Y[k] = 1.0;
  }
  return Y;
}

// =============================================================================
// 1. Mixer delta_h accuracy: known Q → known ΔT
// =============================================================================

TEST(EnergyBoundaryTest, MixerDeltaH_PositiveHeating) {
  // Single air stream at 300K, add 50 kJ/kg → T should rise by ~50
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 300.0;
  st.P_total = 101325.0;
  st.Y = air_Y();

  double delta_h = 50000.0; // 50 kJ/kg
  auto res = sv::mixer_from_streams_and_jacobians({st}, delta_h);

  // Expected: T_mix > 300 K, roughly T_base + delta_h / cp_air
  auto res_base = sv::mixer_from_streams_and_jacobians({st});
  EXPECT_NEAR(res_base.T_mix, 300.0, 0.1);

  // cp_air ~ 1005 J/(kg·K), so ΔT ~ 50000/1005 ~ 49.75 K
  double dT = res.T_mix - res_base.T_mix;
  EXPECT_GT(dT, 40.0);
  EXPECT_LT(dT, 60.0);

  // Verify via enthalpy balance: h(T_mix, X) should equal h(300, X) + delta_h
  std::vector<double> X = mass_to_mole(normalize_fractions(st.Y));
  double h_base = h_mass(300.0, X);
  double h_out = h_mass(res.T_mix, X);
  EXPECT_NEAR(h_out, h_base + delta_h, std::abs(h_base + delta_h) * 1e-6);
}

TEST(EnergyBoundaryTest, MixerDeltaH_NegativeCooling) {
  // Single air stream at 800K, remove 100 kJ/kg → T should drop
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 800.0;
  st.P_total = 101325.0;
  st.Y = air_Y();

  double delta_h = -100000.0; // -100 kJ/kg
  auto res = sv::mixer_from_streams_and_jacobians({st}, delta_h);

  std::vector<double> X = mass_to_mole(normalize_fractions(st.Y));
  double h_base = h_mass(800.0, X);
  double h_out = h_mass(res.T_mix, X);
  EXPECT_NEAR(h_out, h_base + delta_h, std::abs(h_base + delta_h) * 1e-6);
  EXPECT_LT(res.T_mix, 800.0);
}

TEST(EnergyBoundaryTest, MixerDeltaH_ZeroIsNoop) {
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 500.0;
  st.P_total = 200000.0;
  st.Y = air_Y();

  auto res_zero = sv::mixer_from_streams_and_jacobians({st}, 0.0);
  auto res_none = sv::mixer_from_streams_and_jacobians({st});

  EXPECT_DOUBLE_EQ(res_zero.T_mix, res_none.T_mix);
  EXPECT_DOUBLE_EQ(res_zero.P_total_mix, res_none.P_total_mix);
}

TEST(EnergyBoundaryTest, MixerDeltaH_TwoStreams) {
  // Two streams with different T, add delta_h
  sv::Stream s1;
  s1.m_dot = 0.5;
  s1.T = 300.0;
  s1.P_total = 101325.0;
  s1.Y = air_Y();

  sv::Stream s2;
  s2.m_dot = 0.5;
  s2.T = 600.0;
  s2.P_total = 101325.0;
  s2.Y = air_Y();

  double delta_h = 20000.0; // 20 kJ/kg
  auto res = sv::mixer_from_streams_and_jacobians({s1, s2}, delta_h);
  auto res_base = sv::mixer_from_streams_and_jacobians({s1, s2});

  // Mixed T_base should be ~450K (mass-weighted average for same composition)
  EXPECT_NEAR(res_base.T_mix, 450.0, 2.0);

  // delta_h should raise it by ~20 K
  double dT = res.T_mix - res_base.T_mix;
  EXPECT_GT(dT, 15.0);
  EXPECT_LT(dT, 25.0);
}

// =============================================================================
// 2. Mixer dT_mix_d_delta_h Jacobian: finite-difference verification
// =============================================================================

TEST(EnergyBoundaryTest, MixerJacobian_dT_d_delta_h) {
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 500.0;
  st.P_total = 200000.0;
  st.Y = air_Y();

  double delta_h = 10000.0;
  auto res = sv::mixer_from_streams_and_jacobians({st}, delta_h);

  // Finite-difference check
  double eps = 1.0; // 1 J/kg
  auto res_p = sv::mixer_from_streams_and_jacobians({st}, delta_h + eps);
  auto res_m = sv::mixer_from_streams_and_jacobians({st}, delta_h - eps);
  double fd = (res_p.T_mix - res_m.T_mix) / (2.0 * eps);

  EXPECT_NEAR(res.dT_mix_d_delta_h, fd, std::abs(fd) * 1e-4);

  // Also check: dT/d(delta_h) ≈ 1/cp_mix
  std::vector<double> X = mass_to_mole(normalize_fractions(st.Y));
  double cp_val = cp_mass(res.T_mix, X);
  EXPECT_NEAR(res.dT_mix_d_delta_h, 1.0 / cp_val, 1.0 / cp_val * 1e-3);
}

TEST(EnergyBoundaryTest, MixerJacobian_dT_d_delta_h_TwoStreams) {
  sv::Stream s1;
  s1.m_dot = 0.3;
  s1.T = 400.0;
  s1.P_total = 150000.0;
  s1.Y = air_Y();

  sv::Stream s2;
  s2.m_dot = 0.7;
  s2.T = 700.0;
  s2.P_total = 150000.0;
  s2.Y = air_Y();

  double delta_h = -30000.0;
  auto res = sv::mixer_from_streams_and_jacobians({s1, s2}, delta_h);

  double eps = 1.0;
  auto res_p = sv::mixer_from_streams_and_jacobians({s1, s2}, delta_h + eps);
  auto res_m = sv::mixer_from_streams_and_jacobians({s1, s2}, delta_h - eps);
  double fd = (res_p.T_mix - res_m.T_mix) / (2.0 * eps);

  EXPECT_NEAR(res.dT_mix_d_delta_h, fd, std::abs(fd) * 1e-4);
}

// =============================================================================
// 3. Mixer: existing stream Jacobians unaffected by delta_h
// =============================================================================

TEST(EnergyBoundaryTest, MixerStreamJacobians_WithDeltaH) {
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 500.0;
  st.P_total = 200000.0;
  st.Y = air_Y();

  double delta_h = 25000.0;
  auto res = sv::mixer_from_streams_and_jacobians({st}, delta_h);

  // FD: dT_mix/dT_stream
  double eps_T = 0.01;
  sv::Stream st_Tp = st;
  st_Tp.T += eps_T;
  sv::Stream st_Tm = st;
  st_Tm.T -= eps_T;
  auto res_Tp = sv::mixer_from_streams_and_jacobians({st_Tp}, delta_h);
  auto res_Tm = sv::mixer_from_streams_and_jacobians({st_Tm}, delta_h);
  double fd_dT = (res_Tp.T_mix - res_Tm.T_mix) / (2.0 * eps_T);
  EXPECT_NEAR(res.dT_mix_d_stream[0].d_T, fd_dT, std::abs(fd_dT) * 1e-4);

  // FD: dT_mix/dmdot
  double eps_m = 1e-4;
  sv::Stream st_mp = st;
  st_mp.m_dot += eps_m;
  sv::Stream st_mm = st;
  st_mm.m_dot -= eps_m;
  auto res_mp = sv::mixer_from_streams_and_jacobians({st_mp}, delta_h);
  auto res_mm = sv::mixer_from_streams_and_jacobians({st_mm}, delta_h);
  double fd_dmdot = (res_mp.T_mix - res_mm.T_mix) / (2.0 * eps_m);
  EXPECT_NEAR(res.dT_mix_d_stream[0].d_mdot, fd_dmdot,
              std::abs(fd_dmdot) * 1e-3 + 1e-6);
}

// =============================================================================
// 4. Complete combustion with post-combustion delta_h
// =============================================================================

TEST(EnergyBoundaryTest, CompleteCombustion_PostDeltaH) {
  // Combust CH4 + air, then apply heat loss
  sv::Stream fuel;
  fuel.m_dot = 0.05;
  fuel.T = 300.0;
  fuel.P_total = 101325.0;
  fuel.Y = ch4_Y();

  sv::Stream ox;
  ox.m_dot = 1.0;
  ox.T = 600.0;
  ox.P_total = 101325.0;
  ox.Y = air_Y();

  double P = 101325.0;

  auto res_base = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, 0.0);

  // With heat loss of -200 kJ/kg
  double delta_h = -200000.0;
  auto res_cooled = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, delta_h);

  // Cooled T should be lower
  EXPECT_LT(res_cooled.T_mix, res_base.T_mix);

  // Verify enthalpy balance: h(T_cooled, X_products) = h(T_ad, X_products) + delta_h
  std::vector<double> X_b = mass_to_mole(normalize_fractions(res_base.Y_mix));
  double h_ad = h_mass(res_base.T_mix, X_b);
  double h_cooled = h_mass(res_cooled.T_mix, X_b);
  EXPECT_NEAR(h_cooled, h_ad + delta_h, std::abs(h_ad + delta_h) * 1e-5);

  // Products composition should be unchanged (complete combustion)
  for (std::size_t k = 0; k < res_base.Y_mix.size(); ++k) {
    EXPECT_NEAR(res_cooled.Y_mix[k], res_base.Y_mix[k], 1e-10);
  }
}

TEST(EnergyBoundaryTest, CompleteCombustion_DeltaH_Zero_IsNoop) {
  sv::Stream fuel;
  fuel.m_dot = 0.05;
  fuel.T = 300.0;
  fuel.P_total = 101325.0;
  fuel.Y = ch4_Y();

  sv::Stream ox;
  ox.m_dot = 1.0;
  ox.T = 600.0;
  ox.P_total = 101325.0;
  ox.Y = air_Y();

  double P = 101325.0;

  auto res_zero = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, 0.0);
  auto res_none = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P);

  EXPECT_DOUBLE_EQ(res_zero.T_mix, res_none.T_mix);
}

// =============================================================================
// 5. Complete combustion dT_d_delta_h Jacobian: FD verification
// =============================================================================

TEST(EnergyBoundaryTest, CompleteCombustion_Jacobian_dT_d_delta_h) {
  sv::Stream fuel;
  fuel.m_dot = 0.05;
  fuel.T = 300.0;
  fuel.P_total = 101325.0;
  fuel.Y = ch4_Y();

  sv::Stream ox;
  ox.m_dot = 1.0;
  ox.T = 600.0;
  ox.P_total = 101325.0;
  ox.Y = air_Y();

  double P = 101325.0;
  double delta_h = -100000.0;

  auto res = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, delta_h);

  double eps = 1.0;
  auto res_p = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, delta_h + eps);
  auto res_m = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, delta_h - eps);
  double fd = (res_p.T_mix - res_m.T_mix) / (2.0 * eps);

  EXPECT_NEAR(res.dT_mix_d_delta_h, fd, std::abs(fd) * 1e-3);
}

// =============================================================================
// 6. Complete combustion: stream Jacobians still correct with delta_h
// =============================================================================

TEST(EnergyBoundaryTest, CompleteCombustion_StreamJacobians_WithDeltaH) {
  sv::Stream fuel;
  fuel.m_dot = 0.05;
  fuel.T = 300.0;
  fuel.P_total = 101325.0;
  fuel.Y = ch4_Y();

  sv::Stream ox;
  ox.m_dot = 1.0;
  ox.T = 600.0;
  ox.P_total = 101325.0;
  ox.Y = air_Y();

  double P = 101325.0;
  double delta_h = -150000.0;

  auto res = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, delta_h);

  // FD: dT_out/dT_air
  double eps_T = 0.5;
  sv::Stream ox_Tp = ox;
  ox_Tp.T += eps_T;
  sv::Stream ox_Tm = ox;
  ox_Tm.T -= eps_T;
  auto res_Tp = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox_Tp}, P, delta_h);
  auto res_Tm = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox_Tm}, P, delta_h);
  double fd_dT = (res_Tp.T_mix - res_Tm.T_mix) / (2.0 * eps_T);

  // Air is stream index 1
  EXPECT_NEAR(res.dT_mix_d_stream[1].d_T, fd_dT, std::abs(fd_dT) * 5e-3);
}

// =============================================================================
// 7. Stream mix() with delta_h
// =============================================================================

TEST(EnergyBoundaryTest, StreamMix_DeltaH) {
  // Create two streams with State objects
  Stream s1;
  s1.mdot = 0.5;
  s1.state.set_TPX(300.0, 101325.0, standard_dry_air_composition());

  Stream s2;
  s2.mdot = 0.5;
  s2.state.set_TPX(600.0, 101325.0, standard_dry_air_composition());

  // Base mix (adiabatic)
  Stream mixed_base = mix({s1, s2});
  double T_base = mixed_base.state.T;

  // Mix with 30 kJ/kg heat input
  Stream mixed_heated = mix({s1, s2}, -1.0, 30000.0);
  EXPECT_GT(mixed_heated.state.T, T_base);

  // Mix with -30 kJ/kg heat loss
  Stream mixed_cooled = mix({s1, s2}, -1.0, -30000.0);
  EXPECT_LT(mixed_cooled.state.T, T_base);

  // Enthalpy check
  double h_heated = mixed_heated.state.h();
  double h_base = mixed_base.state.h();
  EXPECT_NEAR(h_heated - h_base, 30000.0, 1.0);
}

// =============================================================================
// 8. Equilibrium combustion with delta_h
// =============================================================================

TEST(EnergyBoundaryTest, EquilibriumCombustion_PostDeltaH) {
  sv::Stream fuel;
  fuel.m_dot = 0.05;
  fuel.T = 300.0;
  fuel.P_total = 101325.0;
  fuel.Y = ch4_Y();

  sv::Stream ox;
  ox.m_dot = 1.0;
  ox.T = 600.0;
  ox.P_total = 101325.0;
  ox.Y = air_Y();

  double P = 101325.0;

  auto res_base = sv::adiabatic_T_equilibrium_and_jacobians_from_streams(
      {fuel, ox}, P, 0.0);

  double delta_h = -150000.0;
  auto res_cooled = sv::adiabatic_T_equilibrium_and_jacobians_from_streams(
      {fuel, ox}, P, delta_h);

  EXPECT_LT(res_cooled.T_mix, res_base.T_mix);

  // dT_d_delta_h FD check
  double eps = 1.0;
  auto res_p = sv::adiabatic_T_equilibrium_and_jacobians_from_streams(
      {fuel, ox}, P, delta_h + eps);
  auto res_m = sv::adiabatic_T_equilibrium_and_jacobians_from_streams(
      {fuel, ox}, P, delta_h - eps);
  double fd = (res_p.T_mix - res_m.T_mix) / (2.0 * eps);

  EXPECT_NEAR(res_cooled.dT_mix_d_delta_h, fd, std::abs(fd) * 5e-3);
}
