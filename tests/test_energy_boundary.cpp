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
// 1. Mixer Q parameter accuracy: known Q → known ΔT
// =============================================================================

TEST(EnergyBoundaryTest, MixerQ_PositiveHeating) {
  // Single air stream at 300K, mdot=1 kg/s, Q=50 kW → delta_h=50 kJ/kg
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 300.0;
  st.P_total = 101325.0;
  st.Y = air_Y();

  double Q = 50000.0; // 50 kW
  auto res = sv::mixer_from_streams_and_jacobians({st}, Q, 0.0);

  // Expected: T_mix > 300 K, roughly T_base + Q/(mdot*cp_air)
  auto res_base = sv::mixer_from_streams_and_jacobians({st}, 0.0, 0.0);
  EXPECT_NEAR(res_base.T_mix, 300.0, 0.1);

  // cp_air ~ 1005 J/(kg·K), so ΔT ~ 50000/(1.0*1005) ~ 49.75 K
  double dT = res.T_mix - res_base.T_mix;
  EXPECT_GT(dT, 40.0);
  EXPECT_LT(dT, 60.0);

  // Verify via enthalpy balance: h(T_mix, X) should equal h(300, X) + Q/mdot
  std::vector<double> X = mass_to_mole(normalize_fractions(st.Y));
  double h_base = h_mass(300.0, X);
  double h_out = h_mass(res.T_mix, X);
  double delta_h = Q / st.m_dot;
  EXPECT_NEAR(h_out, h_base + delta_h, std::abs(h_base + delta_h) * 1e-6);
}

TEST(EnergyBoundaryTest, MixerQ_NegativeCooling) {
  // Single air stream at 800K, mdot=1 kg/s, Q=-100 kW → T should drop
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 800.0;
  st.P_total = 101325.0;
  st.Y = air_Y();

  double Q = -100000.0; // -100 kW
  auto res = sv::mixer_from_streams_and_jacobians({st}, Q, 0.0);

  std::vector<double> X = mass_to_mole(normalize_fractions(st.Y));
  double h_base = h_mass(800.0, X);
  double h_out = h_mass(res.T_mix, X);
  double delta_h = Q / st.m_dot;
  EXPECT_NEAR(h_out, h_base + delta_h, std::abs(h_base + delta_h) * 1e-6);
  EXPECT_LT(res.T_mix, 800.0);
}

TEST(EnergyBoundaryTest, MixerQ_ZeroIsNoop) {
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 500.0;
  st.P_total = 200000.0;
  st.Y = air_Y();

  auto res_zero = sv::mixer_from_streams_and_jacobians({st}, 0.0, 0.0);
  auto res_none = sv::mixer_from_streams_and_jacobians({st});

  EXPECT_DOUBLE_EQ(res_zero.T_mix, res_none.T_mix);
  EXPECT_DOUBLE_EQ(res_zero.P_total_mix, res_none.P_total_mix);
}

TEST(EnergyBoundaryTest, MixerQ_TwoStreams) {
  // Two streams with different T, mdot_tot=1.0 kg/s, Q=20 kW
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

  double Q = 20000.0; // 20 kW
  auto res = sv::mixer_from_streams_and_jacobians({s1, s2}, Q, 0.0);
  auto res_base = sv::mixer_from_streams_and_jacobians({s1, s2}, 0.0, 0.0);

  // Mixed T_base should be ~450K (mass-weighted average for same composition)
  EXPECT_NEAR(res_base.T_mix, 450.0, 2.0);

  // Q/(mdot_tot*cp) should raise it by ~20 K
  double dT = res.T_mix - res_base.T_mix;
  EXPECT_GT(dT, 15.0);
  EXPECT_LT(dT, 25.0);
}

// =============================================================================
// 2. Mixer: Fraction parameter accuracy and enthalpy conservation
// =============================================================================

TEST(EnergyBoundaryTest, MixerFraction_NegativeCooling) {
  // Single air stream, fraction=-0.05 → 5% enthalpy loss
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 600.0;
  st.P_total = 101325.0;
  st.Y = air_Y();

  double fraction = -0.05;
  auto res = sv::mixer_from_streams_and_jacobians({st}, 0.0, fraction);
  auto res_base = sv::mixer_from_streams_and_jacobians({st}, 0.0, 0.0);

  // T should drop
  EXPECT_LT(res.T_mix, res_base.T_mix);

  // Enthalpy should be reduced by 5%
  std::vector<double> X = mass_to_mole(normalize_fractions(st.Y));
  double h_base = h_mass(res_base.T_mix, X);
  double h_out = h_mass(res.T_mix, X);
  EXPECT_NEAR(h_out / h_base, 1.0 + fraction, 1e-5);
}

TEST(EnergyBoundaryTest, MixerFraction_PositiveHeating) {
  // Two streams, fraction=0.10 → 10% enthalpy gain
  sv::Stream s1;
  s1.m_dot = 0.5;
  s1.T = 400.0;
  s1.P_total = 101325.0;
  s1.Y = air_Y();

  sv::Stream s2;
  s2.m_dot = 0.5;
  s2.T = 800.0;
  s2.P_total = 101325.0;
  s2.Y = air_Y();

  double fraction = 0.10;
  auto res = sv::mixer_from_streams_and_jacobians({s1, s2}, 0.0, fraction);
  auto res_base = sv::mixer_from_streams_and_jacobians({s1, s2}, 0.0, 0.0);

  // T should rise
  EXPECT_GT(res.T_mix, res_base.T_mix);

  // Enthalpy should increase by 10%
  std::vector<double> X = mass_to_mole(normalize_fractions(s1.Y));
  double h_base = h_mass(res_base.T_mix, X);
  double h_out = h_mass(res.T_mix, X);
  EXPECT_NEAR(h_out / h_base, 1.0 + fraction, 1e-5);
}

TEST(EnergyBoundaryTest, MixerCombined_Q_and_Fraction) {
  // Both Q and fraction set: effective Q = Q + fraction * H_tot
  sv::Stream st;
  st.m_dot = 2.0;
  st.T = 500.0;
  st.P_total = 101325.0;
  st.Y = air_Y();

  double Q = 10000.0; // 10 kW
  double fraction = -0.02; // -2%

  auto res = sv::mixer_from_streams_and_jacobians({st}, Q, fraction);

  // Compute expected delta_h manually
  std::vector<double> X = mass_to_mole(normalize_fractions(st.Y));
  double h_in = h_mass(st.T, X);
  double H_tot = st.m_dot * h_in;
  double Q_eff = Q + fraction * H_tot;
  double delta_h_expected = Q_eff / st.m_dot;

  double h_out = h_mass(res.T_mix, X);
  EXPECT_NEAR(h_out, h_in + delta_h_expected, std::abs(h_in) * 1e-5);
}

// =============================================================================
// 3. Mixer Q Jacobians: FD verification of dT/dmdot with Q correction
// =============================================================================

TEST(EnergyBoundaryTest, MixerJacobian_Q_dT_dmdot_SingleStream) {
  // Verify dT/dmdot includes Q correction: -Q/mdot²
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 500.0;
  st.P_total = 200000.0;
  st.Y = air_Y();

  double Q = 50000.0; // 50 kW
  auto res = sv::mixer_from_streams_and_jacobians({st}, Q, 0.0);

  // FD: dT_mix/dmdot
  double eps_m = 1e-5;
  sv::Stream st_mp = st;
  st_mp.m_dot += eps_m;
  sv::Stream st_mm = st;
  st_mm.m_dot -= eps_m;
  auto res_mp = sv::mixer_from_streams_and_jacobians({st_mp}, Q, 0.0);
  auto res_mm = sv::mixer_from_streams_and_jacobians({st_mm}, Q, 0.0);
  double fd_dmdot = (res_mp.T_mix - res_mm.T_mix) / (2.0 * eps_m);

  EXPECT_NEAR(res.dT_mix_d_stream[0].d_mdot, fd_dmdot,
              std::max(std::abs(fd_dmdot) * 1e-4, 1e-6));
}

TEST(EnergyBoundaryTest, MixerJacobian_Q_dT_dmdot_TwoStreams) {
  // Two streams with Q: verify both dT/dmdot_i derivatives
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

  double Q = -30000.0; // -30 kW
  auto res = sv::mixer_from_streams_and_jacobians({s1, s2}, Q, 0.0);

  // FD: dT_mix/dmdot_1
  double eps_m = 1e-5;
  sv::Stream s1_mp = s1;
  s1_mp.m_dot += eps_m;
  sv::Stream s1_mm = s1;
  s1_mm.m_dot -= eps_m;
  auto res_mp = sv::mixer_from_streams_and_jacobians({s1_mp, s2}, Q, 0.0);
  auto res_mm = sv::mixer_from_streams_and_jacobians({s1_mm, s2}, Q, 0.0);
  double fd_dmdot1 = (res_mp.T_mix - res_mm.T_mix) / (2.0 * eps_m);

  EXPECT_NEAR(res.dT_mix_d_stream[0].d_mdot, fd_dmdot1,
              std::max(std::abs(fd_dmdot1) * 1e-4, 1e-6));

  // FD: dT_mix/dmdot_2
  sv::Stream s2_mp = s2;
  s2_mp.m_dot += eps_m;
  sv::Stream s2_mm = s2;
  s2_mm.m_dot -= eps_m;
  auto res2_mp = sv::mixer_from_streams_and_jacobians({s1, s2_mp}, Q, 0.0);
  auto res2_mm = sv::mixer_from_streams_and_jacobians({s1, s2_mm}, Q, 0.0);
  double fd_dmdot2 = (res2_mp.T_mix - res2_mm.T_mix) / (2.0 * eps_m);

  EXPECT_NEAR(res.dT_mix_d_stream[1].d_mdot, fd_dmdot2,
              std::max(std::abs(fd_dmdot2) * 1e-4, 1e-6));
}

// =============================================================================
// 4. Mixer Fraction Jacobians: FD verification of all derivatives
// =============================================================================

TEST(EnergyBoundaryTest, MixerJacobian_Fraction_dT_dmdot) {
  // Verify dT/dmdot includes fraction correction
  sv::Stream s1;
  s1.m_dot = 0.4;
  s1.T = 400.0;
  s1.P_total = 150000.0;
  s1.Y = air_Y();

  sv::Stream s2;
  s2.m_dot = 0.6;
  s2.T = 800.0;
  s2.P_total = 150000.0;
  s2.Y = air_Y();

  double fraction = -0.03; // -3%
  auto res = sv::mixer_from_streams_and_jacobians({s1, s2}, 0.0, fraction);

  // FD: dT_mix/dmdot_1
  double eps_m = 1e-5;
  sv::Stream s1_mp = s1;
  s1_mp.m_dot += eps_m;
  sv::Stream s1_mm = s1;
  s1_mm.m_dot -= eps_m;
  auto res_mp = sv::mixer_from_streams_and_jacobians({s1_mp, s2}, 0.0, fraction);
  auto res_mm = sv::mixer_from_streams_and_jacobians({s1_mm, s2}, 0.0, fraction);
  double fd_dmdot1 = (res_mp.T_mix - res_mm.T_mix) / (2.0 * eps_m);

  EXPECT_NEAR(res.dT_mix_d_stream[0].d_mdot, fd_dmdot1,
              std::max(std::abs(fd_dmdot1) * 1e-4, 1e-6));
}

TEST(EnergyBoundaryTest, MixerJacobian_Fraction_dT_dT) {
  // Verify dT/dT_i includes fraction correction
  sv::Stream s1;
  s1.m_dot = 0.5;
  s1.T = 500.0;
  s1.P_total = 150000.0;
  s1.Y = air_Y();

  sv::Stream s2;
  s2.m_dot = 0.5;
  s2.T = 700.0;
  s2.P_total = 150000.0;
  s2.Y = air_Y();

  double fraction = 0.05; // +5%
  auto res = sv::mixer_from_streams_and_jacobians({s1, s2}, 0.0, fraction);

  // FD: dT_mix/dT_1
  double eps_T = 0.01;
  sv::Stream s1_Tp = s1;
  s1_Tp.T += eps_T;
  sv::Stream s1_Tm = s1;
  s1_Tm.T -= eps_T;
  auto res_Tp = sv::mixer_from_streams_and_jacobians({s1_Tp, s2}, 0.0, fraction);
  auto res_Tm = sv::mixer_from_streams_and_jacobians({s1_Tm, s2}, 0.0, fraction);
  double fd_dT1 = (res_Tp.T_mix - res_Tm.T_mix) / (2.0 * eps_T);

  EXPECT_NEAR(res.dT_mix_d_stream[0].d_T, fd_dT1,
              std::max(std::abs(fd_dT1) * 1e-4, 1e-6));

  // FD: dT_mix/dT_2
  sv::Stream s2_Tp = s2;
  s2_Tp.T += eps_T;
  sv::Stream s2_Tm = s2;
  s2_Tm.T -= eps_T;
  auto res2_Tp = sv::mixer_from_streams_and_jacobians({s1, s2_Tp}, 0.0, fraction);
  auto res2_Tm = sv::mixer_from_streams_and_jacobians({s1, s2_Tm}, 0.0, fraction);
  double fd_dT2 = (res2_Tp.T_mix - res2_Tm.T_mix) / (2.0 * eps_T);

  EXPECT_NEAR(res.dT_mix_d_stream[1].d_T, fd_dT2,
              std::max(std::abs(fd_dT2) * 1e-4, 1e-6));
}

TEST(EnergyBoundaryTest, MixerJacobian_Fraction_dT_dY) {
  // Verify dT/dY_i includes fraction correction
  sv::Stream st;
  st.m_dot = 1.0;
  st.T = 600.0;
  st.P_total = 150000.0;
  st.Y = air_Y();

  double fraction = -0.04; // -4%
  auto res = sv::mixer_from_streams_and_jacobians({st}, 0.0, fraction);

  // FD: dT_mix/dY_N2 (perturb N2 mass fraction)
  std::size_t ns = num_species();
  std::size_t k_N2 = 0;
  for (std::size_t k = 0; k < ns; ++k) {
    if (species_name(k) == "N2") {
      k_N2 = k;
      break;
    }
  }

  double eps_Y = 1e-5;
  sv::Stream st_Yp = st;
  st_Yp.Y[k_N2] += eps_Y;
  sv::Stream st_Ym = st;
  st_Ym.Y[k_N2] -= eps_Y;
  auto res_Yp = sv::mixer_from_streams_and_jacobians({st_Yp}, 0.0, fraction);
  auto res_Ym = sv::mixer_from_streams_and_jacobians({st_Ym}, 0.0, fraction);
  double fd_dY = (res_Yp.T_mix - res_Ym.T_mix) / (2.0 * eps_Y);

  EXPECT_NEAR(res.dT_mix_d_stream[0].d_Y[k_N2], fd_dY,
              std::max(std::abs(fd_dY) * 0.05, 1.0));
}

// =============================================================================
// 5. Complete combustion with Q and fraction
// =============================================================================

TEST(EnergyBoundaryTest, CompleteCombustion_Q_Cooling) {
  // Combust CH4 + air, then apply heat loss via Q
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
  double mdot_tot = fuel.m_dot + ox.m_dot;

  auto res_base = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, 0.0, 0.0);

  // With heat loss Q = -200 kW
  double Q = -200000.0;
  auto res_cooled = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, Q, 0.0);

  // Cooled T should be lower
  EXPECT_LT(res_cooled.T_mix, res_base.T_mix);

  // Verify enthalpy balance
  std::vector<double> X_b = mass_to_mole(normalize_fractions(res_base.Y_mix));
  double h_ad = h_mass(res_base.T_mix, X_b);
  double h_cooled = h_mass(res_cooled.T_mix, X_b);
  double delta_h = Q / mdot_tot;
  EXPECT_NEAR(h_cooled, h_ad + delta_h, std::abs(h_ad + delta_h) * 1e-5);

  // Products composition should be unchanged (complete combustion)
  for (std::size_t k = 0; k < res_base.Y_mix.size(); ++k) {
    EXPECT_NEAR(res_cooled.Y_mix[k], res_base.Y_mix[k], 1e-10);
  }
}

TEST(EnergyBoundaryTest, CompleteCombustion_Fraction_Cooling) {
  // Combust CH4 + air, then apply heat loss via fraction
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
  double fraction = -0.10; // -10% enthalpy loss

  auto res_base = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, 0.0, 0.0);
  auto res_cooled = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, 0.0, fraction);

  // Cooled T should be lower
  EXPECT_LT(res_cooled.T_mix, res_base.T_mix);

  // Enthalpy should be reduced by fraction
  std::vector<double> X_b = mass_to_mole(normalize_fractions(res_base.Y_mix));
  double h_ad = h_mass(res_base.T_mix, X_b);
  double h_cooled = h_mass(res_cooled.T_mix, X_b);
  EXPECT_NEAR(h_cooled / h_ad, 1.0 + fraction, 1e-4);
}

TEST(EnergyBoundaryTest, CompleteCombustion_Q_Zero_IsNoop) {
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
      {fuel, ox}, P, 0.0, 0.0);
  auto res_none = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P);

  EXPECT_DOUBLE_EQ(res_zero.T_mix, res_none.T_mix);
}

// =============================================================================
// 6. Complete combustion Jacobians: FD verification with Q and fraction
// =============================================================================

TEST(EnergyBoundaryTest, CompleteCombustion_Jacobian_Q_dT_dmdot) {
  // Verify combustion dT/dmdot includes Q correction
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
  double Q = -100000.0; // -100 kW

  auto res = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, Q, 0.0);

  // FD: dT/dmdot_fuel
  double eps_m = 1e-6;
  sv::Stream fuel_mp = fuel;
  fuel_mp.m_dot += eps_m;
  sv::Stream fuel_mm = fuel;
  fuel_mm.m_dot -= eps_m;
  auto res_mp = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel_mp, ox}, P, Q, 0.0);
  auto res_mm = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel_mm, ox}, P, Q, 0.0);
  double fd_dmdot = (res_mp.T_mix - res_mm.T_mix) / (2.0 * eps_m);

  EXPECT_NEAR(res.dT_mix_d_stream[0].d_mdot, fd_dmdot,
              std::max(std::abs(fd_dmdot) * 0.01, 1e-3));
}

TEST(EnergyBoundaryTest, CompleteCombustion_Jacobian_Fraction_dT_dT) {
  // Verify combustion dT/dT_i includes fraction correction
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
  double fraction = -0.05; // -5%

  auto res = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox}, P, 0.0, fraction);

  // FD: dT/dT_air
  double eps_T = 0.5;
  sv::Stream ox_Tp = ox;
  ox_Tp.T += eps_T;
  sv::Stream ox_Tm = ox;
  ox_Tm.T -= eps_T;
  auto res_Tp = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox_Tp}, P, 0.0, fraction);
  auto res_Tm = sv::adiabatic_T_complete_and_jacobian_T_from_streams(
      {fuel, ox_Tm}, P, 0.0, fraction);
  double fd_dT = (res_Tp.T_mix - res_Tm.T_mix) / (2.0 * eps_T);

  EXPECT_NEAR(res.dT_mix_d_stream[1].d_T, fd_dT,
              std::max(std::abs(fd_dT) * 0.06, 1e-4));
}

// =============================================================================
// 7. Equilibrium combustion with Q and fraction
// =============================================================================

TEST(EnergyBoundaryTest, EquilibriumCombustion_Q_Cooling) {
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
  double Q = -150000.0; // -150 kW

  auto res_base = sv::adiabatic_T_equilibrium_and_jacobians_from_streams(
      {fuel, ox}, P, 0.0, 0.0);
  auto res_cooled = sv::adiabatic_T_equilibrium_and_jacobians_from_streams(
      {fuel, ox}, P, Q, 0.0);

  EXPECT_LT(res_cooled.T_mix, res_base.T_mix);
}

TEST(EnergyBoundaryTest, EquilibriumCombustion_Jacobian_Fraction_dT_dmdot) {
  // Verify equilibrium combustion dT/dmdot with fraction
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
  double fraction = -0.08; // -8%

  auto res = sv::adiabatic_T_equilibrium_and_jacobians_from_streams(
      {fuel, ox}, P, 0.0, fraction);

  // FD: dT/dmdot_air
  double eps_m = 1e-6;
  sv::Stream ox_mp = ox;
  ox_mp.m_dot += eps_m;
  sv::Stream ox_mm = ox;
  ox_mm.m_dot -= eps_m;
  auto res_mp = sv::adiabatic_T_equilibrium_and_jacobians_from_streams(
      {fuel, ox_mp}, P, 0.0, fraction);
  auto res_mm = sv::adiabatic_T_equilibrium_and_jacobians_from_streams(
      {fuel, ox_mm}, P, 0.0, fraction);
  double fd_dmdot = (res_mp.T_mix - res_mm.T_mix) / (2.0 * eps_m);

  EXPECT_NEAR(res.dT_mix_d_stream[1].d_mdot, fd_dmdot,
              std::max(std::abs(fd_dmdot) * 0.015, 1e-3));
}

// =============================================================================
// 8. Stream mix() with delta_h (unchanged API)
// =============================================================================

TEST(EnergyBoundaryTest, StreamMix_DeltaH) {
  // Stream mix() keeps delta_h parameter (no Jacobians, simpler use case)
  Stream s1;
  s1.mdot = 0.5;
  s1.state.set_TPX(300.0, 101325.0, dry_air());

  Stream s2;
  s2.mdot = 0.5;
  s2.state.set_TPX(600.0, 101325.0, dry_air());

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
