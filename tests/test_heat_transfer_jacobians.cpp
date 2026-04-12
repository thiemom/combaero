#include "../include/heat_transfer.h"
#include "../include/humidair.h"
#include "../include/thermo.h"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

using namespace combaero;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double cross_section_area(double diameter) {
  return 0.25 * kPi * diameter * diameter;
}

double mass_flow(double rho, double velocity, double diameter) {
  return rho * velocity * cross_section_area(diameter);
}

double velocity_from_mass_flow(double mdot, double rho, double diameter) {
  return mdot / (rho * cross_section_area(diameter));
}

} // namespace

TEST(HeatTransferJacobiansTest, ChannelSmoothFiniteDifference) {
  const auto X = dry_air();

  const double T = 800.0;
  const double P = 2.0e5;
  const double diameter = 0.04;
  const double length = 1.0;
  const double velocity = 25.0;
  const double T_wall = 900.0;

  auto evaluate = [&](double T_eval, double velocity_eval, double T_wall_eval) {
    return channel_smooth(
        T_eval,
        P,
        X,
        velocity_eval,
        diameter,
        length,
        T_wall_eval,
        "gnielinski",
        true,
        1.0,
        0.0,
        1.0,
        1.0);
  };

  ChannelResult base = evaluate(T, velocity, T_wall);
  ASSERT_TRUE(std::isfinite(base.h));
  ASSERT_TRUE(std::isfinite(base.q));

  const double area = cross_section_area(diameter);
  const double rho = density(T, P, X);
  const double mdot = mass_flow(rho, velocity, diameter);

  // --- Mass-flow derivatives ---
  const double eps_mdot = 1e-4 * mdot;
  const double delta_v = eps_mdot / (rho * area);
  ChannelResult plus_m = evaluate(T, velocity + delta_v, T_wall);
  ChannelResult minus_m = evaluate(T, velocity - delta_v, T_wall);

  auto fd_mdot = [&](double ChannelResult::*member) {
    return (plus_m.*member - minus_m.*member) / (2.0 * eps_mdot);
  };

  double fd_dh_dmdot = fd_mdot(&ChannelResult::h);
  EXPECT_NEAR(base.dh_dmdot, fd_dh_dmdot,
              std::max(1e-6, std::abs(fd_dh_dmdot) * 1e-3));

  double fd_ddP_dmdot = fd_mdot(&ChannelResult::dP);
  EXPECT_NEAR(base.ddP_dmdot, fd_ddP_dmdot,
              std::max(1e-6, std::abs(fd_ddP_dmdot) * 1e-3));

  double fd_dq_dmdot = fd_mdot(&ChannelResult::q);
  EXPECT_NEAR(base.dq_dmdot, fd_dq_dmdot,
              std::max(1e-6, std::abs(fd_dq_dmdot) * 1e-3));

  double fd_dT_aw_dmdot = fd_mdot(&ChannelResult::T_aw);
  EXPECT_NEAR(base.dT_aw_dmdot, fd_dT_aw_dmdot,
              std::max(1e-6, std::abs(fd_dT_aw_dmdot) * 1e-3));

  // --- Temperature derivatives (constant mass flow) ---
  const double eps_T = 0.5;
  auto eval_at_T = [&](double T_eval) {
    double rho_eval = density(T_eval, P, X);
    double velocity_eval = velocity_from_mass_flow(mdot, rho_eval, diameter);
    return evaluate(T_eval, velocity_eval, T_wall);
  };

  ChannelResult plus_T = eval_at_T(T + eps_T);
  ChannelResult minus_T = eval_at_T(T - eps_T);

  auto fd_T = [&](double ChannelResult::*member) {
    return (plus_T.*member - minus_T.*member) / (2.0 * eps_T);
  };

  double fd_dh_dT = fd_T(&ChannelResult::h);
  EXPECT_NEAR(base.dh_dT, fd_dh_dT,
              std::max(1e-6, std::abs(fd_dh_dT) * 1e-3));

  double fd_ddP_dT = fd_T(&ChannelResult::dP);
  EXPECT_NEAR(base.ddP_dT, fd_ddP_dT,
              std::max(1e-6, std::abs(fd_ddP_dT) * 1e-3));

  double fd_dq_dT = fd_T(&ChannelResult::q);
  EXPECT_NEAR(base.dq_dT, fd_dq_dT,
              std::max(1e-6, std::abs(fd_dq_dT) * 1e-3));

  double fd_dT_aw_dT = fd_T(&ChannelResult::T_aw);
  EXPECT_NEAR(base.dT_aw_dT, fd_dT_aw_dT,
              std::max(1e-6, std::abs(fd_dT_aw_dT) * 1e-3));

  // --- Wall temperature derivative ---
  const double eps_Twall = 0.5;
  ChannelResult plus_Tw = evaluate(T, velocity, T_wall + eps_Twall);
  ChannelResult minus_Tw = evaluate(T, velocity, T_wall - eps_Twall);
  double fd_dq_dT_wall = (plus_Tw.q - minus_Tw.q) / (2.0 * eps_Twall);
  EXPECT_NEAR(base.dq_dT_wall, fd_dq_dT_wall,
              std::max(1e-6, std::abs(fd_dq_dT_wall) * 1e-3));
}

TEST(HeatTransferJacobiansTest, WallCouplingFiniteDifference) {
  const double h_a = 1500.0;
  const double h_b = 2200.0;
  const double T_aw_a = 850.0;
  const double T_aw_b = 500.0;
  const double t_over_k = 0.002 / 20.0;
  const double A = 0.12;

  auto evaluate = [&](double hA, double hB, double TawA, double TawB) {
    return wall_coupling_and_jacobian(hA, TawA, hB, TawB, t_over_k, A);
  };

  WallCouplingResult base = evaluate(h_a, h_b, T_aw_a, T_aw_b);

  auto fd = [&](auto perturb, double base_value, double eps) {
    auto plus = perturb(base_value + eps);
    auto minus = perturb(base_value - eps);
    return (plus.Q - minus.Q) / (2.0 * eps);
  };

  const double eps_h = 1e-4 * h_a;
  auto perturb_h_a = [&](double value) {
    return evaluate(value, h_b, T_aw_a, T_aw_b);
  };
  double fd_dQ_dh_a = fd(perturb_h_a, h_a, eps_h);
  EXPECT_NEAR(base.dQ_dh_a, fd_dQ_dh_a,
              std::max(1e-6, std::abs(fd_dQ_dh_a) * 1e-6));

  auto perturb_h_b = [&](double value) {
    return evaluate(h_a, value, T_aw_a, T_aw_b);
  };
  double fd_dQ_dh_b = fd(perturb_h_b, h_b, eps_h);
  EXPECT_NEAR(base.dQ_dh_b, fd_dQ_dh_b,
              std::max(1e-6, std::abs(fd_dQ_dh_b) * 1e-6));

  const double eps_Taw = 1e-1;
  auto perturb_Taw_a = [&](double value) {
    return evaluate(h_a, h_b, value, T_aw_b);
  };
  double fd_dQ_dT_aw_a = fd(perturb_Taw_a, T_aw_a, eps_Taw);
  EXPECT_NEAR(base.dQ_dT_aw_a, fd_dQ_dT_aw_a,
              std::max(1e-6, std::abs(fd_dQ_dT_aw_a) * 1e-6));

  auto perturb_Taw_b = [&](double value) {
    return evaluate(h_a, h_b, T_aw_a, value);
  };
  double fd_dQ_dT_aw_b = fd(perturb_Taw_b, T_aw_b, eps_Taw);
  EXPECT_NEAR(base.dQ_dT_aw_b, fd_dQ_dT_aw_b,
              std::max(1e-6, std::abs(fd_dQ_dT_aw_b) * 1e-6));
}

TEST(HeatTransferJacobiansTest, WallCouplingMultiLayerFiniteDifference) {
  const double h_a = 1500.0;
  const double h_b = 2200.0;
  const double T_aw_a = 850.0;
  const double T_aw_b = 500.0;
  const std::vector<double> t_over_k = {0.001 / 50.0, 0.002 / 1.0}; // Steel + Insulation
  const double R_fouling = 0.0001;
  const double A = 0.12;

  auto evaluate = [&](double hA, double hB, double TawA, double TawB) {
    return wall_coupling_and_jacobian(hA, TawA, hB, TawB, t_over_k, A, R_fouling);
  };

  WallCouplingResult base = evaluate(h_a, h_b, T_aw_a, T_aw_b);

  auto fd = [&](auto perturb, double base_value, double eps) {
    auto plus = perturb(base_value + eps);
    auto minus = perturb(base_value - eps);
    return (plus.Q - minus.Q) / (2.0 * eps);
  };

  // dQ/dh_a
  const double eps_h = 1e-4 * h_a;
  auto perturb_h_a = [&](double val) {
    return evaluate(val, h_b, T_aw_a, T_aw_b);
  };
  double fd_dQ_dh_a = fd(perturb_h_a, h_a, eps_h);
  EXPECT_NEAR(base.dQ_dh_a, fd_dQ_dh_a,
              std::max(1e-6, std::abs(fd_dQ_dh_a) * 1e-6));

  // dQ/dh_b
  auto perturb_h_b = [&](double val) {
    return evaluate(h_a, val, T_aw_a, T_aw_b);
  };
  double fd_dQ_dh_b = fd(perturb_h_b, h_b, 1e-4 * h_b);
  EXPECT_NEAR(base.dQ_dh_b, fd_dQ_dh_b,
              std::max(1e-6, std::abs(fd_dQ_dh_b) * 1e-6));

  // dQ/dT_aw_a
  const double eps_T = 0.1;
  auto perturb_Taw_a = [&](double val) {
    return evaluate(h_a, h_b, val, T_aw_b);
  };
  double fd_dQ_dT_aw_a = fd(perturb_Taw_a, T_aw_a, eps_T);
  EXPECT_NEAR(base.dQ_dT_aw_a, fd_dQ_dT_aw_a,
              std::max(1e-6, std::abs(fd_dQ_dT_aw_a) * 1e-6));

  // dQ/dT_aw_b
  auto perturb_Taw_b = [&](double val) {
    return evaluate(h_a, h_b, T_aw_a, val);
  };
  double fd_dQ_dT_aw_b = fd(perturb_Taw_b, T_aw_b, eps_T);
  EXPECT_NEAR(base.dQ_dT_aw_b, fd_dQ_dT_aw_b,
              std::max(1e-6, std::abs(fd_dQ_dT_aw_b) * 1e-6));
}

TEST(HeatTransferJacobiansTest, WallTemperatureProfileTest) {
  const double T_hot = 1000.0;
  const double T_cold = 300.0;
  const double h_hot = 500.0;
  const double h_cold = 50.0;
  const std::vector<double> t_over_k = {0.002 / 50.0, 0.010 / 0.1}; // Steel + Insulation
  const double R_fouling = 0.002;

  double q_calc;
  std::vector<double> profile = wall_temperature_profile(T_hot, T_cold, h_hot, h_cold, t_over_k, R_fouling, q_calc);

  // Profile points: [T_surf_hot, T_interface, T_surf_cold]
  ASSERT_EQ(profile.size(), 3);

  // Monotonicity
  EXPECT_LT(profile[0], T_hot);
  EXPECT_GT(profile[0], profile[1]);
  EXPECT_GT(profile[1], profile[2]);
  EXPECT_GT(profile[2], T_cold);

  // Check consistency: Q = h_hot * (T_hot - T_surf_hot)
  double q_hot = h_hot * (T_hot - profile[0]);
  EXPECT_NEAR(q_calc, q_hot, 1e-7);

  // Check interface drop: delta T = q * R_layer1
  double dT_layer1 = profile[0] - profile[1];
  EXPECT_NEAR(dT_layer1, q_calc * t_over_k[0], 1e-7);

  // Check cold surface consistency: Q = (T_surf_cold - T_cold) / (R_fouling + 1/h_cold)
  double q_cold = (profile[2] - T_cold) / (R_fouling + 1.0 / h_cold);
  EXPECT_NEAR(q_calc, q_cold, 1e-7);
}
