#include <gtest/gtest.h>

#include "../include/state.h"
#include <cmath>

// Helper to check fractional difference
#define EXPECT_REL_ERR(val, ref, tol)                                          \
  EXPECT_LT(std::abs((val) - (ref)) / std::abs(ref), tol)

TEST(StateTest, SettersperfectlyRecoverConditions) {
  State s;
  std::vector<double> X_air = {0.78084, 0.209476, 0.009365, 0.000319, 0.0,
                               0.0,     0.0,      0.0,      0.0,      0.0,
                               0.0,     0.0,      0.0,      0.0};

  double T_target = 800.0;
  double P_target = 2000000.0; // 2 MPa

  s.set_TPX(T_target, P_target, X_air);

  double h_mass = s.h() / s.mw() * 1000.0;
  double s_mass = s.s() / s.mw() * 1000.0;
  double u_mass = s.u() / s.mw() * 1000.0;
  double v_mass = 1.0 / s.rho();

  {
    State s2;
    s2.set_TPX(300, 101325, X_air);
    s2.set_HP_mass(h_mass, P_target);
    EXPECT_REL_ERR(s2.T, T_target, 1e-6);
    EXPECT_REL_ERR(s2.P, P_target, 1e-6);
  }

  {
    State s2;
    s2.set_TPX(300, 101325, X_air);
    s2.set_SP_mass(s_mass, P_target);
    EXPECT_REL_ERR(s2.T, T_target, 1e-6);
    EXPECT_REL_ERR(s2.P, P_target, 1e-6);
  }

  {
    State s2;
    s2.set_TPX(300, 101325, X_air);
    s2.set_UV_mass(u_mass, v_mass);
    EXPECT_REL_ERR(s2.T, T_target, 1e-6);
    EXPECT_REL_ERR(s2.P, P_target, 1e-6);
  }

  {
    State s2;
    s2.set_TPX(300, 101325, X_air);
    s2.set_SV_mass(s_mass, v_mass);
    EXPECT_REL_ERR(s2.T, T_target, 1e-6);
    EXPECT_REL_ERR(s2.P, P_target, 1e-6);
  }

  {
    State s2;
    s2.set_TPX(300, 101325, X_air);
    s2.set_SH_mass(s_mass, h_mass);
    EXPECT_REL_ERR(s2.T, T_target, 1e-6);
    EXPECT_REL_ERR(s2.P, P_target, 1e-6);
  }
}
