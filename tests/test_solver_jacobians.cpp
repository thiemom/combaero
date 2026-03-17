#include <gtest/gtest.h>
#include "solver_interface.h"
#include "thermo.h"
#include <vector>
#include <cmath>

using namespace combaero;
using namespace combaero::solver;

TEST(SolverJacobianTest, OrificeDerivatives) {
    // Boilerplate for species
    std::size_t ns = num_species();
    std::vector<double> Y(ns, 0.0);
    // Assume air-like for testing
    for(size_t i=0; i<ns; ++i) {
        if (species_name(i) == "N2") Y[i] = 0.77;
        else if (species_name(i) == "O2") Y[i] = 0.23;
    }

    double m_dot = 0.1;
    double P_total_up = 200000.0;
    double P_static_up = 190000.0;
    double T_up = 300.0;
    double P_static_down = 180000.0;
    double Cd = 0.62;
    double area = 0.001;

    double beta = 0.5;
    OrificeResult res = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y, P_static_down, Cd, area, beta);

    // FD checks
    double eps = 1.0;

    // d_mdot_dP_total_up
    auto res_p = orifice_residuals_and_jacobian(m_dot, P_total_up + eps, P_static_up, T_up, Y, P_static_down, Cd, area, beta);
    double fd_dP_total_up = (res_p.m_dot_calc - res.m_dot_calc) / eps;
    EXPECT_NEAR(res.d_mdot_dP_total_up, fd_dP_total_up, std::abs(fd_dP_total_up) * 1e-4);

    // d_mdot_dP_static_down
    auto res_down_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y, P_static_down + eps, Cd, area, beta);
    double fd_dP_static_down = (res_down_p.m_dot_calc - res.m_dot_calc) / eps;
    EXPECT_NEAR(res.d_mdot_dP_static_down, fd_dP_static_down, std::abs(fd_dP_static_down) * 1e-4);

    // d_mdot_dP_static_up (via density)
    auto res_up_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up + eps, T_up, Y, P_static_down, Cd, area, beta);
    double fd_dP_static_up = (res_up_p.m_dot_calc - res.m_dot_calc) / eps;
    EXPECT_NEAR(res.d_mdot_dP_static_up, fd_dP_static_up, std::abs(fd_dP_static_up) * 1e-4);

    // d_mdot_dT_up
    double eps_T = 1e-4;
    auto res_T_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up + eps_T, Y, P_static_down, Cd, area, beta);
    auto res_T_m = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up - eps_T, Y, P_static_down, Cd, area, beta);
    double fd_dT_up = (res_T_p.m_dot_calc - res_T_m.m_dot_calc) / (2.0 * eps_T);
    EXPECT_NEAR(res.d_mdot_dT_up, fd_dT_up, std::abs(fd_dT_up) * 1e-5);

    // Y derivatives
    double eps_Y = 1e-6;
    for (size_t i=0; i<ns; ++i) {
        std::vector<double> Y_p = Y;
        Y_p[i] += eps_Y;
        auto res_Y_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y_p, P_static_down, Cd, area, beta);
        double fd_dYi = (res_Y_p.m_dot_calc - res.m_dot_calc) / eps_Y;
        EXPECT_NEAR(res.d_mdot_dY_up[i], fd_dYi, std::abs(fd_dYi) * 1e-3 + 1e-8);
    }
}

TEST(SolverJacobianTest, PipeDerivatives) {
    std::size_t ns = num_species();
    std::vector<double> Y(ns, 0.0);
    for(size_t i=0; i<ns; ++i) {
        if (species_name(i) == "N2") Y[i] = 0.77;
        else if (species_name(i) == "O2") Y[i] = 0.23;
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

    PipeResult res = pipe_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y, P_static_down, L, D, roughness, model);

    double eps_m = 1e-4;
    auto res_m_p = pipe_residuals_and_jacobian(m_dot + eps_m, P_total_up, P_static_up, T_up, Y, P_static_down, L, D, roughness, model);
    double fd_dmdot = (res_m_p.dP_calc - res.dP_calc) / eps_m;
    EXPECT_NEAR(res.d_dP_d_mdot, fd_dmdot, std::abs(fd_dmdot) * 1e-4);

    double eps_P = 1.0;
    auto res_P_p = pipe_residuals_and_jacobian(m_dot, P_total_up, P_static_up + eps_P, T_up, Y, P_static_down, L, D, roughness, model);
    double fd_dP_static = (res_P_p.dP_calc - res.dP_calc) / eps_P;
    EXPECT_NEAR(res.d_dP_dP_static_up, fd_dP_static, std::abs(fd_dP_static) * 1e-4);

    double eps_T = 0.1;
    auto res_T_p = pipe_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up + eps_T, Y, P_static_down, L, D, roughness, model);
    double fd_dT = (res_T_p.dP_calc - res.dP_calc) / eps_T;
    EXPECT_NEAR(res.d_dP_dT_up, fd_dT, std::abs(fd_dT) * 1e-4);

    double eps_Y = 1e-6;
    for (size_t i=0; i<ns; ++i) {
        std::vector<double> Y_p = Y;
        Y_p[i] += eps_Y;
        auto res_Y_p = pipe_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y_p, P_static_down, L, D, roughness, model);
        double fd_dYi = (res_Y_p.dP_calc - res.dP_calc) / eps_Y;
        EXPECT_NEAR(res.d_dP_dY_up[i], fd_dYi, std::abs(fd_dYi) * 1e-2 + 1e-8);
    }
}

TEST(SolverJacobianTest, MomentumChamberDerivatives) {
    // Test momentum chamber residual: P_total = P + 0.5 * rho * v^2
    std::size_t ns = num_species();
    std::vector<double> Y(ns, 0.0);
    for(size_t i=0; i<ns; ++i) {
        Y[i] = (i == species_index_from_name("N2")) ? 0.767 :
               (i == species_index_from_name("O2")) ? 0.233 : 0.0;
    }

    double P = 101325.0;
    double P_total = 105000.0;
    double m_dot = 0.5;
    double T = 600.0;
    double area = 0.01;  // 0.01 m^2

    MomentumChamberResult res = momentum_chamber_residual_and_jacobian(P, P_total, m_dot, T, Y, area);

    // Finite difference checks
    double eps_P = 1.0;
    auto res_P_p = momentum_chamber_residual_and_jacobian(P + eps_P, P_total, m_dot, T, Y, area);
    double fd_dP = (res_P_p.residual - res.residual) / eps_P;
    EXPECT_NEAR(res.d_res_dP, fd_dP, std::abs(fd_dP) * 1e-4 + 1e-8);

    auto res_Pt_p = momentum_chamber_residual_and_jacobian(P, P_total + eps_P, m_dot, T, Y, area);
    double fd_dP_total = (res_Pt_p.residual - res.residual) / eps_P;
    EXPECT_NEAR(res.d_res_dP_total, fd_dP_total, std::abs(fd_dP_total) * 1e-4 + 1e-8);

    double eps_mdot = 1e-4;
    auto res_mdot_p = momentum_chamber_residual_and_jacobian(P, P_total, m_dot + eps_mdot, T, Y, area);
    auto res_mdot_m = momentum_chamber_residual_and_jacobian(P, P_total, m_dot - eps_mdot, T, Y, area);
    double fd_dmdot = (res_mdot_p.residual - res_mdot_m.residual) / (2.0 * eps_mdot);
    EXPECT_NEAR(res.d_res_dmdot, fd_dmdot, std::abs(fd_dmdot) * 1e-4 + 1e-8);

    // Verify residual is correct: P_total - P - 0.5*rho*v^2 = 0
    std::vector<double> X = mass_to_mole(normalize_fractions(Y));
    double rho = density(T, P, X);
    double v = m_dot / (rho * area);
    double q_dynamic = 0.5 * rho * v * v;
    double expected_residual = P_total - P - q_dynamic;
    EXPECT_NEAR(res.residual, expected_residual, 1e-6);
}
