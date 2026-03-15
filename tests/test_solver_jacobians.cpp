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

    OrificeResult res = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y, P_static_down, Cd, area);

    // FD checks
    double eps = 1.0;

    // d_mdot_dP_total_up
    auto res_p = orifice_residuals_and_jacobian(m_dot, P_total_up + eps, P_static_up, T_up, Y, P_static_down, Cd, area);
    double fd_dP_total_up = (res_p.m_dot_calc - res.m_dot_calc) / eps;
    EXPECT_NEAR(res.d_mdot_dP_total_up, fd_dP_total_up, std::abs(fd_dP_total_up) * 1e-4);

    // d_mdot_dP_static_down
    auto res_down_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y, P_static_down + eps, Cd, area);
    double fd_dP_static_down = (res_down_p.m_dot_calc - res.m_dot_calc) / eps;
    EXPECT_NEAR(res.d_mdot_dP_static_down, fd_dP_static_down, std::abs(fd_dP_static_down) * 1e-4);

    // d_mdot_dP_static_up (via density)
    auto res_up_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up + eps, T_up, Y, P_static_down, Cd, area);
    double fd_dP_static_up = (res_up_p.m_dot_calc - res.m_dot_calc) / eps;
    EXPECT_NEAR(res.d_mdot_dP_static_up, fd_dP_static_up, std::abs(fd_dP_static_up) * 1e-4);

    // d_mdot_dT_up
    double eps_T = 1e-4;
    auto res_T_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up + eps_T, Y, P_static_down, Cd, area);
    auto res_T_m = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up - eps_T, Y, P_static_down, Cd, area);
    double fd_dT_up = (res_T_p.m_dot_calc - res_T_m.m_dot_calc) / (2.0 * eps_T);
    EXPECT_NEAR(res.d_mdot_dT_up, fd_dT_up, std::abs(fd_dT_up) * 1e-5);

    // Y derivatives
    double eps_Y = 1e-6;
    for (size_t i=0; i<ns; ++i) {
        std::vector<double> Y_p = Y;
        Y_p[i] += eps_Y;
        auto res_Y_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y_p, P_static_down, Cd, area);
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
