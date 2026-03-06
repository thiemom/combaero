from __future__ import annotations

import numpy as np

import combaero as cb


def central_difference(func, x, h=1e-6, *args):
    """Numerically evaluate df/dx using O(h^2) central difference."""
    return (func(x + h, *args) - func(x - h, *args)) / (2.0 * h)


def test_orifice_jacobian():
    rho = 1.225
    Cd = 0.6
    area = 0.05
    dP_test_points = [10.0, 100.0, 50000.0]

    # We define a pure python wrapper to extract just the value from our fast-path API
    def mdot_func(dP):
        return cb._core.orifice_mdot_and_jacobian(dP, rho, Cd, area)[0]

    for dP in dP_test_points:
        # C++ Analytical Tuple
        mdot_analytic, jac_analytic = cb._core.orifice_mdot_and_jacobian(dP, rho, Cd, area)

        # Python Numerical Gradient
        h = max(1e-5, dP * 1e-6)
        jac_numeric = central_difference(mdot_func, dP, h)

        assert mdot_analytic > 0.0
        assert jac_analytic > 0.0
        # The analytical derivative should perfectly match the central difference
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-6)


def test_pressure_loss_jacobian():
    rho = 10.0
    K = 1.5
    v_test_points = [0.1, 10.0, 300.0]

    def dP_func(v):
        return cb._core.pressure_loss_and_jacobian(v, rho, K)[0]

    for v in v_test_points:
        dP_analytic, jac_analytic = cb._core.pressure_loss_and_jacobian(v, rho, K)
        jac_numeric = central_difference(dP_func, v)

        assert dP_analytic > 0.0
        assert jac_analytic > 0.0
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-6)


def test_nusselt_jacobian():
    Pr_heating = 0.7  # Air
    Pr_cooling = 4.3  # Water

    def nu_func(Re, Pr, heating):
        return cb._core.nusselt_and_jacobian_dittus_boelter(Re, Pr, heating).result[0]

    for Re in [10000.0, 50000.0, 1e6]:
        # Heating
        Nu_analytic, jac_analytic = cb._core.nusselt_and_jacobian_dittus_boelter(
            Re, Pr_heating, True
        ).result
        jac_numeric = central_difference(
            nu_func, Re, 1e-4, Pr_heating, True
        )  # slightly larger step for large Re
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)

        # Cooling
        Nu_analytic, jac_analytic = cb._core.nusselt_and_jacobian_dittus_boelter(
            Re, Pr_cooling, False
        ).result
        jac_numeric = central_difference(nu_func, Re, 1e-4, Pr_cooling, False)
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)


def test_friction_haaland_jacobian():
    e_D = 1e-4

    def f_func(Re):
        return cb._core.friction_and_jacobian_haaland(Re, e_D).result[0]

    for Re in [3000.0, 10000.0, 1e6]:
        f_analytic, jac_analytic = cb._core.friction_and_jacobian_haaland(Re, e_D).result

        # Friction factor gradients are very small, adjust h scalar based on Re magnitude
        h = max(1e-4, Re * 1e-6)
        jac_numeric = central_difference(f_func, Re, h)

        assert f_analytic > 0.0
        # Friction factor decreases with Reynolds number
        assert jac_analytic < 0.0
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)


def test_density_jacobians():
    P = 101325.0
    X = [0.0] * 14
    X[11] = 0.79  # N2
    X[12] = 0.21  # O2

    def rho_T_func(T):
        return cb._core.density_and_jacobians(T, P, X)[0]

    def rho_P_func(P_val):
        return cb._core.density_and_jacobians(300.0, P_val, X)[0]

    for T in [300.0, 1200.0]:
        rho_analytic, d_rho_dT_analytic, d_rho_dP_analytic = cb._core.density_and_jacobians(T, P, X)
        d_rho_dT_numeric = central_difference(rho_T_func, T, 1e-3)
        np.testing.assert_allclose(d_rho_dT_analytic, d_rho_dT_numeric, rtol=1e-5)

    for P_val in [101325.0, 500000.0]:
        rho_analytic, d_rho_dT_analytic, d_rho_dP_analytic = cb._core.density_and_jacobians(
            300.0, P_val, X
        )
        d_rho_dP_numeric = central_difference(rho_P_func, P_val, 1e-1)
        np.testing.assert_allclose(d_rho_dP_analytic, d_rho_dP_numeric, rtol=1e-5, atol=1e-12)


def test_enthalpy_jacobian():
    X = [0.0] * 14
    X[11] = 0.79  # N2
    X[12] = 0.21  # O2

    def h_func(T):
        return cb._core.enthalpy_and_jacobian(T, X)[0]

    for T in [300.0, 1200.0]:
        h_analytic, d_h_dT_analytic = cb._core.enthalpy_and_jacobian(T, X)
        d_h_dT_numeric = central_difference(h_func, T, 1e-2)
        np.testing.assert_allclose(d_h_dT_analytic, d_h_dT_numeric, rtol=1e-5)


def test_viscosity_jacobians():
    P = 101325.0
    X = [0.0] * 14
    X[11] = 0.79  # N2
    X[12] = 0.21  # O2

    def mu_T_func(T):
        return cb._core.viscosity_and_jacobians(T, P, X)[0]

    def mu_P_func(P_val):
        return cb._core.viscosity_and_jacobians(300.0, P_val, X)[0]

    for T in [300.0, 1200.0]:
        mu_analytic, d_mu_dT_analytic, d_mu_dP_analytic = cb._core.viscosity_and_jacobians(T, P, X)
        d_mu_dT_numeric = central_difference(mu_T_func, T, 1e-3)
        np.testing.assert_allclose(d_mu_dT_analytic, d_mu_dT_numeric, rtol=1e-5)

    for P_val in [101325.0, 500000.0]:
        mu_analytic, d_mu_dT_analytic, d_mu_dP_analytic = cb._core.viscosity_and_jacobians(
            300.0, P_val, X
        )
        d_mu_dP_numeric = central_difference(mu_P_func, P_val, 1e-1)
        np.testing.assert_allclose(d_mu_dP_analytic, d_mu_dP_numeric, rtol=1e-5, atol=1e-12)


def test_nusselt_gnielinski_jacobian():
    Pr_val = 0.7
    f_val = 0.05

    def nu_func(Re):
        return cb._core.nusselt_and_jacobian_gnielinski(Re, Pr_val, f_val).result[0]

    for Re in [3000.0, 10000.0, 1e6]:
        Nu_analytic, jac_analytic = cb._core.nusselt_and_jacobian_gnielinski(
            Re, Pr_val, f_val
        ).result
        jac_numeric = central_difference(nu_func, Re, 1e-4)
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)


def test_nusselt_sieder_tate_jacobian():
    Pr_val = 0.7
    mu_ratio = 1.1

    def nu_func(Re):
        return cb._core.nusselt_and_jacobian_sieder_tate(Re, Pr_val, mu_ratio).result[0]

    for Re in [3000.0, 10000.0, 1e6]:
        Nu_analytic, jac_analytic = cb._core.nusselt_and_jacobian_sieder_tate(
            Re, Pr_val, mu_ratio
        ).result
        jac_numeric = central_difference(nu_func, Re, 1e-4)
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)


def test_nusselt_petukhov_jacobian():
    Pr_val = 0.7
    f_val = 0.05

    def nu_func(Re):
        return cb._core.nusselt_and_jacobian_petukhov(Re, Pr_val, f_val).result[0]

    for Re in [10000.0, 50000.0, 1e6]:
        Nu_analytic, jac_analytic = cb._core.nusselt_and_jacobian_petukhov(Re, Pr_val, f_val).result
        jac_numeric = central_difference(nu_func, Re, 1e-4)
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)


def test_pin_fin_jacobians():
    Pr = 0.71
    L_D = 1.5
    S_D = 2.5
    X_D = 2.5

    def nu_func(Re_d):
        return cb._core.pin_fin_nusselt_and_jacobian(Re_d, Pr, L_D, S_D, X_D, True).result[0]

    def f_func(Re_d):
        return cb._core.pin_fin_friction_and_jacobian(Re_d, True).result[0]

    for Re_d in [1000.0, 10000.0, 50000.0]:
        Nu_ana, dNu_ana = cb._core.pin_fin_nusselt_and_jacobian(
            Re_d, Pr, L_D, S_D, X_D, True
        ).result
        dNu_num = central_difference(nu_func, Re_d, max(1e-4, Re_d * 1e-6))
        np.testing.assert_allclose(dNu_ana, dNu_num, rtol=1e-5)

        f_ana, df_ana = cb._core.pin_fin_friction_and_jacobian(Re_d, True).result
        df_num = central_difference(f_func, Re_d, max(1e-4, Re_d * 1e-6))
        np.testing.assert_allclose(df_ana, df_num, rtol=1e-5)


def test_dimple_jacobians():
    d_Dh = 0.2
    h_d = 0.1
    S_d = 2.0  # Must be in [1.5, 3.0]

    def nu_func(Re_Dh):
        return cb._core.dimple_nusselt_enhancement_and_jacobian(Re_Dh, d_Dh, h_d, S_d).result[0]

    def f_func(Re_Dh):
        return cb._core.dimple_friction_multiplier_and_jacobian(Re_Dh, d_Dh, h_d).result[0]

    for Re_Dh in [10000.0, 50000.0]:
        Nu_ana, dNu_ana = cb._core.dimple_nusselt_enhancement_and_jacobian(
            Re_Dh, d_Dh, h_d, S_d
        ).result
        dNu_num = central_difference(nu_func, Re_Dh, max(1e-4, Re_Dh * 1e-6))
        np.testing.assert_allclose(dNu_ana, dNu_num, rtol=1e-5)

        f_ana, df_ana = cb._core.dimple_friction_multiplier_and_jacobian(Re_Dh, d_Dh, h_d).result
        df_num = central_difference(f_func, Re_Dh, max(1e-4, Re_Dh * 1e-6))
        np.testing.assert_allclose(df_ana, df_num, rtol=1e-5)


def test_rib_enhancement_jacobian():
    e_D = 0.05
    pt_h = 10.0
    alpha = 45.0

    def f_func(Re):
        return cb._core.rib_enhancement_factor_high_re_and_jacobian(e_D, pt_h, alpha, Re).result[0]

    for Re in [10000.0, 50000.0]:
        f_ana, df_ana = cb._core.rib_enhancement_factor_high_re_and_jacobian(
            e_D, pt_h, alpha, Re
        ).result
        df_num = central_difference(f_func, Re, max(1e-4, Re * 1e-6))
        np.testing.assert_allclose(df_ana, df_num, rtol=1e-5)


def test_impingement_jacobian():
    Pr = 0.71
    z_D = 3.0

    def f_func(Re_jet):
        return cb._core.impingement_nusselt_and_jacobian(Re_jet, Pr, z_D).result[0]

    for Re_jet in [10000.0, 50000.0]:
        f_ana, df_ana = cb._core.impingement_nusselt_and_jacobian(Re_jet, Pr, z_D).result
        df_num = central_difference(f_func, Re_jet, max(1e-4, Re_jet * 1e-6))
        np.testing.assert_allclose(df_ana, df_num, rtol=1e-5)


def test_film_and_effusion_jacobians():
    x_D = 10.0
    DR = 1.5  # Must be in [1.2, 2.0] for film cooling
    alpha = 30.0
    porosity = 0.05
    s_D = 5.0  # Must be in [4.0, 8.0]

    def film_func(M):
        return cb._core.film_cooling_effectiveness_and_jacobian(x_D, M, DR, alpha).result[0]

    def eff_func(M):
        return cb._core.effusion_effectiveness_and_jacobian(
            x_D, M, DR, porosity, s_D, alpha
        ).result[0]

    for M in [0.5, 1.0, 2.0]:
        f_ana, df_ana = cb._core.film_cooling_effectiveness_and_jacobian(x_D, M, DR, alpha).result
        df_num = central_difference(film_func, M, 1e-5)
        np.testing.assert_allclose(
            df_ana, df_num, rtol=1e-4
        )  # Slightly looser due to M-dependent finite differences

    for M in [1.01, 2.0, 3.5]:
        f_ana, df_ana = cb._core.effusion_effectiveness_and_jacobian(
            x_D, M, DR, porosity, s_D, alpha
        ).result
        df_num = central_difference(
            eff_func, M, 1e-4
        )  # Slightly looser due to M-dependent finite differences
        np.testing.assert_allclose(df_ana, df_num, rtol=1e-3)


def test_stagnation_jacobians():
    T = 300.0
    P = 101325.0
    X = cb.standard_dry_air_composition()

    def mach_func(v):
        return cb._core.mach_number_and_jacobian_v(v, T, X)[0]

    def t_aw_func(v):
        return cb._core.T_adiabatic_wall_and_jacobian_v(T, v, T, P, X)[0]

    for v in [10.0, 100.0, 300.0]:
        ana, dana = cb._core.mach_number_and_jacobian_v(v, T, X)
        dnum = central_difference(mach_func, v, 1e-5)
        np.testing.assert_allclose(dana, dnum, rtol=1e-5)

        ana, dana = cb._core.T_adiabatic_wall_and_jacobian_v(T, v, T, P, X)
        dnum = central_difference(t_aw_func, v, 1e-4)
        np.testing.assert_allclose(dana, dnum, rtol=1e-5)

    def t0_func(M):
        return cb._core.T0_from_static_and_jacobian_M(T, M, X)[0]

    def p0_func(M):
        return cb._core.P0_from_static_and_jacobian_M(P, T, M, X)[0]

    for M in [0.1, 0.5, 1.5]:
        ana, dana = cb._core.T0_from_static_and_jacobian_M(T, M, X)
        dnum = central_difference(t0_func, M, 1e-5)
        np.testing.assert_allclose(dana, dnum, rtol=1e-5)

        ana, dana = cb._core.P0_from_static_and_jacobian_M(P, T, M, X)
        dnum = central_difference(p0_func, M, 1e-5)
        np.testing.assert_allclose(dana, dnum, rtol=1e-4)
