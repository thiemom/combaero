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
        # C++ Analytical Tuple (mdot, d_mdot_ddP, d_mdot_drho)
        mdot_analytic, jac_analytic, _ = cb._core.orifice_mdot_and_jacobian(dP, rho, Cd, area)

        # Python Numerical Gradient
        h = max(1e-5, dP * 1e-6)
        jac_numeric = central_difference(mdot_func, dP, h)

        assert mdot_analytic > 0.0
        assert jac_analytic > 0.0
        # The analytical derivative should perfectly match the central difference
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-6)


def test_orifice_bidirectional_flow():
    """Test that orifice handles bidirectional flow (positive and negative dP)."""
    rho = 1.225
    Cd = 0.6
    area = 0.05

    # Test positive, zero, and negative dP
    dP_test_points = [-50000.0, -100.0, -10.0, 0.0, 10.0, 100.0, 50000.0]

    def mdot_func(dP):
        return cb._core.orifice_mdot_and_jacobian(dP, rho, Cd, area)[0]

    for dP in dP_test_points:
        mdot, d_mdot_ddP, _ = cb._core.orifice_mdot_and_jacobian(dP, rho, Cd, area)

        # Sign consistency: mdot should have same sign as dP
        if dP > 0:
            assert mdot > 0.0, f"Positive dP={dP} should give positive mdot, got {mdot}"
        elif dP < 0:
            assert mdot < 0.0, f"Negative dP={dP} should give negative mdot, got {mdot}"
        else:
            assert abs(mdot) < 1e-10, f"Zero dP should give near-zero mdot, got {mdot}"

        # Jacobian should always be positive (mdot increases with dP)
        assert d_mdot_ddP > 0.0, f"Jacobian should be positive, got {d_mdot_ddP} at dP={dP}"

        # Verify Jacobian accuracy with numerical derivative
        if abs(dP) > 1e-3:  # Skip very small dP where numerical derivative is less accurate
            h = max(1e-5, abs(dP) * 1e-6)
            jac_numeric = central_difference(mdot_func, dP, h)
            np.testing.assert_allclose(
                d_mdot_ddP, jac_numeric, rtol=1e-5, err_msg=f"Jacobian mismatch at dP={dP}"
            )


def test_orifice_smooth_at_zero():
    """Test that orifice flow is smooth through dP = 0 (regularization working)."""
    rho = 1.225
    Cd = 0.6
    area = 0.05

    # Test points very close to zero
    dP_near_zero = [-1.0, -0.1, -0.01, 0.0, 0.01, 0.1, 1.0]

    mdot_values = []
    jac_values = []

    for dP in dP_near_zero:
        mdot, d_mdot_ddP, _ = cb._core.orifice_mdot_and_jacobian(dP, rho, Cd, area)
        mdot_values.append(mdot)
        jac_values.append(d_mdot_ddP)

    # Check mdot is continuous (no jumps)
    mdot_diffs = np.diff(mdot_values)
    max_jump = np.max(np.abs(mdot_diffs))
    assert max_jump < 0.1, f"mdot should be smooth near zero, max jump = {max_jump}"

    # Check Jacobian is continuous and positive everywhere
    for dP, jac in zip(dP_near_zero, jac_values, strict=True):
        assert jac > 0.0, f"Jacobian should be positive at dP={dP}, got {jac}"

    jac_diffs = np.diff(jac_values)
    max_jac_jump = np.max(np.abs(jac_diffs))
    # Jacobian should be relatively smooth (within order of magnitude)
    assert max_jac_jump < 10.0 * np.mean(jac_values), (
        f"Jacobian should be smooth near zero, max jump = {max_jac_jump}"
    )


def test_orifice_symmetry():
    """Test that orifice flow is antisymmetric: mdot(-dP) = -mdot(dP)."""
    rho = 1.225
    Cd = 0.6
    area = 0.05

    dP_test_points = [10.0, 100.0, 1000.0, 50000.0]

    for dP in dP_test_points:
        mdot_pos, jac_pos, drho_pos = cb._core.orifice_mdot_and_jacobian(dP, rho, Cd, area)
        mdot_neg, jac_neg, drho_neg = cb._core.orifice_mdot_and_jacobian(-dP, rho, Cd, area)

        # Antisymmetry: mdot(-dP) = -mdot(dP)
        np.testing.assert_allclose(
            mdot_neg, -mdot_pos, rtol=1e-10, err_msg=f"mdot should be antisymmetric at dP={dP}"
        )

        # Jacobian should be symmetric: d_mdot_ddP is same for +dP and -dP
        np.testing.assert_allclose(
            jac_neg, jac_pos, rtol=1e-10, err_msg=f"Jacobian should be symmetric at dP={dP}"
        )

        # Density derivative should be antisymmetric (like mdot, since d_mdot_drho ∝ mdot)
        np.testing.assert_allclose(
            drho_neg,
            -drho_pos,
            rtol=1e-10,
            err_msg=f"d_mdot_drho should be antisymmetric at dP={dP}",
        )


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


def test_lossless_pressure_jacobian():
    def residual_func(P_out):
        return cb._core.lossless_pressure_and_jacobian(150000.0, P_out)[0]

    for P_out in [100000.0, 150000.0, 200000.0]:
        residual, jac_analytic = cb._core.lossless_pressure_and_jacobian(150000.0, P_out)
        jac_numeric = central_difference(residual_func, P_out, 1e-3)

        assert residual == 150000.0 - P_out
        np.testing.assert_allclose(jac_analytic, -1.0, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-6, atol=1e-9)


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


def test_friction_serghides_jacobian():
    e_D = 1e-4

    def f_func(Re):
        return cb._core.friction_and_jacobian_serghides(Re, e_D).result[0]

    for Re in [3000.0, 10000.0, 1e6]:
        f_analytic, jac_analytic = cb._core.friction_and_jacobian_serghides(Re, e_D).result
        h = max(1e-4, Re * 1e-6)
        jac_numeric = central_difference(f_func, Re, h)

        assert f_analytic > 0.0
        assert jac_analytic < 0.0
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)


def test_friction_colebrook_jacobian():
    e_D = 1e-4

    def f_func(Re):
        return cb._core.friction_and_jacobian_colebrook(Re, e_D).result[0]

    for Re in [3000.0, 10000.0, 1e6]:
        f_analytic, jac_analytic = cb._core.friction_and_jacobian_colebrook(Re, e_D).result
        h = max(1e-4, Re * 1e-6)
        jac_numeric = central_difference(f_func, Re, h)

        assert f_analytic > 0.0
        assert jac_analytic < 0.0
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=2e-5)


def test_friction_petukhov_jacobian():
    def f_func(Re):
        return cb._core.friction_and_jacobian_petukhov(Re).result[0]

    for Re in [5000.0, 10000.0, 1e6]:
        result = cb._core.friction_and_jacobian_petukhov(Re)
        f_analytic, jac_analytic = result.result
        h = max(1e-4, Re * 1e-6)
        jac_numeric = central_difference(f_func, Re, h)

        assert result.status.name in ("VALID", "EXTRAPOLATED")
        assert f_analytic > 0.0
        assert jac_analytic < 0.0
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)

    boundary = cb._core.friction_and_jacobian_petukhov(3000.0)
    assert boundary.status.name == "VALID"


def test_friction_dispatcher_matches_direct_models():
    Re = 50000.0
    e_D = 1e-4

    direct = {
        "haaland": cb._core.friction_and_jacobian_haaland(Re, e_D).result,
        "serghides": cb._core.friction_and_jacobian_serghides(Re, e_D).result,
        "colebrook": cb._core.friction_and_jacobian_colebrook(Re, e_D).result,
        "petukhov": cb._core.friction_and_jacobian_petukhov(Re).result,
    }

    for tag, expected in direct.items():
        got = cb._core.friction_and_jacobian(tag, Re, e_D).result
        np.testing.assert_allclose(got[0], expected[0], rtol=1e-12, atol=0.0)
        np.testing.assert_allclose(got[1], expected[1], rtol=1e-12, atol=0.0)

    invalid = cb._core.friction_and_jacobian("does-not-exist", Re, e_D)
    assert invalid.status.name == "INVALID"


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
    X = cb.species.dry_air()

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
        np.testing.assert_allclose(dana, dnum, rtol=2e-5)

        ana, dana = cb._core.P0_from_static_and_jacobian_M(P, T, M, X)
        dnum = central_difference(p0_func, M, 1e-5)
        np.testing.assert_allclose(dana, dnum, rtol=1e-4)


# =============================================================================
# Compressible Flow Elements Tests
# =============================================================================


def test_orifice_compressible_subsonic():
    """Test compressible orifice in subsonic flow regime."""
    T0, P0, P_back = 300.0, 200000.0, 150000.0
    X = cb.species.dry_air()
    Cd, area, beta = 0.65, 1e-4, 0.5

    mdot, d_P0, d_Pb, d_T0 = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0, P_back, X, Cd, area, beta
    )

    # Verify mdot is positive
    assert mdot > 0, "Mass flow should be positive"

    # Verify Jacobians have correct signs
    assert d_P0 > 0, "d_mdot_dP0 should be positive (higher P0 → more flow)"
    assert d_Pb < 0, "d_mdot_dP_back should be negative (higher P_back → less flow)"
    assert d_T0 < 0, "d_mdot_dT0 should be negative (higher T0 → lower density)"

    # Verify against direct nozzle_flow
    A_eff = Cd * area / np.sqrt(1 - beta**4)
    sol = cb.nozzle_flow(T0, P0, P_back, A_eff, X)
    np.testing.assert_allclose(mdot, sol.mdot, rtol=1e-10)


def test_orifice_compressible_choked():
    """Test compressible orifice in choked flow regime."""
    T0, P0 = 300.0, 200000.0
    X = cb.species.dry_air()
    Cd, area, beta = 0.65, 1e-4, 0.5

    # Get critical pressure ratio
    PR_crit = cb.critical_pressure_ratio(T0, P0, X)
    P_back_choked = 0.5 * PR_crit * P0  # Well below critical

    mdot, d_P0, d_Pb, d_T0 = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0, P_back_choked, X, Cd, area, beta
    )

    # In choked flow, d_mdot_dP_back should be near zero (smoothed)
    assert abs(d_Pb) < 1e-10, "d_mdot_dP_back should be ~0 when well choked"

    # But d_mdot_dP0 should still be positive
    assert d_P0 > 0, "d_mdot_dP0 should still be positive when choked"


def test_orifice_compressible_smooth_transition():
    """Verify d_mdot_dP_back is smooth through choked transition."""
    T0, P0 = 300.0, 200000.0
    X = cb.species.dry_air()
    Cd, area, beta = 0.65, 1e-4, 0.5

    # Sweep through critical pressure ratio
    PR_crit = cb.critical_pressure_ratio(T0, P0, X)
    PR_values = np.linspace(PR_crit - 0.05, PR_crit + 0.05, 20)

    jacobians = []
    for PR in PR_values:
        P_back = PR * P0
        _, _, d_Pb, _ = cb._core.orifice_compressible_mdot_and_jacobian(
            T0, P0, P_back, X, Cd, area, beta
        )
        jacobians.append(d_Pb)

    # Check smoothness: no abrupt jumps
    jac_array = np.array(jacobians)

    # Jacobian should be continuous (no large jumps)
    # Check that it transitions smoothly from ~0 (choked) to negative (unchoked)
    assert jac_array[0] == 0.0 or abs(jac_array[0]) < 1e-10, "Should be ~0 when well choked"
    assert jac_array[-1] < -1e-9, "Should be negative when well unchoked"

    # Check no individual jump is > 50% of the range
    jac_diffs = np.diff(jac_array)
    jac_range = abs(jac_array[-1] - jac_array[0])
    max_jump = np.max(np.abs(jac_diffs))
    assert max_jump < 0.5 * jac_range, "No single jump should be > 50% of total range"


def test_orifice_compressible_reverse_flow():
    """Test compressible orifice with reverse flow (P_back > P0)."""
    T0, P0, P_back = 300.0, 150000.0, 200000.0  # Reversed!
    X = cb.species.dry_air()
    Cd, area, beta = 0.65, 1e-4, 0.5

    mdot, d_P0, d_Pb, d_T0 = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0, P_back, X, Cd, area, beta
    )

    # Mass flow should be negative for reverse flow
    assert mdot < 0, "Mass flow should be negative for reverse flow"

    # Increasing P0 reduces the pressure difference -> less reverse flow (less negative)
    assert d_P0 > 0, "d_mdot_dP0 should be positive in reverse flow"
    # Increasing P_back increases the pressure difference -> more reverse flow (more negative)
    assert d_Pb < 0, "d_mdot_dP_back should be negative in reverse flow"
    # Increasing T0 reduces density -> less reverse flow (less negative)
    assert d_T0 > 0, "d_mdot_dT0 should be positive in reverse flow"

    # Validate all three Jacobians against central FD
    eps_P = 1.0
    mp, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0 + eps_P, P_back, X, Cd, area, beta
    )
    mm, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0 - eps_P, P_back, X, Cd, area, beta
    )
    np.testing.assert_allclose(d_P0, (mp - mm) / (2 * eps_P), rtol=1e-4)

    mp2, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0, P_back + eps_P, X, Cd, area, beta
    )
    mm2, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0, P_back - eps_P, X, Cd, area, beta
    )
    np.testing.assert_allclose(d_Pb, (mp2 - mm2) / (2 * eps_P), rtol=1e-4)

    eps_T = 1e-3
    mp3, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        T0 + eps_T, P0, P_back, X, Cd, area, beta
    )
    mm3, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        T0 - eps_T, P0, P_back, X, Cd, area, beta
    )
    np.testing.assert_allclose(d_T0, (mp3 - mm3) / (2 * eps_T), rtol=1e-4)


def test_orifice_compressible_jacobian_accuracy():
    """Verify Jacobian accuracy via numerical differentiation."""
    T0, P0, P_back = 300.0, 200000.0, 150000.0
    X = cb.species.dry_air()
    Cd, area, beta = 0.65, 1e-4, 0.5

    mdot, d_P0, d_Pb, d_T0 = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0, P_back, X, Cd, area, beta
    )

    # Numerical Jacobian w.r.t. P0
    eps_P0 = max(1e-6, P0 * 1e-6)
    mdot_plus, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0 + eps_P0, P_back, X, Cd, area, beta
    )
    d_P0_num = (mdot_plus - mdot) / eps_P0
    np.testing.assert_allclose(d_P0, d_P0_num, rtol=1e-5)

    # Numerical Jacobian w.r.t. T0
    eps_T0 = max(1e-6, T0 * 1e-6)
    mdot_plus, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        T0 + eps_T0, P0, P_back, X, Cd, area, beta
    )
    d_T0_num = (mdot_plus - mdot) / eps_T0
    np.testing.assert_allclose(d_T0, d_T0_num, rtol=1e-5)


def test_orifice_compressible_matches_nozzle_flow():
    """Verify compressible orifice matches nozzle_flow exactly."""
    T0, P0, P_back = 300.0, 200000.0, 150000.0
    X = cb.species.dry_air()
    Cd, area, beta = 0.65, 1e-4, 0.5

    mdot, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        T0, P0, P_back, X, Cd, area, beta
    )

    # Direct nozzle_flow call
    A_eff = Cd * area / np.sqrt(1 - beta**4)
    sol = cb.nozzle_flow(T0, P0, P_back, A_eff, X)

    # Should match exactly (no smoothing on mdot itself)
    np.testing.assert_allclose(mdot, sol.mdot, rtol=1e-10)


def test_orifice_compressible_smoothing_accuracy_far():
    """Verify smoothing has negligible effect far from critical PR."""
    T0, P0 = 300.0, 200000.0
    X = cb.species.dry_air()
    Cd, area, beta = 0.65, 1e-4, 0.5

    PR_crit = cb.critical_pressure_ratio(T0, P0, X)

    # Test well away from critical (> 3x smoothing width)
    test_cases = [
        ("well unchoked", PR_crit + 0.05),
        ("well choked", PR_crit - 0.05),
    ]

    A_eff = Cd * area / np.sqrt(1 - beta**4)

    for label, PR in test_cases:
        P_back = PR * P0

        # Smoothed implementation
        mdot_smooth, _, d_Pb_smooth, _ = cb._core.orifice_compressible_mdot_and_jacobian(
            T0, P0, P_back, X, Cd, area, beta
        )

        # Exact from nozzle_flow
        sol_exact = cb.nozzle_flow(T0, P0, P_back, A_eff, X)

        # mdot should always match exactly
        np.testing.assert_allclose(
            mdot_smooth, sol_exact.mdot, rtol=1e-10, err_msg=f"mdot error {label}"
        )

        # For well unchoked case, Jacobian should be accurate
        if PR_crit + 0.02 < PR:
            # Compute exact Jacobian numerically
            eps_Pb = max(1e-6, abs(P_back) * 1e-6)
            sol_plus = cb.nozzle_flow(T0, P0, P_back + eps_Pb, A_eff, X)
            sol_minus = cb.nozzle_flow(T0, P0, P_back - eps_Pb, A_eff, X)
            d_Pb_exact = (sol_plus.mdot - sol_minus.mdot) / (2 * eps_Pb)

            # Should be within 0.1% far from transition
            if abs(d_Pb_exact) > 1e-10:
                np.testing.assert_allclose(
                    d_Pb_smooth, d_Pb_exact, rtol=1e-3, err_msg=f"Jacobian error {label}"
                )


def test_channel_compressible_low_mach():
    """Test compressible channel at low Mach number."""
    T_in, P_in, u_in = 400.0, 200000.0, 10.0  # Low velocity
    X = cb.species.dry_air()
    L, D, roughness = 2.0, 0.05, 1e-4

    dP, d_Pin, d_Tin, d_u = cb._core.channel_compressible_mdot_and_jacobian(
        T_in, P_in, u_in, X, L, D, roughness, "haaland"
    )

    # Verify dP is positive
    assert dP > 0, "Pressure drop should be positive"

    # Verify Jacobians have reasonable magnitudes
    assert d_Pin > 0, "d_dP_dP_in should be positive"
    assert abs(d_Tin) > 0, "d_dP_dT_in should be non-zero"
    assert d_u > 0, "d_dP_du_in should be positive (higher velocity → more friction)"


def test_channel_compressible_high_mach():
    """Test compressible channel at higher Mach number."""
    T_in, P_in, u_in = 400.0, 200000.0, 150.0  # Higher velocity
    X = cb.species.dry_air()
    L, D, roughness = 2.0, 0.05, 1e-4

    dP, d_Pin, d_Tin, d_u = cb._core.channel_compressible_mdot_and_jacobian(
        T_in, P_in, u_in, X, L, D, roughness, "haaland"
    )

    # Should still compute successfully
    assert dP > 0, "Pressure drop should be positive"
    assert np.isfinite(dP), "Pressure drop should be finite"


def test_channel_compressible_matches_fanno():
    """Verify compressible channel matches fanno_channel_rough exactly."""
    T_in, P_in, u_in = 400.0, 200000.0, 50.0
    X = cb.species.dry_air()
    L, D, roughness = 2.0, 0.05, 1e-4

    dP, _, _, _ = cb._core.channel_compressible_mdot_and_jacobian(
        T_in, P_in, u_in, X, L, D, roughness, "haaland"
    )

    # Direct fanno_channel_rough call
    sol = cb.fanno_channel_rough(T_in, P_in, u_in, L, D, roughness, X, "haaland")
    dP_direct = P_in - sol.outlet.P

    # Should match exactly
    np.testing.assert_allclose(dP, dP_direct, rtol=1e-10)


def test_channel_compressible_jacobian_accuracy():
    """Verify channel Jacobian accuracy via numerical differentiation."""
    T_in, P_in, u_in = 400.0, 200000.0, 50.0
    X = cb.species.dry_air()
    L, D, roughness = 2.0, 0.05, 1e-4

    dP, d_Pin, d_Tin, d_u = cb._core.channel_compressible_mdot_and_jacobian(
        T_in, P_in, u_in, X, L, D, roughness, "haaland"
    )

    # Numerical Jacobian w.r.t. P_in
    eps_P = max(1e-6, P_in * 1e-6)
    dP_plus, _, _, _ = cb._core.channel_compressible_mdot_and_jacobian(
        T_in, P_in + eps_P, u_in, X, L, D, roughness, "haaland"
    )
    d_Pin_num = (dP_plus - dP) / eps_P
    np.testing.assert_allclose(d_Pin, d_Pin_num, rtol=1e-5)

    # Numerical Jacobian w.r.t. u_in
    eps_u = max(1e-6, u_in * 1e-6)
    dP_plus, _, _, _ = cb._core.channel_compressible_mdot_and_jacobian(
        T_in, P_in, u_in + eps_u, X, L, D, roughness, "haaland"
    )
    d_u_num = (dP_plus - dP) / eps_u
    np.testing.assert_allclose(d_u, d_u_num, rtol=1e-5)


def test_channel_compressible_reverse_flow():
    """Test compressible channel with reverse flow (negative velocity)."""
    T_in, P_in, u_in = 400.0, 200000.0, -50.0  # Negative velocity
    X = cb.species.dry_air()
    L, D, roughness = 2.0, 0.05, 1e-4

    dP, d_Pin, d_Tin, d_u = cb._core.channel_compressible_mdot_and_jacobian(
        T_in, P_in, u_in, X, L, D, roughness, "haaland"
    )

    # Pressure drop should be negative for reverse flow
    assert dP < 0, "Pressure drop should be negative for reverse flow"

    # Jacobians should still be computed
    assert np.isfinite(d_Pin), "d_dP_dP_in should be finite"
    assert np.isfinite(d_u), "d_dP_du_in should be finite"


def test_compressible_network_scenario():
    """Test compressible elements in a realistic network scenario.

    Network: High Pressure -> Channel -> Junction -> Orifice -> Low Pressure
    This mimics the pattern from test_network_scenarios.py but uses compressible models.
    """
    # Use standard air composition (not hardcoded)
    X = cb.species.dry_air()
    Y = cb.mole_to_mass(X)

    # Operating conditions
    T = 300.0  # K
    P_high = 300000.0  # 3 bar
    P_low = 101325.0  # 1 atm

    # Geometry
    L_channel = 2.0  # m
    D_channel = 0.05  # m
    roughness = 1e-4  # m
    A_orifice = 1e-4  # m²
    Cd = 0.65
    beta = 0.5

    # Test channel element
    rho = cb.density(T, P_high, X)
    area_channel = 0.25 * np.pi * D_channel**2
    u_channel = 50.0  # m/s estimate

    dP_channel, d_dP_dP, d_dP_dT, d_dP_du = cb._core.channel_compressible_mdot_and_jacobian(
        T, P_high, u_channel, X, L_channel, D_channel, roughness, "haaland"
    )

    # Verify channel results
    assert dP_channel > 0, "Channel pressure drop should be positive"
    assert d_dP_dP > 0, "d_dP/dP_in should be positive"
    assert d_dP_du > 0, "d_dP/du should be positive (more velocity = more friction)"

    # Compute mass flow through channel
    mdot_channel = rho * area_channel * u_channel

    # Test orifice element at junction
    P_junction = P_high - dP_channel

    mdot_orifice, d_mdot_dP0, d_mdot_dPb, d_mdot_dT0 = (
        cb._core.orifice_compressible_mdot_and_jacobian(
            T, P_junction, P_low, X, Cd, A_orifice, beta
        )
    )

    # Verify orifice results
    assert mdot_orifice > 0, "Orifice mass flow should be positive"
    assert d_mdot_dP0 > 0, "d_mdot/dP0 should be positive"

    # Check pressure ratio to see if we're in compressible regime
    PR = P_low / P_high
    assert PR < 0.8, f"Test should be in compressible regime (PR={PR:.3f} < 0.8)"

    # Check if choked - if so, d_mdot_dP_back should be ~0 (smoothed)
    PR_crit = cb.critical_pressure_ratio(T, P_junction, X)
    if PR_crit > PR:
        # Choked flow - Jacobian should be near zero
        assert abs(d_mdot_dPb) < 1e-8, (
            f"d_mdot/dP_back should be ~0 when choked (PR={PR:.3f} < PR_crit={PR_crit:.3f})"
        )
    else:
        # Unchoked flow - Jacobian should be negative
        assert d_mdot_dPb < 0, "d_mdot/dP_back should be negative when unchoked"

    # Verify mass flow consistency (channel and orifice should be similar order of magnitude)
    # They won't match exactly since we're not solving the full network, but should be close
    flow_ratio = mdot_orifice / mdot_channel
    assert 0.1 < flow_ratio < 10.0, (
        f"Flow rates should be similar order of magnitude (ratio={flow_ratio:.2f})"
    )

    # Test with mass fractions (not mole fractions) for residuals function
    res_orifice = cb._core.orifice_compressible_residuals_and_jacobian(
        mdot_orifice, P_junction, T, Y, P_low, Cd, A_orifice, beta
    )

    # Verify residuals structure
    assert hasattr(res_orifice, "m_dot_calc"), "Should have m_dot_calc attribute"
    assert hasattr(res_orifice, "d_mdot_dP_total_up"), "Should have Jacobian attributes"
    np.testing.assert_allclose(
        res_orifice.m_dot_calc,
        mdot_orifice,
        rtol=1e-10,
        err_msg="Residuals should match direct call",
    )

    # Test channel residuals function
    res_channel = cb._core.channel_compressible_residuals_and_jacobian(
        mdot_channel, P_high, T, Y, P_low, L_channel, D_channel, roughness, "haaland"
    )

    # Verify channel residuals structure
    assert hasattr(res_channel, "dP_calc"), "Should have dP_calc attribute"
    assert hasattr(res_channel, "d_dP_d_mdot"), "Should have Jacobian attributes"
    assert res_channel.dP_calc > 0, "Channel pressure drop should be positive"
