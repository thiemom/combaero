import pytest

import combaero as cb


def test_orifice_jacobian_accuracy():
    m_dot = 0.5
    P_total_up = 2e5
    P_static_up = 1.9e5
    T_up = 300.0
    Y_up = [1.0]  # Pure species or whatever num_species is

    # Ensure num_species matches
    n_sp = cb.num_species()
    Y_up = [0.0] * n_sp
    Y_up[0] = 1.0  # Assume first species is major

    P_static_down = 1e5
    Cd = 0.8
    area = 0.01

    # C++ Calculation
    res_obj = cb.orifice_residuals_and_jacobian(
        m_dot, P_total_up, P_static_up, T_up, Y_up, P_static_down, Cd, area
    )

    # Numerical derivatives (Central Finite Difference)
    eps = 1e-6

    def get_mdot(p_tot, p_stat, t, y, p_down):
        r = cb.orifice_residuals_and_jacobian(m_dot, p_tot, p_stat, t, y, p_down, Cd, area)
        return r.m_dot_calc

    # d_mdot/dP_total_up
    d_tot_num = (
        get_mdot(P_total_up + eps, P_static_up, T_up, Y_up, P_static_down)
        - get_mdot(P_total_up - eps, P_static_up, T_up, Y_up, P_static_down)
    ) / (2 * eps)

    # d_mdot/dP_static_up
    d_stat_num = (
        get_mdot(P_total_up, P_static_up + eps, T_up, Y_up, P_static_down)
        - get_mdot(P_total_up, P_static_up - eps, T_up, Y_up, P_static_down)
    ) / (2 * eps)

    # d_mdot/dT_up
    d_temp_num = (
        get_mdot(P_total_up, P_static_up, T_up + eps, Y_up, P_static_down)
        - get_mdot(P_total_up, P_static_up, T_up - eps, Y_up, P_static_down)
    ) / (2 * eps)

    # d_mdot/dP_static_down
    d_down_num = (
        get_mdot(P_total_up, P_static_up, T_up, Y_up, P_static_down + eps)
        - get_mdot(P_total_up, P_static_up, T_up, Y_up, P_static_down - eps)
    ) / (2 * eps)

    print("\nOrifice Comparison:")
    print(f"Total Up: Analytical={res_obj.d_mdot_dP_total_up:.6e}, Numerical={d_tot_num:.6e}")
    print(f"Static Up: Analytical={res_obj.d_mdot_dP_static_up:.6e}, Numerical={d_stat_num:.6e}")
    print(f"Temp Up: Analytical={res_obj.d_mdot_dT_up:.6e}, Numerical={d_temp_num:.6e}")
    print(
        f"Static Down: Analytical={res_obj.d_mdot_dP_static_down:.6e}, Numerical={d_down_num:.6e}"
    )

    assert res_obj.d_mdot_dP_total_up == pytest.approx(d_tot_num, rel=1e-4)
    assert res_obj.d_mdot_dP_static_up == pytest.approx(d_stat_num, rel=1e-4)
    assert res_obj.d_mdot_dT_up == pytest.approx(d_temp_num, rel=1e-4)
    assert res_obj.d_mdot_dP_static_down == pytest.approx(d_down_num, rel=1e-4)


def test_pipe_jacobian_accuracy():
    m_dot = 0.5
    P_total_up = 2e5
    P_static_up = 1.9e5
    T_up = 300.0
    n_sp = cb.num_species()
    Y_up = [0.0] * n_sp
    Y_up[0] = 1.0
    P_static_down = 1e5
    L = 1.0
    D = 0.1
    roughness = 1e-5
    friction_model = "haaland"

    # C++ Calculation
    res_obj = cb.pipe_residuals_and_jacobian(
        m_dot, P_total_up, P_static_up, T_up, Y_up, P_static_down, L, D, roughness, friction_model
    )

    eps = 1e-6
    if abs(m_dot) < 0.1:
        eps = 1e-7  # Better precision near zero

    def get_dp(m, p_tot, p_stat, t, y, p_down):
        r = cb.pipe_residuals_and_jacobian(
            m, p_tot, p_stat, t, y, p_down, L, D, roughness, friction_model
        )
        return r.dP_calc

    # d_dP/d_mdot
    d_mdot_num = (
        get_dp(m_dot + eps, P_total_up, P_static_up, T_up, Y_up, P_static_down)
        - get_dp(m_dot - eps, P_total_up, P_static_up, T_up, Y_up, P_static_down)
    ) / (2 * eps)

    # d_dP/dP_static_up
    d_stat_num = (
        get_dp(m_dot, P_total_up, P_static_up + eps, T_up, Y_up, P_static_down)
        - get_dp(m_dot, P_total_up, P_static_up - eps, T_up, Y_up, P_static_down)
    ) / (2 * eps)

    # d_dP/dT_up
    d_temp_num = (
        get_dp(m_dot, P_total_up, P_static_up, T_up + eps, Y_up, P_static_down)
        - get_dp(m_dot, P_total_up, P_static_up, T_up - eps, Y_up, P_static_down)
    ) / (2 * eps)

    print("\nPipe Comparison:")
    print(f"m_dot: Analytical={res_obj.d_dP_d_mdot:.6e}, Numerical={d_mdot_num:.6e}")
    print(f"Static Up: Analytical={res_obj.d_dP_dP_static_up:.6e}, Numerical={d_stat_num:.6e}")
    print(f"Temp Up: Analytical={res_obj.d_dP_dT_up:.6e}, Numerical={d_temp_num:.6e}")

    assert res_obj.d_dP_d_mdot == pytest.approx(d_mdot_num, rel=1e-4)
    assert res_obj.d_dP_dP_static_up == pytest.approx(d_stat_num, rel=1e-4)
    assert res_obj.d_dP_dT_up == pytest.approx(d_temp_num, rel=1e-4)
