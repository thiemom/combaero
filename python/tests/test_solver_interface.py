from __future__ import annotations

import pytest
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
        return cb._core.nusselt_and_jacobian_dittus_boelter(Re, Pr, heating)[0]

    for Re in [10000.0, 50000.0, 1e6]:
        # Heating
        Nu_analytic, jac_analytic = cb._core.nusselt_and_jacobian_dittus_boelter(Re, Pr_heating, True)
        jac_numeric = central_difference(nu_func, Re, 1e-4, Pr_heating, True) # slightly larger step for large Re
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)

        # Cooling
        Nu_analytic, jac_analytic = cb._core.nusselt_and_jacobian_dittus_boelter(Re, Pr_cooling, False)
        jac_numeric = central_difference(nu_func, Re, 1e-4, Pr_cooling, False)
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)

def test_friction_haaland_jacobian():
    e_D = 1e-4

    def f_func(Re):
        return cb._core.friction_and_jacobian_haaland(Re, e_D)[0]

    for Re in [3000.0, 10000.0, 1e6]:
        f_analytic, jac_analytic = cb._core.friction_and_jacobian_haaland(Re, e_D)

        # Friction factor gradients are very small, adjust h scalar based on Re magnitude
        h = max(1e-4, Re * 1e-6)
        jac_numeric = central_difference(f_func, Re, h)

        assert f_analytic > 0.0
        # Friction factor decreases with Reynolds number
        assert jac_analytic < 0.0
        np.testing.assert_allclose(jac_analytic, jac_numeric, rtol=1e-5)
