"""Tests for transport properties: reynolds, peclet, thermal_diffusivity."""

import numpy as np

import combaero as cb


def test_reynolds_air() -> None:
    """Reynolds number for air flow."""
    X = cb.standard_dry_air_composition()
    T = 300.0  # K
    P = 101325.0  # Pa
    V = 10.0  # m/s
    L = 0.1  # m (characteristic length)

    Re = cb.reynolds(T, P, X, V, L)

    # Re should be positive and reasonable for air
    assert Re > 0
    # For air at STP: rho ~ 1.2 kg/m^3, mu ~ 1.8e-5 Pa*s
    # Re ~ 1.2 * 10 * 0.1 / 1.8e-5 ~ 67000
    assert 50000 < Re < 100000


def test_peclet_air() -> None:
    """Peclet number for air flow."""
    X = cb.standard_dry_air_composition()
    T = 300.0  # K
    P = 101325.0  # Pa
    V = 10.0  # m/s
    L = 0.1  # m

    Pe = cb.peclet(T, P, X, V, L)

    # Pe should be positive
    assert Pe > 0
    # Pe = Re * Pr, and Pr ~ 0.7 for air
    Re = cb.reynolds(T, P, X, V, L)
    Pr = cb.prandtl(T, P, X)
    assert np.isclose(Pe, Re * Pr, rtol=0.01)


def test_thermal_diffusivity_air() -> None:
    """Thermal diffusivity for air."""
    X = cb.standard_dry_air_composition()
    T = 300.0  # K
    P = 101325.0  # Pa

    alpha = cb.thermal_diffusivity(T, P, X)

    # alpha = k / (rho * cp)
    # For air at STP: alpha ~ 2e-5 m^2/s
    assert 1e-5 < alpha < 5e-5


def test_reynolds_peclet_consistency() -> None:
    """Pe = V*L/alpha and Re = rho*V*L/mu should satisfy Pe = Re*Pr."""
    X = cb.standard_dry_air_composition()
    T = 400.0  # K
    P = 200000.0  # Pa
    V = 50.0  # m/s
    L = 0.5  # m

    Re = cb.reynolds(T, P, X, V, L)
    Pr = cb.prandtl(T, P, X)
    Pe = cb.peclet(T, P, X, V, L)

    # Pe = Re * Pr
    assert np.isclose(Pe, Re * Pr, rtol=1e-6)
