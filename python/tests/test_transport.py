"""Tests for transport properties: reynolds, peclet, thermal_diffusivity."""

import numpy as np

import combaero as cb


def test_reynolds_air() -> None:
    """Reynolds number for air flow."""
    X = cb.species.dry_air()
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
    X = cb.species.dry_air()
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
    X = cb.species.dry_air()
    T = 300.0  # K
    P = 101325.0  # Pa

    alpha = cb.thermal_diffusivity(T, P, X)

    # alpha = k / (rho * cp)
    # For air at STP: alpha ~ 2e-5 m^2/s
    assert 1e-5 < alpha < 5e-5


def test_reynolds_peclet_consistency() -> None:
    """Pe = V*L/alpha and Re = rho*V*L/mu should satisfy Pe = Re*Pr."""
    X = cb.species.dry_air()
    T = 400.0  # K
    P = 200000.0  # Pa
    V = 50.0  # m/s
    L = 0.5  # m

    Re = cb.reynolds(T, P, X, V, L)
    Pr = cb.prandtl(T, P, X)
    Pe = cb.peclet(T, P, X, V, L)

    # Pe = Re * Pr
    assert np.isclose(Pe, Re * Pr, rtol=1e-6)


def _mix(fracs: dict[str, float]) -> list[float]:
    x = np.zeros(cb.num_species())
    for name, v in fracs.items():
        x[cb.species_index_from_name(name)] = v
    return list(x / x.sum())


def test_1000K_polynomial_join_is_value_continuous_and_slopes_nearly_so() -> None:
    """The species polynomials switch range at 1000 K (docs/API_CPP.md,
    Derivatives). Values are continuous; slopes are one-sided but small for
    real mixtures (measured 2026-10-11: dk/dT -0.01% air, -0.10% lean CH4
    products). Bounds pin that, so a data change that widens the kink fails."""
    P, e = 101325.0, 1e-6
    for X, slope_bound in (
        (cb.species.dry_air(), 1e-3),
        (_mix({"N2": 0.73, "O2": 0.11, "H2O": 0.1, "CO2": 0.05, "AR": 0.009}), 1e-2),
    ):
        lo, hi = cb.transport_and_dT(1000 - e, P, X), cb.transport_and_dT(1000 + e, P, X)
        assert abs(hi.k - lo.k) / lo.k < 1e-7
        assert abs(hi.mu - lo.mu) / lo.mu < 1e-7
        assert abs(cb.cp(1000 + e, X) - cb.cp(1000 - e, X)) / cb.cp(1000, X) < 1e-7
        assert abs(hi.dk_dT - lo.dk_dT) / abs(lo.dk_dT) < slope_bound
        dl, dh = cb.dcp_dT(1000 - e, X), cb.dcp_dT(1000 + e, X)
        assert abs(dh - dl) / abs(dl) < slope_bound
