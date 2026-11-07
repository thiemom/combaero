import numpy as np
import pytest

import combaero as cb


def test_combustion_smoothing_impact_and_jacobian():
    """
    Verify smoothing impact on temperature (< 1K) and Jacobian continuity.
    Tests across equivalence ratio range [0, 2.0] including kinks at 0 and 1.
    """
    n_sp = cb.species.num_species
    X_fuel = np.zeros(n_sp)
    X_fuel[cb.species.species_index("CH4")] = 1.0

    X_ox = np.zeros(n_sp)
    X_ox[cb.species.species_index("O2")] = 0.21
    X_ox[cb.species.species_index("N2")] = 0.79

    T_reactants = 300.0
    P = 101325.0

    # Increase resolution around the kink at Phi=1
    phis = np.linspace(1e-9, 2.0, 10001)

    t_no_smooth = []
    t_smooth = []

    for phi in phis:
        X_mix = cb.set_equivalence_ratio_mole(phi, X_fuel, X_ox)
        s_raw = cb.complete_combustion(T_reactants, X_mix, P, smooth_phi0=False, smooth_phi1=False)
        t_no_smooth.append(s_raw.T)

        s_sm = cb.complete_combustion(T_reactants, X_mix, P, smooth_phi0=True, smooth_phi1=True)
        t_smooth.append(s_sm.T)

    t_no_smooth = np.array(t_no_smooth)
    t_smooth = np.array(t_smooth)
    diff = np.abs(t_no_smooth - t_smooth)

    # 1. VERIFY IMPACT: Max temperature difference should be < 1K with k=20000
    max_diff = np.max(diff)
    print(f"\nMax temperature difference with smoothing (k=20000): {max_diff:.4f} K")
    assert max_diff < 1.0, f"Smoothing impact too large: {max_diff:.2f} K > 1.0 K"

    # 2. VERIFY SMOOTHNESS: Max curvature (numerical second derivative)
    # We check both raw and smoothed.
    # dT/dphi jump at Phi=1 in complete combustion is inherent.
    # Smoothing should reduce the magnitude of the discrete jump.
    dT_dphi_raw = np.diff(t_no_smooth) / np.diff(phis)
    dT_dphi_sm = np.diff(t_smooth) / np.diff(phis)

    jump_raw = np.max(np.abs(np.diff(dT_dphi_raw)))
    jump_sm = np.max(np.abs(np.diff(dT_dphi_sm)))

    print(f"Max curvature (jump in dT/dphi) raw: {jump_raw:.4f}")
    print(f"Max curvature (jump in dT/dphi) smoothed: {jump_sm:.4f}")

    # At k=20000, transition is sharp relative to grid, but still finite.
    assert jump_sm < 10000.0, f"Smoothed Jacobian seems to have excessive curvature: {jump_sm:.2f}"
    assert jump_sm < jump_raw, "Smoothing didn't reduce slope discontinuity"


# --- #1037: the zero-fuel smoothing is exact away from zero fuel -------------

_NAMES = [cb.species_name(i) for i in range(cb.num_species())]
_X_AIR = np.array(cb.species.dry_air())
_K0 = 20000.0  # SMOOTHING_K_PHI0 (include/combustion.h), the default k0


def _g(x: np.ndarray) -> np.ndarray:
    """g(x) = ln(1 + e^x) - ln2 sech(x), written without cancellation:
    x/2 + ln cosh(x/2) + ln2 (1 - sech x)."""
    return (
        x / 2 + np.log1p(2 * np.sinh(x / 4) ** 2) + 2 * np.log(2) * np.sinh(x / 2) ** 2 / np.cosh(x)
    )


def _ch4_air(xf: float) -> np.ndarray:
    X = _X_AIR * (1.0 - xf)
    X[_NAMES.index("CH4")] += xf
    return X


@pytest.mark.parametrize("fuel", ["CH4", "H2", "NC7H16", "NH3"])
def test_zero_fuel_smoothing_leaves_no_offset_at_any_phi(fuel: str) -> None:
    """smooth_phi0 on vs off: T_ad agrees to 1e-9 relative wherever k0 a > 50
    (Phi 0.02..2). The softplus it replaced, (ln(1 + e^(k0 a)) - ln2)/k0,
    under-burnt ln2/k0 mol O2 per mol at every phi: CH4 -0.47 K at Phi 0.02,
    -0.34 K at Phi 1 (relative 1.3e-3 at Phi 0.02 on T_ad - T_in)."""
    Xf = np.zeros(len(_NAMES))
    Xf[_NAMES.index(fuel)] = 1.0
    for phi in (0.02, 0.05, 0.1, 0.3, 0.6, 0.9, 1.0, 1.1, 1.5, 2.0):
        X = np.array(cb.set_equivalence_ratio_mole(phi, list(Xf), list(_X_AIR)))
        on = cb.complete_combustion(300.0, X, 101325.0, smooth_phi0=True, smooth_phi1=False)
        off = cb.complete_combustion(300.0, X, 101325.0, smooth_phi0=False, smooth_phi1=False)
        assert pytest.approx(off.T, rel=1e-9, abs=0.0) == on.T, (fuel, phi)
        assert np.allclose(on.X, off.X, rtol=1e-9, atol=1e-15), (fuel, phi)


def test_zero_fuel_smoothing_is_smooth_through_zero_fuel() -> None:
    """The burnt share s = n_req / a follows g(k0 a)/(k0 a) from a = 0 (s = 1/2,
    as before) to 1, monotone in n_req, never above a, with a continuous
    slope: a central difference of n_req(a) across the x = 1e-3 series switch
    and through x = 0..60 matches g'(x)."""
    xs = np.concatenate([np.logspace(-9, 2, 600), np.linspace(9.0e-4, 1.1e-3, 201)])
    xs.sort()
    s = np.empty_like(xs)
    a = np.empty_like(xs)
    for i, x in enumerate(xs):
        X = _ch4_air(x / (2.0 * _K0))
        a[i] = 2.0 * X[_NAMES.index("CH4")]
        s[i] = cb.complete_combustion_to_CO2_H2O_with_fraction(X, True)[1]
    x = _K0 * a
    ref = _g(x) / x
    assert np.abs(s - ref).max() < 1e-12
    assert s.max() <= 1.0
    assert np.all(np.diff(s * a) > 0.0)
    assert s[0] == pytest.approx(0.5, abs=1e-8)
    assert s[x > 40.0].min() == 1.0
    # Zero fuel: unchanged mixture, nothing burnt.
    Xb, f0 = cb.complete_combustion_to_CO2_H2O_with_fraction(_ch4_air(0.0), True)
    assert f0 == 0.0
    assert np.array_equal(np.array(Xb), _ch4_air(0.0))

    def n_req(x: float) -> float:
        X = _ch4_air(x / (2.0 * _K0))
        aa = 2.0 * X[_NAMES.index("CH4")]
        return cb.complete_combustion_to_CO2_H2O_with_fraction(X, True)[1] * aa

    for x0 in (2e-3, 0.999e-3, 1.001e-3, 0.3, 1.0, 5.0, 60.0):
        h = 1e-4 * x0
        fd = (n_req(x0 + h) - n_req(x0 - h)) / (2 * h / _K0)
        t = np.tanh(x0)
        dg = 1 / (1 + np.exp(-x0)) + np.log(2) * t / np.cosh(x0)
        assert fd == pytest.approx(dg, rel=1e-6), x0


@pytest.mark.parametrize("x", [0.0, 1e-4, 1e-3, 1e-2, 0.5, 1.0, 3.0, 10.0, 60.0])
def test_products_jacobian_through_zero_fuel(x: float) -> None:
    """dX_b/dX against differences at k0 a = x near zero fuel (#1037); at
    a = 0 exactly, the one-sided (second-order forward) limit."""
    xf = x / (2.0 * _K0)
    X0 = _ch4_air(xf)
    v = np.zeros(len(_NAMES))
    v[_NAMES.index("CH4")], v[_NAMES.index("N2")] = 1.0, -1.0
    cj = cb.complete_combustion_products_and_jacobian(X0)

    def f(d: float) -> np.ndarray:
        return np.array(cb.complete_combustion_products_and_jacobian(X0 + d * v).X_b)

    if x == 0.0:
        h = 1e-9
        fd = (-3 * f(0.0) + 4 * f(h) - f(2 * h)) / (2 * h)
    else:
        h = min(1e-10, 0.25 * xf)
        fd = (f(h) - f(-h)) / (2 * h)
    an = np.array(cj.dXb_dX) @ v
    assert np.abs(an - fd).max() < 2e-6 * np.abs(fd).max(), x


if __name__ == "__main__":
    test_combustion_smoothing_impact_and_jacobian()
