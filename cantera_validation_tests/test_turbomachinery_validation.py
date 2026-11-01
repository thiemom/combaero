"""#893: cb.compression / cb.expansion against Cantera, on CombAero's NASA-9
database (combaero_nasa9.yaml).

Isentropic efficiency: Cantera's own isentrope (SP) and enthalpy (HP) state
setters. Polytropic efficiency: the path's ODE integrated on Cantera's cp,
dT/d(ln P) = R T / (eta cp) (compression) or eta R T / cp (expansion), at
rtol 1e-12; independent of CombAero's closed form in S(T). Measured: 8e-10
(compression) and 1.6e-9 (expansion). A stepwise integration through
Cantera's HP setter wandered at 1e-6, its own noise.

Tolerance 1e-7: CombAero's default thermo (the C2 refit, #549) sits about
3e-8 from the NASA-9 YAML on the isentropic paths.
"""

from __future__ import annotations

import os

import cantera as ct
import numpy as np
import pytest
from scipy.integrate import solve_ivp

import combaero as cb

CB_SPECIES = list(cb.species.names)


@pytest.fixture(scope="module")
def gas() -> ct.Solution:
    return ct.Solution(os.path.join(os.path.dirname(__file__), "combaero_nasa9.yaml"))


def _ct_X(X: list[float]) -> dict[str, float]:
    return {n: x for n, x in zip(CB_SPECIES, X, strict=True) if x > 0.0}


def _humid_air() -> list[float]:
    return list(cb.mass_to_mole(list(cb.species.humid_air_mass(288.15, 101325.0, 0.6))))


def _products() -> list[float]:
    """Lean CH4 products: complete combustion of air with 2 % fuel by mass."""
    X_air = list(cb.species.dry_air())
    X_f = list(cb.species.pure_species("CH4"))
    s_air, s_f = cb.Stream(), cb.Stream()
    for st, X, m, T in ((s_air, X_air, 1.0, 700.0), (s_f, X_f, 0.02, 300.0)):
        st.set_T(T)
        st.set_P(1.0e6)
        st.set_X(X)
        st.set_mdot(m)
    mixed = cb.mix([s_air, s_f])
    return list(cb.complete_combustion(mixed.T, list(mixed.X), mixed.P).X)


def _ct_isentropic(gas, X, T1, P1, P2, eta, expansion):
    gas.TPX = T1, P1, _ct_X(X)
    h1, s1 = gas.enthalpy_mass, gas.entropy_mass
    gas.SP = s1, P2
    h2s = gas.enthalpy_mass
    h2 = h1 - eta * (h1 - h2s) if expansion else h1 + (h2s - h1) / eta
    gas.HP = h2, P2
    return gas.T


def _ct_polytropic(gas, X, T1, P1, P2, eta, expansion):
    gas.TPX = T1, P1, _ct_X(X)
    R = ct.gas_constant / gas.mean_molecular_weight
    f = eta if expansion else 1.0 / eta

    def rhs(lnP: float, T: np.ndarray) -> list[float]:
        gas.TP = T[0], np.exp(lnP)
        return [f * R * T[0] / gas.cp_mass]

    sol = solve_ivp(rhs, (np.log(P1), np.log(P2)), [T1], rtol=1e-12, atol=1e-10, method="DOP853")
    return float(sol.y[0, -1])


@pytest.mark.parametrize("PR", [2.0, 10.0, 40.0])
@pytest.mark.parametrize("eta", [1.0, 0.88, 0.7])
def test_compression_of_humid_air(gas, PR: float, eta: float) -> None:
    X, T1, P1 = _humid_air(), 288.15, 101325.0
    r = cb.compression(T1, X, PR, eta_is=eta)
    assert pytest.approx(_ct_isentropic(gas, X, T1, P1, PR * P1, eta, False), rel=1e-7) == r.T2
    r = cb.compression(T1, X, PR, eta_poly=eta)
    assert pytest.approx(_ct_polytropic(gas, X, T1, P1, PR * P1, eta, False), rel=1e-7) == r.T2


@pytest.mark.parametrize("PR", [2.0, 5.0, 20.0])
@pytest.mark.parametrize("eta", [1.0, 0.9, 0.75])
def test_expansion_of_combustion_products(gas, PR: float, eta: float) -> None:
    X, T1, P1 = _products(), 1600.0, 2.0e6
    r = cb.expansion(T1, X, PR, eta_is=eta)
    assert pytest.approx(_ct_isentropic(gas, X, T1, P1, P1 / PR, eta, True), rel=1e-7) == r.T2
    r = cb.expansion(T1, X, PR, eta_poly=eta)
    assert pytest.approx(_ct_polytropic(gas, X, T1, P1, P1 / PR, eta, True), rel=1e-7) == r.T2


def test_the_work_is_cantera_enthalpy_change(gas) -> None:
    X, T1, P1, PR = _humid_air(), 288.15, 101325.0, 12.0
    r = cb.compression(T1, X, PR, eta_poly=0.9)
    gas.TPX = T1, P1, _ct_X(X)
    h1 = gas.enthalpy_mass
    gas.TPX = r.T2, PR * P1, _ct_X(X)
    assert r.dh == pytest.approx(gas.enthalpy_mass - h1, rel=1e-7)
