"""Equilibrium inner solves reach their roots (#927).

The WGS extent was a damped Newton with a 1e-6 finite-difference step: with
CO or H2O a trace the admissible interval is narrower than the step and it
stopped at ln(Q/K) = 12 instead of shifting back (0.015 mole fraction short).
The adiabatic temperature was a damped Newton on an S-shaped residual: a rich
burned state at 1500 K ended at 1390 K, the root being ~924 K. A colder
unreacting inlet came back at the former 300 K bracket end.
"""

import math

import numpy as np
import pytest

import combaero as cb
from combaero._core import num_species
from combaero._core import species_index_from_name as ix

P = 101325.0


def _mix(parts: dict[str, float]) -> np.ndarray:
    X = np.zeros(num_species())
    for name, x in parts.items():
        X[ix(name)] = x
    X[ix("N2")] = 1.0 - sum(parts.values())
    return X


def _ln_q(X: list[float]) -> float:
    return math.log(X[ix("CO2")] * X[ix("H2")] / (X[ix("CO")] * X[ix("H2O")]))


@pytest.mark.parametrize("T", [800.0, 1200.0, 1800.0])
@pytest.mark.parametrize(
    "trace",
    [
        {"CO": 1e-8, "H2O": 0.1, "CO2": 0.1, "H2": 0.1},
        {"CO": 0.1, "H2O": 3e-9, "CO2": 0.1, "H2": 0.1},
    ],
    ids=["CO-trace", "H2O-trace"],
)
def test_wgs_with_a_trace_reactant_reaches_ln_k(T: float, trace: dict[str, float]) -> None:
    """ln Q at the result is ln K, as for a well-conditioned start at the same
    T (main: ln Q = 11.7-11.9, CO stuck at ~8e-7)."""
    ref = cb.wgs_equilibrium(T, _mix({"CO": 0.1, "H2O": 0.1, "CO2": 1e-3, "H2": 0.1}), P)
    got = cb.wgs_equilibrium(T, _mix(trace), P)
    assert _ln_q(list(got.state.X)) == pytest.approx(_ln_q(list(ref.state.X)), abs=1e-9)


def test_wgs_with_an_empty_interval_returns_the_input() -> None:
    """No CO and no CO2: the extent cannot move either way."""
    X = _mix({"H2O": 0.2, "H2": 0.1})
    got = np.asarray(cb.wgs_equilibrium(1200.0, X, P).state.X)
    np.testing.assert_allclose(got, X, rtol=0, atol=1e-15)


def test_rich_adiabatic_reforming_reaches_its_energy_balance() -> None:
    """A burned rich CH4 state at 1500 K (from the equilibrium-combustor fuel
    controller, test_fuel_boundary_892): main stopped at 1390 K."""
    X = np.zeros(num_species())
    for name, x in {
        "N2": 0.4857107358380111,
        "AR": 0.005785232893562374,
        "CO2": 0.06541045578043911,
        "H2O": 0.13667079382315606,
        "CH4": 0.3064227816648314,
    }.items():
        X[ix(name)] = x
    T_in, P_in = 1500.0000001870949, 221116.04839715929
    r = cb.reforming_equilibrium_adiabatic(T_in, X, P_in)
    assert r.converged
    Xo = np.asarray(r.state.X)
    H_in = cb.h(T_in, X) / cb.mwmix(X)
    H_out = cb.h(r.state.T, Xo) / cb.mwmix(Xo)
    assert abs(H_out - H_in) <= 1e-9 * abs(H_in)
    # Cantera HP on the model's species: 923.85 K (thermo data differ by ~0.3 K).
    assert pytest.approx(924.18, abs=0.05) == r.state.T


@pytest.mark.parametrize("T", [250.0, 280.0])
def test_a_cold_unreacting_inlet_keeps_its_temperature(T: float) -> None:
    """Air alone: nothing reacts, so the adiabatic T is the inlet T (main: 300 K)."""
    X = np.asarray(cb.species.dry_air())
    assert pytest.approx(T, abs=1e-9) == cb.combustion_equilibrium(T, X, P).state.T
    assert pytest.approx(T, abs=1e-9) == cb.wgs_equilibrium_adiabatic(T, X, P).state.T
