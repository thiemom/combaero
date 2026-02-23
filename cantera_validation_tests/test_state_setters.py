import os

import cantera as ct
import numpy as np
import pytest

import combaero as cb

CB_SPECIES = [
    "N2",
    "O2",
    "AR",
    "CO2",
    "H2O",
    "CH4",
    "C2H6",
    "C3H8",
    "IC4H10",
    "NC5H12",
    "NC6H14",
    "NC7H16",
    "CO",
    "H2",
]


@pytest.fixture(scope="module")
def gas():
    """Returns a Cantera gas object initialized with CombAero's NASA9 database."""
    yaml_path = os.path.join(os.path.dirname(__file__), "combaero_nasa9.yaml")
    gas = ct.Solution(yaml_path)
    return gas


def get_random_state(gas):
    """Generate a random valid state for testing."""
    T = np.random.uniform(300, 2500)
    P = np.random.uniform(0.1 * ct.one_atm, 50 * ct.one_atm)

    cb_X = np.random.rand(len(CB_SPECIES))
    for i, name in enumerate(CB_SPECIES):
        if name not in gas.species_names:
            cb_X[i] = 0.0

    cb_X /= np.sum(cb_X)
    ct_X = {name: cb_X[i] for i, name in enumerate(CB_SPECIES) if cb_X[i] > 0}

    return T, P, cb_X, ct_X


def assert_state_equal(gas, cb_state, tol=1e-5):
    """Assert that the Cantera gas object and combaero state are equivalent."""
    assert np.isclose(gas.T, cb_state.T, rtol=tol), f"T mismatch: {gas.T} vs {cb_state.T}"
    assert np.isclose(gas.P, cb_state.P, rtol=tol), f"P mismatch: {gas.P} vs {cb_state.P}"

    for i, name in enumerate(CB_SPECIES):
        ct_val = gas[name].X[0] if name in gas.species_names else 0.0
        cb_val = cb_state.X[i]
        assert np.isclose(ct_val, cb_val, atol=1e-5), f"X mismatch for {name}: {ct_val} vs {cb_val}"


@pytest.mark.parametrize("i", range(50))
def test_set_TPX(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X

    cb_state = cb.State()
    cb_state.set_TPX(T, P, cb_X)

    assert_state_equal(gas, cb_state)


@pytest.mark.parametrize("i", range(50))
def test_set_TPY(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    # Need to get mass fractions in the format CombAero expects
    cb_Y = np.zeros(len(CB_SPECIES))
    for j, name in enumerate(CB_SPECIES):
        if name in gas.species_names:
            cb_Y[j] = gas[name].Y[0]

    cb_state = cb.State()
    cb_state.set_TPY(T, P, cb_Y)

    assert_state_equal(gas, cb_state, tol=1e-5)


@pytest.mark.parametrize("i", range(50))
def test_property_HP(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    target_h = gas.enthalpy_mass

    cb_state = cb.State()
    cb_state.set_TPX(300, 101325, cb_X)
    cb_state.HP = (target_h, P)

    assert_state_equal(gas, cb_state, tol=1e-3)


@pytest.mark.parametrize("i", range(50))
def test_property_SP(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    target_s = gas.entropy_mass

    cb_state = cb.State()
    cb_state.set_TPX(300, 101325, cb_X)
    cb_state.SP = (target_s, P)

    assert_state_equal(gas, cb_state, tol=1e-3)


@pytest.mark.parametrize("i", range(50))
def test_property_UV(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    target_u = gas.int_energy_mass
    target_v = gas.volume_mass

    cb_state = cb.State()
    cb_state.set_TPX(300, 101325, cb_X)
    cb_state.UV = (target_u, target_v)

    assert_state_equal(gas, cb_state, tol=1e-3)


@pytest.mark.parametrize("i", range(50))
def test_property_SV(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    target_s = gas.entropy_mass
    target_v = gas.volume_mass

    cb_state = cb.State()
    cb_state.set_TPX(300, 101325, cb_X)
    cb_state.SV = (target_s, target_v)

    assert_state_equal(gas, cb_state, tol=1e-3)


@pytest.mark.parametrize("i", range(50))
def test_property_PV(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    target_v = gas.volume_mass

    cb_state = cb.State()
    cb_state.set_TPX(300, 101325, cb_X)
    cb_state.PV = (P, target_v)

    assert_state_equal(gas, cb_state, tol=1e-3)


@pytest.mark.parametrize("i", range(50))
def test_property_UP(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    target_u = gas.int_energy_mass

    cb_state = cb.State()
    cb_state.set_TPX(300, 101325, cb_X)
    cb_state.UP = (target_u, P)

    assert_state_equal(gas, cb_state, tol=1e-3)


@pytest.mark.parametrize("i", range(50))
def test_property_VH(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    target_v = gas.volume_mass
    target_h = gas.enthalpy_mass

    cb_state = cb.State()
    cb_state.set_TPX(300, 101325, cb_X)
    cb_state.VH = (target_v, target_h)

    assert_state_equal(gas, cb_state, tol=1e-3)


@pytest.mark.parametrize("i", range(50))
def test_property_SH(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    target_h = gas.enthalpy_mass
    target_s = gas.entropy_mass

    cb_state = cb.State()
    cb_state.set_TPX(300, 101325, cb_X)
    cb_state.SH = (target_s, target_h)

    assert_state_equal(gas, cb_state, tol=1e-3)


@pytest.mark.parametrize("i", range(50))
def test_property_DP(gas, i):
    T, P, cb_X, ct_X = get_random_state(gas)

    gas.TPX = T, P, ct_X
    target_rho = gas.density

    cb_state = cb.State()
    cb_state.set_TPX(300, 101325, cb_X)
    cb_state.DP = (target_rho, P)

    assert_state_equal(gas, cb_state, tol=1e-3)
