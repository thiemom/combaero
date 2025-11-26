import math

import numpy as np

from combaero import (
    cp,
    h,
    density,
    speed_of_sound,
    specific_gas_constant,
)
from combaero._core import num_species, species_index_from_name


def make_air_mole_fractions() -> np.ndarray:
    """Return a simple dry air-like composition as mole fractions."""

    n = num_species()
    X = np.zeros(n, dtype=float)

    i_n2 = species_index_from_name("N2")
    i_o2 = species_index_from_name("O2")

    X[i_n2] = 0.79
    X[i_o2] = 0.21
    return X


def test_basic_thermo_for_air() -> None:
    T = 300.0  # K
    P = 101_325.0  # Pa
    X = make_air_mole_fractions()

    cp_val = cp(T, X)
    h_val = h(T, X)
    rho_val = density(T, P, X)
    c_val = speed_of_sound(T, X)
    R_s = specific_gas_constant(X)

    # Basic sanity checks: finite and within broad physical ranges.
    assert math.isfinite(cp_val)
    assert 10.0 < cp_val < 100.0

    assert math.isfinite(h_val)

    assert math.isfinite(rho_val)
    assert 0.1 < rho_val < 10.0

    assert math.isfinite(c_val)
    assert 100.0 < c_val < 2000.0

    assert math.isfinite(R_s)
    assert 100.0 < R_s < 500.0
