import math

import numpy as np

from combaero import (equivalence_ratio_mass, equivalence_ratio_mole,
                      mole_to_mass, set_equivalence_ratio_mass,
                      set_equivalence_ratio_mole)
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


def make_fuel_mole_fractions() -> np.ndarray:
    """Pure CH4 fuel stream in mole fractions."""

    n = num_species()
    X = np.zeros(n, dtype=float)

    i_ch4 = species_index_from_name("CH4")
    X[i_ch4] = 1.0
    return X


def test_equivalence_ratio_mole_round_trip() -> None:
    X_fuel = make_fuel_mole_fractions()
    X_ox = make_air_mole_fractions()

    phis = [0.5, 1.0, 2.0]

    for phi_target in phis:
        X_mix = np.array(set_equivalence_ratio_mole(phi_target, X_fuel, X_ox), dtype=float)

        # Normalization check
        assert math.isfinite(X_mix.sum())
        assert abs(X_mix.sum() - 1.0) < 1e-12

        phi_back = equivalence_ratio_mole(X_mix, X_fuel, X_ox)
        assert math.isfinite(phi_back)
        assert abs(phi_back - phi_target) < 1e-10


def test_equivalence_ratio_mass_round_trip() -> None:
    X_fuel = make_fuel_mole_fractions()
    X_ox = make_air_mole_fractions()

    # Convert to mass fractions
    Y_fuel = np.array(mole_to_mass(X_fuel), dtype=float)
    Y_ox = np.array(mole_to_mass(X_ox), dtype=float)

    phis = [0.5, 1.0, 2.0]

    for phi_target in phis:
        Y_mix = np.array(set_equivalence_ratio_mass(phi_target, Y_fuel, Y_ox), dtype=float)

        # Normalization check
        assert math.isfinite(Y_mix.sum())
        assert abs(Y_mix.sum() - 1.0) < 1e-12

        phi_back = equivalence_ratio_mass(Y_mix, Y_fuel, Y_ox)
        assert math.isfinite(phi_back)
        assert abs(phi_back - phi_target) < 1e-10
