import numpy as np
import pytest

import combaero as cb


def test_phi_stoich():
    # Methane (CH4) + Air
    # CH4 + 2 O2 -> CO2 + 2 H2O
    # Stoichiometric air: O2_req = 2.0 mol/mol_fuel

    # Standard Air: ~21% O2, ~79% N2
    # In 1 mol of mixture at phi=1:
    # 1 mol CH4 requires 2 mol O2.
    # Total O2 needed = 2.0.
    # If we have 2 mol O2, we are stoich.

    # Composition mapping
    X_air = cb.species.dry_air()
    X_fuel = cb.species.pure_species("CH4")

    phi_target = 1.0
    # X_mix = (phi * r_st * X_fuel + X_ox) / (phi * r_st + 1)
    X_mix = cb.set_equivalence_ratio_mole(phi_target, X_fuel, X_air)

    phi_calc = cb.equivalence_ratio(X_mix)
    assert pytest.approx(phi_calc, rel=2e-3) == phi_target


def test_phi_rich_lean():
    X_air = cb.species.dry_air()
    X_fuel = cb.species.pure_species("CH4")

    for phi_target in [0.5, 0.8, 1.2, 2.0]:
        X_mix = cb.set_equivalence_ratio_mole(phi_target, X_fuel, X_air)
        phi_calc = cb.equivalence_ratio(X_mix)
        assert pytest.approx(phi_calc, rel=2e-3) == phi_target


def test_phi_mass():
    Y_air = cb.species.dry_air_mass()
    Y_fuel = cb.species.to_mass(cb.species.pure_species("CH4"))

    for phi_target in [0.7, 1.0, 1.5]:
        Y_mix = cb.set_equivalence_ratio_mass(phi_target, Y_fuel, Y_air)
        phi_calc = cb.equivalence_ratio_mass(Y_mix)
        assert pytest.approx(phi_calc, rel=2e-3) == phi_target


def test_phi_burned():
    # Burned stoichiometric methane: CH4 + 2 O2 -> CO2 + 2 H2O
    # This should still return Phi = 1.0
    X_burned = np.zeros(len(cb.species.names))
    X_burned[cb.species.indices["CO2"]] = 1.0 / 3.0
    X_burned[cb.species.indices["H2O"]] = 2.0 / 3.0

    phi_calc = cb.equivalence_ratio(X_burned)
    assert pytest.approx(phi_calc, rel=1e-5) == 1.0


def test_phi_guards():
    # Pure O2 (non-combustible) -> Phi = 0.0
    X_O2 = np.zeros(len(cb.species.names))
    X_O2[cb.species.indices["O2"]] = 1.0
    assert cb.equivalence_ratio(X_O2) == 0.0

    # Pure CH4 (no oxygen) -> Phi = 1e9 (very rich)
    X_CH4 = np.zeros(len(cb.species.names))
    X_CH4[cb.species.indices["CH4"]] = 1.0
    assert cb.equivalence_ratio(X_CH4) == 1.0e9


def test_phi_exceptions():
    # This test is now obsolete as we use guards instead of exceptions
    pass


if __name__ == "__main__":
    # Manually run a quick check
    print(
        "Methane Air Phi=1 (Molar):",
        cb.equivalence_ratio(
            cb.set_equivalence_ratio_mole(1.0, cb.species.pure_species("CH4"), cb.species.dry_air())
        ),
    )
    print(
        "Burned Stoich Methane Phi:",
        cb.equivalence_ratio(np.array([0, 0, 0, 1 / 3, 2 / 3, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    )
