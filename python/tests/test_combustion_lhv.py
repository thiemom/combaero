"""Tests for fuel lower heating value (LHV) helpers."""

from __future__ import annotations

import pytest

import combaero as cb


def _methane_mole_fractions() -> list[float]:
    x = [0.0] * 14
    x[5] = 1.0  # CH4 in thermo_transport_data species ordering
    return x


def test_fuel_lhv_ch4_reasonable_reference_values() -> None:
    """Pure methane LHV should be close to canonical values."""
    x_ch4 = _methane_mole_fractions()

    lhv_molar = cb.fuel_lhv_molar(x_ch4)  # J/mol
    lhv_mass = cb.fuel_lhv_mass(x_ch4)  # J/kg

    assert 7.5e5 < lhv_molar < 8.5e5
    assert 4.8e7 < lhv_mass < 5.2e7


def test_fuel_lhv_mass_molar_consistency() -> None:
    """Mass and molar LHV outputs should be internally consistent."""
    x_ch4 = _methane_mole_fractions()

    lhv_molar = cb.fuel_lhv_molar(x_ch4)
    lhv_mass = cb.fuel_lhv_mass(x_ch4)
    mw = cb.mwmix(x_ch4)  # g/mol

    expected_mass = lhv_molar * 1000.0 / mw
    assert lhv_mass == pytest.approx(expected_mass, rel=1e-12)


def test_fuel_lhv_reference_temperature_argument_changes_result() -> None:
    """Reference temperature input should affect LHV value."""
    x_ch4 = _methane_mole_fractions()

    lhv_298 = cb.fuel_lhv_molar(x_ch4, reference_temperature=298.15)
    lhv_350 = cb.fuel_lhv_molar(x_ch4, reference_temperature=350.0)

    assert lhv_350 != pytest.approx(lhv_298)


def test_fuel_lhv_input_validation() -> None:
    """LHV helpers should validate fuel mole fractions."""
    with pytest.raises(RuntimeError, match="sum"):
        cb.fuel_lhv_molar([0.2] * 14)
