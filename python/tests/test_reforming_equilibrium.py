"""Tests for reforming equilibrium functions."""

from __future__ import annotations

import pytest

import combaero as ca
from combaero.species import SpeciesLocator


@pytest.fixture
def sp() -> SpeciesLocator:
    """Species locator fixture."""
    return SpeciesLocator.from_core()


class TestReformingEquilibrium:
    """Tests for reforming_equilibrium and reforming_equilibrium_adiabatic."""

    def test_reforming_equilibrium_isothermal(self, sp: SpeciesLocator) -> None:
        """Test isothermal reforming equilibrium with CH4."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        T_in = 2000.0
        result = ca.reforming_equilibrium(T_in, X)

        # Temperature should be unchanged (isothermal)
        assert pytest.approx(T_in, rel=1e-6) == result.T

        # CH4 should be mostly reformed at high T
        assert result.X[sp.indices["CH4"]] < X[sp.indices["CH4"]] * 0.5

        # CO and H2 should be produced
        assert result.X[sp.indices["CO"]] > 0.01
        assert result.X[sp.indices["H2"]] > 0.01

    def test_reforming_equilibrium_adiabatic(self, sp: SpeciesLocator) -> None:
        """Test adiabatic reforming equilibrium - temperature should drop."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        T_in = 2200.0
        result = ca.reforming_equilibrium_adiabatic(T_in, X)

        # Temperature should decrease (endothermic reforming)
        assert T_in > result.T
        assert result.T > 1500.0  # But still reasonable

        # CH4 should be reformed
        assert result.X[sp.indices["CH4"]] < X[sp.indices["CH4"]] * 0.5

    def test_reforming_multiple_hydrocarbons(self, sp: SpeciesLocator) -> None:
        """Test reforming with multiple hydrocarbons (natural gas)."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.015
        X[sp.indices["C2H6"]] = 0.003
        X[sp.indices["C3H8"]] = 0.002
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        result = ca.reforming_equilibrium(2000.0, X)

        # All hydrocarbons should be reformed
        assert result.X[sp.indices["CH4"]] < X[sp.indices["CH4"]] * 0.5
        assert result.X[sp.indices["C2H6"]] < X[sp.indices["C2H6"]] * 0.5
        assert result.X[sp.indices["C3H8"]] < X[sp.indices["C3H8"]] * 0.5

        # CO and H2 should be produced
        assert result.X[sp.indices["CO"]] > 0.01
        assert result.X[sp.indices["H2"]] > 0.01

    def test_reforming_c2plus_only(self, sp: SpeciesLocator) -> None:
        """Test reforming with only C2+ hydrocarbons (no CH4)."""
        X = sp.empty()
        X[sp.indices["C2H6"]] = 0.010
        X[sp.indices["C3H8"]] = 0.005
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.085
        X[sp.indices["N2"]] = 0.70

        result = ca.reforming_equilibrium(2000.0, X)

        # C2H6 and C3H8 should be reformed even without CH4
        assert result.X[sp.indices["C2H6"]] < X[sp.indices["C2H6"]] * 0.5
        assert result.X[sp.indices["C3H8"]] < X[sp.indices["C3H8"]] * 0.5

        # CO and H2 should be produced
        assert result.X[sp.indices["CO"]] > 0.01
        assert result.X[sp.indices["H2"]] > 0.01

    def test_element_conservation(self, sp: SpeciesLocator) -> None:
        """Test that element ratios are conserved (using N2 as reference)."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.015
        X[sp.indices["C2H6"]] = 0.003
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70
        X[sp.indices["AR"]] = 0.002

        result = ca.reforming_equilibrium(2000.0, X)

        # N2 is inert, so element ratios per N2 should be conserved
        n2_in = X[sp.indices["N2"]]
        n2_out = result.X[sp.indices["N2"]]

        # C atoms per N2: CH4 + 2*C2H6 + CO + CO2
        C_per_N2_in = (
            X[sp.indices["CH4"]]
            + 2 * X[sp.indices["C2H6"]]
            + X[sp.indices["CO"]]
            + X[sp.indices["CO2"]]
        ) / n2_in
        C_per_N2_out = (
            result.X[sp.indices["CH4"]]
            + 2 * result.X[sp.indices["C2H6"]]
            + result.X[sp.indices["CO"]]
            + result.X[sp.indices["CO2"]]
        ) / n2_out
        assert C_per_N2_in == pytest.approx(C_per_N2_out, rel=1e-5)

        # Ar per N2 should be unchanged
        ar_per_n2_in = X[sp.indices["AR"]] / n2_in
        ar_per_n2_out = result.X[sp.indices["AR"]] / n2_out
        assert ar_per_n2_in == pytest.approx(ar_per_n2_out, rel=1e-10)


class TestSmrWgsEquilibrium:
    """Tests for SMR+WGS equilibrium (CH4 only)."""

    def test_smr_wgs_isothermal(self, sp: SpeciesLocator) -> None:
        """Test isothermal SMR+WGS equilibrium."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        result = ca.smr_wgs_equilibrium(2000.0, X)

        # CH4 should be reformed
        assert result.X[sp.indices["CH4"]] < X[sp.indices["CH4"]] * 0.5

        # CO and H2 should be produced
        assert result.X[sp.indices["CO"]] > 0.01
        assert result.X[sp.indices["H2"]] > 0.01

    def test_smr_wgs_adiabatic(self, sp: SpeciesLocator) -> None:
        """Test adiabatic SMR+WGS equilibrium."""
        X = sp.empty()
        X[sp.indices["CH4"]] = 0.02
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.08
        X[sp.indices["N2"]] = 0.70

        T_in = 2200.0
        result = ca.smr_wgs_equilibrium_adiabatic(T_in, X)

        # Temperature should decrease
        assert T_in > result.T

    def test_smr_wgs_fallback_to_wgs(self, sp: SpeciesLocator) -> None:
        """Test that SMR+WGS falls back to WGS when no CH4 present."""
        X = sp.empty()
        X[sp.indices["CO"]] = 0.05
        X[sp.indices["H2O"]] = 0.20
        X[sp.indices["CO2"]] = 0.05
        X[sp.indices["N2"]] = 0.70

        result = ca.smr_wgs_equilibrium(2000.0, X)

        # Should still work (WGS only)
        assert pytest.approx(2000.0, rel=1e-6) == result.T
        # WGS shifts CO + H2O -> CO2 + H2
        assert result.X[sp.indices["H2"]] > 0.0
