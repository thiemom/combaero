"""
Tests for CombustionState dataclass.

Tests verify:
- CombustionState structure and nested CompleteState objects
- combustion_state() function (phi as input)
- combustion_state_from_streams() function (phi as output)
- Consistency with individual function calls
- Various fuel compositions and conditions
"""

import pytest

import combaero as cb


class TestCombustionStateStructure:
    """Test CombustionState dataclass structure."""

    def test_combustion_state_has_all_attributes(self):
        """Test that CombustionState has all required attributes."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

        # Core attributes
        assert hasattr(state, "phi")
        assert hasattr(state, "fuel_name")
        assert hasattr(state, "reactants")
        assert hasattr(state, "products")
        assert hasattr(state, "mixture_fraction")
        assert hasattr(state, "fuel_burn_fraction")

        # Nested CompleteState objects
        assert hasattr(state.reactants, "thermo")
        assert hasattr(state.reactants, "transport")
        assert hasattr(state.products, "thermo")
        assert hasattr(state.products, "transport")

    def test_combustion_state_repr(self):
        """Test CombustionState __repr__ method."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(
            X_CH4, X_air, phi=1.0, T_reactants=300, P=101325, fuel_name="CH4"
        )

        repr_str = repr(state)
        assert "CombustionState" in repr_str
        assert "phi=" in repr_str
        assert "T_reactants=" in repr_str
        assert "T_products=" in repr_str
        assert "fuel='CH4'" in repr_str


class TestCombustionStateFromPhi:
    """Test combustion_state() function (phi as input)."""

    def test_stoichiometric_methane_air(self):
        """Test stoichiometric methane-air combustion."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

        # Check phi is echoed back
        assert state.phi == pytest.approx(1.0, rel=1e-15)

        # Check reactant temperature
        assert pytest.approx(300, rel=1e-15) == state.reactants.thermo.T

        # Check product temperature (adiabatic flame temperature)
        # For stoichiometric CH4-air at 300K, T_ad should be around 2200-2400 K
        assert 2200 < state.products.thermo.T < 2400

        # Check fuel burn fraction (complete combustion)
        assert state.fuel_burn_fraction == pytest.approx(1.0, rel=1e-15)

    def test_lean_combustion(self):
        """Test lean combustion (phi < 1)."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=0.8, T_reactants=300, P=101325)

        assert state.phi == pytest.approx(0.8, rel=1e-15)

        # Lean combustion should have lower T_ad than stoichiometric
        state_stoich = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)
        assert state.products.thermo.T < state_stoich.products.thermo.T

    def test_rich_combustion(self):
        """Test rich combustion (phi > 1)."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.2, T_reactants=300, P=101325)

        assert state.phi == pytest.approx(1.2, rel=1e-15)

        # Rich combustion should have lower T_ad than stoichiometric
        state_stoich = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)
        assert state.products.thermo.T < state_stoich.products.thermo.T

    def test_preheated_reactants(self):
        """Test combustion with preheated reactants."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state_cold = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)
        state_hot = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=600, P=101325)

        # Preheated reactants should give higher T_ad
        assert state_hot.products.thermo.T > state_cold.products.thermo.T

    def test_elevated_pressure(self):
        """Test combustion at elevated pressure."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=1e6)

        # Check pressure is maintained
        assert pytest.approx(1e6, rel=1e-15) == state.reactants.thermo.P
        assert pytest.approx(1e6, rel=1e-15) == state.products.thermo.P

    def test_fuel_name_storage(self):
        """Test that fuel name is stored correctly."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(
            X_CH4, X_air, phi=1.0, T_reactants=300, P=101325, fuel_name="Methane"
        )

        assert state.fuel_name == "Methane"

    def test_mixture_fraction_range(self):
        """Test that mixture fraction is in valid range [0, 1]."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        for phi in [0.5, 0.8, 1.0, 1.2, 1.5]:
            state = cb.combustion_state(X_CH4, X_air, phi=phi, T_reactants=300, P=101325)
            assert 0.0 <= state.mixture_fraction <= 1.0


class TestCombustionStateFromStreams:
    """Test combustion_state_from_streams() function (phi as output)."""

    def test_stoichiometric_from_streams(self):
        """Test stoichiometric combustion from measured streams."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        # Stoichiometric CH4-air: mdot_air/mdot_fuel = 17.2
        fuel = cb.Stream()
        fuel.set_T(300).set_P(101325).set_X(X_CH4).set_mdot(0.01)

        air = cb.Stream()
        air.set_T(298).set_P(101325).set_X(X_air).set_mdot(0.172)

        state = cb.combustion_state_from_streams(fuel, air, fuel_name="CH4")

        # Phi should be computed to be approximately 1.0
        assert state.phi == pytest.approx(1.0, abs=0.05)

        # Check fuel name
        assert state.fuel_name == "CH4"

    def test_lean_from_streams(self):
        """Test lean combustion from measured streams."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        # Excess air -> lean combustion
        fuel = cb.Stream()
        fuel.set_T(300).set_P(101325).set_X(X_CH4).set_mdot(0.01)

        air = cb.Stream()
        air.set_T(298).set_P(101325).set_X(X_air).set_mdot(0.25)  # More air

        state = cb.combustion_state_from_streams(fuel, air)

        # Phi should be less than 1.0
        assert state.phi < 1.0

    def test_rich_from_streams(self):
        """Test rich combustion from measured streams."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        # Less air -> rich combustion
        fuel = cb.Stream()
        fuel.set_T(300).set_P(101325).set_X(X_CH4).set_mdot(0.01)

        air = cb.Stream()
        air.set_T(298).set_P(101325).set_X(X_air).set_mdot(0.12)  # Less air

        state = cb.combustion_state_from_streams(fuel, air)

        # Phi should be greater than 1.0
        assert state.phi > 1.0


class TestNestedCompleteStateAccess:
    """Test access to nested CompleteState properties."""

    def test_reactant_thermo_properties(self):
        """Test access to reactant thermodynamic properties."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

        # All ThermoState properties should be accessible
        assert state.reactants.thermo.T > 0
        assert state.reactants.thermo.P > 0
        assert state.reactants.thermo.rho > 0
        assert state.reactants.thermo.cp > 0
        assert state.reactants.thermo.cv > 0
        assert state.reactants.thermo.gamma > 1.0
        assert state.reactants.thermo.a > 0
        assert state.reactants.thermo.mw > 0

    def test_reactant_transport_properties(self):
        """Test access to reactant transport properties."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

        # All TransportState properties should be accessible
        assert state.reactants.transport.mu > 0
        assert state.reactants.transport.k > 0
        assert state.reactants.transport.nu > 0
        assert state.reactants.transport.alpha > 0
        assert 0.5 < state.reactants.transport.Pr < 1.5

    def test_product_thermo_properties(self):
        """Test access to product thermodynamic properties."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

        # All ThermoState properties should be accessible
        assert state.products.thermo.T > state.reactants.thermo.T  # Higher T after combustion
        assert state.products.thermo.P > 0
        assert state.products.thermo.rho > 0
        assert state.products.thermo.cp > 0
        assert state.products.thermo.cv > 0
        assert state.products.thermo.gamma > 1.0
        assert state.products.thermo.a > 0
        assert state.products.thermo.mw > 0

    def test_product_transport_properties(self):
        """Test access to product transport properties."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

        # All TransportState properties should be accessible
        assert state.products.transport.mu > state.reactants.transport.mu  # Higher mu at high T
        assert state.products.transport.k > 0
        assert state.products.transport.nu > 0
        assert state.products.transport.alpha > 0
        assert 0.5 < state.products.transport.Pr < 1.5


class TestPhysicalConsistency:
    """Test physical consistency of combustion results."""

    def test_energy_conservation(self):
        """Test that enthalpy is conserved (adiabatic combustion)."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

        # For adiabatic combustion: h_reactants = h_products
        # Allow small tolerance due to numerical precision
        assert state.reactants.thermo.h == pytest.approx(state.products.thermo.h, rel=1e-6)

    def test_pressure_maintained(self):
        """Test that pressure is maintained through combustion."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        P = 5e5  # 5 bar
        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=P)

        assert pytest.approx(P, rel=1e-15) == state.reactants.thermo.P
        assert pytest.approx(P, rel=1e-15) == state.products.thermo.P

    def test_temperature_increases(self):
        """Test that temperature increases during combustion."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

        # Products should be hotter than reactants
        assert state.products.thermo.T > state.reactants.thermo.T

        # Temperature rise should be significant (>1000 K for CH4-air)
        assert (state.products.thermo.T - state.reactants.thermo.T) > 1000

    def test_gamma_decreases_at_high_temperature(self):
        """Test that gamma decreases at high temperature (products)."""
        X_CH4 = [0] * 5 + [1.0] + [0] * 8
        X_air = cb.standard_dry_air_composition()

        state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

        # Gamma should be lower for hot products than cold reactants
        assert state.products.thermo.gamma < state.reactants.thermo.gamma
