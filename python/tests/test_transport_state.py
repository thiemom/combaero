"""
Tests for TransportState dataclass and transport_state() function.

Tests verify:
- All properties match individual function calls within machine precision
- Consistency between properties (e.g., nu = mu/rho, Pr = mu*cp/k)
- Various gas compositions and conditions
- Edge cases and error handling
"""

import numpy as np
import pytest

import combaero as cb


class TestTransportStateStructure:
    """Test TransportState struct behavior."""

    def test_transport_state_has_all_attributes(self):
        """Verify TransportState has all expected attributes."""
        X = cb.standard_dry_air_composition()
        state = cb.transport_state(T=300, P=101325, X=X)

        # Check all attributes exist
        assert hasattr(state, "T")
        assert hasattr(state, "P")
        assert hasattr(state, "rho")
        assert hasattr(state, "mu")
        assert hasattr(state, "k")
        assert hasattr(state, "nu")
        assert hasattr(state, "alpha")
        assert hasattr(state, "Pr")
        assert hasattr(state, "cp")

    def test_transport_state_readonly(self):
        """Verify TransportState attributes are read-only."""
        X = cb.standard_dry_air_composition()
        state = cb.transport_state(T=300, P=101325, X=X)

        # Attempt to modify should raise AttributeError
        with pytest.raises(AttributeError):
            state.T = 400

        with pytest.raises(AttributeError):
            state.mu = 2e-5

    def test_transport_state_repr(self):
        """Verify TransportState has nice __repr__."""
        X = cb.standard_dry_air_composition()
        state = cb.transport_state(T=300, P=101325, X=X)

        repr_str = repr(state)
        assert "TransportState" in repr_str
        assert "T=" in repr_str
        assert "P=" in repr_str
        assert "mu=" in repr_str


class TestInputEchoing:
    """Test that inputs are echoed back correctly."""

    def test_temperature_echoed(self):
        """Verify temperature is echoed back."""
        X = cb.standard_dry_air_composition()
        T = 450.0
        state = cb.transport_state(T=T, P=101325, X=X)
        assert state.T == T

    def test_pressure_echoed(self):
        """Verify pressure is echoed back."""
        X = cb.standard_dry_air_composition()
        P = 200000.0
        state = cb.transport_state(T=300, P=P, X=X)
        assert state.P == P


class TestPropertyConsistency:
    """Test that all properties match individual function calls."""

    def test_all_properties_match_individual_calls(self):
        """Verify all properties match individual function calls."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325

        state = cb.transport_state(T=T, P=P, X=X)

        # Compare with individual function calls
        assert state.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)
        assert state.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert state.k == pytest.approx(cb.thermal_conductivity(T, P, X), rel=1e-15)
        assert state.nu == pytest.approx(cb.kinematic_viscosity(T, P, X), rel=1e-15)
        assert state.alpha == pytest.approx(cb.thermal_diffusivity(T, P, X), rel=1e-15)
        assert state.Pr == pytest.approx(cb.prandtl(T, P, X), rel=1e-15)
        assert state.cp == pytest.approx(cb.cp_mass(T, X), rel=1e-15)


class TestPhysicalRelationships:
    """Test physical relationships between properties."""

    def test_kinematic_viscosity_equals_mu_over_rho(self):
        """Verify nu = mu / rho."""
        X = cb.standard_dry_air_composition()
        state = cb.transport_state(T=300, P=101325, X=X)

        nu_from_ratio = state.mu / state.rho
        assert state.nu == pytest.approx(nu_from_ratio, rel=1e-14)

    def test_thermal_diffusivity_relationship(self):
        """Verify alpha = k / (rho * cp)."""
        X = cb.standard_dry_air_composition()
        state = cb.transport_state(T=300, P=101325, X=X)

        alpha_from_formula = state.k / (state.rho * state.cp)
        assert state.alpha == pytest.approx(alpha_from_formula, rel=1e-14)

    def test_prandtl_number_relationship(self):
        """Verify Pr = mu * cp / k."""
        X = cb.standard_dry_air_composition()
        state = cb.transport_state(T=300, P=101325, X=X)

        Pr_from_formula = state.mu * state.cp / state.k
        assert state.Pr == pytest.approx(Pr_from_formula, rel=1e-14)


class TestVariousCompositions:
    """Test with various gas compositions."""

    def test_dry_air(self):
        """Test with standard dry air."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325
        state = cb.transport_state(T=T, P=P, X=X)

        # Verify all properties match individual calls
        assert state.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert state.k == pytest.approx(cb.thermal_conductivity(T, P, X), rel=1e-15)
        assert state.Pr == pytest.approx(cb.prandtl(T, P, X), rel=1e-15)

    def test_pure_nitrogen(self):
        """Test with pure nitrogen."""
        X = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        T = 300
        P = 101325
        state = cb.transport_state(T=T, P=P, X=X)

        # Verify all properties match individual calls
        assert state.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert state.k == pytest.approx(cb.thermal_conductivity(T, P, X), rel=1e-15)
        assert state.Pr == pytest.approx(cb.prandtl(T, P, X), rel=1e-15)

    def test_pure_methane(self):
        """Test with pure methane."""
        X = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        T = 300
        P = 101325
        state = cb.transport_state(T=T, P=P, X=X)

        # Verify all properties match individual calls
        assert state.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert state.k == pytest.approx(cb.thermal_conductivity(T, P, X), rel=1e-15)
        assert state.Pr == pytest.approx(cb.prandtl(T, P, X), rel=1e-15)

    def test_combustion_products(self):
        """Test with typical combustion products."""
        X = [0.72, 0.03, 0.0, 0.10, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        T = 1500
        P = 101325
        state = cb.transport_state(T=T, P=P, X=X)

        # Verify all properties match individual calls
        assert state.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert state.k == pytest.approx(cb.thermal_conductivity(T, P, X), rel=1e-15)
        assert state.Pr == pytest.approx(cb.prandtl(T, P, X), rel=1e-15)


class TestVariousConditions:
    """Test at various temperature and pressure conditions."""

    def test_ambient_conditions(self):
        """Test at ambient conditions (20°C, 1 atm)."""
        X = cb.standard_dry_air_composition()
        T = 293.15
        P = 101325
        state = cb.transport_state(T=T, P=P, X=X)

        # Verify all properties match individual calls
        assert state.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert state.k == pytest.approx(cb.thermal_conductivity(T, P, X), rel=1e-15)
        assert state.nu == pytest.approx(cb.kinematic_viscosity(T, P, X), rel=1e-15)

    def test_high_temperature(self):
        """Test at high temperature (combustor conditions)."""
        X = cb.standard_dry_air_composition()
        T = 1500
        P = 101325
        state = cb.transport_state(T=T, P=P, X=X)

        # Verify all properties match individual calls
        assert state.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert state.k == pytest.approx(cb.thermal_conductivity(T, P, X), rel=1e-15)
        assert state.alpha == pytest.approx(cb.thermal_diffusivity(T, P, X), rel=1e-15)

    def test_high_pressure(self):
        """Test at high pressure (compressor discharge)."""
        X = cb.standard_dry_air_composition()
        T = 600
        P = 2e6
        state = cb.transport_state(T=T, P=P, X=X)

        # Verify all properties match individual calls
        assert state.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)
        assert state.nu == pytest.approx(cb.kinematic_viscosity(T, P, X), rel=1e-15)
        assert state.Pr == pytest.approx(cb.prandtl(T, P, X), rel=1e-15)

    def test_low_temperature(self):
        """Test at low temperature (cryogenic)."""
        X = cb.standard_dry_air_composition()
        T = 200
        P = 101325
        state = cb.transport_state(T=T, P=P, X=X)

        # Verify all properties match individual calls
        assert state.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert state.k == pytest.approx(cb.thermal_conductivity(T, P, X), rel=1e-15)


class TestTemperatureDependence:
    """Test temperature dependence of transport properties."""

    def test_viscosity_increases_with_temperature(self):
        """Verify viscosity increases with temperature (Sutherland's law)."""
        X = cb.standard_dry_air_composition()
        state1 = cb.transport_state(T=300, P=101325, X=X)
        state2 = cb.transport_state(T=600, P=101325, X=X)

        assert state2.mu > state1.mu

    def test_thermal_conductivity_increases_with_temperature(self):
        """Verify thermal conductivity increases with temperature."""
        X = cb.standard_dry_air_composition()
        state1 = cb.transport_state(T=300, P=101325, X=X)
        state2 = cb.transport_state(T=600, P=101325, X=X)

        assert state2.k > state1.k

    def test_prandtl_number_temperature_dependence(self):
        """Verify Prandtl number is relatively constant with temperature."""
        X = cb.standard_dry_air_composition()
        state1 = cb.transport_state(T=300, P=101325, X=X)
        state2 = cb.transport_state(T=600, P=101325, X=X)

        # Pr should be relatively constant for air (~0.7-0.8)
        assert 0.6 < state1.Pr < 0.9
        assert 0.6 < state2.Pr < 0.9


class TestPressureDependence:
    """Test pressure dependence of transport properties."""

    def test_density_proportional_to_pressure(self):
        """Verify density is proportional to pressure (ideal gas)."""
        X = cb.standard_dry_air_composition()
        T = 300
        state1 = cb.transport_state(T=T, P=101325, X=X)
        state2 = cb.transport_state(T=T, P=202650, X=X)

        # rho2 / rho1 should be approximately P2 / P1
        ratio = state2.rho / state1.rho
        expected_ratio = 202650 / 101325

        assert ratio == pytest.approx(expected_ratio, rel=1e-10)

    def test_kinematic_viscosity_inversely_proportional_to_pressure(self):
        """Verify nu ~ 1/P at constant temperature."""
        X = cb.standard_dry_air_composition()
        T = 300
        state1 = cb.transport_state(T=T, P=101325, X=X)
        state2 = cb.transport_state(T=T, P=202650, X=X)

        # nu2 / nu1 should be approximately P1 / P2
        ratio = state2.nu / state1.nu
        expected_ratio = 101325 / 202650

        assert ratio == pytest.approx(expected_ratio, rel=1e-10)


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_zero_temperature_raises_error(self):
        """Verify zero temperature raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.transport_state(T=0, P=101325, X=X)

    def test_negative_temperature_raises_error(self):
        """Verify negative temperature raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.transport_state(T=-100, P=101325, X=X)

    def test_zero_pressure_raises_error(self):
        """Verify zero pressure raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.transport_state(T=300, P=0, X=X)

    def test_negative_pressure_raises_error(self):
        """Verify negative pressure raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.transport_state(T=300, P=-101325, X=X)

    def test_empty_composition_raises_error(self):
        """Verify empty composition raises error."""
        with pytest.raises(ValueError):
            cb.transport_state(T=300, P=101325, X=[])


class TestTemperatureRange:
    """Test across valid temperature range."""

    def test_temperature_range_200_to_2000K(self):
        """Test across typical temperature range."""
        X = cb.standard_dry_air_composition()
        temperatures = [200, 300, 500, 1000, 1500, 2000]

        for T in temperatures:
            state = cb.transport_state(T=T, P=101325, X=X)

            # All properties should be finite and positive
            assert np.isfinite(state.mu)
            assert np.isfinite(state.k)
            assert np.isfinite(state.Pr)
            assert state.mu > 0
            assert state.k > 0
            assert state.Pr > 0


class TestPressureRange:
    """Test across valid pressure range."""

    def test_pressure_range_1kPa_to_10MPa(self):
        """Test across typical pressure range."""
        X = cb.standard_dry_air_composition()
        pressures = [1e3, 1e4, 1e5, 1e6, 1e7]

        for P in pressures:
            state = cb.transport_state(T=300, P=P, X=X)

            # All properties should be finite and positive
            assert np.isfinite(state.mu)
            assert np.isfinite(state.k)
            assert np.isfinite(state.nu)
            assert state.mu > 0
            assert state.k > 0
            assert state.nu > 0


class TestConsistencyAcrossConditions:
    """Test consistency of relationships across various conditions."""

    def test_prandtl_always_positive(self):
        """Verify Pr > 0 for all conditions."""
        X = cb.standard_dry_air_composition()
        conditions = [
            (200, 1e5),
            (300, 1e5),
            (1000, 1e5),
            (2000, 1e5),
            (300, 1e6),
        ]

        for T, P in conditions:
            state = cb.transport_state(T=T, P=P, X=X)
            assert state.Pr > 0

    def test_all_properties_positive(self):
        """Verify all properties are positive."""
        X = cb.standard_dry_air_composition()
        state = cb.transport_state(T=300, P=101325, X=X)

        assert state.rho > 0
        assert state.mu > 0
        assert state.k > 0
        assert state.nu > 0
        assert state.alpha > 0
        assert state.Pr > 0
        assert state.cp > 0
