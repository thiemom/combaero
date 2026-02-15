"""
Tests for CompleteState dataclass and complete_state() function.

Tests verify:
- Nested thermo and transport states match individual function calls
- All properties are accessible via nested structure
- Various gas compositions and conditions
"""

import pytest

import combaero as cb


class TestCompleteStateStructure:
    """Test CompleteState struct behavior."""

    def test_complete_state_has_nested_states(self):
        """Verify CompleteState has thermo and transport attributes."""
        X = cb.standard_dry_air_composition()
        state = cb.complete_state(T=300, P=101325, X=X)

        # Check nested states exist
        assert hasattr(state, "thermo")
        assert hasattr(state, "transport")

    def test_thermo_state_has_all_attributes(self):
        """Verify nested ThermoState has all 16 properties."""
        X = cb.standard_dry_air_composition()
        state = cb.complete_state(T=300, P=101325, X=X)

        # Check all thermo attributes
        assert hasattr(state.thermo, "T")
        assert hasattr(state.thermo, "P")
        assert hasattr(state.thermo, "rho")
        assert hasattr(state.thermo, "cp")
        assert hasattr(state.thermo, "cv")
        assert hasattr(state.thermo, "h")
        assert hasattr(state.thermo, "s")
        assert hasattr(state.thermo, "u")
        assert hasattr(state.thermo, "gamma")
        assert hasattr(state.thermo, "a")
        assert hasattr(state.thermo, "cp_mass")
        assert hasattr(state.thermo, "cv_mass")
        assert hasattr(state.thermo, "h_mass")
        assert hasattr(state.thermo, "s_mass")
        assert hasattr(state.thermo, "u_mass")
        assert hasattr(state.thermo, "mw")

    def test_transport_state_has_all_attributes(self):
        """Verify nested TransportState has all 9 properties."""
        X = cb.standard_dry_air_composition()
        state = cb.complete_state(T=300, P=101325, X=X)

        # Check all transport attributes
        assert hasattr(state.transport, "T")
        assert hasattr(state.transport, "P")
        assert hasattr(state.transport, "rho")
        assert hasattr(state.transport, "mu")
        assert hasattr(state.transport, "k")
        assert hasattr(state.transport, "nu")
        assert hasattr(state.transport, "alpha")
        assert hasattr(state.transport, "Pr")
        assert hasattr(state.transport, "cp")

    def test_complete_state_readonly(self):
        """Verify CompleteState attributes are read-only."""
        X = cb.standard_dry_air_composition()
        state = cb.complete_state(T=300, P=101325, X=X)

        # Attempt to modify should raise AttributeError
        with pytest.raises(AttributeError):
            state.thermo = None

        with pytest.raises(AttributeError):
            state.transport = None

    def test_complete_state_repr(self):
        """Verify CompleteState has nice __repr__."""
        X = cb.standard_dry_air_composition()
        state = cb.complete_state(T=300, P=101325, X=X)

        repr_str = repr(state)
        assert "CompleteState" in repr_str
        assert "T=" in repr_str
        assert "P=" in repr_str


class TestConsistencyWithIndividualFunctions:
    """Test that nested states match individual function calls."""

    def test_thermo_matches_thermo_state(self):
        """Verify state.thermo matches thermo_state() call."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325

        complete = cb.complete_state(T=T, P=P, X=X)
        thermo = cb.thermo_state(T=T, P=P, X=X)

        # All properties should match
        assert pytest.approx(thermo.T, rel=1e-15) == complete.thermo.T
        assert pytest.approx(thermo.P, rel=1e-15) == complete.thermo.P
        assert complete.thermo.rho == pytest.approx(thermo.rho, rel=1e-15)
        assert complete.thermo.cp == pytest.approx(thermo.cp, rel=1e-15)
        assert complete.thermo.cv == pytest.approx(thermo.cv, rel=1e-15)
        assert complete.thermo.h == pytest.approx(thermo.h, rel=1e-15)
        assert complete.thermo.s == pytest.approx(thermo.s, rel=1e-15)
        assert complete.thermo.u == pytest.approx(thermo.u, rel=1e-15)
        assert complete.thermo.gamma == pytest.approx(thermo.gamma, rel=1e-15)
        assert complete.thermo.a == pytest.approx(thermo.a, rel=1e-15)
        assert complete.thermo.mw == pytest.approx(thermo.mw, rel=1e-15)

    def test_transport_matches_transport_state(self):
        """Verify state.transport matches transport_state() call."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325

        complete = cb.complete_state(T=T, P=P, X=X)
        transport = cb.transport_state(T=T, P=P, X=X)

        # All properties should match
        assert pytest.approx(transport.T, rel=1e-15) == complete.transport.T
        assert pytest.approx(transport.P, rel=1e-15) == complete.transport.P
        assert complete.transport.rho == pytest.approx(transport.rho, rel=1e-15)
        assert complete.transport.mu == pytest.approx(transport.mu, rel=1e-15)
        assert complete.transport.k == pytest.approx(transport.k, rel=1e-15)
        assert complete.transport.nu == pytest.approx(transport.nu, rel=1e-15)
        assert complete.transport.alpha == pytest.approx(transport.alpha, rel=1e-15)
        assert complete.transport.Pr == pytest.approx(transport.Pr, rel=1e-15)
        assert complete.transport.cp == pytest.approx(transport.cp, rel=1e-15)


class TestVariousCompositions:
    """Test with various gas compositions."""

    def test_dry_air(self):
        """Test with standard dry air."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325
        state = cb.complete_state(T=T, P=P, X=X)

        # Verify matches individual calls
        thermo = cb.thermo_state(T=T, P=P, X=X)
        transport = cb.transport_state(T=T, P=P, X=X)

        assert state.thermo.h == pytest.approx(thermo.h, rel=1e-15)
        assert state.transport.mu == pytest.approx(transport.mu, rel=1e-15)

    def test_pure_nitrogen(self):
        """Test with pure nitrogen."""
        X = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        T = 300
        P = 101325
        state = cb.complete_state(T=T, P=P, X=X)

        # Verify matches individual calls
        thermo = cb.thermo_state(T=T, P=P, X=X)
        transport = cb.transport_state(T=T, P=P, X=X)

        assert state.thermo.mw == pytest.approx(thermo.mw, rel=1e-15)
        assert state.transport.Pr == pytest.approx(transport.Pr, rel=1e-15)

    def test_pure_methane(self):
        """Test with pure methane."""
        X = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        T = 300
        P = 101325
        state = cb.complete_state(T=T, P=P, X=X)

        # Verify matches individual calls
        thermo = cb.thermo_state(T=T, P=P, X=X)
        transport = cb.transport_state(T=T, P=P, X=X)

        assert state.thermo.gamma == pytest.approx(thermo.gamma, rel=1e-15)
        assert state.transport.k == pytest.approx(transport.k, rel=1e-15)


class TestVariousConditions:
    """Test at various temperature and pressure conditions."""

    def test_ambient_conditions(self):
        """Test at ambient conditions (20°C, 1 atm)."""
        X = cb.standard_dry_air_composition()
        T = 293.15
        P = 101325
        state = cb.complete_state(T=T, P=P, X=X)

        # Verify matches individual calls
        thermo = cb.thermo_state(T=T, P=P, X=X)
        transport = cb.transport_state(T=T, P=P, X=X)

        assert state.thermo.rho == pytest.approx(thermo.rho, rel=1e-15)
        assert state.transport.nu == pytest.approx(transport.nu, rel=1e-15)

    def test_high_temperature(self):
        """Test at high temperature (combustor conditions)."""
        X = cb.standard_dry_air_composition()
        T = 1500
        P = 101325
        state = cb.complete_state(T=T, P=P, X=X)

        # Verify matches individual calls
        thermo = cb.thermo_state(T=T, P=P, X=X)
        transport = cb.transport_state(T=T, P=P, X=X)

        assert state.thermo.a == pytest.approx(thermo.a, rel=1e-15)
        assert state.transport.mu == pytest.approx(transport.mu, rel=1e-15)

    def test_high_pressure(self):
        """Test at high pressure (compressor discharge)."""
        X = cb.standard_dry_air_composition()
        T = 600
        P = 2e6
        state = cb.complete_state(T=T, P=P, X=X)

        # Verify matches individual calls
        thermo = cb.thermo_state(T=T, P=P, X=X)
        transport = cb.transport_state(T=T, P=P, X=X)

        assert state.thermo.rho == pytest.approx(thermo.rho, rel=1e-15)
        assert state.transport.Pr == pytest.approx(transport.Pr, rel=1e-15)


class TestReferencePresssure:
    """Test entropy with different reference pressures."""

    def test_default_reference_pressure(self):
        """Test with default reference pressure (1 atm)."""
        X = cb.standard_dry_air_composition()
        state = cb.complete_state(T=300, P=101325, X=X)

        # Should match thermo_state with default P_ref
        thermo = cb.thermo_state(T=300, P=101325, X=X)
        assert state.thermo.s == pytest.approx(thermo.s, rel=1e-15)

    def test_custom_reference_pressure(self):
        """Test with custom reference pressure."""
        X = cb.standard_dry_air_composition()
        P_ref = 200000.0
        state = cb.complete_state(T=300, P=101325, X=X, P_ref=P_ref)

        # Should match thermo_state with same P_ref
        thermo = cb.thermo_state(T=300, P=101325, X=X, P_ref=P_ref)
        assert state.thermo.s == pytest.approx(thermo.s, rel=1e-15)


class TestConsistencyBetweenNestedStates:
    """Test consistency between thermo and transport nested states."""

    def test_temperature_matches(self):
        """Verify T is same in both nested states."""
        X = cb.standard_dry_air_composition()
        state = cb.complete_state(T=300, P=101325, X=X)

        assert state.thermo.T == state.transport.T

    def test_pressure_matches(self):
        """Verify P is same in both nested states."""
        X = cb.standard_dry_air_composition()
        state = cb.complete_state(T=300, P=101325, X=X)

        assert state.thermo.P == state.transport.P

    def test_density_matches(self):
        """Verify rho is same in both nested states."""
        X = cb.standard_dry_air_composition()
        state = cb.complete_state(T=300, P=101325, X=X)

        assert state.thermo.rho == pytest.approx(state.transport.rho, rel=1e-15)


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_zero_temperature_raises_error(self):
        """Verify zero temperature raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.complete_state(T=0, P=101325, X=X)

    def test_negative_temperature_raises_error(self):
        """Verify negative temperature raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.complete_state(T=-100, P=101325, X=X)

    def test_zero_pressure_raises_error(self):
        """Verify zero pressure raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.complete_state(T=300, P=0, X=X)

    def test_negative_pressure_raises_error(self):
        """Verify negative pressure raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.complete_state(T=300, P=-101325, X=X)

    def test_empty_composition_raises_error(self):
        """Verify empty composition raises error."""
        with pytest.raises(ValueError):
            cb.complete_state(T=300, P=101325, X=[])
