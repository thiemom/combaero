"""
Tests for ThermoState dataclass and thermo_state() function.

Tests verify:
- All properties match individual function calls within machine precision
- Consistency between molar and mass-specific properties
- Physical relationships (e.g., h = u + P*V_molar)
- Various gas compositions and conditions
- Edge cases and error handling
"""

import numpy as np
import pytest

import combaero as cb


class TestThermoStateStructure:
    """Test ThermoState struct behavior."""

    def test_thermo_state_has_all_attributes(self):
        """Verify ThermoState has all expected attributes."""
        X = cb.standard_dry_air_composition()
        state = cb.thermo_state(T=300, P=101325, X=X)

        # Check all attributes exist
        assert hasattr(state, "T")
        assert hasattr(state, "P")
        assert hasattr(state, "rho")
        assert hasattr(state, "cp")
        assert hasattr(state, "cv")
        assert hasattr(state, "h")
        assert hasattr(state, "s")
        assert hasattr(state, "u")
        assert hasattr(state, "gamma")
        assert hasattr(state, "a")
        assert hasattr(state, "cp_mass")
        assert hasattr(state, "cv_mass")
        assert hasattr(state, "h_mass")
        assert hasattr(state, "s_mass")
        assert hasattr(state, "u_mass")
        assert hasattr(state, "mw")

    def test_thermo_state_readonly(self):
        """Verify ThermoState attributes are read-only."""
        X = cb.standard_dry_air_composition()
        state = cb.thermo_state(T=300, P=101325, X=X)

        # Attempt to modify should raise AttributeError
        with pytest.raises(AttributeError):
            state.T = 400

        with pytest.raises(AttributeError):
            state.rho = 2.0

    def test_thermo_state_repr(self):
        """Verify ThermoState has nice __repr__."""
        X = cb.standard_dry_air_composition()
        state = cb.thermo_state(T=300, P=101325, X=X)

        repr_str = repr(state)
        assert "ThermoState" in repr_str
        assert "T=" in repr_str
        assert "P=" in repr_str
        assert "rho=" in repr_str


class TestInputEchoing:
    """Test that inputs are echoed back correctly."""

    def test_temperature_echoed(self):
        """Verify temperature is echoed back."""
        X = cb.standard_dry_air_composition()
        T = 450.0
        state = cb.thermo_state(T=T, P=101325, X=X)
        assert state.T == T

    def test_pressure_echoed(self):
        """Verify pressure is echoed back."""
        X = cb.standard_dry_air_composition()
        P = 200000.0
        state = cb.thermo_state(T=300, P=P, X=X)
        assert state.P == P


class TestPropertyConsistency:
    """Test that all properties match individual function calls."""

    def test_molar_properties_match_individual_calls(self):
        """Verify molar properties match individual function calls."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325

        state = cb.thermo_state(T=T, P=P, X=X)

        # Compare with individual function calls
        assert state.cp == pytest.approx(cb.cp(T, X), rel=1e-15)
        assert state.cv == pytest.approx(cb.cv(T, X), rel=1e-15)
        assert state.h == pytest.approx(cb.h(T, X), rel=1e-15)
        assert state.s == pytest.approx(cb.s(T, X, P), rel=1e-15)
        assert state.u == pytest.approx(cb.u(T, X), rel=1e-15)
        assert state.gamma == pytest.approx(cb.isentropic_expansion_coefficient(T, X), rel=1e-15)
        assert state.a == pytest.approx(cb.speed_of_sound(T, X), rel=1e-15)

    def test_mass_properties_match_individual_calls(self):
        """Verify mass-specific properties match individual function calls."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325

        state = cb.thermo_state(T=T, P=P, X=X)

        # Compare with individual function calls
        assert state.cp_mass == pytest.approx(cb.cp_mass(T, X), rel=1e-15)
        assert state.cv_mass == pytest.approx(cb.cv_mass(T, X), rel=1e-15)
        assert state.h_mass == pytest.approx(cb.h_mass(T, X), rel=1e-15)
        assert state.s_mass == pytest.approx(cb.s_mass(T, X, P), rel=1e-15)
        assert state.u_mass == pytest.approx(cb.u_mass(T, X), rel=1e-15)

    def test_density_matches_individual_call(self):
        """Verify density matches individual function call."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325

        state = cb.thermo_state(T=T, P=P, X=X)

        assert state.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)

    def test_molecular_weight_matches_individual_call(self):
        """Verify molecular weight matches individual function call."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325

        state = cb.thermo_state(T=T, P=P, X=X)

        assert state.mw == pytest.approx(cb.mwmix(X), rel=1e-15)


class TestPhysicalRelationships:
    """Test physical relationships between properties."""

    def test_gamma_equals_cp_over_cv(self):
        """Verify gamma = cp / cv."""
        X = cb.standard_dry_air_composition()
        state = cb.thermo_state(T=300, P=101325, X=X)

        gamma_from_ratio = state.cp / state.cv
        assert state.gamma == pytest.approx(gamma_from_ratio, rel=1e-14)

    def test_mass_molar_conversion(self):
        """Verify mass-specific = molar / (mw/1000)."""
        X = cb.standard_dry_air_composition()
        state = cb.thermo_state(T=300, P=101325, X=X)

        mw_kg_per_mol = state.mw / 1000.0  # Convert g/mol to kg/mol

        # cp_mass = cp / mw_kg_per_mol
        assert state.cp_mass == pytest.approx(state.cp / mw_kg_per_mol, rel=1e-14)
        assert state.cv_mass == pytest.approx(state.cv / mw_kg_per_mol, rel=1e-14)
        assert state.h_mass == pytest.approx(state.h / mw_kg_per_mol, rel=1e-14)
        assert state.s_mass == pytest.approx(state.s / mw_kg_per_mol, rel=1e-14)
        assert state.u_mass == pytest.approx(state.u / mw_kg_per_mol, rel=1e-14)

    def test_ideal_gas_law(self):
        """Verify ideal gas law: P = rho * R_specific * T."""
        X = cb.standard_dry_air_composition()
        state = cb.thermo_state(T=300, P=101325, X=X)

        R_specific = cb.specific_gas_constant(X)
        P_from_eos = state.rho * R_specific * state.T

        assert pytest.approx(P_from_eos, rel=1e-14) == state.P

    def test_enthalpy_internal_energy_relationship(self):
        """Verify h = u + P*V_molar where V_molar = RT/P."""
        X = cb.standard_dry_air_composition()
        state = cb.thermo_state(T=300, P=101325, X=X)

        R = 8.31446261815324  # Universal gas constant [J/(mol*K)]
        V_molar = R * state.T / state.P
        h_from_u = state.u + state.P * V_molar

        assert state.h == pytest.approx(h_from_u, rel=1e-14)


class TestVariousCompositions:
    """Test with various gas compositions."""

    def test_dry_air(self):
        """Test with standard dry air."""
        X = cb.standard_dry_air_composition()
        T = 300
        P = 101325
        state = cb.thermo_state(T=T, P=P, X=X)

        # Verify properties match individual calls
        assert state.mw == pytest.approx(cb.mwmix(X), rel=1e-15)
        assert state.gamma == pytest.approx(cb.isentropic_expansion_coefficient(T, X), rel=1e-15)

    def test_pure_nitrogen(self):
        """Test with pure nitrogen."""
        # Create pure N2 composition (14 species, N2 is index 0)
        X = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        T = 300
        P = 101325
        state = cb.thermo_state(T=T, P=P, X=X)

        # Verify properties match individual calls
        assert state.mw == pytest.approx(cb.mwmix(X), rel=1e-15)
        assert state.gamma == pytest.approx(cb.isentropic_expansion_coefficient(T, X), rel=1e-15)

    def test_pure_methane(self):
        """Test with pure methane."""
        # Create pure CH4 composition (14 species, CH4 is index 5)
        X = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        T = 300
        P = 101325
        state = cb.thermo_state(T=T, P=P, X=X)

        # Verify properties match individual calls
        assert state.mw == pytest.approx(cb.mwmix(X), rel=1e-15)

    def test_combustion_products(self):
        """Test with typical combustion products."""
        # Approximate post-combustion mixture
        # Species order: N2, O2, AR, CO2, H2O, CH4, C2H6, C3H8, IC4H10, NC5H12, NC6H14, NC7H16, CO, H2
        X = [0.72, 0.03, 0.0, 0.10, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        T = 1500
        P = 101325
        state = cb.thermo_state(T=T, P=P, X=X)

        # Verify properties match individual calls
        assert state.mw == pytest.approx(cb.mwmix(X), rel=1e-15)
        assert state.gamma == pytest.approx(cb.isentropic_expansion_coefficient(T, X), rel=1e-15)


class TestVariousConditions:
    """Test at various temperature and pressure conditions."""

    def test_ambient_conditions(self):
        """Test at ambient conditions (20degC, 1 atm)."""
        X = cb.standard_dry_air_composition()
        T = 293.15
        P = 101325
        state = cb.thermo_state(T=T, P=P, X=X)

        # Verify properties match individual calls
        assert state.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)
        assert state.a == pytest.approx(cb.speed_of_sound(T, X), rel=1e-15)

    def test_high_temperature(self):
        """Test at high temperature (combustor conditions)."""
        X = cb.standard_dry_air_composition()
        T = 1500
        P = 101325
        state = cb.thermo_state(T=T, P=P, X=X)

        # Verify properties match individual calls
        assert state.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)
        assert state.a == pytest.approx(cb.speed_of_sound(T, X), rel=1e-15)

    def test_high_pressure(self):
        """Test at high pressure (compressor discharge)."""
        X = cb.standard_dry_air_composition()
        T = 600
        P = 2e6
        state = cb.thermo_state(T=T, P=P, X=X)

        # Verify properties match individual calls
        assert state.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)

    def test_low_temperature(self):
        """Test at low temperature (cryogenic)."""
        X = cb.standard_dry_air_composition()
        T = 200
        P = 101325
        state = cb.thermo_state(T=T, P=P, X=X)

        # Verify properties match individual calls
        assert state.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)


class TestReferencePresssure:
    """Test entropy with different reference pressures."""

    def test_default_reference_pressure(self):
        """Test with default reference pressure (1 atm)."""
        X = cb.standard_dry_air_composition()
        state = cb.thermo_state(T=300, P=101325, X=X)

        # Should match individual call with default P_ref
        assert state.s == pytest.approx(cb.s(300, X, 101325), rel=1e-15)

    def test_custom_reference_pressure(self):
        """Test with custom reference pressure."""
        X = cb.standard_dry_air_composition()
        P_ref = 200000.0
        state = cb.thermo_state(T=300, P=101325, X=X, P_ref=P_ref)

        # Should match individual call with custom P_ref
        assert state.s == pytest.approx(cb.s(300, X, 101325, P_ref), rel=1e-15)

    def test_entropy_pressure_dependence(self):
        """Test that entropy depends on reference pressure."""
        X = cb.standard_dry_air_composition()
        state1 = cb.thermo_state(T=300, P=101325, X=X, P_ref=101325)
        state2 = cb.thermo_state(T=300, P=101325, X=X, P_ref=200000)

        # Entropy should differ with different reference pressures
        assert state1.s != pytest.approx(state2.s, rel=1e-10)


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_zero_temperature_raises_error(self):
        """Verify zero temperature raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.thermo_state(T=0, P=101325, X=X)

    def test_negative_temperature_raises_error(self):
        """Verify negative temperature raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.thermo_state(T=-100, P=101325, X=X)

    def test_zero_pressure_raises_error(self):
        """Verify zero pressure raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.thermo_state(T=300, P=0, X=X)

    def test_negative_pressure_raises_error(self):
        """Verify negative pressure raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.thermo_state(T=300, P=-101325, X=X)

    def test_empty_composition_raises_error(self):
        """Verify empty composition raises error."""
        with pytest.raises(ValueError):
            cb.thermo_state(T=300, P=101325, X=[])

    def test_zero_reference_pressure_raises_error(self):
        """Verify zero reference pressure raises error."""
        X = cb.standard_dry_air_composition()
        with pytest.raises(ValueError):
            cb.thermo_state(T=300, P=101325, X=X, P_ref=0)


class TestTemperatureRange:
    """Test across valid temperature range."""

    def test_temperature_range_200_to_2000K(self):
        """Test across typical temperature range."""
        X = cb.standard_dry_air_composition()
        temperatures = [200, 300, 500, 1000, 1500, 2000]

        for T in temperatures:
            state = cb.thermo_state(T=T, P=101325, X=X)

            # All properties should be finite and reasonable
            assert np.isfinite(state.rho)
            assert np.isfinite(state.cp)
            assert np.isfinite(state.h)
            assert state.rho > 0
            assert state.cp > 0
            assert state.gamma > 1.0


class TestPressureRange:
    """Test across valid pressure range."""

    def test_pressure_range_1kPa_to_10MPa(self):
        """Test across typical pressure range."""
        X = cb.standard_dry_air_composition()
        pressures = [1e3, 1e4, 1e5, 1e6, 1e7]

        for P in pressures:
            state = cb.thermo_state(T=300, P=P, X=X)

            # Density should scale linearly with pressure (ideal gas)
            assert np.isfinite(state.rho)
            assert state.rho > 0

            # Other properties should be independent of pressure
            assert np.isfinite(state.cp)
            assert np.isfinite(state.h)


class TestConsistencyAcrossConditions:
    """Test consistency of relationships across various conditions."""

    def test_gamma_always_greater_than_one(self):
        """Verify gamma > 1 for all conditions."""
        X = cb.standard_dry_air_composition()
        conditions = [
            (200, 1e5),
            (300, 1e5),
            (1000, 1e5),
            (2000, 1e5),
            (300, 1e6),
        ]

        for T, P in conditions:
            state = cb.thermo_state(T=T, P=P, X=X)
            assert state.gamma > 1.0

    def test_speed_of_sound_increases_with_temperature(self):
        """Verify speed of sound increases with temperature."""
        X = cb.standard_dry_air_composition()
        state1 = cb.thermo_state(T=300, P=101325, X=X)
        state2 = cb.thermo_state(T=600, P=101325, X=X)

        assert state2.a > state1.a

    def test_density_inversely_proportional_to_temperature(self):
        """Verify density ~ 1/T at constant pressure."""
        X = cb.standard_dry_air_composition()
        P = 101325

        state1 = cb.thermo_state(T=300, P=P, X=X)
        state2 = cb.thermo_state(T=600, P=P, X=X)

        # rho1 / rho2 should be approximately T2 / T1
        ratio = state1.rho / state2.rho
        expected_ratio = 600 / 300

        assert ratio == pytest.approx(expected_ratio, rel=1e-10)
