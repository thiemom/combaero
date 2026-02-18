"""Tests for air_properties bundle function.

Verifies the AirProperties dataclass-like struct and air_properties() function.
Tests that all properties match individual function calls with machine precision.
"""

import pytest

import combaero as cb


class TestAirPropertiesStruct:
    """Test AirProperties struct/class behavior."""

    def test_air_properties_returns_struct(self):
        """Test that air_properties returns AirProperties object."""
        T = 300.0  # K
        P = 101325.0  # Pa

        props = cb.air_properties(T, P)

        # Should be AirProperties instance
        assert isinstance(props, cb.AirProperties)

    def test_air_properties_has_all_attributes(self):
        """Test that AirProperties has all expected attributes."""
        T = 300.0  # K
        P = 101325.0  # Pa

        props = cb.air_properties(T, P)

        # Check all 10 properties exist
        assert hasattr(props, "rho")
        assert hasattr(props, "mu")
        assert hasattr(props, "k")
        assert hasattr(props, "cp")
        assert hasattr(props, "cv")
        assert hasattr(props, "Pr")
        assert hasattr(props, "nu")
        assert hasattr(props, "alpha")
        assert hasattr(props, "gamma")
        assert hasattr(props, "a")

    def test_air_properties_readonly(self):
        """Test that AirProperties attributes are read-only."""
        T = 300.0  # K
        P = 101325.0  # Pa

        props = cb.air_properties(T, P)

        # Should not be able to modify properties
        with pytest.raises(AttributeError):
            props.rho = 999.0

    def test_air_properties_repr(self):
        """Test that AirProperties has nice repr."""
        T = 300.0  # K
        P = 101325.0  # Pa

        props = cb.air_properties(T, P)
        repr_str = repr(props)

        # Should contain key properties (AirProperties is now an alias for TransportState)
        assert "AirProperties" in repr_str or "TransportState" in repr_str
        assert "mu=" in repr_str


class TestAirPropertiesDryAir:
    """Test air_properties with dry air (humidity=0)."""

    def test_dry_air_density_matches(self):
        """Test that density matches individual function call."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()

        props = cb.air_properties(T, P, humidity=0.0)
        rho_expected = cb.density(T, P, X)

        # Must match within machine precision
        assert abs(props.rho - rho_expected) < 1e-12

    def test_dry_air_viscosity_matches(self):
        """Test that viscosity matches individual function call."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()

        props = cb.air_properties(T, P, humidity=0.0)
        mu_expected = cb.viscosity(T, P, X)

        # Must match within machine precision
        assert abs(props.mu - mu_expected) < 1e-15

    def test_dry_air_thermal_conductivity_matches(self):
        """Test that thermal conductivity matches individual function call."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()

        props = cb.air_properties(T, P, humidity=0.0)
        k_expected = cb.thermal_conductivity(T, P, X)

        # Must match within machine precision
        assert abs(props.k - k_expected) < 1e-15

    def test_dry_air_cp_matches(self):
        """Test that cp matches individual function call."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()

        props = cb.air_properties(T, P, humidity=0.0)
        cp_expected = cb.cp_mass(T, X)

        # Must match within machine precision
        assert abs(props.cp - cp_expected) < 1e-12

    def test_dry_air_cv_matches(self):
        """Test that cv matches individual function call."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()

        props = cb.air_properties(T, P, humidity=0.0)
        cv_expected = cb.cv_mass(T, X)

        # Must match within machine precision
        assert abs(props.cv - cv_expected) < 1e-12

    def test_dry_air_prandtl_matches(self):
        """Test that Prandtl number matches individual function call."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()

        props = cb.air_properties(T, P, humidity=0.0)
        Pr_expected = cb.prandtl(T, P, X)

        # Must match within machine precision
        assert abs(props.Pr - Pr_expected) < 1e-15

    def test_dry_air_speed_of_sound_matches(self):
        """Test that speed of sound matches individual function call."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()

        props = cb.air_properties(T, P, humidity=0.0)
        a_expected = cb.speed_of_sound(T, X)

        # Must match within machine precision
        assert abs(props.a - a_expected) < 1e-12

    def test_dry_air_derived_properties(self):
        """Test that derived properties are computed correctly."""
        T = 300.0  # K
        P = 101325.0  # Pa

        props = cb.air_properties(T, P, humidity=0.0)

        # nu = mu / rho
        nu_expected = props.mu / props.rho
        assert abs(props.nu - nu_expected) < 1e-15

        # alpha = k / (rho * cp)
        alpha_expected = props.k / (props.rho * props.cp)
        assert abs(props.alpha - alpha_expected) < 1e-15

        # gamma = cp / cv
        gamma_expected = props.cp / props.cv
        assert abs(props.gamma - gamma_expected) < 1e-15


class TestAirPropertiesHumidAir:
    """Test air_properties with humid air (humidity>0)."""

    def test_humid_air_density_matches(self):
        """Test that humid air density matches individual function call."""
        T = 300.0  # K
        P = 101325.0  # Pa
        humidity = 0.5  # 50% RH
        X = cb.humid_air_composition(T, P, humidity)

        props = cb.air_properties(T, P, humidity=humidity)
        rho_expected = cb.density(T, P, X)

        # Must match within machine precision
        assert abs(props.rho - rho_expected) < 1e-12

    def test_humid_air_all_properties_match(self):
        """Test that all humid air properties match individual calls."""
        T = 300.0  # K
        P = 101325.0  # Pa
        humidity = 0.7  # 70% RH
        X = cb.humid_air_composition(T, P, humidity)

        props = cb.air_properties(T, P, humidity=humidity)

        # All properties must match
        assert abs(props.rho - cb.density(T, P, X)) < 1e-12
        assert abs(props.mu - cb.viscosity(T, P, X)) < 1e-15
        assert abs(props.k - cb.thermal_conductivity(T, P, X)) < 1e-15
        assert abs(props.cp - cb.cp_mass(T, X)) < 1e-12
        assert abs(props.cv - cb.cv_mass(T, X)) < 1e-12
        assert abs(props.Pr - cb.prandtl(T, P, X)) < 1e-15
        assert abs(props.a - cb.speed_of_sound(T, X)) < 1e-12

    def test_humid_air_density_lower_than_dry(self):
        """Test that humid air is less dense than dry air."""
        T = 300.0  # K
        P = 101325.0  # Pa

        props_dry = cb.air_properties(T, P, humidity=0.0)
        props_humid = cb.air_properties(T, P, humidity=0.8)

        # Humid air should be less dense (water vapor is lighter)
        assert props_humid.rho < props_dry.rho


class TestAirPropertiesTemperatureEffects:
    """Test air_properties at different temperatures."""

    def test_temperature_range_200_to_600K(self):
        """Test air_properties across temperature range."""
        P = 101325.0  # Pa
        temperatures = [200.0, 300.0, 400.0, 500.0, 600.0]

        for T in temperatures:
            props = cb.air_properties(T, P)

            # All properties should be positive
            assert props.rho > 0
            assert props.mu > 0
            assert props.k > 0
            assert props.cp > 0
            assert props.cv > 0
            assert props.Pr > 0
            assert props.nu > 0
            assert props.alpha > 0
            assert props.gamma > 1.0  # Always > 1 for ideal gas
            assert props.a > 0

    def test_density_decreases_with_temperature(self):
        """Test that density decreases with increasing temperature."""
        P = 101325.0  # Pa

        props_cold = cb.air_properties(T=250.0, P=P)
        props_hot = cb.air_properties(T=400.0, P=P)

        # Density should decrease with temperature (ideal gas law)
        assert props_hot.rho < props_cold.rho

    def test_viscosity_increases_with_temperature(self):
        """Test that viscosity increases with temperature."""
        P = 101325.0  # Pa

        props_cold = cb.air_properties(T=250.0, P=P)
        props_hot = cb.air_properties(T=400.0, P=P)

        # Viscosity increases with temperature (Sutherland's law)
        assert props_hot.mu > props_cold.mu

    def test_speed_of_sound_increases_with_temperature(self):
        """Test that speed of sound increases with temperature."""
        P = 101325.0  # Pa

        props_cold = cb.air_properties(T=250.0, P=P)
        props_hot = cb.air_properties(T=400.0, P=P)

        # Speed of sound ~ sqrt(T)
        assert props_hot.a > props_cold.a


class TestAirPropertiesPressureEffects:
    """Test air_properties at different pressures."""

    def test_pressure_range_50kPa_to_200kPa(self):
        """Test air_properties across pressure range."""
        T = 300.0  # K
        pressures = [50000.0, 101325.0, 150000.0, 200000.0]

        for P in pressures:
            props = cb.air_properties(T, P)

            # All properties should be positive
            assert props.rho > 0
            assert props.mu > 0
            assert props.k > 0
            assert props.cp > 0
            assert props.Pr > 0

    def test_density_increases_with_pressure(self):
        """Test that density increases with pressure."""
        T = 300.0  # K

        props_low = cb.air_properties(T, P=50000.0)
        props_high = cb.air_properties(T, P=200000.0)

        # Density should increase with pressure (ideal gas law)
        assert props_high.rho > props_low.rho

    def test_kinematic_viscosity_decreases_with_pressure(self):
        """Test that kinematic viscosity decreases with pressure."""
        T = 300.0  # K

        props_low = cb.air_properties(T, P=50000.0)
        props_high = cb.air_properties(T, P=200000.0)

        # nu = mu / rho, so nu decreases as rho increases
        assert props_high.nu < props_low.nu


class TestAirPropertiesEdgeCases:
    """Test edge cases and error handling."""

    def test_negative_temperature_raises_error(self):
        """Test that negative temperature raises error."""
        with pytest.raises(ValueError, match="temperature must be positive"):
            cb.air_properties(T=-100.0, P=101325.0)

    def test_zero_temperature_raises_error(self):
        """Test that zero temperature raises error."""
        with pytest.raises(ValueError, match="temperature must be positive"):
            cb.air_properties(T=0.0, P=101325.0)

    def test_negative_pressure_raises_error(self):
        """Test that negative pressure raises error."""
        with pytest.raises(ValueError, match="pressure must be positive"):
            cb.air_properties(T=300.0, P=-101325.0)

    def test_humidity_below_zero_raises_error(self):
        """Test that humidity < 0 raises error."""
        with pytest.raises(ValueError, match="humidity must be in range"):
            cb.air_properties(T=300.0, P=101325.0, humidity=-0.1)

    def test_humidity_above_one_raises_error(self):
        """Test that humidity > 1 raises error."""
        with pytest.raises(ValueError, match="humidity must be in range"):
            cb.air_properties(T=300.0, P=101325.0, humidity=1.5)

    def test_humidity_exactly_zero(self):
        """Test that humidity=0 works (dry air)."""
        props = cb.air_properties(T=300.0, P=101325.0, humidity=0.0)
        assert props.rho > 0

    def test_humidity_exactly_one(self):
        """Test that humidity=1 works (saturated air)."""
        props = cb.air_properties(T=300.0, P=101325.0, humidity=1.0)
        assert props.rho > 0


class TestAirPropertiesTypical:
    """Test air_properties with typical engineering values."""

    def test_standard_conditions(self):
        """Test air properties at standard conditions (15C, 1 atm)."""
        T = 288.15  # K (15C)
        P = 101325.0  # Pa (1 atm)

        props = cb.air_properties(T, P)

        # Verify all properties match individual calls
        X = cb.standard_dry_air_composition()
        assert props.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)
        assert props.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert props.k == pytest.approx(cb.thermal_conductivity(T, P, X), rel=1e-15)
        assert props.cp == pytest.approx(cb.cp_mass(T, X), rel=1e-15)
        assert props.Pr == pytest.approx(cb.prandtl(T, P, X), rel=1e-15)
        assert props.a == pytest.approx(cb.speed_of_sound(T, X), rel=1e-15)

    def test_hot_day(self):
        """Test air properties on a hot day (35C, 1 atm, 60% RH)."""
        T = 308.15  # K (35C)
        P = 101325.0  # Pa
        humidity = 0.6  # 60% RH

        props = cb.air_properties(T, P, humidity=humidity)

        # Verify all properties match individual calls with humid air composition
        X = cb.humid_air_composition(T, P, humidity)
        assert props.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)
        assert props.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert props.Pr == pytest.approx(cb.prandtl(T, P, X), rel=1e-15)

    def test_cold_day(self):
        """Test air properties on a cold day (-10C, 1 atm)."""
        T = 263.15  # K (-10C)
        P = 101325.0  # Pa

        props = cb.air_properties(T, P)

        # Verify all properties match individual calls
        X = cb.standard_dry_air_composition()
        assert props.rho == pytest.approx(cb.density(T, P, X), rel=1e-15)
        assert props.mu == pytest.approx(cb.viscosity(T, P, X), rel=1e-15)
        assert props.a == pytest.approx(cb.speed_of_sound(T, X), rel=1e-15)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
