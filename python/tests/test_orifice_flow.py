"""Tests for orifice flow utilities with real gas support.

Verifies the OrificeFlowResult dataclass and orifice flow utility functions.
Tests real gas corrections via compressibility factor Z.
"""

import numpy as np
import pytest

import combaero as cb


class TestOrificeFlowResult:
    """Test OrificeFlowResult struct/class behavior."""

    def test_orifice_flow_returns_struct(self):
        """Test that orifice_flow returns OrificeFlowResult object."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05  # 50 mm
        geom.D = 0.1  # 100 mm

        result = cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5)

        assert isinstance(result, cb.OrificeFlowResult)

    def test_orifice_flow_has_all_attributes(self):
        """Test that OrificeFlowResult has all expected attributes."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        result = cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5)

        # Check all 7 properties exist
        assert hasattr(result, "mdot")
        assert hasattr(result, "v")
        assert hasattr(result, "Re_D")
        assert hasattr(result, "Re_d")
        assert hasattr(result, "Cd")
        assert hasattr(result, "epsilon")
        assert hasattr(result, "rho_corrected")

    def test_orifice_flow_readonly(self):
        """Test that OrificeFlowResult attributes are read-only."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        result = cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5)

        # Should not be able to modify properties
        with pytest.raises(AttributeError):
            result.mdot = 999.0

    def test_orifice_flow_repr(self):
        """Test that OrificeFlowResult has nice repr."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        result = cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5)
        repr_str = repr(result)

        # Should contain key properties
        assert "OrificeFlowResult" in repr_str
        assert "mdot=" in repr_str
        assert "kg/s" in repr_str


class TestOrificeFlowIdealGas:
    """Test orifice_flow with ideal gas (Z=1.0)."""

    def test_ideal_gas_default_Z(self):
        """Test that default Z=1.0 gives ideal gas behavior."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        result = cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5)

        # Should have positive values
        assert result.mdot > 0
        assert result.v > 0
        assert result.Re_D > 0
        assert result.Re_d > 0
        assert 0.5 < result.Cd < 0.7  # Typical range
        assert result.epsilon == 1.0  # Incompressible
        assert result.rho_corrected > 0

    def test_velocity_matches_formula(self):
        """Test that velocity matches v = mdot / (rho * A)."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        result = cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5, Z=1.0)

        # Calculate expected velocity
        A = np.pi * geom.d**2 / 4.0
        v_expected = result.mdot / (result.rho_corrected * A)

        # Must match within machine precision
        assert abs(result.v - v_expected) < 1e-12

    def test_reynolds_numbers_match_formula(self):
        """Test that Reynolds numbers match formulas."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1
        mu = 1.8e-5

        result = cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=mu, Z=1.0)

        # Re_D = 4 * mdot / (pi * D * μ)
        Re_D_expected = (4.0 * result.mdot) / (np.pi * geom.D * mu)
        assert abs(result.Re_D - Re_D_expected) < 1e-10

        # Re_d = 4 * mdot / (pi * d * μ)
        Re_d_expected = (4.0 * result.mdot) / (np.pi * geom.d * mu)
        assert abs(result.Re_d - Re_d_expected) < 1e-10


class TestOrificeFlowRealGas:
    """Test orifice_flow with real gas (Z != 1.0)."""

    def test_real_gas_density_correction(self):
        """Test that Z factor correctly adjusts density."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1
        T = 300
        P = 10e6  # 100 bar
        mu = 1.1e-5
        Z = 0.85  # Real gas

        result = cb.orifice_flow(geom, dP=50000, T=T, P=P, mu=mu, Z=Z)

        # Calculate ideal gas density
        R_gas = 8.314  # J/(mol*K)
        MW_air = 0.02897  # kg/mol
        rho_ideal = (P * MW_air) / (R_gas * T)

        # Real gas correction
        rho_expected = rho_ideal / Z

        # Must match within machine precision
        assert abs(result.rho_corrected - rho_expected) < 1e-12

    def test_real_gas_higher_mdot_than_ideal(self):
        """Test that real gas (Z<1) gives higher mass flow than ideal gas."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1
        T = 300
        P = 10e6  # 100 bar
        mu = 1.1e-5

        # Ideal gas
        result_ideal = cb.orifice_flow(geom, dP=50000, T=T, P=P, mu=mu, Z=1.0)

        # Real gas (higher density)
        result_real = cb.orifice_flow(geom, dP=50000, T=T, P=P, mu=mu, Z=0.85)

        # Real gas has higher density, so higher mass flow
        assert result_real.rho_corrected > result_ideal.rho_corrected
        assert result_real.mdot > result_ideal.mdot

    def test_Z_range_typical_values(self):
        """Test with typical Z values for real gases."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        Z_values = [0.7, 0.8, 0.9, 0.95, 1.0]
        mdot_values = []

        for Z in Z_values:
            result = cb.orifice_flow(geom, dP=50000, T=300, P=10e6, mu=1.1e-5, Z=Z)
            mdot_values.append(result.mdot)
            assert result.mdot > 0
            assert result.rho_corrected > 0

        # Mass flow should increase as Z decreases (higher density)
        for i in range(len(mdot_values) - 1):
            assert mdot_values[i] > mdot_values[i + 1]


class TestUtilityFunctions:
    """Test utility functions."""

    def test_orifice_velocity_from_mdot_ideal_gas(self):
        """Test velocity calculation with ideal gas."""
        mdot = 1.0  # kg/s
        rho = 1.2  # kg/m^3
        d = 0.05  # m
        Z = 1.0

        v = cb.orifice_velocity_from_mdot(mdot, rho, d, Z)

        # v = mdot / (rho * A)
        A = np.pi * d**2 / 4.0
        v_expected = mdot / (rho * A)

        assert abs(v - v_expected) < 1e-15

    def test_orifice_velocity_from_mdot_real_gas(self):
        """Test velocity calculation with real gas correction."""
        mdot = 1.0  # kg/s
        rho = 1.2  # kg/m^3
        d = 0.05  # m
        Z = 0.85

        v = cb.orifice_velocity_from_mdot(mdot, rho, d, Z)

        # v = mdot / (rho_corrected * A) where rho_corrected = rho / Z
        rho_corrected = rho / Z
        A = np.pi * d**2 / 4.0
        v_expected = mdot / (rho_corrected * A)

        assert abs(v - v_expected) < 1e-15

    def test_orifice_area_from_beta(self):
        """Test orifice area calculation from beta."""
        D = 0.1  # m
        beta = 0.5

        A = cb.orifice_area_from_beta(D, beta)

        # A = pi * (D * beta / 2)^2
        d = D * beta
        A_expected = np.pi * d**2 / 4.0

        assert abs(A - A_expected) < 1e-15

    def test_beta_from_diameters(self):
        """Test beta calculation from diameters."""
        d = 0.05  # m
        D = 0.1  # m

        beta = cb.beta_from_diameters(d, D)

        # beta = d / D
        beta_expected = d / D

        assert abs(beta - beta_expected) < 1e-15

    def test_orifice_Re_d_from_mdot(self):
        """Test Reynolds number calculation from mass flow."""
        mdot = 1.0  # kg/s
        d = 0.05  # m
        mu = 1.8e-5  # Pa*s

        Re_d = cb.orifice_Re_d_from_mdot(mdot, d, mu)

        # Re_d = 4 * mdot / (pi * d * μ)
        Re_d_expected = (4.0 * mdot) / (np.pi * d * mu)

        assert abs(Re_d - Re_d_expected) < 1e-10

    def test_round_trip_beta_area(self):
        """Test round-trip: D, beta → A → beta."""
        D = 0.1
        beta_original = 0.5

        # beta → area
        A = cb.orifice_area_from_beta(D, beta_original)

        # area → d → beta
        d = np.sqrt(4.0 * A / np.pi)
        beta_recovered = cb.beta_from_diameters(d, D)

        assert abs(beta_recovered - beta_original) < 1e-15


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_negative_temperature_raises_error(self):
        """Test that negative temperature raises error."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        with pytest.raises(ValueError, match="temperature must be positive"):
            cb.orifice_flow(geom, dP=10000, T=-100, P=101325, mu=1.8e-5)

    def test_negative_pressure_raises_error(self):
        """Test that negative pressure raises error."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        with pytest.raises(ValueError, match="pressure must be positive"):
            cb.orifice_flow(geom, dP=10000, T=300, P=-101325, mu=1.8e-5)

    def test_negative_Z_raises_error(self):
        """Test that negative Z raises error."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        with pytest.raises(ValueError, match="compressibility factor Z must be positive"):
            cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5, Z=-0.5)

    def test_zero_Z_raises_error(self):
        """Test that zero Z raises error."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        with pytest.raises(ValueError, match="compressibility factor Z must be positive"):
            cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5, Z=0.0)

    def test_invalid_geometry_raises_error(self):
        """Test that invalid geometry raises error."""
        geom = cb.OrificeGeometry()
        geom.d = 0.15  # Larger than D
        geom.D = 0.1

        with pytest.raises(ValueError, match="invalid geometry"):
            cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5)

    def test_zero_dP_gives_zero_mdot(self):
        """Test that zero pressure drop gives zero mass flow."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        result = cb.orifice_flow(geom, dP=0.0, T=300, P=101325, mu=1.8e-5)

        assert result.mdot == 0.0


class TestBetaEffects:
    """Test effects of beta ratio on flow."""

    def test_beta_range(self):
        """Test with various beta ratios."""
        D = 0.1
        beta_values = [0.2, 0.3, 0.5, 0.7]

        for beta in beta_values:
            geom = cb.OrificeGeometry()
            geom.d = D * beta
            geom.D = D

            result = cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5)

            assert result.mdot > 0
            assert 0.5 < result.Cd < 0.7
            assert geom.beta == pytest.approx(beta, rel=1e-10)

    def test_larger_beta_higher_mdot(self):
        """Test that larger beta (larger orifice) gives higher mass flow."""
        D = 0.1
        beta_values = [0.3, 0.5, 0.7]
        mdot_values = []

        for beta in beta_values:
            geom = cb.OrificeGeometry()
            geom.d = D * beta
            geom.D = D

            result = cb.orifice_flow(geom, dP=10000, T=300, P=101325, mu=1.8e-5)
            mdot_values.append(result.mdot)

        # Mass flow should increase with beta
        for i in range(len(mdot_values) - 1):
            assert mdot_values[i + 1] > mdot_values[i]


class TestPressureEffects:
    """Test effects of pressure on flow."""

    def test_pressure_range(self):
        """Test with various pressures."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        pressures = [1e5, 5e5, 1e6, 5e6, 1e7]  # 1 to 100 bar

        for P in pressures:
            result = cb.orifice_flow(geom, dP=10000, T=300, P=P, mu=1.8e-5)

            assert result.mdot > 0
            assert result.rho_corrected > 0

    def test_higher_pressure_higher_density(self):
        """Test that higher pressure gives higher density."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        result_low = cb.orifice_flow(geom, dP=10000, T=300, P=1e5, mu=1.8e-5)
        result_high = cb.orifice_flow(geom, dP=10000, T=300, P=1e7, mu=1.8e-5)

        # Higher pressure → higher density → higher mass flow
        assert result_high.rho_corrected > result_low.rho_corrected
        assert result_high.mdot > result_low.mdot


class TestTypicalApplications:
    """Test with typical engineering applications."""

    def test_natural_gas_metering_low_pressure(self):
        """Test natural gas metering at low pressure (ideal gas)."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05  # 50 mm
        geom.D = 0.1  # 100 mm (beta = 0.5)

        # Low pressure natural gas (approximately ideal)
        result = cb.orifice_flow(geom, dP=5000, T=288, P=2e5, mu=1.1e-5, Z=0.998)

        # Should have reasonable values
        assert result.mdot > 0
        assert result.v > 0
        assert result.Cd > 0.5
        assert result.Re_D > 1000

    def test_natural_gas_metering_high_pressure(self):
        """Test natural gas metering at high pressure (real gas)."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        # High pressure natural gas (significant real gas effects)
        result = cb.orifice_flow(geom, dP=50000, T=300, P=10e6, mu=1.1e-5, Z=0.85)

        # Should have reasonable values
        assert result.mdot > 0
        assert result.v > 0
        assert result.rho_corrected > 10  # High density at 100 bar
        assert result.Cd > 0.5

    def test_air_flow_metering(self):
        """Test air flow metering at standard conditions."""
        geom = cb.OrificeGeometry()
        geom.d = 0.025  # 25 mm
        geom.D = 0.05  # 50 mm (beta = 0.5)

        # Standard air
        result = cb.orifice_flow(geom, dP=2000, T=288, P=101325, mu=1.8e-5, Z=1.0)

        # Should have reasonable values
        assert 0.01 < result.mdot < 1.0  # Reasonable range
        assert 1.0 < result.rho_corrected < 1.5  # Air density
        assert result.Cd > 0.5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
