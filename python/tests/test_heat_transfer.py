"""Tests for heat transfer correlations.

Verifies Python bindings for Nusselt number and LMTD functions.
Units: Nu dimensionless [-], h in [W/(m²·K)], LMTD in [K].
"""

import pytest

import combaero as cb


class TestNusseltCorrelations:
    """Test Nusselt number correlation bindings."""

    def test_nusselt_dittus_boelter_heating(self):
        """Test Dittus-Boelter correlation for heating."""
        Re = 1e5  # Reynolds number [-]
        Pr = 0.7  # Prandtl number (air) [-]

        Nu = cb.nusselt_dittus_boelter(Re, Pr, heating=True)

        # Nusselt number should be dimensionless and positive
        assert Nu > 0
        assert isinstance(Nu, float)

        # For Re=1e5, Pr=0.7, Nu should be O(200-300)
        assert 150 < Nu < 400

    def test_nusselt_dittus_boelter_cooling(self):
        """Test Dittus-Boelter correlation for cooling."""
        Re = 1e5  # Reynolds number [-]
        Pr = 0.7  # Prandtl number [-]

        Nu_heating = cb.nusselt_dittus_boelter(Re, Pr, heating=True)
        Nu_cooling = cb.nusselt_dittus_boelter(Re, Pr, heating=False)

        # For Pr < 1 (like air), cooling (n=0.3) actually gives higher Nu than heating (n=0.4)
        # because Pr^0.3 > Pr^0.4 when Pr < 1
        # Both should be positive and in reasonable range
        assert Nu_cooling > 0
        assert Nu_heating > 0
        assert 150 < Nu_heating < 250
        assert 150 < Nu_cooling < 250

        # They should be close but different
        assert abs(Nu_cooling - Nu_heating) / Nu_heating < 0.1  # Within 10%

    def test_nusselt_gnielinski_with_friction(self):
        """Test Gnielinski correlation with explicit friction factor."""
        Re = 1e5  # Reynolds number [-]
        Pr = 0.7  # Prandtl number [-]
        f = cb.friction_colebrook(Re, 0.0)  # Smooth pipe friction [-]

        Nu = cb.nusselt_gnielinski(Re, Pr, f)

        assert Nu > 0
        assert isinstance(Nu, float)

        # Gnielinski should give similar results to Dittus-Boelter
        Nu_db = cb.nusselt_dittus_boelter(Re, Pr)
        assert abs(Nu - Nu_db) / Nu_db < 0.2  # Within 20%

    def test_nusselt_gnielinski_auto_friction(self):
        """Test Gnielinski correlation with automatic friction."""
        Re = 1e5  # Reynolds number [-]
        Pr = 0.7  # Prandtl number [-]

        Nu_auto = cb.nusselt_gnielinski(Re, Pr)

        # Should match manual friction for smooth pipe
        f = cb.friction_petukhov(Re)  # Smooth pipe friction [-]
        Nu_manual = cb.nusselt_gnielinski(Re, Pr, f)

        # Should be very close
        assert abs(Nu_auto - Nu_manual) / Nu_manual < 0.01  # Within 1%

    def test_nusselt_sieder_tate_no_correction(self):
        """Test Sieder-Tate with no viscosity correction (mu_ratio=1)."""
        Re = 1e5  # Reynolds number [-]
        Pr = 0.7  # Prandtl number [-]

        Nu = cb.nusselt_sieder_tate(Re, Pr, mu_ratio=1.0)

        assert Nu > 0
        assert isinstance(Nu, float)

        # Sieder-Tate uses different constant (0.027) and exponent (Pr^1/3)
        # vs Dittus-Boelter (0.023, Pr^0.4), so they will differ
        # Just verify it's in a reasonable range
        assert 150 < Nu < 300

    def test_nusselt_sieder_tate_viscosity_correction(self):
        """Test Sieder-Tate with viscosity correction."""
        Re = 1e5  # Reynolds number [-]
        Pr = 0.7  # Prandtl number [-]

        # Heating: wall cooler than bulk, mu_wall > mu_bulk
        mu_ratio_heating = 0.8  # μ_bulk / μ_wall [-]
        Nu_heating = cb.nusselt_sieder_tate(Re, Pr, mu_ratio_heating)

        # Cooling: wall hotter than bulk, mu_wall < mu_bulk
        mu_ratio_cooling = 1.2  # μ_bulk / μ_wall [-]
        Nu_cooling = cb.nusselt_sieder_tate(Re, Pr, mu_ratio_cooling)

        # No correction
        Nu_no_corr = cb.nusselt_sieder_tate(Re, Pr, 1.0)

        # Heating should have lower Nu, cooling should have higher Nu
        assert Nu_heating < Nu_no_corr < Nu_cooling

    def test_nusselt_reynolds_number_effect(self):
        """Test Nusselt number increases with Reynolds number."""
        Pr = 0.7  # Fixed Prandtl number [-]

        Re_low = 1e4  # Reynolds number [-]
        Re_high = 1e6  # Reynolds number [-]

        Nu_low = cb.nusselt_dittus_boelter(Re_low, Pr)
        Nu_high = cb.nusselt_dittus_boelter(Re_high, Pr)

        # Higher Re should have higher Nu
        assert Nu_high > Nu_low

    def test_nusselt_prandtl_number_effect(self):
        """Test Nusselt number increases with Prandtl number."""
        Re = 1e5  # Fixed Reynolds number [-]

        Pr_low = 0.7  # Air [-]
        Pr_high = 7.0  # Water [-]

        Nu_low = cb.nusselt_dittus_boelter(Re, Pr_low)
        Nu_high = cb.nusselt_dittus_boelter(Re, Pr_high)

        # Higher Pr should have higher Nu
        assert Nu_high > Nu_low


class TestHeatTransferCoefficient:
    """Test heat transfer coefficient calculation."""

    def test_htc_from_nusselt_basic(self):
        """Test HTC calculation from Nusselt number."""
        Nu = 100.0  # Nusselt number [-]
        k = 0.025  # Thermal conductivity (air) [W/(m·K)]
        L = 0.05  # Characteristic length (pipe diameter) [m]

        h = cb.htc_from_nusselt(Nu, k, L)

        # h = Nu * k / L = 100 * 0.025 / 0.05 = 50 W/(m²·K)
        assert abs(h - 50.0) < 0.01

        # Units check: h should be in W/(m²·K)
        assert h > 0
        assert isinstance(h, float)

    def test_htc_units_verification(self):
        """Verify HTC has correct units [W/(m²·K)]."""
        Nu = 200.0  # Nusselt number [-]
        k = 0.6  # Thermal conductivity (water) [W/(m·K)]
        L = 0.1  # Diameter [m]

        h = cb.htc_from_nusselt(Nu, k, L)

        # h = 200 * 0.6 / 0.1 = 1200 W/(m²·K)
        assert abs(h - 1200.0) < 0.1

        # Typical HTC for water: 500-10,000 W/(m²·K)
        assert 500 < h < 10000

    def test_htc_diameter_effect(self):
        """Test HTC decreases with increasing diameter."""
        Nu = 100.0  # Fixed Nusselt number [-]
        k = 0.025  # Fixed thermal conductivity [W/(m·K)]

        D_small = 0.01  # Small diameter [m]
        D_large = 0.1  # Large diameter [m]

        h_small = cb.htc_from_nusselt(Nu, k, D_small)
        h_large = cb.htc_from_nusselt(Nu, k, D_large)

        # Smaller diameter gives higher HTC
        assert h_small > h_large

        # Ratio should match diameter ratio
        assert abs(h_small / h_large - D_large / D_small) < 0.01


class TestLMTD:
    """Test log mean temperature difference."""

    def test_lmtd_counterflow(self):
        """Test LMTD for counter-flow heat exchanger."""
        # Counter-flow: hot and cold streams flow in opposite directions
        # Hot: 100°C -> 60°C
        # Cold: 20°C -> 50°C
        dT1 = 100 - 50  # Hot in - Cold out = 50 K
        dT2 = 60 - 20  # Hot out - Cold in = 40 K

        lmtd_val = cb.lmtd(dT1, dT2)

        # LMTD = (50 - 40) / ln(50/40) = 10 / 0.223 ≈ 44.8 K
        assert abs(lmtd_val - 44.8) < 0.5

        # LMTD should be between dT1 and dT2
        assert min(dT1, dT2) < lmtd_val < max(dT1, dT2)

    def test_lmtd_parallelflow(self):
        """Test LMTD for parallel-flow heat exchanger."""
        # Parallel-flow: both streams flow in same direction
        # Hot: 100°C -> 60°C
        # Cold: 20°C -> 40°C
        dT1 = 100 - 20  # Hot in - Cold in = 80 K
        dT2 = 60 - 40  # Hot out - Cold out = 20 K

        lmtd_val = cb.lmtd(dT1, dT2)

        # LMTD = (80 - 20) / ln(80/20) = 60 / 1.386 ≈ 43.3 K
        assert abs(lmtd_val - 43.3) < 0.5

        # LMTD should be between dT1 and dT2
        assert min(dT1, dT2) < lmtd_val < max(dT1, dT2)

    def test_lmtd_equal_temperatures(self):
        """Test LMTD when temperature differences are equal."""
        dT1 = 50.0  # Temperature difference [K]
        dT2 = 50.0  # Same temperature difference [K]

        lmtd_val = cb.lmtd(dT1, dT2)

        # When dT1 = dT2, LMTD = dT1 = dT2 (arithmetic mean)
        assert abs(lmtd_val - 50.0) < 0.01

    def test_lmtd_nearly_equal_temperatures(self):
        """Test LMTD when temperature differences are nearly equal."""
        dT1 = 50.0  # Temperature difference [K]
        dT2 = 50.01  # Nearly same [K]

        lmtd_val = cb.lmtd(dT1, dT2)

        # Should return arithmetic mean to avoid 0/0
        assert abs(lmtd_val - 50.005) < 0.01

    def test_lmtd_units_kelvin(self):
        """Verify LMTD returns temperature in Kelvin."""
        dT1 = 100.0  # Temperature difference [K]
        dT2 = 50.0  # Temperature difference [K]

        lmtd_val = cb.lmtd(dT1, dT2)

        # LMTD should be in Kelvin (same units as input)
        # LMTD = (100 - 50) / ln(100/50) = 50 / 0.693 ≈ 72.1 K
        assert abs(lmtd_val - 72.1) < 0.5

        # Should be positive and between the two inputs
        assert 50 < lmtd_val < 100


class TestIntegration:
    """Integration tests combining multiple functions."""

    def test_pipe_flow_heat_transfer(self):
        """Test complete pipe flow heat transfer calculation."""
        # Air flow in a pipe
        Re = 5e4  # Reynolds number [-]
        Pr = 0.7  # Prandtl number (air) [-]
        k = 0.026  # Thermal conductivity (air at 300K) [W/(m·K)]
        D = 0.05  # Pipe diameter [m]

        # Calculate Nusselt number
        Nu = cb.nusselt_dittus_boelter(Re, Pr)

        # Calculate heat transfer coefficient
        h = cb.htc_from_nusselt(Nu, k, D)

        # Typical HTC for air: 10-100 W/(m²·K)
        assert 10 < h < 200

        # Verify units are consistent
        assert h > 0
        assert isinstance(h, float)

    def test_heat_exchanger_design(self):
        """Test heat exchanger design calculation."""
        # Counter-flow heat exchanger
        T_hot_in = 373.15  # 100°C [K]
        T_hot_out = 333.15  # 60°C [K]
        T_cold_in = 293.15  # 20°C [K]
        T_cold_out = 323.15  # 50°C [K]

        # Calculate LMTD
        dT1 = T_hot_in - T_cold_out  # 50 K
        dT2 = T_hot_out - T_cold_in  # 40 K
        lmtd_val = cb.lmtd(dT1, dT2)

        # LMTD should be reasonable
        assert 40 < lmtd_val < 50

        # Heat transfer rate: Q = U * A * LMTD
        # For given U and A, higher LMTD means higher Q
        assert lmtd_val > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
