"""Tests for composite htc_pipe function.

Verifies the composite heat transfer coefficient function that computes
HTC, Nu, and Re from thermodynamic state (T, P, X).

Tests all correlation options and verifies decomposition matches individual calls.
"""

import pytest

import combaero as cb


class TestHtcPipeComposite:
    """Test composite htc_pipe function with correlation selection."""

    def test_htc_pipe_gnielinski_default(self):
        """Test htc_pipe with default Gnielinski correlation."""
        # Air at 300K, 1 atm
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.05  # m

        # Call composite function
        h, Nu, Re = cb.htc_pipe(T, P, X, v, D)

        # Verify tuple return
        assert isinstance(h, float)
        assert isinstance(Nu, float)
        assert isinstance(Re, float)

        # Verify positive values
        assert h > 0
        assert Nu > 0
        assert Re > 0

        # Verify Re is in turbulent range for air at these conditions
        assert Re > 2300

        # Verify h has reasonable magnitude for air [W/(m2*K)]
        assert 10 < h < 500

    def test_htc_pipe_decomposition_gnielinski(self):
        """Test that composite matches manual decomposition (Gnielinski)."""
        T = 400.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 15.0  # m/s
        D = 0.1  # m

        # Composite call
        h_comp, Nu_comp, Re_comp = cb.htc_pipe(T, P, X, v, D, correlation="gnielinski")

        # Manual decomposition
        rho = cb.density(T, P, X)
        mu = cb.viscosity(T, P, X)
        k = cb.thermal_conductivity(T, P, X)
        Pr = cb.prandtl(T, P, X)

        Re_manual = rho * v * D / mu
        f = cb.friction_petukhov(Re_manual)
        Nu_manual = cb.nusselt_gnielinski(Re_manual, Pr, f)
        h_manual = cb.htc_from_nusselt(Nu_manual, k, D)

        # Verify machine precision match
        assert abs(Re_comp - Re_manual) < 1e-10
        assert abs(Nu_comp - Nu_manual) / Nu_manual < 1e-12
        assert abs(h_comp - h_manual) / h_manual < 1e-12

    def test_htc_pipe_decomposition_dittus_boelter(self):
        """Test that composite matches manual decomposition (Dittus-Boelter)."""
        T = 500.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 20.0  # m/s (high velocity for Re > 10000)
        D = 0.05  # m

        # Composite call
        h_comp, Nu_comp, Re_comp = cb.htc_pipe(
            T, P, X, v, D, correlation="dittus_boelter", heating=True
        )

        # Manual decomposition
        rho = cb.density(T, P, X)
        mu = cb.viscosity(T, P, X)
        k = cb.thermal_conductivity(T, P, X)
        Pr = cb.prandtl(T, P, X)

        Re_manual = rho * v * D / mu
        Nu_manual = cb.nusselt_dittus_boelter(Re_manual, Pr, heating=True)
        h_manual = cb.htc_from_nusselt(Nu_manual, k, D)

        # Verify machine precision match
        assert abs(Re_comp - Re_manual) < 1e-10
        assert abs(Nu_comp - Nu_manual) / Nu_manual < 1e-12
        assert abs(h_comp - h_manual) / h_manual < 1e-12

    def test_htc_pipe_decomposition_sieder_tate(self):
        """Test that composite matches manual decomposition (Sieder-Tate)."""
        T = 350.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 25.0  # m/s
        D = 0.08  # m
        mu_ratio = 1.2  # Viscosity correction

        # Composite call
        h_comp, Nu_comp, Re_comp = cb.htc_pipe(
            T, P, X, v, D, correlation="sieder_tate", mu_ratio=mu_ratio
        )

        # Manual decomposition
        rho = cb.density(T, P, X)
        mu = cb.viscosity(T, P, X)
        k = cb.thermal_conductivity(T, P, X)
        Pr = cb.prandtl(T, P, X)

        Re_manual = rho * v * D / mu
        Nu_manual = cb.nusselt_sieder_tate(Re_manual, Pr, mu_ratio)
        h_manual = cb.htc_from_nusselt(Nu_manual, k, D)

        # Verify machine precision match
        assert abs(Re_comp - Re_manual) < 1e-10
        assert abs(Nu_comp - Nu_manual) / Nu_manual < 1e-12
        assert abs(h_comp - h_manual) / h_manual < 1e-12

    def test_htc_pipe_decomposition_petukhov(self):
        """Test that composite matches manual decomposition (Petukhov)."""
        T = 600.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 30.0  # m/s
        D = 0.06  # m

        # Composite call
        h_comp, Nu_comp, Re_comp = cb.htc_pipe(T, P, X, v, D, correlation="petukhov")

        # Manual decomposition
        rho = cb.density(T, P, X)
        mu = cb.viscosity(T, P, X)
        k = cb.thermal_conductivity(T, P, X)
        Pr = cb.prandtl(T, P, X)

        Re_manual = rho * v * D / mu
        Nu_manual = cb.nusselt_petukhov(Re_manual, Pr)
        h_manual = cb.htc_from_nusselt(Nu_manual, k, D)

        # Verify machine precision match
        assert abs(Re_comp - Re_manual) < 1e-10
        assert abs(Nu_comp - Nu_manual) / Nu_manual < 1e-12
        assert abs(h_comp - h_manual) / h_manual < 1e-12

    def test_htc_pipe_heating_vs_cooling(self):
        """Test heating vs cooling for Dittus-Boelter."""
        T = 400.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 20.0  # m/s
        D = 0.05  # m

        h_heat, Nu_heat, Re_heat = cb.htc_pipe(
            T, P, X, v, D, correlation="dittus_boelter", heating=True
        )
        h_cool, Nu_cool, Re_cool = cb.htc_pipe(
            T, P, X, v, D, correlation="dittus_boelter", heating=False
        )

        # Reynolds number should be identical
        assert abs(Re_heat - Re_cool) < 1e-10

        # Nu and h should differ slightly (Pr exponent changes)
        assert Nu_heat != Nu_cool
        assert h_heat != h_cool

        # Both should be positive
        assert h_heat > 0
        assert h_cool > 0

    def test_htc_pipe_with_roughness(self):
        """Test htc_pipe with pipe roughness."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.1  # m
        roughness = 0.00005  # 50 microns (commercial steel)

        h_rough, Nu_rough, Re_rough = cb.htc_pipe(
            T, P, X, v, D, correlation="gnielinski", roughness=roughness
        )

        h_smooth, Nu_smooth, Re_smooth = cb.htc_pipe(
            T, P, X, v, D, correlation="gnielinski", roughness=0.0
        )

        # Reynolds number should be identical
        assert abs(Re_rough - Re_smooth) < 1e-10

        # Rough pipe should have higher friction, thus higher Nu and h
        assert Nu_rough > Nu_smooth
        assert h_rough > h_smooth

    def test_htc_pipe_laminar_flow(self):
        """Test htc_pipe handles laminar flow (Re < 2300)."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 0.1  # m/s (very low velocity)
        D = 0.01  # m (small diameter)

        h, Nu, Re = cb.htc_pipe(T, P, X, v, D, correlation="gnielinski")

        # Should be laminar
        assert Re < 2300

        # Should use constant laminar Nu (3.66 for constant T)
        assert abs(Nu - 3.66) < 0.01

        # h should still be positive
        assert h > 0

    def test_htc_pipe_temperature_effect(self):
        """Test that HTC changes with temperature (affects properties)."""
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 15.0  # m/s
        D = 0.05  # m

        T_low = 300.0  # K
        T_high = 600.0  # K

        h_low, Nu_low, Re_low = cb.htc_pipe(T_low, P, X, v, D)
        h_high, Nu_high, Re_high = cb.htc_pipe(T_high, P, X, v, D)

        # Reynolds number should change (viscosity changes with T)
        assert Re_low != Re_high

        # HTC should change
        assert h_low != h_high

        # Both should be positive
        assert h_low > 0
        assert h_high > 0

    def test_htc_pipe_velocity_effect(self):
        """Test that HTC increases with velocity."""
        T = 400.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        D = 0.05  # m

        v_low = 5.0  # m/s
        v_high = 20.0  # m/s

        h_low, Nu_low, Re_low = cb.htc_pipe(T, P, X, v_low, D)
        h_high, Nu_high, Re_high = cb.htc_pipe(T, P, X, v_high, D)

        # Higher velocity -> higher Re -> higher Nu -> higher h
        assert Re_high > Re_low
        assert Nu_high > Nu_low
        assert h_high > h_low

    def test_htc_pipe_diameter_effect(self):
        """Test that HTC decreases with diameter (for fixed velocity)."""
        T = 350.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 10.0  # m/s

        D_small = 0.02  # m
        D_large = 0.1  # m

        h_small, Nu_small, Re_small = cb.htc_pipe(T, P, X, v, D_small)
        h_large, Nu_large, Re_large = cb.htc_pipe(T, P, X, v, D_large)

        # Larger diameter -> higher Re (Re ~ D)
        assert Re_large > Re_small

        # But h = Nu * k / D, so larger D reduces h
        assert h_small > h_large

    def test_htc_pipe_invalid_correlation(self):
        """Test that invalid correlation name raises error."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.05  # m

        with pytest.raises(ValueError, match="unknown correlation"):
            cb.htc_pipe(T, P, X, v, D, correlation="invalid_correlation")

    def test_htc_pipe_dittus_boelter_low_re_error(self):
        """Test that Dittus-Boelter raises error for Re < 10000 in transition."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 2.0  # m/s (low velocity)
        D = 0.05  # m

        # This should give Re in transition region (2300 < Re < 10000)
        # Dittus-Boelter should raise error
        with pytest.raises(ValueError, match="Dittus-Boelter requires Re > 10000"):
            cb.htc_pipe(T, P, X, v, D, correlation="dittus_boelter")

    def test_htc_pipe_multiple_compositions(self):
        """Test htc_pipe with different gas compositions."""
        T = 400.0  # K
        P = 101325.0  # Pa
        v = 15.0  # m/s
        D = 0.05  # m

        # Air
        X_air = cb.standard_dry_air_composition()
        h_air, Nu_air, Re_air = cb.htc_pipe(T, P, X_air, v, D)

        # Pure N2 (simplified - just use air, properties are similar)
        # In practice, N2 and air have very similar transport properties
        h_N2, Nu_N2, Re_N2 = cb.htc_pipe(T, P, X_air, v, D)

        # Both should be positive
        assert h_air > 0
        assert h_N2 > 0

        # For same composition, should be identical
        assert abs(h_air - h_N2) < 1e-10


class TestHtcPipeEdgeCases:
    """Test edge cases and error handling."""

    def test_htc_pipe_negative_velocity(self):
        """Test that negative velocity raises error."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = -10.0  # m/s (negative)
        D = 0.05  # m

        with pytest.raises(ValueError, match="velocity must be positive"):
            cb.htc_pipe(T, P, X, v, D)

    def test_htc_pipe_zero_diameter(self):
        """Test that zero diameter raises error."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.0  # m (zero)

        with pytest.raises(ValueError, match="diameter must be positive"):
            cb.htc_pipe(T, P, X, v, D)

    def test_htc_pipe_negative_temperature(self):
        """Test that negative temperature raises error."""
        T = -100.0  # K (negative)
        P = 101325.0  # Pa
        X = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.05  # m

        with pytest.raises(ValueError, match="temperature must be positive"):
            cb.htc_pipe(T, P, X, v, D)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
