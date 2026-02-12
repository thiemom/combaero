"""Tests for composite pressure_drop_pipe function.

Verifies the composite function that combines:
1. Thermodynamic property evaluation (density, viscosity)
2. Reynolds number calculation
3. Friction factor correlation
4. Darcy-Weisbach pressure drop equation

All in a single convenient call.
"""

import numpy as np
import pytest

import combaero as cb


class TestPressureDropPipeComposite:
    """Test composite pressure_drop_pipe function (Phase 4)."""

    def test_basic_functionality(self):
        """Test basic pressure drop calculation with air."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.1  # m
        L = 100.0  # m

        dP, Re, f = cb.pressure_drop_pipe(T, P, X_air, v, D, L)

        # All outputs should be positive
        assert dP > 0
        assert Re > 0
        assert f > 0

        # Sanity checks for air at 10 m/s in 0.1m pipe
        assert 50000 < Re < 100000  # Turbulent flow
        assert 0.01 < f < 0.05  # Typical friction factor
        assert 100 < dP < 10000  # Reasonable pressure drop

    def test_matches_manual_calculation(self):
        """Verify composite function matches step-by-step manual calculation."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 15.0  # m/s
        D = 0.05  # m
        L = 50.0  # m
        roughness = 0.000045  # m (commercial steel)

        # Composite function
        dP_comp, Re_comp, f_comp = cb.pressure_drop_pipe(T, P, X_air, v, D, L, roughness, "haaland")

        # Manual calculation
        rho = cb.density(T, P, X_air)
        Re_manual = cb.reynolds(T, P, X_air, v, D)
        e_D = roughness / D
        f_manual = cb.friction_haaland(Re_manual, e_D)
        dP_manual = cb.pipe_dP(v, L, D, f_manual, rho)

        # Must match within machine precision
        assert abs(Re_comp - Re_manual) < 1e-10
        assert abs(f_comp - f_manual) < 1e-14
        assert abs(dP_comp - dP_manual) < 1e-10

    def test_all_correlations(self):
        """Test all friction factor correlations work."""
        T = 350.0  # K
        P = 200000.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 20.0  # m/s
        D = 0.08  # m
        L = 75.0  # m
        roughness = 0.00015  # m

        correlations = ["haaland", "serghides", "colebrook", "petukhov"]
        results = {}

        for corr in correlations:
            dP, Re, f = cb.pressure_drop_pipe(T, P, X_air, v, D, L, roughness, corr)
            results[corr] = (dP, Re, f)

            # All should give valid results
            assert dP > 0
            assert Re > 0
            assert f > 0

        # Reynolds number should be identical for all correlations
        Re_values = [r[1] for r in results.values()]
        for Re_val in Re_values:
            assert abs(Re_val - Re_values[0]) < 1e-10

        # Rough-pipe correlations should agree closely
        rough_correlations = ["haaland", "serghides", "colebrook"]
        f_rough = [results[c][2] for c in rough_correlations]
        f_rough_mean = np.mean(f_rough)
        for f_val in f_rough:
            # Within 1% (these correlations are very accurate)
            assert abs(f_val - f_rough_mean) / f_rough_mean < 0.01

        # Petukhov (smooth-pipe) should give lower f than rough-pipe correlations
        assert results["petukhov"][2] < f_rough_mean

    def test_smooth_vs_rough_pipe(self):
        """Test that rough pipes have higher friction factor."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 12.0  # m/s
        D = 0.1  # m
        L = 100.0  # m

        # Smooth pipe
        dP_smooth, Re_smooth, f_smooth = cb.pressure_drop_pipe(T, P, X_air, v, D, L, roughness=0.0)

        # Rough pipe (commercial steel)
        dP_rough, Re_rough, f_rough = cb.pressure_drop_pipe(
            T, P, X_air, v, D, L, roughness=0.000045
        )

        # Reynolds number should be identical
        assert abs(Re_smooth - Re_rough) < 1e-10

        # Rough pipe should have higher friction factor
        assert f_rough > f_smooth

        # Rough pipe should have higher pressure drop
        assert dP_rough > dP_smooth

    def test_velocity_dependence(self):
        """Test pressure drop scales with velocity squared."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()
        D = 0.1  # m
        L = 100.0  # m

        v1 = 5.0  # m/s
        v2 = 10.0  # m/s (double velocity)

        dP1, Re1, f1 = cb.pressure_drop_pipe(T, P, X_air, v1, D, L)
        dP2, Re2, f2 = cb.pressure_drop_pipe(T, P, X_air, v2, D, L)

        # Reynolds number should double
        assert abs(Re2 / Re1 - 2.0) < 0.01

        # For turbulent flow, f changes with Re (decreases as Re increases)
        # dP ~ f * v^2, so if f decreases and v doubles, dP increases by less than 4x
        # Typical range: 3.2x to 3.8x depending on Re regime
        assert 3.0 < dP2 / dP1 < 4.0

    def test_length_dependence(self):
        """Test pressure drop scales linearly with length."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.1  # m

        L1 = 50.0  # m
        L2 = 100.0  # m (double length)

        dP1, Re1, f1 = cb.pressure_drop_pipe(T, P, X_air, v, D, L1)
        dP2, Re2, f2 = cb.pressure_drop_pipe(T, P, X_air, v, D, L2)

        # Reynolds number should be identical
        assert abs(Re1 - Re2) < 1e-10

        # Friction factor should be identical
        assert abs(f1 - f2) < 1e-14

        # Pressure drop should double
        assert abs(dP2 / dP1 - 2.0) < 1e-12

    def test_temperature_effect(self):
        """Test effect of temperature on pressure drop (via viscosity)."""
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.1  # m
        L = 100.0  # m

        T_cold = 250.0  # K
        T_hot = 400.0  # K

        dP_cold, Re_cold, f_cold = cb.pressure_drop_pipe(T_cold, P, X_air, v, D, L)
        dP_hot, Re_hot, f_hot = cb.pressure_drop_pipe(T_hot, P, X_air, v, D, L)

        # Higher temperature -> lower density, higher viscosity
        # Net effect: lower Re, higher f, but lower dP (density effect dominates)
        assert Re_hot < Re_cold
        assert dP_hot < dP_cold

    def test_pressure_effect(self):
        """Test effect of pressure on pressure drop (via density)."""
        T = 300.0  # K
        X_air = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.1  # m
        L = 100.0  # m

        P_low = 50000.0  # Pa (0.5 atm)
        P_high = 200000.0  # Pa (2 atm)

        dP_low, Re_low, f_low = cb.pressure_drop_pipe(T, P_low, X_air, v, D, L)
        dP_high, Re_high, f_high = cb.pressure_drop_pipe(T, P_high, X_air, v, D, L)

        # Higher pressure -> higher density
        # Higher density -> higher Re, lower f, higher dP
        assert Re_high > Re_low
        assert dP_high > dP_low

    def test_default_parameters(self):
        """Test default parameters (smooth pipe, haaland correlation)."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.1  # m
        L = 100.0  # m

        # Call with defaults
        dP_default, Re_default, f_default = cb.pressure_drop_pipe(T, P, X_air, v, D, L)

        # Call with explicit defaults
        dP_explicit, Re_explicit, f_explicit = cb.pressure_drop_pipe(
            T, P, X_air, v, D, L, roughness=0.0, correlation="haaland"
        )

        # Should be identical
        assert dP_default == dP_explicit
        assert Re_default == Re_explicit
        assert f_default == f_explicit

    def test_invalid_correlation(self):
        """Test error handling for invalid correlation name."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.1  # m
        L = 100.0  # m

        with pytest.raises(ValueError, match="unknown correlation"):
            cb.pressure_drop_pipe(T, P, X_air, v, D, L, correlation="invalid")

    def test_zero_length(self):
        """Test zero length gives zero pressure drop."""
        T = 300.0  # K
        P = 101325.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 10.0  # m/s
        D = 0.1  # m
        L = 0.0  # m

        dP, Re, f = cb.pressure_drop_pipe(T, P, X_air, v, D, L)

        # Zero length -> zero pressure drop
        assert dP == 0.0

        # But Re and f should still be calculated
        assert Re > 0
        assert f > 0

    def test_consistency_with_existing_functions(self):
        """Verify composite matches using existing functions separately."""
        T = 320.0  # K
        P = 150000.0  # Pa
        X_air = cb.standard_dry_air_composition()
        v = 18.0  # m/s
        D = 0.12  # m
        L = 80.0  # m
        roughness = 0.0001  # m

        # Use composite function
        dP_comp, Re_comp, f_comp = cb.pressure_drop_pipe(
            T, P, X_air, v, D, L, roughness, "serghides"
        )

        # Use existing functions step by step
        rho = cb.density(T, P, X_air)
        Re_manual = cb.reynolds(T, P, X_air, v, D)
        e_D = roughness / D
        f_manual = cb.friction_serghides(Re_manual, e_D)
        dP_manual = cb.pipe_dP(v, L, D, f_manual, rho)

        # Perfect match
        assert abs(dP_comp - dP_manual) < 1e-10
        assert abs(Re_comp - Re_manual) < 1e-10
        assert abs(f_comp - f_manual) < 1e-14
