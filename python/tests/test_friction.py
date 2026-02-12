"""Tests for friction factor correlations.

Verifies Python bindings for friction factor functions match C++ implementation.
All friction factors are dimensionless [-].
"""

import pytest

import combaero as cb


class TestFrictionFactors:
    """Test friction factor correlation bindings."""

    def test_friction_haaland_smooth_pipe(self):
        """Test Haaland correlation for smooth pipe (e_D = 0)."""
        Re = 1e5  # Reynolds number [-]
        e_D = 0.0  # Smooth pipe (relative roughness) [-]

        f = cb.friction_haaland(Re, e_D)

        # Friction factor should be dimensionless and positive
        assert f > 0
        assert isinstance(f, float)

        # For smooth pipe at Re=1e5, f should be around 0.018
        # (from Moody diagram)
        assert 0.015 < f < 0.025

    def test_friction_haaland_rough_pipe(self):
        """Test Haaland correlation for rough pipe."""
        Re = 1e5  # Reynolds number [-]
        e_D = 0.001  # Relative roughness [-]

        f = cb.friction_haaland(Re, e_D)

        # Rough pipe should have higher friction than smooth
        f_smooth = cb.friction_haaland(Re, 0.0)
        assert f > f_smooth

        # Reasonable range for this roughness
        assert 0.020 < f < 0.030

    def test_friction_serghides_smooth_pipe(self):
        """Test Serghides correlation for smooth pipe."""
        Re = 1e5  # Reynolds number [-]
        e_D = 0.0  # Smooth pipe [-]

        f = cb.friction_serghides(Re, e_D)

        assert f > 0
        assert isinstance(f, float)

        # Should be very close to Haaland for smooth pipe
        f_haaland = cb.friction_haaland(Re, e_D)
        assert abs(f - f_haaland) / f_haaland < 0.05  # Within 5%

    def test_friction_serghides_vs_haaland(self):
        """Test Serghides is more accurate than Haaland."""
        Re = 5e4  # Reynolds number [-]
        e_D = 0.002  # Relative roughness [-]

        f_serghides = cb.friction_serghides(Re, e_D)
        f_haaland = cb.friction_haaland(Re, e_D)

        # Both should be positive and similar
        assert f_serghides > 0
        assert f_haaland > 0

        # Serghides claims <0.3% accuracy vs Colebrook
        # Haaland claims ~2-3% accuracy
        # They should be close but not identical
        assert abs(f_serghides - f_haaland) / f_haaland < 0.03

    def test_friction_colebrook_smooth_pipe(self):
        """Test Colebrook-White for smooth pipe."""
        Re = 1e5  # Reynolds number [-]
        e_D = 0.0  # Smooth pipe [-]

        f = cb.friction_colebrook(Re, e_D)

        assert f > 0
        assert isinstance(f, float)

        # Should match Serghides very closely (both very accurate)
        f_serghides = cb.friction_serghides(Re, e_D)
        assert abs(f - f_serghides) / f_serghides < 0.005  # Within 0.5%

    def test_friction_colebrook_convergence(self):
        """Test Colebrook-White converges with custom tolerance."""
        Re = 1e5  # Reynolds number [-]
        e_D = 0.001  # Relative roughness [-]

        # Test with different tolerances
        f_tight = cb.friction_colebrook(Re, e_D, tol=1e-12, max_iter=50)
        f_loose = cb.friction_colebrook(Re, e_D, tol=1e-6, max_iter=10)

        # Both should converge to similar values
        assert abs(f_tight - f_loose) / f_tight < 1e-5

    def test_friction_petukhov_smooth_pipe(self):
        """Test Petukhov correlation for smooth pipe."""
        Re = 1e5  # Reynolds number [-]

        f = cb.friction_petukhov(Re)

        assert f > 0
        assert isinstance(f, float)

        # Should match smooth pipe Colebrook
        f_colebrook = cb.friction_colebrook(Re, 0.0)
        assert abs(f - f_colebrook) / f_colebrook < 0.01  # Within 1%

    def test_friction_reynolds_number_effect(self):
        """Test friction factor decreases with Reynolds number."""
        e_D = 0.001  # Fixed roughness [-]

        Re_low = 1e4  # Reynolds number [-]
        Re_high = 1e6  # Reynolds number [-]

        f_low = cb.friction_haaland(Re_low, e_D)
        f_high = cb.friction_haaland(Re_high, e_D)

        # Higher Re should have lower friction (in turbulent regime)
        assert f_low > f_high

    def test_friction_roughness_effect(self):
        """Test friction factor increases with roughness."""
        Re = 1e5  # Fixed Reynolds number [-]

        e_D_smooth = 0.0  # Smooth pipe [-]
        e_D_rough = 0.01  # Rough pipe [-]

        f_smooth = cb.friction_haaland(Re, e_D_smooth)
        f_rough = cb.friction_haaland(Re, e_D_rough)

        # Rougher pipe should have higher friction
        assert f_rough > f_smooth

    def test_friction_valid_range(self):
        """Test friction correlations work across valid Reynolds number range."""
        e_D = 0.0  # Smooth pipe [-]

        # Test at transition region (Re ~ 3000)
        Re_transition = 3000
        f_transition = cb.friction_haaland(Re_transition, e_D)
        assert 0.02 < f_transition < 0.06  # Reasonable range

        # Test at moderate turbulent flow (Re ~ 1e4)
        Re_moderate = 1e4
        f_moderate = cb.friction_haaland(Re_moderate, e_D)
        assert 0.015 < f_moderate < 0.04

        # Test at high turbulent flow (Re ~ 1e6)
        Re_high = 1e6
        f_high = cb.friction_haaland(Re_high, e_D)
        assert 0.008 < f_high < 0.015

        # Verify friction decreases with increasing Re
        assert f_transition > f_moderate > f_high

    def test_friction_all_correlations_consistent(self):
        """Test all correlations give consistent results."""
        Re = 1e5  # Reynolds number [-]
        e_D = 0.0005  # Relative roughness [-]

        f_haaland = cb.friction_haaland(Re, e_D)
        f_serghides = cb.friction_serghides(Re, e_D)
        f_colebrook = cb.friction_colebrook(Re, e_D)

        # All should be within 3% of each other
        mean_f = (f_haaland + f_serghides + f_colebrook) / 3

        assert abs(f_haaland - mean_f) / mean_f < 0.03
        assert abs(f_serghides - mean_f) / mean_f < 0.03
        assert abs(f_colebrook - mean_f) / mean_f < 0.03

    def test_friction_units_dimensionless(self):
        """Verify friction factor is dimensionless."""
        Re = 1e5  # Reynolds number [-]
        e_D = 0.001  # Relative roughness [-]

        f = cb.friction_haaland(Re, e_D)

        # Friction factor should be O(0.01-0.1) for typical cases
        # This is a sanity check on units
        assert 0.001 < f < 1.0

        # Darcy friction factor is used in: DeltaP = f * (L/D) * (ρ*v2/2)
        # So f must be dimensionless for units to work out


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
