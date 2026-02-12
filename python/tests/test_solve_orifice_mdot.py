"""Tests for solve_orifice_mdot iterative solver.

Verifies the iterative solver that accounts for Cd-Re_D coupling
in orifice flow calculations.
"""

import pytest

import combaero as cb


class TestSolveOrificeMdot:
    """Test iterative solver for orifice mass flow with Cd-Re coupling."""

    def test_incompressible_convergence(self):
        """Test basic convergence for incompressible flow."""
        # Create orifice geometry
        geom = cb.OrificeGeometry()
        geom.d = 0.05  # 50 mm orifice
        geom.D = 0.1  # 100 mm pipe

        # Flow conditions (air at standard conditions)
        dP = 10000.0  # Pa (10 kPa)
        rho = 1.2  # kg/m3
        mu = 1.8e-5  # Pa*s

        # Solve
        mdot = cb.solve_orifice_mdot(geom, dP, rho, mu)

        # Should return a positive mass flow
        assert mdot > 0

        # Typical range for these conditions
        assert 0.1 < mdot < 1.0

    def test_zero_pressure_drop(self):
        """Test that zero pressure drop gives zero flow."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        mdot = cb.solve_orifice_mdot(geom, 0.0, 1.2, 1.8e-5)

        assert mdot == 0.0

    def test_compressible_flow(self):
        """Test with compressible flow (kappa > 1)."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        dP = 10000.0  # Pa
        rho = 1.2  # kg/m3
        mu = 1.8e-5  # Pa*s
        P_upstream = 101325.0  # Pa
        kappa = 1.4  # Air

        mdot_comp = cb.solve_orifice_mdot(geom, dP, rho, mu, P_upstream, kappa)

        # Compressible flow should give slightly lower mdot due to expansion
        mdot_incomp = cb.solve_orifice_mdot(geom, dP, rho, mu, P_upstream, 0.0)

        assert mdot_comp > 0
        assert mdot_comp < mdot_incomp

    def test_consistency_with_direct_calculation(self):
        """Verify solver matches direct calculation when Cd is known."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        dP = 5000.0  # Pa
        rho = 1.2  # kg/m3
        mu = 1.8e-5  # Pa*s

        # Solve iteratively
        mdot_solved = cb.solve_orifice_mdot(geom, dP, rho, mu)

        # Calculate Re_D from solved mdot
        Re_D = (4.0 * mdot_solved) / (3.14159 * geom.D * mu)

        # Use OrificeState to calculate Cd
        state = cb.OrificeState()
        state.Re_D = Re_D
        state.dP = dP
        state.rho = rho
        state.mu = mu
        Cd = cb.Cd_sharp_thin_plate(geom, state)

        # Direct calculation with this Cd
        mdot_direct = cb.orifice_mdot_Cd(geom, Cd, dP, rho)

        # Should match within tolerance (relaxed due to iteration convergence)
        # The solver converges mdot, but Cd changes slightly between iterations
        assert abs(mdot_solved - mdot_direct) / mdot_solved < 0.01  # 1%

    def test_different_correlations(self):
        """Test that different Cd correlations work."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        dP = 10000.0
        rho = 1.2
        mu = 1.8e-5

        # Test all supported correlations
        correlations = [
            cb.CdCorrelation.ReaderHarrisGallagher,
            cb.CdCorrelation.Stolz,
            cb.CdCorrelation.Miller,
        ]

        results = {}
        for corr in correlations:
            mdot = cb.solve_orifice_mdot(geom, dP, rho, mu, correlation=corr)
            results[corr] = mdot
            assert mdot > 0

        # All should give similar results (within 10%)
        mdot_values = list(results.values())
        mdot_mean = sum(mdot_values) / len(mdot_values)
        for mdot in mdot_values:
            assert abs(mdot - mdot_mean) / mdot_mean < 0.10

    def test_beta_dependence(self):
        """Test that larger beta gives more flow."""
        dP = 10000.0
        rho = 1.2
        mu = 1.8e-5
        D = 0.1

        # Small beta
        geom_small = cb.OrificeGeometry()
        geom_small.d = 0.03
        geom_small.D = D

        # Large beta
        geom_large = cb.OrificeGeometry()
        geom_large.d = 0.07
        geom_large.D = D

        mdot_small = cb.solve_orifice_mdot(geom_small, dP, rho, mu)
        mdot_large = cb.solve_orifice_mdot(geom_large, dP, rho, mu)

        # Larger orifice -> more flow
        assert mdot_large > mdot_small

    def test_pressure_drop_dependence(self):
        """Test that higher dP gives more flow."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        rho = 1.2
        mu = 1.8e-5

        dP_low = 5000.0  # Pa
        dP_high = 20000.0  # Pa

        mdot_low = cb.solve_orifice_mdot(geom, dP_low, rho, mu)
        mdot_high = cb.solve_orifice_mdot(geom, dP_high, rho, mu)

        # Higher pressure drop -> more flow
        assert mdot_high > mdot_low

        # Should scale roughly as sqrt(dP)
        ratio = mdot_high / mdot_low
        expected_ratio = (dP_high / dP_low) ** 0.5
        assert abs(ratio - expected_ratio) / expected_ratio < 0.05

    def test_custom_tolerance(self):
        """Test that custom tolerance works."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        dP = 10000.0
        rho = 1.2
        mu = 1.8e-5

        # Tight tolerance
        mdot_tight = cb.solve_orifice_mdot(geom, dP, rho, mu, tol=1e-8)

        # Loose tolerance
        mdot_loose = cb.solve_orifice_mdot(geom, dP, rho, mu, tol=1e-4)

        # Both should converge to similar values
        assert abs(mdot_tight - mdot_loose) / mdot_tight < 1e-3

    def test_invalid_geometry(self):
        """Test error handling for invalid geometry."""
        geom = cb.OrificeGeometry()
        geom.d = 0.15  # Larger than pipe!
        geom.D = 0.1

        with pytest.raises(ValueError, match="invalid geometry"):
            cb.solve_orifice_mdot(geom, 10000.0, 1.2, 1.8e-5)

    def test_invalid_parameters(self):
        """Test error handling for invalid parameters."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        # Negative dP
        with pytest.raises(ValueError, match="dP must be non-negative"):
            cb.solve_orifice_mdot(geom, -10000.0, 1.2, 1.8e-5)

        # Zero density
        with pytest.raises(ValueError, match="rho must be positive"):
            cb.solve_orifice_mdot(geom, 10000.0, 0.0, 1.8e-5)

        # Zero viscosity
        with pytest.raises(ValueError, match="mu must be positive"):
            cb.solve_orifice_mdot(geom, 10000.0, 1.2, 0.0)

    def test_high_pressure_ratio_compressible(self):
        """Test compressible flow with significant pressure ratio."""
        geom = cb.OrificeGeometry()
        geom.d = 0.05
        geom.D = 0.1

        dP = 25000.0  # Pa (25% pressure drop)
        P_upstream = 100000.0  # Pa
        rho = 1.2  # kg/m3
        mu = 1.8e-5  # Pa*s
        kappa = 1.4  # Air

        mdot = cb.solve_orifice_mdot(geom, dP, rho, mu, P_upstream, kappa)

        # Should still converge
        assert mdot > 0

        # Expansibility effect should be noticeable
        mdot_incomp = cb.solve_orifice_mdot(geom, dP, rho, mu, P_upstream, 0.0)

        # Compressible should be 5-10% lower
        reduction = (mdot_incomp - mdot) / mdot_incomp
        assert 0.05 < reduction < 0.15

    def test_realistic_natural_gas_metering(self):
        """Test with realistic natural gas metering conditions."""
        # ISO 5167 typical application
        geom = cb.OrificeGeometry()
        geom.d = 0.06  # 60 mm
        geom.D = 0.1  # 100 mm (beta = 0.6)

        dP = 15000.0  # Pa
        P_upstream = 5000000.0  # Pa (50 bar)
        rho = 40.0  # kg/m3 (high pressure gas)
        mu = 1.1e-5  # Pa*s (natural gas)
        kappa = 1.3  # Natural gas

        mdot = cb.solve_orifice_mdot(geom, dP, rho, mu, P_upstream, kappa)

        # Should converge to reasonable value
        assert mdot > 0
        assert 1.0 < mdot < 50.0  # kg/s (typical range)
