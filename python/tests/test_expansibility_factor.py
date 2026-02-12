"""Tests for expansibility_factor function.

Verifies the ISO 5167-2:2003 expansibility factor calculation
for compressible gas flow through orifices.
"""

import pytest

import combaero as cb


class TestExpansibilityFactor:
    """Test expansibility_factor for compressible orifice flow."""

    def test_incompressible_limit(self):
        """Test that epsilon = 1.0 for incompressible flow (dP = 0)."""
        beta = 0.5
        dP = 0.0  # No pressure drop
        P_upstream = 101325.0  # Pa
        kappa = 1.4  # Air

        epsilon = cb.expansibility_factor(beta, dP, P_upstream, kappa)

        assert epsilon == 1.0

    def test_kappa_limit(self):
        """Test that epsilon = 1.0 when kappa <= 1 (incompressible)."""
        beta = 0.5
        dP = 10000.0  # Pa
        P_upstream = 101325.0  # Pa
        kappa = 1.0  # Incompressible limit

        epsilon = cb.expansibility_factor(beta, dP, P_upstream, kappa)

        assert epsilon == 1.0

    def test_typical_air_flow(self):
        """Test with typical air flow conditions."""
        beta = 0.5  # d/D = 0.5
        dP = 10000.0  # Pa (10 kPa)
        P_upstream = 101325.0  # Pa (1 atm)
        kappa = 1.4  # Air

        epsilon = cb.expansibility_factor(beta, dP, P_upstream, kappa)

        # Should be less than 1 (compressible correction)
        assert epsilon < 1.0
        # Should be close to 1 for small pressure ratio
        assert epsilon > 0.95
        # Typical value for these conditions
        assert 0.96 < epsilon < 0.99

    def test_higher_pressure_drop(self):
        """Test with higher pressure drop (larger compressibility effect)."""
        beta = 0.6
        dP = 25000.0  # Pa (25% of upstream pressure)
        P_upstream = 100000.0  # Pa
        kappa = 1.4

        epsilon = cb.expansibility_factor(beta, dP, P_upstream, kappa)

        # Larger pressure drop -> more expansion -> lower epsilon
        assert epsilon < 0.95
        assert epsilon > 0.85

    def test_beta_dependence(self):
        """Test that epsilon decreases with increasing beta."""
        dP = 10000.0  # Pa
        P_upstream = 101325.0  # Pa
        kappa = 1.4

        beta_small = 0.3
        beta_large = 0.7

        eps_small = cb.expansibility_factor(beta_small, dP, P_upstream, kappa)
        eps_large = cb.expansibility_factor(beta_large, dP, P_upstream, kappa)

        # Larger beta -> more expansion -> lower epsilon
        assert eps_large < eps_small

    def test_pressure_ratio_dependence(self):
        """Test that epsilon decreases with increasing pressure ratio."""
        beta = 0.5
        P_upstream = 101325.0  # Pa
        kappa = 1.4

        dP_low = 5000.0  # Pa (5% pressure drop)
        dP_high = 20000.0  # Pa (20% pressure drop)

        eps_low = cb.expansibility_factor(beta, dP_low, P_upstream, kappa)
        eps_high = cb.expansibility_factor(beta, dP_high, P_upstream, kappa)

        # Higher pressure drop -> more expansion -> lower epsilon
        assert eps_high < eps_low

    def test_kappa_dependence(self):
        """Test that epsilon varies with kappa (gas type)."""
        beta = 0.5
        dP = 15000.0  # Pa
        P_upstream = 101325.0  # Pa

        kappa_air = 1.4  # Air
        kappa_co2 = 1.3  # CO2 (more polyatomic)

        eps_air = cb.expansibility_factor(beta, dP, P_upstream, kappa_air)
        eps_co2 = cb.expansibility_factor(beta, dP, P_upstream, kappa_co2)

        # Different kappa -> different expansion behavior
        assert eps_air != eps_co2

    def test_iso_5167_example(self):
        """Test with conditions similar to ISO 5167 examples."""
        # Typical natural gas metering conditions
        beta = 0.6
        dP = 15000.0  # Pa
        P_upstream = 5000000.0  # Pa (50 bar)
        kappa = 1.3  # Natural gas

        epsilon = cb.expansibility_factor(beta, dP, P_upstream, kappa)

        # Small pressure ratio -> epsilon close to 1
        assert 0.99 < epsilon < 1.0

    def test_boundary_beta_values(self):
        """Test with beta near validity boundaries."""
        dP = 10000.0  # Pa
        P_upstream = 101325.0  # Pa
        kappa = 1.4

        # Near lower limit (β = 0.1)
        eps_low = cb.expansibility_factor(0.15, dP, P_upstream, kappa)
        assert 0.95 < eps_low <= 1.0

        # Near upper limit (β = 0.75)
        eps_high = cb.expansibility_factor(0.7, dP, P_upstream, kappa)
        assert 0.90 < eps_high < 1.0

    def test_invalid_beta(self):
        """Test error handling for invalid beta values."""
        dP = 10000.0
        P_upstream = 101325.0
        kappa = 1.4

        # beta <= 0
        with pytest.raises(ValueError, match="beta must be in range"):
            cb.expansibility_factor(0.0, dP, P_upstream, kappa)

        # beta >= 1
        with pytest.raises(ValueError, match="beta must be in range"):
            cb.expansibility_factor(1.0, dP, P_upstream, kappa)

    def test_invalid_pressure(self):
        """Test error handling for invalid pressures."""
        beta = 0.5
        kappa = 1.4

        # Negative upstream pressure
        with pytest.raises(ValueError, match="P_upstream must be positive"):
            cb.expansibility_factor(beta, 10000.0, -101325.0, kappa)

        # Negative dP
        with pytest.raises(ValueError, match="dP must be non-negative"):
            cb.expansibility_factor(beta, -10000.0, 101325.0, kappa)

    def test_consistency_with_formula(self):
        """Verify implementation matches ISO 5167-2 formula exactly."""
        beta = 0.5
        dP = 10000.0  # Pa
        P_upstream = 101325.0  # Pa
        kappa = 1.4

        # Manual calculation
        tau = dP / P_upstream
        beta4 = beta**4
        beta8 = beta**8
        coeff = 0.351 + 0.256 * beta4 + 0.93 * beta8
        expansion_term = 1.0 - (1.0 - tau) ** (1.0 / kappa)
        epsilon_expected = 1.0 - coeff * expansion_term

        epsilon_actual = cb.expansibility_factor(beta, dP, P_upstream, kappa)

        # Should match within machine precision
        assert abs(epsilon_actual - epsilon_expected) < 1e-14

    def test_range_validity(self):
        """Test that function works across full valid range."""
        P_upstream = 101325.0
        kappa = 1.4

        # Test various beta values
        for beta in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
            # Test various pressure ratios
            for tau in [0.05, 0.10, 0.15, 0.20, 0.25]:
                dP = tau * P_upstream
                epsilon = cb.expansibility_factor(beta, dP, P_upstream, kappa)

                # All should be valid
                assert 0.8 < epsilon <= 1.0
                # Epsilon should decrease with tau
                # (more pressure drop -> more expansion)
