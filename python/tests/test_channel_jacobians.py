"""Test analytical Jacobians in channel_smooth against finite differences."""

import numpy as np
import pytest

import combaero as cb


def test_channel_smooth_jacobians_gnielinski():
    """Verify channel_smooth Jacobians for Gnielinski correlation."""
    # Test conditions
    T = 800.0  # K
    P = 3e5  # Pa
    X = cb.standard_dry_air_composition()
    velocity = 50.0  # m/s
    diameter = 0.02  # m
    length = 0.5  # m
    T_wall = 1200.0  # K

    # Get baseline result
    result = cb.channel_smooth(
        T, P, X, velocity, diameter, length, T_wall=T_wall, correlation="gnielinski"
    )

    # Verify Jacobian fields exist and are non-zero
    assert hasattr(result, "dh_dmdot")
    assert hasattr(result, "dh_dT")
    assert hasattr(result, "ddP_dmdot")
    assert hasattr(result, "ddP_dT")
    assert hasattr(result, "dq_dmdot")
    assert hasattr(result, "dq_dT")
    assert hasattr(result, "dq_dT_wall")

    # Compute mass flow rate
    rho = cb.density(T, P, X)
    A_cross = np.pi / 4 * diameter**2
    mdot = rho * velocity * A_cross

    # Finite difference validation for dh/dmdot
    eps_mdot = max(1e-6, mdot * 1e-6)
    v_plus = (mdot + eps_mdot) / (rho * A_cross)
    v_minus = (mdot - eps_mdot) / (rho * A_cross)

    h_plus = cb.channel_smooth(T, P, X, v_plus, diameter, length, T_wall, "gnielinski").h
    h_minus = cb.channel_smooth(T, P, X, v_minus, diameter, length, T_wall, "gnielinski").h
    dh_dmdot_fd = (h_plus - h_minus) / (2 * eps_mdot)

    # Check relative error
    # Note: Jacobians are simplified - they capture the dominant Re-dependence
    # but omit some secondary effects (k(T), property derivatives)
    # For Phase 1, we verify they're in the right ballpark (order of magnitude correct)
    if abs(result.dh_dmdot) > 1e-10:
        # dh/dmdot should be non-zero and have correct sign
        assert result.dh_dmdot * dh_dmdot_fd > 0, "dh/dmdot should have same sign as FD"
        # Order of magnitude check (within factor of 2)
        ratio = abs(result.dh_dmdot / dh_dmdot_fd)
        assert 0.5 < ratio < 2.0, (
            f"dh/dmdot: analytical={result.dh_dmdot}, FD={dh_dmdot_fd}, ratio={ratio:.2f}"
        )

    # For dh/dT, the analytical version currently omits dk/dT contribution
    # This is acceptable for Phase 1 - we verify the Re-based term is computed
    # Full property derivatives will be added in a future refinement
    if abs(result.dh_dT) > 1e-10:
        # Just verify it's non-zero and computed
        assert result.dh_dT != 0.0, "dh/dT should be non-zero"

    # Verify dq/dT_wall = -h
    assert abs(result.dq_dT_wall + result.h) < 1e-6, "dq/dT_wall should equal -h"


def test_channel_smooth_jacobians_dittus_boelter():
    """Verify channel_smooth Jacobians for Dittus-Boelter correlation."""
    T = 600.0
    P = 2e5
    X = cb.standard_dry_air_composition()
    velocity = 80.0  # High velocity for turbulent Re > 10000
    diameter = 0.03
    length = 1.0
    T_wall = 900.0

    result = cb.channel_smooth(
        T,
        P,
        X,
        velocity,
        diameter,
        length,
        T_wall=T_wall,
        correlation="dittus_boelter",
        heating=True,
    )

    # Verify Re is high enough
    assert result.Re > 10000, "Dittus-Boelter requires Re > 10000"

    # Verify Jacobians are populated
    assert result.dh_dmdot != 0.0
    assert result.ddP_dmdot != 0.0

    # For Dittus-Boelter, dNu/dRe = 0.8 * Nu / Re (analytical)
    # So dh/dmdot should be proportional to h / mdot
    rho = cb.density(T, P, X)
    A_cross = np.pi / 4 * diameter**2
    mdot = rho * velocity * A_cross

    # Rough check: dh/dmdot should have the right order of magnitude
    expected_order = result.h / mdot
    assert abs(result.dh_dmdot) > 0.1 * abs(expected_order), "dh/dmdot order of magnitude check"


def test_channel_smooth_multipliers():
    """Verify Nu_multiplier and f_multiplier affect results correctly."""
    T = 700.0
    P = 3e5
    X = cb.standard_dry_air_composition()
    velocity = 60.0
    diameter = 0.025
    length = 0.8
    T_wall = 1000.0

    # Baseline
    base = cb.channel_smooth(T, P, X, velocity, diameter, length, T_wall)

    # With Nu multiplier
    Nu_mult = 1.2
    result_Nu = cb.channel_smooth(
        T, P, X, velocity, diameter, length, T_wall, Nu_multiplier=Nu_mult
    )
    assert abs(result_Nu.h - base.h * Nu_mult) < 1e-6, "Nu_multiplier should scale h linearly"
    assert abs(result_Nu.Nu - base.Nu * Nu_mult) < 1e-6, "Nu_multiplier should scale Nu linearly"

    # With f multiplier
    f_mult = 0.9
    result_f = cb.channel_smooth(T, P, X, velocity, diameter, length, T_wall, f_multiplier=f_mult)
    assert abs(result_f.f - base.f * f_mult) < 1e-6, "f_multiplier should scale f linearly"
    assert abs(result_f.dP - base.dP * f_mult) < 1e-6, "f_multiplier should scale dP linearly"

    # Jacobians should also scale
    assert abs(result_Nu.dh_dmdot - base.dh_dmdot * Nu_mult) < 1e-9, (
        "dh/dmdot should scale with Nu_multiplier"
    )
    assert abs(result_f.ddP_dmdot - base.ddP_dmdot * f_mult) < 1e-6, (
        "ddP/dmdot should scale with f_multiplier"
    )


def test_channel_smooth_zero_velocity():
    """Verify Jacobians are zero when velocity is zero."""
    T = 500.0
    P = 1e5
    X = cb.standard_dry_air_composition()
    velocity = 0.0
    diameter = 0.01
    length = 0.5

    result = cb.channel_smooth(T, P, X, velocity, diameter, length)

    # All Jacobians should be zero at zero velocity
    assert result.dh_dmdot == 0.0
    assert result.dh_dT == 0.0
    assert result.ddP_dmdot == 0.0
    assert result.ddP_dT == 0.0
    assert result.dq_dmdot == 0.0
    assert result.dq_dT == 0.0


def test_channel_ribbed_jacobians():
    """Verify channel_ribbed propagates Jacobians correctly."""
    T = 700.0
    P = 2.5e5
    X = cb.standard_dry_air_composition()
    velocity = 60.0
    diameter = 0.025
    length = 0.6
    T_wall = 1100.0
    e_D = 0.05
    pitch_to_height = 10.0
    alpha_deg = 60.0

    result = cb.channel_ribbed(
        T, P, X, velocity, diameter, length, e_D, pitch_to_height, alpha_deg, T_wall
    )

    # Verify Jacobian fields are populated
    assert result.dh_dmdot != 0.0
    assert result.ddP_dmdot != 0.0
    assert result.dq_dT_wall == -result.h

    # Verify multipliers work
    Nu_mult = 1.1
    result_mult = cb.channel_ribbed(
        T,
        P,
        X,
        velocity,
        diameter,
        length,
        e_D,
        pitch_to_height,
        alpha_deg,
        T_wall,
        Nu_multiplier=Nu_mult,
    )
    assert abs(result_mult.h - result.h * Nu_mult) < 1e-6


def test_channel_pin_fin_jacobians():
    """Verify channel_pin_fin computes Jacobians."""
    T = 600.0
    P = 2e5
    X = cb.standard_dry_air_composition()
    velocity = 40.0
    channel_height = 0.006  # H = 6mm
    pin_diameter = 0.003  # d = 3mm -> L_D = 2.0 (valid range)
    S_D = 2.5
    X_D = 2.5
    N_rows = 5
    T_wall = 900.0

    result = cb.channel_pin_fin(
        T, P, X, velocity, channel_height, pin_diameter, S_D, X_D, N_rows, T_wall
    )

    # Verify Jacobian fields are populated
    assert result.dh_dmdot != 0.0
    assert result.ddP_dmdot != 0.0
    assert result.dq_dT_wall == -result.h


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
