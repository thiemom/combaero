"""Tests for Bassett 2001 tee junction K functions and solver interface."""

import math

import pytest

import combaero._core as _core

# ---------------------------------------------------------------------------
# K function reference values from Bassett 2001
# ---------------------------------------------------------------------------


def test_K5_formula():
    """K5(q) = q^2 - 1.5*q + 0.5. Independent of theta and psi."""
    # K5(0.5) = 0.25 - 0.75 + 0.5 = 0.0
    assert abs(_core.tee_K5(0.5) - 0.0) < 1e-12
    assert abs(_core.tee_K5(0.0) - 0.5) < 1e-12
    assert abs(_core.tee_K5(1.0) - 0.0) < 1e-12


def test_K6_at_90deg_psi1():
    """K6 at theta=90deg, psi=1, q=0.5 should match Bassett Fig. 7a (~0.87)."""
    theta = math.pi / 2.0
    psi = 1.0
    result = _core.tee_K6(0.5, psi, theta)
    assert abs(result - 0.8673) < 5e-4


def test_K12_at_90deg_psi1():
    """K12 at theta=90deg, psi=1, q=0.5 should match Bassett Fig. 10b (~0.35)."""
    theta = math.pi / 2.0
    psi = 1.0
    result = _core.tee_K12(0.5, psi, theta)
    assert abs(result - 0.348) < 5e-3


# ---------------------------------------------------------------------------
# Blend functions
# ---------------------------------------------------------------------------


def test_blend_weight_boundaries():
    """blend_weight(0) = 0.5, (+inf) -> 1, (-inf) -> 0."""
    assert abs(_core.tee_blend_weight(0.0) - 0.5) < 1e-12
    assert _core.tee_blend_weight(10.0) > 0.999
    assert _core.tee_blend_weight(-10.0) < 0.001


def test_blend_agreement_at_zero_flow():
    """At q=0 merging and branching K must be identical (Test 5 from spec)."""
    tol = 1e-10
    for theta in [math.pi / 6, math.pi / 4, math.pi / 3, math.pi / 2]:
        for psi in [0.5, 1.0, 2.0]:
            assert (
                abs(
                    _core.merging_tee_K_straight(0.0, psi, theta)
                    - _core.branching_tee_K_straight(0.0, psi, theta)
                )
                < tol
            )
            assert (
                abs(
                    _core.merging_tee_K_branch(0.0, psi, theta)
                    - _core.branching_tee_K_branch(0.0, psi, theta)
                )
                < tol
            )


# ---------------------------------------------------------------------------
# Validity check
# ---------------------------------------------------------------------------


def test_tee_check_inputs_valid():
    result = _core.tee_check_inputs(0.5, 1.0, math.pi / 2)
    assert result["valid"] is True
    assert result["q_in_range"] is True
    assert result["psi_valid"] is True
    assert result["theta_valid"] is True


def test_tee_check_inputs_invalid_q():
    result = _core.tee_check_inputs(1.5, 1.0, math.pi / 2)
    assert result["valid"] is False
    assert result["q_in_range"] is False


def test_tee_check_inputs_invalid_psi():
    result = _core.tee_check_inputs(0.5, 0.01, math.pi / 2)
    assert result["psi_valid"] is False


def test_tee_check_inputs_invalid_theta():
    result = _core.tee_check_inputs(0.5, 1.0, math.pi)
    assert result["theta_valid"] is False


# ---------------------------------------------------------------------------
# Robustness: no NaN outside valid range
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("psi", [1e-6, 1e-3, 0.04])
def test_small_psi_no_nan(psi):
    """Effective K functions must not produce NaN for tiny psi."""
    theta = math.pi / 2
    for q in [-0.5, 0.0, 0.5, 1.0]:
        assert math.isfinite(_core.merging_tee_K_straight(q, psi, theta))
        assert math.isfinite(_core.merging_tee_K_branch(q, psi, theta))


def test_extreme_theta_no_nan():
    """K functions must not produce NaN for theta outside validated range."""
    for theta in [0.01, math.pi, 2 * math.pi]:
        assert math.isfinite(_core.merging_tee_K_straight(0.5, 1.0, theta))
        assert math.isfinite(_core.merging_tee_K_branch(0.5, 1.0, theta))


# ---------------------------------------------------------------------------
# Solver interface: residuals and Jacobians
# ---------------------------------------------------------------------------


def _dry_air_Y():
    Y = [0.0] * 15
    Y[0] = 0.767  # N2
    Y[1] = 0.233  # O2
    return Y


def test_merging_tee_residuals_basic():
    """Merging tee returns finite residuals and Jacobians for typical conditions."""
    Y = _dry_air_Y()
    res = _core.merging_tee_residuals_and_jacobian(
        m_dot_com=0.5,
        m_dot_branch=0.2,
        dP0_straight=150.0,
        dP0_branch=200.0,
        P_static_com=2e5,
        T_com=600.0,
        Y_com=Y,
        theta=math.pi / 2,
        psi=1.0,
        F_C=0.01,
    )
    assert math.isfinite(res.R_straight)
    assert math.isfinite(res.R_branch)
    assert math.isfinite(res.dR_straight_d_mdot_com)
    assert math.isfinite(res.dR_branch_d_mdot_com)
    assert math.isfinite(res.dR_straight_d_mdot_branch)
    assert math.isfinite(res.dR_branch_d_mdot_branch)
    assert 0.0 <= res.q <= 1.0
    assert res.topology_valid


def test_branching_tee_residuals_basic():
    """Branching tee returns finite residuals for typical conditions."""
    Y = _dry_air_Y()
    res = _core.branching_tee_residuals_and_jacobian(
        m_dot_com=0.5,
        m_dot_branch=0.15,
        dP0_straight=100.0,
        dP0_branch=180.0,
        P_static_com=1.5e5,
        T_com=500.0,
        Y_com=Y,
        theta=math.pi / 3,
        psi=1.2,
        F_C=0.008,
    )
    assert math.isfinite(res.R_straight)
    assert math.isfinite(res.R_branch)
    assert res.topology_valid


def test_zero_flow_no_nan():
    """Zero total flow must return finite residuals (regularisation test)."""
    Y = _dry_air_Y()
    res = _core.merging_tee_residuals_and_jacobian(
        0.0, 0.0, 50.0, 50.0, 1e5, 400.0, Y, math.pi / 2, 1.0, 0.01
    )
    assert math.isfinite(res.R_straight)
    assert math.isfinite(res.R_branch)
    assert math.isfinite(res.dR_straight_d_mdot_com)


def test_topology_valid_flag_reversed():
    """Negative branch flow should flag topology_valid=False."""
    Y = _dry_air_Y()
    res = _core.merging_tee_residuals_and_jacobian(
        m_dot_com=0.5,
        m_dot_branch=-0.2,
        dP0_straight=50.0,
        dP0_branch=50.0,
        P_static_com=1e5,
        T_com=400.0,
        Y_com=Y,
        theta=math.pi / 2,
        psi=1.0,
        F_C=0.01,
    )
    assert not res.topology_valid


def test_status_extrapolated_outside_range():
    """Results outside validated range (q>1 equivalent) should have EXTRAPOLATED status."""
    Y = _dry_air_Y()
    # Very large branch flow relative to common (q > 1 via negative straight flow)
    res = _core.merging_tee_residuals_and_jacobian(
        m_dot_com=0.1,
        m_dot_branch=0.5,
        dP0_straight=50.0,
        dP0_branch=50.0,
        P_static_com=1e5,
        T_com=400.0,
        Y_com=Y,
        theta=math.pi / 2,
        psi=1.0,
        F_C=0.01,
    )
    # q = 0.5/0.1 = 5.0 > 1.0 -> status should be EXTRAPOLATED
    assert res.status == _core.CorrelationValidity.EXTRAPOLATED


def test_m_dot_jacobians_by_fd():
    """Verify m_dot Jacobians against finite differences."""
    Y = _dry_air_Y()
    mc = 0.5
    mb = 0.2
    dP0_s = 150.0
    dP0_b = 200.0
    P = 2e5
    T = 600.0
    theta = math.pi / 2
    psi = 1.0
    F_C = 0.01
    eps = 1e-5

    res0 = _core.merging_tee_residuals_and_jacobian(mc, mb, dP0_s, dP0_b, P, T, Y, theta, psi, F_C)

    # FD for m_dot_com
    r_p = _core.merging_tee_residuals_and_jacobian(
        mc + eps, mb, dP0_s, dP0_b, P, T, Y, theta, psi, F_C
    )
    r_m = _core.merging_tee_residuals_and_jacobian(
        mc - eps, mb, dP0_s, dP0_b, P, T, Y, theta, psi, F_C
    )
    fd_Rs = (r_p.R_straight - r_m.R_straight) / (2 * eps)
    fd_Rb = (r_p.R_branch - r_m.R_branch) / (2 * eps)

    # Allow 1e-4 relative + 0.5 Pa/(kg/s) absolute tolerance
    tol_s = abs(fd_Rs) * 1e-4 + 0.5
    tol_b = abs(fd_Rb) * 1e-4 + 0.5
    assert abs(res0.dR_straight_d_mdot_com - fd_Rs) < tol_s
    assert abs(res0.dR_branch_d_mdot_com - fd_Rb) < tol_b
