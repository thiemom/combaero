"""Tests for sharp-edge and conical area change elements."""

import numpy as np
import pytest

import combaero as cb


# ---------------------------------------------------------------------------
# Test 1: Expansion, forward flow — dP > 0, sign correct
# ---------------------------------------------------------------------------
def test_expansion_forward_positive_dp() -> None:
    """Forward flow through expansion (F0 < F1): dP must be positive."""
    res = cb._core.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04)
    assert res.dP > 0, f"Expected dP > 0 for expansion, got {res.dP}"


# ---------------------------------------------------------------------------
# Test 2: Contraction, forward flow — dP > 0, smaller than expansion
# ---------------------------------------------------------------------------
def test_contraction_less_than_expansion() -> None:
    """Contraction loss should be less than expansion at the same conditions."""
    res_exp = cb._core.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04)
    res_con = cb._core.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.04, F1=0.01)
    assert res_con.dP > 0, f"Expected dP > 0 for contraction, got {res_con.dP}"
    assert res_con.dP < res_exp.dP, (
        f"Contraction loss ({res_con.dP:.4f}) should be less than expansion loss ({res_exp.dP:.4f})"
    )


# ---------------------------------------------------------------------------
# Test 3: Reverse flow gives negative dP
# ---------------------------------------------------------------------------
def test_reverse_flow_negative_dp() -> None:
    """Reverse flow (m_dot < 0) through expansion geometry: dP must be negative."""
    res = cb._core.sharp_area_change(m_dot=-0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04)
    assert res.dP < 0, f"Expected dP < 0 for reverse flow, got {res.dP}"


# ---------------------------------------------------------------------------
# Test 4: Jacobian accuracy — finite difference check on dS/dm and dS/drho
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "m_dot, F0, F1, label",
    [
        (0.1, 0.01, 0.04, "expansion_fwd"),
        (0.1, 0.04, 0.01, "contraction_fwd"),
        (-0.1, 0.01, 0.04, "expansion_rev"),
        (-0.1, 0.04, 0.01, "contraction_rev"),
    ],
)
def test_jacobian_fd_accuracy(m_dot: float, F0: float, F1: float, label: str) -> None:
    """Central FD validation for dS_dm and dS_drho."""
    rho = 1.2
    mu = 1.8e-5

    res = cb._core.sharp_area_change(m_dot=m_dot, rho=rho, mu=mu, F0=F0, F1=F1)

    # dS/dm via central FD
    eps_m = 1e-7
    dp_p = cb._core.sharp_area_change(m_dot=m_dot + eps_m, rho=rho, mu=mu, F0=F0, F1=F1).dP
    dp_m = cb._core.sharp_area_change(m_dot=m_dot - eps_m, rho=rho, mu=mu, F0=F0, F1=F1).dP
    fd_dm = (dp_p - dp_m) / (2 * eps_m)
    np.testing.assert_allclose(res.dS_dm, fd_dm, rtol=1e-4, err_msg=f"{label} dS_dm")

    # dS/drho via central FD
    eps_rho = 1e-6
    dp_p = cb._core.sharp_area_change(m_dot=m_dot, rho=rho + eps_rho, mu=mu, F0=F0, F1=F1).dP
    dp_m = cb._core.sharp_area_change(m_dot=m_dot, rho=rho - eps_rho, mu=mu, F0=F0, F1=F1).dP
    fd_drho = (dp_p - dp_m) / (2 * eps_rho)
    np.testing.assert_allclose(res.dS_drho, fd_drho, rtol=1e-4, err_msg=f"{label} dS_drho")


# ---------------------------------------------------------------------------
# Test 5: Zero-crossing — sign, monotonicity, and smoothness
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "F0, F1, label",
    [
        (0.01, 0.04, "expansion"),
        (0.04, 0.01, "contraction"),
    ],
)
def test_zero_crossing_solver_properties(F0: float, F1: float, label: str) -> None:
    """Sweep m_dot through zero and verify properties the solver relies on.

    1. sign(dP) == sign(m_dot)  — loss always opposes flow
    2. dS_dm > 0 everywhere     — Jacobian has correct sign for convergence
    3. dP, dS_dm are continuous — no jumps in the transition region
    """
    rho, mu, m_scale = 1.2, 1.8e-5, 1e-4
    ms = np.linspace(-5 * m_scale, 5 * m_scale, 201)

    dPs = np.empty_like(ms)
    dS_dms = np.empty_like(ms)
    for i, m in enumerate(ms):
        r = cb._core.sharp_area_change(m, rho, mu, F0, F1, m_scale=m_scale)
        dPs[i] = r.dP
        dS_dms[i] = r.dS_dm

    # All outputs must be finite
    assert np.all(np.isfinite(dPs)), f"{label}: NaN/inf in dP"
    assert np.all(np.isfinite(dS_dms)), f"{label}: NaN/inf in dS_dm"

    # 1. Sign consistency: dP and m_dot must have same sign (or both zero)
    sign_ok = ms * dPs >= -1e-15
    assert np.all(sign_ok), f"{label}: sign(dP) != sign(m_dot) at m_dot = {ms[~sign_ok]}"

    # 2. Monotonicity: dS_dm > 0 everywhere (Jacobian always positive)
    assert np.all(dS_dms >= 0), f"{label}: dS_dm < 0 at m_dot = {ms[dS_dms < 0]}"

    # 3. Smoothness: no jumps relative to the total range of each quantity
    max_dP_jump = np.max(np.abs(np.diff(dPs)))
    max_dS_jump = np.max(np.abs(np.diff(dS_dms)))
    dP_range = np.ptp(dPs) if np.ptp(dPs) > 0 else 1.0
    dS_range = np.ptp(dS_dms) if np.ptp(dS_dms) > 0 else 1.0
    # Adjacent-point jump should be < 5% of total range (201 points → ~0.5% expected)
    assert max_dP_jump < 0.05 * dP_range, (
        f"{label}: discontinuity in dP, max jump = {max_dP_jump:.2e} vs range = {dP_range:.2e}"
    )
    assert max_dS_jump < 0.05 * dS_range, (
        f"{label}: discontinuity in dS_dm, max jump = {max_dS_jump:.2e} vs range = {dS_range:.2e}"
    )


# ---------------------------------------------------------------------------
# Test 6: Equal areas — no loss
# ---------------------------------------------------------------------------
def test_equal_areas_small_loss() -> None:
    """When F0 == F1, ar=1, zeta_con=0 and zeta_exp = alpha-1 ~ small."""
    res = cb._core.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.02, F1=0.02)
    res_exp = cb._core.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04)
    # Equal-area loss should be much smaller than a 1:4 expansion
    assert res.dP < 0.1 * res_exp.dP, (
        f"Equal-area loss ({res.dP:.4f}) should be << expansion loss ({res_exp.dP:.4f})"
    )


# ---------------------------------------------------------------------------
# Test 7: Compressibility correction increases loss
# ---------------------------------------------------------------------------
def test_mach_correction() -> None:
    """k1 = 1 + 0.35*M^2 should increase loss at nonzero Mach."""
    res_m0 = cb._core.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, Mach=0.0)
    res_m5 = cb._core.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, Mach=0.5)
    expected_ratio = 1.0 + 0.35 * 0.5**2  # = 1.0875
    np.testing.assert_allclose(res_m5.dP / res_m0.dP, expected_ratio, rtol=1e-10)

    # mach_clamped should be False for M < MACH_CLAMP (1.5)
    assert not res_m5.mach_clamped


def test_mach_clamp_flag() -> None:
    """Mach numbers above MACH_CLAMP=1.5 are clamped and flagged."""
    res_ok = cb._core.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, Mach=1.0)
    res_clamped = cb._core.sharp_area_change(
        m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, Mach=2.0
    )
    assert not res_ok.mach_clamped
    assert res_clamped.mach_clamped
    # dP at M=2.0 should equal dP at M=1.5 (clamped)
    res_at_clamp = cb._core.sharp_area_change(
        m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, Mach=1.5
    )
    np.testing.assert_allclose(res_clamped.dP, res_at_clamp.dP, rtol=1e-12)


# ---------------------------------------------------------------------------
# Test 8: Network-level wrapper (solver_interface)
# ---------------------------------------------------------------------------
def test_network_level_wrapper() -> None:
    """area_change_residuals_and_jacobian should produce consistent results."""
    Y = cb.species.dry_air()
    m_dot = 0.05
    T_up = 300.0
    P_static_up = 101325.0
    F0 = 0.01
    F1 = 0.04

    res = cb._core.area_change_residuals_and_jacobian(
        m_dot=m_dot,
        P_total_up=P_static_up,
        P_static_up=P_static_up,
        T_up=T_up,
        Y_up=Y,
        P_static_down=P_static_up,
        F0=F0,
        F1=F1,
    )

    # Basic sanity: positive dP for expansion
    assert res.dP_calc > 0
    # d_dP_d_mdot should be positive (more flow -> more loss)
    assert res.d_dP_d_mdot > 0
    # Y derivatives should have correct length
    assert len(res.d_dP_dY_up) == len(Y)

    # FD check on d_dP_d_mdot
    eps = 1e-7
    res_p = cb._core.area_change_residuals_and_jacobian(
        m_dot=m_dot + eps,
        P_total_up=P_static_up,
        P_static_up=P_static_up,
        T_up=T_up,
        Y_up=Y,
        P_static_down=P_static_up,
        F0=F0,
        F1=F1,
    )
    fd_dm = (res_p.dP_calc - res.dP_calc) / eps
    np.testing.assert_allclose(res.d_dP_d_mdot, fd_dm, rtol=1e-3)


# ===========================================================================
# Conical area change tests
# ===========================================================================


# ---------------------------------------------------------------------------
# Conical: diffuser loss > nozzle loss at same geometry
# ---------------------------------------------------------------------------
def test_conical_diffuser_gt_nozzle() -> None:
    """For conical, expansion (diffuser) loss should exceed contraction (nozzle)."""
    L = 0.1
    res_diff = cb._core.conical_area_change(
        m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, length=L
    )
    res_nozz = cb._core.conical_area_change(
        m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.04, F1=0.01, length=L
    )
    assert res_diff.dP > 0
    assert res_nozz.dP > 0
    assert res_diff.dP > res_nozz.dP


# ---------------------------------------------------------------------------
# Conical: length -> 0 recovers sharp-edge behaviour
# ---------------------------------------------------------------------------
def test_conical_short_length_matches_sharp() -> None:
    """At very short length (theta -> pi/2), conical should approach sharp-edge."""
    m_dot, rho, mu = 0.1, 1.2, 1.8e-5
    F0, F1 = 0.01, 0.04
    res_sharp = cb._core.sharp_area_change(m_dot=m_dot, rho=rho, mu=mu, F0=F0, F1=F1)
    res_con = cb._core.conical_area_change(m_dot=m_dot, rho=rho, mu=mu, F0=F0, F1=F1, length=1e-6)
    # Should be within ~10% of sharp (alpha(Re) != 1 causes small difference)
    np.testing.assert_allclose(res_con.dP, res_sharp.dP, rtol=0.15)


# ---------------------------------------------------------------------------
# Conical: longer section -> lower loss (better diffuser efficiency)
# ---------------------------------------------------------------------------
def test_conical_longer_less_loss() -> None:
    """A longer diffuser has a smaller half-angle and lower loss."""
    m_dot, rho, mu = 0.1, 1.2, 1.8e-5
    F0, F1 = 0.01, 0.04
    res_short = cb._core.conical_area_change(m_dot=m_dot, rho=rho, mu=mu, F0=F0, F1=F1, length=0.05)
    res_long = cb._core.conical_area_change(m_dot=m_dot, rho=rho, mu=mu, F0=F0, F1=F1, length=0.5)
    assert res_long.dP < res_short.dP


# ---------------------------------------------------------------------------
# Conical: dS_dmu is always 0 (no Re dependence)
# ---------------------------------------------------------------------------
def test_conical_dS_dmu_zero() -> None:
    """Conical loss has no viscosity dependence."""
    res = cb._core.conical_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, length=0.1)
    assert res.dS_dmu == 0.0


# ---------------------------------------------------------------------------
# Conical: FD Jacobian accuracy
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "m_dot, F0, F1, label",
    [
        (0.1, 0.01, 0.04, "conical_diffuser_fwd"),
        (0.1, 0.04, 0.01, "conical_nozzle_fwd"),
        (-0.1, 0.01, 0.04, "conical_diffuser_rev"),
        (-0.1, 0.04, 0.01, "conical_nozzle_rev"),
    ],
)
def test_conical_jacobian_fd(m_dot: float, F0: float, F1: float, label: str) -> None:
    """Central FD validation for conical dS_dm and dS_drho."""
    rho, mu, L = 1.2, 1.8e-5, 0.1

    res = cb._core.conical_area_change(m_dot=m_dot, rho=rho, mu=mu, F0=F0, F1=F1, length=L)

    eps_m = 1e-7
    dp_p = cb._core.conical_area_change(
        m_dot=m_dot + eps_m, rho=rho, mu=mu, F0=F0, F1=F1, length=L
    ).dP
    dp_m = cb._core.conical_area_change(
        m_dot=m_dot - eps_m, rho=rho, mu=mu, F0=F0, F1=F1, length=L
    ).dP
    fd_dm = (dp_p - dp_m) / (2 * eps_m)
    np.testing.assert_allclose(res.dS_dm, fd_dm, rtol=1e-4, err_msg=f"{label} dS_dm")

    eps_rho = 1e-6
    dp_p = cb._core.conical_area_change(
        m_dot=m_dot, rho=rho + eps_rho, mu=mu, F0=F0, F1=F1, length=L
    ).dP
    dp_m = cb._core.conical_area_change(
        m_dot=m_dot, rho=rho - eps_rho, mu=mu, F0=F0, F1=F1, length=L
    ).dP
    fd_drho = (dp_p - dp_m) / (2 * eps_rho)
    np.testing.assert_allclose(res.dS_drho, fd_drho, rtol=1e-4, err_msg=f"{label} dS_drho")


# ---------------------------------------------------------------------------
# Conical: zero-crossing solver properties
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "F0, F1, label",
    [
        (0.01, 0.04, "conical_expansion"),
        (0.04, 0.01, "conical_contraction"),
    ],
)
def test_conical_zero_crossing(F0: float, F1: float, label: str) -> None:
    """Verify sign consistency, monotonicity, and smoothness through zero."""
    rho, mu, m_scale, L = 1.2, 1.8e-5, 1e-4, 0.1
    ms = np.linspace(-5 * m_scale, 5 * m_scale, 201)
    dPs = np.empty_like(ms)
    dS_dms = np.empty_like(ms)
    for i, m in enumerate(ms):
        r = cb._core.conical_area_change(m, rho, mu, F0, F1, length=L, m_scale=m_scale)
        dPs[i] = r.dP
        dS_dms[i] = r.dS_dm

    assert np.all(np.isfinite(dPs)), f"{label}: NaN/inf in dP"
    assert np.all(np.isfinite(dS_dms)), f"{label}: NaN/inf in dS_dm"
    sign_ok = ms * dPs >= -1e-15
    assert np.all(sign_ok), f"{label}: sign mismatch at m_dot = {ms[~sign_ok]}"
    assert np.all(dS_dms >= 0), f"{label}: dS_dm < 0 at m_dot = {ms[dS_dms < 0]}"

    max_dP_jump = np.max(np.abs(np.diff(dPs)))
    dP_range = np.ptp(dPs) if np.ptp(dPs) > 0 else 1.0
    assert max_dP_jump < 0.05 * dP_range, f"{label}: dP discontinuity"


# ---------------------------------------------------------------------------
# Conical: network-level wrapper
# ---------------------------------------------------------------------------
def test_conical_network_wrapper() -> None:
    """conical_area_change_residuals_and_jacobian basic sanity + FD check."""
    Y = cb.species.dry_air()
    m_dot, T_up, P = 0.05, 300.0, 101325.0
    F0, F1, L = 0.01, 0.04, 0.1

    res = cb._core.conical_area_change_residuals_and_jacobian(
        m_dot=m_dot,
        P_total_up=P,
        P_static_up=P,
        T_up=T_up,
        Y_up=Y,
        P_static_down=P,
        F0=F0,
        F1=F1,
        length=L,
    )
    assert res.dP_calc > 0
    assert res.d_dP_d_mdot > 0
    assert len(res.d_dP_dY_up) == len(Y)

    eps = 1e-7
    res_p = cb._core.conical_area_change_residuals_and_jacobian(
        m_dot=m_dot + eps,
        P_total_up=P,
        P_static_up=P,
        T_up=T_up,
        Y_up=Y,
        P_static_down=P,
        F0=F0,
        F1=F1,
        length=L,
    )
    fd_dm = (res_p.dP_calc - res.dP_calc) / eps
    np.testing.assert_allclose(res.d_dP_d_mdot, fd_dm, rtol=1e-3)


# ---------------------------------------------------------------------------
# Test 9: Public API Export Check
# ---------------------------------------------------------------------------
def test_area_change_exports() -> None:
    """Verify all area change functions and results are exported in the main package."""
    # Functions
    assert hasattr(cb, "sharp_area_change")
    assert hasattr(cb, "conical_area_change")
    assert hasattr(cb, "area_change_residuals_and_jacobian")
    assert hasattr(cb, "conical_area_change_residuals_and_jacobian")

    # Result Structs
    assert hasattr(cb, "AreaChangeResult")
    assert hasattr(cb, "AreaChangeElementResult")

    # Check they match the core versions
    assert cb.sharp_area_change == cb._core.sharp_area_change
    assert cb.AreaChangeResult == cb._core.AreaChangeResult


# ---------------------------------------------------------------------------
# Test 10: Edge Cases - Zero Area and D_h fallback
# ---------------------------------------------------------------------------
def test_zero_area_safety() -> None:
    """Verify that F=0 doesn't crash and returns zero/bounded results."""
    # F0=0 leads to zero velocity fallback in diagnostics, but core handles it with eps
    res = cb.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.0, F1=0.04)
    assert np.isfinite(res.dP)
    assert res.dP > 0  # Should be very high loss but finite


def test_dh_fallback() -> None:
    """Verify that D_h=0 falls back to circular formula."""
    F0, F1 = 0.01, 0.01
    # Circular D = sqrt(4*0.01/pi) approx 0.1128
    res_auto = cb.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=F0, F1=F1, D_h=0.0)

    D_manual = np.sqrt(4.0 * 0.01 / np.pi)
    res_manual = cb.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=F0, F1=F1, D_h=D_manual)

    np.testing.assert_allclose(res_auto.dP, res_manual.dP, rtol=1e-10)


def test_area_change_validation() -> None:
    """Verify that AreaChangeElement.validate() catches bad geometry."""
    # Good geometry
    cb.network.AreaChangeElement("ac1", "n1", "n2", F0=0.01, F1=0.02).validate()

    # Bad F0
    ac_bad = cb.network.AreaChangeElement("ac1", "n1", "n2", F0=0.0, F1=0.02)
    with pytest.raises(ValueError, match="invalid Upstream Area"):
        ac_bad.validate()

    # Bad Dh
    ac_bad_dh = cb.network.AreaChangeElement("ac1", "n1", "n2", F0=0.01, F1=0.02, D_h=-1.0)
    with pytest.raises(ValueError, match="invalid Hydraulic Diameter"):
        ac_bad_dh.validate()
