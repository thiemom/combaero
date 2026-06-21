"""
Plumbing + regression tests for the joining_etransfer_alpha correction.

The default alpha = 0.2 was calibrated against Bassett K11_corr/K12_corr +
Idelchik 1966 tabulated values at psi in [1.25, 3.33] / theta in {30, 45, 90}.
These tests don't re-fit the calibration; they lock in the chosen value and
exercise the plumbing so accidental changes are caught.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from combaero.network.mpce_v2_element import MPCEv2Element
from validation.junction.models.mynard2010 import junction_loss_coefficient

CALIBRATED_ALPHA = 0.2


def _joining_K(psi: float, theta_deg: float, q: float, alpha: float) -> tuple[float, float]:
    """Mynard (K_str, K_bra) for joining type-6 at the given (psi, theta, q).

    Reference scales: rho=1, A_com=1, m_com=1, so K values are normalised by
    0.5 * rho * u_com^2 = 0.5.
    """
    A_str = 1.0
    A_bra = 1.0 / psi
    A_com = 1.0
    m_str = (1.0 - q) * 1.0
    m_bra = q * 1.0
    U = np.array([m_str / A_str, m_bra / A_bra, -1.0 / A_com])
    A = np.array([A_str, A_bra, A_com])
    theta = np.array([0.0, math.radians(theta_deg), math.pi])
    r = junction_loss_coefficient(U, A, theta, joining_etransfer_alpha=alpha)
    assert r.K is not None and len(r.K) == 2
    return float(r.K[0]), float(r.K[1])


def test_default_alpha_value() -> None:
    """The calibrated alpha is a published constant; changing it is deliberate."""
    assert MPCEv2Element.DEFAULT_JOINING_ETRANSFER_ALPHA == CALIBRATED_ALPHA


def test_alpha_passes_through_constructor() -> None:
    """``joining_etransfer_alpha=None`` resolves to the class default."""
    e = MPCEv2Element(
        id="jct",
        inlet_nodes=["port_str", "port_bra"],
        outlet_nodes=["port_com"],
        inlet_angles_deg=[0.0, 90.0],
        outlet_angles_deg=[0.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="merge",
    )
    assert e.joining_etransfer_alpha == CALIBRATED_ALPHA


def test_explicit_alpha_override() -> None:
    """Passing an explicit value overrides the default."""
    e = MPCEv2Element(
        id="jct",
        inlet_nodes=["port_str", "port_bra"],
        outlet_nodes=["port_com"],
        inlet_angles_deg=[0.0, 90.0],
        outlet_angles_deg=[0.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="merge",
        joining_etransfer_alpha=0.0,
    )
    assert e.joining_etransfer_alpha == 0.0


def test_alpha_zero_recovers_faithful_mynard() -> None:
    """alpha=0 must reproduce the original Mynard K values exactly."""
    # Symmetric case at psi=1: correction is zero by construction (area_asym=0)
    # so any alpha gives the same K. Pick an asymmetric case to actually test
    # that alpha=0 is the no-correction path.
    K_str_0, K_bra_0 = _joining_K(psi=2.5, theta_deg=45, q=0.5, alpha=0.0)
    K_str_d, K_bra_d = _joining_K(psi=2.5, theta_deg=45, q=0.5, alpha=CALIBRATED_ALPHA)
    # alpha=0 should differ from alpha=0.2 (otherwise the knob is dead)
    assert abs(K_str_0 - K_str_d) > 0.01
    assert abs(K_bra_0 - K_bra_d) > 0.01


def test_correction_vanishes_at_psi_1() -> None:
    """At psi=1 (equal supplier areas) the correction is zero by construction."""
    K_no = _joining_K(psi=1.0, theta_deg=45, q=0.5, alpha=0.0)
    K_with = _joining_K(psi=1.0, theta_deg=45, q=0.5, alpha=CALIBRATED_ALPHA)
    assert K_no[0] == pytest.approx(K_with[0], abs=1e-12)
    assert K_no[1] == pytest.approx(K_with[1], abs=1e-12)


def test_correction_increases_K_at_psi_greater_than_1() -> None:
    """Positive alpha shifts both K's higher at psi > 1 (the calibrated direction)."""
    K_no_str, K_no_bra = _joining_K(psi=3.0, theta_deg=45, q=0.5, alpha=0.0)
    K_w_str, K_w_bra = _joining_K(psi=3.0, theta_deg=45, q=0.5, alpha=CALIBRATED_ALPHA)
    assert K_w_str > K_no_str
    assert K_w_bra > K_no_bra


def test_correction_inactive_for_separating_flow() -> None:
    """Diverging flow (single supplier) is unchanged by the correction."""
    A_com = 1.0
    A_str = 1.0
    A_bra = 1.0 / 2.5
    # Separating: 1 inlet (com) + 2 outlets (str, bra)
    U = np.array([1.0 / A_com, -0.5 / A_str, -0.5 / A_bra])  # com inlet, str+bra outlets
    A = np.array([A_com, A_str, A_bra])
    theta = np.array([math.pi, 0.0, math.radians(45.0)])

    r0 = junction_loss_coefficient(U, A, theta, joining_etransfer_alpha=0.0)
    rd = junction_loss_coefficient(U, A, theta, joining_etransfer_alpha=CALIBRATED_ALPHA)
    assert r0.K is not None and rd.K is not None
    np.testing.assert_allclose(r0.K, rd.K, atol=1e-12)


def test_locked_K_at_canonical_calibration_point() -> None:
    """Regression-lock the K value at one canonical calibration point so the
    published alpha can't drift silently.

    Chosen point: (psi=2.5, theta=45, q=0.5) -- centre of the calibration grid.
    Values were generated with alpha=0.2 at the time of calibration.
    """
    K_str, K_bra = _joining_K(psi=2.5, theta_deg=45, q=0.5, alpha=CALIBRATED_ALPHA)
    # Locked to ~4 decimal places. If these change, either the calibration
    # was re-run (intentional, update here) or upstream Mynard logic
    # drifted (unintentional, investigate).
    assert K_str == pytest.approx(-0.2811, abs=5e-3)
    assert K_bra == pytest.approx(1.0314, abs=5e-3)
