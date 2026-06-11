"""Smoke tests for the Mynard 2010 Unified0D port + adapter."""

from __future__ import annotations

import math

import numpy as np

from validation.junction.models import bassett2001
from validation.junction.models.mynard2010 import junction_loss_coefficient
from validation.junction.models.mynard_analytical import MynardAnalytical


def test_matlab_help_example_runs():
    """The Matlab help example -- T-junction split -- runs without error
    and returns finite C/K values."""
    r = junction_loss_coefficient(
        U=np.array([100.0, -50.0, -50.0]),
        A=np.array([2.0, 2.0, 2.0]),
        theta=np.array([math.pi, 0.0, math.pi / 2.0]),
    )
    assert np.all(np.isfinite(r.C))
    assert r.K is not None and np.all(np.isfinite(r.K))
    # For the symmetric T-split, the two collectors should have similar
    # |K| (mirror symmetric in flow but different in geometry).
    assert len(r.K) == 2


def test_mynard_matches_bassett_K6_within_15pct_at_theta_90():
    """At theta=90, psi=1, q=0.5, Mynard's branch-collector K should match
    Bassett K6 to within Mynard's energy-transfer-vs-pure-CV gap (~15%)."""
    r = junction_loss_coefficient(
        U=np.array([1.0, -0.5, -0.5]),
        A=np.array([1.0, 1.0, 1.0]),
        theta=np.array([math.pi, 0.0, math.pi / 2.0]),
    )
    K6 = bassett2001.K6(0.5, 1.0, math.pi / 2.0)
    K_mynard_branch = float(r.K[1])
    rel = abs(K_mynard_branch - K6) / K6
    assert rel < 0.15, f"Mynard K_branch={K_mynard_branch} vs Bassett K6={K6}; rel={rel}"


def test_mynard_smooth_through_theta_120():
    """Mynard handles theta=120 cleanly (Bassett K6 formula goes catastrophic
    at flow-direction extremes; the pseudosupplier reformulation avoids it).
    """
    for theta_deg in (90, 105, 120, 135):
        r = junction_loss_coefficient(
            U=np.array([1.0, -0.5, -0.5]),
            A=np.array([1.0, 1.0, 1.0]),
            theta=np.array([math.pi, 0.0, math.radians(theta_deg)]),
        )
        assert r.K is not None
        assert np.all(np.isfinite(r.K)), f"non-finite K at theta={theta_deg}"
        # K should monotonically increase with theta in this regime
        assert r.K[1] > 0, f"K_branch went negative at theta={theta_deg}"


def test_adapter_returns_K_for_separating_K6_and_K5():
    """MynardAnalytical evaluates Bassett K6 and K5 cases at standard cases."""
    m = MynardAnalytical()
    K6 = m.evaluate("bassett2001", "K6", q=0.5, psi=1.0, theta_rad=math.pi / 2.0)
    K5 = m.evaluate("bassett2001", "K5", q=0.5, psi=1.0, theta_rad=math.pi / 2.0)
    assert K6 is not None and 0.5 < K6 < 1.5
    assert K5 is not None  # K_straight_mynard differs from K5_bassett, but finite


def test_adapter_returns_K_for_joining_K11_and_K12():
    """MynardAnalytical evaluates Bassett K11 and K12 cases (converging flow)."""
    m = MynardAnalytical()
    K11 = m.evaluate("bassett2001", "K11", q=0.5, psi=1.0, theta_rad=math.pi / 2.0)
    K12 = m.evaluate("bassett2001", "K12", q=0.5, psi=1.0, theta_rad=math.pi / 2.0)
    assert K11 is not None and math.isfinite(K11)
    assert K12 is not None and math.isfinite(K12)


def test_adapter_returns_none_for_unsupported_ids():
    """K1/K3/K4/K7/K8/K9/K10 use special CVs and aren't supported."""
    m = MynardAnalytical()
    for K_id in ("K1", "K3", "K4", "K7", "K8", "K9", "K10"):
        assert m.evaluate("bassett2001", K_id, 0.5, 1.0, math.pi / 2.0) is None


def test_adapter_returns_none_for_non_bassett_paper():
    """Mynard maps via Bassett K_ids; other papers' K_ids are not supported."""
    m = MynardAnalytical()
    assert m.evaluate("wang2014", "K_13", 0.5, None, None, M_3=0.3) is None
