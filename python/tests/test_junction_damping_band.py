"""FLOW_RATIO_DAMPING: the regulariser's band, and why it is cancelled.

``damping = 1 - exp(-FlowRatio / tau)`` multiplies Mynard's collector
coefficient C. The Matlab reference justifies it as "avoids infinite C when
FlowRatio approaches zero". In the K form MPCEv2 uses (Mynard Eq 18) that
divergence is multiplied by FlowRatio^2, so K reaches the dead-branch limit
with or without the knob. These tests pin that -- exactly, because it is a
limit, not a tuning -- and pin the band over which the knob is inert.

Direct closure calls, no network: the scorecard evidence lives on #271.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from combaero.network import _mynard2010 as mynard

_A = np.array([0.01, 0.01, 0.01])
_THETA = np.array([math.pi, 0.0, math.pi / 2.0])  # common, straight, lateral


def _dividing(fr_lateral: float) -> np.ndarray:
    """Supplier at port 0; lateral collector carries fraction fr of it."""
    return np.array([10.0, -10.0 * (1.0 - fr_lateral), -10.0 * fr_lateral])


@pytest.fixture
def damping(monkeypatch):
    def _set(tau: float):
        monkeypatch.setattr(mynard, "FLOW_RATIO_DAMPING", tau)

    return _set


@pytest.mark.parametrize("tau", [1e-9, 0.005, 0.01, 0.02])
def test_band_is_inert_in_the_interior(damping, tau):
    """Over (0, 0.02] the knob changes nothing where the data lives."""
    U = _dividing(0.5)
    damping(0.02)
    reference = mynard.junction_loss_coefficient(U, _A, _THETA).K
    damping(tau)
    K = mynard.junction_loss_coefficient(U, _A, _THETA).K

    assert reference is not None and K is not None
    np.testing.assert_allclose(K, reference, rtol=1e-9)


def test_above_the_band_it_distorts(damping):
    """The label's upper bound is real: tau = 0.2 moves an interior point."""
    U = _dividing(0.5)
    damping(0.02)
    reference = mynard.junction_loss_coefficient(U, _A, _THETA).K
    damping(0.2)
    K = mynard.junction_loss_coefficient(U, _A, _THETA).K

    assert reference is not None and K is not None
    assert not np.allclose(K, reference, rtol=1e-3)


@pytest.mark.parametrize("fr", [1e-3, 1e-4, 1e-6])
def test_K_reaches_the_dead_branch_limit_with_or_without_damping(damping, fr):
    """Eq 18 multiplies C by FlowRatio^2, cancelling the 1/FlowRatio blow-up.

    K_lat -> 1 (all of the common dynamic head lost into a dead branch) is the
    physical limit, and it is reached identically with the knob on or off --
    the knob protects a quantity the K form never exposes.
    """
    U = _dividing(fr)
    with np.errstate(all="ignore"):
        damping(0.02)
        K_on = mynard.junction_loss_coefficient(U, _A, _THETA).K
        damping(1e-300)
        K_off = mynard.junction_loss_coefficient(U, _A, _THETA).K

    assert K_on is not None and K_off is not None
    assert K_on[1] == pytest.approx(1.0, abs=1e-3)
    assert K_off[1] == pytest.approx(1.0, abs=1e-3)


def test_damping_is_what_bounds_C_and_nothing_else(damping):
    """The knob's one real effect: C stays bounded as FlowRatio -> 0.

    This is the quantity the Matlab comment is about. It would matter for
    Mynard's native C-form residual; it is why the knob is retained.
    """
    with np.errstate(all="ignore"):
        damping(0.02)
        C_on = [
            abs(float(mynard.junction_loss_coefficient(_dividing(fr), _A, _THETA).C[2]))
            for fr in (1e-2, 1e-4, 1e-6)
        ]
        damping(1e-300)
        C_off = [
            abs(float(mynard.junction_loss_coefficient(_dividing(fr), _A, _THETA).C[2]))
            for fr in (1e-2, 1e-4, 1e-6)
        ]

    assert max(C_on) < 20.0, "damping no longer bounds C"
    assert C_off[-1] > 1e4, "C no longer diverges without damping -- the knob has lost its reason"


def test_exact_zero_flow_edge_is_not_guarded_by_damping(damping):
    """Documented edge: at FlowRatio exactly 0, K_lat = -1 either way.

    A sign-flipped dead branch. Pinned so that whoever fixes it as part of
    the degenerate-iterate work (#271 step 2) sees this test fail and removes
    it rather than discovering the edge again.
    """
    with np.errstate(all="ignore"):
        damping(0.02)
        K_on = mynard.junction_loss_coefficient(_dividing(0.0), _A, _THETA).K
        damping(1e-300)
        K_off = mynard.junction_loss_coefficient(_dividing(0.0), _A, _THETA).K

    assert K_on is not None and K_off is not None
    assert K_on[1] == pytest.approx(-1.0)
    assert K_off[1] == pytest.approx(-1.0)
