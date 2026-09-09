"""The soft-barrier weight has units, so it cannot be a constant.

#302 fixed the barrier's fixed point by raising `soft_penalty_alpha` from 1e7
to 1e11. That was the right direction and the wrong kind of answer. The weight
multiplies a squared mass flow to produce a pressure::

    R_i = (Pt_i - Pt_jct) + alpha * max(0, -e_i * mdot_i)^2

so it carries Pa/(kg/s)^2 and is only meaningful against a particular
network's scales. The barrier's fixed point sits at `slack* = sqrt(dP/alpha)`,
and what matters is that as a FRACTION of the flow, `slack*/m_ref`, which a
fixed alpha holds constant for exactly one network size.

Measured on one junction scaled over five decades with every dimensionless
group held fixed -- same pressure, Mach, angle, area ratio, split and loss
coefficients, only bigger or smaller -- the alpha needed to converge follows
`1/m_ref^2` exactly: 1e14, 1e12, 1e10, 1e8, 1e6, 1e4 for size factors 1e-3
through 1e2. Two decades of alpha per decade of size, with the prediction
`alpha_ref * (m_ref/m)^2` landing on the measured value every time. The 1e11
fallback fails on that junction at a hundredth of its size.

`NetworkSolver` now derives the weight from the network's own reference
pressure and mass flow, `alpha = P_ref / (f * m_ref)^2` with `f` =
`BARRIER_SLACK_FRACTION`, and freezes it for the solve -- so it is a constant
in the residual and the Jacobian is unchanged.

Issue #272. Record in `validation/junction/MPCE_CPP_PORT_DESIGN.md` section 7e.
"""

from __future__ import annotations

import dataclasses
import math
import random
import warnings

import pytest

from combaero.network import NetworkSolver
from combaero.network.mpce_element import (
    BARRIER_SLACK_FRACTION,
    DEFAULT_SOFT_PENALTY_ALPHA,
    MultiPortChamberElement,
    scaled_penalty_alpha,
)
from validation.junction import random_robustness as rr


@pytest.fixture(scope="module")
def base_case():
    rng = random.Random(20260906)
    for _ in range(400):
        case = rr.sample(rng)
        if rr.has_root(case) is not True:
            continue
        if (
            case.joining
            and case.drive == "flow_and_pressures"
            and 156.0 < case.theta_deg < 156.5
            and 2.5 < case.psi < 2.6
        ):
            return case
    pytest.skip("the traced case is no longer produced by this seed")


def _converges(case) -> bool:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return rr.classify(rr.build(case)) == "converged"


# ---------------------------------------------------------------------------
# The weight itself
# ---------------------------------------------------------------------------


def test_the_weight_places_the_fixed_point_at_the_intended_fraction():
    """alpha = P_ref / (f*m_ref)^2 must invert to slack* = f*m_ref, taking
    dP = P_ref. Checked by going back through sqrt(dP/alpha) rather than by
    restating the formula."""
    for pressure, mdot in ((1.0e5, 0.1), (6.245e5, 9.2576e-2), (2.0e6, 30.0)):
        alpha = scaled_penalty_alpha(pressure, mdot)
        slack_star = math.sqrt(pressure / alpha)
        assert slack_star == pytest.approx(BARRIER_SLACK_FRACTION * mdot, rel=1e-12)


def test_the_weight_scales_as_one_over_mdot_squared():
    """The measured law: two decades of alpha per decade of size."""
    a1 = scaled_penalty_alpha(1.0e5, 1.0)
    a2 = scaled_penalty_alpha(1.0e5, 0.1)
    a3 = scaled_penalty_alpha(1.0e5, 0.01)

    assert a2 / a1 == pytest.approx(100.0, rel=1e-12)
    assert a3 / a2 == pytest.approx(100.0, rel=1e-12)


@pytest.mark.parametrize(
    ("pressure", "mdot"),
    [(0.0, 1.0), (1.0e5, 0.0), (-1.0e5, 1.0), (1.0e5, -1.0), (float("nan"), 1.0)],
)
def test_a_degenerate_reference_state_falls_back(pressure, mdot):
    """A bad reference must never produce a WEAKER barrier than the fallback,
    which is what a naive divide would do at m_ref -> large or P_ref -> 0."""
    assert scaled_penalty_alpha(pressure, mdot) == DEFAULT_SOFT_PENALTY_ALPHA


# ---------------------------------------------------------------------------
# Which weight an element actually uses
# ---------------------------------------------------------------------------


def _element() -> MultiPortChamberElement:
    return MultiPortChamberElement(
        id="jct",
        inlet_nodes=["a", "b"],
        outlet_nodes=["c"],
        inlet_angles_deg=[0.0, 90.0],
        outlet_angles_deg=[0.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="merge",
        strict=False,
    )


def test_without_a_solver_the_element_uses_the_fallback():
    assert _element().effective_penalty_alpha() == DEFAULT_SOFT_PENALTY_ALPHA


def test_the_solver_supplied_weight_is_used():
    element = _element()
    element._barrier_alpha_scaled = 1.234e13

    assert element.effective_penalty_alpha() == 1.234e13


def test_an_explicit_setting_always_wins():
    """The tuning knob has to keep working: a caller who sets the weight means
    it, and must not have it silently overridden by the network's scales."""
    element = _element()
    element._barrier_alpha_scaled = 1.234e13
    element.soft_penalty_alpha = 5.0e7

    assert element.effective_penalty_alpha() == 5.0e7


def test_the_solver_hands_the_scale_to_every_chamber_element(base_case):
    solver = NetworkSolver(rr.build(base_case))
    elements = [
        e for e in solver.network.elements.values() if isinstance(e, MultiPortChamberElement)
    ]
    assert elements
    assert all(e._barrier_alpha_scaled is None for e in elements)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solver.solve(timeout=30.0)

    for element in elements:
        assert element._barrier_alpha_scaled is not None
        assert element.effective_penalty_alpha() == element._barrier_alpha_scaled


def test_the_weight_is_frozen_for_the_solve_not_recomputed_per_iterate(base_case):
    """If it tracked the state it would enter the Jacobian, which is the one
    thing that would make this a residual-form change."""
    solver = NetworkSolver(rr.build(base_case))
    element = next(
        e for e in solver.network.elements.values() if isinstance(e, MultiPortChamberElement)
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solver.solve(timeout=30.0)
    first = element._barrier_alpha_scaled

    import numpy as np

    solver._residuals(np.array(solver._last_x) * 1.5)
    assert element._barrier_alpha_scaled == first


# ---------------------------------------------------------------------------
# The measurement that motivated it
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("factor", [1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0, 1.0e1, 1.0e2, 1.0e3])
def test_the_same_junction_converges_at_every_size(base_case, factor):
    """Seven decades of size. Scaling the port area at fixed pressure, Mach,
    angle, area ratio, split and loss coefficients leaves every dimensionless
    group untouched: it is the same junction, so it must still solve."""
    assert _converges(dataclasses.replace(base_case, area=base_case.area * factor))


@pytest.mark.parametrize("factor", [1.0e-4, 1.0e-3, 1.0e-2])
def test_the_fixed_fallback_would_have_failed_at_those_sizes(base_case, factor, monkeypatch):
    """The control. Without the hand-off the element falls back to the fixed
    weight, and these are the sizes it cannot solve -- which is the whole
    reason the weight is derived rather than declared. If this ever goes green
    the motivating measurement has changed and the constant should be
    re-derived, not the test relaxed."""
    monkeypatch.setattr(NetworkSolver, "_apply_barrier_scale", lambda self: None)

    assert not _converges(dataclasses.replace(base_case, area=base_case.area * factor))
