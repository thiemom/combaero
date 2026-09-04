"""Coverage for the MPCEv2 FD-fallback guards, which nothing exercised.

`MPCEv2Element` falls back to finite differences whenever the sympy branch
declines, and that loop carries three guards that `continue` past a perturbed
evaluation -- leaving the corresponding Jacobian column at its initialised
ZERO. A zero column is a wrong derivative, not an error: indistinguishable at
the call site from "this unknown has no influence".

An audit (issue #271) found only the sign-flip guard firing across the
validation dataset and the whole test suite, which is *not* evidence the other
two are unreachable -- that dataset holds only 3-port junctions with
well-behaved iterates. These tests pin the preconditions that make each guard
reachable, so the guards surface if the closure beneath them changes.

They deliberately test `_mynard2010`'s contract rather than the guard's effect,
because the contract is what determines reachability.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb
from combaero.network._mynard2010 import junction_loss_coefficient
from combaero.network.components import NetworkMixtureState
from combaero.network.mpce_v2_element import MPCEv2Element

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_A3 = np.array([0.01, 0.01, 0.01])
_THETA3 = np.array([math.pi, 0.0, math.pi / 2.0])


# ---------------------------------------------------------------------------
# Preconditions that make each guard reachable
# ---------------------------------------------------------------------------


def test_mynard_supplies_K_for_three_ports():
    """The supported case: 3 ports, K present and one entry per non-common port."""
    for label, U in (
        ("separating", np.array([-10.0, 6.0, 4.0])),
        ("joining", np.array([-6.0, -4.0, 10.0])),
    ):
        result = junction_loss_coefficient(U, _A3, _THETA3)
        assert result.K is not None, f"{label}: K must be defined for a 3-port junction"
        assert len(result.K) == 2, f"{label}: expected one K per non-common port"


def test_mynard_returns_no_K_beyond_three_ports():
    """Mynard's closure is 3-branch by derivation, so K is None for N > 3.

    This is what makes MPCEv2's `K is None` guard reachable. If the closure ever
    gains N > 3 support, this test fails and the guard -- and the zero-Jacobian
    consequence pinned below -- must be revisited.
    """
    U = np.array([-10.0, 4.0, 3.0, 3.0])
    A = np.array([0.01] * 4)
    theta = np.array([math.pi, 0.0, math.pi / 2.0, math.pi / 4.0])

    result = junction_loss_coefficient(U, A, theta)

    assert result.K is None, "Mynard gained N > 3 support -- revisit the K-shape guard"


@pytest.mark.parametrize(
    "label, U",
    [
        ("every port a supplier", np.array([-4.0, -3.0, -3.0])),
        ("every port a collector", np.array([4.0, 3.0, 3.0])),
    ],
)
def test_mynard_raises_on_a_degenerate_flow_split(label: str, U: np.ndarray):
    """Mynard indexes the collector and supplier masks without checking them.

    Newton can pass through an iterate where every port flows the same way, so
    the `except Exception` guard is reachable. Pinned so that a future closure
    which returns instead of raising surfaces here rather than silently
    changing which branch the Jacobian takes.
    """
    with pytest.raises((IndexError, ValueError)):
        junction_loss_coefficient(U, _A3, _THETA3)


def test_perturbing_a_near_zero_port_flips_its_velocity_sign():
    """The precondition for the sign-flip guard, which the audit saw fire.

    The FD loop perturbs by max(|mdot| * 1e-4, 1e-7); for a port sitting near
    zero flow that step crosses zero, and Mynard's supplier/collector
    classification changes underneath the difference.
    """
    port_mdots = np.array([-0.10, 0.10, -1e-9])
    rho = np.full(3, 1.16)

    U = -port_mdots / (rho * _A3)
    eps = max(abs(port_mdots[2]) * 1e-4, 1e-7)
    perturbed = port_mdots.copy()
    perturbed[2] += eps
    U_pert = -perturbed / (rho * _A3)

    assert (np.sign(U_pert) != np.sign(U)).any(), (
        "the perturbation no longer crosses zero -- the sign-flip guard's "
        "trigger condition has changed"
    )


# ---------------------------------------------------------------------------
# The consequence, pinned as a known defect so a fix surfaces
# ---------------------------------------------------------------------------


def _four_port_element() -> MPCEv2Element:
    element = MPCEv2Element.__new__(MPCEv2Element)
    element.id = "jct"
    element.N = 4
    element.port_nodes = [f"p{i}" for i in range(4)]
    element.port_areas = [0.01] * 4
    element.port_angles_deg = [180.0, 0.0, 90.0, 45.0]
    element._port_signs = [-1.0, 1.0, 1.0, 1.0]
    element._port_element_ids = [f"e{i}" for i in range(4)]
    element.flow_direction = "branch"
    element.strict = False
    element.joining_etransfer_alpha = 0.2
    element.jacobian_method = "sympy"
    element.penalty_alpha = 0.0
    element.eta_scale = 1.0  # faithful-port default; see test_junction_tuned_constants
    return element


def _four_port_states():
    return [
        NetworkMixtureState(P=1.0e5, Pt=100_300.0, T=300.0, Tt=300.5, m_dot=0.0, Y=_Y)
        for _ in range(4)
    ]


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Known defect (issue #271): MPCEv2Element does not restrict N, and "
        "Mynard supplies no K beyond 3 ports, so every FD column is skipped "
        "and the whole dKQ block is silently zero. Either the element must "
        "reject N > 3 at construction -- honest, the closure is 3-branch by "
        "derivation -- or the closure must be extended. Returning a zero "
        "Jacobian with finite residuals is not a third option. Remove this "
        "xfail with whichever fix lands."
    ),
)
def test_four_port_junction_has_a_nonzero_loss_jacobian():
    element = _four_port_element()

    _, jac = element.residuals(_four_port_states(), 100_300.0, [-0.12, 0.05, 0.04, 0.03])

    mdot_entries = sum(1 for row in range(4) for key in jac[row] if key.endswith(".m_dot"))
    assert mdot_entries > 0, (
        "the entire dKQ block is zero: every mass-flow derivative of the loss "
        "term is missing, with no error and no warning"
    )


def test_four_port_junction_still_returns_finite_residuals():
    """Why the defect above is silent rather than loud: nothing looks wrong.

    The residuals are perfectly well-formed; only the Jacobian is empty, so a
    solver has no signal that anything is amiss.
    """
    element = _four_port_element()

    residuals, _ = element.residuals(_four_port_states(), 100_300.0, [-0.12, 0.05, 0.04, 0.03])

    assert len(residuals) == 5  # 4 impulse rows + mass row
    assert all(math.isfinite(r) for r in residuals)
