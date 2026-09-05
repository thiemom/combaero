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


# The N > 3 consequence (a silently all-zero loss Jacobian) is no longer
# reachable: MPCEv2Element refuses more than three ports at construction, and
# that refusal is pinned in test_mpce_v2_degenerate_iterates.py. The closure-
# level precondition above (K is None beyond three ports) is what makes the
# refusal necessary and stays.
