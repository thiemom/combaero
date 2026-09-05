"""Bassett 2001 Fig 7b -- his own measured case -- through the network solver.

theta = 45 deg, psi = 3 (lateral one third of the main area), separating
flow. Bassett measured this curve (Fig 7b); the K6 formula (Eq 27) fits it.
Every point on it sits at M ~ 0.03, deep inside the incompressible regime the
closure assumes -- so a failure here is numerical, not a physics-envelope
limit, and a case the literature documents is the highest-priority thing to
keep passing (issue #271).

Measured 2026-09-05 with the production-faithful adapter (strict=False):

    topology     q=0.2      q=0.4      q=0.6      q=0.8
    imposed_q    ok         ok         ok         ok
    three_pb     FAIL       ok 0.444   ok 1.246   ok 2.765     (K6: 0.444 1.247 2.769)
    mfb_two_pb   FAIL       FAIL       ok 1.403   WRONG ROOT: converged, K 3.44 vs 2.77

Under strict=True the three_pb q=0.4/0.6 cases raised on a transient
wrong-sign iterate. mfb_two_pb q=0.8 converges -- under BOTH strict settings,
run cold -- to a mirror root 24% off Bassett, reported as a success: the port
signs are right, so the direction-only verifier cannot see it. The tests
below pin what passes as passing, and what does not as strict xfails --
targets, not tolerated failures.
"""

from __future__ import annotations

import math

import pytest

from validation.junction.models import bassett2001
from validation.junction.models.mpce_v2_network import MPCEv2Network

_THETA = math.radians(45.0)
_PSI = 3.0


@pytest.fixture
def model():
    # Fresh per test: the same case has been seen to land differently after
    # other solves on one instance, and a target must be measured cold.
    return MPCEv2Network(strict=False)


def _run(model, topology: str, q: float):
    return model.evaluate_network("bassett2001", "K6", q, _PSI, _THETA, topology=topology)


@pytest.mark.parametrize("q", [0.2, 0.4, 0.6, 0.8])
def test_imposed_q_converges_across_the_curve(model, q):
    r = _run(model, "imposed_q", q)
    assert r.converged, r.message


@pytest.mark.parametrize("q", [0.4, 0.6, 0.8])
def test_three_pb_converges_and_lands_on_bassett_K6(model, q):
    """The case that raised under strict=True. With the soft barrier steering
    the reversed first iterate back, it converges to Bassett's own curve."""
    r = _run(model, "three_pb", q)

    assert r.converged, r.message
    assert r.K_lateral == pytest.approx(bassett2001.K6(q, _PSI, _THETA), abs=0.05)


@pytest.mark.xfail(
    strict=False,
    reason=(
        "mfb_two_pb q=0.8 is NOT reliably solved, and the outcome depends on the "
        "process: run cold in one process it converges to K_lat 3.44 against "
        "Bassett's 2.77 and is reported as a success (the port signs are correct, "
        "so the direction-only verifier cannot see the mirror root); under the "
        "parallel suite the same case fails or is demoted instead. Both outcomes "
        "are wrong, in different ways, so a strict marker would flap. This is the "
        "one deliberately non-strict xfail in the junction tests, and the "
        "non-determinism is itself a recorded defect (issue #271). It comes off "
        "when the case converges to Bassett's curve every time."
    ),
)
def test_mfb_two_pb_mirror_root_is_caught(model):
    r = _run(model, "mfb_two_pb", 0.8)

    assert (not r.converged) or r.K_lateral == pytest.approx(
        bassett2001.K6(0.8, _PSI, _THETA), abs=0.15
    ), (
        f"reported converged at K_lat={r.K_lateral:.3f} vs Bassett {bassett2001.K6(0.8, _PSI, _THETA):.3f}"
    )


# ---------------------------------------------------------------------------
# Targets: a documented case must not stay failing quietly
# ---------------------------------------------------------------------------

_TARGET = (
    "Bassett Fig 7b at a low lateral fraction does not converge through the "
    "pressure-driven topologies: 'iteration is not making good progress'. "
    "It sits inside the paper's measured conditions (M ~ 0.03), so this is a "
    "solver-path failure, not physics. Priority target for the step-3 port "
    "and the degenerate-iterate work (issue #271). strict=True so the moment "
    "it starts passing, this marker must come off."
)


@pytest.mark.xfail(strict=True, reason=_TARGET)
def test_three_pb_low_lateral_fraction_converges(model):
    r = _run(model, "three_pb", 0.2)
    assert r.converged, r.message


@pytest.mark.xfail(strict=True, reason=_TARGET)
@pytest.mark.parametrize("q", [0.2, 0.4])
def test_mfb_two_pb_low_lateral_fraction_converges(model, q):
    r = _run(model, "mfb_two_pb", q)
    assert r.converged, r.message
