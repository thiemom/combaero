"""Bassett 2001 Fig 7b -- his own measured case -- through the network solver.

theta = 45 deg, psi = 3 (lateral one third of the main area), separating
flow. Bassett measured this curve (Fig 7b); the K6 formula (Eq 27) fits it.
Every point on it sits at M ~ 0.03, deep inside the incompressible regime the
closure assumes, so a case the literature documents is the highest-priority
thing to keep passing (issue #271).

What each topology can actually prove (docs/archive/JUNCTION_OPERATING_POINT_271.md):

* ``imposed_q`` imposes both mass flows, so q is a boundary condition, the
  root is unique, and the extracted K is the model's answer at that q. This
  is the only topology whose K value means anything.
* ``three_pb`` imposes every total pressure through lossless connections, and
  ``_extract_K`` normalises by the fixed reference flow -- so the extracted K
  is the imposed target BY CONSTRUCTION. Falsified directly: a model with
  eta_scale=3.0, whose imposed_q K at q=0.8 moves 2.88 -> 3.07, still reports
  2.7677 against a target of 2.7689. Asserting that K here (as the first
  version of this file did) asserts nothing. Convergence only.
* ``mfb_two_pb`` fixes the inlet flow and both outlet pressures, which
  constrains only K_lat - K_str. The split is an outcome: q=0.8 converges at
  q=0.857, and its K matches the model's own imposed_q K there to four
  decimals. Comparing it against Bassett at the TARGET q is a comparison at
  the wrong operating point, so this file does not do that either.

Measured 2026-09-05 with the production-faithful adapter (strict=False):

    topology     q=0.2        q=0.4        q=0.6      q=0.8
    imposed_q    ok           ok           ok         ok
    three_pb     no root      ok           ok         ok
    mfb_two_pb   no root      no root      ok         ok  (at q=0.857)

The q=0.2 and q=0.4 failures are NOT solver-path failures. Those boundary-value
problems have no solution: the model's K_lat - K_str has a minimum of 0.551 for
this geometry and the targets sit at 0.122 and 0.385. They are xfail targets
that come off when the K_straight gap (#272) closes, not when the solver or
the seed improves.
"""

from __future__ import annotations

import math

import pytest

from validation.junction.models.mpce_v2_network import MPCEv2Network

_THETA = math.radians(45.0)
_PSI = 3.0


@pytest.fixture
def model():
    # Fresh per test: a target must be measured cold, not after other solves.
    #
    # These solves used to be PYTHONHASHSEED-dependent -- mfb_two_pb q=0.8
    # landed on a different root in 5 of 10 identical processes -- because
    # solver.py's _propagate_pressure_guess seeded its BFS from a set of
    # node-ID strings. Fixed in the same change as this file; the assertions
    # below hold under hash seeds 0, 1, 7 and 13, and would flap without it.
    return MPCEv2Network(strict=False)


def _run(model, topology: str, q: float):
    return model.evaluate_network("bassett2001", "K6", q, _PSI, _THETA, topology=topology)


# ---------------------------------------------------------------------------
# imposed_q: the only topology whose K value is the model's own answer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", [0.2, 0.4, 0.6, 0.8])
def test_imposed_q_converges_across_the_curve(model, q):
    r = _run(model, "imposed_q", q)
    assert r.converged, r.message


@pytest.mark.parametrize(
    "q, expected",
    [(0.2, 0.4534), (0.4, 0.5816), (0.6, 1.3888), (0.8, 2.8842)],
)
def test_imposed_q_reproduces_the_model_curve(model, q, expected):
    """Pins the model's own Fig 7b curve, so a physics change has to state
    itself here. These are NOT Bassett's values -- the gap to his measured
    curve is the accuracy question, tracked on the scorecard, and it is
    largest exactly where K_straight is worst."""
    r = _run(model, "imposed_q", q)

    assert r.converged, r.message
    assert r.K_lateral == pytest.approx(expected, abs=0.01)


# ---------------------------------------------------------------------------
# Pressure-driven: convergence only, for the reasons in the module docstring
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", [0.4, 0.6, 0.8])
def test_three_pb_converges(model, q):
    """The case that raised under strict=True. With the soft barrier steering
    the reversed first iterate back it converges. No K assertion: in this
    topology the extracted K is the imposed target whatever the model does."""
    assert _run(model, "three_pb", q).converged


@pytest.mark.parametrize("q", [0.6, 0.8])
def test_mfb_two_pb_converges(model, q):
    assert _run(model, "mfb_two_pb", q).converged


def test_mfb_two_pb_root_is_the_models_own_answer(model):
    """q=0.8 settles at q=0.857 and must reproduce the model there.

    This was once recorded as an uncaught mirror root, from comparing K
    against Bassett at the target q. The honest check is against the model's
    own imposed_q curve at the q the solve actually reached.
    """
    r = _run(model, "mfb_two_pb", 0.8)
    assert r.converged, r.message

    reference = _run(model, "imposed_q", 0.8569)
    assert reference.converged, reference.message
    assert r.K_lateral == pytest.approx(reference.K_lateral, abs=0.01)


# ---------------------------------------------------------------------------
# Targets: infeasible until the K_straight gap closes
# ---------------------------------------------------------------------------

_TARGET = (
    "Bassett Fig 7b at a low lateral fraction has NO root in the "
    "pressure-driven topologies: they constrain K_lat - K_str, the model's "
    "value of that has a minimum of 0.551 for this geometry, and the targets "
    "are 0.122 (q=0.2) and 0.385 (q=0.4). 'Not making good progress' is the "
    "solver reporting an infeasible system honestly, and no initial guess "
    "changes it. This comes off when #272 closes the K_straight gap, not "
    "when the solver improves. strict=True so it is noticed either way."
)


@pytest.mark.xfail(strict=True, reason=_TARGET)
def test_three_pb_low_lateral_fraction_converges(model):
    assert _run(model, "three_pb", 0.2).converged


@pytest.mark.xfail(strict=True, reason=_TARGET)
@pytest.mark.parametrize("q", [0.2, 0.4])
def test_mfb_two_pb_low_lateral_fraction_converges(model, q):
    assert _run(model, "mfb_two_pb", q).converged
