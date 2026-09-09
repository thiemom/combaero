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

**The 2026-09-05 reading of this case was wrong, and the correction is the
reason most of it now passes.** That reading recorded no root at q = 0.2 and
0.4, on the grounds that the model's ``K_lat - K_str`` bottoms out at 0.551
while the targets sat at 0.122 and 0.385, and filed them as xfails waiting on
the #272 K_straight gap.

The targets were built wrongly. Bassett Table 1 indexes each coefficient on
the mass-flow fraction in ITS OWN leg -- ``K5`` on ``mdot_A/mdot_C`` (the
straight leg), ``K6`` on ``mdot_B/mdot_C`` (the lateral) -- so at one physical
operating point they are ``K5(1 - q)`` and ``K6(q)``. All three network
adapters read both off the same abscissa, which asks the boundary conditions
for a state made of two different operating points. The correct target at
q = 0.2 is ``K6(0.2) - K5(0.8) = 0.422``, not 0.122, and the model hits it
exactly -- it reproduces both of Bassett's coefficients to three decimals.
See ``bassett2001.separating_pair_at``.

Measured 2026-09-08, after the correction:

    topology     q=0.2   q=0.3   q=0.4   q>=0.5
    imposed_q    ok      ok      ok      ok
    three_pb     ok      ok      ok      solver fails (a root EXISTS)
    mfb_two_pb   ok      ok      ok      ok, and it lands ON the asked q

``mfb_two_pb`` is the one that changed character. It used to drift far from
the operating point it was asked for -- q = 0.8 settling at 0.857 -- and now
reaches 0.7970. The old drift was the solve honestly chasing a target built
from a mismatched pair.

``three_pb`` above q = 0.5 is now the only failure, and it is a SOLVER
failure, not an infeasible system: because the model reproduces Bassett's
coefficients, the corrected target equals the model's own value at the asked
q, so a root exists there by construction. That is a different xfail from the
one it replaces.
"""

from __future__ import annotations

import math

import pytest

from validation.junction.models.mpce_network import MPCENetwork

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
    return MPCENetwork(strict=False)


def _run(model, topology: str, q: float):
    return model.evaluate_network("bassett2001", "K6", q, _PSI, _THETA, topology=topology)


# ---------------------------------------------------------------------------
# imposed_q: the only topology whose K value is the model's own answer
# ---------------------------------------------------------------------------


_INADMISSIBLE = (
    "The model is not DISSIPATIVE at a low lateral fraction for this "
    "geometry: its mass-weighted mean K goes negative below q ~ 0.22 "
    "(psi=3, theta=45), so the junction would create flow work. The "
    "post-solve energy check added for #271 defect 10 now refuses those "
    "states, which is correct -- Bassett's own coefficients keep the "
    "weighted mean at +0.19 or better everywhere. The three_pb cases fail "
    "for the same reason at one remove: that topology wanders down to "
    "q ~ 0.02, inside the inadmissible band. These come off when #272 "
    "closes the K_straight gap."
)


@pytest.mark.parametrize("q", [0.2, 0.4, 0.6, 0.8])
def test_imposed_q_converges_across_the_curve(model, q):
    r = _run(model, "imposed_q", q)
    assert r.converged, r.message


@pytest.mark.parametrize(
    "q, expected",
    [(0.2, 0.3623), (0.4, 0.4448), (0.6, 1.2520), (0.8, 2.7937)],
)
def test_imposed_q_reproduces_the_model_curve(model, q, expected):
    """Pins the model's own Fig 7b curve, so a physics change has to state
    itself here.

    Updated 2026-09-05 with the dividing-streamline recovery. They were
    0.4534 / 0.5816 / 1.3888 / 2.8842 and are now within 1% of Bassett's own
    K6 (0.3622 / 0.4445 / 1.2467 / 2.7689) at every point, because with the
    recovery restored and Mynard's fitted transfer off the closure reproduces
    his analytical pair rather than approximating it.
    """
    r = _run(model, "imposed_q", q)

    assert r.converged, r.message
    assert r.K_lateral == pytest.approx(expected, abs=0.01)


# ---------------------------------------------------------------------------
# Pressure-driven: convergence only, for the reasons in the module docstring
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", [0.2, 0.3, 0.4])
def test_three_pb_converges(model, q):
    """The case that raised under strict=True. With the soft barrier steering
    the reversed first iterate back it converges. No K assertion: in this
    topology the extracted K is the imposed target whatever the model does.

    The low-q end used to be an xfail on the grounds that no root existed --
    see the module docstring; the target was built from a mismatched pair.
    q = 0.4 was separately an xfail for wandering off the operating point.
    Both come off with the corrected pairing.
    """
    assert _run(model, "three_pb", q).converged


@pytest.mark.parametrize("q", [0.6, 0.8])
def test_mfb_two_pb_converges(model, q):
    assert _run(model, "mfb_two_pb", q).converged


def test_mfb_two_pb_root_is_the_models_own_answer(model):
    """The K it reports must be the model's own answer at the q it REACHED.

    Once recorded as an uncaught mirror root, from comparing K against Bassett
    at the target q. The honest check is against the model's own imposed_q
    curve at the achieved operating point, whatever that turns out to be --
    written against ``q_converged`` rather than a hard-coded number, because
    that number moved once already when the target pairing was corrected
    (0.857 -> 0.797).
    """
    r = _run(model, "mfb_two_pb", 0.8)
    assert r.converged, r.message
    assert r.q_converged is not None

    reference = _run(model, "imposed_q", r.q_converged)
    assert reference.converged, reference.message
    assert r.K_lateral == pytest.approx(reference.K_lateral, abs=0.01)


def test_mfb_two_pb_now_lands_on_the_operating_point_it_was_asked_for(model):
    """The correction's most visible effect, and the reason the drift was
    never a modelling problem.

    With the target built from a mismatched pair the split was an outcome that
    wandered -- q = 0.8 settling at 0.857. With the pair taken at one
    operating point it settles within 0.005 of what was asked across the
    range. If this regresses, suspect the pairing before the solver.
    """
    for q in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
        r = _run(model, "mfb_two_pb", q)
        assert r.converged, f"q={q}: {r.message}"
        assert r.q_converged == pytest.approx(q, abs=0.005), f"asked q={q}, reached {r.q_converged}"


# ---------------------------------------------------------------------------
# Targets: infeasible until the K_straight gap closes
# ---------------------------------------------------------------------------


_THREE_PB_HIGH_Q = (
    "three_pb does not converge above q ~ 0.5, and unlike the reason this "
    "xfail replaces, the system IS feasible: the model reproduces Bassett's "
    "K5 and K6 to three decimals, so the corrected target equals the model's "
    "own K_lat - K_str at the asked q and a root exists there by "
    "construction. This topology imposes every total pressure and leaves the "
    "flow level free, so the residual constrains only the DIFFERENCE and the "
    "solve has a one-parameter family to wander along. A solver and seeding "
    "problem, not a modelling one -- the opposite of what the xfail it "
    "replaces claimed. strict=True so it is noticed either way."
)


@pytest.mark.parametrize("q", [0.5, 0.6, 0.8])
@pytest.mark.xfail(strict=True, reason=_THREE_PB_HIGH_Q)
def test_three_pb_converges_at_high_lateral_fraction(model, q):
    assert _run(model, "three_pb", q).converged


def test_three_pb_low_lateral_fraction_lands_on_the_asked_point(model):
    """Was a strict xfail on the grounds that no root existed. It did exist --
    the target was built from two different operating points."""
    r = _run(model, "three_pb", 0.2)
    assert r.converged, r.message
    assert r.q_converged == pytest.approx(0.2, abs=0.005)


def test_mfb_two_pb_low_lateral_fraction_converges(model):
    """The other half of the same correction."""
    assert _run(model, "mfb_two_pb", 0.2).converged


def test_mfb_two_pb_mid_lateral_fraction_converges(model):
    """q=0.4 was an xfail target until the dividing-streamline recovery
    landed. It converges now, so it is pinned as passing."""
    assert _run(model, "mfb_two_pb", 0.4).converged
