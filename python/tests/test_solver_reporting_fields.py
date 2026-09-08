"""What a solve reports about itself.

`NetworkSolver.solve` used to answer only two questions: a boolean
`__success__` and a free-text `__message__`. That conflates things a caller
has to tell apart, and forces string matching for the rest:

- `__success__` is False both when Newton never got there AND when it found a
  root that a junction then rejected as unphysical. The second case reports a
  tiny `__final_norm__` alongside `success=False`, which reads as a
  contradiction until the two are reported separately.
- "why did it fail" was only available by matching substrings of
  `__message__` -- and part of that text comes from SciPy, which is free to
  reword it. `validation/junction/random_robustness.classify` did exactly
  that.
- "the elements were checked and passed" and "nothing checked anything" both
  came back as silence.

So: `__converged__` (the root finder's own verdict), `__consistent__`
(True/False/**None for not checked**), `__inconsistent_elements__`,
`__outcome__` (a `SolveOutcome`, always set) and `__worst_residuals__`.
`__success__` keeps its meaning exactly -- converged AND consistent -- because
callers depend on it.

Issue #272, follow-up to #299.
"""

from __future__ import annotations

import math
import random
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    NetworkSolver,
    PressureBoundary,
    SolveOutcome,
)
from validation.junction import random_robustness as rr


def _solve(net, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkSolver(net).solve(timeout=30.0, **kw)


@pytest.fixture(scope="module")
def outcomes():
    """A spread of real solves, so the fields are exercised on both verdicts."""
    rng = random.Random(20260906)
    out = []
    for _ in range(60):
        out.append(_solve(rr.build(rr.sample(rng))))
    assert any(s["__success__"] for s in out)
    assert any(not s["__success__"] for s in out)
    return out


# ---------------------------------------------------------------------------
# Always present, always set
# ---------------------------------------------------------------------------


def test_every_solve_reports_every_field(outcomes):
    for sol in outcomes:
        for key in (
            "__converged__",
            "__consistent__",
            "__inconsistent_elements__",
            "__outcome__",
            "__worst_residuals__",
        ):
            assert key in sol, f"{key} missing"
        assert isinstance(sol["__outcome__"], SolveOutcome)


def test_the_outcome_is_a_plain_string_too():
    """A StrEnum, so an existing `== "converged"` and any JSON encoder keep
    working. Pinned because switching to a bare Enum would break both
    silently."""
    assert SolveOutcome.CONVERGED == "converged"
    assert str(SolveOutcome.NO_PROGRESS) == "no_progress"
    assert f"{SolveOutcome.TIMEOUT}" == "timeout"


def test_a_network_with_no_unknowns_still_reports_an_outcome(monkeypatch):
    """The early return for a fully-constrained network used to hand back a
    three-key dict, so a caller reading `__outcome__` got a KeyError on the
    one case that never fails. Driven by emptying the unknown list rather
    than by an empty network, which the graph validator rejects outright."""
    Y = list(cb.mole_to_mass(cb.species.dry_air()))
    net = FlowNetwork()
    net.add_node(PressureBoundary("inlet", Pt=2.0e5, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("outlet", Pt=1.9e5, Tt=300.0, Y=Y))
    net.add_element(ChannelElement("ch", "inlet", "outlet", length=0.3, diameter=0.05))

    solver = NetworkSolver(net)
    real = solver._build_x0

    def no_unknowns():
        real()
        solver.unknown_names = []
        return np.array([])

    monkeypatch.setattr(solver, "_build_x0", no_unknowns)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = solver.solve(timeout=30.0)

    assert sol["__success__"] is True
    assert sol["__outcome__"] == SolveOutcome.NO_UNKNOWNS
    assert sol["__consistent__"] is None
    assert sol["__worst_residuals__"] == []


# ---------------------------------------------------------------------------
# success = converged AND consistent
# ---------------------------------------------------------------------------


def test_success_still_means_both(outcomes):
    for sol in outcomes:
        expected = bool(sol["__converged__"]) and sol["__consistent__"] is not False
        assert bool(sol["__success__"]) is expected


def test_a_converged_solve_reports_converged_and_the_converged_outcome(outcomes):
    for sol in outcomes:
        if not sol["__success__"]:
            continue
        assert sol["__converged__"] is True
        assert sol["__outcome__"] == SolveOutcome.CONVERGED
        assert sol["__inconsistent_elements__"] == []


def test_a_failed_solve_names_a_reason_other_than_converged(outcomes):
    for sol in outcomes:
        if sol["__success__"]:
            continue
        assert sol["__outcome__"] != SolveOutcome.CONVERGED


# ---------------------------------------------------------------------------
# "not checked" is not "fine"
# ---------------------------------------------------------------------------


def test_not_checked_reports_none_not_true():
    """A network with no element that owns a consistency check must report
    None. Reporting True would let a caller read the absence of any checking
    as a clean bill of health."""
    Y = list(cb.mole_to_mass(cb.species.dry_air()))
    net = FlowNetwork()
    net.add_node(PressureBoundary("inlet", Pt=2.0e5, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("outlet", Pt=1.9e5, Tt=300.0, Y=Y))
    net.add_element(ChannelElement("ch", "inlet", "outlet", length=0.3, diameter=0.05))

    sol = _solve(net)

    assert sol["__success__"] is True
    assert sol["__converged__"] is True
    assert sol["__consistent__"] is None, "nothing checked, so this cannot be True"


def test_a_non_converged_solve_does_not_claim_consistency(outcomes):
    """The checks only run on a converged solve, so an unconverged one has no
    verdict to report -- not a negative one."""
    for sol in outcomes:
        if sol["__converged__"]:
            continue
        assert sol["__consistent__"] is None


def test_a_checked_and_passing_solve_reports_true(outcomes):
    """The other side: where a junction DID verify the root, the field must be
    True rather than None, or it carries no information."""
    assert any(sol["__consistent__"] is True for sol in outcomes)


# ---------------------------------------------------------------------------
# The residual rows
# ---------------------------------------------------------------------------


def test_worst_residuals_are_named_sorted_and_real(outcomes):
    for sol in outcomes:
        worst = sol["__worst_residuals__"]
        if not worst:
            continue
        names = set(sol["__unknown_names__"])
        values = [row["residual"] for row in worst]
        assert all(row["name"] in names for row in worst)
        assert values == sorted(values, reverse=True)
        assert all(v >= 0.0 for v in values)


def test_the_worst_row_is_consistent_with_the_norm(outcomes):
    """The largest row cannot exceed the norm, nor be a vanishing fraction of
    it -- that would mean the two describe different states."""
    for sol in outcomes:
        worst = sol["__worst_residuals__"]
        norm = sol["__final_norm__"]
        if not worst or not math.isfinite(norm) or norm <= 0.0:
            continue
        assert worst[0]["residual"] <= norm * (1.0 + 1e-6)
        assert worst[0]["residual"] >= norm / math.sqrt(len(sol["__unknown_names__"])) * 0.99


# ---------------------------------------------------------------------------
# The consumer that used to string-match
# ---------------------------------------------------------------------------


def test_the_random_harness_classifies_from_the_outcome_not_the_message():
    """The proof that the string matching is gone: rewriting `__message__` to
    something unrecognisable must not change a single bucket."""
    rng = random.Random(4242)
    cases = [rr.sample(rng) for _ in range(25)]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        before = [rr.classify(rr.build(c)) for c in cases]

    original = NetworkSolver.solve

    def scrambled(self, *args, **kwargs):
        sol = original(self, *args, **kwargs)
        sol["__message__"] = "lorem ipsum dolor sit amet"
        return sol

    NetworkSolver.solve = scrambled
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            after = [rr.classify(rr.build(c)) for c in cases]
    finally:
        NetworkSolver.solve = original

    assert before == after, "the buckets still depend on the message text"
    assert set(before) != {"converged"}, "the sample must contain failures to be worth anything"


# ---------------------------------------------------------------------------
# The case the split exists for
# ---------------------------------------------------------------------------


def test_a_rejected_root_reports_converged_true_and_consistent_false():
    """The whole reason these are two fields. A root the junction rejects
    reports `success=False` with a SMALL residual norm, which is unreadable
    from `__success__` alone. Forced rather than sampled: the random harness
    produces no rejection in 80 draws, so this branch would otherwise be
    exercised by nothing.
    """
    rng = random.Random(20260906)
    for _ in range(40):
        net = rr.build(rr.sample(rng))
        solver = NetworkSolver(net)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            clean = solver.solve(timeout=30.0)
        if not clean["__success__"]:
            continue

        element = next(e for e in net.elements.values() if hasattr(e, "verify_solution_consistent"))
        element.verify_solution_consistent = lambda *a, **k: False
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = NetworkSolver(net).solve(timeout=30.0)

        assert sol["__converged__"] is True, "the root finder still got there"
        assert sol["__consistent__"] is False
        assert sol["__success__"] is False
        assert sol["__outcome__"] == SolveOutcome.INCONSISTENT
        assert element.id in sol["__inconsistent_elements__"]
        # The point: a tiny norm alongside success=False.
        assert sol["__final_norm__"] < 1e-3
        return
    pytest.skip("no converging draw available to reject")
