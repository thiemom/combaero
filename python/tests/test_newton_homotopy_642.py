"""A failed cold solve gets a last retry: a Newton homotopy from its cold guess (#642).

H(x, lam) = F(x) - (1 - lam) F(x0) is solved by the guess x0 at lam = 0 -- a
state with flow -- and lam walks to 1, each step warm-started from the last.
Opt-in (``NetworkSolver.homotopy_retry``; the GUI backend and the stress
harness turn it on): every solve that does not converge pays its budget.
Stress random draws 17, 25 and 33 failed every other start within 20 s and
now solve.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

import combaero as cb
from combaero.network import NetworkSolver, OrificeElement, PlenumNode, PressureBoundary
from combaero.network import solver as solver_mod
from combaero.network.graph import FlowNetwork

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402


def _solve(g: FlowNetwork, **kw) -> dict:
    s = NetworkSolver(g)
    s.homotopy_retry = True
    kw.setdefault("timeout", 60)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return s.solve(**kw)


def test_a_draw_every_other_start_failed_is_solved_by_the_homotopy(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The mechanism, not the machine: the 5 s budget is wall-clock, and a CI
    runner ran out of it before the rescue (on a desktop it takes ~1.6 s)."""
    monkeypatch.setattr(solver_mod, "NEWTON_HOMOTOPY_BUDGET", 120.0)
    g, _ = FAMILIES["random"](17)
    r = _solve(g, timeout=None)
    assert r["__success__"], r["__message__"]
    assert "Newton homotopy" in r["__message__"]
    # The physical checks ran on the final, offset-free solve.
    assert abs(r["__ledger__"]["mass_closure_rel"]) < 1e-6


def test_without_the_retry_that_draw_fails() -> None:
    """Off by default: every non-converging solve would pay its budget."""
    g, _ = FAMILIES["random"](17)
    s = NetworkSolver(g)
    assert s.homotopy_retry is False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=60)
    assert not r["__success__"]


def test_a_cold_solve_that_succeeds_never_enters_it() -> None:
    Y = list(cb.mole_to_mass(cb.species.dry_air()))
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2e5, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("p"))
    g.add_node(PressureBoundary("b", Pt=1e5, Tt=300.0, Y=Y))
    g.add_element(OrificeElement("o1", "a", "p", Cd=0.8, diameter=0.01, correlation="fixed"))
    g.add_element(OrificeElement("o2", "p", "b", Cd=0.8, diameter=0.01, correlation="fixed"))
    s = NetworkSolver(g)
    s.homotopy_retry = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"] and "homotopy" not in r["__message__"]
    assert getattr(s, "_homotopy_offset", None) is None


def test_its_budget_is_the_callers(monkeypatch: pytest.MonkeyPatch) -> None:
    """homotopy_budget caps it within the timeout; None (the GUI) gives it all
    the time the timeout leaves."""
    g, _ = FAMILIES["random"](17)
    s = NetworkSolver(g)
    s.homotopy_retry = True
    seen: list[float] = []
    real = s._solve_flow_retries

    def spy(**kw):
        if getattr(s, "_homotopy_offset", None) is not None:
            seen.append(kw["timeout"])
        return real(**kw)

    monkeypatch.setattr(s, "_solve_flow_retries", spy)
    for budget, ceiling in ((1.0, 1.0), (None, 40.0)):
        s.homotopy_budget = budget
        seen.clear()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s.solve(timeout=40.0)
        assert seen and max(seen) <= ceiling
    # None really is more than the default cap: the first step gets the rest
    # of the 40 s, not 5.
    assert max(seen) > solver_mod.NEWTON_HOMOTOPY_BUDGET
