"""The starting-point safety net and the caught-exception record.

The root finder's penalty already caught exceptions at rejected trial points;
the starting point was the one evaluation it did not cover, so a property
refusing a physical-but-extreme state there aborted ``solve()`` with its
retries (#552). Now that attempt ends with outcome ERROR, the retries run, and
every caught exception is returned in ``__caught_exceptions__``. Strict mode
(``raise_internal_errors`` or COMBAERO_SOLVER_STRICT, which this suite sets)
raises as before.
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    WallLayer,
)
from combaero.network import solver as solver_mod

Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _chain() -> FlowNetwork:
    g = FlowNetwork()
    for name, Pt in (("a", 1.1e5), ("b", 1.0e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 300.0, Y
        g.add_node(bnd)
    g.add_node(PlenumNode("p"))
    for eid, a, b in (("o1", "a", "p"), ("o2", "p", "b")):
        g.add_element(OrificeElement(eid, a, b, Cd=0.8, diameter=0.01, correlation="fixed"))
    return g


def _raise_on_calls(
    monkeypatch: pytest.MonkeyPatch, calls: set[int] | None, exc: type[Exception] = RuntimeError
) -> None:
    """Node 'p' raises ``exc`` on the given call numbers (None: always)."""
    original = PlenumNode.compute_derived_state
    n = [0]

    def patched(self: PlenumNode, upstream_states: list) -> tuple:
        n[0] += 1
        if self.id == "p" and (calls is None or n[0] in calls):
            raise exc("test: refused state")
        return original(self, upstream_states)

    monkeypatch.setattr(PlenumNode, "compute_derived_state", patched)


def _lenient(s: NetworkSolver) -> NetworkSolver:
    s.raise_internal_errors = False
    return s


def test_a_starting_point_that_raises_ends_the_attempt_with_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _raise_on_calls(monkeypatch, None)
    s = _lenient(NetworkSolver(_chain()))
    with pytest.warns(RuntimeWarning, match="starting point could not be evaluated at node 'p'"):
        r = s.solve(auto_retry=False)
    assert not r["__success__"] and not r["__converged__"]
    assert r["__outcome__"] == "error"
    assert "node 'p'" in r["__message__"]
    rec = r["__caught_exceptions__"][0]
    assert (rec["phase"], rec["where"], rec["type"]) == ("start", "node 'p'", "RuntimeError")
    assert rec["message"] == "test: refused state"
    assert rec["frames"][-1].startswith("test_solver_safety_net.py:")  # innermost: the raise
    assert r["__caught_exception_count__"] >= 1


def test_strict_mode_raises_from_the_attribute_and_the_environment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _raise_on_calls(monkeypatch, None)
    s = NetworkSolver(_chain())
    s.raise_internal_errors = True
    with pytest.raises(RuntimeError, match="refused state"):
        s.solve(auto_retry=False)
    # None follows the environment: this suite runs strict (conftest) ...
    monkeypatch.setenv(solver_mod.SOLVER_STRICT_ENV, "1")
    with pytest.raises(RuntimeError, match="refused state"):
        NetworkSolver(_chain()).solve(auto_retry=False)
    # ... and "0" is the default users get.
    monkeypatch.setenv(solver_mod.SOLVER_STRICT_ENV, "0")
    with pytest.warns(RuntimeWarning, match="starting point"):
        r = NetworkSolver(_chain()).solve(auto_retry=False)
    assert r["__outcome__"] == "error"


def test_the_retry_runs_after_a_cold_start_that_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """Only the first evaluation raises: the cold attempt ends ERROR and the
    compressible network's warm-start retry converges. Before the net, the
    exception left solve() and no retry ran."""
    _raise_on_calls(monkeypatch, {1})
    s = _lenient(NetworkSolver(_chain()))
    with pytest.warns(RuntimeWarning, match="starting point"):
        r = s.solve()
    assert r["__success__"]
    assert [c["phase"] for c in r["__caught_exceptions__"]] == ["start"]


@pytest.mark.parametrize("bug", [TypeError, AttributeError, KeyError, IndexError])
def test_a_code_bug_always_raises(monkeypatch: pytest.MonkeyPatch, bug: type[Exception]) -> None:
    _raise_on_calls(monkeypatch, None, bug)
    with pytest.raises(bug):
        _lenient(NetworkSolver(_chain())).solve(auto_retry=False)


def test_rejected_probes_are_recorded(monkeypatch: pytest.MonkeyPatch) -> None:
    """A trial point that raises was rejected before, silently; now it is
    in the record, with where and why, while the solve still converges."""
    _raise_on_calls(monkeypatch, {3})
    with pytest.warns(RuntimeWarning, match="unphysical intermediate state"):
        r = NetworkSolver(_chain()).solve(auto_retry=False)
    assert r["__success__"]
    recs = r["__caught_exceptions__"]
    assert len(recs) == 1 and recs[0]["phase"].startswith("probe")
    assert recs[0]["where"] == "node 'p'"


def test_a_clean_solve_records_nothing() -> None:
    r = NetworkSolver(_chain()).solve()
    assert r["__success__"]
    assert r["__caught_exceptions__"] == [] and r["__caught_exception_count__"] == 0


def test_an_unknown_wall_material_raises_at_construction() -> None:
    """It was swallowed at every evaluation, keeping the manual conductivity."""
    with pytest.raises(ValueError, match="unknown material 'inconel_718'"):
        WallLayer(thickness=1e-3, conductivity=20.0, material="inconel_718")
    assert (
        WallLayer(thickness=1e-3, conductivity=20.0, material="inconel718").material == "inconel718"
    )
    assert WallLayer(thickness=1e-3, conductivity=20.0, material="Custom").conductivity == 20.0
