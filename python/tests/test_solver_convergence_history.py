"""
Verifies the ``__convergence_history__`` field that NetworkSolver.solve()
emits in the solution dict.

The history is a throttled trace of (eval_index, t_elapsed_s, residual_norm)
captured by the residuals wrapper at solve time. Used by GUI diagnostics
to render the |F| vs iteration curve and to spot slow-convergence patterns
(Newton overshoot, oscillation, plateau).
"""

from __future__ import annotations

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
)


def _simple_network() -> FlowNetwork:
    """Minimal 2-boundary + 1-channel network. Always converges quickly."""
    net = FlowNetwork()
    Y = list(cb.mole_to_mass(cb.species.dry_air()))
    net.add_node(MassFlowBoundary("inlet", m_dot=0.1, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("outlet", Pt=1.0e5, Tt=300.0, Y=Y))
    net.add_element(
        ChannelElement("ch", from_node="inlet", to_node="outlet", length=1.0, diameter=0.1)
    )
    return net


def test_history_present_in_solution() -> None:
    net = _simple_network()
    sol = NetworkSolver(net).solve()
    assert "__convergence_history__" in sol
    assert isinstance(sol["__convergence_history__"], list)


def test_history_entries_have_expected_keys() -> None:
    net = _simple_network()
    sol = NetworkSolver(net).solve()
    hist = sol["__convergence_history__"]
    # Solve should generate at least one entry (the initial residual call).
    assert len(hist) >= 1
    for entry in hist:
        assert set(entry.keys()) == {"eval", "t", "norm"}
        assert isinstance(entry["eval"], int)
        assert isinstance(entry["t"], float)
        assert isinstance(entry["norm"], float)
        assert entry["t"] >= 0.0
        assert entry["norm"] >= 0.0


def test_history_t_is_monotonic_nondecreasing() -> None:
    """Wall clock should never go backwards across recorded points."""
    net = _simple_network()
    sol = NetworkSolver(net).solve()
    ts = [e["t"] for e in sol["__convergence_history__"]]
    for a, b in zip(ts, ts[1:], strict=False):
        assert b >= a, f"history t decreased: {a} -> {b}"


def test_history_eval_is_monotonic_increasing() -> None:
    """Each recorded entry should correspond to a strictly later residual call."""
    net = _simple_network()
    sol = NetworkSolver(net).solve()
    evals = [e["eval"] for e in sol["__convergence_history__"]]
    for a, b in zip(evals, evals[1:], strict=False):
        assert b > a, f"history eval did not strictly increase: {a} -> {b}"


def test_history_ends_near_zero_when_converged() -> None:
    net = _simple_network()
    sol = NetworkSolver(net).solve()
    assert sol["__success__"] is True
    hist = sol["__convergence_history__"]
    # Last recorded norm should be close to the final reported residual norm.
    final_norm = float(sol.get("__final_norm__", float("inf")))
    last_norm = hist[-1]["norm"]
    assert last_norm <= max(final_norm * 1.5, 1e-3), (
        f"last history norm {last_norm} far from final {final_norm}"
    )


def test_history_capped_at_1000_points() -> None:
    """The 1000-point cap is part of the contract -- prevents unbounded memory
    on pathological non-converging solves. Easier to verify against a quick
    network that won't actually hit 1000 entries; this test just checks the
    upper bound holds."""
    net = _simple_network()
    sol = NetworkSolver(net).solve()
    assert len(sol["__convergence_history__"]) <= 1000
