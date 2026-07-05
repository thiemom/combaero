import time

import pytest

from combaero.network import FlowNetwork, NetworkSolver, OrificeElement, PressureBoundary


def test_solver_no_timeout():
    """Verify solver still works normally when no timeout is hit."""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", Pt=200000, Tt=300)
    outlet = PressureBoundary("outlet", Pt=100000, Tt=300)
    orf = OrificeElement("orf", "inlet", "outlet", Cd=0.6, diameter=0.035682, correlation="fixed")

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(orf)

    solver = NetworkSolver(graph)
    solution = solver.solve(timeout=1.0)  # 1 second is plenty

    assert "orf.m_dot" in solution
    assert solution["orf.m_dot"] > 0


def test_solver_timeout_trigger():
    """Verify solver triggers a timeout and returns the best iterate."""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", Pt=200000, Tt=300)
    outlet = PressureBoundary("outlet", Pt=100000, Tt=300)
    orf = OrificeElement("orf", "inlet", "outlet", Cd=0.6, diameter=0.035682, correlation="fixed")

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(orf)

    solver = NetworkSolver(graph)

    # Inject an artificial delay in _residuals_and_jacobian to force a timeout
    original_raj = solver._residuals_and_jacobian

    def slow_raj(x, **kwargs):
        time.sleep(0.2)
        return original_raj(x, **kwargs)

    solver._residuals_and_jacobian = slow_raj

    with pytest.warns(UserWarning, match="Solver timed out"):
        solution = solver.solve(timeout=0.1)  # 0.1s is shorter than one evaluation

    # It should return the initial guess (the 'best' seen so far)
    assert "orf.m_dot" in solution
    assert solution["orf.m_dot"] == 0.1  # Default initial guess for m_dot


def test_solver_maxfev_limit():
    """Verify solver respects maxfev for hybr method."""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", Pt=500000, Tt=300)
    outlet = PressureBoundary("outlet", Pt=100000, Tt=300)
    orf = OrificeElement("orf", "inlet", "outlet", Cd=0.6, diameter=0.112838, correlation="fixed")

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(orf)

    solver = NetworkSolver(graph)

    # Force failure via low maxfev
    with pytest.warns(UserWarning, match="NetworkSolver did not converge"):
        solution = solver.solve(method="hybr", options={"maxfev": 1})

    assert "orf.m_dot" in solution


def test_solver_unexpected_error_graceful_exit():
    """Verify solver returns best iterate on unexpected residuals error."""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", Pt=200000, Tt=300)
    outlet = PressureBoundary("outlet", Pt=100000, Tt=300)
    orf = OrificeElement("orf", "inlet", "outlet", Cd=0.6, diameter=0.035682, correlation="fixed")

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(orf)

    solver = NetworkSolver(graph)
    call_count = 0
    original_raj = solver._residuals_and_jacobian

    def failing_raj(x, **kwargs):
        nonlocal call_count
        call_count += 1
        if call_count > 1:
            raise RuntimeError("Unexpected failure")
        return original_raj(x, **kwargs)

    solver._residuals_and_jacobian = failing_raj

    with pytest.warns(UserWarning, match="Unexpected error during residual evaluation"):
        solution = solver.solve()

    assert "orf.m_dot" in solution


# ---------------------------------------------------------------------------
# Stall-triggered LM handover: plateau detector unit tests
# ---------------------------------------------------------------------------


def _trace(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Build a running-best (t, best|F|) trace from raw (t, |F|) points."""
    out: list[tuple[float, float]] = []
    best = float("inf")
    for t, norm in points:
        best = min(best, norm)
        out.append((t, best))
    return out


def test_stall_detector_fires_on_far_plateau():
    """A residual parked at |F| ~ 1e4 (far above tolerance) for longer
    than the window is a stall (the 16 kg/s 5-tee signature: hybr
    plateaus, LM-from-x0 converges in seconds)."""
    from combaero.network.solver import _hybr_stall_detected

    trace = _trace([(0.5 * i, 1.0e4) for i in range(30)])
    assert _hybr_stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3)


def test_stall_detector_ignores_steady_improvement():
    """20% improvement per window is progress, not a stall."""
    from combaero.network.solver import _hybr_stall_detected

    trace = _trace([(0.5 * i, 1.0e4 * (0.8 ** (0.5 * i / 5.0))) for i in range(60)])
    assert not _hybr_stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3)


def test_stall_detector_ignores_slow_grind_near_tolerance():
    """Regression (live verification 2026-07-05): the outlet-ref seeded
    retry on the 5-tee at 0.18 kg/s crawls at |F| ~ 1e-2 (12x tol) and
    then converges. Cutting it over to LM turned a success into a
    failure, so plateaus below far_factor*tol must NOT count."""
    from combaero.network.solver import _hybr_stall_detected

    trace = _trace([(0.5 * i, 1.2e-2) for i in range(30)])
    assert not _hybr_stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3, far_factor=100.0)


def test_stall_detector_grace_period():
    """Histories shorter than the window never stall."""
    from combaero.network.solver import _hybr_stall_detected

    trace = _trace([(0.5 * i, 1.0e4) for i in range(8)])  # 3.5 s < 5 s
    assert not _hybr_stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3)
    assert not _hybr_stall_detected([], window=5.0)


def test_stall_error_is_timeout_subclass():
    """SolverStallError must route through every existing timeout handler
    (end the phase, keep the best iterate, reach the LM fallback)."""
    from combaero.network.solver import SolverStallError, SolverTimeoutError

    assert issubclass(SolverStallError, SolverTimeoutError)
