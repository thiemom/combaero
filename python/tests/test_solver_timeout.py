import time

import pytest

from combaero.network import FlowNetwork, NetworkSolver, OrificeElement, PressureBoundary


def test_solver_no_timeout():
    """Verify solver still works normally when no timeout is hit."""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", P_total=200000, T_total=300)
    outlet = PressureBoundary("outlet", P_total=100000, T_total=300)
    orf = OrificeElement("orf", "inlet", "outlet", Cd=0.6, area=0.001)

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
    inlet = PressureBoundary("inlet", P_total=200000, T_total=300)
    outlet = PressureBoundary("outlet", P_total=100000, T_total=300)
    orf = OrificeElement("orf", "inlet", "outlet", Cd=0.6, area=0.001)

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(orf)

    solver = NetworkSolver(graph)

    # Inject an artificial delay in _residuals to force a timeout
    original_residuals = solver._residuals

    def slow_residuals(x):
        time.sleep(0.2)
        return original_residuals(x)

    solver._residuals = slow_residuals

    with pytest.warns(UserWarning, match="Solver timed out"):
        solution = solver.solve(timeout=0.1)  # 0.1s is shorter than one evaluation

    # It should return the initial guess (the 'best' seen so far)
    assert "orf.m_dot" in solution
    assert solution["orf.m_dot"] == 0.1  # Default initial guess for m_dot


def test_solver_maxfev_limit():
    """Verify solver respects maxfev for hybr method."""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", P_total=200000, T_total=300)
    outlet = PressureBoundary("outlet", P_total=100000, T_total=300)
    orf = OrificeElement("orf", "inlet", "outlet", Cd=0.6, area=0.001)

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(orf)

    solver = NetworkSolver(graph)

    # Force failure via low maxfev
    with pytest.warns(UserWarning, match="The number of calls to function has reached maxfev"):
        solution = solver.solve(method="hybr", options={"maxfev": 1})

    assert "orf.m_dot" in solution


def test_solver_unexpected_error_graceful_exit():
    """Verify solver returns best iterate on unexpected residuals error."""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", P_total=200000, T_total=300)
    outlet = PressureBoundary("outlet", P_total=100000, T_total=300)
    orf = OrificeElement("orf", "inlet", "outlet", Cd=0.6, area=0.001)

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(orf)

    solver = NetworkSolver(graph)
    original_residuals = solver._residuals

    call_count = 0

    def failing_residuals(x):
        nonlocal call_count
        call_count += 1
        if call_count > 1:
            raise RuntimeError("Unexpected failure")
        return original_residuals(x)

    solver._residuals = failing_residuals

    with pytest.warns(UserWarning, match="Unexpected error during residual evaluation"):
        solution = solver.solve()

    assert "orf.m_dot" in solution
