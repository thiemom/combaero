import time

import numpy as np
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

    # It should return the initial guess (the 'best' seen so far). Compare
    # against the solver's OWN initial guess rather than a constant: on a
    # network with no MassFlowBoundary the reference flow used to be a hard
    # 0.1 kg/s and is now estimated from the imposed pressure difference.
    assert "orf.m_dot" in solution
    seeded = dict(zip(solver.unknown_names, solver._build_x0(), strict=True))
    assert solution["orf.m_dot"] == pytest.approx(seeded["orf.m_dot"])


def test_the_warm_start_retry_stays_inside_the_timeout():
    """A network without a junction keeps its whole budget for the cold
    attempt, and the retry gets only what is left: none after a timeout. The
    proxy's 10 s floor once ran after a 0.1 s budget had expired (#585)."""
    graph = FlowNetwork()
    graph.add_node(PressureBoundary("inlet", Pt=200000, Tt=300))
    graph.add_node(PressureBoundary("outlet", Pt=100000, Tt=300))
    graph.add_element(
        OrificeElement(
            "orf",
            "inlet",
            "outlet",
            Cd=0.6,
            diameter=0.035682,
            correlation="fixed",
            regime="compressible",
        )
    )
    solver = NetworkSolver(graph)
    original_raj = solver._residuals_and_jacobian
    calls = [0]

    def slow_raj(x, **kwargs):
        calls[0] += 1
        time.sleep(0.2)
        return original_raj(x, **kwargs)

    solver._residuals_and_jacobian = slow_raj
    with pytest.warns(UserWarning, match="Solver timed out"):
        solution = solver.solve(timeout=0.1)
    assert not solution["__success__"]
    # The cold attempt's 3 evaluations; the proxy and the seeded retry after
    # it took 9 in all.
    assert calls[0] <= 4


def test_solver_maxfev_limit():
    """Verify solver respects maxfev for hybr method."""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", Pt=500000, Tt=300)
    outlet = PressureBoundary("outlet", Pt=100000, Tt=300)
    # Incompressible: its pressure-drop row is nonlinear in the flow. The
    # compressible row m - m_calc is linear in it, so one evaluation
    # converges and maxfev=1 cannot force a failure.
    orf = OrificeElement(
        "orf",
        "inlet",
        "outlet",
        Cd=0.6,
        diameter=0.112838,
        correlation="fixed",
        regime="incompressible",
    )

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(orf)

    solver = NetworkSolver(graph)

    # Force failure via low maxfev, from an x0 that is deliberately far from
    # the answer. Relying on a poor default seed to make this fail stopped
    # working when the reference flow started being estimated from the
    # imposed pressure difference rather than fixed at 0.1 kg/s: one
    # evaluation from a good seed can already be close enough.
    solver._build_x0()
    far = np.array(
        [
            1e-6 if name.endswith(".m_dot") else v
            for name, v in zip(solver.unknown_names, solver._build_x0(), strict=True)
        ]
    )
    with pytest.warns(UserWarning, match="NetworkSolver did not converge"):
        solution = solver.solve(method="hybr", x0=far, options={"maxfev": 1})

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
    from combaero.network.solver import _stall_detected

    trace = _trace([(0.5 * i, 1.0e4) for i in range(30)])
    assert _stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3)


def test_stall_detector_ignores_steady_improvement():
    """20% improvement per window is progress, not a stall."""
    from combaero.network.solver import _stall_detected

    trace = _trace([(0.5 * i, 1.0e4 * (0.8 ** (0.5 * i / 5.0))) for i in range(60)])
    assert not _stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3)


def test_stall_detector_ignores_slow_grind_near_tolerance():
    """Regression (live verification 2026-07-05): the outlet-ref seeded
    retry on the 5-tee at 0.18 kg/s crawls at |F| ~ 1e-2 (12x tol) and
    then converges. Cutting it over to LM turned a success into a
    failure, so plateaus below far_factor*tol must NOT count."""
    from combaero.network.solver import _stall_detected

    trace = _trace([(0.5 * i, 1.2e-2) for i in range(30)])
    assert not _stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3, far_factor=100.0)


def test_stall_detector_catches_a_dead_plateau_near_tolerance():
    """#626: a cold attempt's LM fallback parked at |F| = 0.0718, under the
    100x-tol gate, with not one better point in 175 s; a crawl is spared,
    a plateau that never improves over three windows is not."""
    from combaero.network.solver import _STALL_DEAD_WINDOWS, _stall_detected

    t_dead = _STALL_DEAD_WINDOWS * 5.0
    trace = _trace([(0.5 * i, 7.18e-2) for i in range(int(2 * t_dead / 0.5))])
    assert _stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3, far_factor=100.0)
    short = _trace([(0.5 * i, 7.18e-2) for i in range(int(t_dead / 0.5) - 1)])
    assert not _stall_detected(short, window=5.0, min_rel=0.10, tol=1e-3, far_factor=100.0)


def test_stall_detector_spares_a_crawl_near_tolerance():
    """Below the gate, any improvement at all is a crawl, not a stall."""
    from combaero.network.solver import _stall_detected

    trace = _trace([(0.5 * i, 7.18e-2 * (1.0 - 1e-4 * i)) for i in range(120)])
    assert not _stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3, far_factor=100.0)


def test_stall_detector_grace_period():
    """Histories shorter than the window never stall."""
    from combaero.network.solver import _stall_detected

    trace = _trace([(0.5 * i, 1.0e4) for i in range(8)])  # 3.5 s < 5 s
    assert not _stall_detected(trace, window=5.0, min_rel=0.10, tol=1e-3)
    assert not _stall_detected([], window=5.0)


def test_stall_error_is_timeout_subclass():
    """SolverStallError must route through every existing timeout handler
    (end the phase, keep the best iterate, reach the LM fallback)."""
    from combaero.network.solver import SolverStallError, SolverTimeoutError

    assert issubclass(SolverStallError, SolverTimeoutError)
