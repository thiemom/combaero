import time

import pytest

from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    NetworkNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)


class DivergentNode(NetworkNode):
    """A node that has an unsolvable residual equation: x^2 + 1 = 0"""

    def unknowns(self) -> list[str]:
        return [f"{self.id}.P", f"{self.id}.P_total"]

    def residuals(self, state) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Return a residual that can never be zero for any real P
        # We scale it so it doesn't immediately blow up but never converges
        res = [(state.P / 1e5) ** 2 + 1.0]
        # dR/dP = 2 * (P / 1e5) * (1/1e5)
        jac = {0: {f"{self.id}.P": 2 * state.P / 1e10}}
        return res, jac

    def resolve_topology(self, graph) -> None:
        pass


@pytest.mark.parametrize("method", NetworkSolver.SUPPORTED_METHODS)
def test_timeout_all_methods(method):
    """
    Verify that every supported method respects the wall-clock timeout.
    We use a sleep in the residual function to simulate a slow solver.
    """
    graph = FlowNetwork()
    # Inconsistent network: inlet=2bar, outlet=1bar, connected by lossless.
    # Plus a node that has an unsolvable residual: x^2 + 1 = 0
    inlet = PressureBoundary("inlet", P_total=200000.0, T_total=300.0)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0)
    node = DivergentNode("bad_node")

    conn1 = OrificeElement("orf1", "inlet", "bad_node", Cd=0.6, area=0.001)
    conn2 = LosslessConnectionElement("c2", "bad_node", "outlet")

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_node(node)
    graph.add_element(conn1)
    graph.add_element(conn2)

    solver = NetworkSolver(graph)

    original_residuals = solver._residuals

    def slow_residuals(x):
        time.sleep(0.2)  # 200ms per evaluation
        return original_residuals(x)

    solver._residuals = slow_residuals

    # Set a very short timeout.
    # Even if the solver only does 1-2 evaluations, it should hit 0.1s quickly.
    start = time.perf_counter()
    # Accept either "Solver timed out" or other non-convergence messages
    # as long as we don't hang indefinitely.
    with pytest.warns(UserWarning):
        solution = solver.solve(method=method, timeout=0.1)
    duration = time.perf_counter() - start

    # Ensure it didn't take much longer than the timeout + 1 iteration
    # 0.5s is a safe upper bound for a 0.1s timeout given 0.05s steps
    assert duration < 0.6, f"Method {method} took too long: {duration:.2f}s"
    assert "bad_node.P" in solution


def test_unsupported_method():
    """Verify that using an unsupported method raises a ValueError."""
    graph = FlowNetwork()
    solver = NetworkSolver(graph)
    with pytest.raises(ValueError, match="not supported"):
        solver.solve(method="invalid_method")


if __name__ == "__main__":
    # Allow running this script directly for quick verification
    for m in NetworkSolver.SUPPORTED_METHODS:
        print(f"Testing method: {m}...")
        try:
            test_timeout_all_methods(m)
            print(f"  {m}: PASSED")
        except Exception as e:
            print(f"  {m}: FAILED - {e}")
