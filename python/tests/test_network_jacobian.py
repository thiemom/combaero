import numpy as np
from scipy.optimize._numdiff import approx_derivative

from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
)


def _assert_jacobian_matches_fd(solver: NetworkSolver, max_rel_tol: float) -> None:
    x0 = solver._build_x0()

    _, jac_analytical = solver._residuals_and_jacobian(x0)
    jac_analytical = jac_analytical.toarray()

    rel_step = 1e-7
    abs_step = np.maximum(np.abs(x0) * rel_step, 1e-10)
    jac_numerical = approx_derivative(
        lambda x: solver._residuals(x), x0, method="2-point", abs_step=abs_step
    )

    diff = np.abs(jac_analytical - jac_numerical)
    mask = np.abs(jac_numerical) > 1e-5
    if np.any(mask):
        rel_diff = diff[mask] / np.abs(jac_numerical[mask])
        max_rel_diff = np.max(rel_diff)
        assert max_rel_diff < max_rel_tol, f"Jacobian mismatch! Max relative diff: {max_rel_diff}"
    else:
        max_diff = np.max(diff)
        assert max_diff < 1e-5


def test_jacobian_accuracy():
    """Global Jacobian FD check for a representative incompressible network."""
    net = FlowNetwork()

    # Nodes
    n_in = PressureBoundary("in", P_total=200000.0, T_total=300.0)
    n1 = PlenumNode("n1")
    n_out = PressureBoundary("out", P_total=101325.0, T_total=300.0)

    net.add_node(n_in)
    net.add_node(n1)
    net.add_node(n_out)

    # Elements
    orf = OrificeElement("orf", from_node="in", to_node="n1", Cd=0.6, diameter=0.112838)
    pipe = PipeElement(
        "pipe", from_node="n1", to_node="out", length=1.0, diameter=0.1, roughness=1e-5
    )

    net.add_element(orf)
    net.add_element(pipe)

    solver = NetworkSolver(net)
    _assert_jacobian_matches_fd(solver, max_rel_tol=0.01)


def test_jacobian_accuracy_compressible_network() -> None:
    """Global Jacobian FD check for a compressible pipe-orifice network."""
    net = FlowNetwork()

    n_in = PressureBoundary("in", P_total=4.8e5, T_total=540.0)
    n_mid = PlenumNode("mid")
    n_out = PressureBoundary("out", P_total=1.9e5, T_total=520.0)

    net.add_node(n_in)
    net.add_node(n_mid)
    net.add_node(n_out)

    net.add_element(
        PipeElement(
            "pipe",
            from_node="in",
            to_node="mid",
            length=3.2,
            diameter=0.022,
            roughness=8.0e-5,
            regime="compressible",
            friction_model="haaland",
        )
    )
    net.add_element(
        OrificeElement(
            "orf",
            from_node="mid",
            to_node="out",
            Cd=0.69,
            diameter=0.012361,
            regime="compressible",
        )
    )

    solver = NetworkSolver(net)
    _assert_jacobian_matches_fd(solver, max_rel_tol=0.05)


if __name__ == "__main__":
    test_jacobian_accuracy()
