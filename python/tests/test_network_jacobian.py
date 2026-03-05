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


def test_jacobian_accuracy():
    """
    Verifies that the analytical sparse Jacobian assembled by NetworkSolver matches
    numerical finite differences within a tight tolerance.
    """
    # 1. Setup a representative network with various components
    net = FlowNetwork()

    # Nodes
    n_in = PressureBoundary("in", P_total=200000.0, T_total=300.0)
    n1 = PlenumNode("n1")
    n_out = PressureBoundary("out", P_total=101325.0, T_total=300.0)

    net.add_node(n_in)
    net.add_node(n1)
    net.add_node(n_out)

    # Elements
    orf = OrificeElement("orf", from_node="in", to_node="n1", Cd=0.6, area=0.01)
    pipe = PipeElement(
        "pipe", from_node="n1", to_node="out", length=1.0, diameter=0.1, roughness=1e-5
    )

    net.add_element(orf)
    net.add_element(pipe)

    solver = NetworkSolver(net)
    x0 = solver._build_x0()

    # 2. Get analytical Jacobian
    res0, jac_analytical = solver._residuals_and_jacobian(x0)
    jac_analytical = jac_analytical.toarray()

    # 3. Compute numerical Jacobian
    epsilon = 1e-7
    jac_numerical = approx_derivative(
        lambda x: solver._residuals(x), x0, method="2-point", abs_step=epsilon
    )

    # 4. Compare
    diff = np.abs(jac_analytical - jac_numerical)

    # Relative difference where Jacobian is significant
    mask = np.abs(jac_numerical) > 1e-5
    if np.any(mask):
        rel_diff = diff[mask] / np.abs(jac_numerical[mask])
        max_rel_diff = np.max(rel_diff)
        # We expect < 0.1% for most terms, allow slightly more for non-linear density effects
        assert max_rel_diff < 1e-3, f"Jacobian mismatch! Max relative diff: {max_rel_diff}"
    else:
        max_diff = np.max(diff)
        assert max_diff < 1e-5


if __name__ == "__main__":
    test_jacobian_accuracy()
