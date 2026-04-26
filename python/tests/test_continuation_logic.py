import numpy as np
import pytest

from combaero.network import (
    ChannelElement,
    FlowNetwork,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)


def test_continuation_logic():
    # 1. Setup a simple network: P_bound -> Channel -> Plenum -> Channel -> P_bound
    net = FlowNetwork()

    n1 = PressureBoundary("inlet", Pt=200000.0, Tt=300.0)
    n2 = PlenumNode("plenum")
    n3 = PressureBoundary("outlet", Pt=100000.0, Tt=300.0)

    e1 = ChannelElement("p1", "inlet", "plenum", length=1.0, diameter=0.1, roughness=1e-5)
    e2 = ChannelElement("p2", "plenum", "outlet", length=1.0, diameter=0.1, roughness=1e-5)

    net.add_node(n1)
    net.add_node(n2)
    net.add_node(n3)
    net.add_element(e1)
    net.add_element(e2)

    solver = NetworkSolver(net)

    # 3. Cold start solve
    res1 = solver.solve(init_strategy="default")
    assert res1["__success__"] is True
    x1 = res1["__x_solution__"]

    # 4. Continuation solve (no changes)
    # It should converge immediately
    res2 = solver.solve(init_strategy="continuation", x0=x1)
    assert res2["__success__"] is True
    x2 = res2["__x_solution__"]
    assert np.allclose(x1, x2, atol=1e-8)

    # 5. Parameter sweep (change inlet pressure)
    n1.Pt = 250000.0
    res3 = solver.solve(init_strategy="continuation", x0=x2)
    assert res3["__success__"] is True
    x3 = res3["__x_solution__"]
    # Check that it actually changed
    assert not np.allclose(x2, x3, atol=1e-3)

    # 6. Topology mismatch (length mismatch check)
    # Adding a node changes the number of unknowns
    # The solver check is usually a length check on x0
    with pytest.raises(ValueError, match="length.*does not match"):
        # Create a new solver with a different network
        net_new = FlowNetwork()
        net_new.add_node(PressureBoundary("inlet", Pt=200000.0))
        net_new.add_node(PlenumNode("p1"))
        net_new.add_node(PlenumNode("p2"))  # Extra node
        net_new.add_node(PressureBoundary("outlet", Pt=100000.0))
        net_new.add_element(ChannelElement("e1", "inlet", "p1", 1, 0.1, 1e-5))
        net_new.add_element(ChannelElement("e2", "p1", "p2", 1, 0.1, 1e-5))
        net_new.add_element(ChannelElement("e3", "p2", "outlet", 1, 0.1, 1e-5))

        solver_new = NetworkSolver(net_new)
        solver_new.solve(init_strategy="continuation", x0=x1)


def test_continuation_missing_x0():
    net = FlowNetwork()
    net.add_node(PressureBoundary("inlet", Pt=100000.0))
    net.add_node(PressureBoundary("outlet", Pt=90000.0))
    net.add_element(ChannelElement("p1", "inlet", "outlet", 1, 0.1, 1e-5))

    solver = NetworkSolver(net)
    with pytest.raises(ValueError, match="no prior converged solution"):
        solver.solve(init_strategy="continuation", x0=None)


if __name__ == "__main__":
    pytest.main([__file__])
