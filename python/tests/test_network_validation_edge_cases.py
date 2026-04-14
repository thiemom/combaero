import pytest

from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def test_validation_disjoint_boundaries():
    # Two separate subgraphs, each with its own PressureBoundary.
    # The current validation checks if every node is reachable from *some* PressureBoundary.
    # So this might pass validation, but the solver might fail or produce a singular matrix
    # if it doesn't solve them sequentially or if the solver expects a unified connected graph.
    # Actually, the requirement "No isolated subgraphs" implies the graph should be fully connected.
    # Let's see if this raises an error.
    graph = FlowNetwork()

    inlet_1 = PressureBoundary("inlet_1")
    outlet_1 = PressureBoundary("outlet_1")
    plenum_1 = PlenumNode("plenum_1")

    inlet_2 = PressureBoundary("inlet_2")
    outlet_2 = PressureBoundary("outlet_2")
    plenum_2 = PlenumNode("plenum_2")

    graph.add_node(inlet_1)
    graph.add_node(outlet_1)
    graph.add_node(plenum_1)

    graph.add_node(inlet_2)
    graph.add_node(outlet_2)
    graph.add_node(plenum_2)

    graph.add_element(OrificeElement("orf_1", "inlet_1", "plenum_1", 0.6, 0.05))
    graph.add_element(OrificeElement("orf_2", "plenum_1", "outlet_1", 0.6, 0.05))

    graph.add_element(OrificeElement("orf_3", "inlet_2", "plenum_2", 0.6, 0.05))
    graph.add_element(OrificeElement("orf_4", "plenum_2", "outlet_2", 0.6, 0.05))

    # We expect this to either be deemed valid (two separate systems) or invalid (not fully connected).
    # The current requirement says "No isolated subgraphs", but our implementation only checks
    # if all nodes are reachable from *some* BC. Let's see what happens.
    graph.validate()


def test_validation_multiple_edges_between_same_nodes():
    # Adding two elements between the exact same nodes.
    # E.g. parallel channels without intermediate plenums.
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet")
    outlet = PressureBoundary("outlet")

    graph.add_node(inlet)
    graph.add_node(outlet)

    # Element 1
    graph.add_element(OrificeElement("orf_1", "inlet", "outlet", 0.6, 0.05))
    # Element 2
    graph.add_element(OrificeElement("orf_2", "inlet", "outlet", 0.6, 0.02))

    # This should be structurally valid, but let's make sure it doesn't crash
    # the solver or validation.
    graph.validate()

    solver = NetworkSolver(graph)
    sol = solver.solve()
    assert "orf_1.m_dot" in sol
    assert "orf_2.m_dot" in sol


def test_validation_mass_flow_mismatch():
    # Only MassFlow boundaries pumping into a Plenum with no PressureBoundary sink.
    # The pressure should go to infinity.
    graph = FlowNetwork()

    in1 = MassFlowBoundary("in1")
    in2 = MassFlowBoundary("in2")
    plenum = PlenumNode("p1")

    graph.add_node(in1)
    graph.add_node(in2)
    graph.add_node(plenum)

    graph.add_element(OrificeElement("o1", "in1", "p1", 0.6, 0.05))
    graph.add_element(OrificeElement("o2", "in2", "p1", 0.6, 0.05))

    # Wait, our existing validator requires at least 1 pressure boundary, so this will fail early.
    with pytest.raises(ValueError, match="at least one PressureBoundary"):
        graph.validate()


def test_validation_mismatching_mass_flow_with_pressure_node():
    # Two mass flow boundaries forcing flow into a system that only has a pressure boundary.
    # This is solvable: the pressure boundary will act as a sink for both.
    graph = FlowNetwork()
    in1 = MassFlowBoundary("in1")
    in2 = MassFlowBoundary("in2")
    plenum = PlenumNode("p1")
    out = PressureBoundary("out")

    graph.add_node(in1)
    graph.add_node(in2)
    graph.add_node(plenum)
    graph.add_node(out)

    graph.add_element(OrificeElement("o1", "in1", "p1", 0.6, 0.05))
    graph.add_element(OrificeElement("o2", "in2", "p1", 0.6, 0.05))
    graph.add_element(OrificeElement("o3", "p1", "out", 0.6, 0.05))

    graph.validate()


def test_validation_closed_loop():
    # Validates that a closed loop with no boundaries fails,
    # OR a closed loop connected to a single pressure boundary (stagnant) is solvable.
    graph = FlowNetwork()
    p_ref = PressureBoundary("ref")
    p1 = PlenumNode("p1")
    p2 = PlenumNode("p2")
    p3 = PlenumNode("p3")

    graph.add_node(p_ref)
    graph.add_node(p1)
    graph.add_node(p2)
    graph.add_node(p3)

    graph.add_element(OrificeElement("o_ref", "ref", "p1", 0.6, 0.05))
    graph.add_element(OrificeElement("o_1", "p1", "p2", 0.6, 0.05))
    graph.add_element(OrificeElement("o_2", "p2", "p3", 0.6, 0.05))
    graph.add_element(OrificeElement("o_3", "p3", "p1", 0.6, 0.05))  # closed loop

    graph.validate()

    # Should solve to zero flow, since the only boundary is a single uniform pressure reservoir
    solver = NetworkSolver(graph)
    sol = solver.solve()
    assert abs(sol["o_1.m_dot"]) < 1e-5
    assert abs(sol["o_2.m_dot"]) < 1e-5
    assert abs(sol["o_3.m_dot"]) < 1e-5


def test_validation_elements_between_boundaries():
    # A single element connected between two boundaries.
    # While simple, it has no interior nodes. Let's make sure it doesn't fail
    # the generic interior node validation (which should be skipped for boundaries).
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet")
    outlet = PressureBoundary("outlet")

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(OrificeElement("orf_1", "inlet", "outlet", 0.6, 0.05))

    graph.validate()  # Should pass

    solver = NetworkSolver(graph)
    sol = solver.solve()
    assert "orf_1.m_dot" in sol
