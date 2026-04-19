import pytest

from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def test_validation_empty_network():
    graph = FlowNetwork()
    with pytest.raises(ValueError, match="contains no nodes"):
        graph.validate()

    inlet = PressureBoundary("inlet")
    inlet.P_total = 100000.0
    inlet.T_total = 300.0
    graph.add_node(inlet)

    with pytest.raises(ValueError, match="contains no elements"):
        graph.validate()


def test_validation_no_pressure_boundary():
    graph = FlowNetwork()
    node = PlenumNode("plenum1")
    graph.add_node(node)

    # Needs elements to pass the 'no elements' check first
    # But wait, Orifice needs to_node and from_node.
    node2 = PlenumNode("plenum2")
    graph.add_node(node2)
    elem = LosslessConnectionElement("conn", "plenum1", "plenum2")
    graph.add_element(elem)

    with pytest.raises(ValueError, match="must contain at least one PressureBoundary"):
        graph.validate()


def test_validation_self_loop():
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet")
    graph.add_node(inlet)

    elem = LosslessConnectionElement("conn", "inlet", "inlet")
    with pytest.raises(ValueError, match="Self-loops are not permitted"):
        graph.add_element(elem)


def test_validation_isolated_nodes():
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet")
    outlet = PressureBoundary("outlet")
    isolated = PlenumNode("isolated")

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_node(isolated)

    elem = OrificeElement("conn", "inlet", "outlet", Cd=0.6, diameter=1.128379, correlation="fixed")
    graph.add_element(elem)

    with pytest.raises(ValueError, match="Every node must be reachable from a PressureBoundary"):
        graph.validate()


def test_validation_unreachable_from_bc():
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet")
    graph.add_node(inlet)

    # Island subgraph with no BC
    p1 = PlenumNode("p1")
    p2 = PlenumNode("p2")
    conn = LosslessConnectionElement("conn", "p1", "p2")

    graph.add_node(p1)
    graph.add_node(p2)
    graph.add_element(conn)

    # Provide the main graph an element so it doesn't fail on 'no elements'
    # Actually wait, node 'inlet' has no connections. Let's fix that.
    outlet = PressureBoundary("outlet")
    graph.add_node(outlet)
    elem1 = OrificeElement(
        "conn1", "inlet", "outlet", Cd=0.6, diameter=1.128379, correlation="fixed"
    )
    graph.add_element(elem1)

    with pytest.raises(ValueError, match="Every node must be reachable from a PressureBoundary"):
        graph.validate()


def test_validation_interior_node_connections():
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet")
    outlet = PressureBoundary("outlet")
    p1 = PlenumNode("p1")

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_node(p1)

    # p1 is connected to inlet only. It is a dead end.
    conn = LosslessConnectionElement("conn", "inlet", "p1")
    graph.add_element(conn)

    # outlet is connected to inlet
    conn2 = OrificeElement(
        "conn2", "inlet", "outlet", Cd=0.6, diameter=1.128379, correlation="fixed"
    )
    graph.add_element(conn2)

    with pytest.raises(ValueError, match="Interior nodes must have at least 2 connections"):
        graph.validate()
