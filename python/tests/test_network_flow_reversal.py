from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)


def test_flow_reversal_simple_orifice():
    # 1 -> 2 where P2 > P1
    # Expected: m_dot should be negative, representing flow from 2 to 1
    graph = FlowNetwork()
    p1 = PressureBoundary("p1")
    p1.initial_guess = {
        "p1.P": 100000.0,
        "p1.P_total": 100000.0,
        "p1.T": 300.0,
        "p1.T_total": 300.0,
    }

    p2 = PressureBoundary("p2")
    p2.initial_guess = {
        "p2.P": 200000.0,
        "p2.P_total": 200000.0,
        "p2.T": 300.0,
        "p2.T_total": 300.0,
    }

    graph.add_node(p1)
    graph.add_node(p2)
    # The element is defined as flowing from p1 to p2
    graph.add_element(OrificeElement("orf", "p1", "p2", Cd=0.6, diameter=0.252313))

    graph.validate()
    solver = NetworkSolver(graph)
    sol = solver.solve()

    # Since P2 > P1, flow should be forced backwards from p2 to p1, hence negative mass flow
    assert "orf.m_dot" in sol
    assert sol["orf.m_dot"] < 0, f"Expected negative mass flow, got {sol['orf.m_dot']}"


def test_flow_reversal_momentum_chamber():
    # MomentumChamber: P_total > P_static
    # Let's test two cases:
    # 1. Flow entering chamber (should use static pressure of chamber to draw flow, then stagnate it)
    # 2. Flow leaving chamber (should use total pressure of chamber to drive flow)

    # For a momentum chamber where P_total > P_static:
    # Let's set P_total = 200kPa, P_static = 180kPa

    # CASE A: Flow leaving chamber to a boundary of 150kPa
    # Flow direction is chamber -> out
    graph_a = FlowNetwork()
    chamber_a = MomentumChamberNode("chamber")
    chamber_a.initial_guess = {"chamber.P": 180000.0, "chamber.P_total": 200000.0}  # solver guess
    # Flown is out of chamber
    inlet_a = PressureBoundary("inlet_a")
    inlet_a.initial_guess = {
        "inlet_a.P": 250000.0,
        "inlet_a.P_total": 250000.0,
        "inlet_a.T": 300.0,
        "inlet_a.T_total": 300.0,
    }

    outlet_a = PressureBoundary("outlet_a")
    outlet_a.initial_guess = {
        "outlet_a.P": 150000.0,
        "outlet_a.P_total": 150000.0,
        "outlet_a.T": 300.0,
        "outlet_a.T_total": 300.0,
    }

    graph_a.add_node(inlet_a)
    graph_a.add_node(chamber_a)
    graph_a.add_node(outlet_a)

    graph_a.add_element(LosslessConnectionElement("in_channel", "inlet_a", "chamber"))
    graph_a.add_element(
        OrificeElement("out_channel", "chamber", "outlet_a", Cd=0.6, diameter=0.112838)
    )

    graph_a.validate()
    sol_a = NetworkSolver(graph_a).solve()

    # Flow should be positive for both
    assert sol_a["in_channel.m_dot"] > 0
    assert sol_a["out_channel.m_dot"] > 0

    # CASE B: Reversal! The boundary is higher than the chamber.
    graph_b = FlowNetwork()
    chamber_b = MomentumChamberNode("chamber")

    # The outlet is actually at higher pressure than the inlet!
    inlet_b = PressureBoundary("inlet_b")
    inlet_b.initial_guess = {
        "inlet_b.P": 150000.0,
        "inlet_b.P_total": 150000.0,
        "inlet_b.T": 300.0,
        "inlet_b.T_total": 300.0,
    }

    outlet_b = PressureBoundary("outlet_b")
    outlet_b.initial_guess = {
        "outlet_b.P": 250000.0,
        "outlet_b.P_total": 250000.0,
        "outlet_b.T": 300.0,
        "outlet_b.T_total": 300.0,
    }

    graph_b.add_node(inlet_b)
    graph_b.add_node(chamber_b)
    graph_b.add_node(outlet_b)

    # Elements defined as inlet_b -> chamber -> outlet_b
    graph_b.add_element(LosslessConnectionElement("in_channel", "inlet_b", "chamber"))
    graph_b.add_element(
        OrificeElement("out_channel", "chamber", "outlet_b", Cd=0.6, diameter=0.112838)
    )

    graph_b.validate()
    sol_b = NetworkSolver(graph_b).solve()

    # Because outlet > inlet, flow is pushed backwards.
    # The elements are defined 1->2 but flow is 2->1.
    assert sol_b["in_channel.m_dot"] < 0
    assert sol_b["out_channel.m_dot"] < 0
