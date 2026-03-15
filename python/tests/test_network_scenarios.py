import pytest

from combaero.network import (
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
)


def _get_air_Y():
    """Helper to create standard air composition."""
    Y = [0.0] * 14
    Y[0] = 0.79  # N2
    Y[1] = 0.21  # O2
    return Y


def build_pressure_bounds():
    """Create standard pressure boundaries."""
    inlet = PressureBoundary("inlet", P_total=150000.0, T_total=300.0, Y=_get_air_Y())
    outlet = PressureBoundary("outlet", P_total=100000.0, T_total=300.0, Y=_get_air_Y())
    return inlet, outlet


def test_step_1_simple_series():
    """Test: PB -> Pipe -> Plenum -> Pipe -> PB"""
    graph = FlowNetwork()
    inlet, outlet = build_pressure_bounds()

    n1 = PlenumNode("n1")
    p1 = PipeElement("p1", "inlet", "n1", length=2.0, diameter=0.1, roughness=1e-4)
    p2 = PipeElement("p2", "n1", "outlet", length=2.0, diameter=0.1, roughness=1e-4)

    for n in [inlet, outlet, n1]:
        graph.add_node(n)
    for e in [p1, p2]:
        graph.add_element(e)

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    # Verify solution exists and mass is conserved
    assert "p1.m_dot" in sol
    assert "p2.m_dot" in sol
    assert abs(sol["p1.m_dot"] - sol["p2.m_dot"]) < 1e-6
    assert sol["p1.m_dot"] > 0  # Positive flow from high to low pressure


def test_step_2_add_orifices():
    """Test: PB -> Pipe -> n1 -> Orifice -> mc1 -> Pipe -> PB"""
    graph = FlowNetwork()
    inlet, outlet = build_pressure_bounds()

    n1 = PlenumNode("n1")
    mc1 = MomentumChamberNode("mc1")

    p1 = PipeElement("p1", "inlet", "n1", length=2.0, diameter=0.1, roughness=1e-4)
    o1 = OrificeElement("o1", "n1", "mc1", Cd=0.6, area=0.005)
    p2 = PipeElement("p2", "mc1", "outlet", length=2.0, diameter=0.1, roughness=1e-4)

    for n in [inlet, outlet, n1, mc1]:
        graph.add_node(n)
    for e in [p1, o1, p2]:
        graph.add_element(e)

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    # Verify mass conservation across all elements
    assert abs(sol["p1.m_dot"] - sol["o1.m_dot"]) < 1e-6
    assert abs(sol["o1.m_dot"] - sol["p2.m_dot"]) < 1e-6
    assert sol["p1.m_dot"] > 0


def test_step_3_full_series_no_bypass():
    """Test: Full user series without bypass"""
    graph = FlowNetwork()
    inlet, outlet = build_pressure_bounds()

    n1 = PlenumNode("n1")
    mc1 = MomentumChamberNode("mc1")
    n2 = PlenumNode("n2")
    plen = PlenumNode("plenum")
    mc2 = MomentumChamberNode("mc2")

    p1 = PipeElement("p1", "inlet", "n1", length=2.0, diameter=0.1, roughness=1e-4)
    o1 = OrificeElement("o1", "n1", "mc1", Cd=0.6, area=0.008)
    p2 = PipeElement("p2", "mc1", "n2", length=2.0, diameter=0.1, roughness=1e-4)
    o2 = OrificeElement("o2", "n2", "plenum", Cd=0.6, area=0.008)
    p3 = PipeElement("p3", "plenum", "mc2", length=2.0, diameter=0.1, roughness=1e-4)
    p4 = PipeElement("p4", "mc2", "outlet", length=2.0, diameter=0.1, roughness=1e-4)

    for n in [inlet, outlet, n1, mc1, n2, plen, mc2]:
        graph.add_node(n)
    for e in [p1, o1, p2, o2, p3, p4]:
        graph.add_element(e)

    # Set initial guesses for better convergence
    n1.initial_guess = {"n1.P": 140000.0, "n1.P_total": 140000.0}
    mc1.initial_guess = {"mc1.P": 135000.0, "mc1.P_total": 135000.0}
    n2.initial_guess = {"n2.P": 130000.0, "n2.P_total": 130000.0}
    plen.initial_guess = {"plenum.P": 120000.0, "plenum.P_total": 120000.0}
    mc2.initial_guess = {"mc2.P": 110000.0, "mc2.P_total": 110000.0}

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    # Verify mass conservation across entire series
    mass_flows = [sol[f"{e}.m_dot"] for e in ["p1", "o1", "p2", "o2", "p3", "p4"]]
    for i in range(len(mass_flows) - 1):
        assert abs(mass_flows[i] - mass_flows[i + 1]) < 1e-6, (
            f"Mass not conserved between elements {i} and {i + 1}"
        )

    assert all(mf > 0 for mf in mass_flows), "All mass flows should be positive"


def test_step_4_adding_bypass():
    """Test: Adding bypass branch (mc1 -> orifice -> mc2)"""
    graph = FlowNetwork()
    inlet, outlet = build_pressure_bounds()

    n1 = PlenumNode("n1")
    mc1 = MomentumChamberNode("mc1")
    n2 = PlenumNode("n2")
    plen = PlenumNode("plenum")
    mc2 = MomentumChamberNode("mc2")

    # Main series elements
    p1 = PipeElement("p1", "inlet", "n1", length=2.0, diameter=0.1, roughness=1e-4)
    o1 = OrificeElement("o1", "n1", "mc1", Cd=0.6, area=0.008)
    p2 = PipeElement("p2", "mc1", "n2", length=2.0, diameter=0.1, roughness=1e-4)
    o2 = OrificeElement("o2", "n2", "plenum", Cd=0.6, area=0.008)
    p3 = PipeElement("p3", "plenum", "mc2", length=2.0, diameter=0.1, roughness=1e-4)
    p4 = PipeElement("p4", "mc2", "outlet", length=2.0, diameter=0.1, roughness=1e-4)

    # Bypass branch
    n_bypass = PlenumNode("n_bypass")
    p_bypass = PipeElement("p_bypass", "mc1", "n_bypass", length=5.0, diameter=0.08, roughness=1e-4)
    o_bypass = OrificeElement("o_bypass", "n_bypass", "mc2", Cd=0.6, area=0.005)

    # Add all nodes and elements
    for n in [inlet, outlet, n1, mc1, n2, plen, mc2, n_bypass]:
        graph.add_node(n)
    for e in [p1, o1, p2, o2, p3, p4, p_bypass, o_bypass]:
        graph.add_element(e)

    # Set initial guesses
    n1.initial_guess = {"n1.P": 145000.0, "n1.P_total": 145000.0}
    mc1.initial_guess = {"mc1.P": 140000.0, "mc1.P_total": 140000.0}
    n_bypass.initial_guess = {"n_bypass.P": 130000.0, "n_bypass.P_total": 130000.0}
    n2.initial_guess = {"n2.P": 135000.0, "n2.P_total": 135000.0}
    plen.initial_guess = {"plenum.P": 125000.0, "plenum.P_total": 125000.0}
    mc2.initial_guess = {"mc2.P": 115000.0, "mc2.P_total": 115000.0}

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    # Verify mass conservation at junctions
    # At mc1: p1.m_dot = p2.m_dot + p_bypass.m_dot
    assert abs(sol["p1.m_dot"] - (sol["p2.m_dot"] + sol["p_bypass.m_dot"])) < 1e-4

    # At mc2: p3.m_dot + bypass.m_dot = p4.m_dot (bypass joins here)
    assert abs((sol["p3.m_dot"] + sol["p_bypass.m_dot"]) - sol["p4.m_dot"]) < 1e-4

    # All flows should be positive
    assert all(sol[f"{e}.m_dot"] > 0 for e in ["p1", "p2", "p3", "p4", "p_bypass"])

    # Bypass should carry less flow than main branch (smaller diameter)
    assert sol["p_bypass.m_dot"] < sol["p2.m_dot"]


def test_network_validation_edge_cases():
    """Test additional edge cases for network validation"""

    # Test self-loop rejection
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", P_total=150000.0, T_total=300.0, Y=_get_air_Y())
    outlet = PressureBoundary("outlet", P_total=100000.0, T_total=300.0, Y=_get_air_Y())

    graph.add_node(inlet)
    graph.add_node(outlet)

    # Should raise ValueError for self-loop
    with pytest.raises(ValueError, match="Self-loops are not permitted"):
        graph.add_element(OrificeElement("self_loop", "inlet", "inlet", Cd=0.6, area=0.005))

    # Test isolated node detection
    graph2 = FlowNetwork()
    inlet2 = PressureBoundary("inlet2", P_total=150000.0, T_total=300.0, Y=_get_air_Y())
    outlet2 = PressureBoundary("outlet2", P_total=100000.0, T_total=300.0, Y=_get_air_Y())
    isolated = PlenumNode("isolated")

    graph2.add_node(inlet2)
    graph2.add_node(outlet2)
    graph2.add_node(isolated)
    graph2.add_element(OrificeElement("conn", "inlet2", "outlet2", Cd=0.6, area=0.005))

    # Should raise ValueError for isolated node
    with pytest.raises(ValueError, match="isolated node"):
        graph2.validate()


def test_parallel_elements():
    """Test parallel elements between same nodes"""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", P_total=150000.0, T_total=300.0, Y=_get_air_Y())
    outlet = PressureBoundary("outlet", P_total=100000.0, T_total=300.0, Y=_get_air_Y())

    graph.add_node(inlet)
    graph.add_node(outlet)

    # Two parallel orifices with different areas
    orf1 = OrificeElement("orf_1", "inlet", "outlet", Cd=0.6, area=0.005)
    orf2 = OrificeElement("orf_2", "inlet", "outlet", Cd=0.6, area=0.002)

    graph.add_element(orf1)
    graph.add_element(orf2)

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    # Both should have positive flow
    assert sol["orf_1.m_dot"] > 0
    assert sol["orf_2.m_dot"] > 0

    # Larger orifice should carry more flow
    assert sol["orf_1.m_dot"] > sol["orf_2.m_dot"]

    # Total flow should be reasonable
    total_flow = sol["orf_1.m_dot"] + sol["orf_2.m_dot"]
    assert total_flow > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
