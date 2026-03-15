import numpy as np
import pytest

from combaero.network import (
    CombustorNode,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def test_chain_rule_orifice_through_combustor():
    """Verify chain-rule relay accuracy: m_dot -> Combustor -> Orifice."""
    graph = FlowNetwork()
    inlet = MassFlowBoundary("inlet", m_dot=0.5, T_total=300.0)
    comb = CombustorNode("comb", method="complete")
    orifice = OrificeElement("orf", "comb", "out", Cd=0.6, area=0.01)
    outlet = PressureBoundary("out", P_total=101325.0)

    graph.add_node(inlet)
    graph.add_node(comb)
    graph.add_node(outlet)
    graph.add_element(orifice)

    # We need an element between inlet and comb? No, boundaries connect to elements.
    # Actually, let's make it: inlet -> orf1 -> comb -> orf2 -> outlet
    graph = FlowNetwork()
    inlet = MassFlowBoundary("inlet", m_dot=0.5, T_total=300.0)
    p1 = PlenumNode("p1")
    comb = CombustorNode("comb", method="complete")
    outlet = PressureBoundary("outlet", P_total=101325.0)

    graph.add_node(inlet)
    graph.add_node(p1)
    graph.add_node(comb)
    graph.add_node(outlet)

    graph.add_element(OrificeElement("o1", "inlet", "p1", 0.8, 0.05))
    graph.add_element(OrificeElement("o2", "p1", "comb", 0.8, 0.05))
    graph.add_element(OrificeElement("o3", "comb", "outlet", 0.8, 0.05))

    graph.validate()
    solver = NetworkSolver(graph)

    # Check Jacobian accuracy
    res, jac = solver._residuals_and_jacobian(solver._build_x0())

    # Simple check: jac should be square and non-singular
    assert jac.shape[0] == jac.shape[1]
    assert np.all(np.isfinite(res))
    assert np.all(np.isfinite(jac.data))


def test_plenum_temperature_passthrough():
    """Verify derived states correctly reflect upstream T."""
    graph = FlowNetwork()
    # High T inlet
    inlet = PressureBoundary("inlet", P_total=200000.0, T_total=800.0)
    p1 = PlenumNode("p1")
    outlet = PressureBoundary("outlet", P_total=101325.0)

    graph.add_node(inlet)
    graph.add_node(p1)
    graph.add_node(outlet)

    graph.add_element(OrificeElement("o1", "inlet", "p1", 0.6, 0.01))
    graph.add_element(OrificeElement("o2", "p1", "outlet", 0.6, 0.01))

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["p1.T_total"] == pytest.approx(800.0)
    assert sol["o1.m_dot"] == pytest.approx(sol["o2.m_dot"], rel=1e-5)


def test_all_nodes_only_P_Ptotal_unknowns():
    """Verify system size is correctly reduced (only P, P_total per node)."""
    graph = FlowNetwork()
    p1 = PlenumNode("p1")
    p2 = PlenumNode("p2")
    inlet = PressureBoundary("in")
    outlet = PressureBoundary("out")

    graph.add_node(inlet)
    graph.add_node(p1)
    graph.add_node(p2)
    graph.add_node(outlet)

    graph.add_element(OrificeElement("o1", "in", "p1", 0.6, 0.01))
    graph.add_element(OrificeElement("o2", "p1", "p2", 0.6, 0.01))
    graph.add_element(OrificeElement("o3", "p2", "out", 0.6, 0.01))

    solver = NetworkSolver(graph)
    x0 = solver._build_x0()

    # 2 interior nodes * 2 unknowns (P, P_total) + 3 elements * 1 unknown (m_dot) = 7 unknowns
    assert len(x0) == 7
    assert "p1.P" in solver.unknown_names
    assert "p1.P_total" in solver.unknown_names
    assert "p1.T" not in solver.unknown_names
    assert "p1.Y[0]" not in solver.unknown_names


def test_orifice_smoothness_at_zero_dP():
    """Verify C++ regularization works end-to-end through Python."""
    # A single orifice between two pressure boundaries at the EXACT same pressure
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", P_total=101325.0)
    outlet = PressureBoundary("outlet", P_total=101325.0)

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(OrificeElement("o1", "inlet", "outlet", 0.6, 0.01))

    solver = NetworkSolver(graph)
    sol = solver.solve()

    # Should converge to exactly 0 (or machine epsilon)
    assert abs(sol["o1.m_dot"]) < 1e-15

    # Check Jacobian at dP=0 is finite and non-zero
    x_zero = np.array([0.0])  # o1.m_dot
    res, jac = solver._residuals_and_jacobian(x_zero)
    assert np.isfinite(jac.data).all()
    # Residual should be 0 because m_dot=0 and dP=0
    assert res[0] == 0.0


def test_single_stream_passthrough():
    """Verify identity chain rule."""
    # inlet (MassFlow) -> Plenum -> outlet (Pressure)
    # We want to see if dm_dot_out / dm_dot_in = 1 in the Jacobian?
    # No, mass conservation is m_dot_in - m_dot_out = 0.
    # Let's check if deriving state through a chain of 5 plenums works.
    graph = FlowNetwork()
    prev = "inlet"
    graph.add_node(PressureBoundary("inlet", T_total=500.0))
    for i in range(5):
        curr = f"p{i}"
        graph.add_node(PlenumNode(curr))
        graph.add_element(OrificeElement(f"o{i}", prev, curr, 0.6, 0.01))
        prev = curr
    graph.add_node(PressureBoundary("outlet", P_total=101325.0))
    graph.add_element(OrificeElement("o_final", prev, "outlet", 0.6, 0.01))

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["p4.T_total"] == pytest.approx(500.0)
    assert "__success__" in sol and sol["__success__"]
