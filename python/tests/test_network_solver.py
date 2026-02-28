import combaero as cb
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
)


def test_network_solver_simple_orifice():
    """
    Test 1: Simple PressureBoundary -> Orifice -> PressureBoundary.
    Verifies that the solver successfully finds the exact m_dot analytical value.
    """
    graph = FlowNetwork()

    inlet = PressureBoundary("inlet")
    inlet.P_total = 500000.0
    inlet.T_total = 300.0
    inlet.X = [0.0] * 14
    inlet.X[0] = 0.79  # N2
    inlet.X[1] = 0.21  # O2

    outlet = PressureBoundary("outlet")
    outlet.P_total = 101325.0
    outlet.T_total = 300.0
    outlet.X = inlet.X

    orf = OrificeElement("orf_1", "inlet", "outlet", Cd=0.6, area=0.005)

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(orf)

    solver = NetworkSolver(graph)
    solution = solver.solve()

    assert "orf_1.m_dot" in solution

    # Analytical check
    rho_in = cb.density(inlet.T_total, inlet.P_total, inlet.X)
    dP = inlet.P_total - outlet.P_total
    expected_m_dot_incompressible = 0.6 * 0.005 * (2.0 * rho_in * dP) ** 0.5

    # Should be close
    assert abs(solution["orf_1.m_dot"] - expected_m_dot_incompressible) < 1e-4


def test_network_solver_pipe_plenum_orifice():
    """
    Test 2: Multi-node network
    Boundary -> Pipe -> Plenum -> Orifice -> Boundary
    Verifies mass conservation holds across junctions.
    """
    graph = FlowNetwork()

    inlet = PressureBoundary("inlet")
    inlet.P_total = 200000.0
    inlet.T_total = 300.0

    X_air = [0.0] * 14
    X_air[0] = 0.79
    X_air[1] = 0.21
    inlet.X = X_air

    plenum = PlenumNode("plenum")

    outlet = PressureBoundary("outlet")
    outlet.P_total = 101325.0
    outlet.T_total = 300.0
    outlet.X = X_air

    pipe = PipeElement("pipe_1", "inlet", "plenum", length=5.0, diameter=0.1, roughness=1e-5)
    orf = OrificeElement("orf_1", "plenum", "outlet", Cd=0.6, area=0.002)

    graph.add_node(inlet)
    graph.add_node(plenum)
    graph.add_node(outlet)

    graph.add_element(pipe)
    graph.add_element(orf)

    solver = NetworkSolver(graph)
    solution = solver.solve(method="hybr")

    m_dot_pipe = solution["pipe_1.m_dot"]
    m_dot_orf = solution["orf_1.m_dot"]

    assert abs(m_dot_pipe - m_dot_orf) < 1e-6, "Mass not conserved across plenum"
    assert "plenum.P" in solution
    assert "plenum.P_total" in solution


#
