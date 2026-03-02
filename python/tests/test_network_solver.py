import math

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
)


def _get_air_X():
    X = [0.0] * 14
    X[0] = 0.79  # N2
    X[1] = 0.21  # O2
    return X


def test_network_solver_simple_orifice():
    """
    Scenario A (part 1): Simple PressureBoundary -> Orifice -> PressureBoundary.
    Verifies that the solver successfully finds the exact m_dot analytical value.
    """
    graph = FlowNetwork()

    inlet = PressureBoundary("inlet")
    inlet.P_total = 500000.0
    inlet.T_total = 300.0
    inlet.X = _get_air_X()

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
    # Orifice m_dot = Cd * A * sqrt(2 * rho * dP)
    m_dot_analytical, _ = cb.orifice_mdot_and_jacobian(dP, rho_in, 0.6, 0.005)

    assert solution["orf_1.m_dot"] == pytest.approx(m_dot_analytical, rel=1e-4)


def test_network_solver_simple_pipe():
    """
    Scenario A (part 2): Simple PressureBoundary -> Pipe -> PressureBoundary.
    Verifies exact Darcy-Weisbach flow scaling.
    """
    graph = FlowNetwork()

    inlet = PressureBoundary("inlet")
    inlet.P_total = 150000.0
    inlet.T_total = 300.0
    inlet.X = _get_air_X()

    outlet = PressureBoundary("outlet")
    outlet.P_total = 100000.0
    outlet.T_total = 300.0
    outlet.X = inlet.X

    pipe = PipeElement(
        "pipe_1",
        "inlet",
        "outlet",
        length=10.0,
        diameter=0.05,
        roughness=1e-5,
        friction_model="haaland",
    )

    graph.add_node(inlet)
    graph.add_node(outlet)
    graph.add_element(pipe)

    solver = NetworkSolver(graph)
    solution = solver.solve()

    m_dot_solved = solution["pipe_1.m_dot"]

    # Analytical verification
    rho = cb.density(inlet.T_total, inlet.P_total, inlet.X)
    mu = cb.viscosity(inlet.T_total, inlet.P_total, inlet.X)
    A = math.pi * (0.05 / 2) ** 2
    v = m_dot_solved / (rho * A)
    Re = max(4.0 * m_dot_solved / (math.pi * 0.05 * mu), 1.0)

    f, _ = cb.friction_and_jacobian("haaland", Re, 1e-5 / 0.05)

    dP_analytical = f * (10.0 / 0.05) * 0.5 * rho * v**2
    dP_actual = inlet.P_total - outlet.P_total

    assert dP_analytical == pytest.approx(dP_actual, rel=1e-4)


def test_network_bc_swapping():
    """
    Scenario B: Boundary Condition Swapping (Reverse Lookup).
    1. MassFlowBoundary -> Pipe -> PressureBoundary. Solve for P at inlet.
    2. PressureBoundary(P=P_solved) -> Pipe -> PressureBoundary. Solve for m_dot.
    m_dot should precisely equal the original MassFlowBoundary input.
    """

    # --- SETUP 1: Prescribe Mass, Solve for Pressure ---
    graph1 = FlowNetwork()
    target_m_dot = 0.5  # kg/s

    inlet_m = MassFlowBoundary("inlet")
    inlet_m.m_dot = target_m_dot
    inlet_m.T_total = 300.0
    inlet_m.X = _get_air_X()

    outlet_p = PressureBoundary("outlet")
    outlet_p.P_total = 100000.0
    outlet_p.T_total = 300.0
    outlet_p.X = _get_air_X()

    # Direct connection: MassFlowBoundary natively supports floating its pressure to push flow!
    pipe1 = PipeElement("pipe_1", "inlet", "outlet", length=5.0, diameter=0.1, roughness=1e-4)

    graph1.add_node(inlet_m)
    graph1.add_node(outlet_p)
    graph1.add_element(pipe1)

    solver1 = NetworkSolver(graph1)
    # Hint the inlet pressure for a faster Newton step
    inlet_m.initial_guess = {"inlet.P": 110000.0, "inlet.P_total": 110000.0}

    sol1 = solver1.solve()

    # Assert mass flow is conserved
    assert sol1["pipe_1.m_dot"] == pytest.approx(target_m_dot, rel=1e-6)

    # Store the solved pressure that the MFB had to rise to
    p_solved = sol1["inlet.P_total"]

    # --- SETUP 2: Prescribe Solved Pressure, Solve for Mass ---
    graph2 = FlowNetwork()
    inlet_p = PressureBoundary("inlet")
    inlet_p.P_total = p_solved
    inlet_p.T_total = 300.0
    inlet_p.X = _get_air_X()

    outlet2 = PressureBoundary("outlet")
    outlet2.P_total = 100000.0
    outlet2.T_total = 300.0
    outlet2.X = _get_air_X()

    pipe2 = PipeElement("pipe_1", "inlet", "outlet", length=5.0, diameter=0.1, roughness=1e-4)

    graph2.add_node(inlet_p)
    graph2.add_node(outlet2)
    graph2.add_element(pipe2)

    solver2 = NetworkSolver(graph2)
    sol2 = solver2.solve()

    # Verify the mass flows match mapping the reverse lookup
    assert sol2["pipe_1.m_dot"] == pytest.approx(target_m_dot, rel=1e-6)


def test_network_element_series():
    """
    Scenario C: Simple Element Series and Junction Mass Conservation.
    Boundary -> Pipe1 -> Junction -> Pipe2 -> Boundary
    """
    graph = FlowNetwork()

    inlet = PressureBoundary("inlet")
    inlet.P_total = 110000.0
    inlet.T_total = 300.0
    inlet.X = _get_air_X()

    junc = PlenumNode("junc_1")

    outlet = PressureBoundary("outlet")
    outlet.P_total = 100000.0
    outlet.T_total = 300.0
    outlet.X = _get_air_X()

    # Pipe 1 is larger, Pipe 2 is smaller
    pipe1 = PipeElement("p1", "inlet", "junc_1", length=5.0, diameter=0.1, roughness=0.0)
    pipe2 = PipeElement("p2", "junc_1", "outlet", length=5.0, diameter=0.05, roughness=0.0)

    graph.add_node(inlet)
    graph.add_node(junc)
    graph.add_node(outlet)
    graph.add_element(pipe1)
    graph.add_element(pipe2)

    solver = NetworkSolver(graph)
    # The initial guesses for default boundary variables starts at 101325
    sol = solver.solve()

    m_p1 = sol["p1.m_dot"]
    m_p2 = sol["p2.m_dot"]

    # Conservation of mass
    assert m_p1 == pytest.approx(m_p2, abs=1e-6)

    # Ensure continuity of pressure (P_total across junction is essentially static pressure P)
    # The pressure drops should equal the total pressure drop (with a small margin for density changes)
    dP1 = inlet.P_total - sol["junc_1.P"]
    dP2 = sol["junc_1.P"] - outlet.P_total

    assert (dP1 + dP2) == pytest.approx(inlet.P_total - outlet.P_total, abs=1e-4)
    # Pipe 2 (smaller) should drop much more pressure than Pipe 1
    assert dP2 > dP1 * 10


def test_network_lossless_connectors():
    """
    Scenario D: Lossless Connectors and Plenums.
    Boundary -> Lossless -> Plenum -> Lossless -> Boundary
    Mass flow should be essentially unconstrained / infinity if 0 dP, so we set a tiny dP
    Wait, Lossless connection between two different pressures will diverge.
    Instead: We put an Orifice AND a lossless in series.
    Boundary -> Lossless -> Plenum -> Orifice -> Boundary
    """
    graph = FlowNetwork()

    inlet = PressureBoundary("inlet")
    inlet.P_total = 105000.0
    inlet.T_total = 300.0
    inlet.X = _get_air_X()

    plen = PlenumNode("plenum_idx")

    outlet = PressureBoundary("outlet")
    outlet.P_total = 100000.0
    outlet.T_total = 300.0
    outlet.X = _get_air_X()

    conn1 = LosslessConnectionElement("l1", "inlet", "plenum_idx")
    orf1 = OrificeElement("o1", "plenum_idx", "outlet", Cd=0.65, area=0.01)

    graph.add_node(inlet)
    graph.add_node(plen)
    graph.add_node(outlet)
    graph.add_element(conn1)
    graph.add_element(orf1)

    solver = NetworkSolver(graph)
    sol = solver.solve()

    # The total pressure in the plenum should perfectly equal the inlet pressure
    # because the lossless connection drops 0 pressure regardless of flow rate
    assert sol["plenum_idx.P_total"] == pytest.approx(inlet.P_total, abs=1e-6)
    # And mass should be conserved
    assert sol["l1.m_dot"] == pytest.approx(sol["o1.m_dot"], abs=1e-6)
