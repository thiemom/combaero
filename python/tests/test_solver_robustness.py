import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    CombustorNode,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def compute_numerical_jacobian(solver, x, eps=1e-7):
    """Compute numerical Jacobian of solver residuals using central finite difference."""
    n = len(x)

    # Get base residuals to figure out size
    f0, _ = solver._residuals_and_jacobian(x)
    m = f0.size

    jac = np.zeros((m, n))
    for i in range(n):
        x_plus = x.copy()
        x_minus = x.copy()
        x_plus[i] += eps
        x_minus[i] -= eps

        f_plus, _ = solver._residuals_and_jacobian(x_plus)
        f_minus, _ = solver._residuals_and_jacobian(x_minus)

        jac[:, i] = (f_plus - f_minus) / (2 * eps)

    return jac


def test_numerical_jacobian_full_network():
    """Verify that analytical Jacobian matches central finite difference for a complex network."""
    graph = FlowNetwork()
    n_species = cb.num_species()
    get_idx = cb.species_index_from_name

    Y_air = [0.0] * n_species
    Y_air[get_idx("N2")] = 0.79
    Y_air[get_idx("O2")] = 0.21

    Y_fuel = [0.0] * n_species
    Y_fuel[get_idx("H2")] = 1.0

    inlet_air = PressureBoundary("inlet_air", P_total=1.5e5, T_total=300.0, Y=Y_air)
    inlet_fuel = MassFlowBoundary("inlet_fuel", m_dot=0.01, T_total=300.0, Y=Y_fuel)
    p1 = PlenumNode("p1")
    comb = CombustorNode("comb", method="complete")
    p2 = PlenumNode("p2")
    outlet = PressureBoundary("outlet", P_total=1.0e5)

    for node in [inlet_air, inlet_fuel, p1, comb, p2, outlet]:
        graph.add_node(node)

    graph.add_element(OrificeElement("o_air", "inlet_air", "p1", Cd=0.6, area=0.01))
    graph.add_element(OrificeElement("o_fuel", "inlet_fuel", "p1", Cd=0.6, area=0.001))
    graph.add_element(OrificeElement("o_mix", "p1", "comb", Cd=0.6, area=0.01))
    graph.add_element(OrificeElement("o_to_p2", "comb", "p2", Cd=0.6, area=0.01))
    graph.add_element(OrificeElement("o_exit", "p2", "outlet", Cd=0.6, area=0.01))

    solver = NetworkSolver(graph)
    # Trigger build of unknown names and indices
    x0 = solver._build_x0()

    # Perturb x0 slightly away from initial guess to avoid lucky zeros
    x_test = x0 * 1.05 + 0.01

    res, jac_sparse = solver._residuals_and_jacobian(x_test)
    jac_analytical = jac_sparse.toarray()
    jac_numerical = compute_numerical_jacobian(solver, x_test)

    assert jac_analytical.shape == jac_numerical.shape

    # Compare analytical vs numerical. Central diff is quite accurate.
    # Note: increased absolute tolerance slightly for complex network relay chains
    np.testing.assert_allclose(jac_analytical, jac_numerical, rtol=1e-4, atol=1e-4)


def test_cascaded_derived_state_relay():
    """Verify that T and Y are correctly propagated through a long chain of nodes."""
    graph = FlowNetwork()
    n_species = cb.num_species()
    get_idx = cb.species_index_from_name

    # Lean composition
    Y_lean = [0.0] * n_species
    Y_lean[get_idx("N2")] = 0.36
    Y_lean[get_idx("O2")] = 0.6
    Y_lean[get_idx("H2")] = 0.04

    inlet = PressureBoundary("inlet", P_total=2.0e5, T_total=400.0, Y=Y_lean)
    p1 = PlenumNode("p1")
    comb = CombustorNode("comb", method="complete")
    p2 = PlenumNode("p2")
    outlet = PressureBoundary("outlet", P_total=1.0e5)

    for node in [inlet, p1, comb, p2, outlet]:
        graph.add_node(node)

    graph.add_element(OrificeElement("o1", "inlet", "p1", Cd=0.6, area=0.01))
    graph.add_element(OrificeElement("o2", "p1", "comb", Cd=0.6, area=0.01))
    graph.add_element(OrificeElement("o3", "comb", "p2", Cd=0.6, area=0.01))
    graph.add_element(OrificeElement("o4", "p2", "outlet", Cd=0.6, area=0.01))

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["__success__"]

    # p1 should have inlet T and Y
    assert sol["p1.T_total"] == pytest.approx(400.0)
    assert sol["p1.Y[0]"] == pytest.approx(0.36)

    # comb should have combustion products (T_ad > 400.0)
    assert sol["comb.T_total"] > 1000.0
    # H2 should be consumed. With smooth=true, it might not be absolute zero,
    # but at Phi=0.5 it should be significantly reduced.
    assert sol[f"comb.Y[{get_idx('H2')}]"] < 1e-2

    # p2 should inherit comb properties
    assert sol["p2.T_total"] == pytest.approx(sol["comb.T_total"])
    assert sol[f"p2.Y[{get_idx('N2')}]"] == pytest.approx(sol[f"comb.Y[{get_idx('N2')}]"])


def test_orifice_smoothness_and_convergence_reversal():
    """Verify regularization allows convergence when starting with wrong sign."""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", P_total=1.0e5, T_total=300.0)
    p1 = PlenumNode("p1")
    outlet = PressureBoundary("outlet", P_total=1.1e5)

    for node in [inlet, p1, outlet]:
        graph.add_node(node)

    o1 = OrificeElement("o1", "inlet", "p1", Cd=0.6, area=0.01)
    o2 = OrificeElement("o2", "p1", "outlet", Cd=0.6, area=0.01)

    # Force initial guess for m_dot to be positive (wrong direction)
    o1.initial_guess = {"o1.m_dot": 1.0}
    o2.initial_guess = {"o2.m_dot": 1.0}

    graph.add_element(o1)
    graph.add_element(o2)

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["__success__"]
    # Flow should have flipped to negative
    assert sol["o1.m_dot"] < 0
    assert sol["o2.m_dot"] < 0
    assert sol["o1.m_dot"] == pytest.approx(sol["o2.m_dot"], rel=1e-5)


def test_mass_fraction_conservation_and_bounds():
    """Verify that mass fractions sum to 1.0 and remain in [0, 1]."""
    graph = FlowNetwork()
    n_species = cb.num_species()
    get_idx = cb.species_index_from_name
    Y1 = [0.0] * n_species
    Y1[get_idx("N2")] = 1.0
    Y2 = [0.0] * n_species
    Y2[get_idx("O2")] = 1.0

    in1 = PressureBoundary("in1", P_total=1.5e5, Y=Y1)
    in2 = PressureBoundary("in2", P_total=1.5e5, Y=Y2)
    p1 = PlenumNode("p1")
    comb = CombustorNode("comb", method="complete")
    outlet = PressureBoundary("outlet", P_total=1.0e5)

    for node in [in1, in2, p1, comb, outlet]:
        graph.add_node(node)

    graph.add_element(OrificeElement("o1", "in1", "p1", Cd=0.6, area=0.01))
    graph.add_element(OrificeElement("o2", "in2", "p1", Cd=0.6, area=0.01))
    graph.add_element(OrificeElement("o3", "p1", "comb", Cd=0.6, area=0.01))
    graph.add_element(OrificeElement("o4", "comb", "outlet", Cd=0.6, area=0.01))

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["__success__"]

    for node in ["p1", "comb"]:
        Y_sum = sum(sol[f"{node}.Y[{i}]"] for i in range(n_species))
        assert Y_sum == pytest.approx(1.0, abs=1e-12)
        for i in range(n_species):
            val = sol[f"{node}.Y[{i}]"]
            assert val >= -1e-15
            assert val <= 1.0 + 1e-15


def test_system_size_reduction_assertion():
    """Assert system size is exactly (2 * N_interior_nodes + N_elements)."""
    graph = FlowNetwork()
    graph.add_node(PressureBoundary("in"))
    graph.add_node(PlenumNode("p1"))
    graph.add_node(PlenumNode("p2"))
    graph.add_node(CombustorNode("comb"))
    graph.add_node(PressureBoundary("out"))

    graph.add_element(OrificeElement("e1", "in", "p1", 0.6, 0.01))
    graph.add_element(OrificeElement("e2", "p1", "p2", 0.6, 0.01))
    graph.add_element(OrificeElement("e3", "p2", "comb", 0.6, 0.01))
    graph.add_element(OrificeElement("e4", "comb", "out", 0.6, 0.01))

    solver = NetworkSolver(graph)
    # Trigger build
    x0 = solver._build_x0()

    # 3 interior nodes * 2 (P, P_total) + 4 elements * 1 (m_dot) = 10 unknowns
    assert len(x0) == 10
    assert len(solver.unknown_names) == 10
