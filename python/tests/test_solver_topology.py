import numpy as np
import pytest

import combaero as cb
from combaero.network.components import (
    MassFlowBoundary,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.graph import FlowNetwork
from combaero.network.solver import NetworkSolver


def test_all_nodes_only_P_Ptotal_unknowns():
    """
    Assert that the system size matches 2 * num_interior_nodes + num_variable_elements.
    In the new architecture, T and Y are derived, not unknowns.
    """
    net = FlowNetwork()

    # Boundary nodes (no unknowns)
    net.add_node(PressureBoundary("B1", P_total=101325.0, T_total=300.0))
    net.add_node(PressureBoundary("B2", P_total=100000.0, T_total=300.0))

    # Interior nodes (2 unknowns each: P, P_total)
    net.add_node(PlenumNode("N1"))
    net.add_node(PlenumNode("N2"))

    # Elements (1 unknown each: m_dot)
    net.add_element(OrificeElement("E1", "B1", "N1", Cd=0.6, area=0.001))
    net.add_element(OrificeElement("E2", "N1", "N2", Cd=0.6, area=0.001))
    net.add_element(OrificeElement("E3", "N2", "B2", Cd=0.6, area=0.001))

    solver = NetworkSolver(net)
    x0 = solver._build_x0()

    # 2 interior nodes * 2 + 3 elements * 1 = 7 unknowns
    expected_size = 2 * 2 + 3 * 1
    assert len(x0) == expected_size
    assert len(solver.unknown_names) == expected_size

    # Check naming pattern
    assert "N1.P" in solver.unknown_names
    assert "N1.P_total" in solver.unknown_names
    assert "E1.m_dot" in solver.unknown_names
    # T and Y should NOT be in unknown_names
    for name in solver.unknown_names:
        assert ".T" not in name
        assert ".Y[" not in name


def test_single_stream_passthrough():
    """
    Verify T and Y are preserved across a single-stream plenum.
    1 stream -> Plenum -> Orifice -> Boundary
    """
    net = FlowNetwork()
    # Inflow at 500K with specific composition
    # Using mole fractions: 70% N2, 30% O2
    X_in = np.zeros(14)
    X_in[0] = 0.7  # N2
    X_in[1] = 0.3  # O2
    Y_in = list(cb.mole_to_mass(X_in))

    net.add_node(PressureBoundary("IN", P_total=200000.0, T_total=500.0, Y=Y_in))
    net.add_node(PlenumNode("MID"))
    net.add_node(PressureBoundary("OUT", P_total=100000.0, T_total=300.0))

    net.add_element(OrificeElement("E1", "IN", "MID", Cd=0.6, area=0.01))
    net.add_element(OrificeElement("E2", "MID", "OUT", Cd=0.6, area=0.01))

    solver = NetworkSolver(net)
    sol = solver.solve()

    assert sol["__success__"]

    # T should be 500K (adiabatic passthrough)
    assert pytest.approx(sol["MID.T"], abs=1e-3) == 500.0
    # Y should match input
    for i in range(len(Y_in)):
        assert pytest.approx(sol[f"MID.Y[{i}]"], abs=1e-6) == Y_in[i]


def test_plenum_temperature_passthrough():
    """
    Verify that a plenum correctly mixes multiple streams.
    """
    net = FlowNetwork()

    # Stream 1: 1 kg/s at 300K
    net.add_node(MassFlowBoundary("MF1", m_dot=1.0, T_total=300.0))
    # Stream 2: 1 kg/s at 500K
    net.add_node(MassFlowBoundary("MF2", m_dot=1.0, T_total=500.0))

    net.add_node(PlenumNode("MIX"))
    net.add_node(PressureBoundary("OUT", P_total=100000.0, T_total=400.0))

    # Connect both to MIX
    net.add_element(OrificeElement("E1", "MF1", "MIX", Cd=0.8, area=0.1))
    net.add_element(OrificeElement("E2", "MF2", "MIX", Cd=0.8, area=0.1))
    net.add_element(OrificeElement("E3", "MIX", "OUT", Cd=0.8, area=0.1))

    solver = NetworkSolver(net)
    sol = solver.solve()

    assert sol["__success__"]

    # With equal mass flows and similar Cp (air), T should be ~400K
    assert pytest.approx(sol["MIX.T"], abs=2.0) == 400.0
    assert sol["MIX.T"] > 300.0
    assert sol["MIX.T"] < 500.0


if __name__ == "__main__":
    pytest.main([__file__])
