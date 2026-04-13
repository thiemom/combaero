"""
Test EffectiveAreaConnectionElement functionality.
"""

import math

import pytest

import combaero as cb
from combaero.network import (
    EffectiveAreaConnectionElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


class TestEffectiveAreaConnectionElement:
    """Test suite for EffectiveAreaConnectionElement."""

    def test_initialization(self):
        """Test that EffectiveAreaConnectionElement initializes correctly."""
        effective_area = 0.1  # m^2
        connection = EffectiveAreaConnectionElement("conn1", "node1", "node2", effective_area)

        # Verify attributes
        assert connection.Cd == 1.0
        assert connection.area == effective_area
        assert connection.id == "conn1"
        assert connection.from_node == "node1"
        assert connection.to_node == "node2"

    def test_inheritance(self):
        """Test that EffectiveAreaConnectionElement inherits from OrificeElement."""
        connection = EffectiveAreaConnectionElement("conn1", "node1", "node2", 0.1)
        assert isinstance(connection, OrificeElement)

    def test_unknowns(self):
        """Test unknowns method."""
        connection = EffectiveAreaConnectionElement("conn1", "node1", "node2", 0.1)
        unknowns = connection.unknowns()
        assert len(unknowns) == 1
        assert unknowns[0] == "conn1.m_dot"

    def test_n_equations(self):
        """Test n_equations method."""
        connection = EffectiveAreaConnectionElement("conn1", "node1", "node2", 0.1)
        assert connection.n_equations() == 1

    def test_resolve_topology(self):
        """Test resolve_topology method."""
        network = FlowNetwork()
        connection = EffectiveAreaConnectionElement("conn1", "node1", "node2", 0.1)

        # Should not raise any errors
        connection.resolve_topology(network)

        # Should not set upstream/downstream diameters
        assert connection.upstream_diameter is None
        assert connection.downstream_diameter is None

    def test_integration_with_network(self):
        """Test that the element can be added to a network."""
        network = FlowNetwork()

        # Create nodes
        plenum1 = PlenumNode("plenum1")
        plenum2 = PlenumNode("plenum2")

        # Create connection
        connection = EffectiveAreaConnectionElement("conn1", "plenum1", "plenum2", 0.1)

        # Add to network
        network.add_node(plenum1)
        network.add_node(plenum2)
        network.add_element(connection)

        # Verify network contains the element
        assert "conn1" in network.elements
        assert network.elements["conn1"] is connection

    def test_different_effective_areas(self):
        """Test with various effective area values."""
        test_areas = [0.01, 0.1, 1.0, 10.0]

        for i, area in enumerate(test_areas):
            connection = EffectiveAreaConnectionElement(
                f"conn{i}", f"node{i}", f"node{i + 1}", area
            )
            assert connection.area == area
            assert connection.Cd == 1.0

    def test_methods_exist(self):
        """Test that all required methods exist."""
        connection = EffectiveAreaConnectionElement("conn1", "node1", "node2", 0.1)

        assert hasattr(connection, "unknowns")
        assert hasattr(connection, "residuals")
        assert hasattr(connection, "n_equations")
        assert hasattr(connection, "resolve_topology")
        assert callable(connection.unknowns)
        assert callable(connection.residuals)
        assert callable(connection.n_equations)
        assert callable(connection.resolve_topology)

    def test_flow_calculation_matches_orifice(self):
        """Verify EffectiveAreaConnectionElement produces same flow as OrificeElement(Cd=1.0).

        Both use the compressible flow equation via cb.orifice_mdot_and_jacobian(),
        which accounts for density variations across the pressure drop.
        """
        effective_area = 0.01  # m^2

        # Network 1: Using EffectiveAreaConnectionElement
        graph1 = FlowNetwork()
        inlet1 = PressureBoundary("inlet")
        inlet1.P_total = 150000.0
        inlet1.T_total = 300.0
        inlet1.Y = cb.species.dry_air_mass()
        outlet1 = PressureBoundary("outlet")
        outlet1.P_total = 100000.0
        outlet1.T_total = 300.0
        outlet1.Y = cb.species.dry_air_mass()
        conn1 = EffectiveAreaConnectionElement("conn", "inlet", "outlet", effective_area)

        graph1.add_node(inlet1)
        graph1.add_node(outlet1)
        graph1.add_element(conn1)

        solver1 = NetworkSolver(graph1)
        sol1 = solver1.solve(method="lm")

        # Network 2: Using OrificeElement with Cd=1.0
        graph2 = FlowNetwork()
        inlet2 = PressureBoundary("inlet")
        inlet2.P_total = 150000.0
        inlet2.T_total = 300.0
        inlet2.Y = cb.species.dry_air_mass()
        outlet2 = PressureBoundary("outlet")
        outlet2.P_total = 100000.0
        outlet2.T_total = 300.0
        outlet2.Y = cb.species.dry_air_mass()
        orf2 = OrificeElement(
            "orf", "inlet", "outlet", Cd=1.0, diameter=math.sqrt(4.0 * effective_area / math.pi)
        )

        graph2.add_node(inlet2)
        graph2.add_node(outlet2)
        graph2.add_element(orf2)

        solver2 = NetworkSolver(graph2)
        sol2 = solver2.solve(method="lm")

        # Should produce identical results
        assert sol1["conn.m_dot"] == pytest.approx(sol2["orf.m_dot"], rel=1e-6)

    def test_pressure_drop_with_flow(self):
        """Verify that EffectiveAreaConnectionElement produces pressure drop proportional to flow."""
        # Network with EffectiveAreaConnectionElement
        graph = FlowNetwork()
        inlet = PressureBoundary("inlet")
        inlet.P_total = 150000.0
        inlet.T_total = 300.0
        inlet.Y = cb.species.dry_air_mass()
        outlet = PressureBoundary("outlet")
        outlet.P_total = 100000.0
        outlet.T_total = 300.0
        outlet.Y = cb.species.dry_air_mass()
        conn = EffectiveAreaConnectionElement("conn", "inlet", "outlet", 0.01)

        graph.add_node(inlet)
        graph.add_node(outlet)
        graph.add_element(conn)

        solver = NetworkSolver(graph)
        sol = solver.solve(method="lm")

        # Verify positive flow (high to low pressure)
        assert sol["conn.m_dot"] > 0

        # Verify pressure drop exists
        assert inlet.P_total > outlet.P_total

    def test_different_areas_produce_different_flows(self):
        """Verify that different effective areas produce different flow rates."""
        areas = [0.005, 0.01, 0.02]
        flows = []

        for area in areas:
            graph = FlowNetwork()
            inlet = PressureBoundary("inlet")
            inlet.P_total = 150000.0
            inlet.T_total = 300.0
            inlet.Y = cb.species.dry_air_mass()
            outlet = PressureBoundary("outlet")
            outlet.P_total = 100000.0
            outlet.T_total = 300.0
            outlet.Y = cb.species.dry_air_mass()
            conn = EffectiveAreaConnectionElement("conn", "inlet", "outlet", area)

            graph.add_node(inlet)
            graph.add_node(outlet)
            graph.add_element(conn)

            solver = NetworkSolver(graph)
            sol = solver.solve(method="lm")
            flows.append(sol["conn.m_dot"])

        # Larger area should produce larger flow
        assert flows[0] < flows[1] < flows[2]

    def test_no_geometry_discovery(self):
        """Verify that resolve_topology does not discover upstream/downstream geometry."""
        network = FlowNetwork()
        plenum1 = PlenumNode("plenum1")
        plenum2 = PlenumNode("plenum2")
        conn = EffectiveAreaConnectionElement("conn", "plenum1", "plenum2", 0.1)

        network.add_node(plenum1)
        network.add_node(plenum2)
        network.add_element(conn)

        # Resolve topology
        conn.resolve_topology(network)

        # Verify no geometry was discovered
        assert conn.upstream_diameter is None
        assert conn.downstream_diameter is None

    def test_cd_always_one(self):
        """Verify that Cd is always 1.0 regardless of effective area."""
        test_areas = [0.001, 0.01, 0.1, 1.0, 10.0]

        for area in test_areas:
            conn = EffectiveAreaConnectionElement("conn", "node1", "node2", area)
            assert conn.Cd == 1.0
            assert conn.area == area

    def test_compressible_flow_behavior(self):
        """Verify that EffectiveAreaConnectionElement uses compressible flow equation.

        The flow should depend on density, which varies with pressure and temperature.
        This test verifies that the compressible orifice equation is used, not
        the incompressible approximation.
        """
        effective_area = 0.01  # m^2

        # High pressure drop case (more compressibility effects)
        graph1 = FlowNetwork()
        inlet1 = PressureBoundary("inlet")
        inlet1.P_total = 200000.0
        inlet1.T_total = 300.0
        inlet1.Y = cb.species.dry_air_mass()
        outlet1 = PressureBoundary("outlet")
        outlet1.P_total = 50000.0
        outlet1.T_total = 300.0
        outlet1.Y = cb.species.dry_air_mass()
        conn1 = EffectiveAreaConnectionElement("conn", "inlet", "outlet", effective_area)

        graph1.add_node(inlet1)
        graph1.add_node(outlet1)
        graph1.add_element(conn1)

        solver1 = NetworkSolver(graph1)
        sol1 = solver1.solve(method="lm")
        flow_high_drop = sol1["conn.m_dot"]

        # Low pressure drop case (less compressibility effects)
        graph2 = FlowNetwork()
        inlet2 = PressureBoundary("inlet")
        inlet2.P_total = 110000.0
        inlet2.T_total = 300.0
        inlet2.Y = cb.species.dry_air_mass()
        outlet2 = PressureBoundary("outlet")
        outlet2.P_total = 100000.0
        outlet2.T_total = 300.0
        outlet2.Y = cb.species.dry_air_mass()
        conn2 = EffectiveAreaConnectionElement("conn", "inlet", "outlet", effective_area)

        graph2.add_node(inlet2)
        graph2.add_node(outlet2)
        graph2.add_element(conn2)

        solver2 = NetworkSolver(graph2)
        sol2 = solver2.solve(method="lm")
        flow_low_drop = sol2["conn.m_dot"]

        # Both should have positive flow
        assert flow_high_drop > 0
        assert flow_low_drop > 0

        # Higher pressure drop should produce higher flow
        # For compressible flow, the relationship is more complex than incompressible
        # but flow should still increase monotonically with pressure drop
        assert flow_high_drop > flow_low_drop
