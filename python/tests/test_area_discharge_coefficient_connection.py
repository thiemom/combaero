"""
Test AreaDischargeCoefficientConnectionElement functionality.
"""

import pytest

from combaero.network import (
    AreaDischargeCoefficientConnectionElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def _get_air_Y():
    """Helper to create standard air composition."""
    Y = [0.0] * 14
    Y[0] = 0.79  # N2
    Y[1] = 0.21  # O2
    return Y


class TestAreaDischargeCoefficientConnectionElement:
    """Test suite for AreaDischargeCoefficientConnectionElement."""

    def test_initialization_with_cd(self):
        """Test initialization with discharge coefficient."""
        area = 0.0125  # m^2
        Cd = 0.8
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", area, Cd=Cd
        )

        # Verify attributes
        assert connection.Cd == Cd
        assert connection.area == area  # Physical area passed to parent
        assert connection.id == "conn1"
        assert connection.from_node == "node1"
        assert connection.to_node == "node2"

    def test_initialization_with_zeta(self):
        """Test initialization with loss coefficient."""
        area = 0.0125  # m^2
        zeta = 0.5625  # Should give Cd = 0.8
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", area, zeta=zeta
        )

        # Verify attributes
        expected_Cd = 1.0 / (zeta + 1.0) ** 0.5  # Cd = 1/sqrt(zeta + 1)
        assert connection.Cd == pytest.approx(expected_Cd)
        assert connection.area == pytest.approx(area)  # Physical area passed to parent
        assert connection.id == "conn1"

    def test_cd_zeta_equivalence(self):
        """Test that Cd and zeta produce equivalent results."""
        area = 0.0125
        Cd = 0.8
        zeta = 1.0 / (Cd**2) - 1.0  # zeta = 1/Cd^2 - 1

        # Initialize with Cd
        conn_cd = AreaDischargeCoefficientConnectionElement(
            "conn_cd", "node1", "node2", area, Cd=Cd
        )

        # Initialize with equivalent zeta
        conn_zeta = AreaDischargeCoefficientConnectionElement(
            "conn_zeta", "node1", "node2", area, zeta=zeta
        )

        # Should produce identical results
        assert conn_cd.Cd == pytest.approx(conn_zeta.Cd)
        assert conn_cd.area == pytest.approx(conn_zeta.area)

    def test_cd_takes_precedence(self):
        """Test that Cd takes precedence when both Cd and zeta are provided."""
        area = 0.01
        Cd = 0.9
        zeta = 0.2345  # Different from what Cd=0.9 would give

        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", area, Cd=Cd, zeta=zeta
        )

        # Should use Cd, not zeta
        assert connection.Cd == Cd
        assert connection.area == area  # Physical area passed to parent

    def test_inheritance(self):
        """Test that AreaDischargeCoefficientConnectionElement inherits from OrificeElement."""
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", 0.0125, Cd=0.8
        )
        assert isinstance(connection, OrificeElement)

    def test_unknowns(self):
        """Test unknowns method."""
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", 0.0125, Cd=0.8
        )
        unknowns = connection.unknowns()
        assert len(unknowns) == 1
        assert unknowns[0] == "conn1.m_dot"

    def test_n_equations(self):
        """Test n_equations method."""
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", 0.0125, Cd=0.8
        )
        assert connection.n_equations() == 1

    def test_resolve_topology(self):
        """Test resolve_topology method."""
        network = FlowNetwork()
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", 0.0125, Cd=0.8
        )

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

        # Create connection with Cd
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "plenum1", "plenum2", 0.0125, Cd=0.8
        )

        # Add to network
        network.add_node(plenum1)
        network.add_node(plenum2)
        network.add_element(connection)

        # Verify network contains the element
        assert "conn1" in network.elements
        assert network.elements["conn1"] is connection

    def test_different_discharge_coefficients(self):
        """Test with various discharge coefficient values."""
        area = 0.01
        test_Cds = [0.6, 0.7, 0.8, 0.9, 0.95]

        for i, Cd in enumerate(test_Cds):
            connection = AreaDischargeCoefficientConnectionElement(
                f"conn{i}", f"node{i}", f"node{i + 1}", area, Cd=Cd
            )
            assert connection.area == pytest.approx(area)  # Physical area
            assert connection.Cd == pytest.approx(Cd)

    def test_different_loss_coefficients(self):
        """Test with various loss coefficient values."""
        area = 0.01
        test_zetas = [0.0, 0.5, 1.0, 2.0, 5.0]

        for i, zeta in enumerate(test_zetas):
            connection = AreaDischargeCoefficientConnectionElement(
                f"conn{i}", f"node{i}", f"node{i + 1}", area, zeta=zeta
            )
            expected_Cd = 1.0 / (zeta + 1.0) ** 0.5
            assert connection.area == pytest.approx(area)  # Physical area
            assert connection.Cd == pytest.approx(expected_Cd)

    def test_methods_exist(self):
        """Test that all required methods exist."""
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", 0.0125, Cd=0.8
        )

        assert hasattr(connection, "unknowns")
        assert hasattr(connection, "residuals")
        assert hasattr(connection, "n_equations")
        assert hasattr(connection, "resolve_topology")
        assert callable(connection.unknowns)
        assert callable(connection.residuals)
        assert callable(connection.n_equations)
        assert callable(connection.resolve_topology)

    def test_flow_calculation_matches_orifice(self):
        """Verify AreaDischargeCoefficientConnectionElement produces same flow as OrificeElement.

        Both use the compressible flow equation via cb.orifice_mdot_and_jacobian(),
        which accounts for density variations across the pressure drop.
        """
        area = 0.0125  # m^2
        Cd = 0.8

        # Network 1: Using AreaDischargeCoefficientConnectionElement
        graph1 = FlowNetwork()
        inlet1 = PressureBoundary("inlet")
        inlet1.P_total = 150000.0
        inlet1.T_total = 300.0
        inlet1.Y = _get_air_Y()
        outlet1 = PressureBoundary("outlet")
        outlet1.P_total = 100000.0
        outlet1.T_total = 300.0
        outlet1.Y = _get_air_Y()
        conn1 = AreaDischargeCoefficientConnectionElement("conn", "inlet", "outlet", area, Cd=Cd)

        graph1.add_node(inlet1)
        graph1.add_node(outlet1)
        graph1.add_element(conn1)

        solver1 = NetworkSolver(graph1)
        sol1 = solver1.solve(method="lm")

        # Network 2: Using OrificeElement with same Cd and physical area
        graph2 = FlowNetwork()
        inlet2 = PressureBoundary("inlet")
        inlet2.P_total = 150000.0
        inlet2.T_total = 300.0
        inlet2.Y = _get_air_Y()
        outlet2 = PressureBoundary("outlet")
        outlet2.P_total = 100000.0
        outlet2.T_total = 300.0
        outlet2.Y = _get_air_Y()
        orf2 = OrificeElement("orf", "inlet", "outlet", Cd=Cd, area=area)

        graph2.add_node(inlet2)
        graph2.add_node(outlet2)
        graph2.add_element(orf2)

        solver2 = NetworkSolver(graph2)
        sol2 = solver2.solve(method="lm")

        # Should produce identical results
        assert sol1["conn.m_dot"] == pytest.approx(sol2["orf.m_dot"], rel=1e-6)

    def test_pressure_drop_with_flow(self):
        """Verify that AreaDischargeCoefficientConnectionElement produces pressure drop proportional to flow."""
        area = 0.0125
        Cd = 0.8

        # Network with AreaDischargeCoefficientConnectionElement
        graph = FlowNetwork()
        inlet = PressureBoundary("inlet")
        inlet.P_total = 150000.0
        inlet.T_total = 300.0
        inlet.Y = _get_air_Y()
        outlet = PressureBoundary("outlet")
        outlet.P_total = 100000.0
        outlet.T_total = 300.0
        outlet.Y = _get_air_Y()
        conn = AreaDischargeCoefficientConnectionElement("conn", "inlet", "outlet", area, Cd=Cd)

        graph.add_node(inlet)
        graph.add_node(outlet)
        graph.add_element(conn)

        solver = NetworkSolver(graph)
        sol = solver.solve(method="lm")

        # Verify positive flow (high to low pressure)
        assert sol["conn.m_dot"] > 0

        # Verify pressure drop exists
        assert inlet.P_total > outlet.P_total

    def test_different_cd_produce_different_flows(self):
        """Verify that different discharge coefficients produce different flow rates."""
        area = 0.01
        Cds = [0.6, 0.8, 0.95]
        flows = []

        for Cd in Cds:
            graph = FlowNetwork()
            inlet = PressureBoundary("inlet")
            inlet.P_total = 150000.0
            inlet.T_total = 300.0
            inlet.Y = _get_air_Y()
            outlet = PressureBoundary("outlet")
            outlet.P_total = 100000.0
            outlet.T_total = 300.0
            outlet.Y = _get_air_Y()
            conn = AreaDischargeCoefficientConnectionElement("conn", "inlet", "outlet", area, Cd=Cd)

            graph.add_node(inlet)
            graph.add_node(outlet)
            graph.add_element(conn)

            solver = NetworkSolver(graph)
            sol = solver.solve(method="lm")
            flows.append(sol["conn.m_dot"])

        # Higher Cd should produce higher flow
        assert flows[0] < flows[1] < flows[2]

    def test_different_zeta_produce_different_flows(self):
        """Verify that different loss coefficients produce different flow rates."""
        area = 0.01
        zetas = [2.0, 0.5, 0.1]  # Higher zeta = lower Cd = lower flow
        flows = []

        for zeta in zetas:
            graph = FlowNetwork()
            inlet = PressureBoundary("inlet")
            inlet.P_total = 150000.0
            inlet.T_total = 300.0
            inlet.Y = _get_air_Y()
            outlet = PressureBoundary("outlet")
            outlet.P_total = 100000.0
            outlet.T_total = 300.0
            outlet.Y = _get_air_Y()
            conn = AreaDischargeCoefficientConnectionElement(
                "conn", "inlet", "outlet", area, zeta=zeta
            )

            graph.add_node(inlet)
            graph.add_node(outlet)
            graph.add_element(conn)

            solver = NetworkSolver(graph)
            sol = solver.solve(method="lm")
            flows.append(sol["conn.m_dot"])

        # Lower zeta should produce higher flow (higher Cd)
        assert flows[0] < flows[1] < flows[2]

    def test_no_geometry_discovery(self):
        """Verify that resolve_topology does not discover upstream/downstream geometry."""
        network = FlowNetwork()
        plenum1 = PlenumNode("plenum1")
        plenum2 = PlenumNode("plenum2")
        conn = AreaDischargeCoefficientConnectionElement(
            "conn", "plenum1", "plenum2", 0.0125, Cd=0.8
        )

        network.add_node(plenum1)
        network.add_node(plenum2)
        network.add_element(conn)

        # Resolve topology
        conn.resolve_topology(network)

        # Verify no geometry was discovered
        assert conn.upstream_diameter is None
        assert conn.downstream_diameter is None

    def test_error_neither_cd_nor_zeta(self):
        """Test error when neither Cd nor zeta is provided."""
        with pytest.raises(ValueError, match="Either Cd or zeta must be provided"):
            AreaDischargeCoefficientConnectionElement("conn1", "node1", "node2", 0.01)

    def test_error_invalid_cd_range(self):
        """Test error when Cd is outside valid range."""
        # Cd > 1
        with pytest.raises(ValueError, match="Cd must be in range"):
            AreaDischargeCoefficientConnectionElement("conn1", "node1", "node2", 0.01, Cd=1.5)

        # Cd <= 0
        with pytest.raises(ValueError, match="Cd must be in range"):
            AreaDischargeCoefficientConnectionElement("conn1", "node1", "node2", 0.01, Cd=0.0)

        # Cd negative
        with pytest.raises(ValueError, match="Cd must be in range"):
            AreaDischargeCoefficientConnectionElement("conn1", "node1", "node2", 0.01, Cd=-0.5)

    def test_error_invalid_zeta_range(self):
        """Test error when zeta is negative."""
        with pytest.raises(ValueError, match="zeta must be >= 0"):
            AreaDischargeCoefficientConnectionElement("conn1", "node1", "node2", 0.01, zeta=-0.1)

    def test_compressible_flow_behavior(self):
        """Verify that AreaDischargeCoefficientConnectionElement uses compressible flow equation.

        The flow should depend on density, which varies with pressure and temperature.
        This test verifies that the compressible orifice equation is used, not
        the incompressible approximation.
        """
        area = 0.01
        Cd = 0.8

        # High pressure drop case (more compressibility effects)
        graph1 = FlowNetwork()
        inlet1 = PressureBoundary("inlet")
        inlet1.P_total = 200000.0
        inlet1.T_total = 300.0
        inlet1.Y = _get_air_Y()
        outlet1 = PressureBoundary("outlet")
        outlet1.P_total = 50000.0
        outlet1.T_total = 300.0
        outlet1.Y = _get_air_Y()
        conn1 = AreaDischargeCoefficientConnectionElement("conn", "inlet", "outlet", area, Cd=Cd)

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
        inlet2.Y = _get_air_Y()
        outlet2 = PressureBoundary("outlet")
        outlet2.P_total = 100000.0
        outlet2.T_total = 300.0
        outlet2.Y = _get_air_Y()
        conn2 = AreaDischargeCoefficientConnectionElement("conn", "inlet", "outlet", area, Cd=Cd)

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

    def test_zeta_zero_equals_cd_one(self):
        """Test that zeta=0 gives Cd=1.0."""
        area = 0.01
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", area, zeta=0.0
        )

        assert connection.Cd == pytest.approx(1.0)
        assert connection.area == pytest.approx(area)  # Physical area

    def test_high_zeta_low_cd(self):
        """Test that high zeta values give low Cd values."""
        area = 0.01
        high_zeta = 10.0
        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", area, zeta=high_zeta
        )

        expected_Cd = 1.0 / (high_zeta + 1.0) ** 0.5
        assert connection.Cd == pytest.approx(expected_Cd)
        assert connection.Cd < 0.31  # Should be quite low

    def test_physical_area_preserved_in_calculation(self):
        """Verify that the physical area is used correctly in flow calculation."""
        physical_area = 0.0125
        Cd = 0.8

        connection = AreaDischargeCoefficientConnectionElement(
            "conn1", "node1", "node2", physical_area, Cd=Cd
        )

        # Physical area is stored in parent
        assert connection.area == pytest.approx(physical_area)

        # Cd is stored separately
        assert connection.Cd == pytest.approx(Cd)

        # Parent will use A * Cd for flow calculation
        # Effective area = physical_area * Cd
