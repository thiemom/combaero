"""Tests for chamber mixing functionality in network solver."""

from __future__ import annotations

import combaero as cb
from combaero.network.components import (
    CombustorNode,
    MassFlowBoundary,
    MixtureState,
    MomentumChamberNode,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.graph import FlowNetwork


def test_plenum_node_creation():
    """Test PlenumNode creation and basic properties."""
    plenum = PlenumNode("test_plenum")

    # Should have P and P_total unknowns (Derived-State architecture)
    unknowns = plenum.unknowns()
    assert "test_plenum.P" in unknowns
    assert "test_plenum.P_total" in unknowns
    assert "test_plenum.T" not in unknowns
    assert len(unknowns) == 2


def test_momentum_chamber_node_creation():
    """Test MomentumChamberNode creation and basic properties."""
    chamber = MomentumChamberNode("test_chamber", area=0.05)

    # Should have P and P_total unknowns (Derived-State architecture)
    unknowns = chamber.unknowns()
    assert "test_chamber.P" in unknowns
    assert "test_chamber.P_total" in unknowns
    assert "test_chamber.T" not in unknowns
    assert len(unknowns) == 2


def test_combustor_node_creation():
    """Test CombustorNode creation and basic properties."""
    combustor = CombustorNode("test_combustor", method="complete")

    # Should have P and P_total unknowns (Derived-State architecture)
    unknowns = combustor.unknowns()
    assert "test_combustor.P" in unknowns
    assert "test_combustor.P_total" in unknowns
    assert "test_combustor.T" not in unknowns
    assert len(unknowns) == 2

    assert combustor.method == "complete"


def test_combustor_fuel_boundary():
    """Test setting fuel boundary on combustor."""
    combustor = CombustorNode("test_combustor")

    # Create a fuel boundary
    fuel_bc = MassFlowBoundary(
        "fuel_injector",
        m_dot=0.02,
        T_total=300.0,
        Y=cb.species.dry_air_mass(),
    )

    # Set fuel boundary
    combustor.set_fuel_boundary(fuel_bc)

    assert combustor.fuel_boundary is not None
    assert combustor.fuel_boundary.id == "fuel_injector"


def test_chamber_topology_resolution():
    """Test that chambers properly resolve their topology."""
    # Create a network with chambers
    network = FlowNetwork()

    # Add nodes
    plenum = PlenumNode("plenum1")
    chamber = MomentumChamberNode("chamber1", area=0.1)
    combustor = CombustorNode("combustor1")

    network.add_node(plenum)
    network.add_node(chamber)
    network.add_node(combustor)

    # Resolve topology (this should be called automatically, but we'll test it explicitly)
    plenum.resolve_topology(network)
    chamber.resolve_topology(network)
    combustor.resolve_topology(network)

    # Initially, no upstream elements
    assert len(plenum.upstream_elements) == 0
    assert len(chamber.upstream_elements) == 0
    assert len(combustor.upstream_elements) == 0


def test_plenum_with_multiple_upstream():
    """Test that plenum adds species unknowns when it has multiple upstream connections."""
    plenum = PlenumNode("mixing_plenum")

    # Simulate multiple upstream connections by adding to upstream_elements
    plenum.upstream_elements = ["elem1", "elem2"]  # Two upstream elements

    # Now should only have P and P_total unknowns (Pure Pressure-Flow)
    unknowns = plenum.unknowns()
    assert "mixing_plenum.P" in unknowns
    assert "mixing_plenum.P_total" in unknowns
    assert "mixing_plenum.T" not in unknowns
    assert len([u for u in unknowns if "Y[" in u]) == 0
    assert len(unknowns) == 2

    # Verification of species as global unknowns is no longer applicable in Phase 1


def test_chamber_pressure_relationship():
    """Test that chamber nodes enforce proper pressure relationships."""
    # Test PlenumNode: P_total = P_static (Legacy path)
    plenum = PlenumNode("plenum")
    state = MixtureState(
        P=100000.0,
        P_total=100000.0,
        T=300.0,
        T_total=300.0,
        m_dot=0.1,
        Y=cb.species.dry_air_mass(),
    )

    residuals, jac = plenum.residuals(state)

    # Should have one residual: P_total - P = 0 (Legacy path)
    assert len(residuals) == 1
    assert residuals[0] == 0.0  # P_total - P = 100000 - 100000 = 0

    # Test with mismatch (Legacy path)
    state_mismatch = MixtureState(
        P=100000.0,
        P_total=101000.0,
        T=300.0,
        T_total=300.0,
        m_dot=0.1,
        Y=cb.species.dry_air_mass(),
    )
    residuals_mismatch, jac_mismatch = plenum.residuals(state_mismatch)

    assert len(residuals_mismatch) == 1
    assert residuals_mismatch[0] == 1000.0  # P_total - P = 101000 - 100000 = 1000

    # Removed test_momentum_chamber_velocity because 'v' is no longer an unknown


def test_combustor_methods():
    """Test different combustion methods."""
    # Test complete combustion method
    combustor_complete = CombustorNode("combustor1", method="complete")
    assert combustor_complete.method == "complete"

    # Test equilibrium method
    combustor_equilibrium = CombustorNode("combustor2", method="equilibrium")
    assert combustor_equilibrium.method == "equilibrium"

    # Pressure loss fraction tests removed (dead code)


def test_chamber_integration_with_network():
    """Test that chambers can be integrated into a network."""
    network = FlowNetwork()

    # Create nodes
    inlet = PressureBoundary("inlet", P_total=500000.0, T_total=700.0)
    plenum = PlenumNode("plenum")
    chamber = MomentumChamberNode("chamber", area=0.1)
    combustor = CombustorNode("combustor")
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=400.0)

    # Add to network
    network.add_node(inlet)
    network.add_node(plenum)
    network.add_node(chamber)
    network.add_node(combustor)
    network.add_node(outlet)

    # Verify all nodes are in the network
    assert "inlet" in network.nodes
    assert "plenum" in network.nodes
    assert "chamber" in network.nodes
    assert "combustor" in network.nodes
    assert "outlet" in network.nodes

    # Verify node types
    assert isinstance(network.nodes["plenum"], PlenumNode)
    assert isinstance(network.nodes["chamber"], MomentumChamberNode)
    assert isinstance(network.nodes["combustor"], CombustorNode)
