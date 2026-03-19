"""Test FlowNetwork wall management functionality."""

import pytest

from combaero.network import FlowNetwork, PipeElement, PlenumNode, PressureBoundary, WallConnection


def test_flow_network_defaults():
    """Test FlowNetwork default values."""
    network = FlowNetwork()

    # Default values
    assert network.thermal_coupling_enabled is True
    assert network.walls == {}


def test_add_wall_success():
    """Test successful wall addition."""
    network = FlowNetwork()

    # Add nodes and elements
    pb = PressureBoundary("pb", P_total=101325, T_total=300)
    plenum = PlenumNode("plenum")
    pipe1 = PipeElement("pipe1", "pb", "plenum", length=1.0, diameter=0.05, roughness=1e-5)
    pipe2 = PipeElement("pipe2", "pb", "plenum", length=1.0, diameter=0.05, roughness=1e-5)

    network.add_node(pb)
    network.add_node(plenum)
    network.add_element(pipe1)
    network.add_element(pipe2)

    # Add wall
    wall = WallConnection(
        id="wall1",
        element_a="pipe1",
        element_b="pipe2",
        wall_thickness=0.002,
        wall_conductivity=25.0,
    )
    network.add_wall(wall)

    assert len(network.walls) == 1
    assert "wall1" in network.walls
    assert network.walls["wall1"] is wall


def test_add_wall_duplicate_id():
    """Test adding wall with duplicate ID."""
    network = FlowNetwork()

    # Add nodes and elements
    pb = PressureBoundary("pb", P_total=101325, T_total=300)
    plenum = PlenumNode("plenum")
    pipe1 = PipeElement("pipe1", "pb", "plenum", length=1.0, diameter=0.05, roughness=1e-5)

    network.add_node(pb)
    network.add_node(plenum)
    network.add_element(pipe1)

    # Add first wall
    wall1 = WallConnection("wall1", "pipe1", "pipe1", 0.002, 25.0)
    network.add_wall(wall1)

    # Try to add duplicate
    wall1_dup = WallConnection("wall1", "pipe1", "pipe1", 0.002, 25.0)
    with pytest.raises(ValueError, match="Wall 'wall1' already exists in network"):
        network.add_wall(wall1_dup)


def test_add_wall_unknown_element_a():
    """Test adding wall with unknown element_a."""
    network = FlowNetwork()

    # Add nodes and elements
    pb = PressureBoundary("pb", P_total=101325, T_total=300)
    plenum = PlenumNode("plenum")
    pipe1 = PipeElement("pipe1", "pb", "plenum", length=1.0, diameter=0.05, roughness=1e-5)

    network.add_node(pb)
    network.add_node(plenum)
    network.add_element(pipe1)

    # Try to add wall with unknown element_a
    wall = WallConnection("wall2", "unknown_pipe", "pipe1", 0.002, 25.0)
    with pytest.raises(ValueError, match="Unknown element_a 'unknown_pipe' for wall 'wall2'"):
        network.add_wall(wall)


def test_add_wall_unknown_element_b():
    """Test adding wall with unknown element_b."""
    network = FlowNetwork()

    # Add nodes and elements
    pb = PressureBoundary("pb", P_total=101325, T_total=300)
    plenum = PlenumNode("plenum")
    pipe1 = PipeElement("pipe1", "pb", "plenum", length=1.0, diameter=0.05, roughness=1e-5)

    network.add_node(pb)
    network.add_node(plenum)
    network.add_element(pipe1)

    # Try to add wall with unknown element_b
    wall = WallConnection("wall3", "pipe1", "unknown_pipe", 0.002, 25.0)
    with pytest.raises(ValueError, match="Unknown element_b 'unknown_pipe' for wall 'wall3'"):
        network.add_wall(wall)


def test_thermal_coupling_toggle():
    """Test thermal coupling enabled toggle."""
    network = FlowNetwork()

    # Default is True
    assert network.thermal_coupling_enabled is True

    # Can be disabled
    network.thermal_coupling_enabled = False
    assert network.thermal_coupling_enabled is False

    # Can be re-enabled
    network.thermal_coupling_enabled = True
    assert network.thermal_coupling_enabled is True


def test_to_dict_with_walls():
    """Test serialization includes walls and thermal coupling."""
    network = FlowNetwork()

    # Add nodes and elements
    pb = PressureBoundary("pb", P_total=101325, T_total=300)
    plenum = PlenumNode("plenum")
    pipe1 = PipeElement("pipe1", "pb", "plenum", length=1.0, diameter=0.05, roughness=1e-5)
    pipe2 = PipeElement("pipe2", "pb", "plenum", length=1.0, diameter=0.05, roughness=1e-5)

    network.add_node(pb)
    network.add_node(plenum)
    network.add_element(pipe1)
    network.add_element(pipe2)

    # Add wall
    wall = WallConnection(
        id="wall1",
        element_a="pipe1",
        element_b="pipe2",
        wall_thickness=0.002,
        wall_conductivity=25.0,
        contact_area=0.05,
    )
    network.add_wall(wall)

    # Disable thermal coupling
    network.thermal_coupling_enabled = False

    # Serialize
    data = network.to_dict()

    # Check structure
    assert "walls" in data
    assert "thermal_coupling_enabled" in data
    assert data["thermal_coupling_enabled"] is False
    assert len(data["walls"]) == 1
    assert "wall1" in data["walls"]

    # Check wall data
    wall_data = data["walls"]["wall1"]
    assert wall_data["type"] == "WallConnection"
    assert wall_data["kwargs"]["element_a"] == "pipe1"
    assert wall_data["kwargs"]["element_b"] == "pipe2"
    assert wall_data["kwargs"]["wall_thickness"] == 0.002
    assert wall_data["kwargs"]["wall_conductivity"] == 25.0
    assert wall_data["kwargs"]["contact_area"] == 0.05


def test_from_dict_with_walls():
    """Test deserialization restores walls and thermal coupling."""
    # Create original network
    network1 = FlowNetwork()

    # Add nodes and elements
    pb = PressureBoundary("pb", P_total=101325, T_total=300)
    plenum = PlenumNode("plenum")
    pipe1 = PipeElement("pipe1", "pb", "plenum", length=1.0, diameter=0.05, roughness=1e-5)
    pipe2 = PipeElement("pipe2", "pb", "plenum", length=1.0, diameter=0.05, roughness=1e-5)

    network1.add_node(pb)
    network1.add_node(plenum)
    network1.add_element(pipe1)
    network1.add_element(pipe2)

    # Add wall
    wall = WallConnection(
        id="wall1",
        element_a="pipe1",
        element_b="pipe2",
        wall_thickness=0.002,
        wall_conductivity=25.0,
    )
    network1.add_wall(wall)

    # Disable thermal coupling
    network1.thermal_coupling_enabled = False

    # Serialize and deserialize
    data = network1.to_dict()
    network2 = FlowNetwork.from_dict(data)

    # Check restoration
    assert network2.thermal_coupling_enabled is False
    assert len(network2.walls) == 1
    assert "wall1" in network2.walls

    restored_wall = network2.walls["wall1"]
    assert restored_wall.element_a == "pipe1"
    assert restored_wall.element_b == "pipe2"
    assert restored_wall.wall_thickness == 0.002
    assert restored_wall.wall_conductivity == 25.0


def test_from_dict_default_thermal_coupling():
    """Test from_dict uses default thermal coupling when not specified."""
    data = {
        "nodes": {},
        "elements": {},
        "walls": {},
        # No thermal_coupling_enabled key
    }

    network = FlowNetwork.from_dict(data)
    assert network.thermal_coupling_enabled is True  # Default value
