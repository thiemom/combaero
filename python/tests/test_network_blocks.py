import itertools
import os
import random

import pytest

import combaero as cb
from combaero.network.components import (
    ChannelElement,
    LosslessConnectionElement,
    MixtureState,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

# ==============================================================================
# COMBINATORIAL TEST HARNESS CONFIGURATION
# ==============================================================================
NODE_FACTORIES = {
    "plenum": lambda id: PlenumNode(id),
    "momentum_chamber": lambda id: cb.network.MomentumChamberNode(id),
    "boundary": lambda id: PressureBoundary(id),
    "combustor_complete": lambda id: cb.network.CombustorNode(id, method="complete"),
    "combustor_eq": lambda id: cb.network.CombustorNode(id, method="equilibrium"),
}

ELEMENT_FACTORIES = {
    "channel": lambda id, f, t: ChannelElement(id, f, t, length=2.0, diameter=0.15, roughness=1e-5),
    "channel_fanno": lambda id, f, t: ChannelElement(
        id,
        f,
        t,
        length=2.0,
        diameter=0.15,
        roughness=1e-5,
        regime="compressible",
        friction_model="colebrook",
    ),
    "orifice": lambda id, f, t: OrificeElement(id, f, t, Cd=0.6, diameter=0.079788),
    "lossless": lambda id, f, t: LosslessConnectionElement(id, f, t),
}

PERMUTATIONS = list(
    itertools.product(NODE_FACTORIES.keys(), ELEMENT_FACTORIES.keys(), NODE_FACTORIES.keys())
)

MAX_TOPOLOGY_TESTS = int(os.environ.get("MAX_TOPOLOGY_TESTS", "50"))
if len(PERMUTATIONS) > MAX_TOPOLOGY_TESTS:
    random.seed(int(os.environ.get("TOPOLOGY_TEST_SEED", "42")))
    PERMUTATIONS = random.sample(PERMUTATIONS, MAX_TOPOLOGY_TESTS)
# ==============================================================================


def test_atomic_network_components():
    """
    Test 1: Simple Boundary -> Channel -> Orifice -> Sink configuration.
    Currently only validates instantiation and topological property storage.
    """
    # 1. Initialize states
    Y_air = [0.0] * 14
    Y_air[11] = 0.79  # N2
    Y_air[12] = 0.21  # O2

    state_in = MixtureState(
        P=500000.0, P_total=500000.0, T=300.0, T_total=300.0, m_dot=0.1, Y=Y_air
    )

    _state_out = MixtureState(
        P=101325.0, P_total=101325.0, T=300.0, T_total=300.0, m_dot=0.1, Y=Y_air
    )

    # 2. Instantiate nodes (Boundary, intermediate Plenum, Sink)
    inlet = PressureBoundary("inlet")
    node_1 = PlenumNode("node_1")
    _outlet = PressureBoundary("outlet")

    # 3. Instantiate elements
    channel_1 = ChannelElement(
        id="channel_1",
        from_node="inlet",
        to_node="node_1",
        length=2.0,
        diameter=0.05,
        roughness=1e-5,
    )

    orifice_1 = OrificeElement(
        id="orifice_1", from_node="node_1", to_node="outlet", Cd=0.6, diameter=0.079788
    )

    # 4. Verify basic properties
    assert channel_1.area > 0.0
    assert orifice_1.n_equations() == 1
    assert channel_1.n_equations() == 1
    assert len(inlet.unknowns()) == 0
    assert len(node_1.unknowns()) >= 2  # P, P_total, and possibly T or X

    # 5. Verify mixture state wrappers
    rho = state_in.density()
    assert abs(rho - cb.density(300.0, 500000.0, cb.mass_to_mole(Y_air))) < 1e-9


def test_momentum_chamber_network():
    """
    Test 2: Momentum Chamber -> Channel -> Orifice.
    Validates instantiation of the MomentumChamberNode and its topological footprint.
    """
    node_mc = cb.network.MomentumChamberNode("chamber_1", area=0.1)

    # 1. Verify MomentumChamber properties
    unknowns = node_mc.unknowns()
    assert len(unknowns) >= 2
    assert "chamber_1.P" in unknowns
    assert "chamber_1.P_total" in unknowns

    # 1b. Verify directional flow vectors (angles)
    assert node_mc.get_port_angle("channel_in") == 0.0  # default
    node_mc.set_port_angle("channel_in", 45.0)
    assert node_mc.get_port_angle("channel_in") == 45.0

    # Placeholder residuals ensure conservation equations can later be injected
    state = cb.network.MixtureState(1e5, 1.1e5, 300.0, 300.0, 1.0, [0.0] * 14)
    residuals, _ = node_mc.residuals(state)
    assert len(residuals) == 1
    assert residuals[0] == 1.1e5 - 1e5


def test_flownetwork_topology():
    """
    Test 3: Validates the FlowNetwork topological manager and auto-discovery.
    Specifically checks that OrificeElement finds its upstream ChannelElement diameter.
    """
    graph = cb.network.FlowNetwork()

    inlet = cb.network.PressureBoundary("inlet")
    node_1 = cb.network.PlenumNode("node_1")
    outlet = cb.network.PressureBoundary("outlet")

    channel = cb.network.ChannelElement(
        id="channel_1",
        from_node="inlet",
        to_node="node_1",
        length=2.0,
        diameter=0.15,
        roughness=1e-5,
    )

    orifice = cb.network.OrificeElement(
        id="orf_1", from_node="node_1", to_node="outlet", Cd=0.6, diameter=0.079788
    )

    # 1. Register graph
    graph.add_node(inlet)
    graph.add_node(node_1)
    graph.add_node(outlet)

    graph.add_element(channel)
    graph.add_element(orifice)

    # 2. Verify pre-resolution state
    assert orifice.upstream_diameter is None

    # 3. Resolve topology
    graph.resolve_all_topology()

    # 4. Verify post-resolution auto-discovery
    assert orifice.upstream_diameter == 0.15


def test_flownetwork_validation():
    """
    Test 6: Validates that FlowNetwork.validate() catches ill-posed networks.
    """
    # 1. No PressureBoundary -> ValueError
    graph1 = cb.network.FlowNetwork()
    graph1.add_node(cb.network.MassFlowBoundary("bnd_mass"))
    graph1.add_node(cb.network.PlenumNode("plenum_1"))

    # We must add an outlet to plenum_1 otherwise it fails early on interior connection < 2 rule
    graph1.add_node(cb.network.MassFlowBoundary("bnd_mass_out"))

    graph1.add_element(
        cb.network.ChannelElement("channel_1", "bnd_mass", "plenum_1", 1.0, 0.1, 1e-5)
    )
    graph1.add_element(
        cb.network.ChannelElement("channel_2", "plenum_1", "bnd_mass_out", 1.0, 0.1, 1e-5)
    )

    with pytest.raises(ValueError, match="at least one PressureBoundary"):
        graph1.validate()

    # 2. Add PressureBoundary but network has no losses -> ValueError
    graph2 = cb.network.FlowNetwork()
    bnd_press_in = cb.network.PressureBoundary("bnd_press_in")
    plenum = cb.network.PlenumNode("plenum_1")
    bnd_press_out = cb.network.PressureBoundary("bnd_press_out")

    graph2.add_node(plenum)
    graph2.add_node(bnd_press_in)
    graph2.add_node(bnd_press_out)

    graph2.add_element(
        cb.network.LosslessConnectionElement("lossless1", "bnd_press_in", "plenum_1")
    )
    graph2.add_element(
        cb.network.LosslessConnectionElement("lossless2", "plenum_1", "bnd_press_out")
    )

    with pytest.raises(ValueError, match="at least one pressure drop element"):
        graph2.validate()

    # 3. Add an element with loss -> Validates successfully
    graph3 = cb.network.FlowNetwork()
    graph3.add_node(cb.network.PlenumNode("plenum_1"))
    graph3.add_node(cb.network.PressureBoundary("bnd_press_in"))
    graph3.add_node(cb.network.PressureBoundary("bnd_press_out"))

    graph3.add_element(cb.network.OrificeElement("orifice", "bnd_press_in", "plenum_1", 0.6, 0.05))
    graph3.add_element(
        cb.network.OrificeElement("orifice2", "plenum_1", "bnd_press_out", 0.6, 0.05)
    )

    graph3.validate()  # Should not raise
