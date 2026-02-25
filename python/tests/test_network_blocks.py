import itertools

import pytest

import combaero as cb
from combaero.network.components import (
    BoundaryNode,
    LosslessConnectionElement,
    MixtureState,
    OrificeElement,
    PipeElement,
    PlenumNode,
)

# ==============================================================================
# COMBINATORIAL TEST HARNESS CONFIGURATION
# ==============================================================================
NODE_FACTORIES = {
    "plenum": lambda id: PlenumNode(id),
    "momentum_chamber": lambda id: cb.network.MomentumChamberNode(id),
    "boundary": lambda id: BoundaryNode(id),
}

ELEMENT_FACTORIES = {
    "pipe": lambda id, f, t: PipeElement(id, f, t, length=2.0, diameter=0.15, roughness=1e-5),
    "orifice": lambda id, f, t: OrificeElement(id, f, t, Cd=0.6, area=0.005),
    "lossless": lambda id, f, t: LosslessConnectionElement(id, f, t),
}

PERMUTATIONS = list(
    itertools.product(NODE_FACTORIES.keys(), ELEMENT_FACTORIES.keys(), NODE_FACTORIES.keys())
)
# ==============================================================================


def test_atomic_network_components():
    """
    Test 1: Simple Boundary -> Pipe -> Orifice -> Sink configuration.
    Currently only validates instantiation and topological property storage.
    """
    # 1. Initialize states
    X_air = [0.0] * 14
    X_air[11] = 0.79  # N2
    X_air[12] = 0.21  # O2

    state_in = MixtureState(
        P=500000.0, P_total=500000.0, T=300.0, T_total=300.0, m_dot=0.1, X=X_air
    )

    _state_out = MixtureState(
        P=101325.0, P_total=101325.0, T=300.0, T_total=300.0, m_dot=0.1, X=X_air
    )

    # 2. Instantiate nodes (Boundary, intermediate Plenum, Sink)
    inlet = BoundaryNode("inlet")
    node_1 = PlenumNode("node_1")
    _outlet = BoundaryNode("outlet")

    # 3. Instantiate elements
    pipe_1 = PipeElement(
        id="pipe_1", from_node="inlet", to_node="node_1", length=2.0, diameter=0.05, roughness=1e-5
    )

    orifice_1 = OrificeElement(
        id="orifice_1", from_node="node_1", to_node="outlet", Cd=0.6, area=0.005
    )

    # 4. Verify basic properties
    assert pipe_1.area > 0.0
    assert orifice_1.n_equations() == 1
    assert pipe_1.n_equations() == 1
    assert len(inlet.unknowns()) == 0
    assert len(node_1.unknowns()) == 1

    # 5. Verify mixture state wrappers
    rho = state_in.density()
    assert abs(rho - cb.density(300.0, 500000.0, X_air)) < 1e-9


def test_momentum_chamber_network():
    """
    Test 2: Momentum Chamber -> Pipe -> Orifice.
    Validates instantiation of the MomentumChamberNode and its topological footprint.
    """
    node_mc = cb.network.MomentumChamberNode("chamber_1")

    # 1. Verify MomentumChamber properties
    unknowns = node_mc.unknowns()
    assert len(unknowns) == 2
    assert "chamber_1.P" in unknowns
    assert "chamber_1.P_total" in unknowns

    # 1b. Verify directional flow vectors (angles)
    assert node_mc.get_port_angle("pipe_in") == 0.0  # default
    node_mc.set_port_angle("pipe_in", 45.0)
    assert node_mc.get_port_angle("pipe_in") == 45.0

    # Placeholder residuals ensure conservation equations can later be injected
    state = cb.network.MixtureState(1e5, 1.1e5, 300.0, 300.0, 1.0, [0.0] * 14)
    residuals = node_mc.residuals(state)
    assert len(residuals) == 1
    assert residuals[0] == 1.1e5 - 1e5


def test_flownetwork_topology():
    """
    Test 3: Validates the FlowNetwork topological manager and auto-discovery.
    Specifically checks that OrificeElement finds its upstream PipeElement diameter.
    """
    graph = cb.network.FlowNetwork()

    inlet = cb.network.BoundaryNode("inlet")
    node_1 = cb.network.PlenumNode("node_1")
    outlet = cb.network.BoundaryNode("outlet")

    pipe = cb.network.PipeElement(
        id="pipe_1", from_node="inlet", to_node="node_1", length=2.0, diameter=0.15, roughness=1e-5
    )

    orifice = cb.network.OrificeElement(
        id="orf_1", from_node="node_1", to_node="outlet", Cd=0.6, area=0.005
    )

    # 1. Register graph
    graph.add_node(inlet)
    graph.add_node(node_1)
    graph.add_node(outlet)

    graph.add_element(pipe)
    graph.add_element(orifice)

    # 2. Verify pre-resolution state
    assert orifice.upstream_diameter is None

    # 3. Resolve topology
    graph.resolve_all_topology()

    # 4. Verify post-resolution auto-discovery
    assert orifice.upstream_diameter == 0.15


@pytest.mark.parametrize("n_upstream_key, el_type_key, n_downstream_key", PERMUTATIONS)
def test_all_topological_combinations(n_upstream_key, el_type_key, n_downstream_key):
    """
    Test 5: Scalable Matrix Runner.
    Dynamically tests every registered Network Node against every Network Element.
    Ensures safe instantiation and graph topology resolution across all possible architectures.
    """
    graph = cb.network.FlowNetwork()

    # 1. Instantiate endpoints dynamically from factory keys
    n_up = NODE_FACTORIES[n_upstream_key]("node_up")
    n_down = NODE_FACTORIES[n_downstream_key]("node_down")

    # 2. Instantiate intermediate element dynamically
    element = ELEMENT_FACTORIES[el_type_key]("element_mid", "node_up", "node_down")

    # 3. Assemble generic graph
    graph.add_node(n_up)
    graph.add_node(n_down)
    graph.add_element(element)

    # 4. Run global physics interface tests
    assert isinstance(element.n_equations(), int)
    assert isinstance(n_up.unknowns(), list)

    # 5. Verify Topology doesn't crash on unfamiliar combinations
    try:
        graph.resolve_all_topology()
    except Exception as e:
        pytest.fail(
            f"Topology auto-discovery crashed on {n_upstream_key} -> {el_type_key} -> {n_downstream_key}: {str(e)}"
        )


def test_flownetwork_multibranch():
    """
    Test 4: Validates that Plenums/Chambers can accept multiple inlets and outlets.
    Verifies that Orifice auto-discovery safely aborts when upstream geometry is ambiguous.
    """
    graph = cb.network.FlowNetwork()

    # Nodes
    chamber = cb.network.MomentumChamberNode("chamber_main")
    graph.add_node(chamber)
    for i in range(1, 4):
        graph.add_node(cb.network.BoundaryNode(f"bnd_{i}"))

    # Elements: Two pipes entering the chamber
    pipe_a = cb.network.PipeElement("pipe_a", "bnd_1", "chamber_main", 1.0, 0.1, 1e-5)
    pipe_b = cb.network.PipeElement("pipe_b", "bnd_2", "chamber_main", 1.0, 0.2, 1e-5)

    # Elements: One orifice leaving the chamber
    orf_out = cb.network.OrificeElement("orf_out", "chamber_main", "bnd_3", 0.6, 0.05)

    graph.add_element(pipe_a)
    graph.add_element(pipe_b)
    graph.add_element(orf_out)

    # 1. Verify Topology Array Lengths
    upstreams = graph.get_upstream_elements("chamber_main")
    downstreams = graph.get_downstream_elements("chamber_main")

    assert len(upstreams) == 2
    assert len(downstreams) == 1
    assert "pipe_a" in [e.id for e in upstreams]
    assert "pipe_b" in [e.id for e in upstreams]

    # 2. Verify Auto-Discovery Safety
    # The orifice looks upstream, sees TWO pipes feeding the chamber,
    # and should safely abort extracting a diameter instead of crashing or picking wrong.
    assert orf_out.upstream_diameter is None
    graph.resolve_all_topology()
    assert orf_out.upstream_diameter is None
