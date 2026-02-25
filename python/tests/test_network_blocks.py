import combaero as cb
from combaero.network.components import (
    BoundaryNode,
    MixtureState,
    OrificeElement,
    PipeElement,
    PlenumNode,
)


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
