from gui.backend.graph_builder import build_network_from_schema
from gui.backend.schemas import NetworkGraphSchema


def test_plenum_thermal_edges_are_ignored_silently():
    schema = NetworkGraphSchema(
        nodes=[
            {
                "id": "plenum_in",
                "type": "plenum",
                "position": {"x": 0.0, "y": 0.0},
                "data": {},
            },
            {
                "id": "plenum_out",
                "type": "plenum",
                "position": {"x": 400.0, "y": 0.0},
                "data": {},
            },
            {
                "id": "channel_1",
                "type": "channel",
                "position": {"x": 200.0, "y": 0.0},
                "data": {"L": 1.0, "D": 0.1},
            },
        ],
        edges=[
            {
                "id": "edge_flow_1",
                "source": "plenum_in",
                "target": "channel_1",
                "data": {},
            },
            {
                "id": "edge_flow_2",
                "source": "channel_1",
                "target": "plenum_out",
                "data": {},
            },
            {
                "id": "edge_thermal_plenum_channel",
                "source": "plenum_in",
                "target": "channel_1",
                "type": "thermal",
                "data": {"type": "thermal", "thickness": 0.003, "conductivity": 20.0},
            },
        ],
    )

    net = build_network_from_schema(schema)

    assert len(net.walls) == 0


def test_momentum_chamber_thermal_wall_is_built_and_multipliers_forwarded():
    schema = NetworkGraphSchema(
        nodes=[
            {
                "id": "plenum_in",
                "type": "plenum",
                "position": {"x": 0.0, "y": 0.0},
                "data": {},
            },
            {
                "id": "mom_1",
                "type": "momentum_chamber",
                "position": {"x": 200.0, "y": 0.0},
                "data": {"area": 0.2, "Nu_multiplier": 1.3, "f_multiplier": 0.8},
            },
            {
                "id": "channel_1",
                "type": "channel",
                "position": {"x": 400.0, "y": 0.0},
                "data": {"L": 1.0, "D": 0.1, "Nu_multiplier": 1.2, "f_multiplier": 0.9},
            },
        ],
        edges=[
            {
                "id": "edge_flow_1",
                "source": "plenum_in",
                "target": "channel_1",
                "data": {},
            },
            {
                "id": "edge_flow_2",
                "source": "channel_1",
                "target": "mom_1",
                "data": {},
            },
            {
                "id": "edge_thermal_mom_channel",
                "source": "mom_1",
                "target": "channel_1",
                "type": "thermal",
                "data": {"type": "thermal", "thickness": 0.003, "conductivity": 20.0},
            },
        ],
    )

    net = build_network_from_schema(schema)

    assert "edge_thermal_mom_channel" in net.walls
    mom = net.nodes["mom_1"]
    channel = net.elements["channel_1"]

    assert mom.surface.Nu_multiplier == 1.3
    # Momentum chamber is lossless by design; friction multiplier is fixed.
    assert mom.surface.f_multiplier == 1.0
    assert channel.surface.Nu_multiplier == 1.2
    assert channel.surface.f_multiplier == 0.9
