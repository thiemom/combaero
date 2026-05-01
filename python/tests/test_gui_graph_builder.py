from gui.backend.graph_builder import _expand_initial_guess, build_network_from_schema
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
    # Momentum chamber is no longer hardcoded to 1.0; it should forward the UI value.
    assert mom.surface.f_multiplier == 0.8
    assert channel.surface.Nu_multiplier == 1.2
    assert channel.surface.f_multiplier == 0.9


def test_combustor_multipliers_forwarded():
    schema = NetworkGraphSchema(
        nodes=[
            {
                "id": "plenum_in",
                "type": "plenum",
                "position": {"x": 0.0, "y": 0.0},
                "data": {},
            },
            {
                "id": "comb_1",
                "type": "combustor",
                "position": {"x": 200.0, "y": 0.0},
                "data": {
                    "method": "complete",
                    "area": 0.1,
                    "Nu_multiplier": 1.4,
                    "f_multiplier": 0.7,
                },
            },
            {
                "id": "plenum_out",
                "type": "plenum",
                "position": {"x": 400.0, "y": 0.0},
                "data": {},
            },
        ],
        edges=[
            {"id": "e1", "source": "plenum_in", "target": "comb_1", "data": {}},
            {"id": "e2", "source": "comb_1", "target": "plenum_out", "data": {}},
        ],
    )

    net = build_network_from_schema(schema)
    comb = net.nodes["comb_1"]

    assert comb.surface.Nu_multiplier == 1.4
    assert comb.surface.f_multiplier == 0.7


# ---------------------------------------------------------------------------
# _expand_initial_guess: short-key translation for solver unknown names
# ---------------------------------------------------------------------------


def test_expand_initial_guess_empty_returns_empty():
    assert _expand_initial_guess({}, "node1") == {}


def test_expand_initial_guess_P_broadcasts_to_static_and_total():
    out = _expand_initial_guess({"P": 150000.0}, "plenum1")
    assert out == {"plenum1.P": 150000.0, "plenum1.Pt": 150000.0}


def test_expand_initial_guess_T_broadcasts_to_static_and_total():
    out = _expand_initial_guess({"T": 600.0}, "comb")
    assert out == {"comb.T": 600.0, "comb.Tt": 600.0}


def test_expand_initial_guess_m_dot_maps_without_broadcast():
    out = _expand_initial_guess({"m_dot": 1.5}, "loss_elem")
    assert out == {"loss_elem.m_dot": 1.5}


def test_expand_initial_guess_qualified_key_passthrough():
    # Already-qualified keys (hand-authored) must be preserved verbatim.
    out = _expand_initial_guess({"other_node.P": 200000.0}, "this_node")
    assert out == {"other_node.P": 200000.0}


def test_expand_initial_guess_explicit_Ptotal_wins_over_P_broadcast():
    """Regression: ensure explicit `Pt` overrides the `P -> Pt`
    broadcast regardless of dict iteration order."""
    # P first in dict
    out1 = _expand_initial_guess({"P": 100000.0, "Pt": 101325.0}, "n")
    assert out1["n.P"] == 100000.0
    assert out1["n.Pt"] == 101325.0

    # Pt first in dict
    out2 = _expand_initial_guess({"Pt": 101325.0, "P": 100000.0}, "n")
    assert out2["n.P"] == 100000.0
    assert out2["n.Pt"] == 101325.0


def test_expand_initial_guess_explicit_Ttotal_wins_over_T_broadcast():
    out = _expand_initial_guess({"T": 500.0, "Tt": 520.0}, "n")
    assert out["n.T"] == 500.0
    assert out["n.Tt"] == 520.0


def test_expand_initial_guess_mixed_short_and_qualified():
    out = _expand_initial_guess(
        {"P": 150000.0, "m_dot": 1.2, "other.X[0]": 0.767},
        "comb",
    )
    assert out == {
        "comb.P": 150000.0,
        "comb.Pt": 150000.0,
        "comb.m_dot": 1.2,
        "other.X[0]": 0.767,
    }
