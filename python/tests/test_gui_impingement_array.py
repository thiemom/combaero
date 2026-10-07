"""GUI 'impingement_array' node: one configuration, expanded into rows (#465).

The node is the channel-exit configuration only. The backend expands it into
ImpingementArray's row chain, reports the rows back under the node's own id,
and spreads a thermal edge to it over the rows.
"""

from __future__ import annotations

import pytest

from combaero.network import ImpingementCrossflowElement, ImpingementPlateElement, PlenumNode
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import ImpingementArrayData, NetworkGraphSchema, OrificeData


def _node(nid, ntype, data):
    return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}


def _edge(eid, src, tgt, data=None):
    return {"id": eid, "source": src, "target": tgt, "data": data or {}}


def _schema(with_wall: bool = False, upstream_channel: bool = False) -> dict:
    nodes = [
        _node("supply", "pressure_boundary", {"Pt": 1.03e5, "Tt": 300.0}),
        _node("exit", "pressure_boundary", {"Pt": 1.0e5, "Tt": 300.0}),
        _node("arr", "impingement_array", {}),  # every default: the working example
    ]
    edges = [_edge("e2", "arr", "exit")]
    if upstream_channel:
        nodes.append(_node("feed", "channel", {"L": 0.1, "D": 0.1}))
        edges += [_edge("e0", "supply", "feed"), _edge("e1", "feed", "arr")]
    else:
        edges.append(_edge("e1", "supply", "arr"))
    if with_wall:
        nodes += [
            _node("hot_in", "pressure_boundary", {"Pt": 1.02e5, "Tt": 1200.0}),
            _node("hot_out", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1200.0}),
            _node("hot", "channel", {"L": 0.127, "D": 0.05}),
        ]
        edges += [
            _edge("h1", "hot_in", "hot"),
            _edge("h2", "hot", "hot_out"),
            _edge(
                "wall",
                "hot",
                "arr",
                {"type": "thermal", "thickness": 0.001, "conductivity": 20.0},
            ),
        ]
    return {"nodes": nodes, "edges": edges}


def test_the_node_expands_into_the_row_chain() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema()))
    n = ImpingementArrayData().n_rows
    for i in range(1, n + 1):
        assert isinstance(net.elements[f"arr__p{i}"], ImpingementPlateElement)
        assert isinstance(net.elements[f"arr__x{i}"], ImpingementCrossflowElement)
        assert net.elements[f"arr__p{i}"].from_node == "supply"
    assert net.elements[f"arr__x{n}"].to_node == "exit"


def test_the_default_is_a_working_example() -> None:
    """Defaults solve and land inside Florschuetz's range on every row."""
    result = NetworkRunner.from_dict(_schema()).solve(timeout=120.0)
    r = result._element_results["arr"]
    assert r.success
    d = r.model_dump()
    assert d["Re_j_min"] > 2.5e3 and d["Re_j_max"] < 7.0e4
    assert d["Gc_Gj_max"] < 0.8
    assert d["surface_extrapolated"] == 0.0
    assert len(d["rows_Nu"]) == ImpingementArrayData().n_rows


def test_an_upstream_element_feeds_through_a_plenum() -> None:
    """A junction MomentumChamberNode cannot split one stream into N plates."""
    net = build_network_from_schema(NetworkGraphSchema(**_schema(upstream_channel=True)))
    supply = net.elements["arr__p1"].from_node
    assert isinstance(net.nodes[supply], PlenumNode)
    result = NetworkRunner.from_dict(_schema(upstream_channel=True)).solve(timeout=120.0)
    assert result._element_results["arr"].success


def test_one_thermal_edge_heats_every_row_and_reports_as_one() -> None:
    schema = _schema(with_wall=True)
    net = build_network_from_schema(NetworkGraphSchema(**schema))
    n = ImpingementArrayData().n_rows
    assert sorted(net.walls) == sorted(f"wall__r{i}" for i in range(1, n + 1))
    result = NetworkRunner.from_dict(schema).solve(timeout=120.0)
    edge = result._edge_results["wall"]
    assert len(edge["rows_Q"]) == n
    assert edge["Q"] == pytest.approx(sum(edge["rows_Q"]), rel=1e-12)
    assert edge["Q"] > 0.0
    assert edge["T_hot"] == max(edge["rows_T_hot"])


def test_orifice_offers_lichtarowicz() -> None:
    assert OrificeData(correlation="Lichtarowicz").correlation == "Lichtarowicz"
