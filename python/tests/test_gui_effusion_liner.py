"""GUI 'effusion_liner' node (#471): three named ports, N stations behind it.

coolant in (left) -> coolant out (right) along the backside duct; discharge
(top) into the gas, typically a chamber's "s" handle. The default is a
working example: Andrews' plate C holes and pitch over a 30 mm duct. Its
crossflow ratio follows the duct's own drop against the hole drive (not the
duct size); with 200 Pa along the duct and ~3.8 kPa across the holes it peaks
at 0.27, inside the Rohde-scored band.
"""

from __future__ import annotations

import pytest

from combaero.network import EffusionPlateElement, PlenumNode
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import EffusionLinerData, NetworkGraphSchema


def _node(nid, ntype, data):
    return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}


def _edge(eid, src, tgt, src_handle=None, tgt_handle=None):
    e = {"id": eid, "source": src, "target": tgt, "data": {}}
    if src_handle:
        e["sourceHandle"] = src_handle
    if tgt_handle:
        e["targetHandle"] = tgt_handle
    return e


def _schema(liner: dict | None = None, feed_channel: bool = False) -> dict:
    nodes = [
        _node("cin", "pressure_boundary", {"Pt": 1.04e5, "Tt": 600.0}),
        _node("cout", "pressure_boundary", {"Pt": 1.038e5, "Tt": 600.0}),
        _node("gas_in", "pressure_boundary", {"Pt": 1.003e5, "Tt": 1500.0}),
        _node("exit", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1500.0}),
        _node("duct", "channel", {"L": 0.2, "D": 0.08}),
        _node("tail", "channel", {"L": 0.2, "D": 0.08}),
        _node("ch", "momentum_chamber", {}),
        _node("ln", "effusion_liner", liner or {}),
    ]
    edges = [
        _edge("g1", "gas_in", "duct"),
        _edge("g2", "duct", "ch", tgt_handle="flow-target"),
        _edge("g3", "ch", "tail"),
        _edge("g4", "tail", "exit"),
        _edge("l2", "ln", "cout", src_handle="port-coolantout-source"),
        _edge("l3", "ln", "ch", src_handle="port-discharge-source", tgt_handle="side-target"),
    ]
    if feed_channel:
        nodes.append(_node("feed", "channel", {"L": 0.1, "D": 0.06}))
        edges += [
            _edge("f1", "cin", "feed"),
            _edge("l1", "feed", "ln", tgt_handle="port-coolantin-target"),
        ]
    else:
        edges.append(_edge("l1", "cin", "ln", tgt_handle="port-coolantin-target"))
    return {"nodes": nodes, "edges": edges}


def test_the_node_expands_into_stations_between_its_ports() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema()))
    n = EffusionLinerData().n_segments
    assert net.elements["ln__s0"].from_node == "cin"
    assert net.elements[f"ln__s{n}"].to_node == "cout"
    for i in range(1, n + 1):
        p = net.elements[f"ln__p{i}"]
        assert isinstance(p, EffusionPlateElement) and p.to_node == "ch"
    assert net.nodes["ch"].main_inlet == "duct"


def test_the_default_is_a_working_example_inside_the_rohde_band() -> None:
    """Measured 2026-10-08: U1/Vi peaks at 0.27 (station 1). At a 500 Pa duct
    drop it would be 0.41 -- and a taller duct barely helps (0.35 at twice the
    height), which is why the advice is about pressures, not size."""
    result = NetworkRunner.from_dict(_schema()).solve(timeout=120.0)
    r = result._element_results["ln"].model_dump()
    assert r["success"]
    assert 0.0 < r["bleed_fraction"] < 1.0
    assert r["U1_over_Vi_max"] < 0.35
    assert r["crossflow_cd_degraded"] == 0.0
    assert r["is_ingesting"] == 0.0
    assert len(r["rows_T_wall_hot"]) == EffusionLinerData().n_segments
    assert r["m_dot"] == pytest.approx(r["m_dot_out"] + r["m_dot_bleed"], rel=1e-8)


def test_a_channel_feeds_the_liner_through_a_plenum() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema(feed_channel=True)))
    entry = net.elements["ln__s0"].from_node
    assert isinstance(net.nodes[entry], PlenumNode)
    result = NetworkRunner.from_dict(_schema(feed_channel=True)).solve(timeout=120.0)
    assert result._element_results["ln"].success
