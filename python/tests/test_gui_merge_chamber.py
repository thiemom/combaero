"""GUI momentum chamber with side streams (#471).

Inflows on the chamber's "side-target" handle are side streams; the single
inflow on its flow handle is the declared main inlet.
"""

from __future__ import annotations

import pytest

from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import NetworkGraphSchema


def _node(nid, ntype, data):
    return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}


def _edge(eid, src, tgt, tgt_handle=None):
    e = {"id": eid, "source": src, "target": tgt, "data": {}}
    if tgt_handle:
        e["targetHandle"] = tgt_handle
    return e


def _schema(main_handle: str | None = "flow-target") -> dict:
    return {
        "nodes": [
            _node("gas_in", "pressure_boundary", {"Pt": 1.003e5, "Tt": 1200.0}),
            _node("exit", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1200.0}),
            _node("cool", "pressure_boundary", {"Pt": 1.03e5, "Tt": 600.0}),
            _node("duct", "channel", {"L": 0.2, "D": 0.08}),
            _node("tail", "channel", {"L": 0.2, "D": 0.08}),
            _node("ch", "momentum_chamber", {}),
            _node("hole", "orifice", {"diameter": 0.004, "correlation": "fixed", "Cd": 0.7}),
        ],
        "edges": [
            _edge("e1", "gas_in", "duct"),
            _edge("e2", "duct", "ch", main_handle),
            _edge("e3", "ch", "tail"),
            _edge("e4", "tail", "exit"),
            _edge("e5", "cool", "hole"),
            _edge("e6", "hole", "ch", "side-target"),
        ],
    }


def test_the_flow_handle_inflow_is_the_main_inlet() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema()))
    assert net.nodes["ch"].main_inlet == "duct"


def test_an_edge_without_a_handle_counts_as_main() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema(main_handle=None)))
    assert net.nodes["ch"].main_inlet == "duct"


def test_it_solves_and_conserves_mass() -> None:
    result = NetworkRunner.from_dict(_schema()).solve(timeout=120.0)
    e = result._element_results
    assert e["tail"].success
    assert e["tail"].m_dot == pytest.approx(e["duct"].m_dot + e["hole"].m_dot, rel=1e-8)


def test_two_main_inlets_are_refused() -> None:
    s = _schema()
    s["edges"][-1].pop("targetHandle")  # the side stream lands on the flow handle
    s["edges"].append(_edge("e7", "cool", "ch", "side-target"))
    with pytest.raises(ValueError, match="exactly one main inlet"):
        build_network_from_schema(NetworkGraphSchema(**s))
