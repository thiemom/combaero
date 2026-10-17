"""GUI thermal wall orientation (#627, #629).

The wall edge declares how its two elements lie along it; with the flow signs
the solve decides parallel or counter flow, and reports it on the edge. A
saved file without the field keeps "aligned".
"""

from __future__ import annotations

import pytest

from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import NetworkGraphSchema


def _node(nid: str, ntype: str, data: dict) -> dict:
    return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}


def _edge(eid: str, src: str, tgt: str, data: dict | None = None) -> dict:
    return {"id": eid, "source": src, "target": tgt, "data": data or {}}


def _schema(orientation: str | None) -> dict:
    wall = {"type": "thermal", "thickness": 0.0005, "conductivity": 200.0, "area": 0.2}
    if orientation is not None:
        wall["orientation"] = orientation
    return {
        "nodes": [
            _node("hi", "pressure_boundary", {"Pt": 1.02e5, "Tt": 1600.0}),
            _node("ho", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1600.0}),
            _node("ci", "pressure_boundary", {"Pt": 1.02e5, "Tt": 400.0}),
            _node("co", "pressure_boundary", {"Pt": 1.0e5, "Tt": 400.0}),
            _node("hch", "channel", {"L": 1.0, "D": 0.02}),
            _node("cch", "channel", {"L": 1.0, "D": 0.02}),
        ],
        "edges": [
            _edge("h1", "hi", "hch"),
            _edge("h2", "hch", "ho"),
            _edge("c1", "ci", "cch"),
            _edge("c2", "cch", "co"),
            _edge("w", "hch", "cch", wall),
        ],
    }


def test_a_saved_file_without_the_field_builds_an_aligned_wall() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema(None)))
    assert net.walls["w"].orientation == "aligned"


def test_an_unknown_orientation_is_refused_when_the_wall_is_built() -> None:
    with pytest.raises(ValueError, match="orientation"):
        build_network_from_schema(NetworkGraphSchema(**_schema("counter")))


def test_the_declared_orientation_reaches_the_solve_and_comes_back_on_the_edge() -> None:
    """Forward flow on both sides: aligned is parallel, opposed is counter,
    and counterflow moves more heat between the same streams."""
    out = {}
    for orientation in ("aligned", "opposed"):
        result = NetworkRunner.from_dict(_schema(orientation)).solve(timeout=120.0)
        assert result.success, result.message
        out[orientation] = result._edge_results["w"]
    assert out["aligned"]["counterflow"] == 0.0
    assert out["opposed"]["counterflow"] == 1.0
    assert out["opposed"]["Q"] > out["aligned"]["Q"] > 0.0
    for r in out.values():
        assert 0.0 < r["effectiveness"] < 1.0
        assert r["NTU"] > 0.0
