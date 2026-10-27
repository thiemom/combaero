"""The GUI Heat Source node (#792): an EnergyBoundary on the plenum, momentum
chamber or combustor its one edge points at, in any of its three modes; the
solved heat on the source's own result."""

from __future__ import annotations

import math

import pytest
from test_gui_http_routes import _assert_readable_error, _edge, _node, _post

from combaero.network import CombustorNode, MomentumChamberNode, PlenumNode
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.schemas import HeatSourceData, NetworkGraphSchema


def _graph(hs: dict | None = None, target_kind: str = "plenum", wire_to: str = "pl") -> dict:
    """in -> c1 -> pl -> c2 -> out at 300 K, a heat source x wired to `wire_to`."""
    pb = "pressure_boundary"
    nodes = [
        _node("in", pb, {"Pt": 1.05e5, "Tt": 300.0}),
        _node("c1", "channel", {"L": 0.5, "D": 0.05}),
        _node("pl", target_kind, {"area": 0.01} if target_kind == "momentum_chamber" else {}),
        _node("c2", "channel", {"L": 0.5, "D": 0.05}),
        _node("out", pb, {"Pt": 1.0e5, "Tt": 300.0}),
        _node("x", "heat_source", hs if hs is not None else {}),
    ]
    edges = [
        _edge("e1", "in", "c1"),
        _edge("e2", "c1", "pl"),
        _edge("e3", "pl", "c2"),
        _edge("e4", "c2", "out"),
    ]
    if wire_to:
        edges.append(_edge("h", "x", wire_to, sh="heat-out", th="flow-top"))
    return {"nodes": nodes, "edges": edges, "solver_settings": {"timeout": 20.0}}


def _build(graph: dict):
    return build_network_from_schema(NetworkGraphSchema(**graph))


def test_the_defaults_are_a_working_example() -> None:
    d = HeatSourceData()
    assert (d.mode, d.Q, d.fraction, d.Tt_target) == ("Q", 1000.0, 0.01, 600.0)


@pytest.mark.parametrize(
    ("kind", "cls"),
    [
        ("plenum", PlenumNode),
        ("momentum_chamber", MomentumChamberNode),
        ("combustor", CombustorNode),
    ],
)
@pytest.mark.parametrize("mode", ["Q", "fraction", "Tt_target"])
def test_each_mode_puts_one_energy_boundary_on_the_wired_node(
    kind: str, cls: type, mode: str
) -> None:
    net = _build(_graph({"mode": mode, "Q": 2500.0, "fraction": -0.02, "Tt_target": 450.0}, kind))
    node = net.nodes["pl"]
    assert isinstance(node, cls)
    (eb,) = node.energy_boundaries
    assert eb.id == "x"
    assert (eb.Q, eb.fraction, eb.Tt_target) == {
        "Q": (2500.0, 0.0, None),
        "fraction": (0.0, -0.02, None),
        "Tt_target": (0.0, 0.0, 450.0),
    }[mode]
    # the heat edge is not a flow path: no element or node was made of it
    assert "x" not in net.nodes and "x" not in net.elements
    assert not any("x" in eid for eid in net.elements)


def test_a_set_heat_and_its_outlet_temperature_round_trip_through_http() -> None:
    """Q mode heats the plenum; its Tt set as the target gives the Q back,
    reported on the source in both modes beside the plenum's TOTAL Tt."""
    r1 = _post("/solve", _graph({"mode": "Q", "Q": 1000.0}))
    assert r1.status_code == 200, r1.text[:300]
    b1 = r1.json()
    assert b1["success"], b1["message"]
    res1 = b1["node_results"]["x"]
    assert res1["Q"] == pytest.approx(1000.0)
    assert res1["heated_node"] == "pl"
    Tt1 = res1["state"]["Tt"]
    assert Tt1 == pytest.approx(b1["node_results"]["pl"]["state"]["Tt"])
    assert Tt1 > 300.0

    r2 = _post("/solve", _graph({"mode": "Tt_target", "Tt_target": Tt1}))
    b2 = r2.json()
    assert b2["success"], b2["message"]
    res2 = b2["node_results"]["x"]
    assert res2["Q"] == pytest.approx(1000.0, rel=1e-7)
    assert res2["Tt_target"] == pytest.approx(Tt1)
    assert b2["node_results"]["pl"]["state"]["Q"] == pytest.approx(1000.0, rel=1e-7)


def test_a_fraction_reports_the_heat_it_applied() -> None:
    b = _post("/solve", _graph({"mode": "fraction", "fraction": 0.01})).json()
    assert b["success"], b["message"]
    q = b["node_results"]["x"]["Q"]
    assert q > 0.0
    assert q == pytest.approx(b["node_results"]["pl"]["state"]["Q_fraction"], rel=1e-12)


def test_no_through_flow_reports_no_heat() -> None:
    """Equal boundary pressures: no flow, so no Q meets a target. Not a 500 and
    never a garbage number: the solve fails and the heat is null (NaN)."""
    g = _graph({"mode": "Tt_target", "Tt_target": 400.0})
    g["nodes"][0]["data"]["Pt"] = 1.0e5
    r = _post("/solve", g)
    assert r.status_code == 200, r.text[:300]
    b = r.json()
    assert not b["success"]
    assert b["node_results"]["x"]["Q"] is None or math.isnan(b["node_results"]["x"]["Q"])
    assert "no through-flow" in b["message"]


# --- refusals: a readable 4xx naming the source --------------------------------


def test_an_unwired_heat_source_is_refused() -> None:
    _assert_readable_error(_post("/solve", _graph(wire_to="")), "Heat source 'x'", "not wired")


@pytest.mark.parametrize("target", ["c1", "in"])
def test_a_heat_source_on_an_element_or_a_boundary_is_refused(target: str) -> None:
    _assert_readable_error(
        _post("/solve", _graph(wire_to=target)),
        "Heat source 'x'",
        "plenum, momentum chamber or combustor",
    )


def test_a_heat_source_on_two_nodes_is_refused() -> None:
    g = _graph()
    g["nodes"].append(_node("p2", "plenum"))
    g["edges"].append(_edge("h2", "x", "p2", sh="heat-out", th="flow-top"))
    _assert_readable_error(_post("/solve", g), "Heat source 'x'", "2 nodes")


def test_two_outlet_temperatures_on_one_node_are_refused() -> None:
    g = _graph({"mode": "Tt_target", "Tt_target": 400.0})
    g["nodes"].append(_node("y", "heat_source", {"mode": "Tt_target", "Tt_target": 500.0}))
    g["edges"].append(_edge("h2", "y", "pl", sh="heat-out", th="flow-left"))
    _assert_readable_error(_post("/solve", g), "Heat source 'y'", "second time")


def test_two_set_heats_on_one_node_add_up() -> None:
    g = _graph({"mode": "Q", "Q": 600.0})
    g["nodes"].append(_node("y", "heat_source", {"mode": "Q", "Q": 400.0}))
    g["edges"].append(_edge("h2", "y", "pl", sh="heat-out", th="flow-left"))
    b = _post("/solve", g).json()
    assert b["success"], b["message"]
    assert b["node_results"]["pl"]["state"]["Q"] == pytest.approx(1000.0)


def test_a_non_physical_target_is_refused() -> None:
    _assert_readable_error(
        _post("/solve", _graph({"mode": "Tt_target", "Tt_target": -10.0})),
        "Node 'x'",
        "Tt_target",
    )
