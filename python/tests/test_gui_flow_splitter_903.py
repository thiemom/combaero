"""#903: the GUI's flow splitter node -- built from its data and wiring,
reported as one node, exported from a script's network and back."""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import pytest
import test_gui_http_routes as http

from combaero.network import (
    FlowNetwork,
    FlowSplitter,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    SplitterOutletElement,
)
from combaero.network.components import OrificeElement
from gui.backend.export import network_to_gui
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.messages import humanize
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import NetworkGraphSchema

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_flow_splitter_903 as lib  # noqa: E402


def _graph(**data: object) -> dict:
    g = http._palette_network("flow_splitter")
    for n in g["nodes"]:
        if n["id"] == "x":
            n["data"] = {**n["data"], **data}
    return g


def _solve(graph: dict) -> object:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkRunner.from_dict(json.loads(json.dumps(graph))).solve(timeout=60.0)


def test_the_node_builds_one_outlet_element_per_port() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_graph(fractions=[0.3, None])))
    outs = [e for e in net.elements.values() if isinstance(e, SplitterOutletElement)]
    assert [e.id for e in outs] == ["x____tee_jct__x_out0", "x____tee_jct__x_out1"]
    assert [e.fraction for e in outs] == [pytest.approx(0.3), pytest.approx(0.7)]
    assert [e.reference for e in outs] == [False, True]  # the remainder
    assert {e.from_node for e in outs} == {"in"}  # the boundary wired to "in"
    # An element on a port meets it at an inserted plenum (a node at rest).
    assert all(isinstance(net.nodes[e.to_node], PlenumNode) for e in outs)


def test_the_result_is_one_node_with_rows_per_outlet() -> None:
    res = _solve(_graph(fractions=[0.3, None], areas=[1e-4, None]))
    assert res.success, res.message
    r = res._element_results["x"].model_dump()
    assert r["n_outlets"] == 2.0 and r["reference_outlet"] == 1.0
    assert r["rows_fraction"] == [pytest.approx(0.3), pytest.approx(0.7)]
    assert r["rows_fraction_achieved"] == [
        pytest.approx(0.3, abs=1e-12),
        pytest.approx(0.7, abs=1e-12),
    ]
    assert r["fraction_error_max"] < 1e-12
    assert r["m_dot"] == pytest.approx(sum(r["rows_m_dot"]), rel=1e-12)
    # Outlet 1 is metered (a drop, its effective area, Cd and zeta from its
    # area); outlet 2, the reference, is open.
    assert r["rows_valve_dPt"][0] > 0.0 and r["rows_valve_CdA"][0] > 0.0
    assert r["rows_valve_Cd"][0] == pytest.approx(r["rows_valve_CdA"][0] / 1e-4)
    assert r["rows_valve_CdA"][1] is None and r["rows_valve_Cd"][1] is None
    assert r["n_needs_pump"] == 0.0


def test_a_reference_chosen_away_from_the_remainder() -> None:
    """The explicit reference wins; a wrong one (the wide path) is refused,
    naming the outlet the GUI draws, not the inserted node."""
    res = _solve(_graph(fractions=[0.5, None], reference=0))
    assert not res.success
    text = humanize(res.message, NetworkGraphSchema(**_graph()))
    assert "outlet 2" in text and "__tee_jct__" not in text


@pytest.mark.parametrize(
    ("data", "match"),
    [
        ({"fractions": [0.5, 0.4]}, "choose the reference"),
        ({"fractions": [0.5, 0.4], "reference": 0}, "not renormalised"),
        ({"fractions": [None, None]}, "at most one remainder"),
        ({"fractions": [0.5]}, "at least 2"),
        ({"fractions": [0.5, None], "reference": 2}, "not an outlet"),
        ({"fractions": [0.5, None], "areas": [1e-4]}, "areas has 1 entries"),
    ],
)
def test_a_bad_splitter_is_refused_by_name(data: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        build_network_from_schema(NetworkGraphSchema(**_graph(**data)))


def test_an_unwired_outlet_is_refused() -> None:
    g = _graph(fractions=[0.3, 0.3, None])
    with pytest.raises(ValueError, match="port 'out2' is not connected"):
        build_network_from_schema(NetworkGraphSchema(**g))


def test_two_outlets_to_one_node_are_refused() -> None:
    g = http._palette_network("flow_splitter")
    g["edges"] = [e for e in g["edges"] if e["id"] not in ("e2", "e3", "e4", "e5")] + [
        http._edge("d1", "x", "o1", sh="port-out0-source"),
        http._edge("d2", "x", "o1", sh="port-out1-source"),
    ]
    with pytest.raises(ValueError, match="wired to the same node"):
        build_network_from_schema(NetworkGraphSchema(**g))


def test_a_momentum_chamber_at_the_inlet_feeds_an_inserted_plenum() -> None:
    """A chamber carries one bulk velocity and refuses a split: wired straight
    to the inlet, it feeds a plenum through a lossless link."""
    g = http._palette_network("flow_splitter")
    g["nodes"] += [http._node("c0", "channel", {"L": 0.5, "D": 0.05})]
    g["nodes"] += [http._node("mc", "momentum_chamber", {})]
    g["edges"] = [e for e in g["edges"] if e["id"] != "e1"] + [
        http._edge("f1", "in", "c0"),
        http._edge("f2", "c0", "mc"),
        http._edge("f3", "mc", "x", th="port-in-target"),
    ]
    net = build_network_from_schema(NetworkGraphSchema(**g))
    inlet = {e.from_node for e in net.elements.values() if isinstance(e, SplitterOutletElement)}
    assert inlet == {"__tee_jct__x_in"}
    assert isinstance(net.nodes["mc"], MomentumChamberNode)
    res = _solve(g)
    assert res.success, res.message


def test_a_scripted_splitter_exports_and_solves_the_same() -> None:
    """network_to_gui writes the splitter as one node over its outlets; the
    GUI solve gives the direct solve's flows, pressures and valve drops."""
    net, _ = lib._network("pressure")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        doc = network_to_gui(net)
        direct = NetworkSolver(net).solve()
    kinds = {n["id"]: n["type"] for n in doc["nodes"]}
    assert kinds["sp"] == "flow_splitter"
    assert not any(k.startswith("sp__") for k in kinds)
    sp = next(n for n in doc["nodes"] if n["id"] == "sp")["data"]
    assert sp["fractions"] == [0.2, 0.3, pytest.approx(0.5)]
    assert sp["reference"] == 0
    res = _solve(doc)
    assert res.success, res.message
    gui = res._raw
    for o in lib.PATHS:
        assert gui[f"sp__{o}.m_dot"] == pytest.approx(direct[f"sp__{o}.m_dot"], rel=1e-9)
        assert gui[f"{o}.Pt"] == pytest.approx(direct[f"{o}.Pt"], rel=1e-9)
    rows = res._element_results["sp"].model_dump()["rows_valve_dPt"]
    for i, o in enumerate(lib.PATHS):
        want = direct["__element_diag__"][f"sp__{o}"]["valve_dPt"]
        assert rows[i] == pytest.approx(want, rel=1e-6, abs=1e-3)


def test_a_splitter_into_boundaries_and_plenums_wired_directly() -> None:
    """Outlets on nodes (no element between): the outlet elements end there."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("res", Pt=2.0e5, Tt=400.0, Y=lib.AIR))
    net.add_node(PlenumNode("in"))
    net.add_node(PlenumNode("p1"))
    net.add_node(PressureBoundary("b1", Pt=1.0e5, Tt=300.0, Y=lib.AIR))
    net.add_node(PressureBoundary("b2", Pt=1.0e5, Tt=300.0, Y=lib.AIR))
    net.add_element(OrificeElement("feed", "res", "in", diameter=0.02, Cd=0.8, correlation="fixed"))
    net.add_element(
        OrificeElement("h1", "p1", "b1", diameter=0.006, Cd=0.7, correlation="fixed", in_pipe=False)
    )
    FlowSplitter("sp", "in", {"p1": 0.4, "b2": None}, reference="p1").add_to(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        doc = network_to_gui(net)
        direct = NetworkSolver(net).solve()
    assert direct["__success__"], direct["__message__"]
    res = _solve(doc)
    assert res.success, res.message
    for e in ("sp__p1", "sp__b2"):
        assert res._raw[f"{e}.m_dot"] == pytest.approx(direct[f"{e}.m_dot"], rel=1e-9)
