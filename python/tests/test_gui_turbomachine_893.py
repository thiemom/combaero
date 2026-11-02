"""#893: the GUI's compressor and turbine nodes -- built from their data,
with the efficiency type choosing eta_poly or eta_is."""

from __future__ import annotations

import pytest
import test_gui_http_routes as http

from combaero.network import CompressorElement, TurbineElement
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.schemas import NetworkGraphSchema


def _graph(kind: str, **data: object) -> dict:
    g = http._palette_network(kind)
    for n in g["nodes"]:
        if n["id"] == "x":
            n["data"] = {**n["data"], **data}
    return g


@pytest.mark.parametrize(
    ("kind", "cls"), [("compressor", CompressorElement), ("turbine", TurbineElement)]
)
@pytest.mark.parametrize("etype", ["polytropic", "isentropic"])
def test_the_node_builds_its_element(kind: str, cls: type, etype: str) -> None:
    net = build_network_from_schema(
        NetworkGraphSchema(**_graph(kind, efficiency_type=etype, efficiency=0.83))
    )
    el = net.elements["x"]
    assert isinstance(el, cls)
    assert (el.eta_poly, el.eta_is) == ((0.83, None) if etype == "polytropic" else (None, 0.83))


@pytest.mark.parametrize("kind", ["compressor", "turbine"])
def test_the_palette_node_solves_and_reports_its_power(kind: str) -> None:
    body = http._post("/solve", http._palette_network(kind)).json()
    assert body["success"], body["message"]
    res = body["element_results"]["x"]
    assert res["power"] > 0.0
    assert (res["dh"] > 0.0) if kind == "compressor" else (res["dh"] < 0.0)


@pytest.mark.parametrize(("field", "value"), [("pressure_ratio", 0.5), ("efficiency", 1.2)])
def test_invalid_data_is_refused_by_name(field: str, value: float) -> None:
    with pytest.raises(ValueError, match="'x'"):
        build_network_from_schema(NetworkGraphSchema(**_graph("compressor", **{field: value})))


def test_a_default_compressor_burner_turbine_chain_solves() -> None:
    """The palette defaults dropped in a row between two ambient boundaries:
    compressor, channel, combustor on a fuel controller (1500 K), channel,
    turbine. The flow is what the channels and burner lose between the two
    pressure ratios. The defaults leave them 5 % of Pt (PR 4 and 3.8):
    measured, 3 (25 %) and 3.6 (10 %) choke the hot channel and the cold
    start stalls."""
    amb = {"Pt": 101325.0, "Tt": 300.0, "composition": {"source": "humid_air", "mode": "mole"}}
    fuel = {
        **http._drop_defaults("fuel_boundary"),
        "combustor_id": "burner",
        "composition": {"source": "fuel", "mode": "mole"},
    }
    ch = http._drop_defaults("channel")
    g = {
        "nodes": [
            http._node("in", "pressure_boundary", amb),
            http._node("comp", "compressor", http._drop_defaults("compressor")),
            http._node("c1", "channel", dict(ch)),
            http._node("burner", "combustor", http._drop_defaults("combustor")),
            http._node("fuel", "fuel_boundary", fuel),
            http._node("c2", "channel", dict(ch)),
            http._node("turb", "turbine", http._drop_defaults("turbine")),
            http._node("out", "pressure_boundary", amb),
        ],
        "edges": [
            http._edge("e1", "in", "comp"),
            http._edge("e2", "comp", "c1"),
            http._edge("e3", "c1", "burner"),
            http._edge("e4", "fuel", "burner"),
            http._edge("e5", "burner", "c2"),
            http._edge("e6", "c2", "turb"),
            http._edge("e7", "turb", "out"),
        ],
        "solver_settings": {"timeout": 30.0},
    }
    assert g["nodes"][6]["data"]["pressure_ratio"] == 3.8
    body = http._post("/solve", g).json()
    assert body["success"], body["message"]
    assert body["node_results"]["burner"]["state"]["Tt"] == pytest.approx(1500.0, abs=1e-3)
    el = body["element_results"]
    assert el["turb"]["power"] > 0.0 and el["comp"]["power"] > 0.0
