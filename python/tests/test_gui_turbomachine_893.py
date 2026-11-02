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
