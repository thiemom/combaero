"""#892: the GUI's fuel boundary -- built from its node data, refused by name
when it controls nothing or something that is not a combustor."""

from __future__ import annotations

import pytest
import test_gui_http_routes as http

from gui.backend.graph_builder import build_network_from_schema
from gui.backend.schemas import NetworkGraphSchema


def _graph(**fuel: object) -> dict:
    g = http._palette_network("fuel_boundary")
    for n in g["nodes"]:
        if n["id"] == "x":
            n["data"] = {k: v for k, v in {**n["data"], **fuel}.items() if v is not None}
    return g


def test_the_palette_fuel_boundary_meets_its_target() -> None:
    body = http._post("/solve", _graph()).json()
    assert body["success"], body["message"]
    state = body["node_results"]["x"]["state"]
    assert body["node_results"]["cb"]["state"]["Tt"] == pytest.approx(1500.0, abs=1e-3)
    assert state["phi"] < 1.0
    assert state["m_dot"] > 0.0


def test_the_rich_branch_is_taken_from_the_node() -> None:
    body = http._post("/solve", _graph(branch="rich")).json()
    assert body["success"], body["message"]
    assert body["node_results"]["x"]["state"]["phi"] > 1.0


def test_no_combustor_is_refused_by_name() -> None:
    with pytest.raises(ValueError, match="Fuel boundary 'x': no combustor to control"):
        build_network_from_schema(NetworkGraphSchema(**_graph(combustor_id=None)))


def test_a_combustor_id_naming_another_kind_is_refused() -> None:
    g = _graph(combustor_id="c1")
    body = http._post("/solve", g)
    assert "is not a CombustorNode" in body.text
