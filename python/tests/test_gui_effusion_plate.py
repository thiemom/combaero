"""GUI 'effusion_plate' node (#471): the default is a working example.

Discharged into a momentum chamber's side handle the plate takes its gas side
from the chamber; discharged into a plenum it takes the imposed heat flux.
"""

from __future__ import annotations

import pytest

from combaero.network import EffusionPlateElement
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import EffusionPlateData, NetworkGraphSchema


def _node(nid, ntype, data):
    return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}


def _edge(eid, src, tgt, tgt_handle=None):
    e = {"id": eid, "source": src, "target": tgt, "data": {}}
    if tgt_handle:
        e["targetHandle"] = tgt_handle
    return e


def _liner(plate: dict | None = None) -> dict:
    return {
        "nodes": [
            _node("gas_in", "pressure_boundary", {"Pt": 1.003e5, "Tt": 1500.0}),
            _node("exit", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1500.0}),
            _node("cool", "pressure_boundary", {"Pt": 1.03e5, "Tt": 600.0}),
            _node("duct", "channel", {"L": 0.2, "D": 0.08}),
            _node("tail", "channel", {"L": 0.2, "D": 0.08}),
            _node("ch", "momentum_chamber", {}),
            _node("eff", "effusion_plate", plate or {}),
        ],
        "edges": [
            _edge("e1", "gas_in", "duct"),
            _edge("e2", "duct", "ch", "flow-target"),
            _edge("e3", "ch", "tail"),
            _edge("e4", "tail", "exit"),
            _edge("e5", "cool", "eff"),
            _edge("e6", "eff", "ch", "side-target"),
        ],
    }


def test_the_node_builds_andrews_plate_c_by_default() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_liner()))
    e = net.elements["eff"]
    assert isinstance(e, EffusionPlateElement)
    d = EffusionPlateData()
    assert e.hole_diameter == d.hole_diameter == 3.27e-3
    assert e.pitch_x == e.pitch_y == 15.24e-3
    assert e.n_holes == round(0.152 * 0.152 / 15.24e-3**2)
    assert net.nodes["ch"].main_inlet == "duct"


def test_the_default_liner_solves_with_a_chamber_gas_side() -> None:
    result = NetworkRunner.from_dict(_liner()).solve(timeout=120.0)
    r = result._element_results["eff"].model_dump()
    assert r["success"]
    assert r["gas_side_chamber"] == 1.0
    assert 600.0 < r["T_wall_cold"] <= r["T_wall_hot"] < r["T_gas"]
    assert 0.0 < r["eta_overall"] < 1.0


def test_a_plenum_discharge_uses_the_imposed_flux() -> None:
    schema = {
        "nodes": [
            _node("cool", "pressure_boundary", {"Pt": 1.03e5, "Tt": 600.0}),
            _node("burner", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1500.0}),
            _node("eff", "effusion_plate", {"gas_heat_flux": 1.0e5}),
        ],
        "edges": [_edge("e1", "cool", "eff"), _edge("e2", "eff", "burner")],
    }
    result = NetworkRunner.from_dict(schema).solve(timeout=120.0)
    r = result._element_results["eff"].model_dump()
    assert r["gas_side_chamber"] == 0.0
    assert r["q_wall"] == pytest.approx(1.0e5)
    assert r["T_wall_hot"] > r["T_wall_cold"] > 600.0
