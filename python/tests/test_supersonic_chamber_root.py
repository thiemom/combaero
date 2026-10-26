"""A momentum chamber on the supersonic root of its closure is refused.

The chamber closure Pt = P0(P, T0, |m|/A) has two static roots for one
(Pt, m): subsonic, and supersonic below the sonic P*. No element in the
network accelerates a subsonic stream past Mach 1 (a constant-area channel
chokes at 1), so a chamber on the supersonic root is a root-finder artifact.

Reported from the GUI: combustor -> 1 m x 100 mm channel -> 95 mm orifice ->
20 bar, 33 kg/s. The first attempt converged with the channel-to-orifice
junction at Mach 1.46 (static 17.2 bar, Pt 56.7 bar), the channel "leaving"
at Mach 1.67 unchoked and the orifice's static pressure RISING across it --
reported as converged and consistent. The subsonic root has the same flow
and Pt with the junction at Mach 0.63.
"""

from __future__ import annotations

import warnings

import pytest

from combaero.network import NetworkSolver
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import NetworkGraphSchema

JUNCTION = "__junction__ch__ori"


def _node(nid: str, ntype: str, data: dict) -> dict:
    return {"id": nid, "type": ntype, "position": {"x": 0.0, "y": 0.0}, "data": data}


def _edge(src: str, sh: str, tgt: str, th: str) -> dict:
    return {
        "id": f"{src}-{tgt}",
        "source": src,
        "target": tgt,
        "sourceHandle": sh,
        "targetHandle": th,
        "data": {},
    }


def _schema() -> dict:
    air = {"source": "humid_air", "mode": "mole"}
    return {
        "nodes": [
            _node("air", "mass_boundary", {"m_dot": 32.0, "Tt": 897.0, "composition": air}),
            _node(
                "fuel",
                "mass_boundary",
                {"m_dot": 1.0, "Tt": 300.0, "composition": {"source": "fuel", "mode": "mole"}},
            ),
            _node("fuel_ori", "orifice", {"diameter": 0.5, "Cd": 0.8, "correlation": "fixed"}),
            _node("plenum", "plenum", {}),
            _node("comb", "combustor", {"method": "complete"}),
            _node("ch", "channel", {"L": 1.0, "D": 0.1}),
            _node("ori", "orifice", {"diameter": 0.095, "Cd": 0.6}),
            _node("out", "pressure_boundary", {"Pt": 2001325.0, "Tt": 300.0, "composition": air}),
        ],
        "edges": [
            _edge("fuel", "flow-source", "fuel_ori", "flow-target"),
            _edge("fuel_ori", "flow-source", "plenum", "flow-top"),
            _edge("air", "flow-source", "plenum", "flow-left"),
            _edge("plenum", "flow-right", "comb", "flow-target"),
            _edge("comb", "flow-source", "ch", "flow-target"),
            _edge("ch", "flow-source", "ori", "flow-target"),
            _edge("ori", "flow-source", "out", "flow-target"),
        ],
    }


def test_the_supersonic_chamber_root_is_refused() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = NetworkSolver(net).solve(auto_retry=False)
    assert sol["__converged__"] is True
    assert sol["__success__"] is False
    assert sol["__consistent__"] is False
    assert sol["__outcome__"] == "inconsistent"
    assert sol["__supersonic_chambers__"] == [JUNCTION]
    assert "supersonic root" in sol["__message__"]
    assert f"'{JUNCTION}' at Mach 1.46" in sol["__message__"]
    assert "fed by ['ch']" in sol["__message__"]


def test_the_retry_then_finds_the_subsonic_root() -> None:
    result = NetworkRunner.from_dict(_schema()).solve(timeout=60.0)
    assert result.success, result.message
    raw = result._raw
    assert raw["__supersonic_chambers__"] == []
    assert raw[f"{JUNCTION}.mach"] == pytest.approx(0.629, abs=0.01)
    diag = raw["__element_diag__"]
    assert diag["ch"]["mach_out"] < 1.0
    assert diag["ori"]["dP"] > 0.0
    # Same flow and total pressure as the refused root: only the static
    # branch differs.
    assert raw[f"{JUNCTION}.Pt"] == pytest.approx(5.665e6, rel=1e-3)
