"""The GUI combustor -> channel -> orifice net that once landed on a supersonic
chamber root now lands on its one consistent root (#861).

The chamber closure Pt = P0(P, T0, |m|/A) has two static roots for one
(Pt, m): subsonic, and supersonic below the sonic P*. A chamber on the
supersonic root is a root-finder artifact and is refused (exercised by the
junction runner's Wang cases, test_junction_network_runner.py).

Reported from the GUI: combustor -> 1 m x 100 mm channel -> 95 mm orifice ->
20 bar, 33 kg/s. The first attempt "converged" with the channel-to-orifice
junction at Mach 1.46 (Pt 56.7 bar), and the retry with it at Mach 0.63 at
the same Pt. Neither was a root. The orifice (beta 0.95) solves its own flow
through the approach pipe (#651), and its Newton 2-cycled across the pipe's
choke and stopped at a = 40 kg/s, where m_calc = 33. The consistent own flow
at Pt 56.7 bar is about 38 kg/s. Solved on the bracket, the orifice passes
33 kg/s at junction Pt 49.6 bar, Mach 0.86. The cold attempt reaches that
root directly, so no supersonic root is met on the way.
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


def _assert_the_consistent_root(raw: dict) -> None:
    assert raw["__supersonic_chambers__"] == []
    assert raw[f"{JUNCTION}.mach"] == pytest.approx(0.863, abs=0.005)
    assert raw[f"{JUNCTION}.Pt"] == pytest.approx(4.960e6, rel=1e-3)
    diag = raw["__element_diag__"]
    assert diag["ch"]["mach_out"] < 1.0
    assert diag["ori"]["dP"] > 0.0


def test_the_cold_solve_lands_on_the_consistent_root() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = NetworkSolver(net).solve(auto_retry=False)
    assert sol["__converged__"] is True
    assert sol["__success__"] is True
    _assert_the_consistent_root(sol)


def test_the_runner_finds_the_same_root() -> None:
    result = NetworkRunner.from_dict(_schema()).solve(timeout=60.0)
    assert result.success, result.message
    _assert_the_consistent_root(result._raw)
