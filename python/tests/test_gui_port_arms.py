"""A junction's ``port_<i>_*`` result keys belong to the arm the GUI says.

The Inspector's Live Telemetry and the Probe label a tee's or an ejector's
per-port keys with the arm letter the node card draws (C/S/B, P/S/O), from
``PORT_ARMS`` in ``gui/frontend/src/utils/portArms.ts``. That table encodes
graph_builder's port order (inlets first), which nothing checked: this solves
each junction through ``/solve`` and matches every port's flow and pressure to
the element or node wired to that arm.
"""

from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import pytest
from _asgi_client import request

from gui.backend.main import app

_REPO = Path(__file__).resolve().parents[2]
_ARMS_TS = _REPO / "gui" / "frontend" / "src" / "utils" / "portArms.ts"
sys.path.insert(0, str(_REPO / "scripts"))

import gui_result_keys as keys  # noqa: E402
import test_gui_http_routes as http  # noqa: E402


def _arms(name: str) -> list[str]:
    text = _ARMS_TS.read_text(encoding="utf-8")
    m = re.search(r"^\t" + name + r":\s*\[([^\]]*)\]", text, re.M)
    assert m, f"{name} not in {_ARMS_TS.name}"
    arms = re.findall(r'"(\w)"', m.group(1))
    assert len(arms) == 3, (name, arms)
    return arms


def _solve(graph: dict) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = request(app, "POST", "/solve", json_body=graph)
    assert r.status_code == 200, r.text[:300]
    body = r.json()
    assert body["success"], body["message"]
    return body


@pytest.mark.parametrize("model", ["mynard", "constant_k", "ideal_momentum"])
def test_tee_branch_ports_are_common_straight_branch(model: str) -> None:
    # c1 feeds the common port, c2 leaves the straight port, c3 the branch.
    body = _solve(keys._with(http._palette_network("mpce_tee"), "x", junction_model=model))
    x, el = body["element_results"]["x"], body["element_results"]
    wired = {"C": "c1", "S": "c2", "B": "c3"}
    flows = [abs(el[c]["m_dot"]) for c in ("c2", "c3")]
    assert abs(flows[0] - flows[1]) > 1e-3 * max(flows), "arms not distinguishable"
    for i, arm in enumerate(_arms("mpce_tee_branch")):
        assert abs(x[f"port_{i}_m_dot"]) == pytest.approx(abs(el[wired[arm]]["m_dot"]), rel=1e-9)


@pytest.mark.parametrize("model", ["mynard", "constant_k", "ideal_momentum"])
def test_tee_merge_ports_are_straight_branch_common(model: str) -> None:
    # c1 feeds the straight port, c2 the branch port, c3 leaves the common.
    graph = keys._merge_tee({"junction_model": model})
    # Unequal supplies, so the two inlet arms carry different flows.
    graph["nodes"][1]["data"] = {"Pt": 1.048e5, "Tt": 300.0}
    body = _solve(graph)
    x, el = body["element_results"]["x"], body["element_results"]
    wired = {"S": "c1", "B": "c2", "C": "c3"}
    flows = [abs(el[c]["m_dot"]) for c in ("c1", "c2")]
    assert abs(flows[0] - flows[1]) > 1e-3 * max(flows), "arms not distinguishable"
    for i, arm in enumerate(_arms("mpce_tee_merge")):
        assert abs(x[f"port_{i}_m_dot"]) == pytest.approx(abs(el[wired[arm]]["m_dot"]), rel=1e-9)


def test_ejector_ports_are_primary_secondary_outlet() -> None:
    # Each port is wired through a momentum chamber: mp, ms, mo.
    body = _solve(http._palette_network("ejector"))
    x, nodes = body["element_results"]["x"], body["node_results"]
    wired = {"P": "mp", "S": "ms", "O": "mo"}
    for i, arm in enumerate(_arms("ejector")):
        assert x[f"port_{i}_P"] == pytest.approx(nodes[wired[arm]]["state"]["P"], rel=1e-9)
    assert x["m_dot_primary"] == pytest.approx(-x["port_0_m_dot"], rel=1e-12)
