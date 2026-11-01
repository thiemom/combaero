"""The GUI's handle table is what the node components define (#901).

``gui/frontend/src/utils/nodeHandles.generated.json`` (written by
``scripts/gui_node_handles.py``) lists each node type's flow inlets, outlets
and thermal handles. The frontend repairs a flow edge that names a handle its
node lacks from it on load (``utils/edgeHandles.ts``); this file keeps the
table current and the shipped examples clean.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "scripts"))

import gui_node_handles as handles  # noqa: E402

_EXAMPLES = sorted((_REPO / "gui" / "examples").glob("*.json"))


def test_the_table_is_current() -> None:
    """Rerun ``uv run python scripts/gui_node_handles.py`` after changing a
    node component's handles."""
    assert handles.OUT.read_text() == handles.render()


def test_every_node_type_is_listed_with_its_flow_handles() -> None:
    table = json.loads(handles.OUT.read_text())
    assert table["plenum"]["outlets"] == ["flow-right", "flow-bottom"]
    assert table["lossless_connection"]["inlets"] == [None]
    assert table["channel"]["inlets"] == ["flow-target"]
    assert "port-secondary-target" in table["ejector"]["inlets"]
    graph_kinds = {
        "plenum",
        "mass_boundary",
        "pressure_boundary",
        "channel",
        "orifice",
        "combustor",
    }
    assert graph_kinds <= set(table)


@pytest.mark.parametrize("path", _EXAMPLES, ids=lambda p: p.stem)
def test_shipped_example_edges_use_handles_their_nodes_have(path: Path) -> None:
    """An explicit handle id must exist on its node; no id is fine (the canvas
    uses the node's first handle on that side)."""
    table = json.loads(handles.OUT.read_text())
    data = json.loads(path.read_text())
    kind = {n["id"]: n["type"] for n in data["nodes"]}
    for e in data["edges"]:
        if (
            e.get("type") in ("thermal", "probe_link")
            or (e.get("data") or {}).get("type") == "thermal"
        ):
            continue
        for end, key, side in (
            ("source", "sourceHandle", "outlets"),
            ("target", "targetHandle", "inlets"),
        ):
            spec = table.get(kind[e[end]])
            if spec is None or spec["dynamic"] or e.get(key) is None:
                continue
            assert e[key] in spec[side], (path.stem, e["id"], key, e[key], spec[side])
