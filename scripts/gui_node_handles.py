#!/usr/bin/env python3
"""The GUI's node handles, per node type: one table for the frontend, the
backend and the tests (#901).

The canvas draws a flow edge only between handles its two nodes define, and
the backend solves an edge by its endpoints alone -- so an edge that names a
handle its node lacks solves, but is drawn disconnected. This script reads
the node components (``gui/frontend/src/components/flowTypes.ts`` and
``nodes/*.tsx``) and writes, per node type, its flow inlets and outlets (a
handle's ``type="target"`` / ``type="source"``), its thermal handles, and
whether any handle id is computed (``id={...}``, then not checkable).

An inlet or outlet named ``port-*`` is a semantic port (an ejector's
secondary inlet); ``flow-*`` and id-less handles are positional, any of them
the same connection. ``null`` stands for a handle without an id.

Usage (repo root)::

    uv run python scripts/gui_node_handles.py          # rewrite the JSON
    uv run python scripts/gui_node_handles.py --check  # exit 1 when stale
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
COMPONENTS = REPO / "gui" / "frontend" / "src" / "components"
OUT = REPO / "gui" / "frontend" / "src" / "utils" / "nodeHandles.generated.json"


def _node_types() -> dict[str, Path]:
    """Node type -> its component file, from flowTypes.ts."""
    src = (COMPONENTS / "flowTypes.ts").read_text()
    imports = dict(re.findall(r'import (\w+) from "\./nodes/(\w+)(?:\.tsx)?";', src))
    block = src[src.index("export const nodeTypes") :]
    block = block[: block.index("};")]
    return {
        kind: COMPONENTS / "nodes" / f"{imports[comp]}.tsx"
        for kind, comp in re.findall(r"^\s*(\w+): (\w+),", block, re.M)
    }


def handles_of(path: Path) -> dict[str, Any]:
    src = path.read_text()
    inlets: list[str | None] = []
    outlets: list[str | None] = []
    thermal: list[str] = []
    dynamic = False
    for block in re.findall(r"<Handle\b(.*?)/>", src, re.S):
        kind = re.search(r'\btype="(source|target)"', block)
        m = re.search(r'\bid="([^"]+)"', block)
        if m is None and re.search(r"\bid=\{", block):
            dynamic = True
            continue
        hid = m.group(1) if m else None
        if hid is not None and hid.startswith("thermal"):
            thermal.append(hid)
            continue
        if hid is not None and hid.startswith("probe-"):
            continue  # a derived link's anchor (utils/probeLink.ts), never a flow edge
        if kind is None:
            dynamic = True  # a handle whose role is computed
            continue
        (inlets if kind.group(1) == "target" else outlets).append(hid)
    return {"inlets": inlets, "outlets": outlets, "thermal": thermal, "dynamic": dynamic}


def table() -> dict[str, Any]:
    return {kind: handles_of(path) for kind, path in sorted(_node_types().items())}


def render() -> str:
    return json.dumps(table(), indent="\t", sort_keys=True, ensure_ascii=True) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="exit 1 when the JSON is stale")
    args = ap.parse_args()
    text = render()
    if args.check:
        if not OUT.exists() or OUT.read_text() != text:
            sys.exit(f"{OUT.relative_to(REPO)} is stale: run scripts/gui_node_handles.py")
        print(f"{OUT.relative_to(REPO)} is current")
        return
    OUT.write_text(text)
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
