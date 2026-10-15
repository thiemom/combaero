"""The incompressible validity report (#588) in canvas terms.

The solver reports per network element; each flag is attributed to the canvas
node that owns it (a liner for its segments, a tee for its port step), as the
second-law flags are.
"""

from __future__ import annotations

from typing import Any

from .schemas import IncompressibleFlag
from .second_law_view import owner_of


def element_fields(verdict: dict[str, Any] | None) -> dict[str, Any]:
    """The verdict as flat fields on an element's result."""
    if not verdict:
        return {}
    out: dict[str, Any] = {"incompressible_validity": verdict["verdict"]}
    if "mach" in verdict:
        out["incompressible_mach"] = float(verdict["mach"])
    out["incompressible_mach_limit"] = float(verdict["limit"])
    return out


def warnings(result: dict[str, Any], net: Any, gui_ids: set[str]) -> list[IncompressibleFlag]:
    tee_steps = getattr(net, "_gui_tee_steps", {}) or {}
    verdicts = result.get("__incompressible_validity__") or {}
    return [
        IncompressibleFlag(
            element_id=eid,
            owner_id=owner_of(eid, gui_ids, tee_steps),
            mach=float(verdicts[eid]["mach"]),
            limit=float(verdicts[eid]["limit"]),
        )
        for eid in result.get("__incompressible_warnings__") or []
        if eid in verdicts
    ]
