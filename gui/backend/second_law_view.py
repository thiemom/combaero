"""The second-law verdicts and the ledger in canvas terms (#545).

The solver reports per NETWORK element, and some of those the user never
drew: a tee port's sharp step (``__tee_step__<tee>_<port>``), the sub-elements
a composite node builds (``<node>__x1``, ``<node>__s0``, ...). Each verdict is
attributed here to the node on the canvas that owns it, so the Inspector can
flag what the user can actually select.
"""

from __future__ import annotations

from typing import Any

from .schemas import SecondLawFlag

AUTO_LINK = "__auto_link__"


def _auto_link_ends(element_id: str) -> tuple[str, str] | None:
    """(source, target) of the lossless link the builder puts on a drawn
    node-to-node edge, ``__auto_link__<source>__<target>``."""
    if not element_id.startswith(AUTO_LINK):
        return None
    ends = element_id[len(AUTO_LINK) :].split("__")
    return (ends[0], ends[-1]) if len(ends) >= 2 else None


def owner_of(element_id: str, gui_ids: set[str], tee_steps: dict[str, tuple]) -> str:
    """The canvas node an element belongs to (itself when the user drew it).

    A drawn node-to-node edge has no element of its own on the canvas; its
    link belongs to the node it delivers to.
    """
    if element_id in gui_ids:
        return element_id
    if element_id in tee_steps:
        return str(tee_steps[element_id][0])
    ends = _auto_link_ends(element_id)
    if ends:
        return ends[1]
    head = element_id.split("__", 1)[0]
    return head if head in gui_ids else element_id


def element_fields(verdict: dict[str, Any] | None) -> dict[str, Any]:
    """The verdict as flat fields on an element's result."""
    if not verdict:
        return {}
    out: dict[str, Any] = {"second_law": verdict["verdict"]}
    for key in ("sigma_norm_low", "sigma_norm_high", "group"):
        if key in verdict:
            out[f"second_law_{key}" if key == "group" else key] = verdict[key]
    if "reason" in verdict:
        out["second_law_reason"] = verdict["reason"]
    return out


def failures(result: dict[str, Any], net: Any, gui_ids: set[str]) -> list[SecondLawFlag]:
    tee_steps = getattr(net, "_gui_tee_steps", {}) or {}
    verdicts = result.get("__second_law__") or {}
    return [
        SecondLawFlag(
            element_id=eid,
            owner_id=owner_of(eid, gui_ids, tee_steps),
            sigma_norm=float(verdicts[eid]["sigma_norm"]),
            sigma_norm_high=float(verdicts[eid]["sigma_norm_high"]),
        )
        for eid in result.get("__second_law_failures__") or []
        if eid in verdicts
    ]


def ledger(result: dict[str, Any], net: Any, gui_ids: set[str]) -> dict[str, Any]:
    """The ledger with each stream's element attributed to its canvas owner."""
    raw = result.get("__ledger__") or {}
    if not raw:
        return {}
    tee_steps = getattr(net, "_gui_tee_steps", {}) or {}
    streams = []
    for row in raw.get("streams", []):
        via = str(row["via"])
        ends = _auto_link_ends(via)
        if ends:
            # a direct edge: name the node at its far end
            far = ends[1] if row["node"] == ends[0] else ends[0]
            streams.append({**row, "via_owner": far, "direct": True})
        else:
            streams.append({**row, "via_owner": owner_of(via, gui_ids, tee_steps)})
    return {**raw, "streams": streams}
