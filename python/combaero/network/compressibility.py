"""Incompressible validity report (#588).

Every converged solve checks each element running an incompressible closure
against one global Mach limit, ``INCOMPRESSIBLE_MACH_LIMIT`` (0.2 by default,
the solver's ``incompressible_mach_limit``). Above it the element is flagged:
a REPORT, never a change to the solution or its success.

Why one Mach number. The incompressible closures are the low-Mach limits of
their compressible counterparts (the gap closes as Ma^2, test_incompressible_
limit.py), so the error is a function of Mach. Measured at the same flow,
incompressible against compressible total-pressure drop:

* channel: 4.0% at Ma 0.2, 13% at 0.4 -- the worst element so far;
* orifice: 2.9% at a vena-contracta Ma of 0.2, 5% at 0.26.

Ma 0.2 keeps every element measured under ~5%. That is a rule of thumb, not a
guarantee: no generic criterion says when incompressible is "good enough" for
a gas, and the application stays the user's call. The tolerance can be set.

Which Mach. The element's own highest INTERNAL one, named by its
``validity_mach_keys``: an orifice's jet (vena contracta), a channel's
higher end. Port Mach can be low while the jet is not. A momentum chamber
running its incompressible closure is judged at its own Mach (#588).
"""

from __future__ import annotations

from typing import Any

#: Default limit [-]. See the module docstring for the measurement behind it.
INCOMPRESSIBLE_MACH_LIMIT = 0.2


def validity_verdicts(
    network: Any,
    element_diag: dict[str, dict[str, Any]],
    mach_limit: float,
    node_mach: dict[str, float | None] | None = None,
) -> dict[str, dict[str, Any]]:
    """Per incompressible element or chamber: its Mach, the limit, a verdict.

    ``"warn"`` above the limit, ``"ok"`` at or below it, ``"not_checked"``
    (with the reason) when no Mach is reported -- never a silent pass.
    Elements on a compressible closure, or with none, are left out.

    Nodes with a regime (momentum chambers, #588) are judged at their own
    Mach, ``node_mach[id]`` -- the velocity their closure used.
    """
    out: dict[str, dict[str, Any]] = {}
    for nid, node in network.nodes.items():
        if getattr(node, "regime", None) != "incompressible":
            continue
        m = (node_mach or {}).get(nid)
        if m is None:
            out[nid] = {
                "verdict": "not_checked",
                "reason": "no Mach reported (no solved state)",
                "limit": mach_limit,
            }
            continue
        mach = abs(float(m))
        out[nid] = {
            "verdict": "warn" if mach > mach_limit else "ok",
            "mach": mach,
            "limit": mach_limit,
        }
    for eid, element in network.elements.items():
        if getattr(element, "regime", None) != "incompressible":
            continue
        keys = getattr(element, "validity_mach_keys", ())
        diag = element_diag.get(eid) or {}
        machs = [abs(float(diag[k])) for k in keys if k in diag and diag[k] is not None]
        if not machs:
            reason = (
                "no Mach reported (no solved state)"
                if keys
                else "element does not report its internal Mach"
            )
            out[eid] = {"verdict": "not_checked", "reason": reason, "limit": mach_limit}
            continue
        mach = max(machs)
        out[eid] = {
            "verdict": "warn" if mach > mach_limit else "ok",
            "mach": mach,
            "limit": mach_limit,
        }
    return out
