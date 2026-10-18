"""Solver errors in the GUI name what the user drew (#321 item 1).

A tee port reached the user as ``'__tee_jct__node_1788936276012_branch'``,
which nothing on the canvas maps back to. ``gui.backend.messages.humanize``
rewrites it at the HTTP boundary to the node's label and the port letter the
node draws (c / s / b), and adds the hint for the usual cause.
"""

from __future__ import annotations

import asyncio

import pytest
from fastapi import HTTPException
from test_gui_mpce_tee import _mpce_with_channels_schema

from gui.backend.main import solve
from gui.backend.messages import humanize
from gui.backend.schemas import NetworkGraphSchema


def _schema(*nodes: tuple[str, str, dict]) -> NetworkGraphSchema:
    return NetworkGraphSchema(
        nodes=[
            {"id": i, "type": t, "position": {"x": 0.0, "y": 0.0}, "data": d} for i, t, d in nodes
        ],
        edges=[],
    )


def test_tee_port_id_becomes_the_node_label_and_port_letter() -> None:
    # Verbatim from #321.
    msg = (
        "MultiPortChamberElement 'node_1788936276012': port "
        "'__tee_jct__node_1788936276012_branch' is declared as an outlet (sign=+1) but the "
        "connecting element 'node_mttqgstk_0' is wired in the opposite direction (sign=-1)."
    )
    schema = _schema(
        ("node_1788936276012", "mpce_tee", {"label": "T1"}),
        ("node_mttqgstk_0", "channel", {}),
    )
    out = humanize(msg, schema)
    assert "__tee_jct__" not in out
    assert "port 'T1 (Tee) port branch (b)'" in out
    assert "'node_1788936276012' [T1 (Tee)]" in out
    assert "'node_mttqgstk_0' [Channel]" in out
    assert "flow direction (merge / branch)" in out


def test_an_inserted_tee_port_step_is_named_after_its_port() -> None:
    """The builder's sharp step at a port of another bore (#536) and its far
    node are internal ids too."""
    msg = (
        "Element '__tee_step__node_17_straight' failed: state at "
        "'__tee_step__node_17_straight__node' has negative pressure."
    )
    out = humanize(msg, _schema(("node_17", "mpce_tee", {"label": "T1"})))
    assert "__tee_step__" not in out and "__node" not in out
    assert out.count("the area step at T1 (Tee) port straight (s)") == 2


def test_ejector_boundary_port_gets_the_connect_through_an_element_hint() -> None:
    # Verbatim from #321: a boundary wired straight onto an ejector port.
    msg = (
        "MultiPortChamberElement 'node_1788936712695': port 'node_1788936870054' "
        "must be a MomentumChamberNode, got 'PressureBoundary'."
    )
    schema = _schema(
        ("node_1788936712695", "ejector", {}),
        ("node_1788936870054", "pressure_boundary", {"label": "exit"}),
    )
    out = humanize(msg, schema)
    assert "'node_1788936870054' [exit (Pressure Boundary)]" in out
    assert "through an element" in out


def test_a_message_without_internal_ids_is_left_alone() -> None:
    schema = _schema(("ch", "channel", {}))
    assert humanize("The solution converged.", schema) == "The solution converged."
    assert humanize("", schema) == ""


def test_solve_route_reports_a_miswired_tee_in_canvas_terms() -> None:
    """A tee declared as a split whose branch channel feeds INTO it: the route's
    400 names the tee's branch port the way the canvas does."""
    raw = _mpce_with_channels_schema("branch").model_dump()
    raw["edges"] = [e for e in raw["edges"] if e["id"] not in ("e5", "e6")] + [
        {"id": "e5", "source": "p_bra", "target": "ch_bra", "data": {}},
        {
            "id": "e6",
            "source": "ch_bra",
            "target": "mpce",
            "targetHandle": "port-branch-target",
            "data": {},
        },
    ]
    with pytest.raises(HTTPException) as exc:
        asyncio.run(solve(NetworkGraphSchema.model_validate(raw)))
    detail = str(exc.value.detail)
    assert "__tee_jct__" not in detail, detail
    assert "Tee port branch (b)" in detail, detail
    assert "flow direction (merge / branch)" in detail, detail


def test_the_declared_direction_hint_names_the_inspector_choices() -> None:
    """#640's hint names solver classes and attributes; the GUI user picks
    them in the Tee Inspector as a "Junction Model" and a "Flow Direction"."""
    msg = (
        "MultiPortChamberElement 'tee1' (flow_direction='merge') holds port(s) "
        "['__tee_jct__tee1_branch'] at zero or against its declared direction at "
        "the best iterate: a declared-direction junction does not reverse a port, "
        "so this operating point may have no root for it. Check the boundary "
        "pressures, change flow_direction, or use IdealMomentumTeeElement, which "
        "reverses ports."
    )
    out = humanize(msg, _schema(("tee1", "mpce_tee", {"label": "T1"})))
    assert '"Ideal momentum (parameter-free)" junction model' in out
    assert '"Momentum CV (Mynard)"' in out
    assert 'Flow Direction "Merge"' in out
    assert "change the tee's Flow Direction" in out
    assert "T1 (Tee) port branch (b)" in out
    for internal in (
        "IdealMomentumTeeElement",
        "MultiPortChamberElement",
        "flow_direction",
        "__tee_jct__",
    ):
        assert internal not in out
