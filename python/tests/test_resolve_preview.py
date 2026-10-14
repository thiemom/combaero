"""`resolve_preview` reports geometry inheritance without solving (#559).

It runs the authoritative resolve pass and reports, per element/node, each
geometry field's resolved value and provenance, plus whether the item would
refuse (so the GUI can show the inherited placeholder the solver will use, and
flag "you must set this" before Solve rather than as a bare failure).
"""

from __future__ import annotations

import combaero as cb
from combaero.network import (
    AreaChangeElement,
    ChannelElement,
    FlowNetwork,
    PlenumNode,
    PressureBoundary,
)
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.resolve_preview import resolve_preview
from gui.backend.schemas import NetworkGraphSchema

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def test_resolved_channel_reports_value_and_source() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=2e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PlenumNode("p"))
    net.add_node(PressureBoundary("out", Pt=1e5, Tt=300.0, Y=Y_AIR))
    net.add_element(ChannelElement("ch_u", "in", "p", length=1.0, diameter=0.1))
    net.add_element(ChannelElement("ch_inh", "p", "out", length=1.0))

    rep = resolve_preview(net)
    assert rep["ch_u"].status == "resolved"
    assert rep["ch_u"].fields["diameter"].value == 0.1
    assert rep["ch_u"].fields["diameter"].source == "user"
    assert rep["ch_inh"].status == "resolved"
    assert rep["ch_inh"].fields["diameter"].value == 0.1
    assert "ch_u" in rep["ch_inh"].fields["diameter"].source


def test_orphan_channel_is_unresolved_with_a_reason() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("a", Pt=2e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PressureBoundary("b", Pt=1e5, Tt=300.0, Y=Y_AIR))
    net.add_element(ChannelElement("ch_orphan", "a", "b", length=1.0))  # no geometry anywhere

    rep = resolve_preview(net)
    assert rep["ch_orphan"].status == "unresolved"
    assert "diameter" in rep["ch_orphan"].reason
    # The field is still present (value None) so the GUI can mark it required.
    assert rep["ch_orphan"].fields["diameter"].value is None


def test_area_change_reports_each_face() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=2e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PlenumNode("p"))
    net.add_node(PlenumNode("p2"))
    net.add_node(PressureBoundary("out", Pt=1e5, Tt=300.0, Y=Y_AIR))
    net.add_element(ChannelElement("ch_u", "in", "p", length=1.0, diameter=0.1))
    net.add_element(AreaChangeElement("ac", "p", "p2"))
    net.add_element(ChannelElement("ch_d", "p2", "out", length=1.0, diameter=0.08))

    rep = resolve_preview(net)
    assert rep["ac"].status == "resolved"
    assert "ch_u" in rep["ac"].fields["F0"].source
    assert "ch_d" in rep["ac"].fields["F1"].source
    assert rep["ac"].fields["F0"].value is not None


def test_preview_over_a_gui_schema_keeps_only_drawn_items() -> None:
    """The endpoint path: build from a GUI graph and keep only canvas ids."""
    graph = {
        "nodes": [
            {
                "id": "pb_in",
                "type": "pressure_boundary",
                "data": {"id": "pb_in", "Pt": 200000, "Tt": 300},
                "position": {"x": 0, "y": 0},
            },
            {
                "id": "pb_out",
                "type": "pressure_boundary",
                "data": {"id": "pb_out", "Pt": 100000, "Tt": 300},
                "position": {"x": 200, "y": 0},
            },
            {
                "id": "chan",
                "type": "channel",
                "data": {"id": "chan", "L": 1.0, "D": 0.1},
                "position": {"x": 100, "y": 50},
            },
        ],
        "edges": [
            {
                "id": "e1",
                "source": "pb_in",
                "sourceHandle": "flow-source",
                "target": "chan",
                "targetHandle": "flow-target",
                "data": {},
            },
            {
                "id": "e2",
                "source": "chan",
                "sourceHandle": "flow-source",
                "target": "pb_out",
                "targetHandle": "flow-target",
                "data": {},
            },
        ],
    }
    schema = NetworkGraphSchema(**graph)
    net = build_network_from_schema(schema)
    rep = resolve_preview(net)

    gui_ids = {n.id for n in schema.nodes}
    drawn = {k: v for k, v in rep.items() if k in gui_ids}
    assert "chan" in drawn
    assert drawn["chan"].status == "resolved"
    assert drawn["chan"].fields["diameter"].source == "user"


def _channel_graph() -> dict:
    return {
        "nodes": [
            {
                "id": "pb_in",
                "type": "pressure_boundary",
                "data": {"id": "pb_in", "Pt": 200000, "Tt": 300},
                "position": {"x": 0, "y": 0},
            },
            {
                "id": "pb_out",
                "type": "pressure_boundary",
                "data": {"id": "pb_out", "Pt": 100000, "Tt": 300},
                "position": {"x": 200, "y": 0},
            },
            {
                "id": "chan",
                "type": "channel",
                "data": {"id": "chan", "L": 1.0, "D": 0.1},
                "position": {"x": 100, "y": 50},
            },
        ],
        "edges": [
            {
                "id": "e1",
                "source": "pb_in",
                "sourceHandle": "flow-source",
                "target": "chan",
                "targetHandle": "flow-target",
                "data": {},
            },
            {
                "id": "e2",
                "source": "chan",
                "sourceHandle": "flow-source",
                "target": "pb_out",
                "targetHandle": "flow-target",
                "data": {},
            },
        ],
    }


def test_resolve_preview_sync_formats_for_the_api() -> None:
    """The endpoint's sync path: serialized value + source, keyed by canvas id."""
    from gui.backend.main import _resolve_preview_sync

    out = _resolve_preview_sync(NetworkGraphSchema(**_channel_graph()))
    assert out["items"]["chan"]["status"] == "resolved"
    assert out["items"]["chan"]["fields"]["diameter"]["value"] == 0.1
    assert out["items"]["chan"]["fields"]["diameter"]["source"] == "user"
