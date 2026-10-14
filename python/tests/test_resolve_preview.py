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


def test_tee_reports_the_builder_resolved_sources() -> None:
    """A GUI tee's areas are resolved in graph_builder and handed to the
    element as explicit port areas; the preview must report where they really
    came from, not "user"."""
    from test_gui_mpce_tee import _mpce_with_channels_schema, _three_plenum_mpce_schema

    rep = resolve_preview(build_network_from_schema(_mpce_with_channels_schema("branch")))
    assert rep["mpce"].fields["F_C"].source == "channel 'ch_com'"
    assert rep["mpce"].fields["F_branch"].source == "channel 'ch_bra'"

    rep = resolve_preview(
        build_network_from_schema(_mpce_with_channels_schema("branch", F_C_override=0.02))
    )
    assert rep["mpce"].fields["F_C"].source == "user"
    assert rep["mpce"].fields["F_C"].value == 0.02

    # No channels at all: F_C falls back to the default, F_branch from psi.
    rep = resolve_preview(build_network_from_schema(_three_plenum_mpce_schema("branch")))
    assert rep["mpce"].fields["F_C"].source == "default"
    assert rep["mpce"].fields["F_branch"].source == "derived from F_C / psi"


def test_gui_momentum_chamber_with_null_area_inherits_it() -> None:
    """A freshly dropped chamber carries area=None (the GUI no longer seeds
    0.1, which counted as user-set): it inherits the connected channel's flow
    area, and the convective surface follows the inherited area."""
    from gui.backend.graph_builder import build_network_from_schema as build

    def _n(nid: str, ntype: str, data: dict) -> dict:
        return {"id": nid, "type": ntype, "data": {"id": nid, **data}, "position": {"x": 0, "y": 0}}

    def _e(eid: str, src: str, tgt: str) -> dict:
        return {
            "id": eid,
            "source": src,
            "sourceHandle": "flow-source",
            "target": tgt,
            "targetHandle": "flow-target",
            "data": {},
        }

    graph = {
        "nodes": [
            _n("pb_in", "pressure_boundary", {"Pt": 200000, "Tt": 300}),
            _n("ch_u", "channel", {"L": 1.0, "D": 0.1}),
            _n("mc", "momentum_chamber", {"area": None}),
            _n("ch_d", "channel", {"L": 1.0, "D": 0.1}),
            _n("pb_out", "pressure_boundary", {"Pt": 100000, "Tt": 300}),
        ],
        "edges": [
            _e("e1", "pb_in", "ch_u"),
            _e("e2", "ch_u", "mc"),
            _e("e3", "mc", "ch_d"),
            _e("e4", "ch_d", "pb_out"),
        ],
    }
    net = build(NetworkGraphSchema(**graph))
    rep = resolve_preview(net)
    area = rep["mc"].fields["area"]
    assert rep["mc"].status == "resolved"
    assert area.source == "channel 'ch_u'"
    assert rep["mc"].fields["Dh"].source == "channel 'ch_u'"
    assert abs(area.value - 0.25 * 3.141592653589793 * 0.1**2) < 1e-12
    assert net.nodes["mc"].surface.area == area.value

    # A user Dh: the area is derived from it, and says so.
    graph["nodes"][2]["data"]["Dh"] = 0.2
    rep = resolve_preview(build(NetworkGraphSchema(**graph)))
    assert rep["mc"].fields["Dh"].source == "user"
    assert rep["mc"].fields["area"].source == "derived from Dh"


def test_orifice_reports_pipe_diameter_and_beta() -> None:
    """The bore never inherits; the upstream pipe's D does (even when that pipe
    inherits its own D), and the preview reports it with beta = d / D."""
    from combaero.network import OrificeElement

    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=2e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PlenumNode("p1"))
    net.add_node(PlenumNode("p2"))
    net.add_node(PressureBoundary("out", Pt=1e5, Tt=300.0, Y=Y_AIR))
    net.add_element(ChannelElement("ch_u", "in", "p1", length=1.0, diameter=0.1))
    net.add_element(ChannelElement("ch_inh", "p1", "p2", length=1.0))  # inherits 0.1
    net.add_element(OrificeElement("ori", "p2", "out", diameter=0.05))

    f = resolve_preview(net)["ori"].fields
    assert f["diameter"].value == 0.05
    assert f["diameter"].source == "user"
    assert f["pipe_diameter"].value == 0.1
    assert f["pipe_diameter"].source == "channel 'ch_inh'"
    assert abs(f["beta"].value - 0.5) < 1e-12

    # No bore given: the default bore, marked as such.
    net2 = FlowNetwork()
    net2.add_node(PressureBoundary("in", Pt=2e5, Tt=300.0, Y=Y_AIR))
    net2.add_node(PressureBoundary("out", Pt=1e5, Tt=300.0, Y=Y_AIR))
    net2.add_element(OrificeElement("ori", "in", "out"))
    f = resolve_preview(net2)["ori"].fields
    assert f["diameter"].source == "default"
    assert f["pipe_diameter"].value is None
    assert f["beta"].value == 0.0


def test_audit_sources_discrete_loss_area_and_channel_dh() -> None:
    """Closing-audit findings (#559): a discrete loss that takes its area from
    the upstream channel names it (it reported no source), and a channel's
    unset Dh is the circular Dh = D, not "user"."""
    from combaero.network import ConstantFractionLoss, PressureLossElement

    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PlenumNode("p"))
    net.add_node(PressureBoundary("out", Pt=1e5, Tt=300.0, Y=Y_AIR))
    net.add_element(ChannelElement("ch", "in", "p", length=1.0, diameter=0.1))
    net.add_element(PressureLossElement("dl", "p", "out", correlation=ConstantFractionLoss(0.01)))

    rep = resolve_preview(net)
    assert rep["dl"].fields["area"].source == "channel 'ch'"
    assert rep["ch"].fields["Dh"].source == "derived from D (circular)"

    net.add_element(ChannelElement("ch2", "p", "out", length=1.0, diameter=0.1, Dh=0.05))
    assert resolve_preview(net)["ch2"].fields["Dh"].source == "user"
