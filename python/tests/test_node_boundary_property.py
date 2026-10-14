"""Where the network ends is a property of the node, defined once (#545).

Graph validation, the solver and the GUI export each kept their own list of
boundary types: validation matched class names (PressureBoundary,
MassFlowBoundary, WallNode), the export matched GUI kinds and left the wall
out. ``NetworkNode.is_boundary`` (the network may end here) and
``is_flow_boundary`` (mass and enthalpy cross the control volume here) now
answer both questions for every reader, which is what a control-volume
balance over the network needs.
"""

from __future__ import annotations

import pytest

from combaero.network import (
    ChannelElement,
    CombustorNode,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    PlenumNode,
    PressureBoundary,
    WallNode,
)
from gui.backend.runner import NetworkRunner


@pytest.mark.parametrize(
    "node,boundary,flow",
    [
        (PressureBoundary("n"), True, True),
        (MassFlowBoundary("n"), True, True),
        (WallNode("n"), True, False),
        (PlenumNode("n"), False, False),
        (MomentumChamberNode("n"), False, False),
        (CombustorNode("n"), False, False),
    ],
    ids=lambda v: type(v).__name__ if not isinstance(v, bool) else str(v),
)
def test_each_node_says_whether_the_network_ends_there(node, boundary: bool, flow: bool) -> None:
    assert node.is_boundary is boundary
    assert node.is_flow_boundary is flow


class _Sink(PlenumNode):
    """A node type validation has never heard of, declared a boundary."""

    is_boundary = True


def _dead_end(node) -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=2e5, Tt=300.0))
    net.add_node(node)
    net.add_element(ChannelElement("ch", "in", node.id, length=1.0, diameter=0.05))
    return net


def test_validation_reads_the_property_not_the_class_name() -> None:
    """A dead end is refused at an interior node and accepted at any node that
    declares itself a boundary, whatever its class is called."""
    with pytest.raises(ValueError, match="interior node 'end'"):
        _dead_end(PlenumNode("end")).validate()
    _dead_end(_Sink("end")).validate()


def test_the_export_marks_every_boundary_including_a_wall() -> None:
    """The GUI export used its own kind list, which left the wall out."""

    def node(i: str, t: str, **data) -> dict:
        return {"id": i, "type": t, "position": {"x": 0.0, "y": 0.0}, "data": data}

    def edge(i: str, s: str, t: str) -> dict:
        return {
            "id": i,
            "source": s,
            "target": t,
            "sourceHandle": "flow-source",
            "targetHandle": "flow-target",
            "data": {},
        }

    schema = {
        "nodes": [
            node("in", "mass_boundary", m_dot=0.05, Tt=300.0),
            node("p", "plenum"),
            node("ch", "channel", L=1.0, D=0.05),
            node("stub", "channel", L=0.2, D=0.05),
            node("out", "pressure_boundary", Pt=101325.0, Tt=300.0),
            node("end", "wall"),
        ],
        "edges": [
            edge("e1", "in", "p"),
            edge("e2", "p", "ch"),
            edge("e3", "ch", "out"),
            edge("e4", "p", "stub"),
            edge("e5", "stub", "end"),
        ],
    }
    df = NetworkRunner.from_dict(schema).solve().to_dataframe()
    nodes = df[df["type"] == "node"].set_index("id")["is_boundary"]
    assert nodes["in"] and nodes["out"] and nodes["end"]
    assert not nodes["p"]
    assert not df[df["type"] != "node"]["is_boundary"].any()
