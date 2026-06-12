"""
GUI dispatch tests for the mpce_tee element type.

Verifies that the graph_builder converts an `mpce_tee` node to an
`MPCEv2Element` with the correct inlet/outlet wiring per the
``flow_direction`` field, and that the schema accepts both values.
"""

from __future__ import annotations

import pytest

from combaero.network.mpce_v2_element import MPCEv2Element
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.schemas import MPCETeeData, NetworkGraphSchema


def _three_plenum_mpce_schema(flow_direction: str) -> NetworkGraphSchema:
    """Three plenums + an mpce_tee in the middle. Plenum<->tee edges use the
    standard port-{common,straight,branch}-{target,source} handle convention.
    For flow_direction="branch": common is upstream (plenum_com -> tee.common).
    For flow_direction="merge":  common is downstream (tee.common -> plenum_com).
    """
    if flow_direction == "branch":
        edges = [
            {
                "id": "e_com",
                "source": "plenum_com",
                "target": "mpce",
                "targetHandle": "port-common-target",
                "data": {},
            },
            {
                "id": "e_str",
                "source": "mpce",
                "target": "plenum_str",
                "sourceHandle": "port-straight-source",
                "data": {},
            },
            {
                "id": "e_bra",
                "source": "mpce",
                "target": "plenum_bra",
                "sourceHandle": "port-branch-source",
                "data": {},
            },
        ]
    else:
        edges = [
            {
                "id": "e_str",
                "source": "plenum_str",
                "target": "mpce",
                "targetHandle": "port-straight-target",
                "data": {},
            },
            {
                "id": "e_bra",
                "source": "plenum_bra",
                "target": "mpce",
                "targetHandle": "port-branch-target",
                "data": {},
            },
            {
                "id": "e_com",
                "source": "mpce",
                "target": "plenum_com",
                "sourceHandle": "port-common-source",
                "data": {},
            },
        ]

    return NetworkGraphSchema(
        nodes=[
            {"id": "plenum_com", "type": "plenum", "position": {"x": 0, "y": 0}, "data": {}},
            {"id": "plenum_str", "type": "plenum", "position": {"x": 400, "y": 0}, "data": {}},
            {"id": "plenum_bra", "type": "plenum", "position": {"x": 200, "y": 200}, "data": {}},
            {
                "id": "mpce",
                "type": "mpce_tee",
                "position": {"x": 200, "y": 0},
                "data": {"flow_direction": flow_direction, "theta_deg": 45.0},
            },
        ],
        edges=edges,
    )


def _get_element(net, elem_id: str):
    if elem_id not in net.elements:
        raise AssertionError(
            f"element {elem_id!r} not in net.elements (have: {list(net.elements)})"
        )
    return net.elements[elem_id]


def test_mpce_tee_branch_builds_separating_element():
    schema = _three_plenum_mpce_schema("branch")
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    assert isinstance(e, MPCEv2Element)
    assert e.flow_direction == "branch"
    assert len(e.inlet_nodes) == 1
    assert len(e.outlet_nodes) == 2
    assert e.port_angles_deg == [0.0, 0.0, 45.0]


def test_mpce_tee_merge_builds_joining_element():
    schema = _three_plenum_mpce_schema("merge")
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    assert isinstance(e, MPCEv2Element)
    assert e.flow_direction == "merge"
    assert len(e.inlet_nodes) == 2
    assert len(e.outlet_nodes) == 1
    assert e.port_angles_deg == [0.0, 45.0, 0.0]


def test_mpce_tee_data_defaults():
    d = MPCETeeData()
    assert d.flow_direction == "branch"
    assert d.theta_deg == 90.0
    assert d.psi == 1.0


def test_mpce_tee_data_rejects_invalid_flow_direction():
    with pytest.raises(Exception):  # noqa: B017 - pydantic ValidationError
        MPCETeeData(flow_direction="auto")


def test_mpce_tee_psi_from_explicit_areas():
    """When both F_C and F_branch are given, psi is computed from them."""
    schema = _three_plenum_mpce_schema("branch")
    # Override the mpce node's data to set explicit areas
    for n in schema.nodes:
        if n.id == "mpce":
            n.data["F_C"] = 0.02
            n.data["F_branch"] = 0.01
            n.data["psi"] = 999.0  # should be overridden by computed F_C/F_branch
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    # Port areas: [F_C, F_C (straight), A_branch] for branch flow_direction
    assert e.port_areas[0] == 0.02
    assert e.port_areas[1] == 0.02
    assert e.port_areas[2] == 0.01
