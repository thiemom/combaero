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
    assert e.strict is False, (
        "GUI junctions must use soft mode: strict raises on transient "
        "wrong-sign Newton iterates and kills solves that would converge; "
        "the post-solve direction guard rejects artifact roots instead"
    )
    assert len(e.inlet_nodes) == 1
    assert len(e.outlet_nodes) == 2
    assert e.port_angles_deg == [0.0, 0.0, 45.0]


def test_mpce_tee_merge_builds_joining_element():
    schema = _three_plenum_mpce_schema("merge")
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    assert isinstance(e, MPCEv2Element)
    assert e.flow_direction == "merge"
    assert e.strict is False
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


# ---------------------------------------------------------------------------
# Channel-D inheritance (the fix that lets users build geometry by connecting
# channels rather than hand-editing F_C/F_branch in the Inspector).
# ---------------------------------------------------------------------------


def _mpce_with_channels_schema(
    flow_direction: str,
    D_common: float = 0.1,
    D_straight: float = 0.1,
    D_branch: float = 0.05,
    F_C_override: float | None = None,
    F_branch_override: float | None = None,
) -> NetworkGraphSchema:
    """Schema with channels attached to each tee port, optionally with
    explicit area overrides on the tee. Used to test inheritance precedence."""
    nodes = [
        # Outer plenums (terminate channel ends)
        {"id": "p_com", "type": "plenum", "position": {"x": -200, "y": 0}, "data": {}},
        {"id": "p_str", "type": "plenum", "position": {"x": 600, "y": 0}, "data": {}},
        {"id": "p_bra", "type": "plenum", "position": {"x": 200, "y": 300}, "data": {}},
        # Channels attached to each tee port
        {
            "id": "ch_com",
            "type": "channel",
            "position": {"x": -50, "y": 0},
            "data": {"L": 1.0, "D": D_common},
        },
        {
            "id": "ch_str",
            "type": "channel",
            "position": {"x": 450, "y": 0},
            "data": {"L": 1.0, "D": D_straight},
        },
        {
            "id": "ch_bra",
            "type": "channel",
            "position": {"x": 200, "y": 150},
            "data": {"L": 1.0, "D": D_branch},
        },
        # The tee under test
        {
            "id": "mpce",
            "type": "mpce_tee",
            "position": {"x": 200, "y": 0},
            "data": {"flow_direction": flow_direction, "theta_deg": 90.0},
        },
    ]
    if F_C_override is not None:
        nodes[-1]["data"]["F_C"] = F_C_override
    if F_branch_override is not None:
        nodes[-1]["data"]["F_branch"] = F_branch_override

    # Edges: plenum <-> channel <-> tee (channel directly to tee port handle).
    # For branch flow_direction: com is upstream, str+bra are downstream.
    if flow_direction == "branch":
        edges = [
            {"id": "e1", "source": "p_com", "target": "ch_com", "data": {}},
            {
                "id": "e2",
                "source": "ch_com",
                "target": "mpce",
                "targetHandle": "port-common-target",
                "data": {},
            },
            {
                "id": "e3",
                "source": "mpce",
                "target": "ch_str",
                "sourceHandle": "port-straight-source",
                "data": {},
            },
            {"id": "e4", "source": "ch_str", "target": "p_str", "data": {}},
            {
                "id": "e5",
                "source": "mpce",
                "target": "ch_bra",
                "sourceHandle": "port-branch-source",
                "data": {},
            },
            {"id": "e6", "source": "ch_bra", "target": "p_bra", "data": {}},
        ]
    else:  # "merge"
        edges = [
            {"id": "e1", "source": "p_str", "target": "ch_str", "data": {}},
            {
                "id": "e2",
                "source": "ch_str",
                "target": "mpce",
                "targetHandle": "port-straight-target",
                "data": {},
            },
            {"id": "e3", "source": "p_bra", "target": "ch_bra", "data": {}},
            {
                "id": "e4",
                "source": "ch_bra",
                "target": "mpce",
                "targetHandle": "port-branch-target",
                "data": {},
            },
            {
                "id": "e5",
                "source": "mpce",
                "target": "ch_com",
                "sourceHandle": "port-common-source",
                "data": {},
            },
            {"id": "e6", "source": "ch_com", "target": "p_com", "data": {}},
        ]

    return NetworkGraphSchema(nodes=nodes, edges=edges)


def test_inheritance_F_C_from_common_channel():
    """When F_C is unset, derive from common-arm channel diameter."""
    import math

    schema = _mpce_with_channels_schema("branch", D_common=0.1, D_branch=0.05)
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    # Common port area should match pi/4 * 0.1^2
    assert e.port_areas[0] == pytest.approx(math.pi / 4.0 * 0.1**2, rel=1e-9)


def test_inheritance_F_branch_from_branch_channel():
    """When F_branch is unset, derive from branch-arm channel diameter."""
    import math

    schema = _mpce_with_channels_schema("branch", D_common=0.1, D_branch=0.05)
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    # Branch port area should match pi/4 * 0.05^2
    assert e.port_areas[2] == pytest.approx(math.pi / 4.0 * 0.05**2, rel=1e-9)


def test_inheritance_psi_derived_from_channels():
    """psi = F_C / F_branch when both come from inherited channels.

    Reproduces the test-network case: D_common=0.1, D_branch=0.05 should
    give psi = (0.1/0.05)^2 = 4, NOT the schema-default psi=1.
    """
    schema = _mpce_with_channels_schema("branch", D_common=0.1, D_branch=0.05)
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    # port_areas[0] (com) / port_areas[2] (branch) = psi
    psi_implied = e.port_areas[0] / e.port_areas[2]
    assert psi_implied == pytest.approx(4.0, rel=1e-9)


def test_explicit_F_C_overrides_inheritance():
    """Explicit user F_C wins over the inherited channel value."""
    schema = _mpce_with_channels_schema("branch", D_common=0.1, D_branch=0.05, F_C_override=0.05)
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    assert e.port_areas[0] == 0.05  # user value, not inherited 0.00785


def test_inheritance_works_for_merge_direction():
    """Merge flow direction: channels still attach to the same port handles."""
    import math

    schema = _mpce_with_channels_schema("merge", D_common=0.1, D_branch=0.05)
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    # For merge: port_areas = [F_C (str), A_branch, F_C (com)]
    assert e.port_areas[0] == pytest.approx(math.pi / 4.0 * 0.1**2, rel=1e-9)
    assert e.port_areas[1] == pytest.approx(math.pi / 4.0 * 0.05**2, rel=1e-9)
    assert e.port_areas[2] == pytest.approx(math.pi / 4.0 * 0.1**2, rel=1e-9)


def test_inheritance_falls_back_to_default_when_no_channels():
    """No channels attached -> fall back to F_C=0.01 / psi=1 default."""
    schema = _three_plenum_mpce_schema("branch")  # plenums directly on tee
    net = build_network_from_schema(schema)
    e = _get_element(net, "mpce")
    # No channels to inherit from; defaults apply
    assert e.port_areas[0] == 0.01
    assert e.port_areas[2] == 0.01  # F_branch derives from F_C/psi=1
