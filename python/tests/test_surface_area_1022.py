"""A node's wetted surface area is never its flow area (#1022).

``surface.area`` is the A in a wall side's h A; ``area`` is the cross-section
the node's dynamic head is referred to. When a MomentumChamberNode or
CombustorNode inherited its flow area, resolve_topology wrote the
cross-section into ``surface.area`` -- overwriting a user-given wetted area
(356x too small for a 0.05 m pipe feeding a ConvectiveSurface(area=0.7)), and
turning an unset surface (0, no convective surface) into one, so every
inherited-area node silently gained a convective surface. The same write sat
in the tee and MPCE port sizing, the ejector's outlet sizing and the
PressureLossElement's area inheritance.

Decision (issue comment, 2026-11-06): resolve never writes ``surface.area``;
an unset surface stays 0.

Each resolve test below fails with the assignment restored: the given 0.7
becomes the inherited cross-section, and the unset 0 becomes it too.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    CombustorNode,
    ConvectiveSurface,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    MultiPortChamberElement,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    TeeJunctionElement,
    ThermalWall,
    WallLayer,
)
from combaero.network.components import SmoothModel
from combaero.network.ejector_element import EjectorElement
from combaero.network.pressure_loss import ConstantHeadLoss

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_GIVEN = 0.7  # wetted area [m^2]; far from any cross-section below
_D = 0.05
_PIPE_AREA = math.pi * (_D / 2.0) ** 2
_TEE_BRANCH = 0.0012
_MPCE_B = 0.0015


def _surface(given: bool) -> ConvectiveSurface | None:
    return ConvectiveSurface(area=_GIVEN) if given else None


def _ends(net: FlowNetwork) -> None:
    net.add_node(MassFlowBoundary("mfb", m_dot=0.1, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb", Pt=1.0e5, Tt=300.0))


def _chamber_from_channel(given: bool) -> tuple[FlowNetwork, object]:
    """Dh (and so area) inherited from an upstream channel."""
    net = FlowNetwork()
    _ends(net)
    node = MomentumChamberNode("n", surface=_surface(given))
    net.add_node(node)
    net.add_element(ChannelElement("c1", "mfb", "n", length=0.5, diameter=_D))
    net.add_element(ChannelElement("c2", "n", "pb", length=0.5, diameter=_D))
    return net, node


def _chamber_from_dh(given: bool) -> tuple[FlowNetwork, object]:
    """area derived from a user Dh."""
    net = FlowNetwork()
    _ends(net)
    node = MomentumChamberNode("n", Dh=0.2, surface=_surface(given))
    net.add_node(node)
    net.add_element(LosslessConnectionElement("l1", "mfb", "n"))
    net.add_element(LosslessConnectionElement("l2", "n", "pb"))
    return net, node


def _chamber_from_neighbour(given: bool) -> tuple[FlowNetwork, object]:
    """area inferred from a non-channel neighbour (a combustor, #262)."""
    net = FlowNetwork()
    _ends(net)
    net.add_node(CombustorNode("burner", area=0.15))
    node = MomentumChamberNode("n", surface=_surface(given))
    net.add_node(node)
    net.add_element(LosslessConnectionElement("l1", "mfb", "burner"))
    net.add_element(LosslessConnectionElement("l2", "burner", "n"))
    net.add_element(LosslessConnectionElement("l3", "n", "pb"))
    return net, node


def _chamber_default(given: bool) -> tuple[FlowNetwork, object]:
    """No area anywhere: the nominal DEFAULT_CHAMBER_AREA."""
    net = FlowNetwork()
    _ends(net)
    node = MomentumChamberNode("n", surface=_surface(given))
    net.add_node(node)
    net.add_element(LosslessConnectionElement("l1", "mfb", "n"))
    net.add_element(LosslessConnectionElement("l2", "n", "pb"))
    return net, node


def _merge_chamber(given: bool) -> tuple[FlowNetwork, object]:
    """A merge chamber inherits its main inlet's area."""
    net = FlowNetwork()
    _ends(net)
    net.add_node(MassFlowBoundary("side", m_dot=0.01, Tt=300.0, Y=_Y))
    node = MomentumChamberNode("n", surface=_surface(given), main_inlet="main")
    net.add_node(node)
    net.add_element(ChannelElement("main", "mfb", "n", length=0.5, diameter=_D))
    net.add_element(LosslessConnectionElement("jet", "side", "n"))
    net.add_element(ChannelElement("out", "n", "pb", length=0.5, diameter=_D))
    return net, node


def _combustor_from_channel(given: bool) -> tuple[FlowNetwork, object]:
    net = FlowNetwork()
    _ends(net)
    node = CombustorNode("n", surface=_surface(given))
    net.add_node(node)
    net.add_element(ChannelElement("c1", "mfb", "n", length=0.5, diameter=_D))
    net.add_element(ChannelElement("c2", "n", "pb", length=0.5, diameter=_D))
    return net, node


def _combustor_from_neighbour(given: bool) -> tuple[FlowNetwork, object]:
    net = FlowNetwork()
    _ends(net)
    net.add_node(MomentumChamberNode("up", area=0.15))
    node = CombustorNode("n", surface=_surface(given))
    net.add_node(node)
    net.add_element(LosslessConnectionElement("l1", "mfb", "up"))
    net.add_element(LosslessConnectionElement("l2", "up", "n"))
    net.add_element(LosslessConnectionElement("l3", "n", "pb"))
    return net, node


def _combustor_default(given: bool) -> tuple[FlowNetwork, object]:
    net = FlowNetwork()
    _ends(net)
    node = CombustorNode("n", surface=_surface(given))
    net.add_node(node)
    net.add_element(LosslessConnectionElement("l1", "mfb", "n"))
    net.add_element(LosslessConnectionElement("l2", "n", "pb"))
    return net, node


def _tee_port(port: str) -> Callable[[bool], tuple[FlowNetwork, object]]:
    """A tee sizes its auto-area straight/branch port chambers."""

    def build(given: bool) -> tuple[FlowNetwork, object]:
        net = FlowNetwork()
        _ends(net)
        net.add_node(PressureBoundary("pb2", Pt=1.0e5, Tt=300.0))
        # The port chambers resolve before the common chamber has its area, so
        # their own search finds nothing (the nominal default, Dh unset) and
        # the tee's sizing (step 4 of its resolve) is what gives them an area.
        nodes = {}
        for nid in ("straight", "branch"):
            nodes[nid] = MomentumChamberNode(nid, surface=_surface(given and nid == port))
            net.add_node(nodes[nid])
        net.add_node(MomentumChamberNode("common"))
        net.add_element(ChannelElement("feed", "mfb", "common", length=0.5, diameter=_D))
        net.add_element(
            TeeJunctionElement(
                "tee",
                common_node="common",
                straight_node="straight",
                branch_node="branch",
                theta=90.0,
                F_branch=_TEE_BRANCH,
            )
        )
        net.add_element(LosslessConnectionElement("out_s", "straight", "pb"))
        net.add_element(LosslessConnectionElement("out_b", "branch", "pb2"))
        return net, nodes[port]

    build.expected_area = _PIPE_AREA if port == "straight" else _TEE_BRANCH

    return build


def _mpce_port(given: bool) -> tuple[FlowNetwork, object]:
    """An MPCE sizes an auto-area outlet port chamber from port_areas."""
    net = FlowNetwork()
    _ends(net)
    net.add_node(PressureBoundary("pb2", Pt=1.0e5, Tt=300.0))
    # As for the tee: the outlet ports resolve before the inlet chamber has
    # its area, so the MPCE's port sizing (step 3) is what sizes them.
    node = MomentumChamberNode("b", surface=_surface(given))
    net.add_node(node)
    net.add_node(MomentumChamberNode("c"))
    net.add_node(MomentumChamberNode("a"))
    net.add_element(ChannelElement("feed", "mfb", "a", length=0.5, diameter=_D))
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["a"],
            outlet_nodes=["b", "c"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[_PIPE_AREA, _MPCE_B, 0.001],
        )
    )
    net.add_element(LosslessConnectionElement("out_b", "b", "pb"))
    net.add_element(LosslessConnectionElement("out_c", "c", "pb2"))
    return net, node


_mpce_port.expected_area = _MPCE_B


def _ejector_outlet(given: bool) -> tuple[FlowNetwork, object]:
    """An ejector sizes its auto-area outlet chamber to the mixing area A_3."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("pin", Pt=7.0e5, Tt=440.0))
    net.add_node(PressureBoundary("sin", Pt=2.3e4, Tt=280.0))
    net.add_node(PressureBoundary("pb", Pt=4.0e4, Tt=300.0))
    net.add_node(MomentumChamberNode("p", area=_PIPE_AREA))
    net.add_node(MomentumChamberNode("s", area=_PIPE_AREA))
    node = MomentumChamberNode("o", surface=_surface(given))
    net.add_node(node)
    net.add_element(LosslessConnectionElement("lp", "pin", "p"))
    net.add_element(LosslessConnectionElement("ls", "sin", "s"))
    net.add_element(
        EjectorElement(
            id="ej1",
            primary_node="p",
            secondary_node="s",
            outlet_node="o",
            throat_area=3.14e-5,
            nozzle_exit_area=1.0e-4,
            mixing_area=8.0e-4,
            regime="compressible",
        )
    )
    net.add_element(LosslessConnectionElement("lo", "o", "pb"))
    return net, node


def _loss_element(given: bool) -> tuple[FlowNetwork, object]:
    """A PressureLossElement inherits its area from an upstream channel."""
    net = FlowNetwork()
    _ends(net)
    net.add_node(PlenumNode("n"))
    net.add_element(ChannelElement("c1", "mfb", "n", length=0.5, diameter=_D))
    elem = PressureLossElement(
        "loss", "n", "pb", correlation=ConstantHeadLoss(zeta=1.0), surface=_surface(given)
    )
    net.add_element(elem)
    return net, elem


_CASES = {
    "chamber_channel_Dh": _chamber_from_channel,
    "chamber_user_Dh": _chamber_from_dh,
    "chamber_neighbour": _chamber_from_neighbour,
    "chamber_default": _chamber_default,
    "merge_chamber": _merge_chamber,
    "combustor_channel_Dh": _combustor_from_channel,
    "combustor_neighbour": _combustor_from_neighbour,
    "combustor_default": _combustor_default,
    "tee_straight_port": _tee_port("straight"),
    "tee_branch_port": _tee_port("branch"),
    "mpce_port": _mpce_port,
    "ejector_outlet": _ejector_outlet,
    "pressure_loss": _loss_element,
}


@pytest.mark.parametrize("case", sorted(_CASES))
def test_a_given_surface_area_survives_resolve(case: str) -> None:
    net, obj = _CASES[case](True)
    net.resolve_all_topology()
    # The flow area was inherited (the path under test ran) and differs from
    # the wetted area, so a write of one into the other would show.
    assert obj.area and obj.area > 0.0
    assert obj.area != pytest.approx(_GIVEN)
    expected = getattr(_CASES[case], "expected_area", None)
    if expected is not None:
        # The element's own port sizing ran (not the node's default).
        assert obj.area == pytest.approx(expected)
    assert obj.surface.area == _GIVEN
    assert obj.has_convective_surface


@pytest.mark.parametrize("case", sorted(_CASES))
def test_an_unset_surface_stays_unset(case: str) -> None:
    net, obj = _CASES[case](False)
    net.resolve_all_topology()
    assert obj.area and obj.area > 0.0
    assert obj.surface.area == 0.0
    assert not obj.has_convective_surface


@pytest.mark.parametrize("case", sorted(_CASES))
def test_resolve_twice_is_idempotent(case: str) -> None:
    net, obj = _CASES[case](True)
    net.resolve_all_topology()
    net.resolve_all_topology()
    assert obj.surface.area == _GIVEN


# ---------------------------------------------------------------------------
# A wall-coupled network: the heat follows the given wetted area
# ---------------------------------------------------------------------------


def _walled(kind: str, node_area: float | None, wetted: float) -> FlowNetwork:
    """Hot gas through a node walled to a cold channel. The channel's own wetted
    area (pi D L = 0.94 m^2) exceeds every node area used, so the wall's area
    min(A_node, A_channel) is the node's surface."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.3e5, Tt=900.0))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    surface = ConvectiveSurface(model=SmoothModel(), area=wetted)
    if kind == "combustor":
        node = CombustorNode("hn", method="complete", area=node_area, surface=surface)
    else:
        node = MomentumChamberNode("hn", area=node_area, surface=surface)
    g.add_node(node)
    for eid, a, b in (("hin", "hi", "hn"), ("hout", "hn", "ho")):
        g.add_element(ChannelElement(eid, a, b, length=0.3, diameter=_D, regime="compressible"))
    g.add_node(PressureBoundary("ci", Pt=1.3e5, Tt=300.0))
    g.add_node(PressureBoundary("co", Pt=1.0e5, Tt=300.0))
    g.add_element(
        ChannelElement(
            "cch",
            "ci",
            "co",
            length=3.0,
            diameter=0.1,
            regime="compressible",
            surface=ConvectiveSurface(model=SmoothModel()),
        )
    )
    g.add_wall(ThermalWall("w", "hn", "cch", layers=[WallLayer(0.001, 20.0)]))
    return g


def _wall_q(net: FlowNetwork) -> float:
    s = NetworkSolver(net)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    return float(r["w.Q"])


@pytest.mark.parametrize("kind", ["chamber", "combustor"])
def test_wall_heat_follows_the_given_wetted_area(kind: str) -> None:
    # Inherited flow area (the pipe's) vs the same area given explicitly: the
    # flow area's provenance must not change the wall's heat. With the old
    # write the inherited node walled through its 0.002 m^2 cross-section
    # instead of the 0.1 m^2 wetted area -- ~50x less heat.
    q_inherited = _wall_q(_walled(kind, None, 0.1))
    q_explicit = _wall_q(_walled(kind, _PIPE_AREA, 0.1))
    assert q_inherited == pytest.approx(q_explicit, rel=1e-6)
    # And the heat scales with the wetted area: twice the area, more heat
    # (less than twice: the gas cools and the wall warms).
    q_double = _wall_q(_walled(kind, None, 0.2))
    assert abs(q_double) > 1.3 * abs(q_inherited)


def test_a_wall_on_an_unsurfaced_inherited_area_node_is_refused() -> None:
    # Previously the inherited cross-section made a surface appear, so this
    # wall ran on 0.002 m^2 of invented area. Now the node has no convective
    # surface and validation says so by name.
    net = _walled("chamber", None, 0.1)
    net.nodes["hn"].surface = ConvectiveSurface(model=SmoothModel())
    with pytest.raises(ValueError, match="has no convective surface"):
        NetworkSolver(net).solve()


# ---------------------------------------------------------------------------
# The GUI: no wetted-area field on these components; a wall drawn to one
# supplies it
# ---------------------------------------------------------------------------


def _gui_schema(kind: str, data: dict, walled: bool, wall: dict | None = None) -> object:
    from gui.backend.schemas import NetworkGraphSchema

    pos = {"x": 0.0, "y": 0.0}
    nodes = [
        {"id": "pin", "type": "pressure_boundary", "position": pos, "data": {"Pt": 1.2e5}},
        {"id": "pout", "type": "pressure_boundary", "position": pos, "data": {"Pt": 1.0e5}},
        {"id": "ch_in", "type": "channel", "position": pos, "data": {"L": 1.0, "D": _D}},
        {"id": "ch_out", "type": "channel", "position": pos, "data": {"L": 1.0, "D": _D}},
        {"id": "x", "type": kind, "position": pos, "data": data},
    ]
    if kind == "discrete_loss":
        # An element: pipe -> plenum -> loss -> outlet pipe.
        nodes.append({"id": "mid", "type": "plenum", "position": pos, "data": {}})
        chain = ["pin", "ch_in", "mid", "x", "ch_out", "pout"]
    else:
        chain = ["pin", "ch_in", "x", "ch_out", "pout"]
    edges = [
        {"id": f"e{i}", "source": a, "target": b, "data": {}}
        for i, (a, b) in enumerate(zip(chain, chain[1:], strict=False))
    ]
    if walled:
        wd = {"type": "thermal", "thickness": 0.003, "conductivity": 20.0, **(wall or {})}
        edges.append({"id": "w", "source": "x", "target": "ch_in", "type": "thermal", "data": wd})
    return NetworkGraphSchema(nodes=nodes, edges=edges)


def _gui_obj(net: FlowNetwork) -> object:
    return net.nodes.get("x") or net.elements["x"]


_GUI_KINDS = {
    "combustor": {"method": "complete"},
    "momentum_chamber": {},
    "discrete_loss": {"correlation_type": "constant_head"},
}


@pytest.mark.parametrize("kind", sorted(_GUI_KINDS))
def test_gui_wall_supplies_the_wetted_area_of_an_inherited_area_component(kind: str) -> None:
    from gui.backend.graph_builder import build_network_from_schema

    net = build_network_from_schema(_gui_schema(kind, _GUI_KINDS[kind], True, {"area": 0.07}))
    net.resolve_all_topology()
    obj = _gui_obj(net)
    assert obj.area == pytest.approx(_PIPE_AREA)  # inherited from the pipe
    assert obj.surface.area == pytest.approx(0.07)  # the wall's area, not the pipe's
    net.validate()  # walled: the side has a convective surface


@pytest.mark.parametrize("kind", sorted(_GUI_KINDS))
def test_gui_wall_per_side_area_is_that_sides_wetted_area(kind: str) -> None:
    from gui.backend.graph_builder import build_network_from_schema

    wall = {"area": None, "area_a": 0.03, "area_b": 0.05}
    net = build_network_from_schema(_gui_schema(kind, _GUI_KINDS[kind], True, wall))
    net.resolve_all_topology()
    assert _gui_obj(net).surface.area == pytest.approx(0.03)  # side a is the component


@pytest.mark.parametrize("kind", sorted(_GUI_KINDS))
def test_gui_unwalled_inherited_area_component_has_no_surface(kind: str) -> None:
    from gui.backend.graph_builder import build_network_from_schema

    net = build_network_from_schema(_gui_schema(kind, _GUI_KINDS[kind], False))
    net.resolve_all_topology()
    obj = _gui_obj(net)
    assert obj.area == pytest.approx(_PIPE_AREA)
    assert obj.surface.area == 0.0
