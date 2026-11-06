"""#899: the cold start's reference state is per flow component. A part of
the network no element connects to the rest -- a coolant supply coupled to
the gas side only through a wall -- gets the start it would get alone."""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _gas(n: FlowNetwork) -> None:
    """Hot gas, pressure driven: 2.4 bar 808 K through two paths to 1.4 bar."""
    n.add_node(PressureBoundary("in", Pt=2.4e5, Tt=808.0, Y=AIR))
    n.add_node(PressureBoundary("out", Pt=1.4e5, Tt=808.0, Y=AIR))
    for k in range(2):
        n.add_node(MomentumChamberNode(f"m{k}", area=3e-4, regime="compressible"))
        n.add_element(
            ChannelElement(
                f"g{k}",
                "in",
                f"m{k}",
                length=0.05,
                diameter=0.02,
                regime="compressible",
                surface=ConvectiveSurface(model=SmoothModel()),
            )
        )
        n.add_element(
            ChannelElement(
                f"x{k}", f"m{k}", "out", length=1e-3, diameter=0.018, regime="compressible"
            )
        )


def _coolant(n: FlowNetwork) -> None:
    """A small, cold, mass-flow driven duct."""
    n.add_node(MassFlowBoundary("c_in", m_dot=0.0167, Tt=358.0, Y=AIR))
    n.add_node(PlenumNode("c_mid"))
    n.add_node(PressureBoundary("c_out", Pt=2.0e5, Tt=358.0, Y=AIR))
    n.add_element(
        ChannelElement(
            "c",
            "c_in",
            "c_mid",
            length=0.0762,
            diameter=0.0063,
            regime="compressible",
            surface=ConvectiveSurface(model=SmoothModel()),
        )
    )
    n.add_element(ChannelElement("c2", "c_mid", "c_out", length=0.01, diameter=0.0063))


def _x0(n: FlowNetwork) -> dict[str, float]:
    s = NetworkSolver(n)
    n.resolve_all_topology()
    x0 = s._build_x0()
    return dict(zip(s.unknown_names, x0, strict=True))


def _both(wall: bool) -> FlowNetwork:
    n = FlowNetwork()
    _gas(n)
    _coolant(n)
    if wall:
        n.add_wall(
            ThermalWall(
                "w",
                element_a="g0",
                element_b="c",
                layers=[WallLayer(thickness=1e-3, conductivity=20.0)],
                contact_area=1e-3,
            )
        )
    return n


def test_each_component_starts_as_it_would_alone() -> None:
    """Before #899 one global reference state mixed unrelated streams: the
    coolant's 358 K boundaries pulled the reference to 583 K and the gas
    paths started 19% high (0.094 against 0.079 kg/s here); the coolant's
    reference pressure was the gas outlet's 1.4 bar."""
    gas, cool = FlowNetwork(), FlowNetwork()
    _gas(gas)
    _coolant(cool)
    alone = {**_x0(gas), **_x0(cool)}
    joint = _x0(_both(wall=False))
    assert joint.keys() == alone.keys()
    for k in joint:
        assert joint[k] == alone[k], k


GROUP_T = (808.0 + 808.0 + 358.0 + 358.0) / 4


def _parts_at(T: float | None) -> dict[str, float]:
    """Each part's start alone, its boundaries at T (None: their own)."""
    gas, cool = FlowNetwork(), FlowNetwork()
    _gas(gas)
    _coolant(cool)
    if T is not None:
        for net in (gas, cool):
            for node in net.nodes.values():
                if isinstance(node, (PressureBoundary, MassFlowBoundary)):
                    node.Tt = T
    return {**_x0(gas), **_x0(cool)}


def test_a_wall_leaves_the_resistance_seed_at_each_streams_own_temperature() -> None:
    """Joined by a wall, the streams share their thermal group's mean Tt
    (#899), but the resistance-network seed -- here every pressure and flow
    -- runs at each component's own (#990): the joint start is each part's
    start alone. At the group's 583 K the coolant's pressures started 22%
    high (3.68 against 3.02 bar), its kernel's flux too low at that Tt; the
    group's Tt now enters only the seed's first Picard iterate (the reference
    m_dot), hence the tolerance."""
    own = _parts_at(None)
    group = _parts_at(GROUP_T)
    joint = _x0(_both(wall=True))
    assert joint.keys() == own.keys()
    for k in joint:
        assert joint[k] == pytest.approx(own[k], rel=1e-12), k
    assert group["c_in.P"] > 1.2 * own["c_in.P"]


def test_the_group_temperature_stays_the_reference_tt() -> None:
    """Everything else keeps the thermal group's mean (#899): the reference
    state carries it as T, and the component's own boundary mean as T_own."""
    s = NetworkSolver(_both(wall=True))
    comps = s._flow_components()
    group_T = s._thermal_group_T(comps)
    refs = []
    for ci, (nodes, elems) in enumerate(comps):
        s._ref_T_override = group_T[ci]
        try:
            with s._component_view(nodes, elems):
                refs.append(s._infer_reference_state())
        finally:
            s._ref_T_override = None
    assert [r["T"] for r in refs] == [GROUP_T, GROUP_T]
    assert [r["T_own"] for r in refs] == [808.0, 358.0]


def test_a_wall_does_not_join_two_flow_components() -> None:
    s = NetworkSolver(_both(wall=True))
    comps = s._flow_components()
    assert [sorted(c[0]) for c in comps] == [
        ["in", "m0", "m1", "out"],
        ["c_in", "c_mid", "c_out"],
    ]
    assert [c[1] for c in comps] == [["g0", "x0", "g1", "x1"], ["c", "c2"]]


def test_a_network_of_one_component_is_one_component() -> None:
    n = FlowNetwork()
    _gas(n)
    assert len(NetworkSolver(n)._flow_components()) == 1


def test_the_joint_network_solves_to_the_parts_roots() -> None:
    """Without the wall the two parts are independent problems: solved
    together they reach the roots each reaches alone."""
    gas, cool = FlowNetwork(), FlowNetwork()
    _gas(gas)
    _coolant(cool)
    rg = NetworkSolver(gas).solve()
    rc = NetworkSolver(cool).solve()
    r = NetworkSolver(_both(wall=False)).solve()
    assert rg["__success__"] and rc["__success__"] and r["__success__"], r["__message__"]
    for k in ("g0.m_dot", "g1.m_dot", "m0.P", "m1.P"):
        assert r[k] == pytest.approx(rg[k], rel=1e-6), k
    for k in ("c.m_dot", "c_in.Pt", "c_mid.P"):
        assert r[k] == pytest.approx(rc[k], rel=1e-6), k
    assert np.isfinite(r["__final_norm__"])
