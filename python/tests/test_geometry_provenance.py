"""Geometry-inheritance provenance is exposed in element diagnostics (#559).

Every element that inherits or derives its flow area now records where that
value came from ("user" / "inherited from ..." / an `_infer_flow_area` source
string) and surfaces it as `area_source` (or `F0_source`/`F1_source` for an
area change) in `diagnostics()`. This is the authoritative provenance the GUI
reads instead of re-walking the graph itself.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    AreaChangeElement,
    ChannelElement,
    FlowNetwork,
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.mpce_element import MultiPortChamberElement

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _solve(net: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkSolver(net).solve(timeout=20)


def _diag(result: dict, eid: str) -> dict:
    return result.get("__element_diag__", {}).get(eid, {})


def test_channel_reports_user_vs_inherited_source() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=2e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PlenumNode("p"))
    net.add_node(PressureBoundary("out", Pt=1e5, Tt=300.0, Y=Y_AIR))
    net.add_element(ChannelElement("ch_user", "in", "p", length=1.0, diameter=0.1))
    net.add_element(ChannelElement("ch_inh", "p", "out", length=1.0))  # no diameter

    r = _solve(net)
    assert r["__success__"], r.get("__message__")
    assert _diag(r, "ch_user")["area_source"] == "user"
    # The inheriting channel names the element it took its diameter from.
    assert "ch_user" in _diag(r, "ch_inh")["area_source"]
    # And it actually inherited the value.
    assert _diag(r, "ch_inh")["diameter"] == 0.1


def test_area_change_reports_per_face_source() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=2e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PlenumNode("p"))
    net.add_node(PlenumNode("p2"))
    net.add_node(PressureBoundary("out", Pt=1e5, Tt=300.0, Y=Y_AIR))
    net.add_element(ChannelElement("ch_u", "in", "p", length=1.0, diameter=0.1))
    net.add_element(AreaChangeElement("ac", "p", "p2"))  # both faces inherit
    net.add_element(ChannelElement("ch_d", "p2", "out", length=1.0, diameter=0.08))

    r = _solve(net)
    assert r["__success__"], r.get("__message__")
    d = _diag(r, "ac")
    assert "ch_u" in d["F0_source"]
    assert "ch_d" in d["F1_source"]
    assert d["F0"] == pytest.approx(math.pi / 4 * 0.1**2)
    assert d["F1"] == pytest.approx(math.pi / 4 * 0.08**2)


def test_orifice_reports_user_vs_default_bore() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=2e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PlenumNode("p"))
    net.add_node(PressureBoundary("out", Pt=1e5, Tt=300.0, Y=Y_AIR))
    net.add_element(OrificeElement("orf_u", "in", "p", diameter=0.05, Cd=0.6, correlation="fixed"))
    net.add_element(OrificeElement("orf_d", "p", "out", Cd=0.6, correlation="fixed"))

    r = _solve(net)
    assert r["__success__"], r.get("__message__")
    assert _diag(r, "orf_u")["area_source"] == "user"
    # No diameter/area given: the bore falls back to DEFAULT_ORIFICE_BORE (0.08),
    # reported as "default" (it cannot be inferred from a neighbour).
    assert _diag(r, "orf_d")["area_source"] == "default"
    assert _diag(r, "orf_d")["diameter"] == pytest.approx(0.08)


def test_junction_port_reports_inherited_vs_user_area() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_in", Pt=2.1e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PressureBoundary("pb_str", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PressureBoundary("pb_bra", Pt=1.9e5, Tt=300.0, Y=Y_AIR))
    for pid in ("port_com", "port_str", "port_bra"):
        net.add_node(MomentumChamberNode(pid, area=0.01))
    net.add_element(ChannelElement("ch_in", "pb_in", "port_com", length=0.3, diameter=0.05))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    net.add_element(LosslessConnectionElement("lc_bra", "port_bra", "pb_bra"))
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["port_com"],
            outlet_nodes=["port_str", "port_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[None, 0.01, 0.02],  # common inherits; the two outlets are user
            flow_direction="branch",
            strict=False,
        )
    )
    net.resolve_all_topology()
    jct = net.elements["jct"]
    # The idempotent per-pass resolve must label the channel-fed port as inherited
    # (not "user", which an earlier pass's already-set value would wrongly imply).
    assert jct._port_area_sources[0] == "inherited from ch_in"
    assert jct._port_area_sources[1] == "user"
    assert jct._port_area_sources[2] == "user"
