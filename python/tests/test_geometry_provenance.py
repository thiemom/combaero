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
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

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
