"""A channel reports the friction factor its residual uses (#634).

The diagnostics applied the old Re floor sqrt(Re^2 + 64^2) after #634 had
moved the laminar branch to Poiseuille 64/Re with its own tiny re_floor, so
a nearly stagnant duct reported f = 1 while its row used 64/Re."""

from __future__ import annotations

import warnings

import pytest

import combaero as cb
from combaero.network import ChannelElement, FlowNetwork, NetworkSolver, PressureBoundary


def _duct(dp: float) -> dict:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=1e5 + dp, Tt=300.0))
    g.add_node(PressureBoundary("b", Pt=1e5, Tt=300.0))
    g.add_element(
        ChannelElement("c", "a", "b", length=1.0, diameter=0.002, regime="incompressible")
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve()
    assert r["__success__"], r["__message__"]
    return r


def test_a_laminar_duct_reports_poiseuille() -> None:
    r = _duct(100.0)
    assert 1.0 < r["c.Re"] < 1000.0
    assert r["c.f"] == pytest.approx(
        cb._core.friction_and_jacobian("haaland", r["c.Re"], 0.0).result[0], rel=1e-6
    )
    assert r["c.f"] == pytest.approx(64.0 / r["c.Re"], rel=1e-3)


def test_a_stagnant_duct_reports_the_floor_the_row_uses() -> None:
    r = _duct(0.0)
    assert r["c.Re"] == pytest.approx(0.0, abs=1e-6)
    f0 = cb._core.friction_and_jacobian("haaland", r["c.Re"], 0.0).result[0]
    assert r["c.f"] == pytest.approx(f0, rel=1e-9)
    assert r["c.f"] > 1e6
