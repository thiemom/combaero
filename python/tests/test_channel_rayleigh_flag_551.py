"""A wall-heated compressible channel's heat is in its momentum wherever it
discharges (#551, #672), and the omission flag is only a fallback.

Into a node fed by it alone, that node's Tt is the duct's exit Tt (#551). Into
a node that also mixes other streams, or a boundary, the solver forms the
duct's own exit state from its stream and its own wall heat (#672). Where no
exit state can be formed (no flow into the receiving node), the duct is solved
as adiabatic Fanno flow at its inlet Tt. The diagnostic then reports the
first-order Rayleigh change, dPt / Pt = (gamma M^2 / 2) ln(Tt_exit / Tt_in)
(cooling gains it), and the Inspector flags it above 0.1% of the inlet Pt.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _solve(g: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve(timeout=30)
    assert r["__success__"], r["__message__"]
    return r


def _mixing_pair() -> FlowNetwork:
    """A hot and a cold compressible duct on one wall; each discharges into a
    plenum that also takes a small bleed, so neither plenum's Tt is a duct's
    exit Tt."""
    g = FlowNetwork()
    for side, Tt in (("hot", 1500.0), ("cold", 350.0)):
        g.add_node(PressureBoundary(f"{side}_in", Pt=4.0e5, Tt=Tt, Y=Y_AIR))
        g.add_node(PressureBoundary(f"{side}_out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
        g.add_node(PlenumNode(f"{side}_p"))
        g.add_node(MassFlowBoundary(f"{side}_bleed", m_dot=2e-4, Tt=Tt, Y=Y_AIR))
        g.add_element(
            ChannelElement(
                side,
                f"{side}_in",
                f"{side}_p",
                length=0.4,
                diameter=0.01,
                roughness=0.0,
                regime="compressible",
                surface=ConvectiveSurface(area=math.pi * 0.01 * 0.4),
            )
        )
        g.add_element(
            OrificeElement(
                f"{side}_o",
                f"{side}_p",
                f"{side}_out",
                Cd=0.8,
                diameter=0.011,
                correlation="fixed",
                in_pipe=False,
            )
        )
        g.add_element(
            OrificeElement(
                f"{side}_b",
                f"{side}_bleed",
                f"{side}_p",
                Cd=0.8,
                diameter=0.003,
                correlation="fixed",
            )
        )
    g.add_wall(
        ThermalWall(id="w", element_a="hot", element_b="cold", layers=[WallLayer(0.0005, 25.0)])
    )
    return g


def test_a_heated_and_a_cooled_duct_into_mixing_nodes_carry_their_heat() -> None:
    """#672: neither plenum's Tt is a duct's exit Tt, and both ducts still
    carry their heat in the momentum; nothing is flagged."""
    r = _solve(_mixing_pair())
    d = r["__element_diag__"]
    for side in ("hot", "cold"):
        assert d[side]["rayleigh_in_momentum"] == 1.0
        assert "rayleigh_omitted" not in d[side]
    assert r["__second_law_failures__"] == []


def test_the_fallback_estimate_is_the_first_order_rayleigh_formula() -> None:
    """Where no exit state exists the duct reports the omitted change: the
    heated side a loss, the cooled side a gain, on its own Mach numbers."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = NetworkSolver(_mixing_pair())
        r = s.solve(timeout=30)
    assert r["__success__"], r["__message__"]
    for side, sign in (("cold", 1.0), ("hot", -1.0)):
        el = s.network.elements[side]
        st_in = s._get_node_state(s.network.nodes[f"{side}_in"], s._last_x)
        st_out = s._get_node_state(s.network.nodes[f"{side}_p"], s._last_x)
        st_in.m_dot = r[f"{side}.m_dot"]
        _, _, flow = el._compressible_flow(st_in, st_out, 1.0, derivatives=False)
        ref = {"gamma_ref": 1.4}
        est = el._rayleigh_omission(st_in, st_out, flow, ref)
        m2 = 0.5 * (float(flow.M_in) ** 2 + float(flow.M_exit) ** 2)
        expected = 0.5 * ref["gamma_ref"] * m2 * math.log(float(st_out.Tt) / float(st_in.Tt))
        assert est["rayleigh_omitted_frac"] == pytest.approx(expected, rel=1e-12)
        assert sign * est["rayleigh_omitted_frac"] > 1e-3


def test_a_duct_without_a_wall_reports_nothing() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2e5, Tt=300.0))
    g.add_node(PressureBoundary("b", Pt=1e5, Tt=300.0))
    g.add_element(ChannelElement("c", "a", "b", length=1.0, diameter=0.01, regime="compressible"))
    d = _solve(g)["__element_diag__"]["c"]
    assert "rayleigh_omitted" not in d
    assert "rayleigh_in_momentum" not in d
