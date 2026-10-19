"""A wall-heated compressible channel flags the heating its momentum omits (#551).

Where the duct discharges into a node fed by it alone, its heat is in the
momentum (Fanno-Rayleigh, test_channel_fanno_rayleigh_551.py). Into a node
that also mixes other streams, that node's Tt is not the duct's exit Tt and
the duct stays adiabatic Fanno flow at its inlet Tt. To first order Rayleigh
flow loses dPt / Pt = (gamma M^2 / 2) ln(Tt_exit / Tt_in) (cooling gains it),
which the diagnostics then report and the Inspector flags above 0.1% of the
inlet Pt.
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
                f"{side}_o", f"{side}_p", f"{side}_out", Cd=0.8, diameter=0.011, correlation="fixed"
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


def test_a_heated_and_a_cooled_duct_into_mixing_nodes_report_their_omission() -> None:
    """The heated side omits a loss, the cooled side a gain; the estimate is
    the first-order formula on the duct's own Mach numbers and end Tt."""
    r = _solve(_mixing_pair())
    d = r["__element_diag__"]
    assert d["cold"]["rayleigh_omitted_frac"] > 1e-3
    assert d["hot"]["rayleigh_omitted_frac"] < -1e-3
    assert d["cold"]["rayleigh_omitted"] == 1.0 and d["hot"]["rayleigh_omitted"] == 1.0
    assert "rayleigh_in_momentum" not in d["cold"]
    m2 = 0.5 * (d["cold"]["M_in_duct"] ** 2 + d["cold"]["M_exit_duct"] ** 2)
    expected = 0.5 * d["cold"]["gamma_ref"] * m2 * math.log(r["cold_p.T"] / 350.0)
    assert d["cold"]["rayleigh_omitted_frac"] == pytest.approx(expected, rel=1e-9)


def test_a_duct_without_a_wall_reports_nothing() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2e5, Tt=300.0))
    g.add_node(PressureBoundary("b", Pt=1e5, Tt=300.0))
    g.add_element(ChannelElement("c", "a", "b", length=1.0, diameter=0.01, regime="compressible"))
    d = _solve(g)["__element_diag__"]["c"]
    assert "rayleigh_omitted" not in d
    assert "rayleigh_in_momentum" not in d
