"""BorderCarnotLossElement: deprecated in 0.7.0, for removal in 0.8.0.

Audit (2026-10-10, after #272):
- **Not useful.** Its residual, ``Pt_in - Pt_out - L * 0.5 * rho_in * u_in^2``, is
  exactly a constant head loss. The only thing it adds is the formula
  ``L = 4 * (1 - cos(0.75 * delta))^2``, which was derived for a junction's
  lateral branch. Both junction closures now carry that physics, and adding the
  element on top of them double-counts it (Bassett K6 at q = 0.5: 0.867 -> 1.249).
  It was never compared against bend data.
- **Not used in the repo.** No production module, GUI node or example
  instantiates it; only tests did.
- **But public since 0.4.0** and shipped on PyPI, so outside use cannot be ruled
  out. It is therefore deprecated with an exact migration rather than deleted
  outright.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    BorderCarnotLossElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
    PressureLossElement,
)
from combaero.network.pressure_loss import ConstantHeadLoss

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_A = 0.005


def test_constructing_it_warns_with_the_replacement():
    with pytest.warns(DeprecationWarning, match="ConstantHeadLoss"):
        BorderCarnotLossElement("x", "a", "b", delta_geom_deg=90.0, area=_A)


def _solve(kind: str, orifice_scale: float) -> tuple[float, float]:
    zeta = 4.0 * (1.0 - math.cos(0.75 * math.pi / 2.0)) ** 2
    net = FlowNetwork()
    net.add_node(PressureBoundary("src", Pt=2.0e5, Tt=300.0, Y=_Y))
    net.add_node(MomentumChamberNode("mid", area=_A))
    net.add_node(PressureBoundary("snk", Pt=1.0e5, Tt=300.0, Y=_Y))
    if kind == "deprecated":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            loss = BorderCarnotLossElement("x", "src", "mid", delta_geom_deg=90.0, area=_A)
    else:
        loss = PressureLossElement(
            "x", "src", "mid", correlation=ConstantHeadLoss(zeta=zeta, area=_A), area=_A
        )
    net.add_element(loss)
    diameter = math.sqrt(4.0 * _A / math.pi) * orifice_scale
    net.add_element(
        OrificeElement("o", "mid", "snk", Cd=0.8, diameter=diameter, regime="compressible")
    )
    sol = NetworkSolver(net).solve()
    assert sol["__success__"]
    return sol["x.m_dot"], 2.0e5 - sol["mid.Pt"]


@pytest.mark.parametrize("orifice_scale", [0.3, 0.6, 0.9])
def test_the_replacement_is_exact(orifice_scale):
    """The migration in the warning reproduces the element to round-off, from
    a 200 Pa to a 14 kPa drop."""
    m_old, dp_old = _solve("deprecated", orifice_scale)
    m_new, dp_new = _solve("replacement", orifice_scale)
    assert m_new == pytest.approx(m_old, rel=1e-9)
    assert dp_new == pytest.approx(dp_old, rel=1e-9)
