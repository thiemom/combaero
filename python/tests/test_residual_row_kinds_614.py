"""RESIDUAL-ROW-SCALING (#614): elements whose residual is a pressure in Pa
must be scaled as pressure rows.

PressureLossElement and VortexElement write a
total-pressure balance (as AreaChangeElement did until #731 made its row a
mass flow), but the default in
``_build_residual_scales`` filed them as mass-flow rows: weighted
ref_p / ref_mdot (1e5-1e9) too heavily in the norm hybr/LM minimise. On a
stress network the LM fallback then walked every pressure to -4e14 Pa to
shrink those rows, where the Jacobian went non-finite (how #614 surfaced).
"""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    AreaChangeElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    VortexElement,
)
from combaero.network import components as comp
from combaero.network.pressure_loss import ConstantFractionLoss

Y = list(cb.mole_to_mass(cb.species.dry_air()))


@pytest.mark.parametrize("cls", [PressureLossElement, VortexElement])
def test_a_pressure_balance_element_declares_pressure_rows(cls: type) -> None:
    assert cls.residual_scale_kind == "p"


def test_the_area_change_row_is_a_mass_flow() -> None:
    """m - m_calc since #731: its row is in kg/s."""
    assert AreaChangeElement.residual_scale_kind == "mdot"


def _net() -> FlowNetwork:
    g = FlowNetwork()
    for name, Pt in (("a", 1.2e5), ("b", 1.0e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 300.0, Y
        g.add_node(bnd)
    for n in ("p1", "p2", "p3"):
        g.add_node(PlenumNode(n))
    g.add_element(OrificeElement("o1", "a", "p1", Cd=0.8, diameter=0.02, correlation="fixed"))
    g.add_element(
        PressureLossElement("loss", "p1", "p2", correlation=ConstantFractionLoss(xi=0.02))
    )
    g.add_element(AreaChangeElement("ac", "p2", "p3", F0=3e-4, F1=6e-4))
    g.add_element(OrificeElement("o2", "p3", "b", Cd=0.8, diameter=0.02, correlation="fixed"))
    return g


def test_their_rows_get_the_pressure_scale(monkeypatch: pytest.MonkeyPatch) -> None:
    """The loss adds one ref_p row; with the old default it was one more
    ref_mdot row. (The area change's row is a mass flow since #731.)"""
    s = NetworkSolver(_net())
    r = s.solve()
    assert r["__success__"]
    x0 = np.array(r["__x_solution__"])
    ref_p, ref_m = s._reference_scales(x0)
    n_p = int(np.sum(s._build_residual_scales(x0) == ref_p))
    monkeypatch.setattr(comp.PressureLossElement, "residual_scale_kind", "mdot")  # old default
    D_old = s._build_residual_scales(x0)
    assert int(np.sum(D_old == ref_p)) == n_p - 1
    assert ref_p / ref_m > 1e5  # what the misfiling cost, here
