"""RESIDUAL-ROW-SCALING: a strict solve checks every row's declared units
against its Jacobian at the converged solution.

#617 made every element declare its rows' units; this checks the declaration
is right. A Pa row carries a coefficient of order one on a pressure unknown, a
kg/s row on a flow unknown, and a misfiled row shows the other family. The
suite runs strict (conftest.py), so every converged solve in it is checked.
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
from combaero.network.pressure_loss import ConstantFractionLoss
from combaero.network.solver import RowKindMismatchError

Y = list(cb.mole_to_mass(cb.species.dry_air()))


class _LossFiledAsFlow(PressureLossElement):
    residual_scale_kind = "mdot"  # wrong: the residual is a Pt balance in Pa


class _OrificeFiledAsPressure(OrificeElement):
    residual_scale_kind = "p"  # wrong: the residual is m - m_calc in kg/s


def _net(loss_cls: type = PressureLossElement, orifice_cls: type = OrificeElement) -> FlowNetwork:
    g = FlowNetwork()
    for name, Pt in (("a", 1.2e5), ("b", 1.0e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 300.0, Y
        g.add_node(bnd)
    for n in ("p1", "p2", "p3"):
        g.add_node(PlenumNode(n))
    g.add_element(orifice_cls("o1", "a", "p1", Cd=0.8, diameter=0.02, correlation="fixed"))
    g.add_element(loss_cls("loss", "p1", "p2", correlation=ConstantFractionLoss(xi=0.02)))
    g.add_element(AreaChangeElement("ac", "p2", "p3", F0=3e-4, F1=6e-4))
    g.add_element(OrificeElement("o2", "p3", "b", Cd=0.8, diameter=0.02, correlation="fixed"))
    return g


def test_correctly_declared_rows_agree() -> None:
    s = NetworkSolver(_net())
    r = s.solve()
    assert r["__success__"]
    report = s.check_row_kinds(np.array(r["__x_solution__"]))
    assert len(report) == len(r["__x_solution__"])
    verdicts = {(rec["owner"], rec["verdict"]) for rec in report}
    assert ("PressureLossElement 'loss'", "agree") in verdicts
    assert ("AreaChangeElement 'ac'", "agree") in verdicts
    assert ("OrificeElement 'o1'", "agree") in verdicts
    assert ("node 'p1' mass", "agree") in verdicts
    assert not [rec for rec in report if rec["verdict"] == "DISAGREE"]


def test_the_vortex_row_agrees() -> None:
    """The one shipped element no other suite network converges with."""
    g = FlowNetwork()
    for name, Pt in (("a", 1.2e5), ("b", 1.0e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 300.0, Y
        g.add_node(bnd)
    for n in ("p1", "p2"):
        g.add_node(PlenumNode(n))
    g.add_element(OrificeElement("o1", "a", "p1", Cd=0.8, diameter=0.02, correlation="fixed"))
    g.add_element(VortexElement("vx", "p1", "p2", omega_rpm=3000.0))
    g.add_element(OrificeElement("o2", "p2", "b", Cd=0.8, diameter=0.02, correlation="fixed"))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    report = s.check_row_kinds(np.array(r["__x_solution__"]))
    assert [rec["verdict"] for rec in report if rec["owner"] == "VortexElement 'vx'"] == ["agree"]


@pytest.mark.parametrize(
    ("kw", "owner", "declared"),
    [
        ({"loss_cls": _LossFiledAsFlow}, "_LossFiledAsFlow 'loss'", "mdot"),
        ({"orifice_cls": _OrificeFiledAsPressure}, "_OrificeFiledAsPressure 'o1'", "p"),
    ],
    ids=["Pa-row-as-mdot", "kgs-row-as-p"],
)
def test_a_misfiled_row_raises_in_a_strict_solve(kw: dict, owner: str, declared: str) -> None:
    s = NetworkSolver(_net(**kw))
    s.raise_internal_errors = True
    with pytest.raises(RowKindMismatchError, match=f"{owner}: declared '{declared}'"):
        s.solve()


def test_a_non_strict_solve_does_not_check() -> None:
    """Production solves pay nothing; the report stays available on demand."""
    s = NetworkSolver(_net(loss_cls=_LossFiledAsFlow))
    s.raise_internal_errors = False
    r = s.solve()
    assert r["__success__"]
    bad = [
        rec
        for rec in s.check_row_kinds(np.array(r["__x_solution__"]))
        if rec["verdict"] == "DISAGREE"
    ]
    assert [(rec["owner"], rec["measured"]) for rec in bad] == [("_LossFiledAsFlow 'loss'", "p")]


def test_the_check_reads_the_same_kinds_the_scales_do() -> None:
    """One list (_row_kinds) feeds both, so the check tests what is weighted."""
    s = NetworkSolver(_net())
    r = s.solve()
    x = np.array(r["__x_solution__"])
    ref_p, ref_m = s._reference_scales(x)
    expected = [ref_p if kind == "p" else ref_m for _owner, kind in s._row_kinds(x)]
    np.testing.assert_array_equal(s._build_residual_scales(x), expected)
