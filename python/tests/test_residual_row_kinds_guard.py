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
    CrossflowSegmentElement,
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


# Over 20 s on one worker: runs in CI's parallel slow job (see
# validation-tests.yml), still before every merge.
pytestmark = pytest.mark.slow


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


def _bleed_station(d_lat: float) -> FlowNetwork:
    """A compressible bleed station on two 1 um crossflow segments: their
    kg/s rows carry dm/dPt ~ 1.7 kg/(s Pa) and no flow coefficient near 1."""
    g = FlowNetwork()
    for name, Pt in (("in", 5.01e5), ("out", 5.0e5), ("lat", 4.99e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 300.0, Y
        g.add_node(bnd)
    g.add_node(PlenumNode("b"))
    kappa = cb.STATION_KAPPA_BLEED_BASSETT
    g.add_element(
        CrossflowSegmentElement(
            "s0",
            "in",
            "b",
            length=1e-6,
            area=1e-3,
            entry_K=0.0,
            to_kappa=kappa,
            next_seg="s1",
            regime="compressible",
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s1",
            "b",
            "out",
            length=1e-6,
            area=1e-3,
            from_kappa=kappa,
            prev_seg="s0",
            regime="compressible",
        )
    )
    g.add_element(
        OrificeElement(
            "bleed", "b", "lat", diameter=d_lat, Cd=0.6, correlation="fixed", regime="compressible"
        )
    )
    return g


def test_a_gain_near_one_is_not_a_unit_coefficient() -> None:
    """The first rule (any coefficient in [0.5, 2]) read these kg/s rows as
    Pa and a strict solve raised on a correct network (stress sweep,
    2026-10-17). Unit coefficients must be exact."""
    band_flags = 0
    for d_lat in (0.005, 0.01, 0.03):
        s = NetworkSolver(_bleed_station(d_lat))
        s.raise_internal_errors = True
        r = s.solve()  # raised RowKindMismatchError under the band rule
        assert r["__converged__"]
        report = s.check_row_kinds(np.array(r["__x_solution__"]))
        assert not [rec for rec in report if rec["verdict"] == "DISAGREE"]
        _, J = s._residuals_and_jacobian(np.array(r["__x_solution__"]))
        J = J.tocsr()
        names = s.unknown_names
        for rec in report:
            if not rec["owner"].startswith("CrossflowSegmentElement"):
                continue
            row = J.getrow(rec["row"])
            in_band = {
                names[j].rsplit(".", 1)[-1].startswith("m_dot"): True
                for j, v in zip(row.indices, row.data, strict=True)
                if 0.5 <= abs(v) <= 2.0
            }
            band_flags += in_band == {False: True}  # only a pressure in band
    assert band_flags > 0, "the fixture no longer exercises the trap"
