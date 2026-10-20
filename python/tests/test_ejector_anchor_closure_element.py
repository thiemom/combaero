"""EjectorElement(anchor_closure=...): the (P_c*, P_b0) pair the element solves
with, at network level, on the two MEASURED air rigs scored cross-source
(validation/ejector: Henry 2007, read as gauge; Chen 2013's rectangular rig,
read as absolute). Accuracy, labelled cross-source: neither rig is a source of
either closure. The closure-level numbers are in anchor_scorecard.py; these
tests pin that the solved element carries them.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.ejector_element import EjectorElement
from validation.ejector.data import chen2013_figs as chen
from validation.ejector.data import henry2007_fig5 as henry

Y = list(cb.species.dry_air_mass())
ATM = 101325.0
T = 293.15


def _disk(d: float) -> float:
    return math.pi * d * d / 4.0


def _network(
    pp: float, pb: float, a_t: float, a_e: float, a_3: float, d_pipe: float, closure: str
) -> FlowNetwork:
    g = FlowNetwork()
    for nid, pt in (("bp", pp), ("bs", ATM), ("bo", pb)):
        g.add_node(PressureBoundary(nid, Pt=pt, Tt=T, Y=Y))
    for nid in ("mp", "ms", "mo"):
        g.add_node(MomentumChamberNode(nid, area=_disk(d_pipe)))
    for eid, a, b in (("cp", "bp", "mp"), ("cs", "bs", "ms"), ("co", "mo", "bo")):
        g.add_element(
            ChannelElement(eid, a, b, length=0.05, diameter=d_pipe, regime="incompressible")
        )
    g.add_element(
        EjectorElement(
            id="ej",
            primary_node="mp",
            secondary_node="ms",
            outlet_node="mo",
            throat_area=a_t,
            nozzle_exit_area=a_e,
            mixing_area=a_3,
            regime="compressible",
            anchor_closure=closure,
        )
    )
    return g


def _rigs() -> list[tuple[str, float, list[tuple[float, float]], float, float, float, float]]:
    gm = henry.GEOMETRY_MM
    rigs = [
        (
            f"Henry {bar} bar gauge",
            bar * 1e5 + ATM,
            [(x * ATM, w) for x, w in henry.FIG5_OMEGA_VS_PB[bar]],
            _disk(gm["d_throat"] * 1e-3),
            _disk(gm["d_exit"] * 1e-3),
            _disk(gm["D_mix"] * 1e-3),
            0.03,
        )
        for bar in henry.PRIMARY_PRESSURES_BAR
    ]
    a_t = 303.5e-6  # 6.07 mm x 50 mm
    for pp in sorted({r[2] for r in chen.FIG7A}):
        pts = sorted((r[3] * 1e5, r[4]) for r in chen.FIG7A if r[1] == "measured" and r[2] == pp)
        rigs.append((f"Chen rect {pp:g} bar", pp * 1e5, pts, a_t, a_t * 8.35 / 6.07, 1353e-6, 0.15))
    return rigs


def _curve_rms(closure: str, rig: tuple) -> float:
    _, pp, pts, a_t, a_e, a_3, d = rig
    errs = []
    for pb, omega in pts:
        r = NetworkSolver(_network(pp, pb, a_t, a_e, a_3, d, closure)).solve(timeout=30.0)
        assert r["__success__"], (rig[0], pb)
        errs.append(r["cs.m_dot"] / r["cp.m_dot"] - omega)
    return math.sqrt(sum(e * e for e in errs) / len(errs)) / pts[0][1]


@pytest.mark.parametrize("rig", _rigs(), ids=lambda r: r[0])
def test_huang_anchors_solve_every_measured_point_with_half_the_curve_error(rig) -> None:
    """Every measured back pressure solves with either pair. The "huang1999"
    pair's curve RMS (over the plateau) is 12-23%, against 29-48% for
    "kracik_dvorak", on all eight curves."""
    huang = _curve_rms("huang1999", rig)
    kd = _curve_rms("kracik_dvorak", rig)
    assert huang < 0.25, rig[0]
    assert kd > 0.28, rig[0]
    assert huang < kd / 1.4, rig[0]


def test_unknown_anchor_closure_is_refused() -> None:
    with pytest.raises(ValueError, match="anchor_closure"):
        _network(5e5, 1.2e5, 1e-5, 2e-5, 6e-5, 0.03, "lossless")


def test_a_failed_solve_with_the_outlet_above_the_dead_head_says_so() -> None:
    """Above P_b0 a real ejector backflows, outside the model: the solver's
    message names the dead head instead of only reporting no progress."""
    rig = _rigs()[0]
    _, pp, _, a_t, a_e, a_3, d = rig
    g = _network(pp, 3.0 * ATM, a_t, a_e, a_3, d, "huang1999")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve(timeout=30.0)
    assert not r["__success__"]
    assert "above the dead-head back pressure" in r["__message__"]


def test_the_compressible_warm_start_proxy_takes_the_lossless_window() -> None:
    """An incompressible ejector uses no anchors, but its compressible
    warm-start proxy does: the proxy runs with kracik_dvorak (whose wider
    window keeps a forward root where the Huang dead head has none), and the
    element gets its own closure and regime back afterwards."""
    g = _network(5e5, ATM, 1e-5, 2e-5, 6e-5, 0.03, "huang1999")
    ej = g.elements["ej"]
    ej.regime = "incompressible"
    g.resolve_all_topology()
    s = NetworkSolver(g)
    seen = []
    inner = s._solve_impl

    def spy(*args, **kwargs):
        seen.append((ej.regime, ej.anchor_closure))
        return inner(*args, **kwargs)

    s._solve_impl = spy
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s._compressible_seed(timeout=10.0)
    assert seen == [("compressible", "kracik_dvorak")]
    assert (ej.regime, ej.anchor_closure) == ("incompressible", "huang1999")
