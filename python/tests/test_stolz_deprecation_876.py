"""The Stolz metering Cd (ISO 5167:1980) is deprecated in 0.8.0 (#876).

ISO replaced Stolz with Reader-Harris/Gallagher in ISO 5167-1:1991/Amd 1:1998
(and ISO 5167-2:2003). Stolz stays until 0.9.0, with its values unchanged;
every path that selects it warns once per use, naming the replacement:

- ``OrificeElement(correlation="Stolz")``: at construction, never per
  residual (the element evaluates through the quiet ``_metering_*``
  bindings);
- the bindings that take a ``MeteringCdCorrelation`` (``metering_cd_and_dRe``,
  ``metering_cd_in_range``, ``solve_orifice_mdot``, ``orifice_flow``,
  ``orifice_flow_state``);
- the GUI: a saved network with a Stolz orifice still loads and solves, and
  the solve response carries a deprecation note for it.

A cooling plate that refuses the whole ISO 5167 family does not also warn.
"""

from __future__ import annotations

import asyncio
import math
import warnings

import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    EffusionPlateElement,
    FlowNetwork,
    ImpingementPlateElement,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

M = cb.MeteringCdCorrelation
MATCH = r"(?s)Stolz.*deprecated.*0\.9\.0.*ReaderHarrisGallagher"


def _geom(beta: float = 0.5, D: float = 0.1) -> cb.OrificeGeometry:
    g = cb.OrificeGeometry()
    g.d, g.D = beta * D, D
    return g


def _stolz_warnings(record: list[warnings.WarningMessage]) -> list[warnings.WarningMessage]:
    return [
        w
        for w in record
        if issubclass(w.category, DeprecationWarning) and "Stolz" in str(w.message)
    ]


# --- the element ------------------------------------------------------------


def test_the_element_warns_at_construction_and_names_the_replacement() -> None:
    with pytest.warns(DeprecationWarning, match=MATCH) as rec:
        OrificeElement("o", "a", "b", diameter=0.05, correlation="Stolz")
    assert len(_stolz_warnings(list(rec))) == 1
    # attributed to the caller, not to components.py
    assert rec[0].filename == __file__


def test_reader_harris_gallagher_does_not_warn() -> None:
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        OrificeElement("o", "a", "b", diameter=0.05)
        OrificeElement("o", "a", "b", diameter=0.05, correlation="ReaderHarrisGallagher")
        cb.metering_cd_and_dRe(M.ReaderHarrisGallagher, _geom(), 1e5)
    assert _stolz_warnings(rec) == []


def _net(correlation: str) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.5e5, Tt=300.0))
    g.add_node(PlenumNode("hp"))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    from combaero.network import ChannelElement

    g.add_element(ChannelElement("pipe", "hi", "hp", length=0.5, diameter=0.1))
    g.add_element(OrificeElement("o", "hp", "ho", diameter=0.05, correlation=correlation))
    return g


def test_a_stolz_solve_warns_once_and_matches_the_kernel() -> None:
    """One warning for the whole build-and-solve, none per residual; the
    solved Cd is the C++ Stolz kernel's at the solved pipe Reynolds number."""
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        g = _net("Stolz")
        r = NetworkSolver(g).solve()
    assert r["__success__"], r["__message__"]
    assert len(_stolz_warnings(rec)) == 1
    el = g.elements["o"]
    diag = r["__element_diag__"]["o"]
    # The element's Re_D: the pipe's, with mu at the inlet (Tt, P)
    X = list(cb.species.dry_air())
    mu = cb.transport_state(diag["Tt_in"], diag["P_in"], X).mu
    Re_D = 4.0 * r["o.m_dot"] / (math.pi * 0.1 * mu)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        cd_kernel = cb.metering_cd_and_dRe(M.Stolz, el._orifice_geom, Re_D)[0]
    assert diag["Cd"] == pytest.approx(cd_kernel, rel=1e-6)


# --- the quiet path agrees with the public one bit for bit -------------------


@pytest.mark.parametrize("Re_D", [3e3, 1e4, 1e5, 1e7])
@pytest.mark.parametrize("beta", [0.2, 0.5, 0.75])
def test_the_quiet_bindings_agree_bit_for_bit(beta: float, Re_D: float) -> None:
    g = _geom(beta)
    with pytest.warns(DeprecationWarning, match=MATCH):
        loud = cb.metering_cd_and_dRe(M.Stolz, g, Re_D)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        quiet = _solver_tools.metering_cd_and_dRe_quiet(M.Stolz, g, Re_D)
        quiet_range = _solver_tools.metering_cd_in_range_quiet(M.Stolz, g, Re_D)
    assert loud == quiet
    with pytest.warns(DeprecationWarning, match=MATCH):
        assert cb.metering_cd_in_range(M.Stolz, g, Re_D) == quiet_range


# --- the bindings -----------------------------------------------------------


def _calls() -> dict[str, object]:
    g = _geom()
    X = list(cb.species.dry_air())
    return {
        "metering_cd_and_dRe": lambda c: cb.metering_cd_and_dRe(c, g, 1e5),
        "metering_cd_in_range": lambda c: cb.metering_cd_in_range(c, g, 1e5),
        "solve_orifice_mdot": lambda c: cb.solve_orifice_mdot(g, 1e4, 1.2, 1.8e-5, correlation=c),
        "orifice_flow": lambda c: cb.orifice_flow(g, 1e4, 300.0, 1e5, 1.8e-5, correlation=c).mdot,
        "orifice_flow_state": lambda c: (
            cb.orifice_flow_state(g, 1e4, 300.0, 1e5, X, correlation=c).mdot
        ),
    }


@pytest.mark.parametrize("name", sorted(_calls()))
def test_each_binding_warns_on_stolz_only(name: str) -> None:
    call = _calls()[name]
    with pytest.warns(DeprecationWarning, match=MATCH):
        call(M.Stolz)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        call(M.ReaderHarrisGallagher)
    assert _stolz_warnings(rec) == []


def test_the_bindings_still_compute_stolz() -> None:
    """Deprecated, not changed: Stolz's own formula, not a forwarded RHG.
    ISO 5167:1980 at corner taps (L1 = L2' = 0):
    0.5959 + 0.0312 b^2.1 - 0.184 b^8 + 0.0029 b^2.5 (1e6/Re_D)^0.75, with
    0.0029 * 1e6^0.75 carried as 91.71 (orifice.h stolz::C3)."""
    g = _geom(0.5)
    g.tap = cb.OrificeTap.Corner
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        stolz = cb.metering_cd_and_dRe(M.Stolz, g, 1e5)[0]
    b = 0.5
    want = 0.5959 + 0.0312 * b**2.1 - 0.184 * b**8 + 91.71 * b**2.5 * 1e5**-0.75
    assert stolz == pytest.approx(want, rel=1e-12)
    assert stolz != cb.metering_cd_and_dRe(M.ReaderHarrisGallagher, g, 1e5)[0]


# --- plates that refuse the family do not also warn ---------------------------


def test_a_refusing_impingement_plate_refuses_without_the_warning() -> None:
    d = 1e-3
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        p = ImpingementPlateElement(
            "p",
            "s",
            "c",
            d_jet=d,
            xn_d=5.0,
            yn_d=4.0,
            z_d=2.0,
            span=12 * 4.0 * d,
            plate_thickness=d,
            correlation="Stolz",
        )
    assert _stolz_warnings(rec) == []
    g = FlowNetwork()
    g.add_node(PlenumNode("s"))
    g.add_node(PlenumNode("c"))
    g.add_element(p)
    p.resolve_topology(g)
    with pytest.raises(ValueError, match="no\\s+pipe"):
        p.validate()


def test_a_refusing_effusion_plate_refuses_without_the_warning() -> None:
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        e = EffusionPlateElement(
            "p",
            "c",
            "g",
            panel_area=1.0,
            correlation="Stolz",
            hole_diameter=2.16e-3,
            wall_thickness=6.3e-3,
            pitch=15.24e-3,
        )
    assert _stolz_warnings(rec) == []
    with pytest.raises(ValueError, match="NORMED metering correlation"):
        e.validate()


# --- the GUI ----------------------------------------------------------------


def _gui_graph(correlation: str) -> dict:
    def node(i: str, t: str, d: dict, x: int) -> dict:
        return {"id": i, "type": t, "data": {"id": i, **d}, "position": {"x": x, "y": 0}}

    def edge(i: str, s: str, t: str) -> dict:
        return {
            "id": i,
            "source": s,
            "sourceHandle": "flow-source",
            "target": t,
            "targetHandle": "flow-target",
            "data": {},
        }

    return {
        "nodes": [
            node("pb_in", "pressure_boundary", {"Pt": 2e5, "Tt": 300}, 0),
            node(
                "ori",
                "orifice",
                {"diameter": 0.01, "correlation": correlation, "label": "meter"},
                100,
            ),
            node("pb_out", "pressure_boundary", {"Pt": 1e5, "Tt": 300}, 200),
        ],
        "edges": [edge("e1", "pb_in", "ori"), edge("e2", "ori", "pb_out")],
    }


def test_a_saved_gui_network_with_stolz_loads_solves_and_says_so() -> None:
    from gui.backend.main import solve
    from gui.backend.schemas import NetworkGraphSchema, OrificeData

    assert OrificeData(correlation="Stolz").correlation == "Stolz"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        res = asyncio.run(solve(NetworkGraphSchema(**_gui_graph("Stolz"))))
    assert res.success, res.message
    assert len(res.deprecations) == 1
    note = res.deprecations[0]
    assert "meter" in note and "Stolz" in note and "0.9.0" in note
    assert "Reader-Harris/Gallagher" in note


def test_a_gui_network_without_stolz_has_no_deprecation_note() -> None:
    from gui.backend.messages import deprecations
    from gui.backend.schemas import NetworkGraphSchema

    schema = NetworkGraphSchema(**_gui_graph("ReaderHarrisGallagher"))
    assert deprecations(schema) == []
