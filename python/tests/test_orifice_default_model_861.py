"""The compressible orifice model an unset ``compressible_model`` selects (#861).

Selected ONCE at build from the Cd correlation and the geometry, never per
iterate: 'jobson1955' for a Cd referenced to incompressible flow
(Reader-Harris/Gallagher, the Idelchik family, Lichtarowicz,
McGreehan-Schotsch's baseline) on a sharp-edged bore at beta <= 0.75 in both
directions; the nozzle for 'fixed' (the user's measurement basis is unknown),
Stolz (ISO 5167:1980 not held), a rounded or beveled edge, and a wider bore.
An explicit choice always wins. The cooling-plate and connection subclasses
keep the nozzle (McGreehan-Schotsch's crossflow term is Rohde-fitted and
nozzle-referenced).
"""

from __future__ import annotations

import warnings

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    DiameterDischargeCoefficientConnectionElement,
    EffectiveAreaConnectionElement,
    EffusionPlateElement,
    FilmCoolingElement,
    FlowNetwork,
    ImpingementPlateElement,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    SlotFilmElement,
)
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.resolve_preview import resolve_preview
from gui.backend.schemas import NetworkGraphSchema

# Stolz is deprecated (#876): its DeprecationWarning is tested in
# test_stolz_deprecation_876.py; here Stolz is exercised for its values.
pytestmark = pytest.mark.filterwarnings("ignore:.*Stolz.*:DeprecationWarning")

INCOMPRESSIBLE_REFERENCED = {
    "ReaderHarrisGallagher": {},
    "IdelchikThick": {"plate_thickness": 0.01},
    "IdelchikBeveled": {"plate_thickness": 0.002},
    "IdelchikRounded": {"plate_thickness": 0.002},
    "Lichtarowicz": {"plate_thickness": 0.04},
    "McGreehanSchotsch": {"plate_thickness": 0.02},
}


def _net(
    bore: float = 0.02,
    D_up: float | None = None,
    D_down: float | None = None,
    **kw: object,
) -> FlowNetwork:
    """in -> [pipe D_up] -> a -> orifice -> b -> [pipe D_down] -> out; no
    pipe where None (a plenum on that side)."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0))
    n.add_node(PressureBoundary("out", Pt=1.2e5, Tt=300.0))
    n.add_node(PlenumNode("a"))
    n.add_node(PlenumNode("b"))
    for cid, s, t, D in (("c1", "in", "a", D_up), ("c2", "b", "out", D_down)):
        if D is None:
            n.add_element(OrificeElement(cid, s, t, diameter=0.2, correlation="fixed", Cd=1.0))
        else:
            n.add_element(ChannelElement(cid, s, t, length=0.5, diameter=D))
    if kw.get("correlation") == "fixed":
        kw.setdefault("Cd", 0.6)  # 'fixed' needs its Cd (#872)
    n.add_element(OrificeElement("x", "a", "b", diameter=bore, **kw))
    return n


def _model(net: FlowNetwork) -> str:
    net.resolve_all_topology()
    return net.elements["x"].compressible_model


@pytest.mark.parametrize("corr", sorted(INCOMPRESSIBLE_REFERENCED))
def test_an_incompressible_referenced_cd_selects_jobson(corr: str) -> None:
    kw = INCOMPRESSIBLE_REFERENCED[corr]
    o = OrificeElement("x", "a", "b", diameter=0.02, correlation=corr, **kw)
    assert o.compressible_model == "jobson1955"  # before build: the beta = 0 pick
    assert _model(_net(correlation=corr, **kw)) == "jobson1955"  # plenum-fed
    assert _model(_net(D_up=0.04, D_down=0.04, correlation=corr, **kw)) == "jobson1955"


@pytest.mark.parametrize("corr", ["fixed", "Stolz"])
def test_fixed_and_stolz_keep_the_nozzle(corr: str) -> None:
    assert _model(_net(correlation=corr)) == "nozzle"
    assert _model(_net(D_up=0.04, D_down=0.04, correlation=corr)) == "nozzle"


@pytest.mark.parametrize(
    ("corr", "edge"),
    [
        ("ReaderHarrisGallagher", {"edge_radius": 0.001}),
        ("IdelchikRounded", {"edge_radius": 0.001, "plate_thickness": 0.002}),
        ("IdelchikBeveled", {"bevel_depth": 0.001, "plate_thickness": 0.002}),
        ("McGreehanSchotsch", {"edge_radius": 0.001, "plate_thickness": 0.02}),
    ],
)
def test_a_rounded_or_beveled_edge_keeps_the_nozzle(corr: str, edge: dict) -> None:
    assert _model(_net(correlation=corr, **edge)) == "nozzle"


def test_beta_above_jobsons_range_keeps_the_nozzle_in_either_direction() -> None:
    assert cb.JOBSON_1955_BETA_MAX == 0.75
    assert _model(_net(0.0749, D_up=0.1, D_down=0.1)) == "jobson1955"
    assert _model(_net(0.0751, D_up=0.1)) == "nozzle"  # forward
    assert _model(_net(0.0751, D_down=0.1)) == "nozzle"  # reverse approach
    # A lumped restriction (in_pipe=False) has beta = 0 whatever adjoins it.
    assert _model(_net(0.0751, D_up=0.1, D_down=0.1, in_pipe=False)) == "jobson1955"


@pytest.mark.parametrize("model", ["nozzle", "jobson1955", "iso5167"])
def test_an_explicit_choice_always_wins(model: str) -> None:
    assert _model(_net(compressible_model=model)) == model
    assert _model(_net(0.03, D_up=0.03 / 0.9, compressible_model=model)) == model
    if model != "iso5167":
        assert _model(_net(correlation="fixed", compressible_model=model)) == model


def test_cv_needs_jobson_set_explicitly() -> None:
    with pytest.raises(ValueError, match="explicitly"):
        OrificeElement("x", "a", "b", diameter=0.01, Cv=0.98)
    o = OrificeElement("x", "a", "b", diameter=0.01, compressible_model="jobson1955", Cv=0.98)
    assert o.Cv == 0.98


@pytest.mark.parametrize(
    "cls",
    [
        EffusionPlateElement,
        FilmCoolingElement,
        SlotFilmElement,
        ImpingementPlateElement,
        EffectiveAreaConnectionElement,
        DiameterDischargeCoefficientConnectionElement,
    ],
)
def test_the_cooling_and_connection_subclasses_keep_the_nozzle(cls: type) -> None:
    assert cls._selects_compressible_model is False


def test_connection_elements_are_nozzle_instances() -> None:
    e = EffectiveAreaConnectionElement("x", "a", "b", effective_area=1e-4)
    assert e.compressible_model == "nozzle"
    e = DiameterDischargeCoefficientConnectionElement("x", "a", "b", diameter=0.01, Cd=0.6)
    assert e.compressible_model == "nozzle"


def test_the_model_is_fixed_through_a_solve_above_ci_0_7() -> None:
    """A long hole (Lichtarowicz, l/d 2) solves at Ci ~0.8, above Jobson's
    0.7: the model is not switched on Ci, per iterate or at all. It stays
    jobson1955 and model_in_range flags it."""
    net = _net(correlation="Lichtarowicz", plate_thickness=0.04)
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    assert net.elements["x"].compressible_model == "jobson1955"
    assert r["x.compressible_model"] == "jobson1955"
    assert r["x.Ci"] > 0.7 and r["x.model_in_range"] == 0.0


# --- GUI --------------------------------------------------------------------


def _gui(**orf: object) -> dict:
    def node(i: str, t: str, **data: object) -> dict:
        return {"id": i, "type": t, "position": {"x": 0.0, "y": 0.0}, "data": data}

    def edge(i: str, s: str, t: str) -> dict:
        return {"id": i, "source": s, "target": t, "data": {}}

    return {
        "nodes": [
            node("in", "pressure_boundary", Pt=2.0e5, Tt=300.0),
            node("c1", "channel", L=1.0, D=0.1),
            node("orf", "orifice", diameter=0.05, **orf),
            node("c2", "channel", L=1.0, D=0.1),
            node("out", "pressure_boundary", Pt=1.0e5, Tt=300.0),
        ],
        "edges": [
            edge("e1", "in", "c1"),
            edge("e2", "c1", "orf"),
            edge("e3", "orf", "c2"),
            edge("e4", "c2", "out"),
        ],
    }


@pytest.mark.parametrize(
    ("data", "model", "source"),
    [
        ({}, "jobson1955", "selected at build"),  # a file saved without the key
        ({"compressible_model": "default"}, "jobson1955", "selected at build"),
        ({"compressible_model": "nozzle"}, "nozzle", "user"),  # an explicit pick is kept
        ({"compressible_model": "jobson1955"}, "jobson1955", "user"),
        ({"correlation": "fixed", "Cd": 0.6}, "nozzle", "selected at build"),
    ],
)
def test_the_gui_default_selects_at_build_and_keeps_an_explicit_model(
    data: dict, model: str, source: str
) -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_gui(**data)))
    rep = resolve_preview(net)
    assert net.elements["orf"].compressible_model == model
    f = rep["orf"].fields["jobson1955"]
    assert f.value == (1.0 if model == "jobson1955" else 0.0) and f.source == source
