"""Geometry and input checks (#1024).

A zero-area restriction is CLOSED (a wall, row m = 0); zero or negative area on
a duct or node is refused; holes must leave material between them; a hydraulic
diameter must satisfy the isoperimetric bound; a turbomachine's pressure ratio
must be finite; a combustor's method must be one it has; an energy boundary
may not remove all the sensible enthalpy or more.

Each check has its own test, and each test fails with its check removed.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.components import (
    CombustorNode,
    CompressorElement,
    DiameterDischargeCoefficientConnectionElement,
    EffectiveAreaConnectionElement,
    EffusionPlateElement,
    EnergyBoundary,
    FilmCoolingElement,
    ImpingementPlateElement,
    MomentumChamberNode,
    NetworkMixtureState,
    OrificeElement,
    PressureLossElement,
    TurbineElement,
    WallNode,
)
from combaero.network.pressure_loss import ConstantFractionLoss
from gui.backend.resolve_preview import resolve_preview

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _closed(kind: str, eid: str = "v", a: str = "a", b: str = "b") -> OrificeElement:
    if kind == "orifice":
        return OrificeElement(eid, a, b, diameter=0.0, correlation="fixed", Cd=0.6)
    if kind == "orifice_rhg":
        return OrificeElement(eid, a, b, area=0.0)
    if kind == "effective_area":
        return EffectiveAreaConnectionElement(eid, a, b, 0.0)
    return DiameterDischargeCoefficientConnectionElement(eid, a, b, diameter=0.0, Cd=0.6)


def _line(mid: object, mid_nodes: tuple[str, ...] = ("a", "b"), nodes: dict | None = None):
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    for nid in mid_nodes:
        n.add_node((nodes or {}).get(nid) or PlenumNode(nid))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    n.add_element(ChannelElement("c1", "in", mid_nodes[0], length=1.0, diameter=0.05))
    n.add_element(mid)
    n.add_element(ChannelElement("c2", mid_nodes[-1], "out", length=1.0, diameter=0.05))
    return n


def _built(n: FlowNetwork) -> FlowNetwork:
    n.resolve_all_topology()
    n.validate()
    return n


# --- 1. A zero-area restriction is closed ------------------------------------


@pytest.mark.parametrize("kind", ["orifice", "orifice_rhg", "effective_area", "dcd"])
def test_zero_area_restriction_is_a_wall(kind: str) -> None:
    """The row is m = 0 with the constant Jacobian 1 in its own flow, at a
    state where the kernels would refuse A = 0 ("nozzle_flow: A_eff must be
    positive" before #1024)."""
    el = _closed(kind)
    assert el.closed
    n = _built(_line(el))
    s_in = NetworkMixtureState(2e5, 2e5, Tt=300.0, m_dot=0.37, Y=Y_AIR)
    s_out = NetworkMixtureState(1e5, 1e5, Tt=300.0, m_dot=0.37, Y=Y_AIR)
    res, jac = n.elements["v"].residuals(s_in, s_out)
    assert res == [0.37]
    assert jac == {0: {"v.m_dot": 1.0}}


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize("kind", ["orifice", "effective_area", "dcd"])
def test_closed_restriction_in_series_solves_to_zero_flow(kind: str, regime: str) -> None:
    n = _line(_closed(kind))
    for e in n.elements.values():
        e.regime = regime
    r = NetworkSolver(n).solve()
    assert r["__success__"], r["__message__"]
    # On the first attempt: the start seeds the closed flow at zero, not a
    # share of the network's (it took the warm-start retry before).
    assert r["__message__"] == "The solution converged."
    for eid in ("c1", "v", "c2"):
        assert abs(r[f"{eid}.m_dot"]) < 1e-12
    # Each side sits at its own boundary's pressure.
    assert r["a.P"] == pytest.approx(2.0e5, rel=1e-9)
    assert r["b.P"] == pytest.approx(1.0e5, rel=1e-9)
    diag = r["__element_diag__"]["v"]
    assert diag["closed"] == 1.0


def test_the_propagated_start_gives_a_closed_restriction_no_flow() -> None:
    """The fallback start (networks the resistance seed does not cover)
    gives a closed element zero flow, as it does a duct into a WallNode,
    not an equal share of its node's."""
    n = _built(_line(_closed("orifice")))
    s = NetworkSolver(n)
    guess = s._propagate_mdot_guess(s._infer_reference_state())
    assert guess["v"] == 0.0
    assert guess["c1"] != 0.0


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_closed_bypass_leaves_the_open_path_alone(regime: str) -> None:
    """A closed valve in parallel with an open orifice carries nothing; the
    open path's flow is that of the network without the valve's branch."""

    def net(with_closed: bool) -> FlowNetwork:
        n = FlowNetwork()
        n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
        n.add_node(PlenumNode("b"))
        n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
        n.add_element(OrificeElement("o", "in", "b", diameter=0.02, correlation="fixed", Cd=0.6))
        n.add_element(ChannelElement("c2", "b", "out", length=1.0, diameter=0.05))
        if with_closed:
            n.add_node(PlenumNode("a"))
            n.add_element(ChannelElement("c1", "in", "a", length=1.0, diameter=0.05))
            n.add_element(_closed("orifice"))
        for e in n.elements.values():
            e.regime = regime
        return n

    ref = NetworkSolver(net(False)).solve()
    r = NetworkSolver(net(True)).solve()
    assert r["__success__"], r["__message__"]
    assert abs(r["v.m_dot"]) < 1e-12
    assert r["o.m_dot"] == pytest.approx(ref["o.m_dot"], rel=1e-8)


def test_part_cut_off_by_a_closed_restriction_is_refused() -> None:
    """Beyond the closed valve a dead branch with no PressureBoundary: its
    pressure is set by nothing. Refused at validate(), naming the valve."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    n.add_node(PlenumNode("a"))
    n.add_node(PlenumNode("b"))
    n.add_node(WallNode("w"))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    n.add_element(ChannelElement("c1", "in", "a", length=1.0, diameter=0.05))
    n.add_element(ChannelElement("c2", "a", "out", length=1.0, diameter=0.05))
    n.add_element(_closed("orifice"))
    n.add_element(ChannelElement("c3", "b", "w", length=1.0, diameter=0.05))
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=r"\['b', 'w'\].*only through closed restriction.*\['v'\]"):
        n.validate()


def test_fraction_loss_behind_a_closed_valve_is_a_dead_end() -> None:
    """The #1004 drive rule: a fraction loss without m_ref at zero flow in a
    driven network is flagged, unless it is a dead end. Behind a closed valve
    it is one: no flag."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    n.add_node(PlenumNode("a"))
    n.add_node(PlenumNode("b"))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    n.add_element(PressureLossElement("f", "in", "a", correlation=ConstantFractionLoss(xi=0.05)))
    n.add_element(_closed("orifice"))
    n.add_element(ChannelElement("c2", "b", "out", length=1.0, diameter=0.05))
    r = NetworkSolver(n).solve()
    assert r["__success__"], r["__message__"]
    assert r["__fraction_unanchored__"] == []
    assert r["f.fraction_unanchored"] == 0.0


@pytest.mark.parametrize(
    "build",
    [
        lambda: OrificeElement("o", "a", "b", diameter=-0.01, correlation="fixed", Cd=0.6),
        lambda: OrificeElement("o", "a", "b", area=-1e-4, correlation="fixed", Cd=0.6),
        lambda: OrificeElement("o", "a", "b", diameter=math.nan, correlation="fixed", Cd=0.6),
        lambda: EffectiveAreaConnectionElement("o", "a", "b", -1e-4),
        lambda: DiameterDischargeCoefficientConnectionElement(
            "o", "a", "b", diameter=-0.01, Cd=0.6
        ),
    ],
)
def test_negative_restriction_area_is_refused(build) -> None:
    with pytest.raises(ValueError, match="non-negative"):
        build()


def test_resolve_preview_reports_a_closed_restriction() -> None:
    rep = resolve_preview(_line(_closed("orifice")))
    assert rep["v"].fields["closed"].value == 1.0


# --- 2. Ducts and nodes: zero or negative area refused -----------------------


@pytest.mark.parametrize("D", [0.0, -0.05])
def test_channel_without_area_is_refused(D: float) -> None:
    n = _line(ChannelElement("m", "a", "b", length=1.0, diameter=D))
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=r"ChannelElement 'm': diameter D = .* positive"):
        n.validate()


@pytest.mark.parametrize("area", [0.0, -1.0])
def test_chamber_without_area_is_refused(area: float) -> None:
    n = _line(
        ChannelElement("m", "a", "b", length=1.0, diameter=0.05),
        nodes={"a": MomentumChamberNode("a", area=area)},
    )
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=r"MomentumChamberNode 'a': area = .* positive"):
        n.validate()


def test_chamber_negative_Dh_is_refused() -> None:
    n = _line(
        ChannelElement("m", "a", "b", length=1.0, diameter=0.05),
        nodes={"a": MomentumChamberNode("a", Dh=-0.1)},
    )
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=r"MomentumChamberNode 'a': hydraulic diameter Dh"):
        n.validate()


@pytest.mark.parametrize(
    ("kw", "what"),
    [({"area": 0.0}, "area"), ({"area": -1.0}, "area"), ({"Dh": -0.1}, "hydraulic diameter")],
)
def test_combustor_without_area_or_Dh_is_refused(kw: dict, what: str) -> None:
    n = _line(
        ChannelElement("m", "a", "b", length=1.0, diameter=0.05),
        nodes={"a": CombustorNode("a", **kw)},
    )
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=rf"CombustorNode 'a': {what}"):
        n.validate()


# --- 3. Hydraulic diameter: Dh <= 2 sqrt(A/pi) --------------------------------


def test_channel_Dh_must_be_positive_but_may_exceed_its_diameter() -> None:
    """A channel standing for a share of a wider duct keeps that duct's Dh,
    above its own area's 2 sqrt(A/pi) (the C3X example's hole segments):
    not bounded. A Dh of zero is refused."""
    n = _line(ChannelElement("m", "a", "b", length=1.0, diameter=0.05, Dh=0.0))
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=r"ChannelElement 'm': hydraulic diameter Dh = 0.0 m"):
        n.validate()
    _built(_line(ChannelElement("m", "a", "b", length=1.0, diameter=0.05, Dh=0.08)))


def test_chamber_Dh_above_isoperimetric_bound_is_refused() -> None:
    n = _line(
        ChannelElement("m", "a", "b", length=1.0, diameter=0.05),
        nodes={"a": MomentumChamberNode("a", area=1e-3, Dh=0.5)},
    )
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=r"MomentumChamberNode 'a': hydraulic diameter Dh = 0.5"):
        n.validate()
    D = 2.0 * math.sqrt(1e-3 / math.pi)
    _built(
        _line(
            ChannelElement("m", "a", "b", length=1.0, diameter=0.05),
            nodes={"a": MomentumChamberNode("a", area=1e-3, Dh=D)},
        )
    )


def test_combustor_Dh_above_isoperimetric_bound_is_refused() -> None:
    n = _line(
        ChannelElement("m", "a", "b", length=1.0, diameter=0.05),
        nodes={"a": CombustorNode("a", area=1e-3, Dh=0.5)},
    )
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=r"CombustorNode 'a': hydraulic diameter Dh = 0.5"):
        n.validate()


# --- 4. Turbomachines: a finite pressure ratio -------------------------------


@pytest.mark.parametrize("cls", [CompressorElement, TurbineElement])
def test_infinite_pressure_ratio_is_refused(cls: type) -> None:
    with pytest.raises(ValueError, match="must be a finite number"):
        cls("t", "a", "b", pressure_ratio=math.inf, eta_is=0.9)
    cls("t", "a", "b", pressure_ratio=1.0, eta_is=0.9)  # a pass-through stays valid


# --- 5. Combustor method ------------------------------------------------------


def test_combustor_method_typo_is_refused() -> None:
    with pytest.raises(ValueError, match=r"method 'equilibirum' is not one of"):
        CombustorNode("a", method="equilibirum")
    node = CombustorNode("a")
    node.method = "nonsense"  # set after construction: validate() catches it
    n = _line(ChannelElement("m", "a", "b", length=1.0, diameter=0.05), nodes={"a": node})
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=r"method 'nonsense' is not one of"):
        n.validate()


# --- 6. Energy boundary: fraction > -1 ----------------------------------------


@pytest.mark.parametrize("fraction", [-1.0, -1.5])
def test_fraction_removing_all_sensible_enthalpy_is_refused(fraction: float) -> None:
    with pytest.raises(ValueError, match="sensible enthalpy"):
        EnergyBoundary("q", fraction=fraction)
    EnergyBoundary("q", fraction=-0.99)


def test_fractions_summing_to_minus_one_on_a_node_are_refused() -> None:
    node = PlenumNode("a")
    node.add_energy_boundary(EnergyBoundary("q1", fraction=-0.6))
    node.add_energy_boundary(EnergyBoundary("q2", fraction=-0.5))
    n = _line(ChannelElement("m", "a", "b", length=1.0, diameter=0.05), nodes={"a": node})
    n.resolve_all_topology()
    with pytest.raises(ValueError, match=r"Node 'a' \('q1', 'q2'\): fraction -1.1"):
        n.validate()


# --- 7. Holes leave material between them -------------------------------------

D_HOLE = 1e-3


def _effusion(**kw) -> EffusionPlateElement:
    args = {"hole_diameter": D_HOLE, "wall_thickness": 2e-3, "panel_area": 1e-2}
    args.update(kw)
    return EffusionPlateElement("e", "a", "b", **args)


def test_effusion_spanwise_pitch_must_exceed_d() -> None:
    with pytest.raises(ValueError, match=r"spanwise pitch_y = 0.001 m is not larger"):
        _effusion(pitch_x=5 * D_HOLE, pitch_y=D_HOLE).validate()  # equal: touching
    _effusion(pitch_x=5 * D_HOLE, pitch_y=1.01 * D_HOLE).validate()


def test_effusion_streamwise_pitch_must_exceed_the_footprint() -> None:
    # 30 deg: the footprint is d / sin(30) = 2 d; porosity 0.35 passes its own check.
    with pytest.raises(
        ValueError, match=r"streamwise pitch_x = 0.0015 m .* d / sin\(angle\) = 0.002"
    ):
        _effusion(pitch=1.5 * D_HOLE, angle_deg=30.0).validate()
    _effusion(pitch_x=2.01 * D_HOLE, pitch_y=1.5 * D_HOLE, angle_deg=30.0).validate()


def _film(**kw) -> FilmCoolingElement:
    args = {
        "hole_diameter": D_HOLE,
        "holes_per_row": 10,
        "lateral_pitch": 3 * D_HOLE,
        "wall_thickness": 2e-3,
    }
    args.update(kw)
    return FilmCoolingElement("f", "a", "b", **args)


def test_film_lateral_pitch_must_exceed_d() -> None:
    with pytest.raises(ValueError, match=r"spanwise lateral_pitch = 0.0005 m"):
        _film(lateral_pitch=0.5 * D_HOLE).validate()


def test_film_row_pitch_must_exceed_the_footprint() -> None:
    # 30 deg default: footprint 2 d.
    with pytest.raises(ValueError, match=r"streamwise row_pitch = 0.002 m"):
        _film(n_rows=3, row_pitch=2 * D_HOLE).validate()
    _film(n_rows=3, row_pitch=2.01 * D_HOLE).validate()


def _jets(**kw) -> ImpingementPlateElement:
    args = {"d_jet": D_HOLE, "xn_d": 5.0, "yn_d": 4.0, "z_d": 2.0, "span": 0.05}
    args.update(kw)
    return ImpingementPlateElement("i", "a", "b", plate_thickness=1e-3, **args)


def test_impingement_open_area_must_be_below_the_plate() -> None:
    with pytest.raises(ValueError, match=r"jet holes' area is 3.14 times the plate's"):
        _jets(xn_d=0.5, yn_d=0.5).validate()


def test_impingement_spanwise_pitch_must_exceed_d() -> None:
    with pytest.raises(ValueError, match=r"spanwise yn_d \* d_jet = 0.001 m"):
        _jets(yn_d=1.0).validate()


def test_impingement_inline_streamwise_pitch_must_exceed_d() -> None:
    with pytest.raises(ValueError, match=r"streamwise xn_d \* d_jet = 0.0008 m"):
        _jets(xn_d=0.8, yn_d=3.0).validate()


def test_impingement_staggered_uses_the_nearest_neighbour() -> None:
    # Rows 0.8 d apart clear each other when the next row is half a pitch over.
    _jets(xn_d=0.8, yn_d=3.0, pattern="staggered").validate()
    # Diagonal neighbour too close: hypot(0.6, 0.75) = 0.96 d.
    with pytest.raises(ValueError, match=r"staggered .* = 0.9605 d"):
        _jets(xn_d=0.6, yn_d=1.5, pattern="staggered").validate()
    # Same column two rows on: 2 xn = 0.9 d.
    with pytest.raises(ValueError, match=r"staggered .* = 0.9 d"):
        _jets(xn_d=0.45, yn_d=4.0, pattern="staggered").validate()


def test_resolve_preview_reports_a_geometry_refusal() -> None:
    n = _line(
        ChannelElement("m", "a", "b", length=1.0, diameter=0.05),
        nodes={"a": MomentumChamberNode("a", area=1e-3, Dh=0.5)},
    )
    rep = resolve_preview(n)
    assert rep["a"].status == "resolved"
    assert "2 sqrt(A/pi)" in rep["a"].refusal
    assert rep["c1"].refusal == ""
