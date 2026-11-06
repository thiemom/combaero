"""An effusion plate with no coolant flow (#998): flagged, never silent.

* Hot gas, no coolant (a chamber with gas flow): the adiabatic limit, both
  faces at the gas T_aw, and ``wall_above_limit`` when that is above the
  plate material's #882 envelope hot end (``T_wall_limit``).
* No coolant and no gas-side flow (a plenum, a stagnant chamber): no cooling
  path, wall T undetermined: ``no_cooling_path`` = 1 and the solve refused
  with a message naming the plate.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EffusionLiner,
    EffusionPlateElement,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
    SolveOutcome,
)

Y_AIR = cb.mole_to_mass(cb.species.dry_air())
MESSAGE = "no cooling path: wall T undetermined (heats until failure)"


def _pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _plate(to_node: str = "ch", **kw) -> EffusionPlateElement:
    args = {
        "hole_diameter": 0.8e-3,
        "wall_thickness": 2e-3,
        "pitch": 6e-3,
        "panel_length": 0.1,
        "panel_width": 0.1,
        "angle_deg": 30.0,
    }
    args.update(kw)
    return EffusionPlateElement("eff", "cool", to_node, **args)


def _chamber(P_cool: float, T_gas: float, P_gas_in: float = 1.003e5, **kw) -> dict:
    """The plate into a momentum chamber on a gas duct; P_cool equal to the
    chamber's pressure gives no coolant flow, P_gas_in = 1e5 no gas flow."""
    g = FlowNetwork()
    g.add_node(_pb("gas_in", P_gas_in, T_gas))
    g.add_node(_pb("exit", 1.0e5, T_gas))
    g.add_node(_pb("cool", P_cool, 600.0))
    g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
    g.add_element(ChannelElement("duct", "gas_in", "ch", length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(ChannelElement("tail", "ch", "exit", length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(_plate(**kw))
    return NetworkSolver(g).solve()


def _no_coolant_in_hot_gas(T_gas: float, **kw) -> tuple[dict, dict]:
    """The plate in flowing gas, fed by a mass-flow boundary at zero flow."""
    g = FlowNetwork()
    g.add_node(_pb("gas_in", 1.003e5, T_gas))
    g.add_node(_pb("exit", 1.0e5, T_gas))
    g.add_node(MassFlowBoundary("cool", m_dot=0.0, Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
    g.add_element(ChannelElement("duct", "gas_in", "ch", length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(ChannelElement("tail", "ch", "exit", length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(_plate(**kw))
    r = NetworkSolver(g).solve()
    return r, r["__element_diag__"]["eff"]


def test_the_limit_is_the_materials_envelope_hot_end() -> None:
    info = cb.get_material_info("haynes230")
    assert _plate(material="haynes230").T_wall_limit == info["T_envelope_max"]
    assert _plate().T_wall_limit == cb.GENERIC_METAL_SOLIDUS == 1533.15
    with pytest.raises(ValueError, match="unknown material"):
        _plate(material="unobtainium")


def test_hot_gas_with_no_coolant_is_adiabatic_and_flagged_above_the_limit() -> None:
    r, d = _no_coolant_in_hot_gas(1800.0)
    assert r["__success__"], r.get("__message__")
    # A zero-flow root: |m_dot| ~ 1e-35, below the no-through-flow floor.
    assert abs(d["m_dot"]) < cb.MIXER_HEAT_MDOT_FLOOR
    assert d["Re_hole"] == 0.0
    assert d["no_cooling_path"] == 0.0
    # The adiabatic limit: no heat taken, both faces at the gas T_aw.
    assert d["q_wall"] == 0.0
    assert d["T_wall_hot"] == d["T_wall_cold"] == d["T_adiabatic_wall"]
    assert d["T_wall_hot"] > cb.GENERIC_METAL_SOLIDUS
    assert d["T_wall_limit"] == cb.GENERIC_METAL_SOLIDUS
    assert d["wall_above_limit"] == 1.0
    assert r["__wall_above_limit__"] == ["eff"]


def test_hot_gas_with_no_coolant_below_the_limit_is_not_flagged() -> None:
    r, d = _no_coolant_in_hot_gas(1000.0)
    assert r["__success__"], r.get("__message__")
    assert d["q_wall"] == 0.0
    assert d["T_wall_hot"] < cb.GENERIC_METAL_SOLIDUS
    assert d["wall_above_limit"] == 0.0
    assert r["__wall_above_limit__"] == []


def test_the_material_moves_the_limit() -> None:
    """haynes230's solidus (1574.15 K) is above the generic cap (1533.15 K):
    a wall between them flags without a material and not with it."""
    T_gas = 1555.0
    _, generic = _no_coolant_in_hot_gas(T_gas)
    _, h230 = _no_coolant_in_hot_gas(T_gas, material="haynes230")
    assert cb.GENERIC_METAL_SOLIDUS < generic["T_wall_hot"] < h230["T_wall_limit"]
    assert generic["wall_above_limit"] == 1.0
    assert h230["wall_above_limit"] == 0.0


def test_a_cooled_plate_below_the_limit_is_not_flagged() -> None:
    r = _chamber(1.03e5, 1500.0)
    assert r["__success__"], r.get("__message__")
    d = r["__element_diag__"]["eff"]
    assert d["no_cooling_path"] == 0.0
    assert d["wall_above_limit"] == 0.0
    assert r["__wall_undetermined__"] == []


def _plenum(P_cool: float, **kw) -> dict:
    g = FlowNetwork()
    g.add_node(_pb("cool", P_cool, 600.0))
    g.add_node(_pb("gas", 1.0e5, 1500.0))
    g.add_element(_plate("gas", **kw))
    return NetworkSolver(g).solve()


@pytest.mark.parametrize("gas_heat_flux", [0.0, 1.0e5])
def test_a_plenum_with_no_coolant_is_refused_naming_the_plate(gas_heat_flux: float) -> None:
    with pytest.warns(UserWarning, match="no cooling path"):
        r = _plenum(1.0e5, gas_heat_flux=gas_heat_flux)
    assert r["__converged__"]
    assert not r["__success__"]
    assert r["__outcome__"] == SolveOutcome.INCONSISTENT
    assert f"'eff': {MESSAGE}" in r["__message__"]
    assert r["__wall_undetermined__"] == ["eff"]
    d = r["__element_diag__"]["eff"]
    assert d["no_cooling_path"] == 1.0
    assert math.isnan(d["T_wall_hot"]) and math.isnan(d["wall_above_limit"])
    # The keys stay: the limit is the plate's, whatever the operating point.
    assert d["T_wall_limit"] == cb.GENERIC_METAL_SOLIDUS


def test_a_plenum_with_coolant_flowing_solves() -> None:
    r = _plenum(1.03e5)
    assert r["__success__"], r.get("__message__")
    assert r["__element_diag__"]["eff"]["no_cooling_path"] == 0.0
    assert r["__wall_undetermined__"] == []


def test_a_stagnant_chamber_with_no_coolant_is_refused() -> None:
    with pytest.warns(UserWarning, match="no cooling path"):
        r = _chamber(1.0e5, 1500.0, P_gas_in=1.0e5)
    assert not r["__success__"]
    assert MESSAGE in r["__message__"]
    assert r["__element_diag__"]["eff"]["no_cooling_path"] == 1.0
    assert r["__element_diag__"]["eff"]["gas_side_chamber"] == 1.0
    # Converged, so nothing may call it unconverged: the auto-retry's
    # "skipped: neither ... proxy solve converged" note was appended here.
    assert r["__converged__"]
    assert "converge" not in r["__message__"].replace("Converged, but", "")


def _liner_net(P_cin: float, P_cout: float, P_gas_in: float) -> tuple[FlowNetwork, EffusionLiner]:
    g = FlowNetwork()
    g.add_node(_pb("cin", P_cin, 600.0))
    g.add_node(_pb("cout", P_cout, 600.0))
    g.add_node(_pb("gas_in", P_gas_in, 1500.0))
    g.add_node(_pb("exit", 1.0e5, 1500.0))
    g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
    g.add_element(ChannelElement("duct", "gas_in", "ch", length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(ChannelElement("tail", "ch", "exit", length=0.2, diameter=0.08, roughness=0.0))
    liner = EffusionLiner(
        "ln",
        n_segments=2,
        length=0.152,
        width=0.152,
        duct_height=0.03,
        hole_diameter=3.27e-3,
        wall_thickness=6.3e-3,
        pitch_x=15.24e-3,
        pitch_y=15.24e-3,
        material="haynes230",
    )
    liner.add_to(g, "cin", "cout", "ch")
    return g, liner


def test_a_liner_reports_its_panels_flags() -> None:
    g, liner = _liner_net(1.03e5, 1.025e5, 1.003e5)
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    s = liner.summarize(r["__element_diag__"])
    assert s["no_cooling_path"] == 0.0
    assert s["wall_above_limit"] == 0.0
    assert s["T_wall_limit"] == cb.get_material_info("haynes230")["T_envelope_max"]


def test_a_liner_at_zero_flow_is_refused() -> None:
    """Every boundary at one pressure: the zero-flow root leaves ~1e-13 kg/s
    in the chamber, below the no-through-flow floor, so no panel has a
    cooling path (it was reported as an adiabatic 1500 K wall)."""
    g, liner = _liner_net(1.0e5, 1.0e5, 1.0e5)
    with pytest.warns(UserWarning, match="no cooling path"):
        r = NetworkSolver(g).solve()
    assert not r["__success__"]
    assert r["__wall_undetermined__"] == [f"ln__p{i}" for i in (1, 2)]
    assert "'ln__p1', 'ln__p2': " + MESSAGE in r["__message__"]
    assert liner.summarize(r["__element_diag__"])["no_cooling_path"] == 1.0
