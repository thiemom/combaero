"""An element discharging into a PressureBoundary reports the gas arriving.

A PressureBoundary's state is its set Tt/Y. No residual reads it for an
outflow (the solved pressures and flows do not move with the setting), but the
element's exit diagnostics did: a 600 K stream into an outlet set to 300 K
reported Tt_out = T_out = 300 K and a mach_out from 300 K gas.
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _chamber_to_outlet(T_outlet: float, regime: str) -> dict:
    net = FlowNetwork()
    m_dot = 1.2 if regime == "compressible" else 0.3
    net.add_node(MassFlowBoundary("inl", m_dot=m_dot, Tt=600.0, Y=Y_AIR))
    net.add_node(MomentumChamberNode("ch", area=0.003))
    net.add_node(PressureBoundary("out", Pt=2.0e5, Tt=T_outlet, Y=Y_AIR))
    net.add_element(ChannelElement("c1", "inl", "ch", length=0.2, diameter=0.0618, regime=regime))
    net.add_element(ChannelElement("c2", "ch", "out", length=0.2, diameter=0.0618, regime=regime))
    sol = NetworkSolver(net).solve()
    assert sol["__success__"]
    return sol


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
def test_exit_diagnostics_do_not_follow_the_outlet_setting(regime: str) -> None:
    cold = _chamber_to_outlet(300.0, regime)
    hot = _chamber_to_outlet(900.0, regime)
    # The solve itself never depended on the outlet's Tt.
    assert cold["ch.P"] == pytest.approx(hot["ch.P"], rel=1e-12)
    for sol in (cold, hot):
        assert sol["c2.Tt_out"] == pytest.approx(600.0, rel=1e-9)
        assert sol["c2.T_out"] == pytest.approx(600.0, rel=1e-9)
    assert cold["c2.mach_out"] == pytest.approx(hot["c2.mach_out"], rel=1e-9)
    # The boundary's own result stays its setting.
    assert cold["out.Tt"] == 300.0


def test_wall_heat_reaches_the_reported_exit() -> None:
    """A wall-coupled pair, each side discharging into an outlet: the exit Tt
    carries the element's own wall heat, Q = m cp dTt on both sides."""
    for T_set in (300.0, 1000.0):
        net = FlowNetwork()
        net.add_node(MassFlowBoundary("hin", m_dot=0.2, Tt=800.0, Y=Y_AIR))
        net.add_node(MassFlowBoundary("cin", m_dot=0.2, Tt=300.0, Y=Y_AIR))
        net.add_node(PressureBoundary("hout", Pt=101325.0, Tt=T_set, Y=Y_AIR))
        net.add_node(PressureBoundary("cout", Pt=101325.0, Tt=T_set, Y=Y_AIR))
        net.add_element(ChannelElement("h", "hin", "hout", length=2.0, diameter=0.05))
        net.add_element(ChannelElement("c", "cin", "cout", length=2.0, diameter=0.05))
        net.add_wall(ThermalWall("w", "h", "c", [WallLayer(thickness=0.002, conductivity=25.0)]))
        sol = NetworkSolver(net).solve()
        assert sol["__success__"]
        Q = sol["w.Q"]
        assert Q > 1.0e3
        X = list(cb.species.dry_air())
        h_drop = cb.h_mass(800.0, X) - cb.h_mass(sol["h.Tt_out"], X)
        c_rise = cb.h_mass(sol["c.Tt_out"], X) - cb.h_mass(300.0, X)
        assert 0.2 * h_drop == pytest.approx(Q, rel=1e-3)
        assert 0.2 * c_rise == pytest.approx(Q, rel=1e-3)


def test_a_reversed_element_reports_the_gas_it_delivers_to_its_from_boundary() -> None:
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("src", m_dot=0.3, Tt=600.0, Y=Y_AIR))
    net.add_node(PlenumNode("m"))
    net.add_node(PressureBoundary("a", Pt=101325.0, Tt=300.0, Y=Y_AIR))
    net.add_node(PressureBoundary("b", Pt=101325.0, Tt=300.0, Y=Y_AIR))
    net.add_element(ChannelElement("feed", "src", "m", length=0.5, diameter=0.05))
    net.add_element(ChannelElement("e", "a", "m", length=1.0, diameter=0.05))
    net.add_element(ChannelElement("f", "m", "b", length=1.0, diameter=0.05))
    sol = NetworkSolver(net).solve()
    assert sol["__success__"]
    assert sol["e.m_dot"] < 0.0
    assert sol["e.Tt_in"] == pytest.approx(600.0, rel=1e-9)
    assert sol["f.Tt_out"] == pytest.approx(600.0, rel=1e-9)


def test_a_pressure_boundary_supplying_flow_keeps_its_setting() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("a", Pt=1.2e5, Tt=450.0, Y=Y_AIR))
    net.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    net.add_element(ChannelElement("e", "a", "b", length=1.0, diameter=0.05))
    sol = NetworkSolver(net).solve()
    assert sol["__success__"]
    assert sol["e.Tt_in"] == pytest.approx(450.0, rel=1e-12)
    assert sol["e.Tt_out"] == pytest.approx(450.0, rel=1e-9)
