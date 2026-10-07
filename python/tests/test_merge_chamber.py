"""Merge chamber: one main inlet, side streams, one outlet (#471).

A MomentumChamberNode with ``main_inlet`` declared closes several inflows by
the axial momentum balance over a constant-area control volume:

    P_face A + m_main u_main + sum(m_s u_jet cos theta) = P A + m_out u_out

The node's (P, Pt) is the chamber/outlet state; the main-inlet element sees
the face state. The case here is an effusion panel discharging into a
liner's hot-gas chamber at a realistic 3% liner pressure drop.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize._numdiff import approx_derivative

import combaero as cb
from combaero.network import (
    ChannelElement,
    EffusionPlateElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkMixtureState,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)

Y_AIR = cb.mole_to_mass(cb.species.dry_air())
X_AIR = cb.species.dry_air()


def _pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _liner(angle_deg: float = 30.0, correlation: str = "IdelchikThick", **chamber) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(_pb("gas_in", 1.003e5, 1200.0))
    g.add_node(_pb("exit", 1.0e5, 1200.0))
    g.add_node(_pb("cool", 1.03e5, 600.0))
    g.add_node(MomentumChamberNode("ch", main_inlet="duct", **chamber))
    g.add_element(ChannelElement("duct", "gas_in", "ch", length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(ChannelElement("tail", "ch", "exit", length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(
        EffusionPlateElement(
            "eff",
            "cool",
            "ch",
            hole_diameter=0.8e-3,
            wall_thickness=2e-3,
            pitch=6e-3,
            panel_area=0.01,
            angle_deg=angle_deg,
            correlation=correlation,
        )
    )
    return g


def _solve(net: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


@pytest.mark.parametrize("angle", [90.0, 30.0])
def test_the_solution_satisfies_the_control_volume_momentum_balance(angle: float) -> None:
    s, r = _solve(_liner(angle))
    d = r["__element_diag__"]
    ch = s.network.nodes["ch"]
    P, T = r["ch.P"], s._derived_states["ch"][0]
    rho = cb.density(T, P, X_AIR)
    A = ch.area
    m_main, m_out = d["duct"]["m_dot"], d["tail"]["m_dot"]
    P_face = d["duct"]["P_out"]  # the main inlet reports the face it sees
    eff = s.network.elements["eff"]
    supply = NetworkMixtureState(
        T=600.0, P=1.03e5, Pt=1.03e5, Tt=600.0, m_dot=d["eff"]["m_dot"], Y=Y_AIR
    )
    chamber = NetworkMixtureState(T=T, P=P, Pt=r["ch.Pt"], Tt=T, m_dot=m_out, Y=Y_AIR)
    S, _ = eff.injection_momentum(supply, chamber)
    lhs = P_face * A + m_main**2 / (rho * A) + S
    rhs = P * A + m_out**2 / (rho * A)
    # Relative to the momentum fluxes, not to P A, which would hide an error.
    assert lhs - rhs == pytest.approx(0.0, abs=1e-6 * (m_out**2 / (rho * A) + abs(S)))
    assert m_out == pytest.approx(m_main + d["eff"]["m_dot"], rel=1e-9)
    if angle == 90.0:
        assert S == 0.0


def test_axial_injection_pumps_the_main_flow() -> None:
    """More axial jet momentum, more main flow: the jets work like an ejector."""
    flows = []
    for angle in (90.0, 60.0, 30.0, 20.0):
        _, r = _solve(_liner(angle))
        flows.append(r["__element_diag__"]["duct"]["m_dot"])
    assert all(b > a for a, b in zip(flows, flows[1:], strict=False))


def test_normal_injection_costs_the_main_stream_its_mixing_loss() -> None:
    """At 90 deg the main face sits above the chamber in both P and Pt: the
    side stream must be accelerated by the main flow."""
    _, r = _solve(_liner(90.0))
    d = r["__element_diag__"]
    assert d["duct"]["P_out"] > r["ch.P"]
    assert d["duct"]["Pt_out"] > r["ch.Pt"]


@pytest.mark.parametrize("correlation", ["fixed", "IdelchikThick"])
def test_global_jacobian_matches_finite_differences(correlation: str) -> None:
    """Exact, including the discharge-hole Cd's dependence on the flow, which
    the orifice residual and the jet momentum both carry now."""
    s, r = _solve(_liner(30.0, correlation=correlation))
    x = np.array(r["__x_solution__"], dtype=float)
    _, jac = s._residuals_and_jacobian(x)
    jac = jac.toarray()
    fd = approx_derivative(
        lambda v: s._residuals(v), x, method="3-point", abs_step=np.maximum(np.abs(x) * 1e-7, 1e-10)
    )
    mask = np.abs(fd) > 1e-3
    assert (np.abs(jac - fd)[mask] / np.abs(fd[mask])).max() < 1e-4


def test_the_area_is_inherited_from_the_main_inlet() -> None:
    s, _ = _solve(_liner(30.0))
    ch = s.network.nodes["ch"]
    assert ch.area == pytest.approx(s.network.elements["duct"].area)
    assert ch._area_source == "inherited from duct"


def test_a_mismatched_area_is_refused() -> None:
    with pytest.raises(ValueError, match="differs from its main inlet 'duct'"):
        NetworkSolver(_liner(30.0, area=0.01)).solve()


def test_two_inflows_without_a_main_inlet_are_refused() -> None:
    net = _liner(30.0)
    net.nodes["ch"].main_inlet = None
    with pytest.raises(ValueError, match="Declare main_inlet"):
        NetworkSolver(net).solve()


def test_a_main_inlet_that_does_not_flow_in_is_refused() -> None:
    net = _liner(30.0)
    net.nodes["ch"].main_inlet = "tail"
    with pytest.raises(ValueError, match="does not flow into it"):
        NetworkSolver(net).solve()


def test_a_split_is_still_refused() -> None:
    net = _liner(30.0)
    net.add_node(_pb("bleed", 0.99e5, 300.0))
    net.add_element(OrificeElement("bl", "ch", "bleed", diameter=0.01, correlation="fixed"))
    with pytest.raises(ValueError, match="multiple outgoing"):
        NetworkSolver(net).solve()


def test_a_single_inlet_chamber_is_unchanged_by_declaring_it_main() -> None:
    """No side streams: the face offset is zero, the old closure exactly."""

    def solo(main: str | None) -> dict:
        g = FlowNetwork()
        g.add_node(_pb("a", 1.01e5, 300.0))
        g.add_node(_pb("b", 1.0e5, 300.0))
        g.add_node(MomentumChamberNode("ch", main_inlet=main, area=np.pi * 0.04**2))
        g.add_element(ChannelElement("c1", "a", "ch", length=1.0, diameter=0.08))
        g.add_element(ChannelElement("c2", "ch", "b", length=1.0, diameter=0.08))
        return _solve(g)[1]

    a, b = solo(None), solo("c1")
    assert b["c1.m_dot"] == pytest.approx(a["c1.m_dot"], rel=1e-10)
    assert b["ch.P"] == pytest.approx(a["ch.P"], rel=1e-12)
