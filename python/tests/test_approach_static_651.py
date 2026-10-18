"""The orifice and the incompressible channel read their approach pipe's static
state, not their node's (#651).

A PlenumNode has no velocity by definition, so between a pipe and an orifice
it reported the pipe's STAGNATION pressure as its static one: the orifice's
p1 and rho1, and the next pipe's friction density, came out a dynamic head
too high (1.5% flow at beta 0.7, pipe Mach 0.18). Both now take the static
state of the flow at their pipe's mass flux from the node's (Pt, Tt), via
C++'s approach_static_pressure, so the node type between pipes no longer
moves the flow.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


@pytest.mark.parametrize("incompressible", [True, False])
@pytest.mark.parametrize("G", [0.0, 20.0, 150.0])
def test_the_approach_state_is_the_momentum_chambers_closure(
    incompressible: bool, G: float
) -> None:
    """The same state a MomentumChamberNode of the pipe's area would carry:
    its closure's residual is zero there, in both regimes."""
    Pt, Tt, A = 1.5e5, 350.0, 0.01
    a = _solver_tools.approach_static_pressure(Pt, Tt, G, Y_AIR, incompressible)
    assert not a.choked
    r = _solver_tools.momentum_chamber_residual_and_jacobian(
        a.P, Pt, G * A, Tt, Y_AIR, A, incompressible
    )
    assert abs(r.residual) <= 1e-7 * Pt
    if G == 0.0:
        assert pytest.approx(Pt, rel=1e-14) == a.P


@pytest.mark.parametrize("incompressible", [True, False])
def test_the_approach_state_derivatives_are_exact(incompressible: bool) -> None:
    x0 = [1.5e5, 350.0, 120.0]

    def P(x: list[float]) -> float:
        return _solver_tools.approach_static_pressure(x[0], x[1], x[2], Y_AIR, incompressible).P

    a = _solver_tools.approach_static_pressure(*x0, Y_AIR, incompressible)
    for i, name in enumerate(("dP_dPt", "dP_dTt", "dP_dG")):
        h = 1e-6 * x0[i]
        up, dn = list(x0), list(x0)
        up[i] += h
        dn[i] -= h
        fd = (P(up) - P(dn)) / (2.0 * h)
        assert getattr(a, name) == pytest.approx(fd, rel=1e-6), name


def test_past_its_maximum_flux_the_incompressible_state_says_so() -> None:
    """Pt = P + G^2 R T / (2 P) has no root above G = Pt / sqrt(2 R T)."""
    a = _solver_tools.approach_static_pressure(1.0e5, 300.0, 300.0, Y_AIR, True)
    assert a.choked and pytest.approx(0.5e5) == a.P


def _pipe_orifice_sink(node_cls: type, regime: str, beta: float) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.2e5, Tt=300.0))
    g.add_node(node_cls("a", regime=regime) if node_cls is MomentumChamberNode else node_cls("a"))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    g.add_element(ChannelElement("p1", "hi", "a", length=0.5, diameter=0.1, regime=regime))
    g.add_element(
        OrificeElement(
            "o", "a", "ho", Cd=0.6, diameter=beta * 0.1, regime=regime, correlation="fixed"
        )
    )
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    return s, r


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize("beta", [0.5, 0.7])
def test_into_a_sink_the_node_type_does_not_move_the_flow(regime: str, beta: float) -> None:
    """pipe -> node -> orifice -> sink, the topology with no recovery: the
    node type was 1.5% of flow at beta 0.7 before."""
    _, rp = _solve(_pipe_orifice_sink(PlenumNode, regime, beta))
    _, rc = _solve(_pipe_orifice_sink(MomentumChamberNode, regime, beta))
    assert rc["o.m_dot"] == pytest.approx(rp["o.m_dot"], rel=1e-6)


def test_a_plenum_fed_orifice_reads_its_node_unchanged() -> None:
    """No pipe upstream, so no E and no approach state: the kernel reads the
    node's stagnation state as before."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.2e5, Tt=300.0))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    g.add_element(OrificeElement("o", "hi", "ho", Cd=0.6, diameter=0.02, correlation="fixed"))
    _, r = _solve(g)
    o = g.elements["o"]
    assert o.beta == 0.0 and not o._approach_reads_flow()
    assert "o.P_approach" not in r


@pytest.mark.parametrize("node_cls", [PlenumNode, MomentumChamberNode])
def test_the_incompressible_channels_density_columns_are_exact(node_cls: type) -> None:
    """The channel out of node a reads its density at a's (Pt, T) and its own
    flux: the Pt, T and flow columns against central differences."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.3e5, Tt=300.0))
    g.add_node(
        node_cls("a", regime="incompressible") if node_cls is MomentumChamberNode else node_cls("a")
    )
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    g.add_element(
        OrificeElement(
            "o", "hi", "a", Cd=0.6, diameter=0.06, regime="incompressible", correlation="fixed"
        )
    )
    g.add_element(
        ChannelElement("p", "a", "ho", length=2.0, diameter=0.08, regime="incompressible")
    )
    s, r = _solve(g)
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    rows = s._row_kinds(x)
    i = next(k for k, (owner, _) in enumerate(rows) if "'p'" in owner)
    scale = (np.abs(J[i]) * np.abs(x)).max()
    for j, name in enumerate(s.unknown_names):
        h = 1e-6 * max(abs(x[j]), 1e-3)
        e = np.zeros_like(x)
        e[j] = h
        fd = (F(x + e)[i] - F(x - e)[i]) / (2 * h)
        assert J[i, j] * abs(x[j]) == pytest.approx(fd * abs(x[j]), rel=1e-5, abs=1e-7 * scale), (
            name
        )


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_the_approach_orifice_row_is_exact_into_a_sink(regime: str) -> None:
    """The orifice's own consistency solve (m* = m_calc(m*, x)) and its
    implicit-function columns, no recovery involved."""
    s, r = _solve(_pipe_orifice_sink(PlenumNode, regime, 0.7))
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    rows = s._row_kinds(x)
    i = next(k for k, (owner, _) in enumerate(rows) if "'o'" in owner)
    scale = (np.abs(J[i]) * np.abs(x)).max()
    for j, name in enumerate(s.unknown_names):
        h = 1e-6 * max(abs(x[j]), 1e-3)
        e = np.zeros_like(x)
        e[j] = h
        fd = (F(x + e)[i] - F(x - e)[i]) / (2 * h)
        assert J[i, j] * abs(x[j]) == pytest.approx(fd * abs(x[j]), rel=1e-5, abs=1e-7 * scale), (
            name
        )
    # The row reads the plenum's Pt (the sink is a boundary, no unknown).
    assert J[i, s.unknown_names.index("a.Pt")] != 0.0
    # Consistent: the row's own flow column is exactly 1 (fixed Cd).
    assert J[i, s.unknown_names.index("o.m_dot")] == 1.0


def test_between_internal_nodes_the_stall_is_gone() -> None:
    """test_second_law's heated pair (beta 0.89 orifices behind plenums into a
    sink at 2 bar): with p1 read at the iterate's flow the forward kernel ran
    backwards for 0 < Pt - p2 < q_pipe and the solve stalled at |F| = 0.12.
    The consistency solve gives m* the sign of Pt - p2."""
    from combaero.network import MassFlowBoundary

    g = FlowNetwork()
    g.add_node(MassFlowBoundary("in", m_dot=0.1, Tt=800.0, Y=Y_AIR))
    g.add_node(PlenumNode("p"))
    g.add_node(PressureBoundary("out", Pt=2e5, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("c", "in", "p", length=1.0, diameter=0.04, roughness=0.0))
    g.add_element(OrificeElement("o", "p", "out", Cd=0.8, diameter=0.0357, correlation="fixed"))
    _, r = _solve(g)
    q_pipe = (0.1 / (math.pi * 0.02**2)) ** 2 / (2.0 * r["o.P_approach"] / (287.0 * 800.0))
    # The plenum sits above the sink by about the pipe's dynamic head: E and
    # beta 0.89 make the orifice's own drop small.
    assert r["p.Pt"] - 2e5 > 0.5 * q_pipe
