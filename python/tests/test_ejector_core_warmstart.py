"""The compressible ejector seeds its own cold start (EjectorElement.
seed_warm_start), so a network built through the Python API converges like
a GUI-built one. The seed used to live in the GUI backend only: cold, the
stress draws converged 2 of 22 times, seeded 8 of 22 (every critical draw).
"""

from __future__ import annotations

import math

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

Y = list(cb.species.dry_air_mass())
D = 0.05
A = math.pi * D * D / 4


def _net(regime: str = "compressible", P: tuple = (3.0e5, 10865.4, 12253.4)) -> FlowNetwork:
    """A critical-mode stress draw (gadgets.ejector, seed 4), no guesses."""
    g = FlowNetwork()
    for nid, Pt, Tt in zip(("bp", "bs", "bo"), P, (440.0, 280.0, 300.0), strict=True):
        g.add_node(PressureBoundary(nid, Pt=Pt, Tt=Tt, Y=Y))
    for nid in ("mp", "ms", "mo"):
        g.add_node(MomentumChamberNode(nid, area=A))
    for eid, a, b, L in (
        ("cp", "bp", "mp", 0.05),
        ("cs", "bs", "ms", 0.05),
        ("co", "mo", "bo", 1.0),
    ):
        g.add_element(ChannelElement(eid, a, b, length=L, diameter=D, regime="incompressible"))
    g.add_element(
        EjectorElement(
            id="ej",
            primary_node="mp",
            secondary_node="ms",
            outlet_node="mo",
            throat_area=3.14e-5,
            nozzle_exit_area=4.675e-5,
            mixing_area=3.069e-4,
            regime=regime,
        )
    )
    return g


def test_an_api_built_ejector_converges_cold_to_the_critical_point() -> None:
    s = NetworkSolver(_net())
    s.raise_internal_errors = False
    r = s.solve(timeout=30.0)
    assert r["__success__"]
    assert r["ej.critical_mode"] == 1.0
    assert r["ej.P_jct"] == pytest.approx(r["ej.P_c_star"], rel=1e-3)
    assert r["cp.m_dot"] > 0.0 and r["cs.m_dot"] > 0.0


def test_the_seed_comes_from_the_closed_forms_and_spares_user_guesses() -> None:
    g = _net()
    g.nodes["ms"].initial_guess = {"ms.P": 1.2e4, "ms.Pt": 1.2e4}
    g.resolve_all_topology()
    assert g.nodes["mp"].initial_guess == {"mp.P": 3.0e5, "mp.Pt": 3.0e5}
    assert g.nodes["ms"].initial_guess == {"ms.P": 1.2e4, "ms.Pt": 1.2e4}  # the user's
    p_jct = g.elements["ej"].initial_guess["ej.P_jct"]
    assert 10865.4 < p_jct < 3.0e5  # P_c*, between the supplies
    assert g.elements["cp"].initial_guess["cp.m_dot"] > 0.0  # choked primary flow


def test_a_changed_boundary_is_seeded_afresh() -> None:
    """The seed re-runs at every resolve_topology and replaces its own
    values: a sweep re-solving one network object is never seeded from the
    previous operating point."""
    g = _net()
    g.resolve_all_topology()
    p1 = g.elements["ej"].initial_guess["ej.P_jct"]
    g.nodes["bp"].Pt = 6.0e5
    g.resolve_all_topology()
    assert g.nodes["mp"].initial_guess["mp.Pt"] == 6.0e5
    assert g.elements["ej"].initial_guess["ej.P_jct"] > p1


def test_the_incompressible_regime_is_not_seeded() -> None:
    """Its P_jct means something else; the solver's incompressible proxy
    flips a compressible ejector there, and the seed must withdraw."""
    g = _net(regime="incompressible")
    g.resolve_all_topology()
    assert "ej.P_jct" not in getattr(g.elements["ej"], "initial_guess", {})
    g = _net()
    g.resolve_all_topology()
    g.elements["ej"].regime = "incompressible"
    g.resolve_all_topology()
    assert "ej.P_jct" not in getattr(g.elements["ej"], "initial_guess", {})
    assert "mp.Pt" not in getattr(g.nodes["mp"], "initial_guess", {})
