"""An ejector root with a port flowing backwards is not a solution (#619).

Every regime the element models has its streams going forward. Past the jet
pump's shut-off head the suction reverses, and the rows' C1 continuation there
has roots no ejector produces: the stress family's seed 36, driven to an
outlet of 171.6 kPa, "converged" at omega = -169 with 0.47 kg/s entering the
1.7 bar outlet, rising to a 4 bar mixing plane and leaving through the 1 bar
suction.
"""

from __future__ import annotations

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

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
D = 0.05
A = 3.141592653589793 * D * D / 4.0


def _jet_pump(p_out: float) -> FlowNetwork:
    """Seed 36's geometry: incompressible jet pump, 4 bar primary, 1 bar
    suction; shut-off near 1.25 bar at the outlet."""
    g = FlowNetwork()
    for nid, p, tt in (("bp", 4e5, 300.0), ("bs", 1e5, 280.0), ("bo", p_out, 300.0)):
        g.add_node(PressureBoundary(nid, Pt=p, Tt=tt, Y=Y_AIR))
    for nid, p in (("mp", 4e5), ("ms", 1e5), ("mo", p_out)):
        mc = MomentumChamberNode(nid, area=A, regime="incompressible")
        mc.initial_guess = {f"{nid}.P": p, f"{nid}.Pt": p}
        g.add_node(mc)
    for eid, a, b, L in (
        ("cp", "bp", "mp", 1.0),
        ("cs", "bs", "ms", 1.0),
        ("co", "mo", "bo", 0.05),
    ):
        g.add_element(ChannelElement(eid, a, b, length=L, diameter=D, regime="incompressible"))
    g.add_element(
        EjectorElement(
            id="ej",
            primary_node="mp",
            secondary_node="ms",
            outlet_node="mo",
            throat_area=1.5e-5,
            nozzle_exit_area=5.929498449620368e-05,
            mixing_area=6.897920631859684e-04,
            regime="incompressible",
        )
    )
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=30)
    return s, r


def test_below_shut_off_the_jet_pump_solves_forward() -> None:
    s, r = _solve(_jet_pump(1.159e5))
    assert r["__success__"], r["__message__"]
    assert r["__consistent__"] is True
    assert r["ej.omega"] > 0.0


def test_past_shut_off_no_backwards_root_is_reported_as_a_solution() -> None:
    """At 171.6 kPa the solver used to return omega = -169 as converged and
    consistent. Whatever it finds now, a reported solution has every port
    going forward."""
    _, r = _solve(_jet_pump(171555.0))
    if r["__success__"]:
        assert r["ej.omega"] >= 0.0
        assert r["co.m_dot"] >= 0.0


@pytest.mark.parametrize(
    ("port", "factor", "ok"),
    [("cs", -1.0, False), ("co", -1.0, False), ("cp", -1.0, False), ("cs", 1.0, True)],
)
def test_the_verifier_rejects_a_backwards_port(port: str, factor: float, ok: bool) -> None:
    s, r = _solve(_jet_pump(1.159e5))
    assert r["__success__"]
    sol = dict(zip(s.unknown_names, r["__x_solution__"], strict=True))
    sol[f"{port}.m_dot"] *= factor
    assert s.network.elements["ej"].verify_solution_consistent(sol) is ok


def test_a_flow_at_solver_precision_of_zero_is_not_backwards() -> None:
    """Exactly at shut-off the suction flow is zero to solver precision."""
    s, r = _solve(_jet_pump(1.159e5))
    sol = dict(zip(s.unknown_names, r["__x_solution__"], strict=True))
    sol["cs.m_dot"] = -1e-9 * abs(sol["cp.m_dot"])
    assert s.network.elements["ej"].verify_solution_consistent(sol)
