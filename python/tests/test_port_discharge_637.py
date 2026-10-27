"""A junction port's chamber discharging into a pressure boundary (#637).

The GUI joins a port chamber wired straight to a boundary with a link. A
Pt-equal link made the reservoir recover the jet's whole head, which has no
root once the chamber is sonic: a critical-mode ejector's outlet stalled
below its P_c*. PortDischargeElement takes the boundary at the chamber's
STATIC pressure (#360's exit-head rule) and, choked, holds the chamber at the
sonic static pressure P* of its flux (under-expanded discharge).
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
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    PortDischargeElement,
    PressureBoundary,
)
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.schemas import NetworkGraphSchema

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=60)
    assert r["__success__"], r["__message__"]
    return s, r


# --- the sonic static pressure (C++) ------------------------------------------


@pytest.mark.parametrize(("T0", "G"), [(300.0, 250.0), (1600.0, 80.0), (900.0, 1.0)])
def test_the_sonic_static_pressure_is_where_the_flux_is_sonic(T0: float, G: float) -> None:
    """The chamber closure at P* reads Mach 1, and the T0 and G derivatives
    are exact."""
    s = _solver_tools.sonic_static_pressure(T0, G, Y_AIR)
    A = 1e-3
    c = _solver_tools.momentum_chamber_residual_and_jacobian(
        s.P, 2.0 * s.P, G * A, T0, Y_AIR, A, False
    )
    assert pytest.approx(1.0, abs=1e-10) == c.M
    hT, hG = 1e-4 * T0, 1e-5 * G
    fdT = (
        _solver_tools.sonic_static_pressure(T0 + hT, G, Y_AIR).P
        - _solver_tools.sonic_static_pressure(T0 - hT, G, Y_AIR).P
    ) / (2 * hT)
    fdG = (
        _solver_tools.sonic_static_pressure(T0, G + hG, Y_AIR).P
        - _solver_tools.sonic_static_pressure(T0, G - hG, Y_AIR).P
    ) / (2 * hG)
    assert s.dP_dTt == pytest.approx(fdT, rel=1e-7)
    assert s.dP_dG == pytest.approx(fdG, rel=1e-7)


# --- the element ----------------------------------------------------------------


def _chamber_to(
    p_out: float, coupling: str = "auto", Pt_in: float = 5e5, lossy_inlet: bool = False
) -> FlowNetwork:
    """reservoir -> inlet -> chamber (5e-4 m^2) -> discharge -> boundary. The
    inlet is lossless (the chamber then carries the reservoir's Pt), or a
    pipe of the chamber's bore where the discharge itself is Pt-equal: two
    Pt-equal links in series would have no root."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=Pt_in, Tt=300.0))
    g.add_node(MomentumChamberNode("c", area=5e-4))
    g.add_node(PressureBoundary("lo", Pt=p_out, Tt=300.0, coupling=coupling))
    if lossy_inlet:
        D = math.sqrt(4.0 * 5e-4 / math.pi)
        g.add_element(ChannelElement("l", "hi", "c", length=2.0, diameter=D))
    else:
        g.add_element(LosslessConnectionElement("l", "hi", "c"))
    g.add_element(PortDischargeElement("d", "c", "lo"))
    return g


def test_subsonic_the_chamber_static_is_the_boundary_pressure() -> None:
    _, r = _solve(_chamber_to(4.0e5))
    assert r["c.P"] == pytest.approx(4.0e5, rel=1e-7)
    assert r["c.Pt"] == pytest.approx(5.0e5, rel=1e-7)
    assert r["d.exit_choked"] == 0.0


@pytest.mark.parametrize("p_out", [0.5e5, 1.0e5, 2.0e5])
def test_choked_the_discharge_is_under_expanded_at_the_sonic_pressure(p_out: float) -> None:
    """Below the critical ratio the chamber sits at Mach 1 and P*, above the
    boundary, and the flow no longer depends on the boundary: the choked flux
    of the chamber's area at its Pt."""
    _, r = _solve(_chamber_to(p_out))
    assert r["d.exit_choked"] == 1.0
    assert r["d.mach_in"] == pytest.approx(1.0, abs=1e-4)
    assert r["c.P"] == pytest.approx(r["d.P_sonic"], rel=1e-6)
    assert r["c.P"] > p_out
    pr = cb.critical_pressure_ratio(300.0, 5e5, cb.species.dry_air())
    G_star = cb.mass_flux_isentropic(300.0, 5e5, pr * 5e5, cb.species.dry_air())
    assert r["d.m_dot"] == pytest.approx(G_star * 5e-4, rel=1e-4)


def test_a_total_coupled_boundary_keeps_the_pt_equality() -> None:
    _, r = _solve(_chamber_to(4.0e5, coupling="total", lossy_inlet=True))
    assert r["c.Pt"] == pytest.approx(4.0e5, rel=1e-7)


def test_flow_into_the_chamber_takes_the_boundary_as_a_stagnation_source() -> None:
    """Reversed (boundary above the reservoir) the link is Pt-equal."""
    _, r = _solve(_chamber_to(4.0e5, Pt_in=3.0e5, lossy_inlet=True))
    assert r["d.m_dot"] < 0.0
    assert r["c.Pt"] == pytest.approx(4.0e5, rel=1e-7)


@pytest.mark.parametrize("p_out", [4.0e5, 1.0e5, 2.6e5])
def test_the_discharge_rows_jacobian_is_exact(p_out: float) -> None:
    """Subsonic, choked, and at the choke blend: every column against central
    differences in the J_ij |x_j| metric."""
    s, r = _solve(_chamber_to(p_out))
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    rows = s._row_kinds(x)
    i = next(k for k, (owner, _) in enumerate(rows) if "'d'" in owner)
    scale = np.abs(J[i]) * np.abs(x)
    for j, name in enumerate(s.unknown_names):
        h = 1e-7 * max(abs(x[j]), 1e-3)
        e = np.zeros_like(x)
        e[j] = h
        fd = (F(x + e)[i] - F(x - e)[i]) / (2 * h)
        assert J[i, j] * abs(x[j]) == pytest.approx(
            fd * abs(x[j]), rel=1e-5, abs=1e-7 * scale.max()
        ), name


# --- the GUI network of #637 ----------------------------------------------------


def _node(i: str, t: str, d: dict, x: int = 0) -> dict:
    return {"id": i, "type": t, "position": {"x": x, "y": 0}, "data": {"id": i, **d}}


def _edge(i: str, s: str, t: str, sh: str | None = None, th: str | None = None) -> dict:
    e = {"id": i, "source": s, "target": t, "data": {}}
    if sh:
        e["sourceHandle"] = sh
    if th:
        e["targetHandle"] = th
    return e


def _ejector_into(p_out: float) -> dict:
    """The issue's network: three pressure boundaries wired straight onto the
    ejector's ports, mixing area unset (the 8e-4 m^2 default)."""
    return {
        "nodes": [
            _node("pp", "pressure_boundary", {"Pt": 9.01325e5, "Tt": 300}),
            _node("ps", "pressure_boundary", {"Pt": 1.02325e5, "Tt": 300}),
            _node(
                "ej",
                "ejector",
                {
                    "throat_area": 3.14e-5,
                    "nozzle_exit_area": 1e-4,
                    "recovery_efficiency": 1,
                    # K&D's anchors (P_c* ~1.2 bar): the critical-mode test
                    # below places its outlets in that window; the default
                    # huang1999 pair is tested in test_ejector_huang_anchors.py.
                    "anchor_closure": "kracik_dvorak",
                },
                200,
            ),
            _node("po", "pressure_boundary", {"Pt": p_out, "Tt": 300}, 400),
        ],
        "edges": [
            _edge("e1", "pp", "ej", th="port-primary-target"),
            _edge("e2", "ps", "ej", th="port-secondary-target"),
            _edge("e3", "ej", "po", sh="port-outlet-source", th="flow-target"),
        ],
        "solver_settings": {"global_regime": "compressible"},
    }


def test_the_gui_discharges_an_outlet_port_into_a_boundary_through_the_new_link() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_ejector_into(1.013e5)))
    kinds = {type(e).__name__ for e in net.elements.values() if e.to_node == "po"}
    assert kinds == {"PortDischargeElement"}
    # The inlet ports stay Pt-equal: a reservoir feeding a chamber.
    inlet_links = [e for e in net.elements.values() if e.from_node in ("pp", "ps")]
    assert {type(e).__name__ for e in inlet_links} == {"LosslessConnectionElement"}


@pytest.mark.slow
def test_a_critical_mode_ejector_into_a_low_outlet_converges_and_decouples() -> None:
    """Below P_c* (about 1.2 bar here) the ejector is double-choked and its
    entrainment no longer depends on the outlet; it used to stall there (|F|
    2.2e3 at 1.013 bar). Above it, entrainment falls with back pressure."""
    omega = {}
    for p_out in (0.9e5, 1.013e5, 1.5e5, 2.0e5):
        _, r = _solve(build_network_from_schema(NetworkGraphSchema(**_ejector_into(p_out))))
        omega[p_out] = r["ej.entrainment_ratio"]
        if p_out < 1.2e5:
            assert r["ej.critical_mode"] == 1.0
        ch = "__tee_jct__ej_outlet"
        if r["__auto_link__" + ch + "__po.exit_choked"] == 0.0:
            assert r[ch + ".P"] == pytest.approx(p_out, rel=1e-6)
    assert omega[0.9e5] == pytest.approx(omega[1.013e5], rel=1e-6)
    assert omega[1.013e5] > omega[1.5e5] > omega[2.0e5]
    assert math.isfinite(omega[2.0e5])
