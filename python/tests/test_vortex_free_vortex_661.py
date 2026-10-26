"""VortexElement is a free vortex (#661).

It used to add the static radial-equilibrium rise of a Vatistas vortex to the
TOTAL pressure at constant Tt -- neither a free vortex (no Pt change) nor
rotor work (Tt rises with Pt) -- and so destroyed entropy for every omega > 0
(19 of the 33 sigma failures in the #598 audit). A free vortex conserves Pt
and Tt and carries r * V_theta from r_in to r_out; only the static pressure
redistributes, and that is reported.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    VortexElement,
)
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.schemas import NetworkGraphSchema

Y = list(cb.mole_to_mass(cb.species.dry_air()))
X = list(cb.species.dry_air())


def _net(**vortex_kw: float) -> FlowNetwork:
    g = FlowNetwork()
    for name, Pt in (("a", 2.4e5), ("b", 2.0e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 600.0, Y
        g.add_node(bnd)
    for n in ("p1", "p2"):
        g.add_node(PlenumNode(n))
    g.add_element(OrificeElement("o1", "a", "p1", Cd=0.8, diameter=0.02, correlation="fixed"))
    g.add_element(VortexElement("vx", "p1", "p2", **vortex_kw))
    g.add_element(OrificeElement("o2", "p2", "b", Cd=0.8, diameter=0.02, correlation="fixed"))
    return g


def _solve(g: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkSolver(g).solve()


@pytest.mark.parametrize("rpm", [3000.0, 20000.0])
def test_total_pressure_and_temperature_are_conserved(rpm: float) -> None:
    r = _solve(_net(shaft_speed_rpm=rpm, r_in=0.05, r_out=0.10))
    assert r["__converged__"], r["__message__"]
    assert r["p2.Pt"] == pytest.approx(r["p1.Pt"], rel=1e-9)
    assert r["p2.Tt"] == pytest.approx(r["p1.Tt"], rel=1e-9)
    assert r["__second_law_failures__"] == []
    assert r["__second_law__"]["vx"]["verdict"] != "fail"


def test_angular_momentum_is_carried_from_r_in_to_r_out() -> None:
    el = VortexElement("v", "a", "b", shaft_speed_rpm=6000.0, r_in=0.04, r_out=0.10)
    v_in, v_out = el.swirl()
    assert v_in == pytest.approx(6000.0 * 2.0 * math.pi / 60.0 * 0.04)
    assert v_out * 0.10 == pytest.approx(v_in * 0.04)
    given = VortexElement(
        "v", "a", "b", shaft_speed_rpm=6000.0, r_in=0.04, r_out=0.02, v_theta_in=30.0
    )
    assert given.swirl() == pytest.approx((30.0, 60.0))


def test_the_static_rise_is_the_swirl_head_given_up() -> None:
    """Outward the swirl slows and its dynamic head becomes static pressure;
    inward the reverse. At low swirl Mach it is 1/2 rho (V_in^2 - V_out^2)."""
    for r_in, r_out in ((0.05, 0.10), (0.10, 0.05)):
        r = _solve(_net(r_in=r_in, r_out=r_out, v_theta_in=40.0))
        assert r["__converged__"], r["__message__"]
        d = r["vx.v_theta_in"], r["vx.v_theta_out"]
        rho = cb.density(r["p1.T"], r["vx.P_swirl_in"], X)
        expect = 0.5 * rho * (d[0] ** 2 - d[1] ** 2)
        assert r["vx.dP_swirl_static"] == pytest.approx(expect, rel=0.01)
        assert math.copysign(1.0, r["vx.dP_swirl_static"]) == math.copysign(1.0, r_out - r_in)


def test_the_swirl_static_state_is_isentropic_from_the_stagnation_state() -> None:
    """The static P from h = h_t - V^2/2 agrees with C++'s own isentropic
    T_from_stagnation at the swirl Mach, so the two routes are one state."""
    r = _solve(_net(r_in=0.05, r_out=0.02, v_theta_in=120.0))
    assert r["__converged__"], r["__message__"]
    Tt, M = r["p1.Tt"], r["vx.M_theta_out"]
    T = cb.T_from_stagnation(Tt, M, X)
    assert M * cb.speed_of_sound(T, X) == pytest.approx(r["vx.v_theta_out"], rel=1e-6)


def test_no_swirl_at_the_axis() -> None:
    el = VortexElement("v", "a", "b", shaft_speed_rpm=3000.0, r_in=0.0)
    assert el.swirl() == (0.0, 0.0)
    with pytest.raises(ValueError, match="no angular momentum"):
        VortexElement("v", "a", "b", r_in=0.0, v_theta_in=10.0)


def test_the_old_vatistas_parameters_are_deprecated_and_ignored() -> None:
    with pytest.warns(DeprecationWarning, match="r_c and n are ignored"):
        old = VortexElement("v", "a", "b", shaft_speed_rpm=3000.0, r_c=0.02, n=2.0)
    assert old.swirl() == VortexElement("v", "a", "b", shaft_speed_rpm=3000.0).swirl()
    with pytest.raises(TypeError):
        VortexElement("v", "a", "b", 0.02, 0.1)  # the old positional r_c, r_out


def test_a_gui_file_saved_before_661_loads_as_a_free_vortex() -> None:
    """The legacy r_c and n are dropped without a warning; the GUI shows the
    note. Its r_in = 0 (the old drop default) means no swirl."""
    nodes = [
        {"id": "a", "type": "pressure_boundary", "position": {"x": 0, "y": 0}, "data": {}},
        {"id": "b", "type": "pressure_boundary", "position": {"x": 400, "y": 0}, "data": {}},
        {
            "id": "vx",
            "type": "vortex",
            "position": {"x": 200, "y": 0},
            "data": {"r_c": 0.02, "r_out": 0.1, "r_in": 0.0, "omega_rpm": 3000.0, "n": 2.0},
        },
        {
            "id": "vn",
            "type": "vortex",
            "position": {"x": 200, "y": 0},
            "data": {"omega_rpm": 3000.0, "v_theta_in": 25.0},
        },
    ]
    edges = [
        {"id": "e1", "source": "a", "target": "vx", "data": {}},
        {"id": "e2", "source": "vx", "target": "vn", "data": {}},
        {"id": "e3", "source": "vn", "target": "b", "data": {}},
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        net = build_network_from_schema(NetworkGraphSchema(nodes=nodes, edges=edges))
    assert net.elements["vx"].swirl() == (0.0, 0.0)
    assert net.elements["vn"].r_in == 0.05
    assert net.elements["vn"].swirl()[0] == 25.0


def test_omega_rpm_is_a_deprecated_alias_of_shaft_speed_rpm() -> None:
    """#773: omega is 2 pi f, so the shaft speed in rpm is shaft_speed_rpm.
    The old keyword still works, with a DeprecationWarning, until 0.8.0."""
    with pytest.warns(DeprecationWarning, match="removed in 0.8.0"):
        old = VortexElement("v", "a", "b", omega_rpm=3000.0)
    new = VortexElement("v", "a", "b", shaft_speed_rpm=3000.0)
    assert old.shaft_speed_rpm == 3000.0
    assert old.swirl() == new.swirl()
    with pytest.warns(DeprecationWarning, match="removed in 0.8.0"):
        assert new.omega_rpm == 3000.0
    with pytest.warns(DeprecationWarning), pytest.raises(TypeError, match="not both"):
        VortexElement("v", "a", "b", shaft_speed_rpm=1.0, omega_rpm=2.0)


def test_a_gui_file_saved_with_omega_rpm_still_loads() -> None:
    """A network saved before #773 names the shaft speed omega_rpm, on the
    node and in the solver settings; both still reach the element."""

    def graph(node_key: str, settings_key: str) -> dict:
        def node(nid: str, ntype: str, data: dict) -> dict:
            return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}

        return {
            "nodes": [
                node("a", "pressure_boundary", {}),
                node("b", "pressure_boundary", {}),
                node("local", "vortex", {node_key: 3000.0, "r_in": 0.05}),
                node("glob", "vortex", {}),
            ],
            "edges": [
                {"id": "e1", "source": "a", "target": "local", "data": {}},
                {"id": "e2", "source": "local", "target": "glob", "data": {}},
                {"id": "e3", "source": "glob", "target": "b", "data": {}},
            ],
            "solver_settings": {settings_key: 1500.0},
        }

    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        old = build_network_from_schema(NetworkGraphSchema(**graph("omega_rpm", "omega_rpm")))
        new = build_network_from_schema(
            NetworkGraphSchema(**graph("shaft_speed_rpm", "shaft_speed_rpm"))
        )
    for net in (old, new):
        assert net.elements["local"].shaft_speed_rpm == 3000.0
        assert net.elements["glob"].shaft_speed_rpm == 1500.0
    assert old.elements["local"].swirl() == new.elements["local"].swirl()
    v_in = old.elements["local"].swirl()[0]
    assert v_in == pytest.approx(3000.0 * 2.0 * math.pi / 60.0 * 0.05)


def test_the_diagnostics_report_shaft_speed_rpm() -> None:
    r = _solve(_net(shaft_speed_rpm=3000.0, r_in=0.05, r_out=0.10))
    assert r["__converged__"], r["__message__"]
    assert r["vx.shaft_speed_rpm"] == 3000.0
    assert "vx.omega_rpm" not in r
