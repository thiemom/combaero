"""Energy accounting semantics (#481): an injecting MassFlowBoundary mixes its
own stream in, 'fraction' scales sensible enthalpy, heat is spread over the
real flow, and heat a node cannot take is reported, not lost.
"""

from __future__ import annotations

import pytest
from _closure_check import closure

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    EnergyBoundary,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    WallNode,
)

X_AIR = cb.species.dry_air()
Y_AIR = list(cb.mole_to_mass(X_AIR))
I_CH4 = cb.species_index_from_name("CH4")
I_O2 = cb.species_index_from_name("O2")


def _pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _orifice(eid: str, a: str, b: str, d: float) -> OrificeElement:
    return OrificeElement(eid, a, b, Cd=0.7, diameter=d, correlation="fixed")


def _assert_closed(s: NetworkSolver, r: dict) -> None:
    c = closure(s, r)
    assert c is not None and abs(c["energy"]) < c["energy_tol"], c


def _sensible(T: float, X) -> float:
    return cb.h_mass(T, X) - cb.h_mass(cb.SENSIBLE_ENTHALPY_REF_T, X)


# --- an injecting MassFlowBoundary -------------------------------------------


def test_an_injection_mixes_its_own_stream_in() -> None:
    """50 g/s of hot methane-laden gas injected between two orifices: the
    node takes it at its own Tt and Y (it stayed at the inflow's 300 K and
    air composition)."""
    Y_inj = list(Y_AIR)
    Y_inj[I_CH4] += 0.2
    Y_inj = [y / sum(Y_inj) for y in Y_inj]
    g = FlowNetwork()
    g.add_node(_pb("a", 2e5, 300.0))
    g.add_node(_pb("out", 1e5, 300.0))
    g.add_node(MassFlowBoundary("inj", m_dot=0.05, Tt=800.0, Y=Y_inj))
    g.add_element(_orifice("e1", "a", "inj", 0.02))
    g.add_element(_orifice("e2", "inj", "out", 0.03))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    m1 = r["e1.m_dot"]
    assert r["e2.m_dot"] == pytest.approx(m1 + 0.05, rel=1e-9)
    T, Y, _ = s._derived_states["inj"]
    H = m1 * cb.h_mass(300.0, X_AIR) + 0.05 * cb.h_mass(800.0, cb.mass_to_mole(Y_inj))
    assert (m1 + 0.05) * cb.h_mass(T, cb.mass_to_mole(list(Y))) == pytest.approx(H, rel=1e-9)
    assert Y[I_CH4] == pytest.approx(0.05 * Y_inj[I_CH4] / (m1 + 0.05), rel=1e-9)
    _assert_closed(s, r)


# --- 'fraction' is a fraction of the sensible enthalpy ------------------------


def _fraction_net(T_in: float, fraction: float) -> tuple[NetworkSolver, dict]:
    g = FlowNetwork()
    g.add_node(_pb("a", 2e5, T_in))
    g.add_node(_pb("out", 1e5, 300.0))
    p = PlenumNode("p")
    p.add_energy_boundary(EnergyBoundary("loss", fraction=fraction))
    g.add_node(p)
    g.add_element(_orifice("e1", "a", "p", 0.02))
    g.add_element(_orifice("e2", "p", "out", 0.02))
    s = NetworkSolver(g)
    return s, s.solve()


@pytest.mark.parametrize("T_in", [600.0, 1500.0])
def test_a_five_percent_loss_removes_five_percent_of_the_sensible_heat(T_in: float) -> None:
    s, r = _fraction_net(T_in, -0.05)
    assert r["__success__"]
    T = s._derived_states["p"][0]
    assert _sensible(T, X_AIR) == pytest.approx(0.95 * _sensible(T_in, X_AIR), rel=1e-9)
    m = r["e1.m_dot"]
    assert r["p.Q_fraction"] == pytest.approx(-0.05 * m * _sensible(T_in, X_AIR), rel=1e-9)
    _assert_closed(s, r)


def test_a_loss_does_nothing_to_gas_at_the_reference_temperature() -> None:
    """On absolute h a '-5%' HEATED air at 298 K (formation enthalpy)."""
    s, r = _fraction_net(cb.SENSIBLE_ENTHALPY_REF_T, -0.05)
    assert s._derived_states["p"][0] == pytest.approx(cb.SENSIBLE_ENTHALPY_REF_T, abs=1e-9)


def test_a_combustor_loss_is_a_fraction_of_its_products_sensible_heat() -> None:
    Yf = [0.0] * len(Y_AIR)
    Yf[I_CH4] = 1.0

    def burn(f: float):
        st = [cb.MassStream(1.0, 600.0, 1e5, Y_AIR), cb.MassStream(0.03, 300.0, 1e5, Yf)]
        return _solver_tools.adiabatic_T_complete_and_jacobian_T_from_streams(st, 1e5, 0.0, f)

    r0, r1 = burn(0.0), burn(-0.05)
    Xb = cb.mass_to_mole(list(r1.Y_mix))
    assert _sensible(r1.T_mix, Xb) == pytest.approx(0.95 * _sensible(r0.T_mix, Xb), rel=1e-9)


# --- the mixer's Jacobian with heat and 'fraction' ---------------------------


@pytest.mark.parametrize(
    ("Q", "fraction"), [(0.0, 0.0), (3000.0, 0.0), (0.0, -0.05), (3000.0, -0.05)]
)
def test_the_mixer_jacobian_holds_with_heat_and_fraction(Q: float, fraction: float) -> None:
    """d/dm carried (1+f) on the composition term, and d/dY used the base
    enthalpy where the final one belongs: 46% off with Q, 280% with f."""
    Yb = list(Y_AIR)
    Yb[I_CH4] += 0.05
    Yb = [y / sum(Yb) for y in Yb]
    ms, Ts, Ys = [0.1, 0.05], [400.0, 900.0], [Y_AIR, Yb]

    def mix(ms_, Ts_, Ys_):
        st = [cb.MassStream(m, T, 1e5, Y) for m, T, Y in zip(ms_, Ts_, Ys_, strict=True)]
        return _solver_tools.mixer_from_streams_and_jacobians(st, Q=Q, fraction=fraction)

    r = mix(ms, Ts, Ys)
    for i in range(2):
        mp, mm = list(ms), list(ms)
        mp[i] += 1e-7
        mm[i] -= 1e-7
        fd = (mix(mp, Ts, Ys).T_mix - mix(mm, Ts, Ys).T_mix) / 2e-7
        assert r.dT_mix_d_stream[i].d_mdot == pytest.approx(fd, rel=1e-6)
        tp, tm = list(Ts), list(Ts)
        tp[i] += 1e-3
        tm[i] -= 1e-3
        fd = (mix(ms, tp, Ys).T_mix - mix(ms, tm, Ys).T_mix) / 2e-3
        assert r.dT_mix_d_stream[i].d_T == pytest.approx(fd, rel=1e-6)
        for k in (I_CH4, I_O2):
            yp, ym = [list(y) for y in Ys], [list(y) for y in Ys]
            yp[i][k] += 1e-6
            ym[i][k] -= 1e-6
            fd = (mix(ms, Ts, yp).T_mix - mix(ms, Ts, ym).T_mix) / 2e-6
            assert r.dT_mix_d_stream[i].d_Y[k] == pytest.approx(fd, rel=1e-6)


# --- heat over the real flow; what a stagnant node cannot take ---------------


@pytest.mark.parametrize("m", [1e-6, 1e-4, 1e-3])
def test_heat_is_spread_over_the_real_flow_down_to_a_milligram_per_second(m: float) -> None:
    """The floor was 2 g/s (20% of Q withheld at 1 g/s); it is 1 mg/s now."""
    assert cb.MIXER_HEAT_MDOT_FLOOR == 1e-6
    q = 0.5 * m * 1000.0  # a 500 K rise at cp ~ 1 kJ/kg/K
    mix = _solver_tools.mixer_from_streams_and_jacobians(
        [cb.MassStream(m, 400.0, 1e5, Y_AIR)], Q=q, fraction=0.0
    )
    dh = cb.h_mass(mix.T_mix, X_AIR) - cb.h_mass(400.0, X_AIR)
    assert m * dh == pytest.approx(q, rel=1e-9)


def test_heat_given_to_a_stagnant_node_is_reported_as_withheld() -> None:
    """A plenum in a dead-end branch (closed by a WallNode) with a 100 W source: no flow carries it anywhere,
    so the solver reports it and the network balance still closes."""
    g = FlowNetwork()
    g.add_node(_pb("a", 2e5, 300.0))
    g.add_node(_pb("out", 1e5, 300.0))
    g.add_node(PlenumNode("p"))
    dead = PlenumNode("dead")
    dead.add_energy_boundary(EnergyBoundary("heater", Q=100.0))
    g.add_node(dead)
    g.add_element(_orifice("e1", "a", "p", 0.02))
    g.add_element(_orifice("e2", "p", "out", 0.02))
    g.add_node(WallNode("end"))
    g.add_element(_orifice("e3", "p", "dead", 0.01))
    g.add_element(_orifice("e4", "dead", "end", 0.01))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    assert abs(r["e3.m_dot"]) < 1e-9
    assert r["dead.Q_withheld"] == pytest.approx(100.0, rel=1e-6)
    _assert_closed(s, r)


def test_the_dead_fuel_boundary_api_is_gone() -> None:
    """It stored a boundary nothing read: fuel counted only if wired through
    an element. An injecting MassFlowBoundary is the way now."""
    from combaero.network import CombustorNode

    assert not hasattr(CombustorNode("c"), "set_fuel_boundary")
