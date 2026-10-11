"""Global energy balance of the cooling configurations (#471).

The network is the control volume. Enthalpy enters with the streams from the
boundaries and leaves with the streams into them; walls move heat INSIDE the
network, so they cancel -- except heat a wall injects into a node that cannot
take it, a BOUNDARY: that heat leaves with the stream crossing the boundary,
and the solver reports it as the node's ``Q_wall_out`` (and the wall's
``Q_to_boundary``). With it the balance closes to round-off:

    sum_out (m h(T_upstream)) + sum_boundaries Q_wall_out - sum_in m h(T_in) = 0
"""

from __future__ import annotations

import io
from contextlib import redirect_stdout

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

with redirect_stdout(io.StringIO()):
    import test_effusion_liner as liner_cases
    import test_effusion_wall as plate_cases
    import test_impingement_array as array_cases

X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)


def _balance(net: FlowNetwork, s: NetworkSolver, r: dict, with_boundary_heat: bool = True):
    """(out - in, in) enthalpy flow [W] across the network's boundary."""
    bnd = {n for n, o in net.nodes.items() if isinstance(o, (PressureBoundary, MassFlowBoundary))}

    def h(node: str) -> float:
        T = net.nodes[node].Tt if node in bnd else s._derived_states[node][0]
        return cb.h_mass(T, X_AIR)

    e_in = e_out = 0.0
    for eid, e in net.elements.items():
        m = r["__element_diag__"][eid]["m_dot"]
        f, t = e.from_node, e.to_node
        if m < 0.0:
            f, t, m = t, f, -m
        if f in bnd:
            e_in += m * h(f)
        if t in bnd:
            e_out += m * h(f)
    if with_boundary_heat:
        e_out += sum(r.get(f"{b}.Q_wall_out", 0.0) for b in bnd)
    return e_out - e_in, e_in


def _assert_closed(net, s, r) -> None:
    assert r["__success__"], r.get("__message__")
    d, e_in = _balance(net, s, r)
    assert abs(d) < 1e-9 * e_in, d


def test_effusion_plate_into_a_merge_chamber() -> None:
    s, r = plate_cases._liner()
    _assert_closed(s.network, s, r)


@pytest.mark.parametrize("chamber", [False, True])
def test_effusion_liner_with_its_bypass(chamber: bool) -> None:
    g = liner_cases._net(liner_cases._liner(), chamber=chamber)
    s = NetworkSolver(g)
    _assert_closed(g, s, s.solve())


def test_impingement_array_whose_hot_side_dumps_into_a_boundary() -> None:
    """The hot duct runs boundary to boundary, so the heat its walls remove
    lands on the outlet boundary: without Q_wall_out the balance is open by
    exactly that heat (the wall's whole Q)."""
    g = array_cases._net(array_cases._array(5), hot=True)
    s = NetworkSolver(g)
    r = s.solve()
    _assert_closed(g, s, r)
    open_by, _ = _balance(g, s, r, with_boundary_heat=False)
    q_walls = sum(r[f"w__r{i}.Q"] for i in range(1, 6))
    assert open_by == pytest.approx(q_walls, rel=1e-6)
    assert sum(r[f"w__r{i}.Q_to_boundary"] for i in range(1, 6)) == pytest.approx(-q_walls)
    assert r["hot_out.Q_wall_out"] == pytest.approx(-q_walls)


def test_a_wall_between_internal_nodes_needs_no_boundary_term() -> None:
    """Both streams pass a plenum before their outlet: the wall's heat stays
    inside, nothing reaches a boundary, and the balance closes as is."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("hot_in", m_dot=0.1, Tt=800.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("cold_in", m_dot=0.05, Tt=400.0, Y=Y_AIR))
    for n in ("hot_out", "cold_out"):
        g.add_node(PressureBoundary(n, Pt=2e5, Tt=300.0, Y=Y_AIR))
    for n in ("hot_p", "cold_p"):
        g.add_node(PlenumNode(n))
    area = 3.14159 * 0.04
    g.add_element(
        ChannelElement(
            "hot",
            "hot_in",
            "hot_p",
            length=1.0,
            diameter=0.04,
            roughness=0.0,
            surface=ConvectiveSurface(area=area),
        )
    )
    g.add_element(
        ChannelElement(
            "cold",
            "cold_in",
            "cold_p",
            length=1.0,
            diameter=0.04,
            roughness=0.0,
            surface=ConvectiveSurface(area=area),
        )
    )
    g.add_element(
        OrificeElement("ho", "hot_p", "hot_out", Cd=0.8, diameter=0.035682, correlation="fixed")
    )
    g.add_element(
        OrificeElement("co", "cold_p", "cold_out", Cd=0.8, diameter=0.035682, correlation="fixed")
    )
    g.add_wall(
        ThermalWall(id="w", element_a="hot", element_b="cold", layers=[WallLayer(0.002, 25.0)])
    )
    s = NetworkSolver(g)
    r = s.solve()
    _assert_closed(g, s, r)
    assert r["w.Q"] > 0.0
    assert r["w.Q_to_boundary"] == 0.0


@pytest.mark.parametrize("m", [2e-3, 0.01, 0.1, 1.0])
def test_the_mixer_deposits_the_whole_heat(m: float) -> None:
    """Q / m exactly at and above kMixerHeatMdotFloor (2 g/s): the former
    sqrt(m^2 + 1e-6) withheld 5e-5 of Q at 0.1 kg/s and 11% at 2 g/s."""
    from combaero import _solver_tools

    for q in (2000.0, -50.0 * m):
        mix = _solver_tools.mixer_from_streams_and_jacobians(
            [cb.MassStream(m, 800.0, 1e5, Y_AIR)], Q=q, fraction=0.0
        )
        dh = cb.h_mass(mix.T_mix, X_AIR) - cb.h_mass(800.0, X_AIR)
        assert m * dh == pytest.approx(q, rel=1e-10, abs=1e-6)


@pytest.mark.parametrize("m", [5e-4, 1.9e-3, 2.1e-3, 0.05])
def test_the_mixer_heat_term_derivative_across_the_floor(m: float) -> None:
    """dT_mix/dm through Q/m_eff on both sides of the floor's knee."""
    from combaero import _solver_tools

    def t_mix(mm: float) -> float:
        return _solver_tools.mixer_from_streams_and_jacobians(
            [cb.MassStream(mm, 800.0, 1e5, Y_AIR)], Q=5.0, fraction=0.0
        ).T_mix

    mix = _solver_tools.mixer_from_streams_and_jacobians(
        [cb.MassStream(m, 800.0, 1e5, Y_AIR)], Q=5.0, fraction=0.0
    )
    step = 1e-6 * m
    fd = (t_mix(m + step) - t_mix(m - step)) / (2.0 * step)
    assert mix.dT_mix_d_stream[0].d_mdot == pytest.approx(fd, rel=1e-5)


def test_solving_the_same_network_again_does_not_stack_wall_heat() -> None:
    """Each solver hangs a wall EnergyBoundary on the heated nodes. A second
    solver (or a retry) used to add another beside the first, which kept the
    previous solve's heat -- counted twice."""
    g = array_cases._net(array_cases._array(3), hot=True)
    first = NetworkSolver(g).solve()
    s = NetworkSolver(g)
    again = s.solve()
    _assert_closed(g, s, again)
    for i in (1, 2, 3):
        assert again[f"w__r{i}.Q"] == pytest.approx(first[f"w__r{i}.Q"], rel=1e-9)
    for node in g.nodes.values():
        ids = [eb.id for eb in getattr(node, "energy_boundaries", [])]
        assert ids.count(f"_wall_{node.id}") <= 1


def test_a_failed_cold_solve_is_retried_from_the_wall_free_flows(monkeypatch) -> None:
    """The wall-free solution seeds the retry; the result is the coupled one
    and its balance closes (no heat left over from the failed attempt)."""
    g = array_cases._net(array_cases._array(3), hot=True)
    reference = NetworkSolver(g).solve()
    assert "wall-free" not in reference["__message__"]

    original = NetworkSolver._solve_flow_retries

    def fail_once(self, **kw):
        r = original(self, **kw)
        if self.network.thermal_coupling_enabled:
            r = dict(r, __success__=False, __message__="forced", __final_norm__=1.0)
        return r

    monkeypatch.setattr(NetworkSolver, "_solve_flow_retries", fail_once)
    s = NetworkSolver(g)
    r = s.solve()
    assert "wall-free" in r["__message__"]
    assert g.thermal_coupling_enabled
    _assert_closed(g, s, r)
    for i in (1, 2, 3):
        assert r[f"w__r{i}.Q"] == pytest.approx(reference[f"w__r{i}.Q"], rel=1e-6)
