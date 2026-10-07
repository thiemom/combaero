"""ImpingementPlateElement and ImpingementCrossflowElement, unit level (#465).

The array-level evidence against Florschuetz et al. (1981) is in
test_impingement_plate_validation.py. Here: geometry, the Gc/Gj definition,
and every derivative the solver relies on -- element-level against central
differences, and the global Jacobian of a wall-coupled chain against finite
differences, which is the only place the neighbour-flow relay is visible.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.optimize._numdiff import approx_derivative

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    EffusionPlateElement,
    FlowNetwork,
    ImpingementCrossflowElement,
    ImpingementPlateElement,
    NetworkMixtureState,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)
from combaero.network.schema import validate_diagnostics

D = 0.00254
Y_AIR = cb.mole_to_mass(cb.species.dry_air())


def _plate(**kw) -> ImpingementPlateElement:
    args = {
        "d_jet": D,
        "xn_d": 5.0,
        "yn_d": 4.0,
        "z_d": 2.0,
        "span": 12 * 4.0 * D,
        "plate_thickness": D,
    }
    args.update(kw)
    return ImpingementPlateElement("p", "s", "c", **args)


def _state(m_dot: float, T: float = 300.0, P: float = 1.03e5) -> NetworkMixtureState:
    return NetworkMixtureState(T=T, P=P, Pt=P, Tt=T, m_dot=m_dot, Y=Y_AIR)


def test_geometry_follows_the_plate_design() -> None:
    p = _plate(span=0.1)
    assert p.n_holes == round(0.1 / (4.0 * D))
    assert p.hole_count_exact == pytest.approx(0.1 / (4.0 * D))
    assert p.area == pytest.approx(p.n_holes * math.pi * D * D / 4.0)
    assert p.surface.area == pytest.approx(p.n_holes * 5.0 * 4.0 * D * D)
    assert p.Cd == cb.FLORSCHUETZ_1981_DEFAULT_CD


def test_gc_gj_is_florschuetz_mass_velocity_ratio() -> None:
    """Gc on the channel cross-section z * W, Gj on the hole area n pi d^2/4."""
    p = _plate()
    p._crossflow_sources = ["x0"]
    m_j, m_c = 0.004, 0.012
    r = p.htc_and_T(_state(m_j), flows={"x0": m_c})
    W = p.n_holes * p.yn_d * D
    Gc = m_c / (p.z_d * D * W)
    Gj = m_j / (p.n_holes * math.pi * D * D / 4.0)
    assert r.Gc_Gj == pytest.approx(Gc / Gj, rel=1e-12)


def test_h_is_the_correlation_at_the_rows_own_re_and_gc() -> None:
    p = _plate()
    p._crossflow_sources = ["x0"]
    r = p.htc_and_T(_state(0.004), flows={"x0": 0.012})
    tr = cb.complete_state(300.0, 1.03e5, cb.species.dry_air()).transport
    jet = cb.jet_array_impingement_nu(
        cb.florschuetz_1981_inline(), r.Re, r.Gc_Gj, tr.Pr, 5.0, 4.0, 2.0
    )
    assert r.Re == pytest.approx(4.0 * 0.004 / (p.n_holes * math.pi * D * tr.mu), rel=1e-9)
    assert r.h == pytest.approx(jet.Nu * tr.k / D, rel=1e-9)
    assert r.T_aw == 300.0


def test_staggered_pattern_selects_the_staggered_set() -> None:
    a = _plate(pattern="inline").htc_and_T(_state(0.004)).Nu
    b = _plate(pattern="staggered").htc_and_T(_state(0.004)).Nu
    assert a != b


@pytest.mark.parametrize("m_c", [0.0, 0.003, 0.012])
def test_h_derivatives_match_central_differences(m_c: float) -> None:
    p = _plate()
    p._crossflow_sources = ["x0", "x1"]
    flows = {"x0": m_c, "x1": 0.5 * m_c}
    m_j = 0.004
    r = p.htc_and_T(_state(m_j), flows=flows)

    e = 1e-7
    fd_j = (
        p.htc_and_T(_state(m_j + e), flows=flows).h - p.htc_and_T(_state(m_j - e), flows=flows).h
    ) / (2 * e)
    assert r.dh_dmdot == pytest.approx(fd_j, rel=1e-6)
    for src in ("x0", "x1"):
        up = dict(flows, **{src: flows[src] + e})
        dn = dict(flows, **{src: flows[src] - e})
        fd = (p.htc_and_T(_state(m_j), flows=up).h - p.htc_and_T(_state(m_j), flows=dn).h) / (2 * e)
        assert r.dh_dsources[src] == pytest.approx(fd, rel=1e-6, abs=1e-6)
    if m_c > 0:
        assert r.dh_dsources["x0"] < 0.0  # crossflow degrades the target h


def test_lichtarowicz_refuses_a_one_diameter_plate() -> None:
    g = FlowNetwork()
    for n in ("s", "c"):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = (1.03e5 if n == "s" else 1.0e5), 300.0, Y_AIR
        g.add_node(b)
    g.add_element(_plate(correlation="Lichtarowicz"))
    with pytest.raises(ValueError, match="'p'"):
        NetworkSolver(g).solve()


def test_metering_correlations_are_refused() -> None:
    p = _plate(correlation="Stolz")
    g = FlowNetwork()
    g.add_node(PlenumNode("s"))
    g.add_node(PlenumNode("c"))
    g.add_element(p)
    p.resolve_topology(g)
    with pytest.raises(ValueError, match="no\\s+pipe"):
        p.validate()


def _segment_with_neighbours() -> tuple[ImpingementCrossflowElement, list[str]]:
    """x2 between c2 and c3; x1 arrives at c2, plate p3 merges at c3."""
    g = FlowNetwork()
    for n in ("s", "c2", "c3", "c1"):
        g.add_node(PlenumNode(n))
    seg = ImpingementCrossflowElement("x2", "c2", "c3", length=5 * D, height=2 * D, span=0.1)
    g.add_element(
        ImpingementCrossflowElement("x1", "c1", "c2", length=5 * D, height=2 * D, span=0.1)
    )
    g.add_element(seg)
    for row in (2, 3):
        g.add_element(
            ImpingementPlateElement(
                f"p{row}",
                "s",
                f"c{row}",
                d_jet=D,
                xn_d=5.0,
                yn_d=4.0,
                z_d=2.0,
                span=0.1,
                plate_thickness=D,
            )
        )
    seg.resolve_topology(g)
    return seg, [eid for eid, _ in seg.network_flow_inputs()]


def test_the_segment_reads_the_right_neighbours() -> None:
    seg, inputs = _segment_with_neighbours()
    assert seg.network_flow_inputs() == [("p3", "c3"), ("x1", "c2")]
    assert sorted(inputs) == ["p3", "x1"]


def test_the_momentum_term_is_half_of_both_adjacent_merges() -> None:
    seg, _ = _segment_with_neighbours()
    m, m_arr, m_jet = 0.02, 0.015, 0.006
    st = _state(m, P=1.0e5)
    dP = seg._momentum_drop(st, {"x1": m_arr, "p3": m_jet})[0].dP
    rho = cb.density(300.0, 1.0e5, cb.species.dry_air())
    A = 2 * D * 0.1
    merge_up = (m**2 - m_arr**2) / (rho * A * A)  # row 2, arriving -> this segment
    merge_dn = ((m + m_jet) ** 2 - m**2) / (rho * A * A)  # row 3
    assert dP == pytest.approx(0.5 * (merge_up + merge_dn), rel=1e-9)


def test_segment_residual_jacobian_matches_central_differences() -> None:
    seg, _ = _segment_with_neighbours()
    m, flows = 0.02, {"x1": 0.015, "p3": 0.006}
    st_out = _state(m, P=0.99e5)
    _, jac = seg.residuals(_state(m, P=1.0e5), st_out, flows=flows)
    row = jac[0]

    def r(m_=m, fl=flows, T=300.0, P=1.0e5):
        return seg.residuals(_state(m_, T=T, P=P), st_out, flows=fl)[0][0]

    e = 1e-7
    assert row["x2.m_dot"] == pytest.approx((r(m_=m + e) - r(m_=m - e)) / (2 * e), rel=1e-6)
    for src in ("x1", "p3"):
        up, dn = dict(flows, **{src: flows[src] + e}), dict(flows, **{src: flows[src] - e})
        assert row[f"{src}.m_dot"] == pytest.approx((r(fl=up) - r(fl=dn)) / (2 * e), rel=1e-6)


def _wall_coupled_chain(n_rows: int = 3) -> NetworkSolver:
    """A 3-row array whose targets are heated by a hot-gas duct, one wall per row."""
    span = 12 * 4.0 * D
    g = FlowNetwork()

    def pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
        b = PressureBoundary(name)
        b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
        return b

    g.add_node(pb("supply", 1.03e5, 300.0))
    g.add_node(pb("exit", 1.0e5, 300.0))
    g.add_node(pb("hot_in", 1.02e5, 900.0))
    g.add_node(pb("hot_out", 1.0e5, 900.0))
    for i in range(1, n_rows + 1):
        g.add_node(PlenumNode(f"c{i}"))
        if i < n_rows:
            g.add_node(PlenumNode(f"h{i}"))
    for i in range(1, n_rows + 1):
        g.add_element(
            ImpingementPlateElement(
                f"p{i}",
                "supply",
                f"c{i}",
                d_jet=D,
                xn_d=5.0,
                yn_d=4.0,
                z_d=1.0,
                span=span,
                plate_thickness=D,
            )
        )
        to = f"c{i + 1}" if i < n_rows else "exit"
        g.add_element(
            ImpingementCrossflowElement(f"x{i}", f"c{i}", to, length=5 * D, height=D, span=span)
        )
        h_from = "hot_in" if i == 1 else f"h{i - 1}"
        h_to = f"h{i}" if i < n_rows else "hot_out"
        g.add_element(
            ChannelElement(
                f"g{i}",
                h_from,
                h_to,
                length=5 * D,
                diameter=0.02,
                roughness=0.0,
                surface=ConvectiveSurface(area=5 * D * span),
            )
        )
        g.add_wall(
            ThermalWall(
                id=f"w{i}",
                element_a=f"g{i}",
                element_b=f"p{i}",
                layers=[WallLayer(thickness=0.001, conductivity=20.0)],
            )
        )
    return NetworkSolver(g)


def test_the_wall_coupled_chain_solves_and_heats_the_crossflow() -> None:
    solver = _wall_coupled_chain()
    res = solver.solve()
    assert res["__success__"], res.get("__message__")
    # Every crossflow node is heated; row 2 mixes row 1's spent air with its
    # own, so the temperatures need not rise monotonically.
    T = [solver._derived_states[f"c{i}"][0] for i in (1, 2, 3)]
    assert all(t > 300.0 for t in T)
    d = res["__element_diag__"]
    assert d["p1"]["Gc_Gj"] == 0.0 < d["p2"]["Gc_Gj"] < d["p3"]["Gc_Gj"]
    for i in (1, 2, 3):
        assert not validate_diagnostics("ImpingementPlateElement", d[f"p{i}"])
        assert not validate_diagnostics("ImpingementCrossflowElement", d[f"x{i}"])


def test_global_jacobian_carries_the_neighbour_flow_sensitivity() -> None:
    """At the solved point, the analytic Jacobian against central differences.

    Two claims, kept apart:

    * The whole matrix meets the standard the existing wall-coupled check
      holds (test_wall_coupling_integration: max 10%, mean 2% over
      significant entries). Entries a few percent off exist BEFORE this
      element, on the hot channels' own columns too -- the wall relay's
      existing approximations, untouched here.
    * The crossflow-source columns are where dh_dsources lives. Wherever the
      relay changes an entry it must move it toward the finite difference,
      and an entry made of nothing BUT the relay (zero without it) must
      match tightly. Disabling the relay fails this (falsified, see the PR).
    """
    solver = _wall_coupled_chain()
    res = solver.solve()
    assert res["__success__"], res.get("__message__")
    x = np.array(res["__x_solution__"], dtype=float)
    step = np.maximum(np.abs(x) * 1e-7, 1e-10)
    fd = approx_derivative(lambda v: solver._residuals(v), x, method="3-point", abs_step=step)
    with_relay = solver._residuals_and_jacobian(x)[1].toarray()
    solver._relay_flow_inputs = lambda *a, **k: None
    without = solver._residuals_and_jacobian(x)[1].toarray()

    sig = np.abs(fd) > 1e-3
    rel = np.abs(with_relay - fd)[sig] / np.abs(fd[sig])
    assert rel.max() < 0.10, rel.max()
    assert rel.mean() < 0.02, rel.mean()

    pure = 0
    for j, name in enumerate(solver.unknown_names):
        if name not in ("x1.m_dot", "x2.m_dot"):
            continue
        for i in np.where(np.abs(with_relay[:, j] - without[:, j]) > 1e-6 * np.abs(fd[:, j]))[0]:
            err_with = abs(with_relay[i, j] - fd[i, j])
            err_without = abs(without[i, j] - fd[i, j])
            assert err_with < 0.5 * err_without, (
                name,
                i,
                with_relay[i, j],
                without[i, j],
                fd[i, j],
            )
            if without[i, j] == 0.0:
                pure += 1
                assert with_relay[i, j] == pytest.approx(fd[i, j], rel=1e-3)
    assert pure > 0


def _hand_wired(extra) -> FlowNetwork:
    """Two rows, channel exit, plus whatever ``extra`` adds."""
    g = FlowNetwork()
    for n, pt in (("supply", 1.03e5), ("exit", 1.0e5)):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)
    for n in ("c1", "c2"):
        g.add_node(PlenumNode(n))
    span = 12 * 4.0 * D
    for i, to in ((1, "c2"), (2, "exit")):
        g.add_element(
            ImpingementPlateElement(
                f"p{i}",
                "supply",
                f"c{i}",
                d_jet=D,
                xn_d=5.0,
                yn_d=4.0,
                z_d=2.0,
                span=span,
                plate_thickness=D,
            )
        )
        g.add_element(
            ImpingementCrossflowElement(f"x{i}", f"c{i}", to, length=5 * D, height=2 * D, span=span)
        )
    extra(g, span)
    return g


def test_the_channel_exit_configuration_is_accepted() -> None:
    res = NetworkSolver(_hand_wired(lambda g, span: None)).solve()
    assert res["__success__"], res.get("__message__")


def test_a_bypass_crossflow_is_refused() -> None:
    """External crossflow into the chain is #467's configuration, not this one."""

    def bypass(g, span):
        b = PressureBoundary("bypass")
        b.Pt, b.Tt, b.Y = 1.03e5, 300.0, Y_AIR
        g.add_node(b)
        g.add_element(ChannelElement("byp", "bypass", "c1", length=0.01, diameter=0.01))

    with pytest.raises(ValueError, match="#467"):
        NetworkSolver(_hand_wired(bypass)).solve()


def test_spent_air_through_the_target_is_refused() -> None:
    """A gap that drains through effusion holes is #468's configuration."""
    g = FlowNetwork()
    for n, pt in (("supply", 1.03e5), ("gas", 1.0e5)):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)
    g.add_node(PlenumNode("gap"))
    g.add_element(
        ImpingementPlateElement(
            "p",
            "supply",
            "gap",
            d_jet=D,
            xn_d=5.0,
            yn_d=4.0,
            z_d=2.0,
            span=0.1,
            plate_thickness=D,
        )
    )
    g.add_element(
        EffusionPlateElement(
            "eff",
            "gap",
            "gas",
            hole_diameter=D,
            wall_thickness=2 * D,
            pitch=8 * D,
            panel_area=0.01,
        )
    )
    with pytest.raises(ValueError, match="#468"):
        NetworkSolver(g).solve()


def test_two_plates_on_one_node_are_refused() -> None:
    def twin(g, span):
        g.add_element(
            ImpingementPlateElement(
                "p1b",
                "supply",
                "c1",
                d_jet=D,
                xn_d=5.0,
                yn_d=4.0,
                z_d=2.0,
                span=span,
                plate_thickness=D,
            )
        )

    with pytest.raises(ValueError, match="more than one plate"):
        NetworkSolver(_hand_wired(twin)).solve()
