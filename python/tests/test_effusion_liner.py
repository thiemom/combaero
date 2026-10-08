"""EffusionLiner: a duct-fed effusion wall, N stations along the backside (#471).

Each station is a Bassett bleed (the segments either side carry its momentum)
feeding a McGreehan-Schotsch panel whose Cd sees the duct's crossflow,
U1/Vi. The station closure is scored in test_bleed_station_bassett.py and the
crossflow term on Rohde in the orifice harness; here: the assembly, the
coupling between them, and the Jacobian.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.optimize._numdiff import approx_derivative

import combaero as cb
from combaero.network import (
    ChannelElement,
    CrossflowSegmentElement,
    EffusionLiner,
    EffusionPlateElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkMixtureState,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

Y_AIR = cb.mole_to_mass(cb.species.dry_air())
X_AIR = cb.species.dry_air()


def _pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _liner(n: int = 4, **kw) -> EffusionLiner:
    args = {
        "length": 0.152,
        "width": 0.152,
        "duct_height": 0.03,
        "hole_diameter": 3.27e-3,
        "wall_thickness": 6.3e-3,
        "pitch_x": 15.24e-3,
        "pitch_y": 15.24e-3,
    }
    args.update(kw)
    return EffusionLiner("ln", n_segments=n, **args)


def _net(liner: EffusionLiner, p_in=1.04e5, p_out=1.035e5, chamber: bool = False):
    g = FlowNetwork()
    g.add_node(_pb("cin", p_in, 600.0))
    g.add_node(_pb("cout", p_out, 600.0))
    if chamber:
        g.add_node(_pb("gas_in", 1.003e5, 1500.0))
        g.add_node(_pb("exit", 1.0e5, 1500.0))
        g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
        g.add_element(ChannelElement("duct", "gas_in", "ch", length=0.2, diameter=0.08))
        g.add_element(ChannelElement("tail", "ch", "exit", length=0.2, diameter=0.08))
        liner.add_to(g, "cin", "cout", "ch")
    else:
        g.add_node(_pb("gas", 1.0e5, 1500.0))
        liner.add_to(g, "cin", "cout", "gas")
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


def test_the_stations_are_wired_as_bleeds_with_duct_fed_panels() -> None:
    liner = _liner(3)
    g = _net(liner)
    kappa = cb.STATION_KAPPA_BLEED_BASSETT
    s0 = g.elements["ln__s0"]
    assert isinstance(s0, CrossflowSegmentElement) and s0.entry_K == 0.5
    assert (s0.to_kappa, s0.next_seg) == (kappa, "ln__s1")
    last = g.elements["ln__s3"]
    assert (last.from_kappa, last.prev_seg, last.to_kappa) == (kappa, "ln__s2", None)
    assert last.to_node == "cout"
    for i in (1, 2, 3):
        p = g.elements[f"ln__p{i}"]
        assert isinstance(p, EffusionPlateElement)
        assert isinstance(g.nodes[f"ln__b{i}"], PlenumNode)
        assert p.correlation == "McGreehanSchotsch"
        assert p.crossflow_segments == (f"ln__s{i - 1}", f"ln__s{i}")
        assert p.crossflow_area == pytest.approx(0.03 * 0.152)


def test_each_panel_cd_sees_its_stations_crossflow() -> None:
    """The panel's Cd is McGreehan-Schotsch at the U1/Vi of the SOLVED duct
    flows either side of its station, Vi from the station's static drive."""
    s, r = _solve(_net(_liner()))
    d = r["__element_diag__"]
    for i in (1, 2, 3, 4):
        m_a, m_b = d[f"ln__s{i - 1}"]["m_dot"], d[f"ln__s{i}"]["m_dot"]
        P, T = r[f"ln__b{i}.P"], s._derived_states[f"ln__b{i}"][0]
        ratio = cb.crossflow_velocity_ratio(m_a, m_b, P, T, X_AIR, 1.0e5, 0.03 * 0.152)
        assert d[f"ln__p{i}"]["U1_over_Vi"] == pytest.approx(ratio.U1_over_Vi, rel=1e-9)
        panel = s.network.elements[f"ln__p{i}"]
        st = NetworkMixtureState(T=T, P=P, Pt=P, Tt=T, m_dot=d[f"ln__p{i}"]["m_dot"], Y=Y_AIR)
        hole = cb.DischargeHoleGeometry(d=3.27e-3, L=panel.hole_length, r=0.0)
        cd = cb.discharge_cd(
            cb.DischargeCdCorrelation.McGreehanSchotsch1988,
            hole,
            cb.DischargeHoleState(Re=panel._hole_reynolds(st), U1_over_Vi=ratio.U1_over_Vi),
        )
        assert d[f"ln__p{i}"]["Cd"] == pytest.approx(cd, rel=1e-9)


def test_the_bleed_regains_static_pressure_along_the_duct() -> None:
    """A distributing duct: each Bassett station recovers static pressure as
    the duct slows, so the downstream panels see more drive and pass more."""
    _, r = _solve(_net(_liner()))
    d = r["__element_diag__"]
    P = [r[f"ln__b{i}.P"] for i in (1, 2, 3, 4)]
    bleed = [d[f"ln__p{i}"]["m_dot"] for i in (1, 2, 3, 4)]
    assert all(b > a for a, b in zip(P, P[1:], strict=False))
    assert all(b > a for a, b in zip(bleed, bleed[1:], strict=False))


def test_crossflow_lowers_the_cd() -> None:
    """Upstream panels sit in faster crossflow and so have the lower Cd."""
    _, r = _solve(_net(_liner()))
    d = r["__element_diag__"]
    u = [d[f"ln__p{i}"]["U1_over_Vi"] for i in (1, 2, 3, 4)]
    cd = [d[f"ln__p{i}"]["Cd"] for i in (1, 2, 3, 4)]
    assert all(b < a for a, b in zip(u, u[1:], strict=False))
    assert all(b > a for a, b in zip(cd, cd[1:], strict=False))


def test_the_summary_is_the_stations_own_numbers() -> None:
    liner = _liner()
    _, r = _solve(_net(liner))
    d = r["__element_diag__"]
    s = liner.summarize(d)
    assert s["m_dot"] == pytest.approx(d["ln__s0"]["m_dot"])
    assert s["m_dot"] == pytest.approx(s["m_dot_out"] + s["m_dot_bleed"], rel=1e-8)
    assert s["rows_m_dot_bleed"] == [d[f"ln__p{i}"]["m_dot"] for i in (1, 2, 3, 4)]
    hole_area = math.pi * 3.27e-3**2 / 4.0
    assert s["station_psi_min"] == pytest.approx(
        0.03 * 0.152 / (d["ln__p1"]["n_holes"] * hole_area)
    )
    assert s["station_psi_beyond_bassett"] == 1.0
    assert s["coolant_ht_plenum_assumed"] == 1.0


def test_into_a_merge_chamber_it_solves_with_a_chamber_gas_side() -> None:
    liner = _liner()
    _, r = _solve(_net(liner, chamber=True))
    s = liner.summarize(r["__element_diag__"])
    assert len(s["rows_T_wall_hot"]) == 4
    assert all(600.0 < t < 1500.0 for t in s["rows_T_wall_hot"])


def test_global_jacobian_matches_finite_differences() -> None:
    """Above McGreehan-Schotsch's Re floor (1e4), where the correlation's own
    derivative is the value's: exact, including the crossflow chain into the
    neighbour segments and the stations' momentum."""
    s, r = _solve(_net(_liner(3), p_in=1.12e5, p_out=1.11e5))
    d = r["__element_diag__"]
    for i in (1, 2, 3):
        T, P = s._derived_states[f"ln__b{i}"][0], r[f"ln__b{i}.P"]
        st = NetworkMixtureState(T=T, P=P, Pt=P, Tt=T, m_dot=d[f"ln__p{i}"]["m_dot"], Y=Y_AIR)
        assert s.network.elements[f"ln__p{i}"]._hole_reynolds(st) > 1.0e4
    x = np.array(r["__x_solution__"], dtype=float)
    _, jac = s._residuals_and_jacobian(x)
    fd = approx_derivative(
        lambda v: s._residuals(v), x, method="3-point", abs_step=np.maximum(np.abs(x) * 1e-6, 1e-9)
    )
    mask = np.abs(fd) > 1e-3
    assert (np.abs(jac.toarray() - fd)[mask] / np.abs(fd[mask])).max() < 1e-3


def test_a_duct_fed_panel_must_use_the_crossflow_correlation() -> None:
    g = FlowNetwork()
    g.add_node(PlenumNode("b"))
    g.add_node(_pb("gas", 1.0e5, 1500.0))
    p = EffusionPlateElement(
        "p",
        "b",
        "gas",
        hole_diameter=1e-3,
        wall_thickness=2e-3,
        pitch=6e-3,
        panel_area=0.01,
        correlation="IdelchikThick",
        crossflow_segments=(None, None),
        crossflow_area=1e-3,
    )
    g.add_element(p)
    p.resolve_topology(g)
    with pytest.raises(ValueError, match="McGreehanSchotsch"):
        p.validate()
