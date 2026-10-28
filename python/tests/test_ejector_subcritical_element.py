"""The compressible EjectorElement solves the choked SUBCRITICAL regime.

Between the critical back pressure P_c* and the dead-head P_b0 a choked
ejector's entrainment droops (Phase 1A: the Tier-1 chord from (P_c*,
omega_crit) to (P_b0, 0), OPERATING_REGIMES_DESIGN.md sec 8.1) and its
outlet floats. The element used to keep omega_crit there and pin the outlet
port at recovery * P_c*, below the outlet boundary, so a choked ejector with
its outlet in that window had no root (every subcritical stress draw failed,
#619). Fidelity is judged against the Python reference closures at the
solved port states; the droop's SHAPE is the closure's, validated on Henry
et al. 2007 (test_ejector_subcritical.py), not re-scored here. Each anchor
closure is run inside its own window.
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
from combaero.network._ejector_huang1999 import (
    EjectorGeometry,
    critical_back_pressure,
    dead_head_back_pressure,
    entrainment_ratio,
    huang_critical_back_pressure,
    subcritical_entrainment_ratio,
    sudden_expansion_dead_head,
)
from combaero.network.ejector_element import EjectorElement, choke_plane_gamma

Y = list(cb.species.dry_air_mass())
D = 0.05
A = math.pi * D * D / 4
# The reference double-choked ejector. Its (P_c*, P_b0): kracik_dvorak
# ~ (58, 176) kPa, huang1999 (the default) ~ (37, 59) kPa.
A_T, A_E, A_MIX = 3.14e-5, 1.0e-4, 8.0e-4


def _solve(p_out: float, closure: str) -> tuple[dict, NetworkSolver]:
    g = FlowNetwork()
    for nid, Pt, Tt in (("bp", 7.0e5, 440.0), ("bs", 22850.0, 280.0), ("bo", p_out, 300.0)):
        g.add_node(PressureBoundary(nid, Pt=Pt, Tt=Tt, Y=Y))
    for nid in ("mp", "ms", "mo"):
        g.add_node(MomentumChamberNode(nid, area=A))
    for eid, a, b in (("cp", "bp", "mp"), ("cs", "bs", "ms"), ("co", "mo", "bo")):
        g.add_element(ChannelElement(eid, a, b, length=0.2, diameter=D, regime="incompressible"))
    g.add_element(
        EjectorElement(
            id="ej",
            primary_node="mp",
            secondary_node="ms",
            outlet_node="mo",
            throat_area=A_T,
            nozzle_exit_area=A_E,
            mixing_area=A_MIX,
            regime="compressible",
            anchor_closure=closure,
        )
    )
    s = NetworkSolver(g)
    return s.solve(timeout=60.0), s


def _reference(r: dict, closure: str) -> tuple[float, float, float, float]:
    """(omega_crit, P_c*, P_b0, omega_chord) at the solved port states."""
    pg, tg, pe, te, po = r["mp.Pt"], r["mp.Tt"], r["ms.Pt"], r["ms.Tt"], r["mo.Pt"]
    X = [r[f"ms.X[{i}]"] for i in range(len(cb.species.dry_air()))]
    gamma, r_gas = choke_plane_gamma(te, X), float(cb.specific_gas_constant(X))
    geom = EjectorGeometry(A_E / A_T, A_MIX / A_T)
    oc = entrainment_ratio(pg, tg, pe, te, geom, gamma).omega
    if closure == "huang1999":
        pc = huang_critical_back_pressure(pg, tg, pe, te, geom, gamma, r_gas).p_c_pa
        pb0 = sudden_expansion_dead_head(pg, tg, pe, te, geom, gamma, r_gas).p_b0_pa
    else:
        pc = critical_back_pressure(pg, tg, pe, te, geom, gamma, r_gas, 1.0).p_c_pa
        pb0 = dead_head_back_pressure(pg, tg, pe, te, geom, gamma, r_gas, 1.0).p_b0_pa
    return oc, pc, pb0, subcritical_entrainment_ratio(oc, pc, pb0, po)


DROOP = [("kracik_dvorak", p) for p in (8.0e4, 1.0e5, 1.3e5, 1.6e5)] + [
    ("huang1999", p) for p in (4.0e4, 4.5e4, 5.0e4, 5.5e4)
]
CRITICAL = [("kracik_dvorak", p) for p in (3.0e4, 5.0e4)] + [
    ("huang1999", p) for p in (2.0e4, 3.0e4)
]
SWEEP = {
    "kracik_dvorak": (3.0e4, 6.0e4, 1.0e5, 1.4e5, 1.7e5),
    "huang1999": (1.5e4, 3.0e4, 4.5e4, 5.3e4, 5.8e4),
}


@pytest.mark.parametrize(("closure", "p_out"), DROOP)
def test_the_subcritical_entrainment_is_the_droop(closure: str, p_out: float) -> None:
    r, _ = _solve(p_out, closure)
    assert r["__success__"]
    oc, pc, pb0, chord = _reference(r, closure)
    assert 1.02 * pc < r["mo.Pt"] < 0.98 * pb0, "not in the subcritical window"
    assert r["ej.critical_mode"] == 0.0
    omega = r["cs.m_dot"] / r["cp.m_dot"]
    assert 0.0 < omega < oc
    assert omega == pytest.approx(chord, rel=1e-8)


@pytest.mark.parametrize(("closure", "p_out"), CRITICAL)
def test_critical_mode_keeps_the_plateau_exactly(closure: str, p_out: float) -> None:
    """Below P_c* the droop is switched out exactly (no smooth-min shave)."""
    r, _ = _solve(p_out, closure)
    assert r["__success__"]
    oc, pc, _pb0, _ = _reference(r, closure)
    assert r["mo.Pt"] < 0.98 * pc
    assert r["ej.critical_mode"] == 1.0
    assert r["cs.m_dot"] / r["cp.m_dot"] == pytest.approx(oc, rel=1e-10)


@pytest.mark.parametrize("closure", list(SWEEP))
def test_entrainment_falls_monotonically_with_back_pressure(closure: str) -> None:
    omegas = []
    for p_out in SWEEP[closure]:
        r, _ = _solve(p_out, closure)
        assert r["__success__"]
        omegas.append(r["cs.m_dot"] / r["cp.m_dot"])
    assert all(b < a for a, b in zip(omegas[1:], omegas[2:], strict=False))
    assert omegas[0] == pytest.approx(omegas[1], rel=0.05)  # both below P_c*: ~plateau
    assert omegas[-1] < 0.1 * omegas[0]  # near dead-head
