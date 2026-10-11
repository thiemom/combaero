"""Bleed through a duct wall: the station is Bassett's separating straight run (#471).

Bassett, Winterbone & Pearson (2001), separating flow, straight-through loss
K2/K5 (Eq. 15): K = q^2 - 1.5 q + 0.5, q = m_out/m_in, theta- and
psi-independent. The station momentum balance with the bleed carrying
kappa = 0.75 of the arriving axial velocity reproduces it for every split, so
the closure is Bassett's and is scored here against his own measured K5
(Fig. 7c, theta 45/60/90, psi 1/3).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

DATA = Path(__file__).resolve().parents[2] / "validation" / "junction" / "data" / "bassett2001"
X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)


def _k5(q: float) -> float:
    """Bassett Eq. (15), written out here so the check is independent of C++."""
    return q * q - 1.5 * q + 0.5


def _station_K(q: float, kappa: float) -> float:
    """Total-pressure loss coefficient of one whole station from C++, on the
    arriving dynamic head, at constant density."""
    P, T, A, m_a = 1.2e5, 500.0, 4e-3, 0.4
    rho = cb.density(T, P, X_AIR)
    half = cb.station_half_drop(m_a, q * m_a, P, T, X_AIR, A, kappa)
    ua, ub = m_a / (rho * A), q * m_a / (rho * A)
    return (2.0 * half.dP + 0.5 * rho * (ua**2 - ub**2)) / (0.5 * rho * ua**2)


@pytest.mark.parametrize("q", [0.0, 0.3, 0.6, 0.9, 0.97, 1.0])
def test_the_bleed_station_is_bassett_k5(q: float) -> None:
    assert _station_K(q, cb.STATION_KAPPA_BLEED_BASSETT) == pytest.approx(_k5(q), abs=1e-12)


def test_normal_injection_is_not_bassett() -> None:
    """kappa = 0 (the impingement merge) is a different closure: it would
    over-state the bleed's regain."""
    assert abs(_station_K(0.9, cb.STATION_KAPPA_MERGE_NORMAL) - _k5(0.9)) > 0.1


def _fig7c() -> list[tuple[float, float]]:
    pts = []
    for f in sorted(DATA.glob("bassett_fig07c_K5_*_measured.csv")):
        for line in f.read_text().splitlines()[1:]:
            q, k = map(float, line.split(","))
            pts.append((q, k))
    return pts


def test_scored_against_bassett_fig7c() -> None:
    """Fidelity on the author's own data, all four series (theta 45/60/90,
    psi 1/3: theta- and psi-independent as Bassett states). Absolute error in
    K, since K crosses zero. The bleed-relevant end is q >= 0.6 (an effusion
    station bleeds a few percent)."""
    pts = _fig7c()
    assert len(pts) == 43
    err = np.array([_station_K(q, cb.STATION_KAPPA_BLEED_BASSETT) - k for q, k in pts])
    hi = np.array([_station_K(q, cb.STATION_KAPPA_BLEED_BASSETT) - k for q, k in pts if q >= 0.6])
    print(
        f"\nall {len(err)}: bias {err.mean():+.4f}, MAE {np.abs(err).mean():.4f}, "
        f"max {np.abs(err).max():.4f}"
        f"\nq>=0.6 {len(hi)}: bias {hi.mean():+.4f}, MAE {np.abs(hi).mean():.4f}, "
        f"max {np.abs(hi).max():.4f}"
    )
    # Measured 2026-10-08: all 43 bias +0.017, MAE 0.039, max 0.158 (q ~ 0,
    # full diversion, formula 0.49 vs 0.35 measured); q >= 0.6: bias -0.019,
    # MAE 0.021, max 0.046. Reported, not tuned. Near q -> 1 the misses sit at
    # q -> 1, where the data carry a +0.04 offset the formula (K5(1) = 0)
    # cannot have: a measurement datum, most likely the rig's own friction.
    assert np.abs(err).mean() < 0.045
    assert np.abs(hi).mean() < 0.025


def test_a_solved_bleed_network_satisfies_bassett() -> None:
    """Reservoir -> duct -> bleed station -> duct -> outlet, near-frictionless
    and nearly isothermal: the outlet static must sit below the reservoir by
    exactly the entry dynamic head plus Bassett's station drop at the SOLVED
    split -- computed here from K5 alone, independently of the C++ station."""
    A = 1e-3
    g = FlowNetwork()
    for n, pt in (("in", 1.010e5), ("out", 1.000e5), ("lat", 0.990e5)):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)
    g.add_node(PlenumNode("b"))
    kappa = cb.STATION_KAPPA_BLEED_BASSETT
    g.add_element(
        CrossflowSegmentElement(
            "s0", "in", "b", length=1e-6, area=A, entry_K=0.0, to_kappa=kappa, next_seg="s1"
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s1", "b", "out", length=1e-6, area=A, from_kappa=kappa, prev_seg="s0"
        )
    )
    g.add_element(OrificeElement("bleed", "b", "lat", diameter=0.01, Cd=0.6, correlation="fixed"))
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    d = r["__element_diag__"]
    m_a, m_b = d["s0"]["m_dot"], d["s1"]["m_dot"]
    q = m_b / m_a
    assert 0.5 < q < 0.99
    rho = cb.density(300.0, r["b.P"], X_AIR)
    qa, qb = 0.5 * m_a**2 / (rho * A * A), 0.5 * m_b**2 / (rho * A * A)
    # Static: P_a = Pt_in - qa (entry, K_in = 0); station Pt loss K5 qa.
    P_a = 1.010e5 - qa
    P_b = P_a - (_k5(q) * qa - (qa - qb))
    # The outlet dumps (static continuity).
    assert P_b == pytest.approx(1.000e5, abs=2e-3 * (qa + qb))
    assert math.isfinite(r["b.P"])
