"""Mass-flow derivatives of a non-circular channel (#463), and the relay
through a wall that feeds a node's own temperature back to it.

channel_smooth took its mass flow on pi D^2/4; a channel with a separate
hydraulic diameter has another flow area, so dh/dmdot, ddP/dmdot and
dT_aw/dmdot were off by (Dh/D)^2. Values were never affected.

A wall that couples a node's outflow to its inflow (two walls in series)
feeds the node's temperature back into its own heat input: T' = a + b T'.
The relay read that node's half-built sensitivities and dropped the loop --
19% of dT/d(cooling flow).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import test_solver_robustness_481 as tsr

import combaero as cb
from combaero.network import NetworkSolver

AIR = np.array(cb.species.dry_air())


@pytest.mark.parametrize("corr", ["gnielinski", "dittus_boelter"])
def test_channel_smooth_takes_the_flow_area_for_its_mass_flow_derivatives(corr: str) -> None:
    """A flat duct: Dh = 12 mm, flow area 4x that of a 12 mm circle."""
    T, P, Dh, L = 500.0, 2e5, 0.012, 0.5
    A = 4.0 * math.pi * Dh * Dh / 4.0
    rho = cb.density(T, P, list(AIR))
    v = 30.0

    def run(vel: float) -> cb.ChannelResult:
        return cb.channel_smooth(T, P, AIR, vel, Dh, L, 300.0, corr, True, 1.0, 0.0, 1.0, 1.0, A)

    r = run(v)
    hv = 1e-6 * v
    dm = 2 * hv * rho * A  # d(mdot) = rho A dv on THIS duct's area
    assert r.dh_dmdot == pytest.approx((run(v + hv).h - run(v - hv).h) / dm, rel=1e-5)
    assert r.ddP_dmdot == pytest.approx((run(v + hv).dP - run(v - hv).dP) / dm, rel=1e-5)
    assert r.dT_aw_dmdot == pytest.approx((run(v + hv).T_aw - run(v - hv).T_aw) / dm, rel=1e-5)
    # The default is the circular area; values do not depend on it.
    circ = cb.channel_smooth(T, P, AIR, v, Dh, L, 300.0, corr, True, 1.0, 0.0, 1.0, 1.0)
    assert circ.h == r.h and circ.dP == r.dP
    assert circ.dh_dmdot == pytest.approx(4.0 * r.dh_dmdot, rel=1e-12)


def _relay_vs_fd(s: NetworkSolver, x: np.ndarray, var: str) -> tuple[dict, dict]:
    j = s.unknown_names.index(var)
    relay = s._propagate_states(x)
    an = {n: relay[n].get(j, {}).get("T", 0.0) for n in ("n2", "n3", "cp")}
    h = abs(x[j]) * 1e-5
    T = []
    for sgn in (1.0, -1.0):
        xp = x.copy()
        xp[j] += sgn * h
        s._propagate_states(xp)
        T.append({n: s._derived_states[n][0] for n in an})
    return an, {n: (T[0][n] - T[1][n]) / (2 * h) for n in an}


@pytest.mark.parametrize("Dh", [None, 0.012], ids=["circular", "flat"])
@pytest.mark.parametrize("var", ["cc.m_dot", "e1.m_dot", "e2.m_dot"])
def test_the_relay_closes_a_walls_self_loop(Dh: float | None, var: str) -> None:
    g = tsr._series_walls()
    if Dh is not None:
        for eid in ("e1", "e2", "cc"):
            g.elements[eid].Dh = Dh
    s = NetworkSolver(g)
    x = np.array(s.solve()["__x_solution__"])
    an, fd = _relay_vs_fd(s, x, var)
    for n in an:
        assert an[n] == pytest.approx(fd[n], rel=1e-5, abs=1e-6 * max(abs(v) for v in fd.values()))
