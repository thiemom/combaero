"""Fanno flow integrated in Mach number (compressible-channel redesign).

At fixed mass flux and stagnation temperature the station state is algebraic
in M, so the duct length is a smooth integral of dx/dM, which vanishes at
sonic. Checked against the x-march (a different integration of the same
equations), the exact sonic point, the choke root, and the smoothness the
network Jacobians will rely on.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb
from combaero import _core

X = cb.species.dry_air()
D, L, ROUGH = 0.02, 1.0, 1e-5
PT, TT = 2.0e5, 300.0


def test_sonic_is_exactly_mach_one() -> None:
    """dx/dM > 0 below sonic and 0 at M = 1 for the thermally perfect gas."""
    for Tt in (300.0, 1500.0):
        for M in (0.2, 0.6, 0.9, 0.99, 0.9999):
            assert _core.fanno_dx_dmach(300.0, Tt, M, X, D, ROUGH, "haaland") > 0.0
        at_sonic = _core.fanno_dx_dmach(300.0, Tt, 1.0, X, D, ROUGH, "haaland")
        near = _core.fanno_dx_dmach(300.0, Tt, 0.999, X, D, ROUGH, "haaland")
        assert abs(at_sonic) < 1e-12 * near


@pytest.mark.parametrize(
    ("Tt", "model", "fm"),
    [(300.0, "haaland", 1.0), (1500.0, "haaland", 1.0), (300.0, "fixed", 0.02)],
)
@pytest.mark.parametrize("frac", [0.2, 0.6, 0.9, 0.99])
def test_the_exit_state_is_the_x_marchs(Tt: float, model: str, fm: float, frac: float) -> None:
    """Same equations, integrated in x at 20000 steps: 1e-9 at constant f;
    5e-8 once f depends on Re, the floor kinks in mu(T) set."""
    G = frac * _core.fanno_choked_mass_flux(PT, Tt, X, L, D, ROUGH, model, fm)
    r = _core.fanno_duct(PT, Tt, G, X, L, D, ROUGH, model, fm)
    assert not r.choked and r.L_star > L
    if model == "fixed":
        m = cb.fanno_channel(r.inlet.T, r.inlet.P, r.inlet.u, L, D, fm, X, 20000)
    else:
        m = cb.fanno_channel_rough(
            r.inlet.T, r.inlet.P, r.inlet.u, L, D, ROUGH, X, model, fm, 20000
        )
    tol = 1e-9 if model == "fixed" else 5e-8
    assert pytest.approx(m.outlet.P, rel=tol) == r.exit.P
    assert pytest.approx(m.outlet.T, rel=tol) == r.exit.T


def test_the_inlet_is_isentropic_from_the_stagnation_state() -> None:
    G = 250.0
    r = _core.fanno_duct(PT, TT, G, X, L, D, ROUGH, "haaland")
    s = r.inlet
    assert s.rho * s.u == pytest.approx(G, rel=1e-12)
    assert cb.h_mass(s.T, X) + 0.5 * s.u**2 == pytest.approx(cb.h_mass(TT, X), rel=1e-12)
    assert cb.P0_from_static(s.P, s.T, s.M, X) == pytest.approx(PT, rel=1e-9)


def test_the_choked_flux_chokes_exactly_at_the_exit() -> None:
    G_ch = _core.fanno_choked_mass_flux(PT, TT, X, L, D, ROUGH, "haaland")
    r = _core.fanno_duct(PT, TT, G_ch, X, L, D, ROUGH, "haaland")
    assert r.L_star == pytest.approx(L, rel=1e-9)
    assert pytest.approx(1.0, abs=1e-5) == r.exit.M
    assert not _core.fanno_duct(PT, TT, G_ch * (1 - 1e-6), X, L, D, ROUGH, "haaland").choked
    assert _core.fanno_duct(PT, TT, G_ch * (1 + 1e-6), X, L, D, ROUGH, "haaland").choked
    # The x-march's own choke flag, bisected (to the 6 digits recorded),
    # lands on the same flow.
    area = math.pi * D * D / 4
    assert G_ch * area == pytest.approx(0.113271, rel=1e-5)
    assert G_ch < _core.fanno_sonic_mass_flux(PT, TT, X)


def test_the_exit_pressure_falls_monotonically_with_the_flux() -> None:
    G_ch = _core.fanno_choked_mass_flux(PT, TT, X, L, D, ROUGH, "haaland")
    P = [
        _core.fanno_duct(PT, TT, f * G_ch, X, L, D, ROUGH, "haaland").exit.P
        for f in np.linspace(0.05, 0.999, 40)
    ]
    assert all(b < a for a, b in zip(P, P[1:], strict=False))


def test_the_choking_length_is_smooth_in_the_flux() -> None:
    """A fixed quadrature rule: L*(G) differences cleanly at two step sizes
    (an adaptive rule switches its panels and does not)."""
    for G in (100.0, 300.0, 400.0):

        def Ls(g: float) -> float:
            M = _core.fanno_inlet_mach(PT, TT, g, X)
            return _core.fanno_length_between(g, TT, M, 1.0, X, D, ROUGH, "haaland")

        d = [(Ls(G * (1 + h)) - Ls(G * (1 - h))) / (2 * G * h) for h in (1e-4, 1e-5)]
        assert d[0] < 0.0
        assert d[1] == pytest.approx(d[0], rel=1e-5)


def test_an_inlet_at_or_above_the_sonic_flux_is_choked_at_the_inlet() -> None:
    G_max = _core.fanno_sonic_mass_flux(PT, TT, X)
    assert _core.fanno_inlet_mach(PT, TT, G_max * 1.01, X) == 1.0
    r = _core.fanno_duct(PT, TT, G_max * 1.01, X, L, D, ROUGH, "haaland")
    assert r.choked and r.L_star == 0.0


def test_invalid_ducts_are_refused() -> None:
    with pytest.raises(ValueError):
        _core.fanno_duct(PT, TT, 100.0, X, 0.0, D, ROUGH, "haaland")
    with pytest.raises(ValueError):
        _core.fanno_length_between(100.0, TT, 0.5, 0.4, X, D, ROUGH, "haaland")
