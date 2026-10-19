"""Fanno-Rayleigh duct flow: friction AND heat transfer (#551).

The compressible channel kernel is adiabatic Fanno flow; a wall's heat only
raised the downstream node's Tt, so the momentum omitted the Rayleigh change
(heating costs total pressure, cooling gains it). fanno_rayleigh_channel_flow
marches the duct with Tt rising linearly from Tt0 to Tt_exit, exactly for a
thermally perfect gas (see fanno_mach.h).

Fidelity is checked against relations the march does not use: at zero heat
the Mach-quadrature Fanno kernel, at zero friction the algebra of Rayleigh
flow (impulse P + G u and h_t conserved), and the Rayleigh thermal choke.
"""

from __future__ import annotations

import math

import pytest
from scipy.optimize import brentq

import combaero as cb

X = list(cb.species.dry_air())
PT, TT, L, D = 2.0e5, 600.0, 0.5, 0.01
CASES = [("fixed", 0.02), ("haaland", 1.0)]


def _station(G: float, Tt: float, M: float) -> tuple[float, float]:
    s = cb.fanno_state_at_mach(G, Tt, M, X)
    return s.P, s.u


def _flow(P: float, Tte: float, et: bool, model: str, fm: float, **kw) -> cb.FannoChannelFlow:
    return cb.fanno_rayleigh_channel_flow(PT, TT, Tte, X, P, et, L, D, 0.0, model, fm, **kw)


@pytest.mark.parametrize(("model", "fm"), CASES)
@pytest.mark.parametrize("et", [False, True])
@pytest.mark.parametrize("P", [1.99e5, 1.8e5, 1.3e5, 0.5e5, PT - 1.5])
def test_without_heat_it_is_the_fanno_flow(model: str, fm: float, et: bool, P: float) -> None:
    a = cb.fanno_channel_flow(PT, TT, X, P, et, L, D, 0.0, model, fm)
    b = _flow(P, TT, et, model, fm)
    assert b.choked == a.choked
    assert pytest.approx(a.G, rel=2e-7) == b.G
    tol = 1e-6 if PT - P > 100 else 2e-4  # both kernels' differences at 1.5 Pa
    assert b.dG_dPt0 == pytest.approx(a.dG_dPt0, rel=tol)
    assert b.dG_dP_target == pytest.approx(a.dG_dP_target, rel=tol, abs=1e-12)
    # Fanno's Tt0 moves the whole duct's Tt: here that is Tt0 and Tt_exit.
    assert b.dG_dTt0 + b.dG_dTt_exit == pytest.approx(a.dG_dTt0, rel=tol)


@pytest.mark.parametrize(("G", "Tte"), [(150.0, 900.0), (150.0, 450.0), (250.0, 700.0)])
def test_without_friction_it_is_rayleigh_flow(G: float, Tte: float) -> None:
    """f = 0: the exit Mach is the one with the inlet's impulse at Tt_exit."""
    m = cb.fanno_rayleigh_march(G, PT, TT, Tte, X, L, D, 0.0, "fixed", 0.0)
    P_in, u_in = _station(G, TT, m.M_in)

    def impulse(M: float) -> float:
        P, u = _station(G, Tte, M)
        return P + G * u - (P_in + G * u_in)

    assert m.M_exit == pytest.approx(brentq(impulse, 1e-4, 1.0 - 1e-12), rel=1e-8)
    assert m.Pt_exit < PT if Tte > TT else m.Pt_exit > PT


def test_the_rayleigh_thermal_choke() -> None:
    """f = 0, heated to 900 K: the choked flux is sonic at the exit with the
    inlet's impulse."""
    Tte = 900.0
    Gc = cb.fanno_rayleigh_choked_mass_flux(PT, TT, Tte, X, L, D, 0.0, "fixed", 0.0)

    def excess(G: float) -> float:
        P_in, u_in = _station(G, TT, cb.fanno_inlet_mach(PT, TT, G, X))
        P_s, u_s = _station(G, Tte, 1.0)
        return (P_in + G * u_in) - (P_s + G * u_s)

    G_sonic = cb.fanno_sonic_mass_flux(PT, TT, X)
    assert Gc == pytest.approx(brentq(excess, 1.0, G_sonic * (1 - 1e-12)), rel=1e-7)


@pytest.mark.parametrize(("model", "fm"), CASES)
def test_heating_chokes_earlier_and_costs_pressure(model: str, fm: float) -> None:
    G = {Tte: _flow(1.6e5, Tte, False, model, fm).G for Tte in (450.0, 600.0, 900.0)}
    assert G[900.0] < G[600.0] < G[450.0]
    Gc = [
        cb.fanno_rayleigh_choked_mass_flux(PT, TT, Tte, X, L, D, 0.0, model, fm)
        for Tte in (450.0, 600.0, 900.0)
    ]
    assert Gc[2] < Gc[1] < Gc[0]


def _fd(fn, x0: float, rel: float = 1e-5) -> float:
    h = rel * abs(x0)
    return (fn(x0 + h) - fn(x0 - h)) / (2 * h)


@pytest.mark.parametrize(("model", "fm"), CASES)
@pytest.mark.parametrize(("P", "et"), [(1.6e5, False), (1.9e5, True), (1.3e5, True)])
@pytest.mark.parametrize("Tte", [400.0, 600.0, 900.0])
def test_derivatives_against_differences(model: str, fm: float, P: float, et: bool, Tte: float):
    """Unchoked, choked at the exit, and (haaland, cooled to 400 K) choked at
    an interior throat."""
    r = _flow(P, Tte, et, model, fm)

    def G(Pt0: float = PT, Tt0: float = TT, Ttx: float = Tte, Px: float = P) -> float:
        return cb.fanno_rayleigh_channel_flow(
            Pt0, Tt0, Ttx, X, Px, et, L, D, 0.0, model, fm, with_derivatives=False
        ).G

    tol = 1e-5 if not r.choked else 1e-6
    assert r.dG_dPt0 == pytest.approx(_fd(lambda v: G(Pt0=v), PT), rel=tol)
    assert r.dG_dTt0 == pytest.approx(_fd(lambda v: G(Tt0=v), TT), rel=tol)
    assert r.dG_dTt_exit == pytest.approx(_fd(lambda v: G(Ttx=v), Tte), rel=tol)
    if r.choked:
        assert r.dG_dP_target == 0.0
    else:
        assert r.dG_dP_target == pytest.approx(_fd(lambda v: G(Px=v), P), rel=tol)


def test_an_interior_throat() -> None:
    """Cooled from 600 to 400 K with Haaland friction, the choked duct is sonic
    inside and decelerates after it: the throat sets the flux."""
    Gc = cb.fanno_rayleigh_choked_mass_flux(PT, TT, 400.0, X, L, D, 0.0, "haaland", 1.0)
    m = cb.fanno_rayleigh_march(Gc, PT, TT, 400.0, X, L, D, 0.0, "haaland", 1.0)
    assert m.zeta_peak == pytest.approx(0.0, abs=1e-9)
    assert m.M_exit < 0.95


def test_composition_derivatives() -> None:
    Y = cb.mole_to_mass(X)
    k = cb.species_index_from_name("O2")
    v = [0.0] * len(Y)
    v[k] = 1.0
    r = _flow(1.6e5, 800.0, False, "haaland", 1.0, dY_directions=[v])

    def G(s: float) -> float:
        Ys = list(Y)
        Ys[k] += s
        Xs = list(cb.mass_to_mole(Ys))
        return cb.fanno_rayleigh_channel_flow(
            PT, TT, 800.0, Xs, 1.6e5, False, L, D, 0.0, "haaland", 1.0, with_derivatives=False
        ).G

    assert r.dG_ddir[0] == pytest.approx((G(1e-5) - G(-1e-5)) / 2e-5, rel=1e-4)


def test_the_warm_start_only_changes_the_cost() -> None:
    cold = _flow(1.7e5, 800.0, False, "haaland", 1.0)
    for g in (0.5, 0.99, 1.001, 1.3):
        warm = _flow(1.7e5, 800.0, False, "haaland", 1.0, G_guess=cold.G * g)
        assert pytest.approx(cold.G, rel=1e-11) == warm.G


def test_low_drive_is_continuous_through_the_floor_and_vanishes() -> None:
    f = cb.FANNO_FLOW_DRIVE_FLOOR
    lo = _flow(PT - f * (1 - 1e-9), 800.0, False, "haaland", 1.0).G
    hi = _flow(PT - f * (1 + 1e-9), 800.0, False, "haaland", 1.0).G
    assert lo == pytest.approx(hi, rel=1e-6)
    assert _flow(PT, 800.0, False, "haaland", 1.0).G == 0.0
    assert 1e-4 * hi > _flow(PT - 1e-6, 800.0, False, "haaland", 1.0).G


def test_a_vacuum_back_pressure_chokes() -> None:
    r = _flow(0.0, 800.0, False, "haaland", 1.0)
    assert r.choked
    Gc = cb.fanno_rayleigh_choked_mass_flux(PT, TT, 800.0, X, L, D, 0.0, "haaland", 1.0)
    assert pytest.approx(Gc, rel=1e-12) == r.G
    assert math.isfinite(r.dG_dTt_exit)
