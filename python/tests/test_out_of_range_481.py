"""Out-of-range hardening (#481 B6-B8): no throw, no NaN, no jump where a
solver probe can reach -- and exact wherever each closure is defined.

Found by probing every converged solve of the suite through reversal, zero
flow and out-of-range pressures (485 solves): turbulent-only Nusselt
correlations threw in 2300 <= Re < 1e4; the Fanno flow threw at a non-positive
back pressure; the ejector went NaN over most of its low-primary range; the
legacy merging tee went NaN with an idle branch; the junction kernel refused a
state its guards let through; and a node with no inflow jumped to 300 K air
(a momentum chamber to 1 kg/s).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)

AIR = list(cb.species.dry_air())
Y_AIR = list(cb.mole_to_mass(AIR))

# ---------------------------------------------------------------------------
# Turbulent-only Nusselt correlations: the VDI Heat Atlas transition
# ---------------------------------------------------------------------------

T0, P0, D0, L0 = 400.0, 1.5e5, 0.02, 1.0
TURB_ONLY = ["dittus_boelter", "sieder_tate", "petukhov"]


def _velocity(Re: float) -> float:
    rho = cb.density(T0, P0, AIR)
    mu = cb.viscosity(T0, P0, AIR)
    return Re * mu / (rho * D0)


def _smooth(corr: str, Re: float, T: float = T0) -> cb.ChannelResult:
    v = _velocity(Re) * (cb.density(T0, P0, AIR) / cb.density(T, P0, AIR))
    return cb.channel_smooth(T, P0, np.array(AIR), v, D0, L0, 300.0, corr, True, 1.0, 0.0, 1.0, 1.0)


@pytest.mark.parametrize("corr", TURB_ONLY)
def test_turbulent_only_correlations_cover_the_transition(corr: str) -> None:
    """2300 <= Re < 1e4 is the VDI interpolation, linear in Re from laminar Nu
    at 2300 to the correlation at 1e4; continuous at both ends. They threw."""
    nu_lam = 3.66  # NU_LAMINAR_CONST_T: heating
    lo = _smooth(corr, 2300.0 * (1 + 1e-9)).Nu
    hi_in, hi_out = _smooth(corr, 1e4 * (1 - 1e-9)).Nu, _smooth(corr, 1e4 * (1 + 1e-9)).Nu
    assert lo == pytest.approx(nu_lam, rel=1e-6)
    assert hi_in == pytest.approx(hi_out, rel=1e-6)
    mid = _smooth(corr, 6150.0)
    g = (mid.Re - 2300.0) / 7700.0
    assert mid.Nu == pytest.approx((1 - g) * nu_lam + g * hi_out, rel=1e-6)


@pytest.mark.parametrize("corr", TURB_ONLY)
@pytest.mark.parametrize("Re", [5000.0, 3.0e4])
def test_turbulent_only_derivatives_match_differences(corr: str, Re: float) -> None:
    """dh/dmdot and dh/dT, in the band (exact linear interpolation) and above."""
    A = math.pi * D0 * D0 / 4
    r = _smooth(corr, Re)
    rho = cb.density(T0, P0, AIR)
    v = _velocity(Re)
    hv = 1e-6 * v
    hp = cb.channel_smooth(
        T0, P0, np.array(AIR), v + hv, D0, L0, 300.0, corr, True, 1.0, 0.0, 1.0, 1.0
    )
    hm = cb.channel_smooth(
        T0, P0, np.array(AIR), v - hv, D0, L0, 300.0, corr, True, 1.0, 0.0, 1.0, 1.0
    )
    fd_m = (hp.h - hm.h) / (2 * hv * rho * A)
    assert r.dh_dmdot == pytest.approx(fd_m, rel=1e-5)
    hT = 1e-3
    fd_T = (_smooth(corr, Re, T0 + hT).h - _smooth(corr, Re, T0 - hT).h) / (2 * hT)
    assert r.dh_dT == pytest.approx(fd_T, rel=1e-4)


# ---------------------------------------------------------------------------
# Fanno flow against a non-positive back pressure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("P_target", [0.0, -1.0e5])
def test_a_duct_into_zero_or_negative_pressure_is_choked(P_target: float) -> None:
    """It threw (fanno_length_between: need 0 < M1). Continuous with 0+."""
    args = (2e5, 700.0, AIR)
    f = cb.fanno_channel_flow(*args, P_target, True, 1.0, 0.02, 1e-5, "haaland")
    ref = cb.fanno_channel_flow(*args, 1e-6, True, 1.0, 0.02, 1e-5, "haaland")
    assert f.choked and pytest.approx(ref.G, rel=1e-12) == f.G and f.dG_dP_target == 0.0


# ---------------------------------------------------------------------------
# Ejector outside its closures' domains
# ---------------------------------------------------------------------------


def test_the_jet_pump_discharge_is_finite_past_its_domain() -> None:
    """P_py above the secondary supply (no expansion, lambda^2 < 0) gave NaN;
    the C1 floor continues it, monotone, to the secondary stagnation pressure.
    Inside the domain the value is unchanged (pinned by the ejector suite)."""
    Pg, Tg, Pe, Te, om, g = 7.0e5, 440.0, 2.28e4, 280.0, 0.79, 1.401
    p03 = [
        _solver_tools.ejector_jetpump_discharge_and_jacobian(Pg, Tg, Pe, Te, f * Pe, om, g, 1.0).p03
        for f in np.linspace(0.3, 3.0, 55)
    ]
    assert np.all(np.isfinite(p03))
    assert p03[-1] == pytest.approx(Pe, rel=1e-3)


def test_the_ejector_residual_is_finite_down_to_zero_primary_flow() -> None:
    """Row 3 was NaN for every primary flow below 0.9 of the operating one."""
    import test_ejector_element as tej

    s = NetworkSolver(tej._build_network())
    x = np.array(s.solve()["__x_solution__"])
    j = s.unknown_names.index("ch_primary.m_dot")
    for frac in (0.9, 0.5, 0.0, -0.3):
        xp = x.copy()
        xp[j] = frac * x[j]
        f, J = s._residuals_and_jacobian(xp)
        assert np.all(np.isfinite(f)) and np.all(np.isfinite(J.data))


# ---------------------------------------------------------------------------
# Junctions at zero branch flow
# ---------------------------------------------------------------------------


def test_the_merging_tee_with_an_idle_branch_is_finite() -> None:
    """The datum IS the collector there (x = 1, phi = 0): K -> 0 but the closed
    form evaluated 0 * inf."""
    import test_tee_network as ttn

    s = NetworkSolver(ttn._merging_net())
    x = np.array(s.solve()["__x_solution__"])
    xp = x.copy()
    xp[s.unknown_names.index("tee.m_dot_branch")] = 0.0
    f, J = s._residuals_and_jacobian(xp)
    assert np.all(np.isfinite(f)) and np.all(np.isfinite(J.data))


def test_the_merging_tee_with_no_net_supply_is_finite() -> None:
    """Every flow zero (or the suppliers' sum reversed): the pseudodatum is
    stagnant -- dividing by the eps floor gave a negative datum temperature."""
    import test_tee_network as ttn

    s = NetworkSolver(ttn._merging_net())
    x = np.array(s.solve()["__x_solution__"])
    for name in ("tee.m_dot_com", "tee.m_dot_branch"):
        x[s.unknown_names.index(name)] = 0.0
    f, J = s._residuals_and_jacobian(x)
    assert np.all(np.isfinite(f)) and np.all(np.isfinite(J.data))


def test_a_junction_inlet_just_reversed_takes_the_soft_barrier(monkeypatch) -> None:
    """An inlet carrying 1e-9 kg/s outward passed the mass-flow tolerance as
    'not wrong' but was a collector to the velocity masks: no supplier, and
    the kernel refused it with a RuntimeError (Bassett Fig. 7b network)."""
    from validation.junction.models.mpce_network import MPCENetwork

    solvers: list[NetworkSolver] = []
    solve = NetworkSolver.solve

    def spy(self: NetworkSolver, *a, **kw) -> dict:
        solvers.append(self)
        return solve(self, *a, **kw)

    monkeypatch.setattr(NetworkSolver, "solve", spy)
    MPCENetwork(strict=False).evaluate_network(
        "bassett2001", "K6", 0.2, 3.0, math.radians(45.0), topology="imposed_q"
    )
    s = solvers[-1]
    x = np.array(s._last_x if hasattr(s, "_last_x") else s.solve()["__x_solution__"])
    j = s.unknown_names.index("lc_in.m_dot")
    xp = x.copy()
    xp[j] = -1e-8 * x[j]
    f, _ = s._residuals_and_jacobian(xp)
    assert np.all(np.isfinite(f))


# ---------------------------------------------------------------------------
# A node with no inflow
# ---------------------------------------------------------------------------


def _pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _feed_line() -> NetworkSolver:
    """A (600 K) -> momentum chamber -> B (300 K)."""
    g = FlowNetwork()
    g.add_node(_pb("A", 1.05e5, 600.0))
    g.add_node(_pb("B", 1.0e5, 300.0))
    g.add_node(MomentumChamberNode("p"))
    g.add_element(OrificeElement("in", "A", "p", Cd=0.8, diameter=0.01, correlation="fixed"))
    g.add_element(OrificeElement("out", "p", "B", Cd=0.8, diameter=0.01, correlation="fixed"))
    return NetworkSolver(g)


def test_a_chamber_losing_its_last_inflow_carries_no_flow() -> None:
    """With its feed reversed and its outlet still flowing out, a momentum
    chamber has no inflow. It took 1 kg/s there -- a 4.3 kPa jump in its
    dynamic-head row on a Bassett junction port. It now carries none, so its
    row is continuous through the feed's zero."""
    s = _feed_line()
    x = np.array(s.solve()["__x_solution__"])
    j = s.unknown_names.index("in.m_dot")
    node = s.network.nodes["p"]

    def row(m: float) -> float:
        xp = x.copy()
        xp[j] = m
        s._residuals_and_jacobian(xp, compute_jacobian=False)
        return node.residuals(s._get_node_state(node, xp))[0][0]

    r_pos, r_neg = row(1e-12), row(-1e-12)
    assert node._total_m_dot == 0.0
    assert r_neg == pytest.approx(r_pos, abs=1e-6)
