"""#552: the solver does not accept, or crash on, a root whose balance cannot close.

1. A converged root whose network ledger is open in mass or energy is demoted
   to INCONSISTENT (LEDGER_CLOSURE_TOL), as the conftest closure check already
   demands of every test solve.
2. Withdrawing more heat than a stream holds: the mixer used to raise from the
   T(h) inverse (#607 made it refuse a target no positive T reaches), which
   aborted the solve at its starting point. It now continues below the thermo
   floor, so the node lands on the [200, 5000] K clamp and the solve is
   reported INCONSISTENT with the node named.
3. That clamp then exposed a spurious zero-flow root: the isentropic solvers
   floored the static T at 200 K, so a compressible orifice fed at Tt <= 200 K
   returned zero flow under any drive.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    EnergyBoundary,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network import second_law as second_law_mod

Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _cooled_chain(Q: float) -> FlowNetwork:
    g = FlowNetwork()
    for name, Pt in (("a", 1.1e5), ("b", 1.0e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 300.0, Y
        g.add_node(bnd)
    p = PlenumNode("p")
    g.add_node(p)
    g.add_node(PlenumNode("q"))
    for eid, a, b in (("o1", "a", "p"), ("o2", "p", "q"), ("o3", "q", "b")):
        g.add_element(OrificeElement(eid, a, b, Cd=0.8, diameter=0.01, correlation="fixed"))
    p.add_energy_boundary(EnergyBoundary(id="cool", Q=Q))
    return g


@pytest.mark.no_closure_check  # the clamp discards energy: that is the point
def test_cooling_a_stream_below_zero_enthalpy_is_refused_not_raised() -> None:
    """20 kW out of ~6 g/s of 300 K air: no admissible root. main raised
    RuntimeError from calc_T_from_h at the starting point; after that was
    caught, it converged to zero flow under a 10 kPa drive (item 3)."""
    s = NetworkSolver(_cooled_chain(-2e4))
    with pytest.warns(UserWarning, match="clamp"):
        r = s.solve(auto_retry=False)
    assert r["__converged__"] and not r["__success__"]
    assert r["__outcome__"] == "inconsistent"
    assert r["__T_clamped__"] == ["p"]
    assert r["o1.m_dot"] > 1e-3  # flows, not the spurious zero-flow root


def test_admissible_cooling_still_solves() -> None:
    """200 W out of the same stream: 34 K of cooling, a normal root."""
    r = NetworkSolver(_cooled_chain(-200.0)).solve(auto_retry=False)
    assert r["__success__"]
    assert 260.0 < r["p.T"] < 270.0


@pytest.mark.parametrize("Tt", [120.0, 150.0, 199.0, 200.0, 201.0])
def test_the_compressible_orifice_flows_from_a_stagnation_state_below_200_k(Tt: float) -> None:
    """Flow rises as Tt falls (mdot ~ 1/sqrt(Tt)), through 200 K: it was 0 at
    and below 200 K, where the isentropic static T was held at 200 K."""
    area = math.pi * 0.01**2 / 4

    def mdot(T: float) -> float:
        return _solver_tools.orifice_compressible_residuals_and_jacobian(
            0.0, 1.06e5, T, Y, 1.0e5, 0.8, area, 0.0
        ).m_dot_calc

    m = mdot(Tt)
    assert m > 0.0
    # within 1% of the 1/sqrt(Tt) scaling from 300 K (ideal-gas, low Mach)
    assert m == pytest.approx(mdot(300.0) * math.sqrt(300.0 / Tt), rel=1e-2)


def test_the_mixer_continues_below_the_thermo_floor() -> None:
    """A heat withdrawal past what the stream holds returns a temperature
    (linear in Q below the floor, C1 there) instead of raising."""
    st = [cb.MassStream(0.01, 300.0, 1e5, Y)]

    def mix(Q: float) -> _solver_tools.MixerResult:
        return _solver_tools.mixer_from_streams_and_jacobians(st, Q=Q, fraction=0.0)

    Ts = [mix(Q).T_mix for Q in (-1.5e3, -1.7e3, -3e3, -1e4)]
    assert Ts == sorted(Ts, reverse=True)
    assert all(math.isfinite(T) for T in Ts)
    # slope continuous where the continuation starts (between -1.6 and -1.8 kW)
    dQ = 1.0
    slopes = [(mix(Q + dQ).T_mix - mix(Q - dQ).T_mix) / (2 * dQ) for Q in (-1.6e3, -1.8e3)]
    assert slopes[0] == pytest.approx(slopes[1], rel=2e-2)
    assert mix(-1e4).dT_mix_dQ == pytest.approx(slopes[1], rel=1e-6)


@pytest.mark.no_closure_check  # the ledger is forced open: that is the point
def test_a_root_whose_ledger_does_not_close_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    """The guard, by itself: the same converged root with its mass closure
    forced to 1e-4 (the LM roots #552 found were open by 4e-3 to 8e-3)."""
    original = second_law_mod.SecondLawReport.ledger

    def open_ledger(self: second_law_mod.SecondLawReport) -> dict:
        led = original(self)
        led["mass_closure_rel"] = 1e-4
        return led

    monkeypatch.setattr(second_law_mod.SecondLawReport, "ledger", open_ledger)
    with pytest.warns(UserWarning, match="balance does not close"):
        r = NetworkSolver(_cooled_chain(-200.0)).solve(auto_retry=False)
    assert r["__converged__"] and not r["__success__"]
    assert r["__outcome__"] == "inconsistent"


def test_the_ledger_bound_sits_above_the_suites_worst_closure() -> None:
    """Measured over the suite's 25 147 accepted solves: 4.2e-8 (mass) and
    1.1e-7 (energy). The bound keeps a 9x margin over the worse."""
    assert second_law_mod.LEDGER_CLOSURE_TOL >= 9 * 1.1e-7
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        r = NetworkSolver(_cooled_chain(0.0)).solve(auto_retry=False)
    led = r["__ledger__"]
    assert abs(led["mass_closure_rel"]) < 1e-9 and abs(led["energy_closure_rel"]) < 1e-7
