"""A solve that ends on a state diagnostics cannot evaluate reports its outcome (#810).

Found on #802's base: a 6-row impingement array at 0.041 kg/s "converged" (LM,
|F| 7e-4) on a crossflow segment's Fanno continuation with the free exit node
at -113 kPa. The element diagnostics, which raise on a converged state they
cannot evaluate, ran before the second-law check that would have demoted the
root, so ``ChannelElement.diagnostics`` raised out of ``solve()``.

The network here closes its rows exactly at a negative pressure on main: a
constant-fraction loss of 1.5 between a 1 bar supply and a fixed exit flow
puts the exit at -0.5 bar. Its element's diagnostics evaluate the outlet gas
state, as ``ChannelElement``'s do; ``PressureLossElement``'s own guard P <= 0.
"""

from __future__ import annotations

import warnings
from typing import Any

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
    SolveOutcome,
)
from combaero.network.components import PressureLossElement
from combaero.network.pressure_loss import ConstantFractionLoss

X_AIR = list(cb.species.dry_air())
Y_AIR = list(cb.mole_to_mass(X_AIR))


class _GasStateDiagLoss(PressureLossElement):
    """Diagnostics that evaluate the outlet gas state unguarded."""

    def diagnostics(self, state_in: Any, state_out: Any, **_: Any) -> dict[str, float]:
        cs = cb.complete_state(state_out.T, state_out.P, state_out.X)
        return {"rho_out": float(cs.thermo.rho)}


class _BrokenDiagLoss(PressureLossElement):
    """Diagnostics that always fail, as a buggy element's would."""

    def diagnostics(self, state_in: Any, state_out: Any, **_: Any) -> dict[str, float]:
        raise ValueError("diagnostics bug")


def _net(cls: type[PressureLossElement], xi: float) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("supply", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("exit", m_dot=0.01, Tt=300.0, Y=Y_AIR))
    g.add_element(cls("loss", "supply", "exit", correlation=ConstantFractionLoss(xi=xi)))
    return g


def _solve(g: FlowNetwork) -> dict[str, Any]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkSolver(g).solve()


def test_a_root_at_negative_pressure_is_a_failed_solve_not_a_raise() -> None:
    """On main this raised ValueError (pressure must be positive) out of solve()."""
    r = _solve(_net(_GasStateDiagLoss, 1.5))
    assert r["exit.P"] == pytest.approx(-0.5e5, rel=1e-6)  # the rows did close
    assert r["__converged__"] is True
    assert r["__success__"] is False
    assert r["__outcome__"] == SolveOutcome.INCONSISTENT
    assert r["__P_nonpositive__"] == ["exit"]
    assert "at or below zero" in r["__message__"]
    # The element could not be evaluated there: no diagnostics, and the
    # failure is recorded, not dropped.
    assert "rho_out" not in r["__element_diag__"]["loss"]
    recs = [c for c in r["__caught_exceptions__"] if c["phase"] == "diagnostics"]
    assert [c["where"] for c in recs] == ["element 'loss'"]
    assert recs[0]["type"] == "ValueError"
    assert "pressure must be positive" in recs[0]["message"]


def test_a_diagnostics_error_on_a_failed_solve_is_recorded() -> None:
    r = _solve(_net(_BrokenDiagLoss, 1.5))
    assert r["__success__"] is False
    recs = [c for c in r["__caught_exceptions__"] if c["phase"] == "diagnostics"]
    assert recs and recs[0]["message"] == "diagnostics bug"
    assert r["__caught_exception_count__"] >= len(recs)


def test_a_diagnostics_error_on_a_successful_solve_still_raises() -> None:
    """A converged, admissible state the element cannot evaluate is a bug."""
    with pytest.raises(ValueError, match="diagnostics bug"):
        _solve(_net(_BrokenDiagLoss, 0.5))


def test_a_physical_root_is_untouched() -> None:
    r = _solve(_net(_GasStateDiagLoss, 0.5))
    assert r["__success__"] is True
    assert r["__P_nonpositive__"] == []
    assert r["exit.P"] == pytest.approx(0.5e5, rel=1e-6)
    assert r["__element_diag__"]["loss"]["rho_out"] > 0.0
    assert not [c for c in r["__caught_exceptions__"] if c["phase"] == "diagnostics"]


def test_a_channel_pair_past_its_flow_returns_a_failed_result() -> None:
    """Shipped elements only, on main: two 1 m x 10 mm channels from 2 bar
    into a fixed 0.02 kg/s. An auto-retry attempt "converged" with the exit
    below zero pressure and ``ChannelElement.diagnostics`` raised out of
    solve(). Now that attempt is demoted and the solve reports its outcome."""
    from combaero.network import PlenumNode
    from combaero.network.components import ChannelElement

    g = FlowNetwork()
    g.add_node(PressureBoundary("supply", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("mid"))
    g.add_node(MassFlowBoundary("exit", m_dot=0.02, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("c1", "supply", "mid", length=1.0, diameter=0.01))
    g.add_element(ChannelElement("c2", "mid", "exit", length=1.0, diameter=0.01))
    r = _solve(g)
    assert r["__success__"] is False
    assert r["__outcome__"] != SolveOutcome.CONVERGED
    assert r["__message__"]
