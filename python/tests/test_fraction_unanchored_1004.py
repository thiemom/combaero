"""A fraction loss without m_ref on its zero-flow blend is flagged (#1004).

``ConstantFractionLoss`` / ``LinearThetaFractionLoss`` without ``m_ref`` hold
dPt = xi Pt at every flow except zero. Where the network forces a smaller drop
across one, its only root is the blend of its forward and reverse forms at
|m_dot| < FLOW_REVERSAL_M_SCALE, and the solve converged there reporting
success with nothing to say the loss law is not met.

``fraction_unanchored`` flags it (1) when its flow component is DRIVEN,
judged from the inputs alone (boundary Pt spread, a fed flow, a fuel
controller's target, a compressor or turbine), and it is not a dead end by
topology; ``__fraction_unanchored__`` lists it. Reported only: the solve is
not demoted. Zero drive, a dead end, a loss with m_ref or one that flows is
not flagged. On a failed solve the losses are named in the message. The
check runs once, on the outermost solve's result, after every retry has
restored the boundaries.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    CombustorNode,
    CompressorElement,
    EffusionPlateElement,
    FlowNetwork,
    FuelBoundary,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
)
from combaero.network.components import FLOW_REVERSAL_M_SCALE, NetworkMixtureState
from combaero.network.pressure_loss import (
    ConstantFractionLoss,
    ConstantHeadLoss,
    LinearThetaFractionLoss,
)

Y = list(cb.species.dry_air_mass())
CH4 = list(cb.species.to_mass(cb.species.pure_species("CH4")))
P_LO = 1.0e5


def _head(eid: str, a: str, b: str) -> PressureLossElement:
    return PressureLossElement(
        eid, a, b, correlation=ConstantHeadLoss(zeta=2.0, area=1e-3), area=1e-3
    )


def _parallel(P_a: float, correlation: object) -> FlowNetwork:
    """a -> main (head loss) -> b, and beside it a -> loss -> p -> part -> b."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=P_a, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("b", Pt=P_LO, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("p"))
    g.add_element(_head("main", "a", "b"))
    g.add_element(PressureLossElement("loss", "a", "p", correlation=correlation, area=1e-3))
    g.add_element(_head("part", "p", "b"))
    return g


def _series(P_a: float, correlation: object) -> FlowNetwork:
    """a -> loss -> p -> part -> b."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=P_a, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("b", Pt=P_LO, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("p"))
    g.add_element(PressureLossElement("loss", "a", "p", correlation=correlation, area=1e-3))
    g.add_element(_head("part", "p", "b"))
    return g


def _solve(g: FlowNetwork, **kw: object) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkSolver(g).solve(**kw)


#: 10 % drive across a 30 % loss: the branch can only sit on its blend.
BLOCKED = 1.1e5


@pytest.mark.parametrize(
    "correlation",
    [ConstantFractionLoss(xi=0.3), LinearThetaFractionLoss(k=0.5, xi0=0.3)],
    ids=["constant", "linear_theta"],
)
def test_a_blocked_fraction_loss_in_a_flowing_network_is_flagged(correlation: object) -> None:
    r = _solve(_parallel(BLOCKED, correlation))
    assert r["__success__"], r["__message__"]
    assert str(r["__outcome__"]) == "converged"  # reported, never demoted
    assert abs(r["loss.m_dot"]) < FLOW_REVERSAL_M_SCALE
    assert r["main.m_dot"] > 0.05
    assert r["loss.fraction_unanchored"] == 1.0
    assert r["__element_diag__"]["loss"]["fraction_unanchored"] == 1.0
    assert r["__fraction_unanchored__"] == ["loss"]
    # The head losses beside it are anchored by construction.
    assert r["main.fraction_unanchored"] == 0.0
    assert r["part.fraction_unanchored"] == 0.0


def test_the_same_loss_with_m_ref_is_silent() -> None:
    r = _solve(_parallel(BLOCKED, ConstantFractionLoss(xi=0.3, m_ref=0.05)))
    assert r["__success__"], r["__message__"]
    assert r["loss.m_dot"] > FLOW_REVERSAL_M_SCALE  # it flows, its loss faded
    assert r["loss.fraction_unanchored"] == 0.0
    assert r["__fraction_unanchored__"] == []


def test_a_flowing_fraction_loss_is_silent() -> None:
    r = _solve(_parallel(1.5e5, ConstantFractionLoss(xi=0.3)))
    assert r["__success__"], r["__message__"]
    assert r["loss.m_dot"] > 0.05
    assert r["loss.fraction_unanchored"] == 0.0
    assert r["__fraction_unanchored__"] == []


@pytest.mark.parametrize("build", [_parallel, _series], ids=["parallel", "series"])
def test_zero_drive_is_an_answer_and_silent(build: object) -> None:
    r = _solve(build(P_LO, ConstantFractionLoss(xi=0.3)))
    assert r["__success__"], r["__message__"]
    assert r["loss.m_dot"] == pytest.approx(0.0, abs=1e-9)
    assert r["loss.fraction_unanchored"] == 0.0
    assert r["__fraction_unanchored__"] == []


def test_a_blocked_loss_in_series_is_flagged() -> None:
    """The whole component sits on the blend, but it is driven (10 % across
    a 30 % loss): the loss law is not met."""
    r = _solve(_series(BLOCKED, ConstantFractionLoss(xi=0.3)))
    assert r["__success__"], r["__message__"]
    assert abs(r["loss.m_dot"]) < FLOW_REVERSAL_M_SCALE
    assert r["loss.fraction_unanchored"] == 1.0
    assert r["__fraction_unanchored__"] == ["loss"]


def test_zero_fed_flow_is_silent() -> None:
    """A mass-flow boundary feeding 0 kg/s into a loss to an exit: no drive."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("a", m_dot=0.0, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("b", Pt=P_LO, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("p"))
    g.add_element(
        PressureLossElement("loss", "a", "p", correlation=ConstantFractionLoss(xi=0.3), area=1e-3)
    )
    g.add_element(_head("part", "p", "b"))
    r = _solve(g)
    assert r["__success__"], r["__message__"]
    assert r["loss.fraction_unanchored"] == 0.0
    assert r["__fraction_unanchored__"] == []


def test_a_dead_end_behind_a_fraction_loss_is_silent() -> None:
    """a -> main -> b flows; a -> loss -> d -> stub -> e, a mass-flow
    boundary at 0 kg/s (a closed end): zero flow through the loss by mass
    balance, not a blocked loss."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=BLOCKED, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("b", Pt=P_LO, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("d"))
    g.add_node(MassFlowBoundary("e", m_dot=0.0, Tt=300.0, Y=Y))
    g.add_element(_head("main", "a", "b"))
    g.add_element(
        PressureLossElement("loss", "a", "d", correlation=ConstantFractionLoss(xi=0.3), area=1e-3)
    )
    g.add_element(_head("stub", "d", "e"))
    s = NetworkSolver(g)
    assert s._dead_end("loss") and s._dead_end("stub") and not s._dead_end("main")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    assert abs(r["loss.m_dot"]) < FLOW_REVERSAL_M_SCALE
    assert r["loss.fraction_unanchored"] == 0.0
    assert r["__fraction_unanchored__"] == []


def test_a_closed_pocket_is_a_dead_end_and_a_fed_one_is_not() -> None:
    """Topology only: a pocket of plenums joined to the rest by the loss
    alone is a dead end; a nonzero feed into it opens it."""

    def pocket(m_feed: float | None) -> FlowNetwork:
        g = FlowNetwork()
        g.add_node(PressureBoundary("a", Pt=BLOCKED, Tt=300.0, Y=Y))
        g.add_node(PressureBoundary("b", Pt=P_LO, Tt=300.0, Y=Y))
        for n in ("d", "e"):
            g.add_node(PlenumNode(n))
        g.add_element(_head("main", "a", "b"))
        g.add_element(
            PressureLossElement(
                "loss", "a", "d", correlation=ConstantFractionLoss(xi=0.3), area=1e-3
            )
        )
        g.add_element(_head("de1", "d", "e"))
        g.add_element(_head("de2", "e", "d"))
        if m_feed is not None:
            g.add_node(MassFlowBoundary("m", m_dot=m_feed, Tt=300.0, Y=Y))
            g.add_element(_head("feed", "m", "e"))
        return g

    assert NetworkSolver(pocket(None))._dead_end("loss")
    assert NetworkSolver(pocket(0.0))._dead_end("loss")
    assert not NetworkSolver(pocket(0.01))._dead_end("loss")


def test_a_compressor_between_equal_boundaries_drives_a_blocked_loss() -> None:
    """a -> comp (PR 1.05) -> p; p -> head -> b flows; p -> loss (30 %) -> q
    -> head -> b is blocked. The boundaries are equal: only (d) drives."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=P_LO, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("b", Pt=P_LO, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("p"))
    g.add_node(PlenumNode("q"))
    g.add_element(CompressorElement("comp", "a", "p", 1.05, eta_is=0.85))
    g.add_element(_head("main", "p", "b"))
    g.add_element(
        PressureLossElement("loss", "p", "q", correlation=ConstantFractionLoss(xi=0.3), area=1e-3)
    )
    g.add_element(_head("part", "q", "b"))
    s = NetworkSolver(g)
    assert "CompressorElement 'comp'" in s._component_drive(list(g.nodes), list(g.elements))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    assert r["main.m_dot"] > 0.01
    assert abs(r["loss.m_dot"]) < FLOW_REVERSAL_M_SCALE
    assert r["__fraction_unanchored__"] == ["loss"]


def test_the_drive_is_read_from_the_inputs() -> None:
    """(a) the Pt tolerance, (b) a fed flow, (c) a fuel controller's target,
    with no solve: the fuel controller's flow is an unknown, so its target is
    what is read."""

    def drive(g: FlowNetwork) -> str:
        return NetworkSolver(g)._component_drive(list(g.nodes), list(g.elements))

    for dP, driven in ((1e-12 * P_LO, False), (1e-6 * P_LO, True)):
        g = _series(P_LO + dP, ConstantFractionLoss(xi=0.3))
        assert bool(drive(g)) is driven, dP
    g = _series(P_LO, ConstantFractionLoss(xi=0.3))
    g.add_node(MassFlowBoundary("m", m_dot=0.01, Tt=300.0, Y=Y))
    g.add_element(_head("feed", "m", "p"))
    assert "MassFlowBoundary 'm'" in drive(g)

    def fuel_only(target: float) -> FlowNetwork:
        n = FlowNetwork()
        n.add_node(PressureBoundary("air", Pt=P_LO, Tt=400.0, Y=Y))
        n.add_node(PressureBoundary("out", Pt=P_LO, Tt=400.0, Y=Y))
        n.add_node(FuelBoundary("fuel", combustor="burner", Tt_target=target, Y=CH4))
        n.add_node(CombustorNode("burner"))
        n.add_element(_head("inlet", "air", "burner"))
        n.add_element(LosslessConnectionElement("inj", "fuel", "burner"))
        n.add_element(_head("exit", "burner", "out"))
        return n

    assert "FuelBoundary 'fuel' targets" in drive(fuel_only(1500.0))
    assert drive(fuel_only(350.0)) == ""


def test_the_check_runs_once_on_the_restored_boundaries(monkeypatch: pytest.MonkeyPatch) -> None:
    """The boundary ramp (#642) scales the boundaries for its sub-solves and
    restores them; the check runs once, after it, at the user's Pt."""
    g = _series(BLOCKED, ConstantFractionLoss(xi=0.3))
    s = NetworkSolver(g)
    s.bc_ramp_retry = True
    s.homotopy_retry = False
    seen_pt: list[float] = []
    marks: list[float] = []
    real_retries = NetworkSolver._solve_flow_retries
    real_mark = NetworkSolver._mark_fraction_unanchored
    calls = {"n": 0}

    def retries(self: NetworkSolver, **kw: object) -> dict:
        calls["n"] += 1
        seen_pt.append(float(self.network.nodes["a"].Pt))
        if calls["n"] == 1:
            return {"__success__": False, "__converged__": False, "__message__": "forced"}
        return real_retries(self, **kw)

    def mark(self: NetworkSolver, sol: dict) -> dict:
        marks.append(float(self.network.nodes["a"].Pt))
        return real_mark(self, sol)

    monkeypatch.setattr(NetworkSolver, "_solve_flow_retries", retries)
    monkeypatch.setattr(NetworkSolver, "_mark_fraction_unanchored", mark)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=30.0)
    assert any(pt != BLOCKED for pt in seen_pt), "the ramp never scaled the boundaries"
    assert marks == [BLOCKED]
    assert r["__success__"], r["__message__"]
    assert r["__fraction_unanchored__"] == ["loss"]


def test_the_key_is_nan_until_the_solver_sets_it() -> None:
    """The element alone cannot tell (it needs its component's flow)."""
    g = _parallel(BLOCKED, ConstantFractionLoss(xi=0.3))
    el = g.elements["loss"]
    assert el.unanchored_fraction
    assert not g.elements["main"].unanchored_fraction
    assert not PressureLossElement(
        "z", "a", "p", correlation=ConstantFractionLoss(xi=0.0)
    ).unanchored_fraction
    assert not PressureLossElement(
        "m", "a", "p", correlation=ConstantFractionLoss(xi=0.3, m_ref=0.1)
    ).unanchored_fraction
    st = NetworkMixtureState(P=1.1e5, Pt=1.1e5, Tt=300.0, m_dot=0.0, Y=Y)
    d = el.diagnostics(st, st)
    assert math.isnan(d["fraction_unanchored"])


def test_a_failed_solve_names_the_blocked_loss() -> None:
    """Part C: the solve is refused elsewhere (an effusion plate with no
    cooling path, #998, in a second component); the blocked loss is named in
    the message, and still listed."""
    g = _parallel(BLOCKED, ConstantFractionLoss(xi=0.3))
    for n in ("gas_in", "exit", "cool"):
        g.add_node(PressureBoundary(n, Pt=P_LO, Tt=600.0, Y=Y))
    g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
    for eid, a, b in (("duct", "gas_in", "ch"), ("tail", "ch", "exit")):
        g.add_element(ChannelElement(eid, a, b, length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(
        EffusionPlateElement(
            "eff",
            "cool",
            "ch",
            hole_diameter=0.8e-3,
            wall_thickness=2e-3,
            pitch=6e-3,
            panel_length=0.1,
            panel_width=0.1,
            angle_deg=30.0,
        )
    )
    r = _solve(g)
    assert not r["__success__"]
    assert "no cooling path" in r["__message__"]
    assert r["__fraction_unanchored__"] == ["loss"]
    assert "PressureLossElement 'loss': a fraction loss without m_ref" in r["__message__"]


def test_a_converged_solve_does_not_name_it_in_the_message() -> None:
    r = _solve(_parallel(BLOCKED, ConstantFractionLoss(xi=0.3)))
    assert r["__fraction_unanchored__"] == ["loss"]
    assert "without m_ref" not in r["__message__"]
