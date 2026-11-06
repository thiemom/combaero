"""A fraction loss without m_ref on its zero-flow blend is flagged (#1004).

``ConstantFractionLoss`` / ``LinearThetaFractionLoss`` without ``m_ref`` hold
dPt = xi Pt at every flow except zero. Where the network forces a smaller drop
across one, its only root is the blend of its forward and reverse forms at
|m_dot| < FLOW_REVERSAL_M_SCALE, and the solve converged there reporting
success with nothing to say the loss law is not met.

``fraction_unanchored`` flags it (1) when its flow component carries flow;
``__fraction_unanchored__`` lists it. Reported only: the solve is not
demoted. A whole component at ~zero flow (zero drive, or the loss in series
with nothing else to carry the flow) is not flagged, nor is a loss with
m_ref or one that flows. On a failed solve the losses are named in the
message.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EffusionPlateElement,
    FlowNetwork,
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


def test_a_blocked_loss_in_series_is_a_whole_component_at_zero_flow() -> None:
    """Nothing else in its component flows: the rule's other side, by
    decision (a component at ~zero flow is not flagged)."""
    r = _solve(_series(BLOCKED, ConstantFractionLoss(xi=0.3)))
    assert r["__success__"], r["__message__"]
    assert abs(r["loss.m_dot"]) < FLOW_REVERSAL_M_SCALE
    assert r["loss.fraction_unanchored"] == 0.0
    assert r["__fraction_unanchored__"] == []


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
