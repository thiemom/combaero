"""The retry chain's budget split when the last retries are on (#1018).

At the stress harness's 20 s budget the Newton homotopy took its 5 s before
the boundary-condition ramp, which is what solves random 81 and 119, and
rescued neither; its own stress rescues (random 60, 91, 136) take under 2 s.
The cap is 3 s, not 2: each step's timeout is what the cap leaves, and a
solve runs its LM fallback only with more than 2 s left (random 91's last
homotopy step needs it; at 2 s the draw was lost).
So with the ramp to follow it gets at most NEWTON_HOMOTOPY_BUDGET_BEFORE_RAMP.
A network that prefers the warm start (#706) no longer solves from the same
outlet-referenced seed twice.
"""

from __future__ import annotations

from typing import Any

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    ImpingementArray,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network import solver as solver_mod

AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _net() -> FlowNetwork:
    n = FlowNetwork()
    n.add_node(PressureBoundary("a", Pt=1.2e5, Tt=300.0, Y=AIR))
    n.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0, Y=AIR))
    n.add_element(ChannelElement("c", "a", "b", length=0.1, diameter=0.02))
    return n


def _stages(
    s: NetworkSolver, monkeypatch: pytest.MonkeyPatch, *, ramp: bool, cap: float | None
) -> list[Any]:
    """Every stage of ``s`` fails; returns the stages it ran, in order, with
    the homotopy's budget."""
    calls: list[Any] = []

    def chain(**kw: Any) -> dict[str, Any]:
        calls.append("chain")
        return {"__success__": False, "__converged__": False, "__final_norm__": 1.0}

    def homotopy(failed: dict[str, Any], *a: Any) -> dict[str, Any]:
        calls.append(("homotopy", s.homotopy_budget))
        return failed

    def bc_ramp(failed: dict[str, Any], *a: Any) -> dict[str, Any]:
        calls.append("ramp")
        return failed

    monkeypatch.setattr(s, "_solve_flow_retries", chain)
    monkeypatch.setattr(s, "_newton_homotopy_retry", homotopy)
    monkeypatch.setattr(s, "_bc_ramp_retry", bc_ramp)
    s.homotopy_retry = True
    s.homotopy_budget = cap
    s.bc_ramp_retry = ramp
    assert not s.solve(timeout=20.0)["__success__"]
    assert s.homotopy_budget == cap  # restored
    return calls


def test_before_the_ramp_the_homotopy_gets_at_most_three_seconds(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Before #1018 it ran with its full 5 s."""
    calls = _stages(
        NetworkSolver(_net()), monkeypatch, ramp=True, cap=solver_mod.NEWTON_HOMOTOPY_BUDGET
    )
    assert solver_mod.NEWTON_HOMOTOPY_BUDGET_BEFORE_RAMP == 3.0
    assert calls == ["chain", ("homotopy", 3.0), "ramp"]


def test_a_smaller_cap_is_kept(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _stages(NetworkSolver(_net()), monkeypatch, ramp=True, cap=0.5)
    assert calls == ["chain", ("homotopy", 0.5), "ramp"]


def test_without_the_ramp_the_homotopy_keeps_its_budget(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = _stages(
        NetworkSolver(_net()), monkeypatch, ramp=False, cap=solver_mod.NEWTON_HOMOTOPY_BUDGET
    )
    assert calls == ["chain", ("homotopy", solver_mod.NEWTON_HOMOTOPY_BUDGET)]


def test_uncapped_it_still_takes_half_of_the_time_left(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The GUI's homotopy_budget = None: half of what the chain left."""
    calls = _stages(NetworkSolver(_net()), monkeypatch, ramp=True, cap=None)
    assert calls[0] == "chain" and calls[2] == "ramp"
    assert calls[1][0] == "homotopy" and 9.0 < calls[1][1] <= 10.0


def test_the_warm_start_seed_is_not_solved_from_twice(monkeypatch: pytest.MonkeyPatch) -> None:
    """An impingement array prefers the outlet-referenced warm start (#706).
    When that solve fails and the cold one too, the auto-retry used to build
    the same outlet-referenced seed and solve from it again (the CR-3630
    Table 4.1 case: |F| 0.017108 twice, 7 s each). It now goes on to the
    inlet-referenced seed."""
    net = FlowNetwork()
    for name, Pt in (("supply", 1.03e5), ("exit", 1.0e5)):
        net.add_node(PressureBoundary(name, Pt=Pt, Tt=300.0, Y=AIR))
    ImpingementArray(
        "ia",
        n_rows=4,
        d_jet=0.00254,
        xn_d=5.0,
        yn_d=4.0,
        z_d=2.0,
        span=0.122,
        plate_thickness=0.00254,
    ).add_to(net, "supply", "exit")
    s = NetworkSolver(net)
    real_impl, real_seed = s._solve_impl, s._outlet_ref_incompressible_seed
    seeds: list[str] = []

    def impl(**kw: Any) -> dict[str, Any]:
        r = real_impl(**kw)
        if kw.get("x0") is None and not getattr(s, "_in_warmstart_proxy", False):
            # The warm-started and the cold compressible solves fail.
            return {**r, "__success__": False, "__converged__": False, "__final_norm__": 1.0}
        return r

    def seed(**kw: Any) -> Any:
        seeds.append(kw.get("p_ref", "outlet"))
        return real_seed(**kw)

    monkeypatch.setattr(s, "_solve_impl", impl)
    monkeypatch.setattr(s, "_outlet_ref_incompressible_seed", seed)
    res = s.solve(timeout=30.0)
    assert seeds == ["outlet", "inlet"]
    assert res["__success__"]
    assert s._diagnostic_data["solver_settings_used"]["init_strategy"] == "inletref_warmstart"
