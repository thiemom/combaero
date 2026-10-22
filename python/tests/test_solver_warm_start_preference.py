"""Elements that converge slowly from the cold start can ask for the warm one
(#706). ``ImpingementCrossflowElement`` does: a chain of compressible momentum
stations converged in 24-145 residual evaluations cold and 16-17 from the
outlet-referenced incompressible proxy. On the stress suite's impingement
family (seeds 0-39, serial) it took the family from 37/40 converged in 167 s
to 40/40 in 33 s, every previously converged root unchanged to 1.6e-7.
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    ImpingementArray,
    ImpingementCrossflowElement,
    NetworkSolver,
    PressureBoundary,
)

D = 0.00254
Y_AIR = cb.mole_to_mass(cb.species.dry_air())


def _pb(name: str, Pt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, 300.0, Y_AIR
    return b


def _array_net() -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    ImpingementArray(
        "ia", n_rows=10, d_jet=D, xn_d=5.0, yn_d=4.0, z_d=2.0, span=0.122, plate_thickness=D
    ).add_to(net, "supply", "exit")
    return net


def _used(solver: NetworkSolver) -> str:
    return solver._diagnostic_data["solver_settings_used"]["init_strategy"]


def test_only_the_declaring_element_asks() -> None:
    assert ImpingementCrossflowElement.prefers_warm_start
    assert not ChannelElement.prefers_warm_start


def test_an_impingement_network_starts_warm() -> None:
    s = NetworkSolver(_array_net())
    res = s.solve()
    assert res["__success__"]
    assert _used(s) == "outletref_warmstart"


def test_a_network_without_one_starts_cold() -> None:
    net = FlowNetwork()
    net.add_node(_pb("a", 1.03e5))
    net.add_node(_pb("b", 1.0e5))
    net.add_element(ChannelElement("c", "a", "b", length=1.0, diameter=0.05))
    s = NetworkSolver(net)
    assert s.solve()["__success__"]
    assert _used(s) == "default"


def test_an_explicit_strategy_is_tried_first(monkeypatch) -> None:
    """The preference reorders only the "default" start. An explicit one runs
    first; its own auto-retry (which on this array is the same warm start)
    is unchanged."""
    s = NetworkSolver(_array_net())
    real, tried = s._solve_impl, []

    def spy(**kw):
        tried.append(kw.get("init_strategy"))
        return real(**kw)

    monkeypatch.setattr(s, "_solve_impl", spy)
    assert s.solve(init_strategy="analytical_pt_prop")["__success__"]
    assert tried[0] == "analytical_pt_prop"


def test_a_failed_warm_start_falls_back_to_the_cold_path(monkeypatch) -> None:
    """The preference only reorders the attempts: if the warm start fails,
    the cold path runs as before and converges the network."""
    s = NetworkSolver(_array_net())
    real = s._solve_impl

    def warm_fails(**kw):
        if kw.get("init_strategy") == "outletref_warmstart" and kw.get("x0") is None:
            return {"__success__": False, "__message__": "forced", "__final_norm__": 1.0}
        return real(**kw)

    monkeypatch.setattr(s, "_solve_impl", warm_fails)
    res = s.solve()
    assert res["__success__"]
    assert _used(s) == "default"


@pytest.mark.parametrize("auto_retry", [False])
def test_no_auto_retry_means_no_reordering(auto_retry: bool) -> None:
    s = NetworkSolver(_array_net())
    assert s.solve(auto_retry=auto_retry)["__success__"]
    assert _used(s) == "default"
