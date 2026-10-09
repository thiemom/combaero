"""Lazy Jacobian assembly (#489).

The root finder's residual callback evaluates F only, skipping the state
sensitivity relay (dT/dx, dY/dx); the Jacobian is assembled only when MINPACK
asks for it -- at the start and on restarts for hybr, 18% of evaluations on a
36-case mixing study. When it asks at the point F was just evaluated, the
deferred assembly reuses that evaluation's element derivatives and rebuilds
the relay at the same x, so the answer must be the eager Jacobian (to the
round-off of the elements' warm-started inner solves).

Measured against eager assembly on 292 networks: identical solutions
(bitwise), identical evaluation counts; impingement arrays 2.1x faster,
random junctions 24% faster.
"""

from __future__ import annotations

import random

import numpy as np
import pytest
import test_impingement_array as tia
import test_solver_robustness_481 as tsr

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from validation.junction import random_robustness as rr


def _pb(name: str, Pt: float, Tt: float, Y: list[float]) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y
    return b


def _mixing() -> FlowNetwork:
    """Air and CO2 mix in a plenum feeding a compressible duct."""
    air = list(cb.mole_to_mass(cb.species.dry_air()))
    X = [0.0] * cb.num_species()
    X[cb.species_index_from_name("CO2")] = 1.0
    g = FlowNetwork()
    g.add_node(_pb("air", 2.0e5, 300.0, air))
    g.add_node(_pb("co2", 2.05e5, 500.0, list(cb.mole_to_mass(X))))
    g.add_node(_pb("out", 1.0e5, 300.0, air))
    g.add_node(PlenumNode("p"))
    g.add_element(OrificeElement("oa", "air", "p", Cd=0.8, diameter=0.015, correlation="fixed"))
    g.add_element(OrificeElement("oc", "co2", "p", Cd=0.8, diameter=0.012, correlation="fixed"))
    g.add_element(
        ChannelElement(
            "c", "p", "out", length=1.0, diameter=0.02, roughness=1e-5, regime="compressible"
        )
    )
    return g


def _reversed() -> FlowNetwork:
    g = FlowNetwork()
    air = list(cb.mole_to_mass(cb.species.dry_air()))
    g.add_node(_pb("A", 1.02e5, 400.0, air))
    g.add_node(_pb("B", 1.03e5, 600.0, air))
    g.add_node(PlenumNode("p"))
    g.add_element(OrificeElement("e1", "A", "p", Cd=0.8, diameter=0.01, correlation="fixed"))
    g.add_element(OrificeElement("e2", "p", "B", Cd=0.8, diameter=0.01, correlation="fixed"))
    return g


def _junction(seed: int) -> FlowNetwork:
    return rr.build(rr.sample(random.Random(seed)))


NETWORKS = {
    "mixing": _mixing,
    "walls": tsr._series_walls,
    "reversed": _reversed,
    "impingement": lambda: tia._net(tia._array(5), hot=True),
    "junction_a": lambda: _junction(20260906),
    "junction_b": lambda: _junction(7),
}


def _points(s: NetworkSolver) -> list[np.ndarray]:
    """The solution and two perturbed points around it."""
    x = np.array(s.solve()["__x_solution__"])
    rng = np.random.default_rng(3)
    return [x] + [x * (1.0 + 0.02 * rng.standard_normal(x.size)) for _ in range(2)]


@pytest.mark.parametrize("name", sorted(NETWORKS))
def test_the_deferred_jacobian_is_the_eager_one(name: str) -> None:
    s = NetworkSolver(NETWORKS[name]())
    for x in _points(s):
        res_eager, J_eager = s._residuals_and_jacobian(x)
        res_lazy, none = s._residuals_and_jacobian(x, compute_jacobian=False)
        assert none is None
        # Two evaluations at one x agree to the elements' inner solves
        # (warm-started roots, the wall k(T) fixed point), not bitwise.
        np.testing.assert_allclose(res_lazy, res_eager, rtol=1e-9, atol=1e-9)
        J_lazy, J_ref = s._jacobian_at_last(x).toarray(), J_eager.toarray()
        np.testing.assert_allclose(J_lazy, J_ref, rtol=1e-9, atol=1e-12 * np.abs(J_ref).max())


def test_a_deferred_jacobian_is_only_given_for_its_own_point() -> None:
    s = NetworkSolver(_mixing())
    x = np.array(s.solve()["__x_solution__"])
    s._residuals_and_jacobian(x, compute_jacobian=False)
    assert s._jacobian_at_last(x * 1.001) is None
    s._residuals_and_jacobian(x)  # an eager evaluation clears it
    assert s._jacobian_at_last(x) is None


def test_a_solve_assembles_far_fewer_jacobians_than_it_evaluates() -> None:
    """hybr wants J at the start and on restarts, not every evaluation."""
    s = NetworkSolver(tia._net(tia._array(8), hot=True))
    calls = {"F": 0, "relay": 0}
    evaluate, propagate = s._residuals_and_jacobian, s._propagate_states

    def counting_evaluate(x: np.ndarray, **kw: bool) -> tuple:
        calls["F"] += 1
        return evaluate(x, **kw)

    def counting_propagate(x: np.ndarray, with_relay: bool = True) -> dict:
        calls["relay"] += with_relay
        return propagate(x, with_relay)

    s._residuals_and_jacobian = counting_evaluate
    s._propagate_states = counting_propagate
    r = s.solve(auto_retry=False)
    assert r["__success__"]
    assert calls["relay"] < calls["F"] / 3
