"""What a failed solve hands back.

`NetworkSolver.solve` warns "Returning best iterate" when it fails, and the
returned dictionary is the same shape as a successful one: every unknown by
name, derived node states, per-element diagnostics, and the bookkeeping keys
`__success__`, `__message__`, `__final_norm__`, `__x_solution__`,
`__unknown_names__` and `__convergence_history__`.

It was not doing what it said. Three defects, all in what a FAILED solve
reports rather than in whether it succeeds:

1. `final_x` and `final_norm` are captured immediately after the primary
   root() call. Later phases -- above all the LM fallback -- keep evaluating
   through the same wrapper and improve the tracked best, but only re-point
   `final_x` when they reach the convergence tolerance. An improvement that
   fell short was thrown away. Measured over 38 non-converged junction solves:
   30 returned a state worse than the best they had evaluated, by a median of
   5.8x and up to 2e5x.
2. When the automatic retry ran and also failed, its result was returned
   unconditionally, even when the primary attempt had got closer.
3. `_network_builder` read `__residual_norm__`, a key the solver never sets, so
   every validation record's residual norm was silently infinity.

None of the three changes whether a solve converges. The junction scorecard
(1708 of 2073) and the random-boundary sweep (92.0%, 98.5% inside the
documented Mach range) are identical either side of the fix, which is the point:
this corrects the answer that comes back, not the decision about success.
"""

from __future__ import annotations

import math
import random
import warnings

import numpy as np
import pytest

from combaero.network import NetworkSolver
from validation.junction import random_robustness as rr


def _a_failing_solve():
    """A junction network that does not converge, from the random sweep."""
    rng = random.Random(20260906)
    for _ in range(400):
        case = rr.sample(rng)
        if rr.has_root(case) is not True:
            continue
        solver = NetworkSolver(rr.build(case))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solver.solve(timeout=20.0)
        if not sol["__success__"] and sol["__convergence_history__"]:
            return solver, sol
    pytest.skip("no failing case found in the sampled draws")


@pytest.fixture(scope="module")
def failed():
    return _a_failing_solve()


# ---------------------------------------------------------------------------
# What comes back at all
# ---------------------------------------------------------------------------


def test_a_failed_solve_returns_a_full_solution(failed):
    _, sol = failed

    assert sol["__success__"] is False
    assert sol["__message__"]
    physical = [k for k in sol if not k.startswith("__")]
    assert len(physical) > 10, "a failed solve must still describe the state it reached"
    assert "lc_com.m_dot" in sol


def test_the_bookkeeping_keys_are_the_ones_callers_read(failed):
    """Pinned because a consumer read `__residual_norm__` for months and got
    infinity every time: the solver has never set that name."""
    _, sol = failed

    for key in (
        "__success__",
        "__message__",
        "__final_norm__",
        "__x_solution__",
        "__unknown_names__",
        "__convergence_history__",
    ):
        assert key in sol, f"{key} missing from the returned solution"
    assert "__residual_norm__" not in sol


def test_a_failed_solve_warns():
    """A caller who ignores `__success__` must at least see a warning, since
    the state handed back is the best iterate and not a solution."""
    rng = random.Random(20260906)
    for _ in range(400):
        case = rr.sample(rng)
        if rr.has_root(case) is not True:
            continue
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            sol = NetworkSolver(rr.build(case)).solve(timeout=20.0)
        if sol["__success__"]:
            continue
        assert any("did not converge" in str(w.message) for w in caught), (
            f"failed silently: {sol['__message__'][:80]}"
        )
        return
    pytest.skip("no failing case found in the sampled draws")


# ---------------------------------------------------------------------------
# It is the BEST iterate
# ---------------------------------------------------------------------------


def test_the_returned_norm_is_the_best_one_evaluated(failed):
    _, sol = failed
    best = min(h["norm"] for h in sol["__convergence_history__"])

    assert sol["__final_norm__"] == pytest.approx(best, rel=1e-9), (
        "the solve evaluated a better point than the one it returned"
    )


def test_the_returned_state_really_has_that_residual(failed):
    """The norm and the state must describe the same point, or a caller
    warm-starting from the result gets something other than it was told."""
    solver, sol = failed
    x = np.array(sol["__x_solution__"])

    assert float(np.linalg.norm(solver._residuals(x))) == pytest.approx(
        sol["__final_norm__"], rel=1e-6
    )


def test_the_warm_start_point_matches_what_was_returned(failed):
    solver, sol = failed

    assert np.allclose(solver._last_x, np.array(sol["__x_solution__"]))
    assert np.allclose(solver._diagnostic_data["solution"], np.array(sol["__x_solution__"]))


@pytest.mark.parametrize("seed", [20260906, 4242])
def test_no_solve_returns_worse_than_it_reached(seed):
    """The whole class, swept. Before the fix this failed on 30 of 38."""
    rng = random.Random(seed)
    checked = 0
    for _ in range(250):
        case = rr.sample(rng)
        if rr.has_root(case) is not True:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = NetworkSolver(rr.build(case)).solve(timeout=20.0)
        history = sol["__convergence_history__"]
        if sol["__success__"] or not history:
            continue
        checked += 1
        best = min(h["norm"] for h in history)
        assert sol["__final_norm__"] <= best * 1.01 + 1e-12, (
            f"returned |F|={sol['__final_norm__']:.4e} against a best of {best:.4e}"
        )
        if checked >= 12:
            break
    assert checked > 0, "no failing solves sampled"


# ---------------------------------------------------------------------------
# The validation harness records a real number again
# ---------------------------------------------------------------------------


def test_the_junction_harness_records_a_finite_residual_norm():
    from validation.junction.models.mpce_v2_network import MPCEv2Network

    model = MPCEv2Network(strict=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        converged = model.evaluate_network("bassett2001", "K6", 0.6, 3.0, math.radians(45.0))

    assert converged.converged
    assert math.isfinite(converged.residual_norm), (
        "the harness is reading a key the solver does not set"
    )
    assert converged.residual_norm < 1e-3


# ---------------------------------------------------------------------------
# The automatic retry must not replace a closer result with a worse one
# ---------------------------------------------------------------------------


def _retryable_network():
    """A junction network with a compressible element, so the auto-retry path
    is applicable at all (it needs an MPCE element AND something compressible)."""
    import combaero as cb
    from combaero.network import (
        ChannelElement,
        FlowNetwork,
        LosslessConnectionElement,
        MomentumChamberNode,
        PressureBoundary,
    )
    from combaero.network.mpce_v2_element import MPCEv2Element

    Y = list(cb.mole_to_mass(cb.species.dry_air()))
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_in", Pt=2.1e5, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=2.05e5, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_bra", Pt=2.0e5, Tt=300.0, Y=Y))
    for pid in ("port_com", "port_str", "port_bra"):
        net.add_node(MomentumChamberNode(pid, area=0.01))
    net.add_element(
        ChannelElement(
            "ch_in", "pb_in", "port_com", length=0.3, diameter=0.05, regime="compressible"
        )
    )
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    net.add_element(LosslessConnectionElement("lc_bra", "port_bra", "pb_bra"))
    net.add_element(
        MPCEv2Element(
            id="jct",
            inlet_nodes=["port_com"],
            outlet_nodes=["port_str", "port_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[0.01, 0.01, 0.01],
            flow_direction="branch",
            strict=False,
        )
    )
    return net


def test_the_retry_is_only_returned_when_it_got_closer(monkeypatch):
    """Both attempts fail and the retry lands FURTHER away.

    The retry used to be returned unconditionally, so a near-miss could be
    replaced by something far worse while the warning still promised the best
    iterate. Driven through `_solve_impl` rather than a real stiff network, so
    the comparison is tested rather than a particular network's luck.
    """
    solver = NetworkSolver(_retryable_network())
    calls = {"n": 0}

    def fake_impl(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            return {
                "__success__": False,
                "__message__": "primary",
                "__final_norm__": 1.0,
                "marker": "primary",
            }
        return {
            "__success__": False,
            "__message__": "retry",
            "__final_norm__": 100.0,
            "marker": "retry",
        }

    monkeypatch.setattr(solver, "_solve_impl", fake_impl)
    monkeypatch.setattr(solver, "_outlet_ref_incompressible_seed", lambda **kw: np.zeros(3))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = solver.solve(timeout=20.0)

    assert calls["n"] == 2, "the auto-retry did not run, so nothing was compared"
    assert sol["marker"] == "primary", "a worse retry replaced the closer primary"
    assert "auto-retry also" in sol["__message__"]


def test_a_closer_retry_is_returned(monkeypatch):
    """The other side of the comparison, so the guard cannot just always
    return the primary."""
    solver = NetworkSolver(_retryable_network())
    calls = {"n": 0}

    def fake_impl(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            return {
                "__success__": False,
                "__message__": "primary",
                "__final_norm__": 100.0,
                "marker": "primary",
            }
        return {
            "__success__": False,
            "__message__": "retry",
            "__final_norm__": 1.0,
            "marker": "retry",
        }

    monkeypatch.setattr(solver, "_solve_impl", fake_impl)
    monkeypatch.setattr(solver, "_outlet_ref_incompressible_seed", lambda **kw: np.zeros(3))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = solver.solve(timeout=20.0)

    assert calls["n"] == 2
    assert sol["marker"] == "retry"
