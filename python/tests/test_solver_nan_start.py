"""A solve whose starting residual is NaN must not report convergence.

|F| = NaN fails every "<" in the best-iterate tracking and the "> tol"
success guard alike, so such an attempt came back ``__success__ = True`` at a
point whose residual was NaN and -87 kPa. Found by the stress families
(2026-10-17): a compressible ejector whose cold solve fails is retried from
its incompressible proxy, whose P_jct (-4e4 Pa, a jet-pump quantity) gives
the compressible element a NaN row at the seeded start.
"""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.ejector_element import EjectorElement

Y = list(cb.species.dry_air_mass())
D = 0.05
A = np.pi * D * D / 4


def _ejector_net() -> FlowNetwork:
    """The stress draw (gadgets.ejector, seed 4)."""
    g = FlowNetwork()
    for nid, Pt, Tt in (("bp", 3.0e5, 440.0), ("bs", 10865.4, 280.0), ("bo", 12253.4, 300.0)):
        g.add_node(PressureBoundary(nid, Pt=Pt, Tt=Tt, Y=Y))
    for nid, p in (("mp", 3.0e5), ("ms", 10865.4), ("mo", 12253.4)):
        mc = MomentumChamberNode(nid, area=A)
        mc.initial_guess = {f"{nid}.P": p, f"{nid}.Pt": p}
        g.add_node(mc)
    for eid, a, b, L in (
        ("cp", "bp", "mp", 0.05),
        ("cs", "bs", "ms", 0.05),
        ("co", "mo", "bo", 1.0),
    ):
        g.add_element(ChannelElement(eid, a, b, length=L, diameter=D, regime="incompressible"))
    ej = EjectorElement(
        id="ej",
        primary_node="mp",
        secondary_node="ms",
        outlet_node="mo",
        throat_area=3.14e-5,
        nozzle_exit_area=4.675e-5,
        mixing_area=3.069e-4,
        regime="compressible",
    )
    g.add_element(ej)
    return g


@pytest.fixture
def no_ejector_warm_start(monkeypatch: pytest.MonkeyPatch) -> None:
    """With its own warm start (EjectorElement.seed_warm_start) this cold
    solve converges; the NaN path needs the cold attempt to fail and the
    incompressible proxy to seed the retry, so the seed is switched off."""
    monkeypatch.setattr(EjectorElement, "seed_warm_start", lambda self, graph: None)


@pytest.mark.usefixtures("no_ejector_warm_start")
def test_a_nan_start_is_not_reported_converged() -> None:
    s = NetworkSolver(_ejector_net())
    s.raise_internal_errors = False
    with pytest.warns(RuntimeWarning):
        r = s.solve(timeout=30.0)
    if r["__success__"]:
        res = s._residuals(np.array(r["__x_solution__"]))
        assert np.all(np.isfinite(res)), "converged at a NaN residual"
        assert np.max(np.abs(res)) < 1e-3
    starts = [c for c in r["__caught_exceptions__"] if c["phase"] == "start"]
    assert any(
        "non-finite residual" in c["message"] and "EjectorElement 'ej'" in c["where"]
        for c in starts
    )


@pytest.mark.usefixtures("no_ejector_warm_start")
def test_a_strict_solve_raises_on_a_nan_start() -> None:
    s = NetworkSolver(_ejector_net())
    s.raise_internal_errors = True
    with (
        pytest.raises(FloatingPointError, match="non-finite residual"),
        pytest.warns(RuntimeWarning),
    ):
        s.solve(timeout=30.0)


def test_an_explicit_nan_start_ends_with_outcome_error() -> None:
    s = NetworkSolver(_ejector_net())
    s.raise_internal_errors = False
    x0 = np.array(s._build_x0(), dtype=float)
    x0[s.unknown_names.index("ej.P_jct")] = -4.0e4
    with pytest.warns(RuntimeWarning, match="starting point could not be evaluated"):
        r = s.solve(x0=x0, timeout=10.0)
    assert not r["__success__"]
    assert str(r["__outcome__"]) == str(cb.network.SolveOutcome.ERROR)
