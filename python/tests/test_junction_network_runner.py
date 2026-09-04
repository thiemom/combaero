"""Smoke tests for the network-mode validation runner."""

from __future__ import annotations

import math

import pytest

from validation.junction.models.mpce_v1_network import MPCEv1Network
from validation.junction.models.tee_junction_element_network import (
    TeeJunctionElementNetwork,
)
from validation.junction.network_runner import (
    build_network_cells,
    iter_network_records,
    run_network,
)
from validation.junction.schema import load_dataset


def test_tee_network_adapter_evaluates_K6_case():
    """Tee adapter converges at a known mid-range K6 case (psi=1, theta=90, q=0.5)."""
    r = TeeJunctionElementNetwork().evaluate_network("bassett2001", "K6", 0.5, 1.0, math.pi / 2.0)
    assert r.converged, f"Tee network must converge at K6/psi=1/theta=90/q=0.5; msg={r.message}"
    assert r.K_lateral is not None
    # Loose bound: K should be within an order of magnitude of Bassett K6 ~ 0.87.
    assert 0.0 < r.K_lateral < 5.0, f"K_lateral {r.K_lateral} out of plausible range"


def test_mpce_network_adapter_evaluates_K6_case():
    """MPCE adapter converges at the same K6 case."""
    r = MPCEv1Network().evaluate_network("bassett2001", "K6", 0.5, 1.0, math.pi / 2.0)
    assert r.converged, f"MPCE network must converge at K6/psi=1/theta=90/q=0.5; msg={r.message}"
    assert r.K_lateral is not None
    assert 0.0 < r.K_lateral < 5.0


def test_mpce_K5_is_evaluated_not_declared_unsupported():
    """K_straight is a real prediction and must be scored, not skipped.

    MPCE-v1 used to hard-code ``converged=False`` for K5/K2 on the reasoning
    that ``sin^2(theta=0)`` loses the axial coupling. That assumption was never
    measured, so the straight-flow coefficient carried no error number at all
    for the whole of the junction arc. Whatever the model predicts, the harness
    has to produce a value the scorecard can hold against Bassett Fig 7c.
    """
    r = MPCEv1Network().evaluate_network("bassett2001", "K5", 0.5, 1.0, math.pi / 2.0)
    assert r.converged, f"K5 must be evaluated, not skipped: {r.message}"
    assert r.K_straight is not None


def test_mpce_straight_K_uses_the_straight_leg_fraction():
    """K5's q is the straight fraction; K6's is the lateral one.

    Bassett indexes each coefficient by the mass-flow fraction in its own leg,
    so building the network from a K5 file's q without the 1-q transform picks
    a mirrored operating point. Pinned by symmetry: evaluating K5 at q and K6
    at 1-q must land on the SAME network, hence the same extracted pair.
    """
    model = MPCEv1Network()
    straight = model.evaluate_network("bassett2001", "K5", 0.3, 1.0, math.pi / 2.0)
    lateral = model.evaluate_network("bassett2001", "K6", 0.7, 1.0, math.pi / 2.0)

    assert straight.converged and lateral.converged
    assert straight.K_straight == pytest.approx(lateral.K_straight, rel=1e-9)
    assert straight.K_lateral == pytest.approx(lateral.K_lateral, rel=1e-9)


def test_network_runner_produces_records_and_scorecard():
    """End-to-end: TeeJunctionElement run produces records across multiple
    topologies and the network-cell scorecard groups by them."""
    records = run_network(TeeJunctionElementNetwork())
    assert len(records) > 300, f"too few records: {len(records)}"
    topologies = {r.topology for r in records}
    assert topologies == {"imposed_q", "three_pb", "mfb_two_pb"}, (
        f"expected all 3 topologies, got {topologies}"
    )

    cells = build_network_cells(records)
    # Cells are now grouped by topology too.
    assert any(c.topology == "three_pb" for c in cells)
    # K_lateral_sep at psi=1, theta=45 should fully converge under imposed_q.
    matches = [
        c
        for c in cells
        if c.canonical_K == "K_lateral_sep"
        and c.psi_bin == "1"
        and c.theta_bin == "45"
        and c.topology == "imposed_q"
    ]
    assert matches and matches[0].pct_converged > 0.9, (
        "K_lateral_sep imposed_q psi=1 theta=45 should converge cleanly"
    )


def test_topology_filter_in_run_network():
    """Caller can restrict topologies via the `topologies` kwarg."""
    records = run_network(TeeJunctionElementNetwork(), topologies=("imposed_q",))
    topologies = {r.topology for r in records}
    assert topologies == {"imposed_q"}


def test_iter_skips_x_axis_M3_files():
    """Wang files (x_axis=M_3) should be skipped -- network mode uses q only."""
    ds = load_dataset()
    tee = TeeJunctionElementNetwork()
    records = list(iter_network_records(tee, ds, topologies=("imposed_q",)))
    papers = {r.paper for r in records}
    assert "wang2014" not in papers, "wang files should be skipped (M_3 axis)"
