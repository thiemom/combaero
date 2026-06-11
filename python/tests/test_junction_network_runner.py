"""Smoke tests for the network-mode validation runner."""

from __future__ import annotations

import math

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


def test_mpce_K5_returns_unconverged_with_diagnostic():
    """MPCE-v1's documented straight-K limitation is reported via converged=False."""
    r = MPCEv1Network().evaluate_network("bassett2001", "K5", 0.5, 1.0, math.pi / 2.0)
    assert not r.converged
    assert "MPCE-v1" in r.message or "K_straight" in r.message


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
