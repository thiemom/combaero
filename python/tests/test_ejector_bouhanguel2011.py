"""Bouhanguel, Desevaux & Gavignet 2011: the anchor pairs against the flow
regime the authors saw (validation/ejector/bouhanguel2011_regime_check.py).
Accuracy, cross-source: neither closure was built on this rig."""

from __future__ import annotations

from validation.ejector.bouhanguel2011_regime_check import regime_rows


def test_huangs_critical_pressure_agrees_with_every_observed_regime() -> None:
    rows = [r for r in regime_rows() if r.closure == "huang1999"]
    assert len(rows) == 6
    assert all(r.consistent for r in rows)
    # Not at an edge: "not choked" holds with P_c* clearly below Pa at P2_min.
    assert all(r.pc_over_pa_lo < 0.8 for r in rows if r.stated == "subcritical")


def test_the_lossless_critical_pressure_needs_an_edge_case_reading() -> None:
    """K&D contradicts "not choked" under the gauge reading, and under the
    absolute one sits on the edge at P1 = 5 bar (P_c* = Pa at P2_min)."""
    rows = {
        (r.reading, r.p1_bar): r
        for r in regime_rows()
        if r.closure == "kracik_dvorak" and r.stated == "subcritical"
    }
    assert not rows[("gauge", 5.0)].consistent
    assert not rows[("gauge", 4.5)].consistent
    assert abs(rows[("absolute", 5.0)].pc_over_pa_lo - 1.0) < 0.01


def test_both_pairs_call_the_small_tube_choked() -> None:
    assert all(r.consistent for r in regime_rows() if r.stated == "critical")
