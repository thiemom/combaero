"""Cross-source accuracy of the ejector's anchors (validation/ejector/anchor_scorecard.py).

The scorecard is a REPORT: the numbers are the model's limitation on rigs that
are not its sources, never a target. These tests pin what it reports, so that
a change to the closures that moves it is seen, and the README's statement of
the limitation is updated with it.
"""

from __future__ import annotations

import pytest

from validation.ejector.anchor_scorecard import report, score


@pytest.fixture(scope="module")
def rows():
    return score()


def test_both_rigs_are_scored_and_labelled(rows) -> None:
    kinds = {(r.source.split(" ")[0], r.kind) for r in rows}
    assert kinds == {("Henry", "measured"), ("Akbarnejad", "CFD")}
    assert len(rows) == 6 + 10  # Henry 3 pressures x 2 P1 readings; 10 CFD rows
    assert report(rows).startswith("CROSS-SOURCE ACCURACY")


def test_the_dead_head_anchor_runs_far_above_every_reference_breakdown(rows) -> None:
    """The documented limitation (README "Cross-source accuracy of the
    anchors"): P_b0, Kracik & Dvorak's mixing at omega = 0, is lossless and
    independent of A_3/A_t, so the subcritical droop window P_b0 / P_c* is
    1.4-1.8 where the rigs show 1.06-1.32, and 3.2-3.4 at A_3/A_t = 33. On
    every row P_b0 is at least 78% high."""
    for r in rows:
        assert r.pb0_err > 0.75, (r.source, r.case)
        assert r.window_model > r.window_ref + 0.3, (r.source, r.case)
    big = [r for r in rows if r.area_ratio_mix > 30.0]
    assert big and all(r.window_model > 3.0 for r in big)


def test_the_critical_anchor_is_high_but_far_closer(rows) -> None:
    """P_c* (the plateau's end) lands 18-63% high: above, as a lossless
    closure must, and well short of P_b0's error."""
    for r in rows:
        assert 0.15 < r.pc_err < 0.65, (r.source, r.case)
        assert r.pc_err < r.pb0_err


def test_the_sudden_expansion_candidate_lands_near_every_breakdown(rows) -> None:
    """The candidate dead-head (not in the model): the primary's nozzle exit
    into the mixing tube, the dead-headed suction annulus at P_e. It is 3-39%
    high on every row (today's anchor: 78-274%). Its window P_b0 / P_c* is
    0.98-1.06, because P_c* itself is high: the README's reason not to adopt
    it alone."""
    for r in rows:
        assert 0.0 < r.pb0_candidate_err < 0.40, (r.source, r.case)
        assert r.pb0_candidate_err < r.pb0_err / 2.0
        assert 0.95 < r.pb0_candidate / r.pc_model < 1.08
