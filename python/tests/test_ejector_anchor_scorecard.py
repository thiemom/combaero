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
    assert kinds == {("Henry", "measured"), ("Chen", "measured"), ("Akbarnejad", "CFD")}
    # Henry 3 pressures x 2 P1 readings; Chen 5 pressures; 10 CFD rows
    assert len(rows) == 6 + 5 + 10
    assert report(rows).startswith("CROSS-SOURCE ACCURACY")


def test_the_dead_head_anchor_runs_far_above_every_reference_breakdown(rows) -> None:
    """The documented limitation (README "Cross-source accuracy of the
    anchors"): P_b0, Kracik & Dvorak's mixing at omega = 0, is lossless and
    independent of A_3/A_t, so the subcritical droop window P_b0 / P_c* is
    1.4-1.8 where the rigs show 1.06-1.32, and 3.2-3.4 at A_3/A_t = 33. On
    every row P_b0 is at least 62% high (Chen's rectangular rig: 62-68%)."""
    for r in rows:
        assert r.pb0_err > 0.60, (r.source, r.case)
        assert r.window_model > r.window_ref + 0.3, (r.source, r.case)
    big = [r for r in rows if r.area_ratio_mix > 30.0]
    assert big and all(r.window_model > 3.0 for r in big)


def test_the_critical_anchor_is_high_but_far_closer(rows) -> None:
    """P_c* (the plateau's end) lands 14-63% high: above, as a lossless
    closure must, and well short of P_b0's error."""
    for r in rows:
        assert 0.12 < r.pc_err < 0.65, (r.source, r.case)
        assert r.pc_err < r.pb0_err


def test_the_sudden_expansion_candidate_lands_near_every_breakdown(rows) -> None:
    """The candidate dead-head (not in the model): the primary's nozzle exit
    into the mixing tube, the dead-headed suction annulus at P_e. It is 3-39%
    high on every row (today's anchor: 62-274%); on Chen's rectangular rig,
    scored after the candidate was fixed, 3-5%. Its window P_b0 / P_c* is
    0.98-1.06, because P_c* itself is high: the README's reason not to adopt
    it alone."""
    for r in rows:
        assert 0.0 < r.pb0_candidate_err < 0.40, (r.source, r.case)
        assert r.pb0_candidate_err < r.pb0_err / 2.0
        assert 0.95 < r.pb0_candidate / r.pc_model < 1.08


def test_the_candidate_dead_head_holds_on_the_rig_it_was_not_built_on(rows) -> None:
    """Chen's rectangular rig was added after the candidate's form was fixed:
    it lands 3-5% above the extrapolated breakdown on all five curves."""
    chen_rows = [r for r in rows if r.source.startswith("Chen")]
    assert len(chen_rows) == 5
    for r in chen_rows:
        assert 0.0 < r.pb0_candidate_err < 0.07, r.case


def test_the_candidate_pair_keeps_a_positive_window_and_halves_the_curve_error(rows) -> None:
    """The candidate pair (Huang's printed P_c* chain with the sudden-expansion
    P_b0; README "The candidate pair"): its window P_b0/P_c* is 1.37-1.65 on
    every row, A_3/A_t = 33 included, its P_c* lands -19..+17%, and on every
    measured curve the plateau-and-chord RMS is 12-23% against 29-59% today."""
    for r in rows:
        assert 1.3 < r.pb0_candidate / r.pc_candidate < 1.7, (r.source, r.case)
        assert -0.20 < r.pc_candidate_err < 0.20, (r.source, r.case)
    measured = [r for r in rows if r.kind == "measured"]
    assert len(measured) == 11
    for r in measured:
        assert r.curve_rms_candidate < 0.25, (r.source, r.case)
        assert r.curve_rms_candidate < r.curve_rms / 1.4, (r.source, r.case)
