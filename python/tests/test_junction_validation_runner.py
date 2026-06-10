"""
Smoke tests for the junction-validation runner + scorecard.

These guard the dataset loader, the PaperCeiling adapter, the
TeeJunctionRaw adapter, and the end-to-end record/scorecard pipeline
against silent breakage. They do NOT enforce specific accuracy
targets -- that's the job of the scorecard's headline metrics, which a
user reviews manually.
"""

from __future__ import annotations

import math

from validation.junction.models import bassett2001
from validation.junction.models.paper_ceiling import PaperCeiling
from validation.junction.models.tee_junction_raw import TeeJunctionRaw
from validation.junction.runner import run
from validation.junction.schema import load_dataset
from validation.junction.scorecard import build_cells, headline


def test_dataset_loads():
    ds = load_dataset()
    assert {"bassett2001", "hager1984", "perez_garcia2010", "wang2014"} <= ds.papers.keys()
    assert len(ds.files) >= 80, f"expected ~89 files, got {len(ds.files)}"
    assert ds.papers["bassett2001"].q_convention == "bassett"
    assert ds.papers["hager1984"].q_convention == "hager"


def test_bassett_K12_correct_at_known_reference_point():
    """K12 at (q=0.5, psi=1, theta=pi/2) should equal the C++ gtest reference 0.348.

    This pins the K12 transcription -- catches a regression to the historical
    psi^2 typo that gave -0.72 at psi=3 and 0.348 at psi=1 (hiding the bug)."""
    val = bassett2001.K12_corr(0.5, 1.0, math.pi / 2.0)
    assert abs(val - 0.348) < 5e-3, f"K12_corr at (0.5, 1, pi/2) = {val}, expected ~0.348"


def test_bassett_K12_correct_at_psi3():
    """K12 at (q=0.5, psi=3, theta=pi/4) should be ~1.47, matching Fig 10c measured (1.46).

    The psi^2 transcription would give -0.72; the psi (correct) form gives ~1.47."""
    val = bassett2001.K12_corr(0.5, 3.0, math.pi / 4.0)
    assert 1.3 < val < 1.6, f"K12_corr at psi=3 = {val}, expected ~1.47"


def test_runner_produces_scorecard():
    """End-to-end smoke: TeeJunctionRaw vs PaperCeiling, scored against Bassett data.

    Confirms the records pipeline does not crash and produces a non-empty
    set of evaluable cells. Numeric thresholds are deliberately loose."""
    records = run(TeeJunctionRaw(), PaperCeiling())
    assert len(records) > 500, f"expected ~600+ records, got {len(records)}"

    # At least 200 records have both K_model and K_ceiling (the four bound K's
    # cover plenty of Bassett's dataset).
    n_paired = sum(1 for r in records if r.K_model is not None and r.K_ceiling is not None)
    assert n_paired > 100, f"only {n_paired} records with both model+ceiling values"

    cells = build_cells(records)
    h = headline("tee_junction_raw_cpp", cells)
    assert h.n_supported > 100
    # Core RMSE should be modest -- the K's the C++ supports (K5/K6/K11/K12)
    # are Bassett's own analytical correlations, so model should approximately
    # match the ceiling -> low Delta_vs_ceiling.
    assert h.core_rmse < 1.0, f"core_RMSE = {h.core_rmse}, unexpectedly high"
