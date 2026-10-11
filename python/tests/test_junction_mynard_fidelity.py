"""The junction closure against its own paper, Mynard & Valen-Sendstad 2015.

FIDELITY (a miss is our bug): the port in Mynard's configuration must draw his
red Unified0D curves, Figs 4 and 6-11, 31 panels. Measured 2026-10-10 with
the model evaluated at his CFD samples and splined between them, as his plots
are: 29 panels within 2.7 px of every visible red point (<= 0.032 in K), Fig 4
4.2 px (0.021) and Fig 6d 4.1 px (1.5 on its -20..40 axis; 1.75 px with
pchip). Fig 4's dashed eta_j = 0 curve: 2.2 px. The scans are 150 dpi, so a
pixel is 0.005-0.017 in K on the ordinary panels.

Each production departure from his model shows on this metric, which is what
makes it a check (falsify the metric): the dividing-streamline recovery moves
5 panels by up to 31 px, eta off 12 panels by up to 79 px, the joining alpha
exactly the two psi != 1 joining panels by 19 px.

REF3D (Mynard's own laminar CFD, Re 350-2400) is reported, not asserted
against a target: our production closure departs from his model on purpose
and nothing is tuned on his envelope. The numbers are pinned as regression
values so a change to the closure shows up here as a diff.

Also here: Bassett 2001's own calculated curves against our transcription of
his Table 2 -- the fidelity half of the Bassett source, which had none.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from validation.junction import mynard_fidelity as mf
from validation.junction.models import bassett2001 as bassett

_SPEC = {s.key: s for s in mf.PANELS}
_LOOSE = {"fig04": 4.5, "fig06d": 4.5}  # see module docstring


@pytest.fixture(scope="module")
def faithful() -> dict[str, mf.PanelFidelity]:
    return {s.key: mf.panel_fidelity(s) for s in mf.PANELS}


def test_every_panel_is_digitised(faithful):
    """Marker counts read off the figures by eye; red points must exist on
    every panel (the fewest are on Re panels where Ref3D hides the line)."""
    counts = {p.panel: 0 for p in mf.samples()}
    for p in mf.samples():
        counts[p.panel] += 1
    assert counts["fig04"] == 7
    assert counts["fig08c"] == 10
    assert counts["fig11c"] == 8
    assert sum(counts.values()) == 158
    assert min(f.n_red for f in faithful.values()) >= 15


def test_markers_sit_on_design_abscissae():
    """Every marker lands within 2.5 px of a design value after the marker
    recalibration -- if not, the snap would pick the wrong sample."""
    assert max(p.snap_px for p in mf.samples()) < 2.5


@pytest.mark.parametrize("key", list(_SPEC))
def test_port_draws_mynards_unified0d(faithful, key):
    assert faithful[key].max_px < _LOOSE.get(key, 3.0)


def test_port_draws_the_eta_zero_curve():
    """Fig 4's dashed line is his model with the energy-transfer factor off:
    a direct check of Eq 36 and of eta_scale."""
    f = mf.panel_fidelity(_SPEC["fig04"], "unified0d_eta0", dict(mf.FAITHFUL, eta_scale=0.0))
    assert f.max_px < 3.0


@pytest.mark.parametrize(
    ("knob", "panels", "at_least_px"),
    [
        (
            {"dividing_streamline_recovery": 0.5},
            {"fig04", "fig06e", "fig06f", "fig06g", "fig06h"},
            15.0,
        ),
        (
            {"eta_scale": 0.0},
            {"fig04", "fig06e", "fig06f", "fig06g", "fig06h", "fig07c", "fig08c"},
            10.0,
        ),
        ({"joining_etransfer_alpha": 0.2}, {"fig09d", "fig09h"}, 15.0),
    ],
)
def test_each_departure_is_visible(faithful, knob, panels, at_least_px):
    """Falsify the metric: each knob, flipped to production, moves the panels
    it acts on far outside the fidelity bound and leaves every other one
    where it was."""
    cfg = dict(mf.FAITHFUL, **knob)
    for spec in mf.PANELS:
        moved = mf.panel_fidelity(spec, config=cfg).max_px
        if spec.key in panels:
            assert moved > at_least_px, spec.key
        elif "joining_etransfer_alpha" in knob:
            assert moved == pytest.approx(faithful[spec.key].max_px, abs=1e-9), spec.key


def test_fig9_rows_are_swapped_against_the_caption(faithful):
    """The caption puts the side inlet on the top row. Read that way the port
    misses by tens of pixels; read the other way it draws both rows. The
    physics agrees with the swap: a branch whose flow vanishes has K -> -1,
    and on the figure that is the bottom row as lam -> 0."""
    as_captioned = mf.PanelSpec("fig09f", 4, "lam", "straight")
    assert mf.panel_fidelity(as_captioned).max_px > 30.0
    assert faithful["fig09f"].max_px < 3.0
    side = _SPEC["fig09f"]
    assert side.branch == "side"
    assert mf.model_K(side, 0.05) < -0.7


def test_ref3d_scores_are_pinned():
    """CROSS-CHECK on the author's own CFD, laminar Re 350-2400. Measured
    2026-10-10 (MAE faithful / production): all 158 samples 0.212 / 0.219;
    Fig 4 0.022 / 0.126 (his eta was fitted there); Y diverging 0.083 / 0.124;
    the psi != 1 joining samples 0.198 / 0.202 (bias -0.065 / +0.077).
    His data do not discriminate the joining alpha, and nothing is tuned on
    them."""
    groups = dict(mf.ref3d_groups(mf.samples()))

    def mae(sel: list, attr: str) -> float:
        return float(np.mean([abs(getattr(p, attr) - p.K_ref3d) for p in sel]))

    assert mae(groups["all"], "K_faithful") == pytest.approx(0.212, abs=0.005)
    assert mae(groups["all"], "K_production") == pytest.approx(0.219, abs=0.005)
    assert mae(groups["Y diverging, Fig 4"], "K_faithful") == pytest.approx(0.022, abs=0.005)
    assert mae(groups["Y diverging, Fig 4"], "K_production") == pytest.approx(0.126, abs=0.005)
    assert mae(groups["joining, psi != 1"], "K_production") == pytest.approx(0.202, abs=0.005)


# ---------------------------------------------------------------------------
# Bassett 2001: his own calculated curves against our Table 2
# ---------------------------------------------------------------------------

_BASSETT = Path(__file__).resolve().parents[2] / "validation" / "junction" / "data" / "bassett2001"


@pytest.mark.parametrize(
    ("name", "fn", "psi", "theta", "bound"),
    [
        ("bassett_fig10a_K12_theta=30_psi=1_calc.csv", bassett.K12_raw, 1, 30, 0.03),
        ("bassett_fig10a_K12_theta=45_psi=1_calc.csv", bassett.K12_raw, 1, 45, 0.03),
        ("bassett_fig10a_K12_theta=60_psi=1_calc.csv", bassett.K12_raw, 1, 60, 0.03),
        ("bassett_fig10a_K12_theta=90_psi=1_calc.csv", bassett.K12_raw, 1, 90, 0.03),
        ("bassett_fig10b_K12_theta=90_psi=1_calc.csv", bassett.K12_corr, 1, 90, 0.03),
        ("bassett_fig10c_K12_theta=45_psi=1_calc.csv", bassett.K12_corr, 1, 45, 0.03),
        ("bassett_fig10c_K12_theta=45_psi=3_calc.csv", bassett.K12_corr, 3, 45, 0.11),
        ("bassett_fig10c_K12_theta=45_psi=4_calc.csv", bassett.K12_corr, 4, 45, 0.11),
        ("bassett_fig14_K9_theta=90_psi=1_calc_raw.csv", bassett.K9_raw, 1, 90, 0.05),
        ("bassett_fig14_K9_theta=90_psi=1_calc_corr.csv", bassett.K9_corr, 1, 90, 0.03),
        ("bassett_fig07a_K6_theta=120_psi=1_calc.csv", bassett.K6, 1, 120, 0.035),
        ("bassett_fig07a_K6_theta=90_psi=1_calc.csv", bassett.K6, 1, 90, 0.035),
        ("bassett_fig07a_K6_theta=60_psi=1_calc.csv", bassett.K6, 1, 60, 0.035),
        ("bassett_fig07a_K6_theta=45_psi=1_calc.csv", bassett.K6, 1, 45, 0.035),
        ("bassett_fig07b_K6_theta=45_psi=1_calc.csv", bassett.K6, 1, 45, 0.035),
        ("bassett_fig07b_K6_theta=45_psi=3_calc.csv", bassett.K6, 3, 45, 0.035),
        ("bassett_fig07c_K5_calc.csv", bassett.K5, 1, 45, 0.035),
        ("bassett_fig08_K4_theta=135_psi=1_calc.csv", bassett.K4, 1, 135, 0.035),
        ("bassett_fig08_K4_theta=90_psi=1_calc.csv", bassett.K4, 1, 90, 0.035),
        ("bassett_fig08_K4_theta=45_psi=1_calc.csv", bassett.K4, 1, 45, 0.035),
        ("bassett_fig11a_K11_theta=90_psi=1_calc.csv", bassett.K11_corr, 1, 90, 0.035),
        ("bassett_fig11a_K11_theta=60_psi=1_calc.csv", bassett.K11_corr, 1, 60, 0.035),
        ("bassett_fig11a_K11_theta=45_psi=1_calc.csv", bassett.K11_corr, 1, 45, 0.035),
        ("bassett_fig11a_K11_theta=30_psi=1_calc.csv", bassett.K11_corr, 1, 30, 0.035),
        ("bassett_fig11b_K11_theta=45_psi=1_calc.csv", bassett.K11_corr, 1, 45, 0.035),
        ("bassett_fig11b_K11_theta=45_psi=2_calc.csv", bassett.K11_corr, 2, 45, 0.035),
        ("bassett_fig11b_K11_theta=45_psi=3_calc.csv", bassett.K11_corr, 3, 45, 0.035),
        ("bassett_fig11b_K11_theta=45_psi=4_calc.csv", bassett.K11_corr, 4, 45, 0.035),
        ("bassett_fig12a_K8_theta=90_psi=1_calc.csv", bassett.K8_corr, 1, 90, 0.035),
        ("bassett_fig12a_K8_theta=60_psi=1_calc.csv", bassett.K8_corr, 1, 60, 0.035),
        ("bassett_fig12a_K8_theta=45_psi=1_calc.csv", bassett.K8_corr, 1, 45, 0.035),
        ("bassett_fig12a_K8_theta=30_psi=1_calc.csv", bassett.K8_corr, 1, 30, 0.035),
        ("bassett_fig12b_K8_theta=45_psi=1_calc.csv", bassett.K8_corr, 1, 45, 0.035),
        ("bassett_fig12b_K8_theta=45_psi=2_calc.csv", bassett.K8_corr, 2, 45, 0.035),
        ("bassett_fig12b_K8_theta=45_psi=4_calc.csv", bassett.K8_corr, 4, 45, 0.035),
        ("bassett_fig13a_K7_theta=90_psi=1_calc.csv", bassett.K7_corr, 1, 90, 0.035),
        ("bassett_fig13a_K7_theta=60_psi=1_calc.csv", bassett.K7_corr, 1, 60, 0.035),
        ("bassett_fig13a_K7_theta=45_psi=1_calc.csv", bassett.K7_corr, 1, 45, 0.035),
        ("bassett_fig13a_K7_theta=30_psi=1_calc.csv", bassett.K7_corr, 1, 30, 0.035),
        ("bassett_fig13b_K7_theta=45_psi=1_calc.csv", bassett.K7_corr, 1, 45, 0.035),
        ("bassett_fig13b_K7_theta=45_psi=2_calc.csv", bassett.K7_corr, 2, 45, 0.035),
        ("bassett_fig13b_K7_theta=45_psi=4_calc.csv", bassett.K7_corr, 4, 45, 0.035),
    ],
)
def test_bassett_table_2_draws_his_curves(name, fn, psi, theta, bound):
    """Measured 2026-10-10: max |error| 0.012-0.024 at psi = 1 (K9 raw
    0.041), 0.097 / 0.082 at psi = 3 / 4, where the Fig 10c axis is
    coarser. Digitised curves, so an upper bound on the transcription error.

    K4-K8 and K11 (Figs 7, 8, 11-13, traced 2026-10-11): the trace snaps to
    ink within 8 px of our formula, so its max error (<= 0.031) is bounded by
    construction; the evidence that the formulas are his is the coverage and
    the unbiased ink distance (median 0.5-3 px), against 6-21 px for a
    re-planted K8 /psi typo or a raw-for-corrected K11. As committed data the
    files still catch any later change to these formulas."""
    rows = np.loadtxt(_BASSETT / name, delimiter=",", comments="#", ndmin=2, skiprows=1)
    err = [fn(q, psi, math.radians(theta)) - k for q, k in rows]
    assert max(map(abs, err)) < bound
