"""Scoring effusion internal heat transfer against Andrews' Figure 8.

Ground truth for the correlations themselves is in
`tests/test_effusion_internal.cpp`, against the paper's algebra. What this
file guards is the bookkeeping between the data and those correlations --
the two unit conversions the source's own axes demand, which is where a
scoring harness quietly lies.
"""

from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import combaero as cb  # noqa: E402
from validation.cooling import effusion_internal_runner as er  # noqa: E402
from validation.cooling.schema import load_dataset  # noqa: E402


@pytest.fixture(scope="module")
def andrews():
    series = [
        s for s in load_dataset() if s.label.startswith("andrews1988/") and s.kind == "measured"
    ]
    assert series, "andrews1988 is not in the dataset"
    return series


@pytest.fixture(scope="module")
def plate_c(andrews):
    return next(s for s in andrews if s.label.endswith("fig8_h_effusionC"))


def test_owns_only_the_internal_h_series(plate_c) -> None:
    """Dispatch is an explicit predicate, so it must be exactly right."""
    assert er.owns(plate_c)
    claimed = [
        s.label for s in load_dataset() if not s.label.startswith("andrews1988/") and er.owns(s)
    ]
    assert not claimed, f"runner claims series it does not score: {claimed}"


def test_hole_length_is_derived_from_the_ejection_angle(plate_c) -> None:
    """L is the plate thickness only for a normal hole.

    Andrews drills at 90 degrees so the two coincide here, which is exactly
    why the distinction is worth pinning: an angled plate would be silently
    wrong by 1/sin(alpha), and nothing in this dataset would catch it.
    """
    D, X, L = er._geometry(plate_c)
    assert plate_c.alpha_deg == 90.0
    assert pytest.approx(float(plate_c.geometry["thickness_mm"]) * 1e-3) == L

    # At 30 degrees the same plate would give a hole twice as long.
    angled = dataclasses.replace(plate_c, alpha_deg=30.0)
    _, _, L30 = er._geometry(angled)
    assert pytest.approx(L / math.sin(math.radians(30.0))) == L30


def test_the_two_conversions_the_axes_demand(plate_c) -> None:
    """Nu is on the hole internal area; Fig. 8's h is on the plate area.

    `A_h/A` is 1/3.46 for this plate, so dropping it would overstate h by
    that factor -- a definitional error, not a modelling one, and the class
    of mistake #389 exists to catch. Asserted by reconstructing the
    prediction from the correlation plus geometry alone.
    """
    D, X, L = er._geometry(plate_c)
    air = cb.standard_dry_air_composition()
    T, P = er.ASSUMED_T, er.ASSUMED_P
    k = cb.thermal_conductivity(T, P, air)
    mu = cb.viscosity(T, P, air)
    pr = cb.prandtl(T, P, air)

    G = 0.5
    Re = 4.0 * (G * X * X) / (math.pi * D * mu)
    nu = cb.effusion_internal_nusselt(Re, pr, X / L, L / D)
    area_plate = X * X - math.pi * D * D / 4.0
    area_hole = math.pi * D * L

    assert er.predict(plate_c, G) == pytest.approx(nu * k / D * area_hole / area_plate)
    # And the factor really is ~3.46, so its omission is not marginal.
    assert area_plate / area_hole == pytest.approx(3.46, abs=0.02)


def test_G_is_per_gross_plate_area(plate_c) -> None:
    """The coolant per hole is G X^2, not G (X^2 - pi D^2/4).

    The paper's nomenclature separates "plate area" (for G) from "hole
    approach surface area A" (for h), and the hole count settles it:
    N = 4306 per square metre against 1/X^2 = 4305.8.
    """
    _, X, _ = er._geometry(plate_c)
    assert pytest.approx(float(plate_c.geometry["N_per_m2"]), rel=1e-3) == 1.0 / (X * X)


def test_reproduces_andrews_figure_eight_within_the_sources_own_scatter(
    plate_c,
) -> None:
    """The acceptance check, and it is a CROSS-SOURCE one.

    The correlations are Andrews 86-GT-225 and the data is 88-GT-290 --
    same lab, different study, which the validation policy counts as
    cross-source, so a miss is the model's limitation rather than
    automatically our bug.

    Measured: bias -13.5%, RMSE 14.7% over 10 points. One-sided, so it is a
    real systematic under-prediction rather than scatter -- but it sits
    inside the source's own +/-16% spread about its power-law fit, and the
    two worst points are the two the metadata already flags as off-trend in
    the source.
    """
    records = [r for r in er.run_series(plate_c) if r.rel_error is not None]
    assert len(records) == 10

    errs = [r.rel_error for r in records]
    bias = sum(errs) / len(errs)
    rmse = math.sqrt(sum(e * e for e in errs) / len(errs))

    assert -0.20 < bias < -0.05, f"bias moved to {bias:.1%}"
    assert rmse < 0.20, f"RMSE {rmse:.1%}"
    # One-sided: every point under-predicts, which is the shape being
    # recorded. If it ever straddles zero the systematic has gone and this
    # test should be revisited rather than widened.
    assert all(e < 0 for e in errs), "the under-prediction is no longer one-sided"


def test_the_prediction_rises_with_coolant_flow(plate_c) -> None:
    """Monotone in G, which both Nusselt terms require."""
    previous = None
    for G in (0.1, 0.3, 0.6, 1.0, 1.8):
        v = er.predict(plate_c, G)
        assert v is not None and v > 0
        if previous is not None:
            assert v > previous, f"G = {G}"
        previous = v
    assert er.predict(plate_c, 0.0) is None


def test_temperature_assumption_is_reported_not_minimised(plate_c) -> None:
    """The paper states no coolant temperature, so its cost must be visible.

    Reporting it is the point: choosing the temperature that scores best
    would be fitting an unmeasured input to the metric that judges it.
    Across 280-320 K the bias moves -15.3% to -11.7%, so the assumption is
    worth about 4 points and cannot explain the 13.5% gap.
    """
    sens = er.temperature_sensitivity(plate_c)
    assert set(sens) == {280.0, 300.0, 320.0}
    assert er.ASSUMED_T in sens
    spread = max(sens.values()) - min(sens.values())
    assert spread < 0.06, f"temperature moves the result by {spread:.1%}"
    # Warmer coolant helps, but nowhere near enough to close the gap.
    assert all(v < 0 for v in sens.values())


def test_the_runner_does_not_reimplement_the_correlation(plate_c) -> None:
    """Predictions must move when the real correlation moves.

    A runner that inlined the algebra would keep answering after the
    library changed underneath it. R_Nu depends on L/D, so a different
    plate thickness must change the answer through the correlation.
    """
    thicker = dataclasses.replace(plate_c, geometry={**plate_c.geometry, "thickness_mm": 20.0})
    assert er.predict(thicker, 0.5) != er.predict(plate_c, 0.5)
    # A longer hole has a lower entry-length enhancement.
    assert cb.mills_entry_length_factor(20.0 / 3.27) < cb.mills_entry_length_factor(6.3 / 3.27)
