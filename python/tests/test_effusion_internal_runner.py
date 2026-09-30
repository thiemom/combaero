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
from validation.cooling.schema import load_dataset, load_points  # noqa: E402


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


def test_owns_exactly_the_two_effusion_h_figures(plate_c) -> None:
    """Dispatch is an explicit predicate, so it must be exactly right.

    Two figures, two h conventions, one runner -- so the claim list is
    named rather than pattern-matched. Anything else appearing here means
    a series has picked up one of the two `y_axis` values by accident.
    """
    assert er.owns(plate_c)
    claimed = {s.label for s in load_dataset() if er.owns(s)}
    assert claimed == {
        "andrews1988/fig8_h_effusionC",
        "andrews1986/fig8_hm_plate_a",
        "andrews1986/fig8_hm_plate_b",
        "andrews1986/fig8_hm_plate_c",
        "andrews1986/fig8_hm_plate_d",
    }, f"runner claims {sorted(claimed)}"


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


# ---------------------------------------------------------------------------
# Andrews 86-GT-225 Figure 10: the author's own evaluation of his own
# equations. This is an IMPLEMENTATION check -- did we transcribe Eqs. (12)
# to (19) correctly from a poor scan -- which is a different and stronger
# question than Fig. 8's model-vs-data check. The two together separate the
# failure modes #389 exists to keep apart.
# ---------------------------------------------------------------------------

FIG10_RE = 2200.0
FIG10_PR = 0.72
FIG10_L_MM = 6.35
FIG10_X_MM = 6.11  # the 25 x 25 plates, and what Eq. (19) is printed for


def _fig10(name):
    series = next(s for s in load_dataset() if s.label == f"andrews1986/{name}")
    return sorted((p.x, p.y) for p in load_points(series))


def _nu_inf():
    return 0.023 * FIG10_RE**0.8 * FIG10_PR ** (1 / 3)


def test_fig10_throat_series_reproduces_our_entry_length_factor() -> None:
    """The triangles ARE R_Nu, so this checks Eqs. (13) and (14) directly.

    Dividing the throat Nusselt number by the fully developed value leaves
    the entry-length factor alone, so the paper has effectively plotted
    `mills_entry_length_factor` for us. That is the most direct check
    available on two polynomials transcribed from a scan that renders them
    as "RNu 0.13(n) - 0*75(r) + 1.04(r) + 2.24".

    WHAT THIS TEST DOES AND DOES NOT CONSTRAIN, measured by perturbation.
    Fig. 10 spans L/D 4.8 to 9.9, so `w = D/L` runs 0.10 to 0.21 and the
    high-order terms of Eq. (14) barely contribute: changing the `w^3`
    coefficient from 58.6 to 56.6 moves R_Nu under 1% and this test stays
    green. Changing the LINEAR 7.48 to 7.00 turns it red.

    The high-order terms are pinned elsewhere, by
    `MillsBranchesAgreeWhereTheyMeet` in tests/test_effusion_internal.cpp:
    at the L/D = 2 junction `w` is 0.5, where they dominate, and that test
    does catch the 58.6 change. Neither test alone constrains the whole
    polynomial; together they do.
    """
    errs = []
    for LD, y in _fig10("fig10_eq12_throat_only"):
        ours = cb.mills_entry_length_factor(LD)
        errs.append(abs(ours / y - 1.0))
        assert ours == pytest.approx(y, rel=0.03), f"L/D = {LD}"
    assert sum(errs) / len(errs) < 0.02, "mean disagreement grew"

    # And the same numbers must come back through the throat correlation,
    # which is what actually ships.
    for LD, y in _fig10("fig10_eq12_throat_only"):
        assert cb.effusion_throat_nusselt(FIG10_RE, FIG10_PR, LD) / _nu_inf() == (
            pytest.approx(y, rel=0.03)
        )


def test_fig10_eq19_series_reproduces_our_summed_correlation() -> None:
    """Our Eq. (19) against the author's own Eq. (19).

    This is the check that our assembly -- approach rebased by X/(pi L),
    plus throat, both on the hole internal area -- is the assembly the
    paper means. Agreement is 1.5% mean, 2.5% worst, which is digitisation
    precision plus the marker-centre read.
    """
    errs = []
    for LD, y in _fig10("fig10_eq19_summed"):
        ours = (
            cb.effusion_internal_nusselt(FIG10_RE, FIG10_PR, FIG10_X_MM / FIG10_L_MM, LD)
            / _nu_inf()
        )
        errs.append(ours / y - 1.0)
        assert ours == pytest.approx(y, rel=0.04), f"L/D = {LD}"
    assert abs(sum(errs) / len(errs)) < 0.025


def test_fig10_plate_a_is_a_pitch_test_not_an_outlier() -> None:
    """The point that looks 50% wrong is the only one that tests the pitch.

    Plate a is a 10 x 10 array on a 152 mm plate, so its pitch is 15.2 mm
    against 6.08 mm for the 25 x 25 plates -- and Eq. (18)'s approach term
    scales with X/(pi L). Evaluated at plate a's own pitch, Eq. (19)
    reproduces its measurement to under 5%; evaluated at the others' pitch
    it misses by a third.

    The other three points share one X and so cannot distinguish the pitch
    dependence at all, which is what makes this single point worth naming.
    """
    measured = _fig10("fig10_measured")
    plate_a = next(p for p in measured if 5.0 < p[0] < 5.8)
    LD, y = plate_a
    assert y > 4.0, "plate a should be the high point"

    at_own_pitch = (
        cb.effusion_internal_nusselt(FIG10_RE, FIG10_PR, 15.2 / FIG10_L_MM, LD) / _nu_inf()
    )
    at_other_pitch = (
        cb.effusion_internal_nusselt(FIG10_RE, FIG10_PR, FIG10_X_MM / FIG10_L_MM, LD) / _nu_inf()
    )

    assert at_own_pitch == pytest.approx(y, rel=0.05)
    assert at_other_pitch < 0.75 * y, "the wrong pitch should miss badly"


def test_fig10_shows_why_the_approach_term_is_needed() -> None:
    """Throat-only sits far below the measurement, which is the paper's
    whole argument: "the differences emphasise the importance of the hole
    approach flow heat transfer, which is the dominant process at all L/D
    of relevance to gas turbine film cooling applications"."""
    throat = dict(_fig10("fig10_eq12_throat_only"))
    summed = dict(_fig10("fig10_eq19_summed"))
    assert max(throat.values()) < 2.0
    assert min(summed.values()) > 2.3
    for LD_t, y_t in throat.items():
        nearest = min(summed, key=lambda k: abs(k - LD_t))
        assert summed[nearest] > y_t, f"L/D ~ {LD_t}"


def test_fig10_calibration_ticks_land_on_round_numbers() -> None:
    """The digitisation's own evidence, as for murray2018.

    Not declared as a series: `kind: frame` means a straight panel EDGE
    fitted as a power law, and this is a traversal -- left origin, x ticks
    5 to 10, top-right corner, top of the y axis, then y ticks 4 down to
    1.5. It is therefore NOT sorted, and that order is its meaning.
    """
    import csv

    path = (
        Path(__file__).resolve().parents[2]
        / "validation"
        / "cooling"
        / "data"
        / "andrews1986"
        / "fig10_calibration.csv"
    )
    rows = [
        (float(r[0]), float(r[1]))
        for r in csv.reader(path.open())
        if not r[0].strip().startswith("x")
    ]
    assert len(rows) == 13

    origin, x_ticks, top_right, y_top, y_ticks = (rows[0], rows[1:7], rows[7], rows[8], rows[9:])
    x_span = top_right[0] - origin[0]
    y_span = y_top[1] - origin[1]

    worst_x = max(abs(p[0] - t) for p, t in zip(x_ticks, [5, 6, 7, 8, 9, 10], strict=True))
    worst_y = max(abs(p[1] - t) for p, t in zip(y_ticks, [4, 3, 2, 1.5], strict=True))
    assert worst_x / x_span < 0.01, f"x ticks off by {worst_x:.3f}"
    assert worst_y / y_span < 0.01, f"y ticks off by {worst_y:.4f}"

    # The ticks must sit on their axes, which is what makes them ticks.
    assert max(abs(p[1] - origin[1]) for p in x_ticks) < 0.01
    assert max(abs(p[0] - origin[0]) for p in y_ticks) < 0.02
    # And the frame must close.
    assert abs(top_right[1] - y_top[1]) < 0.01


# ---------------------------------------------------------------------------
# Andrews 86-GT-225 Figure 8: h_m against G for the four test plates.
#
# The paper that OWNS the correlations, so this is FIDELITY where
# `andrews1988` Fig. 8 is accuracy -- and the two report as separate rows.
#
# The trap this block exists to hold shut: 86-GT-225's `h_m` is "over the
# hole length" and 88-GT-290's `h` is "based on the surface area, A". Same
# author, same group, same figure number, different quantity.
# ---------------------------------------------------------------------------

FIG8_PLATES = ("a", "b", "c", "d")

# Table 5, with D recovered as X/(X/D) because the scan loses the D column.
# (D_mm, X_mm, L_over_D, X_over_D, technique)
TABLE_5 = {
    "a": (1.178, 15.20, 5.38, 12.9, "Drilled"),
    "b": (0.642, 6.10, 9.92, 9.5, "Spark Eroded"),
    "c": (0.897, 6.10, 7.05, 6.8, "Drilled"),
    "d": (1.298, 6.10, 4.85, 4.7, "Drilled"),
}


@pytest.fixture(scope="module")
def fig8_1986():
    series = {
        s.path.stem[-1]: s
        for s in load_dataset()
        if s.label.startswith("andrews1986/fig8_hm_plate_")
    }
    assert set(series) == set(FIG8_PLATES)
    return series


def test_the_two_papers_h_conventions_are_different_quantities(fig8_1986) -> None:
    """The one error this dataset can make that looks like a model result.

    86-GT-225 nomenclature: "h_m Average heat transfer coefficient, W/m2K,
    over the hole length" -- the hole internal area, pi D L, which is
    already the basis of Eq. (19)'s Nusselt number.
    88-GT-290 nomenclature: "h Convective heat transfer coefficient based
    on the surface area, A", with "A = X^2 - pi/4 D^2" -- the plate area.

    The ratio is A/A_h, between 1.4 and 9.8 over these four plates. Scored
    on the wrong convention the same points read -60%; on the right one
    they read -10%. Nothing in the numbers themselves says which is meant,
    which is why the conventions are distinct `y_axis` values and not a
    flag.
    """
    for tag, series in fig8_1986.items():
        assert series.y_axis == er.H_HOLE_LENGTH
        D, X, L = er._geometry(series)
        ratio = (X * X - math.pi * D * D / 4.0) / (math.pi * D * L)

        on_hole = er.predict(series, 0.5)
        on_plate = er.predict(dataclasses.replace(series, y_axis=er.H_PLATE_AREA), 0.5)
        assert on_plate / on_hole == pytest.approx(1.0 / ratio), tag

    ratios = []
    for series in fig8_1986.values():
        D, X, L = er._geometry(series)
        ratios.append((X * X - math.pi * D * D / 4.0) / (math.pi * D * L))
    assert min(ratios) == pytest.approx(1.4, abs=0.1)
    assert max(ratios) == pytest.approx(9.8, abs=0.2)

    # And on the wrong convention the whole set collapses, which is the
    # number that makes the distinction impossible to wave away.
    wrong = [
        er.predict(dataclasses.replace(s, y_axis=er.H_PLATE_AREA), p.x) / p.y - 1.0
        for s in fig8_1986.values()
        for p in load_points(s)
    ]
    assert sum(wrong) / len(wrong) < -0.5, "the wrong convention no longer collapses"


def test_table_5_geometry_is_internally_consistent(fig8_1986) -> None:
    """The scan loses Table 5's n and D columns; X/D recovers them.

    Two independent checks that the recovery is right, neither of which
    uses the correlation: the recovered D must reproduce the table's own
    printed L/D at the common 6.35 mm thickness, and the four L/D values
    must land on Figure 10's four abscissae -- a different figure, read in
    a different session.
    """
    fig10_abscissae = [x for x, _ in _fig10("fig10_measured")]
    assert len(fig10_abscissae) == 4

    for tag, (D_mm, X_mm, L_over_D, X_over_D, _) in TABLE_5.items():
        series = fig8_1986[tag]
        assert series.geometry["D_mm"] == pytest.approx(D_mm)
        assert series.geometry["X_mm"] == pytest.approx(X_mm)

        # D recovered from the pitch column alone.
        assert X_mm / X_over_D == pytest.approx(D_mm, rel=0.01), tag
        # which then reproduces the printed L/D.
        assert 6.35 / D_mm == pytest.approx(L_over_D, rel=0.01), tag
        # and lands on a Figure 10 point.
        nearest = min(fig10_abscissae, key=lambda v: abs(v - L_over_D))
        assert nearest == pytest.approx(L_over_D, rel=0.02), tag

    # Plate b is the geometry Eq. (19)'s printed 0.27 is evaluated for.
    assert 0.881 * TABLE_5["b"][1] / (math.pi * 6.35) == pytest.approx(0.27, abs=0.001)


def test_reproduces_andrews_own_figure_eight(fig8_1986) -> None:
    """The FIDELITY acceptance check: the correlations against the data in
    the paper that publishes them.

    Measured: -10.4% bias, 11.1% MAE over 41 points. Per plate,
    a -6.2%, b -19.6%, c -5.5%, d -12.5%.

    NOT one-sided -- the worst positive point is +3.4% -- so unlike the
    1988 comparison this is a strong negative bias with a little scatter
    rather than a pure offset. The bounds below are deliberately loose
    enough to be a binding check on the CHAIN (geometry, the two
    conversions, the correlation call) and not a pin on the digitisation's
    last digit.

    WHAT THIS TEST DOES NOT CATCH, measured rather than assumed: moving
    Sparrow's 0.881 to 0.960 (+9%) leaves it GREEN, because the band is
    wide by design. The coefficient is pinned instead by
    `test_the_fidelity_miss_is_not_a_transcription_error` and
    `test_fig10_eq19_series_reproduces_our_summed_correlation`, both of
    which go red on that change. Neither test constrains the other's
    question; both are needed.
    """
    per_plate = {}
    pooled = []
    for tag, series in fig8_1986.items():
        errs = [r.rel_error for r in er.run_series(series) if r.rel_error is not None]
        assert errs, tag
        per_plate[tag] = sum(errs) / len(errs)
        pooled += errs

    assert len(pooled) == 41
    bias = sum(pooled) / len(pooled)
    mae = sum(abs(e) for e in pooled) / len(pooled)
    assert -0.16 < bias < -0.05, f"pooled bias moved to {bias:.1%}"
    assert mae < 0.16, f"pooled MAE {mae:.1%}"

    for tag, expected in (("a", -0.062), ("b", -0.196), ("c", -0.055), ("d", -0.125)):
        assert per_plate[tag] == pytest.approx(expected, abs=0.04), tag

    # Every plate under-predicts on average, and b is the worst. Recorded
    # because b is also the only SPARK ERODED plate in Table 5 -- an
    # observation, not a demonstrated cause. One plate is not evidence, and
    # nothing here corrects for it.
    assert all(v < 0 for v in per_plate.values())
    assert min(per_plate, key=per_plate.get) == "b"
    assert TABLE_5["b"][4] == "Spark Eroded"
    assert [TABLE_5[t][4] for t in ("a", "c", "d")] == ["Drilled"] * 3


def test_the_fidelity_miss_is_not_a_transcription_error(fig8_1986) -> None:
    """The reason this dataset was worth digitising at all.

    A 10% miss on the author's OWN data normally means our transcription is
    wrong. Here it cannot: Figure 10 of the SAME paper plots Andrews' own
    evaluation of his own Eq. (19), and combaero reproduces it to 1.5%. So
    the -10.4% is the correlation missing its author's measurements, not us
    missing the correlation.

    That separation is the whole point of the fidelity/accuracy split, and
    it only holds if Fig. 10's single Reynolds number sits INSIDE Fig. 8's
    range rather than at an end of it -- so that is asserted, not assumed.
    """
    # The transcription bound, restated here as the premise this rests on.
    errs = [
        cb.effusion_internal_nusselt(FIG10_RE, FIG10_PR, FIG10_X_MM / FIG10_L_MM, LD)
        / _nu_inf()
        / y
        - 1.0
        for LD, y in _fig10("fig10_eq19_summed")
    ]
    assert abs(sum(errs) / len(errs)) < 0.025

    air = cb.standard_dry_air_composition()
    mu = cb.viscosity(er.ASSUMED_T, er.ASSUMED_P, air)
    below = above = 0
    for series in fig8_1986.values():
        D, X, _ = er._geometry(series)
        for p in load_points(series):
            Re = 4.0 * p.x * X * X / (math.pi * D * mu)
            if Re < FIG10_RE:
                below += 1
            else:
                above += 1
    assert below >= 10 and above >= 10, (
        f"Fig. 10's Re = {FIG10_RE:.0f} does not bracket Fig. 8's range "
        f"({below} below, {above} above); the transcription check no longer "
        "covers the regime the fidelity miss is measured in"
    )


def test_fig8_calibration_proves_the_log_axes_were_read_linearly(fig8_1986) -> None:
    """The digitisation's own evidence, kept UNCONVERTED on purpose.

    `fig8_calibration_as_read.csv` holds the frame traversal in the
    digitiser's linear coordinates, taken before the log axes were noticed.
    It is committed raw because it is the proof: a linear pick of a log
    axis puts the 0.5 and 1.0 x ticks at 1.1218 and 1.5619, and they read
    1.1299 and 1.5689.

    Both directions are asserted -- that reading it linearly FAILS and the
    log inverse SUCCEEDS -- because only the pair distinguishes "the
    transform is right" from "the ticks are near the anchors anyway".
    """
    import csv

    path = (
        Path(__file__).resolve().parents[2]
        / "validation"
        / "cooling"
        / "data"
        / "andrews1986"
        / "fig8_calibration_as_read.csv"
    )
    rows = [
        (float(r[0]), float(r[1]))
        for r in csv.reader(path.open())
        if not r[0].strip().startswith("x")
    ]
    assert len(rows) == 9
    origin, x01, x05, x10, x20, top_right, top_left, y1000, y500 = rows

    def undo(v, l0, l1, t0, t1):
        return t0 * (t1 / t0) ** ((v - l0) / (l1 - l0))

    # Anchored on the frame; the two INTERIOR ticks are the free check.
    for read, true in ((x05[0], 0.5), (x10[0], 1.0)):
        assert undo(read, x01[0], x20[0], 0.1, 2.0) == pytest.approx(true, rel=0.02)
    for read, true in ((y500[1], 500.0), (y1000[1], 1000.0)):
        assert undo(read, origin[1], top_right[1], 100.0, 4000.0) == pytest.approx(true, rel=0.02)

    # Read linearly, the SAME ticks land nowhere near their labels -- which
    # is what makes the log reading a finding rather than a preference.
    def linear(v, l0, l1, t0, t1):
        return t0 + (t1 - t0) * (v - l0) / (l1 - l0)

    assert linear(x05[0], x01[0], x20[0], 0.1, 2.0) == pytest.approx(1.13, abs=0.02)
    assert linear(y500[1], origin[1], top_right[1], 100.0, 4000.0) == pytest.approx(
        1790.0, abs=20.0
    )

    # Page skew, measured and recorded rather than corrected: 0.21% of
    # value between the two top corners, which moves a mid-range point by
    # 0.39% -- an order of magnitude below the residual being reported.
    assert abs(top_left[1] / top_right[1] - 1.0) < 0.005

    # The committed series must be the CONVERTED ones, not the raw picks.
    spans = [p.x for s in fig8_1986.values() for p in load_points(s)]
    assert 0.08 < min(spans) < 0.1 and 1.7 < max(spans) < 1.85, (
        "the committed abscissae are not on the converted (log) scale; "
        'the paper states G ran "from 0.1 to 1.7 kg/sm2"'
    )
