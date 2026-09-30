"""Overall effusion cooling effectiveness against Andrews 88-GT-290 Fig. 10.

This is the last piece of #387, and its result is partly negative. The
recorded plan was internal convection plus an external film, superposed
into an overall eta. Andrews' own data refuses the film half, and these
tests are what hold that refusal in place -- so that a future change which
re-adds a film term goes red rather than quietly looking better on one
plate.
"""

from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import combaero as cb  # noqa: E402
from validation.cooling import effusion_overall_runner as er  # noqa: E402
from validation.cooling.schema import load_dataset, load_points  # noqa: E402


@pytest.fixture(scope="module")
def plates():
    s = {
        (x.geometry or {}).get("plate"): x
        for x in load_dataset()
        if x.label.startswith("andrews1988/fig10_eta_effusion_")
    }
    assert set(s) == {"B", "C"}, f"got {sorted(s)}"
    return s


def _errors(series):
    return [r.rel_error for r in er.run_series(series) if r.rel_error is not None]


def test_owns_exactly_the_overall_eta_series(plates) -> None:
    """Dispatch is an explicit predicate, so it must be exactly right.

    Fig. 10's other six curves carry the same `y_axis` but `scores: null`,
    so they are OWNED and returned unscored -- that is the designed
    behaviour and is asserted here rather than left to chance.
    """
    claimed = {s.label for s in load_dataset() if er.owns(s)}
    assert claimed == {
        "andrews1988/fig10_eta_effusion_B",
        "andrews1988/fig10_eta_effusion_C",
        "andrews1988/fig10_eta_imp_effusion_AB",
        "andrews1988/fig10_eta_imp_effusion_AC",
        "andrews1988/fig10_eta_impingement_A",
        "andrews1988/fig10_eta_porous_wall",
        "andrews1988/fig10_eta_lamilloy",
        "andrews1988/fig10_eta_transply",
    }, f"runner claims {sorted(claimed)}"

    unscored = [
        s for s in load_dataset() if s.label.startswith("andrews1988/fig10_") and s.scores is None
    ]
    assert len(unscored) == 6
    for s in unscored:
        recs = er.run_series(s)
        assert recs and all(r.predicted is None for r in recs)
        assert all("not scored" in (r.reason or "") for r in recs)


def test_the_gas_side_is_the_rig_not_a_fitted_number(plates) -> None:
    """h_g comes from the duct the paper describes, reproducibly.

    Page 3: "the wall of a 76 mm by 152 mm wide air cooled duct through
    which the product gases from a propane preheater flowed at 750K and a
    Mach number of 0.05". Every number below follows from that sentence
    plus combaero's properties -- nothing is chosen.
    """
    g = er.gas_side()
    assert g["D_h"] == pytest.approx(4 * (0.076 * 0.152) / (2 * (0.076 + 0.152)))
    assert g["U_g"] == pytest.approx(26.8, abs=0.3)
    assert g["Re_duct"] == pytest.approx(38400, rel=0.02)
    assert g["h_g"] == pytest.approx(46.1, rel=0.02)
    # Turbulent, so Dittus-Boelter is at least the right family of
    # correlation -- a laminar duct would need a different one entirely.
    assert g["Re_duct"] > 10000

    # The paper states a coolant-to-gas density ratio of "approximately
    # 2.5"; it falls out of the two temperatures, which is a check that
    # the stated conditions are mutually consistent.
    dr = er.jet_regime(plates["C"], 0.5)["density_ratio"]
    assert dr == pytest.approx(2.5, abs=0.1)


def test_scores_both_plates_and_reports_them_separately(plates) -> None:
    """The acceptance check, and the whole point is that it is TWO numbers.

    plate C  +3.1% bias,  4.6% MAE, 72 points
    plate B +23.7% bias, 23.7% MAE, 66 points

    A pooled +13% would describe neither. The runner's reporting group is
    per plate so the rollup cannot form it.
    """
    c, b = _errors(plates["C"]), _errors(plates["B"])
    assert len(c) == 72 and len(b) == 66

    bias_c, bias_b = sum(c) / len(c), sum(b) / len(b)
    assert 0.0 < bias_c < 0.08, f"plate C bias moved to {bias_c:.1%}"
    assert 0.17 < bias_b < 0.30, f"plate B bias moved to {bias_b:.1%}"
    # The gap between them is the finding; it must stay large.
    assert bias_b - bias_c > 0.15

    assert er._group_of(plates["B"]) != er._group_of(plates["C"])
    assert "plate B" in er._group_of(plates["B"])


def test_a_velocity_ratio_split_would_misreport_the_plate_effect(plates) -> None:
    """Why the reporting group is per plate and NOT per jet regime.

    Splitting the points by velocity ratio looks informative -- roughly
    +7% below VR = 1 against +19% above it. It is an artefact of which
    plate supplies the points: both plates cross VR = 1 inside their own G
    range, and VR < 1 is mostly plate C while VR >= 1 is mostly plate B.

    Holding the regime fixed and varying the plate shows it: the plate
    moves the error by a factor of ten, the regime label by a few points.
    A VR-grouped row would report the plate mix and call it physics.
    """
    cells: dict[tuple[str, bool], list[float]] = {}
    for tag, series in plates.items():
        for p in load_points(series):
            vr = er.jet_regime(series, p.x)["velocity_ratio"]
            err = er.predict(series, p.x) / p.y - 1.0
            cells.setdefault((tag, vr >= 1.0), []).append(err)

    assert set(cells) == {("B", False), ("B", True), ("C", False), ("C", True)}
    mean = {k: sum(v) / len(v) for k, v in cells.items()}

    # Same regime, different plate: a factor of ten.
    for regime in (False, True):
        assert mean[("B", regime)] > 3.0 * mean[("C", regime)], (
            f"VR>={regime}: B {mean[('B', regime)]:.1%} vs C {mean[('C', regime)]:.1%}"
        )
    # Same plate, different regime: a few points.
    for tag in ("B", "C"):
        assert abs(mean[(tag, True)] - mean[(tag, False)]) < 0.06, (
            f"plate {tag} now depends strongly on the regime label"
        )


def test_a_film_term_alone_is_refused_without_a_matching_augmentation(
    plates,
) -> None:
    """The negative result this module records -- stated precisely.

    NOT "there is no film". Baldauf's adiabatic effectiveness already
    contains the counter-rotating vortex pair entraining hot gas under a
    lifted jet: that changes the adiabatic wall temperature, which is
    exactly what an IR-on-insulated-wall measurement records, and Eq. (38)
    decays as steeply as `mu^-4.57` for plate B's spacing.

    What an adiabatic experiment CANNOT yield is the coefficient. The
    two-temperature form needs both:

        q = h_f (T_aw - T_w),   T_aw = T_g - eta_f (T_g - T_c)

    and an adiabatic wall has `q = 0` by construction, so it determines
    `T_aw` and says nothing about `h_f`. Offering `eta_f` with `h_f` left
    at its smooth-duct value therefore over-predicts, and that is what
    this test pins:

        eta_f admissible at h_f = h_0 = (eta_meas (h_i + h_0) - h_i)/h_0

        plate B   0 of 66 points admit any; the largest is -0.147
        plate C   13 of 72 admit some, all below G = 0.45; largest +0.105

    against 0.27 to 0.58 from Baldauf superposed over the ten rows. Admit
    the film and the same measurement demands a gas-side augmentation of
    F = 1.6 to 4.3 to go with it -- which is consistent, and is why the
    shipped closure lumps both into one supplied `h_gas` rather than
    modelling half a pair.
    """
    h_g = er.gas_side()["h_g"]
    admissible = {}
    for tag, series in plates.items():
        admissible[tag] = [
            (p.x, (p.y * (er.internal_h(series, p.x) + h_g) - er.internal_h(series, p.x)) / h_g)
            for p in load_points(series)
        ]

    assert all(v < 0.0 for _, v in admissible["B"]), (
        "plate B now admits a film term at the smooth-duct h_g somewhere; "
        "it did not at any of 66 points"
    )
    assert max(v for _, v in admissible["B"]) < -0.10

    # Plate C does admit a little, at the low-G end only. Recorded rather
    # than rounded away, because it is the honest limit of the claim.
    positive_C = [(g, v) for g, v in admissible["C"] if v > 0.0]
    assert positive_C, "plate C's low-G points used to admit a small film"
    assert max(g for g, _ in positive_C) < 0.45, (
        "plate C now admits film at high G too, not just the low-G end"
    )
    worst_case = max(v for _, v in admissible["C"])
    assert worst_case == pytest.approx(0.105, abs=0.03)

    # What the shipped film model supplies at the same conditions, through
    # the real library. Compared PER PLATE against that plate's own
    # admissible maximum -- comparing plate B's prediction against plate
    # C's headroom would be meaningless.
    gas = er.gas_side()
    rho_c = er.P_AMBIENT / (er.R_AIR * er.T_COOLANT)
    for tag, series in plates.items():
        D = series.geometry["D_mm"] * 1e-3
        X = series.geometry["X_mm"] * 1e-3
        headroom = max(v for _, v in admissible[tag])
        for G in (0.2, 0.6, 1.0):
            M = (G * X * X / (math.pi * D * D / 4.0)) / (gas["rho_g"] * gas["U_g"])
            rows = [
                cb.film_effectiveness_baldauf_2002(
                    n * X / D, M, rho_c / gas["rho_g"], 90.0, X / D, 0.05
                )
                for n in range(1, 11)
            ]
            predicted = cb.film_superposition_sellers(rows)
            assert predicted > 0.2, f"plate {tag} at G={G}: {predicted:.3f}"
            if headroom <= 0.0:
                # Plate B: at the smooth-duct h_g the prediction and the
                # data disagree in SIGN, which is a stronger statement
                # than any ratio.
                continue
            assert predicted > 3.0 * headroom, (
                f"plate {tag} at G={G}: Baldauf+Sellers predicts "
                f"eta_f = {predicted:.3f} against the {headroom:.3f} the data "
                "admits at the smooth-duct h_g -- no longer far enough above"
            )


def test_admitting_the_film_demands_a_gas_side_augmentation(plates) -> None:
    """The other half of the pair, and the number a caller actually needs.

    Inverting the two-temperature form for the augmentation instead of
    refusing the film:

        F = h_f/h_0 = h_i (1 - eta) / ((eta - eta_f) h_0)

    gives F = 1.6 to 4.3 across the two plates with Baldauf's own film.
    Those are high, but the absolute level is confounded: the plate sits
    only x/D_h = 1.5 into the duct, deep in the entry region where
    Dittus-Boelter's fully developed value understates h_0. A standard
    entry correction raises h_0 by 1.75x and brings plate C to 0.9-1.6 --
    ordinary film-cooling augmentation.

    THE CONFOUND CANCELS IN THE B/C RATIO, which is why the finding is
    stated that way: B needs about 1.7x C's augmentation either way, and
    admitting the film barely moves it (1.86 to 1.67 at G = 1.0). So the
    split is NOT explained by the film being absent.
    """
    h_0 = er.gas_side()["h_g"]
    got: dict[str, dict[float, float]] = {}
    for tag, series in plates.items():
        got[tag] = {}
        for G in (0.6, 1.0, 1.4):
            p = min(load_points(series), key=lambda q: abs(q.x - G))
            eta_f = er.film_effectiveness(series, G)
            got[tag][G] = er.implied_gas_side_h(series, p.x, p.y, eta_f) / h_0

    for tag in ("B", "C"):
        for G, f in got[tag].items():
            assert 1.5 < f < 4.5, f"plate {tag} at G={G}: F = {f:.2f}"

    # The ratio is the robust part -- it needs no assumption about h_0 at
    # all, since h_0 divides out.
    for G in (0.6, 1.0, 1.4):
        ratio = got["B"][G] / got["C"][G]
        assert 1.5 < ratio < 2.0, f"G={G}: B/C augmentation ratio {ratio:.2f}"

    # And the film does not collapse the plates: without it the ratio is
    # of the same size, so the split survives either modelling choice.
    for G in (0.6, 1.0, 1.4):
        p_b = min(load_points(plates["B"]), key=lambda q: abs(q.x - G))
        p_c = min(load_points(plates["C"]), key=lambda q: abs(q.x - G))
        bare = er.implied_gas_side_h(plates["B"], p_b.x, p_b.y) / er.implied_gas_side_h(
            plates["C"], p_c.x, p_c.y
        )
        assert abs(bare - got["B"][G] / got["C"][G]) < 0.35, (
            f"G={G}: admitting the film now changes the B/C split materially "
            f"({bare:.2f} against {got['B'][G] / got['C'][G]:.2f})"
        )


def test_the_coolant_heat_up_is_already_inside_the_measured_h(plates) -> None:
    """Why no heat-up term, on evidence rather than on preference.

    The coolant gains 12 K at G = 1.4 and 60 K at G = 0.2 passing through
    the wall, so an explicit heat-up looks obviously missing. It is not:
    86-GT-225's h is fitted from the plate's transient cooling rate
    against the SUPPLY temperature, so it is inside the correlation.

    The test is that adding it does not FLATTEN the error against G -- a
    genuinely missing G-dependent term would. It tilts it instead, and on
    plate C it turns a small positive bias into a large negative one at
    low G.
    """
    series = plates["C"]
    h_g = er.gas_side()["h_g"]
    cp = 1005.0
    grid = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]

    without, with_heatup, rise = [], [], []
    for G in grid:
        measured = min(load_points(series), key=lambda p: abs(p.x - G))
        h_i = er.internal_h(series, G)
        without.append(h_i / (h_i + h_g) / measured.y - 1.0)
        # epsilon-NTU on the coolant passage: the capacity rate caps what
        # the coolant can absorb, so the effective coefficient saturates.
        h_eff = G * cp * (1.0 - math.exp(-h_i / (G * cp)))
        with_heatup.append(h_eff / (h_eff + h_g) / measured.y - 1.0)
        rise.append(h_i * (1.0 - h_i / (h_i + h_g)) * (750.0 - 295.0) / (G * cp))

    # The heat-up really is large at low G, so this is not a small term
    # being waved away.
    assert rise[0] > 40.0 and rise[-1] < 20.0

    span_without = max(without) - min(without)
    span_with = max(with_heatup) - min(with_heatup)
    assert span_with > span_without, (
        f"adding the heat-up now FLATTENS the error ({span_with:.3f} against "
        f"{span_without:.3f}); the term may belong after all"
    )
    # And it is the low-G end that it breaks.
    assert with_heatup[0] < without[0] - 0.04


def test_the_implied_gas_side_h_needs_no_assumption(plates) -> None:
    """The B/C finding, stated the way that does not rest on the duct.

    Inverting the measurement for the gas-side coefficient uses only
    combaero's internal h and the digitised eta. Plate C lands on the
    smooth-duct value; plate B needs up to 2.5 times it.
    """
    duct = er.gas_side()["h_g"]
    got = {}
    for tag, series in plates.items():
        got[tag] = [
            er.implied_gas_side_h(series, p.x, p.y) / duct
            for G in (0.2, 1.0, 1.4)
            for p in [min(load_points(series), key=lambda q: abs(q.x - G))]
        ]
    assert got["C"] == pytest.approx([0.94, 1.27, 1.31], abs=0.06)
    assert got["B"] == pytest.approx([1.44, 2.35, 2.56], abs=0.10)
    # The inversion must actually invert: feeding it back reproduces eta.
    s = plates["B"]
    p = load_points(s)[10]
    assert er.predict(s, p.x, er.implied_gas_side_h(s, p.x, p.y)) == pytest.approx(p.y)


def test_no_jet_parameter_explains_the_split(plates) -> None:
    """Why no threshold is offered, measured rather than asserted.

    The enhancement collapses on none of the velocity ratio, the blowing
    ratio or the momentum flux ratio: at equal values of any of them the
    two plates still differ by well over 25%. Two plates cannot separate a
    jet parameter from a geometry one, so `jet_regime` reports and fits
    nothing.
    """
    duct = er.gas_side()["h_g"]
    samples = {}
    for tag, series in plates.items():
        samples[tag] = [
            (er.jet_regime(series, p.x), er.implied_gas_side_h(series, p.x, p.y) / duct)
            for p in load_points(series)
        ]

    for key in ("velocity_ratio", "blowing_ratio", "momentum_flux_ratio"):
        for target in (1.0, 1.5, 2.0):
            near = {}
            for tag, rows in samples.items():
                hit = min(rows, key=lambda r: abs(r[0][key] - target))
                if abs(hit[0][key] - target) < 0.3:
                    near[tag] = hit[1]
            if len(near) == 2:
                ratio = near["B"] / near["C"]
                assert ratio > 1.25, (
                    f"{key} ~ {target}: B/C is {ratio:.2f}; the plates now "
                    f"nearly collapse on {key}, so a threshold may be "
                    "supportable after all -- revisit the header"
                )


def test_fig10_reproduces_the_papers_own_table_3() -> None:
    """The third channel, and it validates the curve LABELS not just the
    coordinates.

    Page 7 tabulates the ratios BETWEEN these curves at G = 0.2, 0.5 and
    1.0 under "Relative Cooling Effectiveness". The table does not say
    what the ratio is; reading it as `(eta_num - eta_den)/eta_den`
    reproduces all nine printed cells to under 2.1 percentage points, and
    0.57 on average. Nine cells agreeing is not chance, so the reading is
    confirmed as well as the digitisation.
    """
    import bisect

    def curve(name):
        s = next(x for x in load_dataset() if x.label == f"andrews1988/fig10_eta_{name}")
        return sorted((p.x, p.y) for p in load_points(s))

    def at(pts, x):
        xs = [p[0] for p in pts]
        if not xs[0] <= x <= xs[-1]:
            return None
        i = max(1, bisect.bisect_left(xs, x))
        (x0, y0), (x1, y1) = pts[i - 1], pts[i]
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

    data = {
        k: curve(v)
        for k, v in {
            "B": "effusion_B",
            "C": "effusion_C",
            "A/B": "imp_effusion_AB",
            "A/C": "imp_effusion_AC",
            "A": "impingement_A",
        }.items()
    }

    TABLE_3 = [
        (0.2, "A/B", "A", None),
        (0.2, "A/C", "A", None),
        (0.2, "A/B", "B", 13),
        (0.2, "A/C", "C", 21),
        (0.5, "A/B", "A", 6),
        (0.5, "A/C", "A", 17),
        (0.5, "A/B", "B", 24),
        (0.5, "A/C", "C", 20),
        (1.0, "A/B", "A", 5),
        (1.0, "A/C", "A", 13),
        (1.0, "A/B", "B", 26),
        (1.0, "A/C", "C", 19),
    ]
    errs: list[float] = []
    missing: list[tuple] = []
    for G, num, den, printed in TABLE_3:
        a, b = at(data[num], G), at(data[den], G)
        if printed is None:
            # The table prints a dash, and the curve must be absent --
            # impingement A alone starts at G = 0.487. That the two agree
            # is itself a confirmation that this curve is the one the
            # table's A_imp column means.
            assert a is None or b is None, (
                f"Table 3 dashes {num}/{den} at G={G} but both curves reach it"
            )
            continue
        if a is None or b is None:
            # One printed cell is not evaluable: the A/B curve's first
            # digitised point is at G = 0.203, so G = 0.2 would need
            # extrapolation. Recorded rather than clamped -- 0.003 is
            # small, but inventing a point to fill a validation cell is
            # exactly the shortcut this file exists to refuse.
            missing.append((G, num, den))
            continue
        errs.append((a / b - 1.0) * 100 - printed)

    assert missing == [(0.2, "A/B", "B")], f"unexpectedly unevaluable: {missing}"
    assert len(errs) == 9
    assert max(abs(e) for e in errs) < 2.1, f"worst cell off by {max(abs(e) for e in errs):.1f}"
    assert sum(abs(e) for e in errs) / len(errs) < 0.9


def test_fig10_calibration_and_why_the_skew_is_not_corrected() -> None:
    """The page IS rotated, and the raw picks are committed anyway.

    The frame corners give 0.94% skew read as y-per-x-span and 0.96% read
    as x-per-y-span -- two independent measurements of one real rotation.
    But the digitiser's own affine already absorbed it: applying the
    four-corner inverse makes the y ticks WORSE. Both halves are asserted,
    because only the pair justifies leaving the data alone.
    """
    import csv

    path = (
        Path(__file__).resolve().parents[2]
        / "validation"
        / "cooling"
        / "data"
        / "andrews1988"
        / "fig10_calibration.csv"
    )
    rows = [
        (float(r[0]), float(r[1]))
        for r in csv.reader(path.open())
        if not r[0].strip().startswith("x")
    ]
    assert len(rows) == 16
    x_ticks, top_right, top_left, y_ticks = rows[:9], rows[9], rows[10], rows[11:]
    bottom_left, bottom_right = x_ticks[0], x_ticks[-1]

    # Raw: the ticks land on their labels.
    worst_x = max(
        abs(p[0] - t)
        for p, t in zip(x_ticks, [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6], strict=True)
    )
    worst_y = max(abs(p[1] - t) for p, t in zip(y_ticks, [0.9, 0.8, 0.7, 0.6, 0.5], strict=True))
    assert worst_x / 1.6 < 0.006, f"x ticks off by {worst_x:.4f}"
    assert worst_y / 0.6 < 0.005, f"y ticks off by {worst_y:.5f}"

    # The rotation is real, and two independent reads of it agree.
    skew_y = (bottom_right[1] - bottom_left[1]) / 0.6
    skew_x = (bottom_left[0] - top_left[0]) / 1.6
    assert skew_y == pytest.approx(skew_x, rel=0.10), (
        f"the two skew reads disagree ({skew_y:.4f} vs {skew_x:.4f}); "
        "the frame may not be a rigid rotation"
    )
    assert 0.005 < skew_y < 0.02

    # Correcting it makes the y ticks WORSE, which is why it is not applied.
    def unskew_y(x, y):
        yb = bottom_left[1] + (bottom_right[1] - bottom_left[1]) * (
            (x - bottom_left[0]) / (bottom_right[0] - bottom_left[0])
        )
        yt = top_left[1] + (top_right[1] - top_left[1]) * (
            (x - top_left[0]) / (top_right[0] - top_left[0])
        )
        return 0.4 + 0.6 * (y - yb) / (yt - yb)

    worst_corrected = max(
        abs(unskew_y(x, y) - t)
        for (x, y), t in zip(y_ticks, [0.9, 0.8, 0.7, 0.6, 0.5], strict=True)
    )
    assert worst_corrected > worst_y, (
        f"the skew correction now IMPROVES the y ticks ({worst_corrected:.5f} "
        f"against {worst_y:.5f}); apply it and re-derive the data"
    )

    # And the committed DATA must still be the raw picks, not a corrected
    # rebase. Found by falsification: applying the correction to the two
    # scored curves is caught by the Table 3 check only because the other
    # six stay raw -- a UNIFORM correction would slip through, since it
    # moves eta by at most 0.005 and the acceptance bands are wider than
    # that. So two exact coordinates are pinned as the literal record of
    # what was digitised.
    for name, first, last in (
        ("effusion_B", (0.15308, 0.51094), (1.55143, 0.69329)),
        ("effusion_C", (0.09499, 0.51304), (1.57724, 0.79076)),
    ):
        series = next(x for x in load_dataset() if x.label == f"andrews1988/fig10_eta_{name}")
        pts = sorted((q.x, q.y) for q in load_points(series))
        assert pts[0] == pytest.approx(first, abs=5e-5), f"{name} first point"
        assert pts[-1] == pytest.approx(last, abs=5e-5), f"{name} last point"


def test_the_runner_does_not_reimplement_the_correlation(plates) -> None:
    """Predictions must move when the real correlation moves."""
    c = plates["C"]
    thicker = dataclasses.replace(c, geometry={**c.geometry, "thickness_mm": 20.0})
    assert er.predict(thicker, 0.5) != er.predict(c, 0.5)
    assert er.internal_h(c, 0.5) == pytest.approx(
        cb.effusion_internal_nusselt(
            4.0
            * 0.5
            * (15.24e-3) ** 2
            / (
                math.pi
                * 3.27e-3
                * cb.viscosity(er.T_COOLANT, er.P_AMBIENT, cb.standard_dry_air_composition())
            ),
            cb.prandtl(er.T_COOLANT, er.P_AMBIENT, cb.standard_dry_air_composition()),
            15.24 / 6.3,
            6.3 / 3.27,
        )
        * cb.thermal_conductivity(er.T_COOLANT, er.P_AMBIENT, cb.standard_dry_air_composition())
        / 3.27e-3
        * (math.pi * 3.27e-3 * 6.3e-3)
        / ((15.24e-3) ** 2 - math.pi * (3.27e-3) ** 2 / 4.0)
    )


def test_the_gas_side_assumption_is_reported_not_minimised(plates) -> None:
    """Its cost must be visible, and on plate B it is the finding itself.

    Doubling the gas-side h takes plate B from +23.7% to about +1.6% -- so
    a runner that "calibrated" h_g would report plate B as a good fit and
    erase the jet-stirring result entirely. That is precisely why the
    multiplier is reported and never applied.
    """
    sens = er.gas_side_sensitivity(plates["B"])
    assert set(sens) == {0.75, 1.0, 1.5, 2.0}
    assert sens[1.0] > 0.17
    assert abs(sens[2.0]) < 0.05, (
        "doubling h_g no longer nearly closes plate B; the claim in the docstring needs remeasuring"
    )
    # Monotone decreasing: more gas-side h, lower eta, lower over-prediction.
    assert sens[0.75] > sens[1.0] > sens[1.5] > sens[2.0]


def test_no_film_correlation_of_any_magnitude_closes_the_split(plates) -> None:
    """The result that says where the missing physics is NOT.

    A better lift-off film correlation is the obvious next thing to reach
    for -- Badal et al.'s three-regime method (fully attached, transition,
    fully lifted) is a named candidate. This test says it would not close
    this gap, and says so before the effort is spent.

    Scaling Baldauf's `eta_f` from 0 to 1.25x and re-solving for the
    augmentation each plate then requires, the ABSOLUTE F moves a lot --
    plate B's at G = 0.6 runs 2.07 to 5.99 -- but the B/C RATIO barely
    moves at all:

        scale     0.00   0.25   0.50   0.75   1.00   1.25
        B/C       1.85   1.81   1.77   1.72   1.66   1.60   (at G = 1.0)

    Any film correlation can only place `eta_f` somewhere in that span, so
    the split is bounded away from 1 whatever film model is chosen. The
    residual is on the GAS SIDE, and a film correlation -- which is an
    adiabatic measurement and cannot yield a coefficient -- cannot supply
    it. See the module docstring.
    """
    h_0 = er.gas_side()["h_g"]
    ratios = []
    for scale in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25):
        for G in (0.6, 1.0, 1.4):
            per_plate = {}
            for tag, series in plates.items():
                p = min(load_points(series), key=lambda q: abs(q.x - G))
                eta_f = er.film_effectiveness(series, G) * scale
                h_f = er.implied_gas_side_h(series, p.x, p.y, eta_f)
                assert h_f is not None, f"{tag} scale={scale} G={G}"
                per_plate[tag] = h_f / h_0
            ratios.append(per_plate["B"] / per_plate["C"])

    assert min(ratios) > 1.5, (
        f"the B/C split now reaches {min(ratios):.2f} for some film "
        "magnitude; a film correlation may close it after all"
    )
    assert max(ratios) < 2.0
    # And it is the RATIO that is insensitive, not the problem -- the
    # absolute augmentation really does depend strongly on the film, which
    # is why only the ratio is quoted as a finding.
    spread = max(ratios) - min(ratios)
    assert spread < 0.4, f"the ratio moved by {spread:.2f} across film magnitudes"

    plate_b_absolute = [
        er.implied_gas_side_h(plates["B"], p.x, p.y, er.film_effectiveness(plates["B"], 0.6) * s)
        / h_0
        for s in (0.0, 1.0)
        for p in [min(load_points(plates["B"]), key=lambda q: abs(q.x - 0.6))]
    ]
    assert plate_b_absolute[1] / plate_b_absolute[0] > 1.8, (
        "the absolute augmentation is no longer strongly film-dependent; "
        "the reason for quoting only the ratio would need restating"
    )
