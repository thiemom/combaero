"""Evaluate an implemented correlation set against the digitised series.

The runner never computes G from a set's coefficients itself. It calls
cb.evaluate_rib and bisects the Reynolds number until the chain's own e+
lands on the digitised abscissa, so the whole chain -- friction factor,
e+, then G -- is what gets scored. Reimplementing the formula here would
score the model against a second copy of itself, which is the failure
mode #333 exists to prevent.
"""

from __future__ import annotations

from dataclasses import dataclass

import combaero as cb

from validation.cooling.schema import Point, SeriesMetadata, load_points

# Han states G_bar (the ribbed/smooth average) as 1.2 G. Confirmed as item
# 10 of validation/cooling/extractions/han_ribbed.md, from the label on
# Figure 4.47, and independently by Figure 4.46's printed pair 3.7 and 4.5
# (ratio 1.216).
#
# IT IS NOT UNIVERSAL. Digitising both panels of Figure 4.51 measures
# G_bar/G for seven rib configurations directly: 90 deg gives 1.2193,
# matching the printed 1.2162 to 0.3%, while the angled configurations
# average 1.155 -- 4.5% below it, with a spread of 1.153 to 1.220.
#
# So this applies to 90 deg ribs, where it was established. Applying it to
# an angled rib is wrong by about 4.5%, which is larger than the 2.6%
# scatter of the measurement itself. Anything off 90 deg is refused rather
# than scaled; see _gbar_reason.
G_BAR_OVER_G = 1.2
G_BAR_VALID_ALPHA = 90.0

# Bracket for the Re bisection. Deliberately far wider than any set's
# stated validity: a series may sit outside it, and reporting that as
# extrapolated is the point rather than something to avoid.
RE_LO = 1.0e3
RE_HI = 1.0e8
BISECT_STEPS = 200


@dataclass(frozen=True)
class Record:
    """One digitised point, evaluated."""

    series: SeriesMetadata
    x: float  # e+
    measured: float
    predicted: float | None  # None when the set has no binding
    extrapolated: bool
    re_used: float | None
    reason: str | None = None  # why there is no prediction

    @property
    def rel_error(self) -> float | None:
        if self.predicted is None or self.measured == 0.0:
            return None
        return self.predicted / self.measured - 1.0

    @property
    def within_uncertainty(self) -> bool:
        band = self.series.uncertainty
        err = self.rel_error
        if band is None or err is None:
            return False
        return abs(err) <= band


SETS = {
    "han_1988_orthogonal": cb.han_1988_orthogonal,
    "rallabandi_2009_high_re": cb.rallabandi_2009_high_re,
    "han_park_1988_angled": cb.han_park_1988_angled,
}


def _binding_reason(
    rib_set: "cb.RibCorrelationSet", series: SeriesMetadata
) -> str | None:
    """Why this set cannot answer for this series, or None if it can.

    A set outside its NUMERIC validity is still answering the question --
    that is an extrapolation, and the scorecard reports it as one. A set
    asked about a CONFIGURATION it has no binding for is not: scoring
    han_1988_orthogonal, whose valid_alpha is [90, 90] and whose G carries
    no alpha term, against 45 deg rib data produced a 26.7% bias that reads
    as model error when it is really the wrong correlation entirely. That
    is the shape of misleading metric this harness exists to prevent, so
    the configuration mismatch is refused rather than scored.
    """
    alpha = series.alpha_deg
    rng = rib_set.valid_alpha
    if alpha is not None and rng.hi >= rng.lo:
        if not (rng.lo <= alpha <= rng.hi):
            return (
                f"{alpha:g} deg ribs; {rib_set.name} binds "
                f"{rng.lo:g}-{rng.hi:g} deg only"
            )
    return None


def _probe_geometry(rib_set: "cb.RibCorrelationSet") -> cb.RibGeometry:
    """A geometry at the centre of the set's stated validity.

    Only legitimate when the set's G carries no geometry dependence -- see
    _assert_geometry_free. The probe then fixes only WHICH Reynolds number
    reaches a given e+, not the G reported there.
    """
    geom = cb.RibGeometry()
    geom.e_D = _mid(rib_set.valid_eD, 0.0625)
    geom.p_e = _mid(rib_set.valid_pe, 10.0)
    geom.W_H = _mid(rib_set.valid_WH, 1.0)
    geom.alpha_deg = 90.0
    return geom


def _mid(rng: "cb.RibRange", fallback: float) -> float:
    return 0.5 * (rng.lo + rng.hi) if rng.hi > rng.lo else fallback


def _gbar_reason(series: SeriesMetadata) -> str | None:
    """Why a G_bar series cannot be converted, or None if it can.

    G_bar/G is a per-configuration quantity, not a constant: figure 4.51
    ranks nine configurations in both and the lowest differs between the
    panels, which a constant multiplier cannot do. The 1.2 is a 90 degree
    result. Applying it to an angled rib would silently manufacture a
    number, so it is refused instead.
    """
    if series.y_axis != "G_bar":
        return None
    alpha = series.alpha_deg
    if alpha is None or alpha == G_BAR_VALID_ALPHA:
        return None
    return (
        f"G_bar at {alpha:g} deg: the 1.2 ratio is a 90 deg result and "
        "varies by configuration (figure 4.51 ranks G and G_bar "
        "differently), so it is not applied here"
    )


def _binds_geometry(rib_set: "cb.RibCorrelationSet") -> bool:
    """True when the set's G actually depends on the rig geometry."""
    return any(
        t.exponent != 0.0
        for t in (rib_set.G_eD, rib_set.G_pe, rib_set.G_WH, rib_set.G_alpha)
    )


def _assert_geometry_free(rib_set: "cb.RibCorrelationSet") -> None:
    """Refuse to score a geometry-dependent set against geometry-less data.

    Applies only to series with NO recorded rig geometry -- the Han and
    Zhang (1992) figures, whose rig the extracted text never states. For a
    set whose G has zero exponents on e/D, P/e, W/H and alpha this costs
    nothing. For any other set it would mean substituting a default and
    reporting the result as a measurement, so the runner stops instead.

    Figure 4.46 states a geometry per class in its legend, so those series
    are exempt and reach the disputed-label check above instead.
    """
    terms = {
        "e/D": rib_set.G_eD,
        "P/e": rib_set.G_pe,
        "W/H": rib_set.G_WH,
        "alpha": rib_set.G_alpha,
    }
    bound = [name for name, t in terms.items() if t.exponent != 0.0]
    if bound:
        raise ValueError(
            f"{rib_set.name} binds G to {', '.join(bound)}, but the digitised "
            "series carries no rig geometry. Supply the geometry in "
            "metadata.yaml from the primary paper before scoring this set."
        )


def _normalised_R(
    rib_set: "cb.RibCorrelationSet", geom: cb.RibGeometry, raw_R: float
) -> float:
    """Raw R divided back down to what a figure's y-axis actually plots.

    Both R shapes in this project print their correlation as R divided by
    its own normalisers (Han 1988: R/(p/e/10)^0.35; Han and Park 1988:
    R/[(p/e/10)^0.35 (W/H)^m]) -- never raw R. Recomputing the divisor from
    the SAME exposed fields evaluate_rib used keeps this from silently
    drifting out of step with a future change to either set's shape.
    """
    pe_term = rib_set.R_pe
    pe_factor = (
        (geom.p_e / pe_term.reference) ** pe_term.exponent
        if pe_term.reference and pe_term.exponent
        else 1.0
    )
    WH_factor = 1.0
    if rib_set.R_alpha_shape == cb.RAlphaShape.QuadraticAlpha:
        at_90 = abs(geom.alpha_deg - 90.0) < 1e-9
        m = (
            rib_set.R_quad_WH_exponent_at_90
            if at_90
            else rib_set.R_quad_WH_exponent_off_90
        )
        if m:
            W_H = geom.W_H
            if rib_set.R_quad_WH_cap > 0.0:
                W_H = min(W_H, rib_set.R_quad_WH_cap)
            WH_factor = W_H**m
    divisor = pe_factor * WH_factor
    return raw_R / divisor if divisor else raw_R


def _g_at_eplus(
    rib_set: "cb.RibCorrelationSet", geom: cb.RibGeometry, target: float
) -> tuple[float, float, bool, float] | None:
    """G and the normalised R at a target e+, via the real chain.

    The lower panel of Figure 4.46 plots R/(P/e/10)^0.35, which is what
    the set's C_R is defined as, so the normalised value is returned
    alongside G rather than being recomputed by the caller.
    """
    lo, hi = RE_LO, RE_HI
    if cb.evaluate_rib(rib_set, geom, lo).e_plus > target:
        return None
    if cb.evaluate_rib(rib_set, geom, hi).e_plus < target:
        return None
    for _ in range(BISECT_STEPS):
        mid = 0.5 * (lo + hi)
        if cb.evaluate_rib(rib_set, geom, mid).e_plus < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9 * max(1.0, hi):
            break
    re = 0.5 * (lo + hi)
    res = cb.evaluate_rib(rib_set, geom, re)
    pe_term = rib_set.R_pe
    norm = (geom.p_e / pe_term.reference) ** pe_term.exponent if pe_term.reference else 1.0
    r_norm = res.R / norm if norm else res.R
    return res.G, r_norm, res.extrapolated, re


def _r_at_alpha(
    rib_set: "cb.RibCorrelationSet", geom: cb.RibGeometry, alpha_deg: float
) -> float:
    """R at a given rib angle. No Re bisection needed.

    Every set in this project has R carrying no e+ term (see
    rib_correlation.h's header comment: "R CARRIES NO e+ TERM by
    construction"), so R does not depend on Reynolds number at all. Any Re
    gives the same R -- ARBITRARY_RE below is not a choice that matters,
    only a value evaluate_rib needs to run.
    """
    geom.alpha_deg = alpha_deg
    return cb.evaluate_rib(rib_set, geom, ARBITRARY_RE).R


ARBITRARY_RE = 30000.0


def run_series(
    series: SeriesMetadata, points: list[Point] | None = None
) -> list[Record]:
    """Evaluate one series. Unscored series yield records with no prediction.

    `points` overrides the series' own CSV. Its only use is asking the set
    what it predicts at abscissae the digitisation could not pick -- see
    `recovery.py` -- so the recovered bounds are checked against the same
    chain as everything else rather than a second implementation of it.
    """
    if points is None:
        points = load_points(series)

    if series.scores is None:
        return [
            Record(series, p.x, p.y, None, False, None, "not scored by any set")
            for p in points
        ]
    if series.x_axis not in ("e_plus", "alpha_deg"):
        return [
            Record(series, p.x, p.y, None, False, None, f"x axis is {series.x_axis}")
            for p in points
        ]
    if series.x_axis == "alpha_deg" and not series.y_axis.startswith("R"):
        # The alpha-indexed path below only computes R (R, R_normalised,
        # R_normalised_angled -- every R-family name in this dataset starts
        # with "R"), the only quantity independent of e+/Re in every set
        # implemented so far.
        return [
            Record(
                series, p.x, p.y, None, False, None,
                f"alpha-indexed scoring is only implemented for an R "
                f"quantity, not {series.y_axis}",
            )
            for p in points
        ]

    rib_set = SETS[series.scores]()

    gbar_block = _gbar_reason(series)
    if gbar_block is not None:
        return [
            Record(series, p.x, p.y, None, False, None, gbar_block)
            for p in points
        ]

    if series.class_confidence == "disputed" and _binds_geometry(rib_set):
        # A disputed label is usable while nothing depends on it. The
        # moment a set binds geometry, the label IS the input, and a
        # disputed input must not be fed in silently.
        return [
            Record(
                series, p.x, p.y, None, False, None,
                f"class label disputed and {rib_set.name} binds geometry",
            )
            for p in points
        ]

    if not series.geometry:
        # Only series with no recorded rig geometry need this. The figure
        # 4.46 legend states one per class, so those are exempt.
        _assert_geometry_free(rib_set)

    unsupported = _binding_reason(rib_set, series)
    if unsupported is not None:
        return [
            Record(series, p.x, p.y, None, False, None, unsupported) for p in points
        ]

    geom = _probe_geometry(rib_set)
    if series.geometry:
        # The figure 4.46 legend states each class's geometry, so use it.
        # It matters for R, whose plotted ordinate is normalised by P/e.
        geom.e_D = float(series.geometry.get("e_D", geom.e_D))
        geom.p_e = float(series.geometry.get("p_e", geom.p_e))
        geom.W_H = float(series.geometry.get("W_H", geom.W_H))
    # Rib angle, which _probe_geometry defaults to 90. Omitting it scored
    # every angled series as though its ribs were transverse: figure 4.51's
    # 45 and 60 deg classes were evaluated at 90 against
    # han_park_1988_angled, the one set whose whole subject is rib angle.
    # `alpha_deg` is a top-level field, not part of the geometry mapping,
    # which is how it was missed when e_D/p_e/W_H were wired through.
    if series.alpha_deg is not None:
        geom.alpha_deg = float(series.alpha_deg)

    if series.x_axis == "alpha_deg":
        # R vs alpha (figure 4.47's own axis): no e+/Re bisection needed --
        # R does not depend on Re for any set in this project. alpha itself
        # is the point's x-value, not a fixed series field, so extrapolation
        # is read per point from evaluate_rib's own flag rather than the
        # series-level _binding_reason check above (which only applies when
        # alpha is fixed for the whole series).
        records = []
        for pt in points:
            geom.alpha_deg = pt.x
            res = cb.evaluate_rib(rib_set, geom, ARBITRARY_RE)
            predicted = res.R
            if series.y_axis in ("R_normalised", "R_normalised_angled"):
                # The figure plots R divided by its own normalisers, not raw
                # R -- dividing them back out here, from the SAME fields
                # evaluate_rib used, rather than duplicating Eq. 4.17's
                # formula. Mirrors what _g_at_eplus does for han_1988's
                # R_pe division below.
                predicted = _normalised_R(rib_set, geom, res.R)
            records.append(
                Record(series, pt.x, pt.y, predicted, res.extrapolated, ARBITRARY_RE)
            )
        return records

    records: list[Record] = []
    for p in points:
        found = _g_at_eplus(rib_set, geom, p.x)
        if found is None:
            records.append(
                Record(series, p.x, p.y, None, True, None, "e+ unreachable")
            )
            continue
        g, r_norm, extrapolated, re = found
        if series.y_axis == "G_bar":
            predicted = g * G_BAR_OVER_G
        elif series.y_axis == "R_normalised":
            predicted = r_norm
        else:
            predicted = g
        records.append(Record(series, p.x, p.y, predicted, extrapolated, re))
    return records


def run_all(dataset: list[SeriesMetadata]) -> list[Record]:
    out: list[Record] = []
    for series in dataset:
        out.extend(run_series(series))
    return out
