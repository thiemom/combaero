"""Score multi-row film/effusion effectiveness against adiabatic plate data.

The chain under test is Baldauf (2002) for one row, superposed over rows by
Sellers (Gao Eq. 1) or Gao's corrected form (Eq. 7):

    eta_row_i(x) = film_effectiveness_baldauf_2002((x - x_i) * s_x/d, ...)
    eta(x)       = 1 - prod_i (1 - eta_row_i(x))        rows i upstream of x

Both are called through the real implementation, never reimplemented here --
the same discipline as `jet_array_runner`. The only arithmetic this module
owns is bookkeeping: which rows are upstream, how far, and how to average.

WHY ROW MEANS. The measured curve carries a sawtooth, one tooth per hole row,
because eta jumps at each injection and decays between. Superposing
row-resolved Baldauf curves reproduces that structure in principle, but its
phase and amplitude depend on the within-pitch detail of a single row's
footprint, which a laterally averaged correlation does not claim to give. So
both sides are reduced to per-row means over each pitch interval before
comparison, as `andrei2014`'s metadata header specifies. Scoring the raw
sawtooth would report a phase disagreement as a magnitude error.

WHAT A MISS MEANS HERE. `andrei2014` sits outside Baldauf's published
envelope on every curve -- s/D 7.37 against a maximum of 5 throughout, plus
P = 1.0 and M = 3 on some -- so this is a cross-source ACCURACY measurement
on an extrapolated model, never a fidelity check. Records carry the
extrapolation flag from the correlation itself rather than a judgement made
here. See `tests/test_film_effectiveness.cpp`'s
`AndreiTwentyFourteenRigIsOutsideTheEnvelope`.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import combaero as cb

from validation.cooling.schema import Point, SeriesMetadata, load_points

# Mainstream turbulence intensity. Andrei et al. do not state one, and it is a
# required Baldauf input, so this is an ASSUMPTION and is reported as such
# rather than hidden: 5% is a representative combustor-rig value and sits
# inside Baldauf's 0.35-7.5% envelope, so it contributes no extrapolation of
# its own. Sensitivity is reported by `turbulence_sensitivity` below instead of
# being argued about.
ASSUMED_TU = 0.05

# Superposition with alpha = 1 IS Sellers (Gao Eq. 7 collapses onto Eq. 1),
# so this is the uncorrected baseline, not a tuned one. Gao never publishes
# fitted (a, b) for alpha, so nothing here may invent them.
UNCORRECTED_ALPHA = 1.0


@dataclass(frozen=True)
class Record:
    """One row-mean comparison."""

    series: SeriesMetadata
    x: float  # pitch-interval centre, x/s_x
    measured: float  # mean measured eta over that interval
    predicted: float | None
    extrapolated: bool
    reason: str | None = None
    n_samples: int = 0  # digitised samples behind `measured`

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


def owns(series: SeriesMetadata) -> bool:
    """Whether this runner should score ``series``.

    Explicit predicate, because `run_series` returns reason-carrying records
    for series it cannot score, which makes truthiness useless for dispatch.
    """
    return series.x_axis == "x_over_sx" and series.y_axis == "eta_adiabatic"


def _blowing_and_density(series: SeriesMetadata) -> tuple[float, float] | None:
    """Read BR and DR out of the panel/series labels.

    `andrei2014` encodes them in the filename and repeats them in `series`;
    they are conditions of the run, not geometry, so they do not live in
    `geometry`.
    """
    stem = series.label.split("/")[-1]
    br = dr = None
    for token in stem.split("_"):
        if token.startswith("BR"):
            br = float(token[2:].replace("p", "."))
        elif token.startswith("DR"):
            dr = float(token[2:].replace("p", "."))
    if br is None or dr is None:
        return None
    return br, dr


def row_positions(series: SeriesMetadata) -> list[float]:
    """Row abscissae in x/s_x: phase, then unit pitch spacing.

    `row1_x_over_sx` is an extraction fact recovered from the data's own
    injection points and cross-checked three ways -- see the metadata header.
    Unit spacing is what `x/s_x` means.
    """
    g = series.geometry
    n = int(g["rows"])
    x1 = float(g["row1_x_over_sx"])
    return [x1 + i for i in range(n)]


def predict(
    series: SeriesMetadata, x: float, Tu: float = ASSUMED_TU
) -> tuple[float | None, bool]:
    """Superposed eta at x/s_x, and whether any row's call extrapolated."""
    cond = _blowing_and_density(series)
    if cond is None:
        return None, False
    M, P = cond

    g = series.geometry
    sx_over_d = float(g["sx_over_d"])
    sz_over_d = float(g["sz_over_d"])
    alpha_deg = float(g["angle_deg"])

    per_row: list[float] = []
    for x_row in row_positions(series):
        if x_row >= x:
            break  # downstream rows cannot cool an upstream point
        x_over_D = (x - x_row) * sx_over_d
        per_row.append(
            cb.film_effectiveness_baldauf_2002(
                x_over_D, M, P, alpha_deg, sz_over_d, Tu
            )
        )
    if not per_row:
        return None, False

    # alpha = 1 collapses Eq. (7) onto Sellers, so this is the uncorrected
    # baseline; the corrected form needs Gao's unpublished (a, b).
    eta = cb.film_superposition_sellers(per_row)

    # The correlation owns the extrapolation verdict, and s/D, M, P, alpha
    # and Tu are all fixed across a series -- so it is asked ONCE per
    # (series, Tu) and cached, not once per sample. That is not only
    # wasteful: each probe swaps the global warning handler and restores it,
    # and `get_warning_handler` hands back a fresh wrapper around the current
    # one, so restoring nests. Doing it ~6000 times per scorecard run nested
    # deeply enough to overflow the stack and segfault pytest.
    extrapolated = _is_extrapolated(M, P, alpha_deg, sz_over_d, Tu)
    return eta, extrapolated


@lru_cache(maxsize=None)
def _is_extrapolated(
    M: float, P: float, alpha_deg: float, s_over_D: float, Tu: float
) -> bool:
    """Ask the correlation, rather than restating its envelope here.

    Python gets the envelope verdict as warnings -- the `CorrelationStatus`
    out-parameter is C++-only -- so this captures them through the handler.
    Restating the bounds in this file would let them drift from the header
    that owns them.

    CACHED, and that is a correctness requirement rather than an
    optimisation. `get_warning_handler` returns a fresh wrapper around the
    current handler, so save-and-restore NESTS one layer per call; several
    thousand calls nest deeply enough that the next `warn()` overflows the
    stack. The inputs are constant across a series, so the cache also makes
    the call count O(series) instead of O(samples).
    """
    seen: list[str] = []
    previous = cb.get_warning_handler()
    try:
        cb.set_warning_handler(seen.append)
        cb.film_effectiveness_baldauf_2002(10.0, M, P, alpha_deg, s_over_D, Tu)
    finally:
        cb.set_warning_handler(previous)
    return any("outside validated range" in m for m in seen)


def row_intervals(
    series: SeriesMetadata, points: list[Point]
) -> list[tuple[float, float, list[float], list[float]]]:
    """Group samples by pitch interval: (centre, end, abscissae, eta values).

    One interval per row, from that row to the next. The abscissae are handed
    back rather than just their mean because the prediction must be averaged
    over THE SAME points -- see `run_series`.

    A row whose interval the digitisation does not fully cover is still
    reported, with its sample count, so a thin interval is visible rather
    than silently carrying the weight of a full one.
    """
    rows = row_positions(series)
    out: list[tuple[float, float, list[float], list[float]]] = []
    for i, x_row in enumerate(rows):
        x_end = rows[i + 1] if i + 1 < len(rows) else x_row + 1.0
        xs = [p.x for p in points if x_row <= p.x < x_end]
        ys = [p.y for p in points if x_row <= p.x < x_end]
        if not xs:
            continue
        out.append((0.5 * (x_row + x_end), x_end, xs, ys))
    return out


def row_means(
    series: SeriesMetadata, points: list[Point]
) -> list[tuple[float, float, int]]:
    """(interval centre, mean measured eta, sample count) per row."""
    return [
        (x_mid, sum(ys) / len(ys), len(ys))
        for x_mid, _end, _xs, ys in row_intervals(series, points)
    ]


def run_series(series: SeriesMetadata) -> list[Record]:
    """Evaluate one series as per-row means."""
    points = load_points(series)

    if series.scores is None:
        return [Record(series, p.x, p.y, None, False, "not scored by any set") for p in points]
    if not owns(series):
        return [
            Record(
                series,
                p.x,
                p.y,
                None,
                False,
                f"x_axis={series.x_axis}, y_axis={series.y_axis} not implemented",
            )
            for p in points
        ]

    needed = ("rows", "row1_x_over_sx", "sx_over_d", "sz_over_d", "angle_deg")
    missing = [k for k in needed if k not in (series.geometry or {})]
    if missing:
        return [
            Record(series, p.x, p.y, None, False, f"geometry lacks {', '.join(missing)}")
            for p in points
        ]
    if _blowing_and_density(series) is None:
        return [
            Record(series, p.x, p.y, None, False, "cannot read BR/DR from the label")
            for p in points
        ]

    records: list[Record] = []
    for x_mid, _end, xs, ys in row_intervals(series, points):
        measured = sum(ys) / len(ys)

        # Average the PREDICTION over the same abscissae as the measurement.
        # Comparing an interval mean against a single point value at the
        # interval centre is not comparing like with like: eta rises steeply
        # after each injection, so over row 1 -- where it climbs from zero --
        # the centre value overstates the interval mean badly. Reusing the
        # measured abscissae also makes the two means share a sampling
        # distribution, so an uneven digitisation cannot bias one side only.
        preds = []
        extrapolated = False
        for x in xs:
            value, flag = predict(series, x)
            if value is None:
                continue
            preds.append(value)
            extrapolated = extrapolated or flag

        predicted = sum(preds) / len(preds) if preds else None
        reason = None if predicted is not None else "no row upstream of this interval"
        records.append(
            Record(series, x_mid, measured, predicted, extrapolated, reason, len(ys))
        )
    return records


def turbulence_sensitivity(
    series: SeriesMetadata, levels: tuple[float, ...] = (0.01, 0.05, 0.075)
) -> dict[float, float]:
    """Mean |rel error| against Tu, since Andrei state no turbulence intensity.

    Reported, never minimised: picking the Tu that scores best would be
    fitting an unmeasured input to the metric it is judged by.
    """
    points = load_points(series)
    intervals = row_intervals(series, points)
    out: dict[float, float] = {}
    for Tu in levels:
        errs = []
        for _x_mid, _end, xs, ys in intervals:
            measured = sum(ys) / len(ys)
            preds = [v for v, _ in (predict(series, x, Tu) for x in xs) if v is not None]
            if preds and measured:
                errs.append(abs(sum(preds) / len(preds) / measured - 1.0))
        if errs:
            out[Tu] = sum(errs) / len(errs)
    return out


def run_all(dataset: list[SeriesMetadata]) -> list[Record]:
    out: list[Record] = []
    for series in dataset:
        out.extend(run_series(series))
    return out
