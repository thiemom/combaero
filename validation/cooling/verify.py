"""Check digitised series against a figure card read off the printed page.

The card and the points are two independent channels. The card records
what the AXES and the PRINTED EQUATIONS say -- read from the page, without
looking at where the digitiser put anything. The points record where the
marks are. A digitisation bug shows up as disagreement between them.

This is the same discipline validation/cooling/extractions/README.md
applies to reading equations, moved to reading plots. Neither channel is
authoritative; the disagreement is the signal.

What the card catches, with no second digitisation:

  - a dropped axis multiplier (the figure 4.54 bug: two decades)
  - a mis-calibrated axis origin or span
  - points read off the wrong panel
  - a drawn correlation line that does not reproduce its own printed
    equation, which means the calibration is wrong however plausible the
    numbers look
  - double-picked or missed points

What it cannot catch, and still needs a human read:

  - which symbol belongs to which geometry in the legend
  - whether a point was assigned to the right series
  - a figure whose axes are unlabelled in the first place

Run:  uv run python -m validation.cooling.verify
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import yaml

from validation.cooling.schema import (
    Point,
    SeriesMetadata,
    SourceMetadata,
    load_dataset,
    load_points,
)

# A digitised point may sit slightly outside the outermost tick: the axis
# usually runs a little past it, and the mark has width. Wide enough to
# accept honest reading, far too tight to accept a decade.
SPAN_MARGIN = 0.25

# Two points closer than this in both coordinates, relative to the span,
# are almost certainly one mark picked twice.
DUPLICATE_TOL = 1e-3


@dataclass
class Finding:
    series: str
    check: str
    ok: bool
    detail: str


def _span(ticks: list[float], multiplier: float) -> tuple[float, float]:
    lo, hi = min(ticks) * multiplier, max(ticks) * multiplier
    # Log axes are the norm here; widen geometrically so the margin means
    # the same thing at both ends.
    if lo > 0.0:
        factor = (hi / lo) ** SPAN_MARGIN
        return lo / factor, hi * factor
    pad = (hi - lo) * SPAN_MARGIN
    return lo - pad, hi + pad


def _power_slope(points: list[Point]) -> float:
    """Least-squares exponent of a power law through the points."""
    n = len(points)
    sx = sum(math.log(p.x) for p in points)
    sy = sum(math.log(p.y) for p in points)
    sxx = sum(math.log(p.x) ** 2 for p in points)
    sxy = sum(math.log(p.x) * math.log(p.y) for p in points)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def _evaluate_printed(spec: dict, x: float) -> float:
    """The curve the figure prints, at x."""
    if "constant" in spec:
        return float(spec["constant"])
    if "power_law" in spec:
        pl = spec["power_law"]
        return float(pl["C"]) * x ** float(pl["n"])
    if "polynomial" in spec:
        poly = spec["polynomial"]
        scale = float(poly.get("variable_scale", 1.0))
        u = x / scale
        return sum(float(c) * u**i for i, c in enumerate(poly["coeffs"]))
    raise ValueError(f"unrecognised printed-curve spec: {sorted(spec)}")


def check_series(series: SeriesMetadata) -> list[Finding]:
    card = series.verification
    label = series.label
    if card is None:
        return [
            Finding(label, "card", False, "no verification card; series unchecked")
        ]

    points: list[Point] = load_points(series)
    out: list[Finding] = []

    # 1. Axis span. This is what catches a dropped multiplier.
    #
    # `<axis>_limits` overrides the tick span, and either entry may be null
    # for a bound that cannot be read off the scan -- a log axis whose
    # minor ticks continue past the last LABELLED one is the usual case,
    # and Figure 4.193c is exactly that below G = 20. A null bound is not
    # checked rather than being invented, which keeps the card a record of
    # what the page shows.
    vertical_frame = (
        series.kind == "frame"
        and card.get("frame_orientation", "horizontal") == "vertical"
    )
    for axis, values, ticks_key, mult_key, lim_key in (
        ("x", [p.x for p in points], "x_ticks", "x_multiplier", "x_limits"),
        ("y", [p.y for p in points], "y_ticks", "y_multiplier", "y_limits"),
    ):
        ticks = card.get(ticks_key)
        if not ticks:
            continue
        if vertical_frame and axis == "x":
            # A vertical frame sits AT one abscissa; it has no span to check.
            continue
        mult = float(card.get(mult_key, 1.0))
        limits = card.get(lim_key)
        if limits is not None:
            lo = float(limits[0]) * mult if limits[0] is not None else None
            hi = float(limits[1]) * mult if limits[1] is not None else None
        else:
            lo, hi = _span([float(t) for t in ticks], mult)

        outside = [
            v
            for v in values
            if (lo is not None and v < lo) or (hi is not None and v > hi)
        ]
        bound = (
            f"{'unread' if lo is None else format(lo, 'g')} to "
            f"{'unread' if hi is None else format(hi, 'g')}"
        )
        out.append(
            Finding(
                label,
                f"{axis}-span",
                not outside,
                f"axis bounds {bound}; "
                + (
                    f"{len(outside)} of {len(values)} points outside "
                    f"(e.g. {outside[0]:g})"
                    if outside
                    else f"all {len(values)} points inside"
                ),
            )
        )

    # 2. Printed equation. The strongest check available: the equation
    #    comes from the page text, the points from the plot area, and a
    #    mis-calibrated axis cannot reproduce it by accident.
    printed = card.get("printed_curve")
    if printed is not None:
        tol = float(card.get("curve_tolerance", 0.02))
        errs = [_evaluate_printed(printed, p.x) / p.y - 1.0 for p in points]
        rms = math.sqrt(sum(e * e for e in errs) / len(errs))
        out.append(
            Finding(
                label,
                "printed-curve",
                rms <= tol,
                f"RMS against the equation printed on the figure = "
                f"{rms * 100:.2f}% (tolerance {tol * 100:.1f}%)",
            )
        )

    # 3. Printed exponent, where the figure prints a power law whose
    #    COEFFICIENT depends on geometry the curve was drawn at but whose
    #    EXPONENT does not. A free fit must recover it. This is what
    #    separates two curves on one figure -- the 4.193c pair differ only
    #    as 0.35 against 0.42 -- so it catches a swap, which a span check
    #    cannot. It is deliberately NOT a decade check: rescaling x moves
    #    the coefficient and leaves the exponent alone.
    expected_exp = card.get("printed_exponent")
    if expected_exp is not None:
        tol = float(card.get("exponent_tolerance", 0.03))
        fitted = _power_slope(points)
        out.append(
            Finding(
                label,
                "printed-exponent",
                abs(fitted - float(expected_exp)) <= tol,
                f"figure prints {float(expected_exp):g}; free fit gives "
                f"{fitted:.4f} (tolerance {tol:g})",
            )
        )

    # 4. Point count, against what the reader counted on the page.
    expected = card.get("expected_points")
    if expected is not None:
        out.append(
            Finding(
                label,
                "count",
                len(points) == int(expected),
                f"card says {expected} marks, CSV has {len(points)}",
            )
        )

    # 5. Abscissa span, against what a fitted slope needs.
    #
    # A power-law exponent fitted over a short span is dominated by
    # picking noise, however clean the points look. Figure 4.51's 60 deg
    # crossed series is the case: three marks over 0.30 decades, the rest
    # obscured behind other symbols, giving a slope that inverts the
    # ordering every other class shows. The marks are real; the SLOPE is
    # not a measurement, and saying so here keeps a later reader from
    # treating it as one.
    min_decades = card.get("min_decades_for_slope", 0.5)
    if series.kind != "frame" and card.get("x_axis_type", "log") == "log":
        xs_all = [p.x for p in points]
        decades = math.log10(max(xs_all) / min(xs_all)) if min(xs_all) > 0 else 0.0
        wide = decades >= float(min_decades)
        # A short span is not a defect -- marks hide behind each other and
        # three legible points are then the right answer. What must not
        # happen is a short span passing unnoticed and its slope being read
        # as physics. So a narrow series passes only once the metadata
        # ACKNOWLEDGES it, which puts the limitation where a later reader
        # will find it.
        acknowledged = bool(card.get("slope_unreliable", False))
        detail = f"spans {decades:.2f} decades"
        if wide:
            pass
        elif acknowledged:
            detail += (
                f" (< {float(min_decades):g}), acknowledged: usable as data, "
                "its fitted exponent is not a measurement"
            )
        else:
            detail += (
                f" (< {float(min_decades):g}) and not acknowledged; set "
                "slope_unreliable: true if the remaining marks are obscured"
            )
        out.append(Finding(label, "slope-span", wide or acknowledged, detail))

    # 6. Panel distortion, measured on a frame line.
    #
    # A plot frame is horizontal BY CONSTRUCTION, so its fitted slope is
    # the panel's distortion and nothing else. This is a better yardstick
    # than a printed equation, which assumes the draftsman drew the
    # equation faithfully -- exactly what is in question when a line and
    # its label disagree. Figure 4.46 is the case: its R line rises 2.7%
    # against a label reading "= 3.2", and the frame of the panel it is
    # drawn in is flat to -0.00036, so the rise is in the drawing rather
    # than the scan.
    if card.get("kind") == "frame" or series.kind == "frame":
        tol = float(card.get("frame_tolerance", 0.002))
        # A horizontal frame is fitted y against x; a VERTICAL one has x
        # constant and y varying, so it must be fitted the other way round
        # or the slope diverges. A vertical frame measures shear, which a
        # horizontal one cannot see, and it also pins the panel's edge --
        # which is how figure 4.51's two panels were shown aligned before
        # its unlabelled upper abscissa was transferred.
        vertical = card.get("frame_orientation", "horizontal") == "vertical"
        pts = [Point(p.y, p.x) for p in points] if vertical else points
        fitted = _power_slope(pts)
        which = "vertical" if vertical else "horizontal"
        out.append(
            Finding(
                label,
                "frame-slope",
                abs(fitted) <= tol,
                f"frame is {which} by construction; measured slope "
                f"{fitted:+.5f} (tolerance {tol:g}) -- this is the panel's "
                f"distortion, use it to judge every other line in the panel",
            )
        )

    # 7. Double-picked marks.
    xs = [p.x for p in points]
    ys = [p.y for p in points]
    x_span = max(xs) - min(xs) or 1.0
    y_span = max(ys) - min(ys) or 1.0
    dupes = [
        (a, b)
        for i, a in enumerate(points)
        for b in points[i + 1 :]
        if abs(a.x - b.x) / x_span < DUPLICATE_TOL
        and abs(a.y - b.y) / y_span < DUPLICATE_TOL
    ]
    out.append(
        Finding(
            label,
            "distinct",
            not dupes,
            f"{len(dupes)} coincident pairs" if dupes else "no coincident marks",
        )
    )

    # 8. Declared monotonicity, where the physics or the figure demands it.
    trend = card.get("monotonic")
    if trend in ("increasing", "decreasing"):
        ordered = sorted(points, key=lambda p: p.x)
        if trend == "increasing":
            bad = sum(
                1
                for a, b in zip(ordered, ordered[1:], strict=False)
                if b.y < a.y
            )
        else:
            bad = sum(
                1
                for a, b in zip(ordered, ordered[1:], strict=False)
                if b.y > a.y
            )
        out.append(
            Finding(
                label,
                "monotonic",
                bad == 0,
                f"{bad} steps against the declared {trend} trend",
            )
        )

    return out


def check_all() -> list[Finding]:
    out: list[Finding] = []
    for series in load_dataset():
        out.extend(check_series(series))
        if series.class_confidence == "disputed":
            # Reported every run, deliberately. A disputed label that lives
            # only in a metadata comment is a fact with a half-life; one
            # that prints on every verification is not.
            out.append(
                Finding(
                    series.label,
                    "class-label",
                    True,
                    "DISPUTED -- pooled use only, never per-class; "
                    "see this series' cross_check for the evidence",
                )
            )
    return out


def render(findings: list[Finding]) -> str:
    lines = []
    current = None
    for f in findings:
        if f.series != current:
            lines.append("")
            lines.append(f.series)
            current = f.series
        mark = "ok  " if f.ok else "FAIL"
        lines.append(f"  [{mark}] {f.check:<16} {f.detail}")
    failed = [f for f in findings if not f.ok]
    lines.append("")
    lines.append(f"{len(findings) - len(failed)} passed, {len(failed)} failed")
    return "\n".join(lines).lstrip("\n")


def check_panel_pairing(
    directory: Path, pattern_a: str, pattern_b: str, tol_pct: float = 3.0
) -> list[Finding]:
    """Pair one symbol class across two panels of the same figure.

    Where a figure stacks two panels over one abscissa, each experimental
    run is plotted ONCE IN EACH PANEL at the same abscissa. So a class's
    two files must hold the same number of points, pairing to within a
    per cent or two.

    This assumes nothing -- no card, no model, no legend reading -- which
    makes it the sharpest check available on a scatter pass. It catches
    what a span check cannot see: a missed mark, a mark picked twice, and
    a series filed under the wrong class.
    """
    suffix_a = pattern_a.rsplit("*", 1)[-1]
    suffix_b = pattern_b.rsplit("*", 1)[-1]

    out: list[Finding] = []
    for path_a in sorted(directory.glob(pattern_a)):
        stem = path_a.name
        label = stem[: -len(suffix_a)] if suffix_a else stem
        path_b = directory / (label + suffix_b)
        if not path_b.exists():
            out.append(Finding(label, "pairing", False, f"no counterpart {path_b.name}"))
            continue
        a = sorted(load_points_from(path_a), key=lambda p: p.x)
        b = sorted(load_points_from(path_b), key=lambda p: p.x)
        if len(a) != len(b):
            out.append(
                Finding(label, "pair-count", False,
                        f"{len(a)} points here against {len(b)} in "
                        f"{path_b.name}; every run appears in both panels")
            )
        matched = 0
        worst = 0.0
        for pa in a:
            nearest = min(b, key=lambda pb: abs(math.log(pb.x / pa.x)))
            d = (nearest.x / pa.x - 1.0) * 100.0
            if abs(d) < tol_pct:
                matched += 1
            worst = max(worst, abs(d))
        out.append(
            Finding(label, "pair-abscissa", matched == len(a),
                    f"{matched}/{len(a)} paired within {tol_pct:g}%"
                    f" (worst {worst:.1f}%)")
        )
    return out


def load_points_from(path: Path) -> list[Point]:
    """Read a bare x,y CSV with no metadata around it."""
    import csv as _csv

    pts: list[Point] = []
    with open(path, newline="") as fh:
        for row in _csv.reader(fh):
            if not row or not row[0].strip() or row[0].strip() == "x":
                continue
            pts.append(Point(float(row[0]), float(row[1])))
    return pts


def check_candidate(csv_path: Path, card_path: Path) -> list[Finding]:
    """Check a freshly digitised CSV before it joins the dataset.

    The workflow the cards are for is digitise -> check -> commit, not
    commit -> discover. This runs the same checks against a loose file and
    a loose card, so a bad calibration is caught while the figure is still
    open rather than after it is in the tree.
    """
    card = yaml.safe_load(card_path.read_text())
    series = SeriesMetadata(
        path=csv_path,
        source=SourceMetadata(name=csv_path.parent.name or "candidate",
                              citation="(not yet filed)", secondary=True),
        after=None,
        page=None,
        item=card.get("item", "(candidate)"),
        series=card.get("series", csv_path.stem),
        geometry=None,
        alpha_deg=None,
        x_axis=card.get("x_axis", "e_plus"),
        x_scale=float(card.get("x_scale", 1.0)),
        y_axis=card.get("y_axis", "G"),
        kind=card.get("kind", "correlation"),
        extraction="figure-digitised",
        confidence="band",
        uncertainty=card.get("uncertainty"),
        cross_check="(candidate; not yet recorded)",
        scores=None,
        verification=card.get("verification", card),
    )
    return check_series(series)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check digitised series against their figure cards."
    )
    parser.add_argument(
        "--candidate",
        type=Path,
        help="a freshly digitised CSV, not yet in the dataset",
    )
    parser.add_argument(
        "--candidate-dir",
        type=Path,
        help="a directory of freshly digitised CSVs, all from one panel",
    )
    parser.add_argument(
        "--glob",
        default="*.csv",
        help="which files in --candidate-dir to check (default: *.csv)",
    )
    parser.add_argument(
        "--card",
        type=Path,
        help="the figure card for the candidate(s), as a YAML file",
    )
    parser.add_argument(
        "--pair",
        nargs=2,
        metavar=("GLOB_A", "GLOB_B"),
        help="pair one class across two panels, e.g. '*_R.csv' '*_G.csv'",
    )
    args = parser.parse_args()

    if args.pair:
        if not args.candidate_dir:
            parser.error("--pair needs --candidate-dir")
        print(render(check_panel_pairing(args.candidate_dir, *args.pair)))
        return

    if args.candidate or args.candidate_dir:
        if not args.card:
            parser.error("a candidate needs --card")
        if args.candidate:
            print(render(check_candidate(args.candidate, args.card)))
            return
        paths = sorted(args.candidate_dir.glob(args.glob))
        if not paths:
            parser.error(f"no files matching {args.glob} in {args.candidate_dir}")
        findings = []
        for path in paths:
            findings.extend(check_candidate(path, args.card))
        print(render(findings))
        return
    print(render(check_all()))


if __name__ == "__main__":
    main()
