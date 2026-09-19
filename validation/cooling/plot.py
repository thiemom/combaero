"""Redraw a digitised figure from the committed data.

A numeric check tells you a series is self-consistent. It does not tell
you the cloud has the shape the page shows, that a class sits where your
eye says it should, or that two series were swapped in a way the
arithmetic happens to tolerate. Putting the data back on axes and holding
it next to the scan catches that in a second.

This reads the SAME metadata the verifier and scorecard use -- axis type,
tick span and plot limits all come from each series' figure card -- so a
plot cannot drift from what the checks believe.

    uv run python -m validation.cooling.plot --figure 4.46
    uv run python -m validation.cooling.plot --list

Needs the `examples` extra for matplotlib:  uv pip install -e ".[examples]"
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from validation.cooling.schema import SeriesMetadata, load_dataset, load_points

# Distinct markers so ten overlapping symbol classes stay separable -- the
# same problem the source figure has, which is why classes are hard to
# attribute in the first place.
MARKERS = ["o", "s", "^", "v", "D", "<", ">", "P", "X", "*", "h", "p"]

PANEL_ORDER = {"upper": 0, "single": 1, "lower": 2}

# Default output. Gitignored: a redrawn figure is a derived artefact, and
# regenerating it is one command. It carries its source citation on the
# image itself, so committing one later is a policy choice rather than a
# rework -- a reproduction from our own measurements is ours to publish
# provided the source stays named on it.
PLOT_DIR = Path(__file__).parent / "plots"

# What the figures actually print on their axes, rather than the terse
# keys the metadata uses.
AXIS_NAMES = {
    "e_plus": "e+ = (e/D) Re (f/2)^0.5",
    "G": "G",
    "G_bar": "G_bar",
    "R_normalised": "R / (P/e/10)^0.35",
    "f_ratio": "f / f_0",
    "Nu_ratio": "Nu_r / Nu_0",
}


def _curve(spec: dict, xs: list[float]) -> list[float]:
    """Evaluate a printed_curve spec across xs."""
    if "constant" in spec:
        return [float(spec["constant"])] * len(xs)
    if "power_law" in spec:
        pl = spec["power_law"]
        return [float(pl["C"]) * x ** float(pl["n"]) for x in xs]
    if "polynomial" in spec:
        poly = spec["polynomial"]
        scale = float(poly.get("variable_scale", 1.0))
        return [
            sum(float(c) * (x / scale) ** i for i, c in enumerate(poly["coeffs"]))
            for x in xs
        ]
    raise ValueError(f"unrecognised printed_curve: {sorted(spec)}")


def _label(s: SeriesMetadata) -> str:
    geom = s.geometry or {}
    if geom:
        base = (
            f"e/D {geom.get('e_D')}  P/e {geom.get('p_e')}  W/H {geom.get('W_H')}"
        )
    else:
        base = s.series[:46]
    if s.class_confidence == "disputed":
        base += "  [DISPUTED]"
    return base


def plot_figure(
    figure: str, out: Path, dataset=None, show_frames: bool = False
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dataset = dataset or load_dataset()
    chosen = [s for s in dataset if s.figure == figure]
    if not chosen:
        raise SystemExit(f"no series for figure {figure}")

    panels: dict[str, list[SeriesMetadata]] = defaultdict(list)
    for s in chosen:
        panels[s.panel or "single"].append(s)
    names = sorted(panels, key=lambda p: PANEL_ORDER.get(p, 1))

    fig, axes = plt.subplots(
        len(names), 1, figsize=(9.5, 5.2 * len(names)), squeeze=False
    )

    for ax, panel in zip(axes[:, 0], names, strict=False):
        marker_i = 0
        card0 = None
        frame_y: list[float] = []
        for s in sorted(panels[panel], key=lambda s: (s.kind, s.path.name)):
            card = s.verification or {}
            card0 = card0 or card
            pts = load_points(s)
            xs = [p.x for p in pts]
            ys = [p.y for p in pts]

            if s.kind == "frame":
                # Off by default. A frame sits outside the labelled ticks,
                # so drawing it forces the axis open and squashes the data
                # into a strip -- and its job (measuring distortion) is
                # already done numerically by the verifier.
                if not show_frames:
                    continue
                ax.plot(xs, ys, ":", color="0.45", lw=1.2,
                        label="panel frame (distortion standard)")
                frame_y.extend(ys)
                continue

            if s.kind == "correlation":
                ax.plot(xs, ys, "--", lw=1.4, alpha=0.9,
                        label=f"drawn: {s.series[:44]}")
                printed = card.get("printed_curve")
                if printed:
                    grid = sorted(xs)
                    ax.plot(grid, _curve(printed, grid), "-", lw=1.1,
                            color="k", alpha=0.55, label="printed equation")
                continue

            ax.plot(xs, ys, MARKERS[marker_i % len(MARKERS)], ms=6,
                    mfc="none", mew=1.3, label=_label(s))
            marker_i += 1

        card0 = card0 or {}
        ax.set_xscale(card0.get("x_axis_type", "log"))
        ax.set_yscale(card0.get("y_axis_type", "log"))
        for setter, key in ((ax.set_xlim, "plot_x_limits"),
                            (ax.set_ylim, "plot_y_limits")):
            lim = card0.get(key)
            if lim:
                lo, hi = lim
                if key == "plot_y_limits" and frame_y:
                    lo = min(lo, min(frame_y) * 0.97)
                    hi = max(hi, max(frame_y) * 1.03)
                setter(lo, hi)
        real = [s for s in panels[panel] if s.kind != "frame"]
        sample = (real or panels[panel])[0]
        ax.set_xlabel(AXIS_NAMES.get(sample.x_axis, sample.x_axis))
        ylabels = {s.y_axis for s in real}
        ax.set_ylabel(
            " or ".join(AXIS_NAMES.get(y, y) for y in sorted(ylabels))
            if ylabels
            else sample.y_axis
        )
        ax.set_title(f"Figure {figure} -- {panel} panel", fontsize=10)
        ax.grid(True, which="both", alpha=0.25, lw=0.5)
        handles, labels = ax.get_legend_handles_labels()
        seen: dict[str, object] = {}
        for h, lab in zip(handles, labels, strict=False):
            seen.setdefault(lab, h)
        ax.legend(seen.values(), seen.keys(), fontsize=6.5, loc="best",
                  ncol=2, framealpha=0.9)

    # Source citation travels with the image, always. Wrapped, because a
    # citation clipped at the page edge is not a citation.
    import textwrap

    lines = ["Redrawn from digitised measurements."]
    lines += textwrap.wrap(f"Source: {chosen[0].source.citation}", 110)
    if chosen[0].after:
        lines += textwrap.wrap(f"After: {chosen[0].after}", 110)
    fig.text(0.5, 0.004, "\n".join(lines), ha="center", va="bottom",
             fontsize=6, color="0.35", linespacing=1.5)

    fig.tight_layout(rect=(0, 0.013 * len(lines) + 0.01, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figure", help='which figure, e.g. "4.46"')
    ap.add_argument("--out", type=Path, help="output image path")
    ap.add_argument("--list", action="store_true", help="list known figures")
    ap.add_argument("--all", action="store_true",
                    help="redraw every figure in the dataset")
    ap.add_argument("--frames", action="store_true",
                    help="draw panel frames too (opens the axes wide)")
    args = ap.parse_args()

    dataset = load_dataset()

    if args.all:
        figures = sorted({s.figure for s in dataset if s.figure})
        for figure in figures:
            out = PLOT_DIR / f"fig{figure}_redrawn.png"
            print(f"wrote {plot_figure(figure, out, dataset, args.frames)}")
        return

    if args.list or not args.figure:
        seen: dict[str, set] = defaultdict(set)
        for s in dataset:
            if s.figure:
                seen[s.figure].add(s.panel)
        for figure in sorted(seen):
            n = sum(1 for s in dataset if s.figure == figure)
            print(f"  {figure:<8} {n:>3} series, panels: "
                  f"{', '.join(sorted(seen[figure]))}")
        return

    out = args.out or PLOT_DIR / f"fig{args.figure}_redrawn.png"
    print(f"wrote {plot_figure(args.figure, out, dataset, args.frames)}")


if __name__ == "__main__":
    main()
