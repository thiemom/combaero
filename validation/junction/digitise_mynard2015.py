"""Digitise Mynard & Valen-Sendstad 2015, Figs 4 and 6-11, by colour.

The figures are 150-dpi raster images in the accepted manuscript
(docs/junction/a-unified-method-...pdf, gitignored: copyright). Each panel
plots Mynard's own model ("Unified0D", a red line) and his 3D CFD reference
("Ref3D", black dots on a black line). Both are separable by colour, so the
extraction is automatic and repeatable:

- **Calibration** from the tick LABELS: their glyph centroids, fitted linearly
  against the values printed on the axis. The labels are right-aligned to the
  y axis and centred under the x ticks; stray glyphs (the rotated axis title,
  a neighbouring panel's "0") are dropped by keeping the most evenly spaced
  subset. The residual of that fit is recorded per panel.
- **Unified0D** as centre points of the red band: column scans where the curve
  is shallow, row scans where it is steep. A band touching a black pixel is
  dropped, because the Ref3D line is drawn over it and a half-hidden band
  biases its centre. Fig 4 carries two red curves (Mynard's model and the same
  with eta_j = 0, dashed); they are followed by continuity from the right edge
  and dropped where they cross.
- **Ref3D** as the centroids of the filled markers: a 5x5 erosion removes the
  2-px line and leaves the ~9-px dots.

Run (needs poppler's ``pdfimages``):

    uv run python validation/junction/digitise_mynard2015.py [--overlay DIR]

It rewrites ``data/mynard2015/*.csv`` and the calibration block of
``data/mynard2015/calibration.yaml``. ``--overlay`` writes review images (the
scan with the calibrated ticks, the red points and the dots marked) -- those
contain the copyrighted scan and must never be committed.
"""

from __future__ import annotations

import argparse
import itertools
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml
from PIL import Image, ImageDraw
from scipy import ndimage

_REPO = Path(__file__).resolve().parents[2]
_PDF = _REPO / "docs/junction/a-unified-method-for-estimating-pressure-losses-at-vascular-3elbvn7sno.pdf"
_OUT = Path(__file__).resolve().parent / "data/mynard2015"


@dataclass(frozen=True)
class Panel:
    key: str
    image: int  # pdfimages global image number
    y_axis_col: int  # approximate pixel column of the y axis
    x_axis_row: int  # approximate pixel row of the x axis
    xticks: tuple[float, ...]
    yticks: tuple[float, ...]
    n_red: int = 1
    exclude: tuple[tuple[int, int, int, int], ...] = ()  # legend boxes, x0 y0 x1 y1
    ylab_w: int = 45
    extra: dict = field(default_factory=dict)


def _r(a: float, b: float, s: float) -> tuple[float, ...]:
    return tuple(round(a + i * s, 6) for i in range(round((b - a) / s) + 1))


_RE2000, _RE2500, _FR, _ANG = _r(0, 2000, 500), _r(0, 2500, 500), _r(0, 1, 0.2), (0.0, 50.0, 100.0, 150.0)

PANELS: tuple[Panel, ...] = (
    Panel("fig04", 3, 54, 287, _r(0, 90, 15), _r(-0.2, 1.2, 0.2), n_red=2, ylab_w=30, exclude=((55, 0, 200, 90),)),
    Panel("fig06a", 5, 59, 209, _RE2000, _r(0, 2.5, 0.5), exclude=((62, 30, 205, 100),)),
    Panel("fig06b", 5, 252, 209, _FR, _r(0, 2.5, 0.5)),
    Panel("fig06c", 5, 446, 209, _ANG, _r(0, 2.5, 0.5)),
    Panel("fig06d", 5, 640, 209, _FR, _r(-20, 40, 10)),
    Panel("fig06e", 5, 61, 464, _RE2000, _r(0, 2, 0.5)),
    Panel("fig06f", 5, 255, 464, _FR, _r(0, 2, 0.5)),
    Panel("fig06g", 5, 448, 464, _ANG, _r(0, 2, 0.5)),
    Panel("fig06h", 5, 642, 464, _FR, _r(0, 2, 0.5)),
    Panel("fig07a", 6, 43, 242, _RE2000, _r(-0.2, 1.8, 0.2), exclude=((46, 28, 200, 110),)),
    Panel("fig07b", 6, 242, 242, _FR, _r(-0.2, 1.8, 0.2)),
    Panel("fig07c", 6, 442, 242, _ANG, _r(-0.2, 1.8, 0.2)),
    Panel("fig07d", 6, 641, 242, _FR, _r(-0.2, 1.8, 0.2)),
    Panel("fig08a", 7, 39, 238, _RE2000, _r(0, 1.6, 0.2), exclude=((42, 28, 220, 125),)),
    Panel("fig08b", 7, 306, 238, _FR, _r(0, 1.6, 0.2)),
    Panel("fig08c", 7, 573, 238, _r(0, 120, 20), _r(0, 1.6, 0.2)),
    Panel("fig09a", 8, 65, 234, _RE2000, _r(-1, 1.5, 0.5), exclude=((68, 30, 205, 97),)),
    Panel("fig09b", 8, 256, 234, _FR, _r(-1, 1.5, 0.5)),
    Panel("fig09c", 8, 447, 234, _ANG, _r(-1, 1.5, 0.5)),
    Panel("fig09d", 8, 637, 234, _FR, _r(-1, 1.5, 0.5)),
    Panel("fig09e", 8, 65, 511, _RE2000, _r(-1, 1.5, 0.5)),
    Panel("fig09f", 8, 256, 511, _FR, _r(-1, 1.5, 0.5)),
    Panel("fig09g", 8, 447, 511, _ANG, _r(-1, 1.5, 0.5)),
    Panel("fig09h", 8, 637, 511, _FR, _r(-1, 1.5, 0.5)),
    Panel("fig10a", 9, 40, 242, _RE2500, _r(-1, 2.5, 0.5), exclude=((43, 30, 200, 100),)),
    Panel("fig10b", 9, 238, 242, _FR, _r(-1, 2.5, 0.5)),
    Panel("fig10c", 9, 436, 242, _ANG, _r(-1, 2.5, 0.5)),
    Panel("fig10d", 9, 634, 242, _FR, _r(-1, 2.5, 0.5)),
    Panel("fig11a", 10, 40, 233, _RE2500, _r(-1, 1.5, 0.5), exclude=((43, 25, 215, 100),)),
    Panel("fig11b", 10, 307, 233, _FR, _r(-1, 1.5, 0.5)),
    Panel("fig11c", 10, 575, 233, _r(0, 120, 15), _r(-1, 1.5, 0.5)),
)


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------


def _axis_mask(a: np.ndarray) -> np.ndarray:
    return (a.mean(axis=2) < 235) & ((a.max(axis=2) - a.min(axis=2)) < 30)


def _axes(a: np.ndarray, p: Panel) -> dict:
    """The axis lines: the y-axis run that reaches the x axis, and the x-axis
    run that starts at the y axis. Gaps of a few pixels (a marker or a curve
    crossing the axis) are bridged."""
    m = _axis_mask(a)
    g = a.mean(axis=2)
    best = None
    for c in range(p.y_axis_col - 4, p.y_axis_col + 5):
        col = m[:, c]
        r = p.x_axis_row + 2
        while r > p.x_axis_row - 6 and not col[r]:
            r -= 1
        if not col[r]:
            continue
        end = r
        while r - 9 >= 0 and col[r - 9 : r].any():
            r -= 1
        if best is None or end - r > best[0]:
            best = (end - r, c, r)
    _, ycol, top = best
    rows = slice(top + 5, p.x_axis_row - 5)
    w = np.array([(255 - g[rows, c]).sum() for c in range(ycol - 2, ycol + 3)])
    y_axis = float(np.dot(w, np.arange(ycol - 2, ycol + 3)) / w.sum())
    best = None
    for r in range(p.x_axis_row - 4, p.x_axis_row + 5):
        c = p.y_axis_col - 3
        while c < p.y_axis_col + 12 and not m[r, c]:
            c += 1
        start = c
        while c + 4 < m.shape[1] and m[r, c + 1 : c + 5].any():
            c += 1
        if best is None or c - start > best[0]:
            best = (c - start, r, c)
    _, xrow, right = best
    cols = slice(ycol + 5, right - 5)
    w = np.array([(255 - g[r, cols]).sum() for r in range(xrow - 2, xrow + 3)])
    x_axis = float(np.dot(w, np.arange(xrow - 2, xrow + 3)) / w.sum())
    return dict(y_axis=y_axis, ycol=ycol, top=top, x_axis=x_axis, xrow=xrow, right=right)


def _clusters(mask: np.ndarray, axis: int, gap: int) -> list[float]:
    lab, n = ndimage.label(mask)
    spans = sorted((sl[axis].start, sl[axis].stop) for sl in ndimage.find_objects(lab))
    groups: list[list[int]] = []
    for lo, hi in spans:
        if groups and lo <= groups[-1][1] + gap:
            groups[-1][1] = max(groups[-1][1], hi)
        else:
            groups.append([lo, hi])
    return [(lo + hi - 1) / 2.0 for lo, hi in groups]


def _evenly_spaced(v: list[float], n: int) -> list[float]:
    if len(v) <= n:
        return v
    best = None
    for keep in itertools.combinations(range(len(v)), n):
        y = np.array([v[i] for i in keep])
        res = y - np.polyval(np.polyfit(np.arange(n), y, 1), np.arange(n))
        if best is None or np.abs(res).max() < best[0]:
            best = (float(np.abs(res).max()), list(y))
    return best[1]


def calibrate(a: np.ndarray, p: Panel) -> dict:
    ax = _axes(a, p)
    dark = (a.max(axis=2) < 150) & ((a.max(axis=2) - a.min(axis=2)) < 40)
    y0, y1 = max(0, ax["top"] - 8), ax["xrow"] + 8
    x0, x1 = max(0, ax["ycol"] - p.ylab_w), ax["ycol"] - 3
    win = dark[y0:y1, x0:x1].copy()
    lab, _ = ndimage.label(win)
    for i, sl in enumerate(ndimage.find_objects(lab)):
        if sl[1].stop < win.shape[1] - 7:  # y labels are right-aligned to the axis
            win[lab == i + 1] = False
    ylab = _evenly_spaced([y0 + v for v in _clusters(win, 0, 2)], len(p.yticks))
    r0, c0 = ax["xrow"] + 3, max(0, ax["ycol"] - 12)
    xlab = _clusters(dark[r0 : ax["xrow"] + 18, c0 : ax["right"] + 12], 1, 5)
    xlab = _evenly_spaced([c0 + v for v in xlab], len(p.xticks))
    if len(ylab) != len(p.yticks) or len(xlab) != len(p.xticks):
        raise RuntimeError(f"{p.key}: {len(xlab)} x labels / {len(ylab)} y labels found")
    yv = sorted(p.yticks, reverse=True)  # top to bottom
    ky, by = np.polyfit(ylab, yv, 1)
    kx, bx = np.polyfit(xlab, p.xticks, 1)
    res_y = np.abs(np.array(yv) - (ky * np.array(ylab) + by)).max() / abs(ky)
    res_x = np.abs(np.array(p.xticks) - (kx * np.array(xlab) + bx)).max() / abs(kx)
    return dict(ax, kx=float(kx), bx=float(bx), ky=float(ky), by=float(by), res_x=float(res_x), res_y=float(res_y))


# ---------------------------------------------------------------------------
# Curves and markers
# ---------------------------------------------------------------------------


def _runs(idx: list[int]) -> list[list[int]]:
    runs: list[list[int]] = []
    for r in idx:
        if runs and r == runs[-1][1] + 1:
            runs[-1][1] = r
        else:
            runs.append([r, r])
    return runs


def red_tracks(a: np.ndarray, p: Panel, cal: dict) -> tuple[list[list[tuple[float, float]]], int]:
    R, G, B = a[..., 0], a[..., 1], a[..., 2]
    red = (R > 170) & (G < 120) & (B < 120) & (R - G > 90)
    black = a.max(axis=2) < 110
    for x0, y0, x1, y1 in p.exclude:
        red[y0 : y1 + 1, x0 : x1 + 1] = False
    c0, c1 = int(cal["y_axis"]) + 2, cal["right"]
    r0, r1 = cal["top"], int(cal["x_axis"]) - 1
    cols = {c: _runs([r for r in range(r0, r1) if red[r, c]]) for c in range(c0, c1 + 1)}
    thick = [hi - lo + 1 for runs in cols.values() for lo, hi in runs if hi > lo]
    nominal = int(np.median(thick)) if thick else 0
    pts_col: dict[int, list[float]] = {}
    steep = np.zeros_like(red)
    for c, runs in cols.items():
        good = []
        for lo, hi in runs:
            t = hi - lo + 1
            if t < 2:
                continue
            if t > nominal + 1:
                steep[lo : hi + 1, c] = True
                continue
            if black[lo - 1, c] or black[hi + 1, c]:
                continue
            good.append((lo + hi) / 2.0)
        pts_col[c] = good
    pts_row = []
    for r in range(r0, r1):
        for lo, hi in _runs([c for c in range(c0, c1 + 1) if red[r, c]]):
            if not 2 <= hi - lo + 1 <= nominal + 1 or not steep[r, lo : hi + 1].any():
                continue
            if black[r, lo - 1] or black[r, hi + 1]:
                continue
            pts_row.append(((lo + hi) / 2.0, float(r)))
    if p.n_red == 1:
        pts = [(float(c), ys[0]) for c, ys in pts_col.items() if len(ys) == 1] + pts_row
        return [sorted(pts)], nominal
    tracks: list[list[tuple[float, float]]] = [[], []]
    last: list[float | None] = [None, None]
    sep = 4 * nominal + 4  # closer than this, which curve is which is a guess
    for c in sorted(pts_col, reverse=True):
        ys = pts_col[c]
        if len(ys) == 2 and abs(ys[0] - ys[1]) > sep:
            pair = ys
            if last[0] is not None and abs(ys[1] - last[0]) + abs(ys[0] - last[1]) < abs(ys[0] - last[0]) + abs(
                ys[1] - last[1]
            ):
                pair = ys[::-1]
            for k in (0, 1):
                tracks[k].append((float(c), pair[k]))
                last[k] = pair[k]
        elif len(ys) == 1 and last[0] is not None:
            d = [abs(ys[0] - last[0]), abs(ys[0] - last[1])]
            k = int(np.argmin(d))
            if d[k] <= 2.0 and d[1 - k] > sep:
                tracks[k].append((float(c), ys[0]))
                last[k] = ys[0]
    return tracks, nominal


def dots(a: np.ndarray, p: Panel, cal: dict) -> list[tuple[float, float]]:
    black = a.max(axis=2) < 110
    c0, c1 = int(cal["y_axis"]) - 6, cal["right"] + 6
    r0, r1 = cal["top"] - 6, int(cal["x_axis"]) + 6
    st = np.ones((5, 5), bool)
    st[0, 0] = st[0, 4] = st[4, 0] = st[4, 4] = False
    core = ndimage.binary_erosion(black[r0:r1, c0:c1], structure=st)
    lab, n = ndimage.label(core)
    found: list[tuple[float, float]] = []
    for i in range(1, n + 1):
        ys, xs = np.nonzero(lab == i)
        cx, cy = float(xs.mean() + c0), float(ys.mean() + r0)
        if any(x0 <= cx <= x1 and y0 <= cy <= y1 for x0, y0, x1, y1 in p.exclude):
            continue
        near = [j for j, (fx, fy) in enumerate(found) if (fx - cx) ** 2 + (fy - cy) ** 2 < 25.0]
        if near:  # one marker eroded into two cores
            fx, fy = found[near[0]]
            found[near[0]] = ((fx + cx) / 2.0, (fy + cy) / 2.0)
        else:
            found.append((cx, cy))
    return sorted(found)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _images(workdir: Path) -> dict[int, Path]:
    subprocess.run(["pdfimages", "-png", str(_PDF), str(workdir / "img")], check=True)
    return {int(f.stem.split("-")[1]): f for f in workdir.glob("img-*.png")}


def _write_csv(path: Path, header: str, rows: list[tuple[float, float]]) -> None:
    lines = [header, "x, K"] + [f"{x:.5f}, {k:.5f}" for x, k in rows]
    path.write_text("\n".join(lines) + "\n")


def run(overlay_dir: Path | None) -> None:
    calib: dict[str, dict] = {}
    with tempfile.TemporaryDirectory() as tmp:
        imgs = _images(Path(tmp))
        for p in PANELS:
            a = np.asarray(Image.open(imgs[p.image]).convert("RGB")).astype(int)
            cal = calibrate(a, p)
            tracks, nominal = red_tracks(a, p, cal)
            ds = dots(a, p, cal)

            def data(c: float, r: float, cal: dict = cal) -> tuple[float, float]:
                return cal["kx"] * c + cal["bx"], cal["ky"] * r + cal["by"]

            names = ["unified0d"]
            if p.n_red == 2:
                # Name the two curves by how they are DRAWN, not by the model:
                # the dashed one (eta_j = 0 in the legend) repeats gaps of the
                # dash period, 5-12 px; the solid one is only interrupted
                # where the Ref3D line hides it (2-3 px, or one long stretch).
                def dash_gaps(t: list[tuple[float, float]]) -> int:
                    d = np.diff(sorted({round(c) for c, _ in t}))
                    return int(((d >= 5) & (d <= 12)).sum())

                gaps = [dash_gaps(t) for t in tracks]
                tracks = [tracks[int(np.argmin(gaps))], tracks[int(np.argmax(gaps))]]
                names = ["unified0d", "unified0d_eta0"]
                print(f"  {p.key}: dash-period gaps per track {gaps}; the solid one has fewer")
            for name, t in zip(names, tracks, strict=True):
                _write_csv(
                    _OUT / f"mynard_{p.key}_{name}.csv",
                    f"# Mynard 2015 {p.key}: red {name} curve, digitise_mynard2015.py",
                    [data(c, r) for c, r in t],
                )
            _write_csv(
                _OUT / f"mynard_{p.key}_ref3d.csv",
                f"# Mynard 2015 {p.key}: Ref3D markers, digitise_mynard2015.py",
                [data(c, r) for c, r in ds],
            )
            calib[p.key] = dict(
                px_per_x=round(1.0 / abs(cal["kx"]), 4),
                px_per_K=round(1.0 / abs(cal["ky"]), 4),
                calib_residual_px=[round(cal["res_x"], 2), round(cal["res_y"], 2)],
                red_band_px=nominal,
                n_unified0d=[len(t) for t in tracks],
                n_ref3d=len(ds),
            )
            print(f"{p.key}: {calib[p.key]}")
            if overlay_dir is not None:
                _overlay(imgs[p.image], p, cal, tracks, ds, overlay_dir)
    header = (
        "# Written by validation/junction/digitise_mynard2015.py -- do not edit.\n"
        "# px_per_x / px_per_K: the scan's resolution per data unit, used to state\n"
        "# fidelity agreement in pixels. calib_residual_px: worst tick-label\n"
        "# misfit of the linear calibration (x, y).\n"
    )
    (_OUT / "calibration.yaml").write_text(header + yaml.safe_dump(calib, sort_keys=True))


def _overlay(img: Path, p: Panel, cal: dict, tracks: list, ds: list, out: Path) -> None:
    S = 4
    im = Image.open(img).convert("RGB")
    box = (int(cal["y_axis"]) - 50, max(0, cal["top"] - 12), cal["right"] + 12, int(cal["x_axis"]) + 22)
    im = im.crop(box).resize(((box[2] - box[0]) * S, (box[3] - box[1]) * S), Image.NEAREST)
    d = ImageDraw.Draw(im)

    def px(c: float, r: float) -> tuple[float, float]:
        return (c - box[0] + 0.5) * S, (r - box[1] + 0.5) * S

    for v in p.yticks:
        x, y = px(cal["y_axis"], (v - cal["by"]) / cal["ky"])
        d.line([(x - 16, y), (x + 16, y)], fill=(0, 140, 255), width=2)
    for v in p.xticks:
        x, y = px((v - cal["bx"]) / cal["kx"], cal["x_axis"])
        d.line([(x, y - 16), (x, y + 16)], fill=(0, 140, 255), width=2)
    for k, t in enumerate(tracks):
        col = [(0, 230, 230), (255, 0, 255)][k]
        for c, r in t:
            x, y = px(c, r)
            d.ellipse([x - 2, y - 2, x + 2, y + 2], fill=col)
    for c, r in ds:
        x, y = px(c, r)
        d.line([(x - 10, y - 10), (x + 10, y + 10)], fill=(0, 200, 0), width=3)
        d.line([(x - 10, y + 10), (x + 10, y - 10)], fill=(0, 200, 0), width=3)
    im.save(out / f"{p.key}.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlay", type=Path, default=None)
    args = ap.parse_args()
    if args.overlay is not None:
        args.overlay.mkdir(parents=True, exist_ok=True)
    _OUT.mkdir(parents=True, exist_ok=True)
    run(args.overlay)
