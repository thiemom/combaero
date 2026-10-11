"""Score the joining closures on Blaisdell & Manson 1963 Table 12 (centre laterals).

Pre-registered in data/blaisdell_manson1963/PREREGISTRATION.md; this script
only evaluates what that file froze. Metric: per run |K - zeta|, MAE per
(theta, psi) cell, headline = mean of cell MAEs. No band.

Angle used for a run: the built angle from Table 12 footnote 2 (primary) and
the nominal one (sensitivity); the pre-registration does not say which.

    uv run python validation/junction/bm1963_score.py [--nominal]
"""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validation.junction.joining_variants import (  # noqa: E402
    MODELS,
    bassett_cv,
    bm_theory,
    mynard,
    tee,
)

DATA = Path(__file__).parent / "data" / "blaisdell_manson1963" / "table12_centre.csv"

# Development cells (shown in Bassett 2001), keyed by psi rounded to 2 decimals.
DEV = {1.0: {30, 45, 60, 90, 120, 135, 150}, 2.12: {45, 135}, 4.0: {45, 135}}


def load(nominal: bool) -> list[dict]:
    runs = []
    with DATA.open() as f:
        for r in csv.DictReader(f):
            if "" in (r["Qb_over_Qd"], r["zeta_u"], r["zeta_b"]):
                continue
            th = int(r["theta_deg"])
            psi = float(r["Ad_over_Ab"])
            runs.append(
                {
                    "theta": th,
                    "angle": float(th if nominal else r["theta_actual"]),
                    "psi": psi,
                    "q": float(r["Qb_over_Qd"]),
                    "zu": float(r["zeta_u"]),
                    "zb": float(r["zeta_b"]),
                    "dev": th in DEV.get(psi, set()),
                }
            )
    return runs


def predict(name: str, q: float, psi: float, angle: float) -> tuple[float, float]:
    """(K straight -> collector, K lateral -> collector) for a model."""
    if name == "BM theory (Eqs 12/19)":
        return bm_theory(q, psi, angle)
    U, A, th = tee(q, psi, angle)
    # The kernels need both inlets supplying; at q = 0 or 1 nudge by 1e-9.
    eps = 1e-9
    if q <= 0.0:
        U, A, th = tee(eps, psi, angle)
    elif q >= 1.0:
        U, A, th = tee(1.0 - eps, psi, angle)
    K = MODELS[name](U, A, th)
    return K[1], K[2]


def cell_mae(runs: list[dict], name: str) -> dict[tuple[int, float], tuple[float, float, int]]:
    acc: dict[tuple[int, float], list[tuple[float, float]]] = defaultdict(list)
    for r in runs:
        ku, kb = predict(name, r["q"], r["psi"], r["angle"])
        acc[(r["theta"], r["psi"])].append((ku - r["zu"], kb - r["zb"]))
    return {
        k: (float(np.mean([abs(e[0]) for e in v])), float(np.mean([abs(e[1]) for e in v])), len(v))
        for k, v in acc.items()
    }


def cell_bias(runs: list[dict], name: str) -> dict[tuple[int, float], tuple[float, float]]:
    acc: dict[tuple[int, float], list[tuple[float, float]]] = defaultdict(list)
    for r in runs:
        ku, kb = predict(name, r["q"], r["psi"], r["angle"])
        acc[(r["theta"], r["psi"])].append((ku - r["zu"], kb - r["zb"]))
    return {
        k: (float(np.mean([e[0] for e in v])), float(np.mean([e[1] for e in v])))
        for k, v in acc.items()
    }


def headline(cells: dict, keys: list) -> tuple[float, float, int]:
    ks = [k for k in keys if k in cells]
    if not ks:
        return math.nan, math.nan, 0
    return (
        float(np.mean([cells[k][0] for k in ks])),
        float(np.mean([cells[k][1] for k in ks])),
        len(ks),
    )


def implied_factor(runs: list[dict]) -> float:
    """H3: deflection factor f (d_eff = f * theta) minimising zeta_u MAE under the P2 form."""
    fs = np.linspace(0.0, 1.3, 131)
    best = (math.inf, math.nan)
    for f in fs:
        err = 0.0
        for r in runs:
            q = min(max(r["q"], 1e-9), 1.0 - 1e-9)
            U, A, th = tee(q, r["psi"], r["angle"])
            K = _p2_factor(U, A, th, f)
            err += abs(K[1] - r["zu"])
        if err < best[0]:
            best = (err, float(f))
    return best[1]


def _p2_factor(U, A, th, f: float) -> dict[int, float]:
    """P2 with the jet deflection scaled by f (wall rule from the built geometry, as in P2c)."""
    if f == 1.0:
        return bassett_cv(U, A, th)
    # bassett_cv(correct=True) hard-codes 3/4; reproduce it for any f.
    U, A, th = (np.asarray(x, float) for x in (U, A, th))
    c, sup = 0, [1, 2]
    uc = abs(U[c])
    phi = {j: abs(math.remainder(th[j] - th[c], 2 * math.pi)) for j in sup}
    d = {j: f * (math.pi - phi[j]) for j in sup}
    collinear = {j: phi[j] > math.pi - 1e-9 for j in sup}
    with_flow = {j: math.cos(math.pi - phi[j]) > 1e-12 for j in sup}
    D = A[c] + sum(0.5 * A[j] * math.cos(d[j]) for j in sup if not collinear[j] and with_flow[j])
    M = sum((U[j] / uc) ** 2 * A[j] * math.cos(d[j]) for j in sup)
    dp = 2.0 * (A[c] - M) / D
    return {j: dp + (U[j] / uc) ** 2 - 1.0 for j in sup}


def implied_alpha(runs: list[dict]) -> tuple[float, float]:
    """Exploratory: Mynard alpha (eta 0) minimising zeta_u / zeta_b MAE in a cell. Described, never applied."""
    out = []
    for comp in (1, 2):
        best = (math.inf, math.nan)
        for a in np.linspace(0.0, 1.0, 51):
            err = 0.0
            for r in runs:
                q = min(max(r["q"], 1e-9), 1.0 - 1e-9)
                U, A, th = tee(q, r["psi"], r["angle"])
                K = mynard(U, A, th, float(a), 0.0)
                err += abs(K[comp] - (r["zu"] if comp == 1 else r["zb"]))
            if err < best[0]:
                best = (err, float(a))
        out.append(best[1])
    return out[0], out[1]


def main() -> None:
    nominal = "--nominal" in sys.argv
    runs = load(nominal)
    names = ["shipped", "faithful", "P2 raw", "P2c", "BM theory (Eqs 12/19)"]
    cells = {n: cell_mae(runs, n) for n in names}
    keys = sorted({(r["theta"], r["psi"]) for r in runs})
    held = [k for k in keys if k[0] not in DEV.get(k[1], set())]
    dev = [k for k in keys if k not in held]
    print(
        f"Blaisdell & Manson 1963 Table 12, centre laterals: {len(runs)} runs, {len(keys)} cells "
        f"({len(dev)} dev, {len(held)} same-rig held-out); angle = {'nominal' if nominal else 'built (footnote 2)'}"
    )
    print("headline = mean of per-cell MAE (zeta_u / zeta_b)\n")

    def table(title: str, sel: list) -> None:
        print(f"{title}  [{len(sel)} cells]")
        for n in names:
            u, b, _ = headline(cells[n], sel)
            print(f"  {n:24s} {u:7.3f} / {b:7.3f}")
        print()

    table("All cells", keys)
    table("Development (Bassett-shown)", dev)
    table("Held-out", held)
    h1 = [k for k in held if k[1] >= 4.0 and k[0] <= 90]
    table("H1: held-out, psi >= 4, theta <= 90", h1)
    su, sb, _ = headline(cells["shipped"], h1)
    pu, pb, _ = headline(cells["P2 raw"], h1)
    print(
        f"  H1 (P2 raw < shipped on both): {'SUPPORTED' if pu < su and pb < sb else 'NOT SUPPORTED'} "
        f"(zeta_u {pu:.3f} vs {su:.3f}, zeta_b {pb:.3f} vs {sb:.3f})\n"
    )
    for title, sel in (
        ("H2: psi = 1, all cells", [k for k in keys if k[1] == 1.0]),
        ("H2: psi = 1, held-out cells", [k for k in held if k[1] == 1.0]),
    ):
        table(title, sel)
        su, sb, _ = headline(cells["shipped"], sel)
        pu, pb, _ = headline(cells["P2 raw"], sel)
        print(
            f"  H2 (shipped <= P2 raw on both): {'SUPPORTED' if su <= pu and sb <= pb else 'NOT SUPPORTED'} "
            f"(zeta_u {su:.3f} vs {pu:.3f}, zeta_b {sb:.3f} vs {pb:.3f})\n"
        )

    print("H3: implied deflection factor d_eff/theta (zeta_u, P2 form), theta <= 90")
    psis = sorted({k[1] for k in keys})
    print("  theta " + "".join(f"{p:>8.2f}" for p in psis))
    by = defaultdict(list)
    for r in runs:
        by[(r["theta"], r["psi"])].append(r)
    for th in sorted({k[0] for k in keys if k[0] <= 90}):
        row = [implied_factor(by[(th, p)]) if (th, p) in by else math.nan for p in psis]
        print(f"  {th:5d} " + "".join(f"{v:8.2f}" for v in row))
    print()

    bias = cell_bias(runs, "BM theory (Eqs 12/19)")
    bu = float(np.mean([bias[k][0] for k in keys]))
    bb = float(np.mean([bias[k][1] for k in keys]))
    print(
        f"H4: BM theory mean cell bias (model - data): zeta_u {bu:+.3f}, zeta_b {bb:+.3f} -> "
        f"{'SUPPORTED (biased high)' if bu > 0 and bb > 0 else 'NOT SUPPORTED'}\n"
    )

    print("Per-cell MAE zeta_u / zeta_b (H = held-out)")
    print("  theta   psi  n    " + "".join(f"{n[:14]:>18s}" for n in names))
    for k in keys:
        n_ = cells[names[0]][k][2]
        tag = "H" if k in held else " "
        print(
            f"  {k[0]:5d} {k[1]:5.2f} {n_:3d} {tag} "
            + "".join(f"{cells[n][k][0]:8.3f}/{cells[n][k][1]:<8.3f}  " for n in names)
        )
    if "--alpha" in sys.argv:
        print(
            "\nExploratory: implied Mynard alpha (eta 0) per cell, zeta_u / zeta_b -- described, never applied"
        )
        for k in keys:
            au, ab = implied_alpha(by[k])
            print(f"  {k[0]:5d} {k[1]:5.2f}  {au:5.2f} / {ab:5.2f}")


if __name__ == "__main__":
    main()
