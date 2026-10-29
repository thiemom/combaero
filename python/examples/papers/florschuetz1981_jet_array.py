"""Florschuetz, Truman and Metzger (1981): a jet-array impingement rig, row by row.

Source
    validation/sources.yaml entry
    ``sources/heat_transfer/flohrshuetz/flohrshuetz_truman_metzger_1981.pdf``:
    Florschuetz, L. W., Truman, C. R., Metzger, D. E. (1981). Streamwise flow
    and heat transfer distributions for jet array impingement with crossflow.
    ASME J. Heat Transfer 103, 337-342. Issues #380 (Fig. 6) and #465.

What is built
    The paper's configuration: a jet plate fed from one supply plenum, ten
    spanwise rows of 0.1 in holes impinging on a target, the spent air
    leaving in ONE direction down the gap, so the crossflow builds up row by
    row. ``ImpingementArray`` adds one ImpingementPlateElement and one
    ImpingementCrossflowElement per row; the solver finds the jet flow
    distribution itself, and each row reads its own Gc/Gj from it.

What is reproduced
    1. The crossflow build-up: the network's Gc/Gj at every row beside the
       paper's own closed form, Eq. 8 (crossflow_to_jet_ratio_at_row).
    2. The row-by-row heat transfer: the network's Nu/Nu1 at every row
       beside the measured Nu/Nu1 of Fig. 6, panel B(5,4)I, one series per
       channel height z/d = 1, 2, 3. Each digitised point is matched to the
       row whose Eq. 8 abscissa is nearest (the paper plots one point per
       row, rows 2-10). ``--all`` runs all nine panels (27 geometries).

What is NOT reproduced
    The staggered arrays (Fig. 6 is inline only), the absolute Nu1 of Fig. 5
    (scored in python/tests/test_impingement_plate_validation.py), and the
    paper's rig conditions: Fig. 6 is Nu/Nu1, which the correlation makes
    independent of Re_j, so the supply and exit pressures here only need to
    put Re_j inside the correlation's 2,500-70,000.

Validation label
    FIDELITY (docs/VALIDATION_POLICY.md): the measured points are the
    authors' own data for their own correlation. A miss is ours to explain.
    The z/d = 1 series of this panel is the known weak spot: it sits at the
    correlation's own validity corner (xn/d 5, yn/d 4, z/d 1 all at their
    lower bounds), where the fit under-predicts the crossflow degradation
    (+31% bias on the correlation alone, metadata.yaml). The band on the
    data is +/-5.1% (the paper's +/-5% Nu uncertainty and 0.97% digitising),
    a measurement band, not the model's error.

Data
    The digitised points are committed under
    validation/cooling/data/florschuetz1981/ (rc_<row><col>_*.csv; the scan
    stays out of git). This script READS those files, so it runs from a
    repo checkout, not from an installed wheel; it imports nothing from
    validation/.

The same network, panel B(5,4)I at z/d = 2, opens in the GUI from
Examples -> "Florschuetz 1981 jet array" (gui/examples/florschuetz1981_jet_array.json).
"""

from __future__ import annotations

import argparse
import csv
import sys
from functools import cache
from pathlib import Path
from typing import Any

import combaero as cb
from combaero.network import FlowNetwork, ImpingementArray, NetworkSolver, PressureBoundary

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "validation" / "cooling" / "data" / "florschuetz1981"

D_JET = 0.00254  # Florschuetz's 0.1 in holes
N_ROWS = 10
CD = cb.FLORSCHUETZ_1981_DEFAULT_CD
PR = 0.7  # cancels in Nu/Nu1; only needs to be the same in both calls
# 6 kPa across the array keeps every row of every Fig. 6 geometry inside the
# correlation's Re_j >= 2,500 (3 kPa leaves z/d = 1's first row at ~2,000).
P_SUPPLY = 1.06e5  # [Pa] total, supply plenum
P_EXIT = 1.00e5  # [Pa] total, gap exit
T_SUPPLY = 300.0  # [K]

# Fig. 6 panels by file prefix (validation/cooling/data/florschuetz1981/
# metadata.yaml): rows are xn/d 5, 10, 15; columns yn/d 4, 6, 8.
PANELS = {
    f"rc_{r}{c}": (xn, yn)
    for r, xn in enumerate((5.0, 10.0, 15.0))
    for c, yn in enumerate((4.0, 6.0, 8.0))
}
SERIES = {1.0: "circles_zd_1", 2.0: "squares_zd_2", 3.0: "triangles_zd_3"}


def _boundary(name: str, Pt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, T_SUPPLY, cb.mole_to_mass(cb.species.dry_air())
    return b


@cache
def solve_array(xn_d: float, yn_d: float, z_d: float) -> dict[str, Any]:
    """Build and solve one ten-row inline array; return the per-row summary."""
    arr = ImpingementArray(
        "jets",
        n_rows=N_ROWS,
        d_jet=D_JET,
        xn_d=xn_d,
        yn_d=yn_d,
        z_d=z_d,
        span=12 * yn_d * D_JET,  # twelve holes across, as the GUI default
        plate_thickness=D_JET,
    )
    net = FlowNetwork()
    net.add_node(_boundary("supply", P_SUPPLY))
    net.add_node(_boundary("exit", P_EXIT))
    arr.add_to(net, "supply", "exit")
    res = NetworkSolver(net).solve()
    if not res["__success__"]:
        raise RuntimeError(f"solve failed for ({xn_d}, {yn_d}, {z_d}): {res['__message__']}")
    return arr.summarize(res["__element_diag__"])


def nu_over_nu1(xn_d: float, yn_d: float, z_d: float, re_j: float, gc_gj: float) -> float:
    """Nu at (Re_j, Gc/Gj) over Nu at the same Re_j with no crossflow: the
    quantity Fig. 6 plots. Both through the shipped correlation."""
    s = cb.florschuetz_1981_inline()
    at = cb.jet_array_impingement_nu(s, re_j, gc_gj, PR, xn_d, yn_d, z_d).Nu
    at0 = cb.jet_array_impingement_nu(s, re_j, 0.0, PR, xn_d, yn_d, z_d).Nu
    return at / at0


def measured(prefix: str, z_d: float) -> list[tuple[float, float]]:
    """(Gc/Gj, Nu/Nu1) of one digitised Fig. 6 series."""
    path = DATA / f"{prefix}_{SERIES[z_d]}.csv"
    with path.open() as f:
        rows = list(csv.reader(f))[1:]
    return [(float(x), float(y)) for x, y in rows]


def compare(prefix: str, z_d: float, verbose: bool) -> list[float]:
    """Solve one geometry, print the row table, return per-point errors."""
    xn, yn = PANELS[prefix]
    summ = solve_array(xn, yn, z_d)
    gc, re_j, nu = summ["rows_Gc_Gj"], summ["rows_Re_j"], summ["rows_Nu"]
    eq8 = [cb.crossflow_to_jet_ratio_at_row(yn, z_d, CD, i) for i in range(1, N_ROWS + 1)]
    net_ratio = [nu_over_nu1(xn, yn, z_d, re_j[i], gc[i]) for i in range(N_ROWS)]

    # Each digitised point -> the row whose Eq. 8 abscissa is nearest. Where
    # Eq. 8 flattens (z/d = 1, last rows) two points can share a row.
    by_row: dict[int, list[float]] = {i: [] for i in range(N_ROWS)}
    for x, y in measured(prefix, z_d):
        by_row[min(range(N_ROWS), key=lambda i: abs(eq8[i] - x))].append(y)
    errors = [net_ratio[r] / y - 1.0 for r, ys in by_row.items() for y in ys]

    if verbose:
        print(
            f"\n(xn/d, yn/d, z/d) = ({xn:g}, {yn:g}, {z_d:g})   "
            f"Re_j {min(re_j):.0f}-{max(re_j):.0f}"
        )
        print("row  Gc/Gj net  Gc/Gj Eq.8   Nu net  Nu/Nu1 net  Nu/Nu1 Fig.6   error")
        for i in range(N_ROWS):
            head = f"{i + 1:3d}  {gc[i]:9.3f}  {eq8[i]:10.3f}  {nu[i]:7.1f}  {net_ratio[i]:10.3f}"
            if not by_row[i]:
                print(head)
            for k, y in enumerate(by_row[i]):
                lead = head if k == 0 else " " * len(head)
                print(f"{lead}  {y:12.3f}  {net_ratio[i] / y - 1.0:+6.1%}")
        bias = sum(errors) / len(errors)
        mae = sum(abs(e) for e in errors) / len(errors)
        print(f"     {len(errors)} points: bias {bias:+.1%}, MAE {mae:.1%}")
    return errors


def _plot(prefix: str) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("(matplotlib not installed: pip install combaero[examples] for the plot)")
        return
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from plot_utils import show_or_save

    xn, yn = PANELS[prefix]
    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    for z_d, marker in zip(SERIES, ("o", "s", "^"), strict=True):
        summ = solve_array(xn, yn, z_d)
        gc, re_j = summ["rows_Gc_Gj"], summ["rows_Re_j"]
        ratio = [nu_over_nu1(xn, yn, z_d, re_j[i], gc[i]) for i in range(N_ROWS)]
        (line,) = ax.plot(gc, ratio, ".-", label=f"network, z/d = {z_d:g}")
        xs, ys = zip(*measured(prefix, z_d), strict=True)
        ax.plot(xs, ys, marker, mfc="none", color=line.get_color(), label=f"Fig. 6, z/d = {z_d:g}")
    ax.set_xlabel("Gc/Gj")
    ax.set_ylabel("Nu/Nu1")
    ax.set_title(f"Florschuetz et al. 1981, Fig. 6 panel B({xn:g},{yn:g})I -- fidelity")
    ax.set_xlim(0.0, 0.8)
    ax.set_ylim(0.4, 1.05)
    ax.legend(fontsize=8)
    show_or_save(fig, "florschuetz1981_jet_array.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--all", action="store_true", help="all nine Fig. 6 panels")
    parser.add_argument("--no-plot", action="store_true", help="skip the panel B(5,4)I plot")
    args = parser.parse_args()

    if not DATA.is_dir():
        sys.exit(f"Digitised data not found at {DATA}: run from a repo checkout.")

    print("Florschuetz, Truman and Metzger (1981), Fig. 6 -- FIDELITY (the authors' own data)")
    print(f"Inline arrays, {N_ROWS} rows, d = {D_JET * 1000:.2f} mm, Cd = {CD}")

    panel = "rc_00"  # B(5,4)I
    errors = [e for z in SERIES for e in compare(panel, z, verbose=True)]
    n = len(errors)
    print(
        f"\nPanel B(5,4)I, {n} points: bias {sum(errors) / n:+.1%}, "
        f"MAE {sum(abs(e) for e in errors) / n:.1%}"
    )

    if args.all:
        pooled = [e for p in PANELS for z in SERIES for e in compare(p, z, verbose=False)]
        n = len(pooled)
        print(
            f"All nine panels, {n} points: bias {sum(pooled) / n:+.1%}, "
            f"MAE {sum(abs(e) for e in pooled) / n:.1%}"
        )

    # The crossflow build-up against Eq. 8 for the GUI example's geometry.
    summ = solve_array(5.0, 4.0, 2.0)
    worst = max(
        abs(summ["rows_Gc_Gj"][i - 1] / cb.crossflow_to_jet_ratio_at_row(4.0, 2.0, CD, i) - 1.0)
        for i in range(2, N_ROWS + 1)
    )
    print(f"Gc/Gj at (5, 4, 2), rows 2-10: worst |network/Eq. 8 - 1| = {worst:.1%}")
    m = summ["rows_m_dot"]
    print(f"Jet flow last/first row: {m[-1] / m[0]:.3f} (crossflow pulls flow to the exit end)")

    if not args.no_plot:
        _plot(panel)


if __name__ == "__main__":
    main()
