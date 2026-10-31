"""Hylton et al. (1983): the C3X internally cooled turbine vane, station by station.

Source
    validation/sources.yaml entry
    ``sources/heat_transfer/vanes/hylton1983_NASA_CR-168015_C3X_MarkII_vane_heat_transfer.pdf``:
    Hylton, L. D., Mihelc, M. S., Turner, E. R., Nealy, D. A., York, R. E.
    (1983). Analytical and experimental evaluation of the heat transfer
    distribution over the surfaces of turbine vanes. NASA CR-168015.

What is built
    A linear cascade of C3X vanes (ASTM 310 stainless steel), cooled by air
    through ten radial holes, in hot air at about 800 K and exit Mach 0.9.
    The vane surface is split into one segment per hole: each surface sample
    belongs to the hole nearest to it, giving 19 segments (hole 1 has no
    pressure-side segment). Per segment:

    - Gas side: a ChannelElement from the inlet total state (Pt1, Tt1) to the
      MEASURED local static pressure, with the passage's hydraulic diameter.
      The pressure distribution is imposed, not predicted: this network has
      no airfoil boundary-layer model, and a passage chain runs too close to
      sonic on the suction side to stand in for one.
    - Coolant side: the hole's share of its flow and perimeter (f = L / L_hole
      of the segment), as a ChannelElement of the same velocity and Dh
      (diameter D sqrt(f), Dh = D). The inlet temperature is set so that the
      midspan temperature is the report's average hole temperature.
    - Wall: a ThermalWall of ss310 (Gen3 CSP, #881), as thick as the mean
      distance from the segment to the hole's perimeter. A ThermalWall has
      one area -- the segment's -- so the hole side's smaller area enters as
      the factor f pi D / L on its h (``Nu_multiplier`` of the coolant
      sub-channel). That factor is geometry, not a tuner (per-side wall
      areas would carry it directly: #887).

What is reproduced
    1. Untuned: the measured midspan wall temperature Tw of each segment
       against the network's, with smooth-channel (Gnielinski) heat transfer
       on both sides. Run 158 (Tu 8.3 %) and run 149 (Tu 6.5 %).
    2. Tuned: one gas-side ``Nu_multiplier`` per segment, fitted so that the
       network reproduces run 158's measured Tw exactly -- 19 multipliers for
       19 targets, no degrees of freedom left. This is the workflow a user
       follows to match a rig.
    3. Transfer: run 149 solved with run 158's multipliers, which is what the
       tuning is worth on a second condition of the same rig.

What is NOT reproduced
    - **The report's external h is not a target.** It is inferred from the
      measured Tw by a 2-D finite-element conduction model with a
      correlation on the hole side (Nu = Cr 0.022 Pr^0.5 Re^0.8), so
      matching it would be circular. It is printed for comparison only.
    - The leading-edge stagnation flow, transition and free-stream
      turbulence: none is in a channel correlation. The tuners absorb them,
      which is why they differ between runs of different Tu.
    - 2-D conduction around each hole, and conduction between segments.
    - The spanwise variation: the data are midspan.
    - A segment's duct starts at a reservoir, and its heat-transfer reference
      is that reservoir's state, so the recovery temperature on the
      near-sonic suction side is about 5 K high (#888). That is a few kelvin
      of Tw, small beside the untuned error, and the tuners absorb it.

Validation label
    TUNED (docs/VALIDATION_POLICY.md: matching a specific rig is the user's
    job, with ``Nu_multiplier``). The untuned errors are CROSS-SOURCE
    accuracy: a smooth-channel correlation that never saw a vane. Run 149
    with run 158's multipliers is a same-rig transfer, not an independent
    test.

Data
    The transcribed tables are committed under
    validation/cooling/data/hylton1983/ (extraction record:
    validation/cooling/extractions/hylton_1983_cr168015_c3x.md; the report
    itself stays out of git). This script READS those files, so it runs
    from a repo checkout, not from an installed wheel; it imports nothing
    from validation/.

The tuned run-158 network opens in the GUI from Examples -> "Hylton 1983 C3X
vane" (gui/examples/hylton1983_c3x_vane.json); ``--write-gui`` regenerates
that file from this script.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import warnings
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "validation" / "cooling" / "data" / "hylton1983"
GUI_FILE = REPO / "gui" / "examples" / "hylton1983_c3x_vane.json"

SPAN = 0.0762  # [m] vane height, Table IV
ARC = {"SS": 0.17782, "PS": 0.13723}  # [m] surface arcs, Table IV
N_SAMPLES = 400  # surface samples per side for the hole assignment
P_COOLANT_EXIT = 2.0e5  # [Pa] hole exit plenum; the hole flows are imposed
MATERIAL = "ss310"
TUNE_RUN = 158

# Gas-side Nu multipliers fitted to run 158's measured Tw (``--tune``
# regenerates them; test_paper_examples pins that they still reproduce it).
TUNED = {
    "PS0_h2": 2.2380,
    "PS1_h3": 1.8266,
    "PS2_h4": 1.4553,
    "PS3_h5": 1.4799,
    "PS4_h6": 1.2201,
    "PS5_h7": 1.4412,
    "PS6_h8": 1.4107,
    "PS7_h9": 1.0367,
    "PS8_h10": 0.9023,
    "SS0_h2": 2.1936,
    "SS1_h1": 0.8686,
    "SS2_h3": 1.1570,
    "SS3_h4": 1.0934,
    "SS4_h5": 1.1658,
    "SS5_h6": 1.0407,
    "SS6_h7": 1.1157,
    "SS7_h8": 1.1328,
    "SS8_h9": 0.9416,
    "SS9_h10": 0.8034,
}


def _rows(name: str) -> list[dict[str, str]]:
    with open(DATA / name) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def _air() -> list[float]:
    return list(cb.mole_to_mass(cb.species.dry_air()))


def _sides() -> dict[str, np.ndarray]:
    """Each surface as a polyline from the leading-edge origin (minimum x),
    Table III: points 1-29 run over the suction side, 78 back to 40 the
    pressure side."""
    xy = np.array([[float(p["x_cm"]), float(p["y_cm"])] for p in _rows("c3x_coordinates.csv")])
    xy = xy / 100.0
    origin = xy[int(np.argmin(xy[:, 0]))]
    return {"SS": np.vstack([origin, xy[0:29]]), "PS": np.vstack([origin, xy[77:38:-1]])}


def _resample(poly: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    seg = np.hypot(*np.diff(poly, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    t = np.linspace(0.0, s[-1], n)
    return np.column_stack([np.interp(t, s, poly[:, 0]), np.interp(t, s, poly[:, 1])]), t


def run_conditions(run: int) -> dict[str, float]:
    """Table IX: inlet total pressure (from the psia column) and temperature, Tu."""
    r = {int(c["run"]): c for c in _rows("c3x_run_conditions.csv")}[run]
    return {
        "Pt1": float(r["PT1_kPa_from_psia"]) * 1e3,
        "Tt1": float(r["TT1_K"]),
        "Tu": float(r["Tu_pct"]),
        "M2": float(r["M2"]),
    }


def _measured(table: list[dict[str, str]], key: str, side: str, frac: np.ndarray) -> np.ndarray:
    """A measured column of one surface (the stagnation tap included),
    interpolated linearly in s/arc; a station a run does not list is skipped."""
    p = sorted(
        (float(r["s_over_arc"]), float(r[key]))
        for r in table
        if r["surface"] in (side, "LE") and r[key].strip()
    )
    return np.interp(frac, [q[0] for q in p], [q[1] for q in p])


@cache
def segments(run: int) -> tuple[dict[str, Any], ...]:
    """The hole-assigned surface segments with their measured averages."""
    holes = _rows("c3x_cooling_holes_xy.csv")
    centre = np.array([[float(h["x_cm"]), float(h["y_cm"])] for h in holes]) / 100.0
    radius = np.array([float(h["D_cm"]) for h in holes]) / 200.0
    tw = _rows(f"c3x_run{run}_Tw_h.csv")
    ps = _rows(f"c3x_run{run}_Ps.csv")
    passage = _rows("c3x_passage_geometry.csv")
    dh_s = [float(r["s_over_arc_ss"]) for r in passage]
    dh = [float(r["Dh_cm"]) / 100.0 for r in passage]
    out = []
    for side, poly in _sides().items():
        pts, s = _resample(poly, N_SAMPLES)
        frac = s / ARC[side]
        dist = (
            np.hypot(pts[:, None, 0] - centre[None, :, 0], pts[:, None, 1] - centre[None, :, 1])
            - radius[None, :]
        )
        nearest = np.argmin(dist, axis=1)

        Tw = _measured(tw, "Tw_K", side, frac)
        h_rep = _measured(tw, "h_W_m2K", side, frac)
        ps_ratio = _measured(ps, "Ps_over_Pt1", side, frac)
        start, j = 0, 0
        for k in range(1, N_SAMPLES + 1):
            if k < N_SAMPLES and nearest[k] == nearest[start]:
                continue
            idx = slice(start, k)
            i = int(nearest[start])
            mid = 0.5 * (frac[start] + frac[k - 1])
            out.append(
                {
                    "id": f"{side}{j}_h{i + 1}",
                    "side": side,
                    "hole": i + 1,
                    "s0": float(frac[start]),
                    "s1": float(frac[k - 1]),
                    "L": float(s[k - 1] - s[start] + s[1] - s[0]),
                    "thickness": float(np.mean(dist[idx, i])),
                    "Dh": float(np.interp(mid, dh_s, dh)),
                    "Ps_over_Pt1": float(np.interp(mid, frac, ps_ratio)),
                    "Tw_measured": float(np.mean(Tw[idx])),
                    "h_report": float(np.mean(h_rep[idx])),
                }
            )
            start, j = k, j + 1
    # The share of each hole's flow and perimeter that faces each segment.
    for seg in out:
        same = [g for g in out if g["hole"] == seg["hole"]]
        seg["f"] = seg["L"] / sum(g["L"] for g in same)
    return tuple(out)


def coolant(run: int) -> dict[int, dict[str, float]]:
    """Appendix A coolant table: flow [kg/s], D [m], midspan mean T [K]."""
    return {
        int(c["hole"]): {
            "m_dot": float(c["mdot_kg_s"]),
            "D": float(c["D_cm"]) / 100.0,
            "T_avg": float(c["T_avg_K"]),
        }
        for c in _rows(f"c3x_run{run}_coolant.csv")
    }


def build(run: int, multipliers: dict[str, float], T_in: dict[int, float]) -> FlowNetwork:
    """The network of one run, with the given gas-side Nu multipliers and
    hole inlet temperatures."""
    cond = run_conditions(run)
    holes = coolant(run)
    Y = _air()
    net = FlowNetwork()
    for g in segments(run):
        e = g["id"]
        net.add_node(PressureBoundary(f"{e}_in", Pt=cond["Pt1"], Tt=cond["Tt1"], Y=Y))
        net.add_node(
            PressureBoundary(f"{e}_out", Pt=g["Ps_over_Pt1"] * cond["Pt1"], Tt=cond["Tt1"], Y=Y)
        )
        net.add_element(
            ChannelElement(
                e,
                from_node=f"{e}_in",
                to_node=f"{e}_out",
                length=g["L"],
                diameter=g["Dh"],
                roughness=0.0,
                regime="compressible",
                surface=ConvectiveSurface(
                    model=SmoothModel(), Nu_multiplier=multipliers.get(e, 1.0)
                ),
            )
        )
        hole = holes[g["hole"]]
        c = f"c{g['hole']}_{e}"
        net.add_node(
            MassFlowBoundary(
                f"{c}_in", m_dot=g["f"] * hole["m_dot"], Tt=T_in.get(g["hole"], hole["T_avg"]), Y=Y
            )
        )
        net.add_node(PressureBoundary(f"{c}_out", Pt=P_COOLANT_EXIT, Tt=hole["T_avg"], Y=Y))
        net.add_element(
            ChannelElement(
                c,
                from_node=f"{c}_in",
                to_node=f"{c}_out",
                length=SPAN,
                diameter=hole["D"] * math.sqrt(g["f"]),
                Dh=hole["D"],
                roughness=0.0,
                regime="compressible",
                surface=ConvectiveSurface(
                    model=SmoothModel(), Nu_multiplier=area_ratio(g, hole["D"])
                ),
            )
        )
        net.add_wall(
            ThermalWall(
                f"w_{e}",
                element_a=e,
                element_b=c,
                layers=[WallLayer(thickness=g["thickness"], conductivity=20.0, material=MATERIAL)],
                contact_area=g["L"] * SPAN,
            )
        )
    return net


def area_ratio(g: dict[str, Any], D: float) -> float:
    """Hole-side area over the segment's: the wall has one area (the segment's)."""
    return g["f"] * math.pi * D / g["L"]


def solve(
    run: int, multipliers: dict[str, float] | None = None, T_in: dict[int, float] | None = None
) -> tuple[dict[str, float], dict[int, float]]:
    """Solve one run. The hole inlet temperatures are iterated until each
    hole's midspan temperature is the report's average (the report assumes a
    linear rise along the span)."""
    multipliers = multipliers or {}
    holes = coolant(run)
    T_in = dict(T_in or {})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for _ in range(6):
            sol = NetworkSolver(build(run, multipliers, T_in)).solve()
            if not sol.get("__success__", True):
                raise RuntimeError(f"run {run}: the network did not solve")
            change = 0.0
            for i, hole in holes.items():
                subs = [g for g in segments(run) if g["hole"] == i]
                if not subs:
                    continue
                rise = sum(
                    g["f"] * (sol[f"c{i}_{g['id']}.Tt_out"] - sol[f"c{i}_{g['id']}.Tt_in"])
                    for g in subs
                )
                new = hole["T_avg"] - 0.5 * rise
                change = max(change, abs(new - T_in.get(i, hole["T_avg"])))
                T_in[i] = new
            if change < 0.01:
                break
    return sol, T_in


def wall_temperatures(run: int, sol: dict[str, float]) -> list[dict[str, Any]]:
    """Per segment: measured and network Tw, gas-side h and the report's h."""
    return [
        {
            **g,
            "Tw_network": sol[f"w_{g['id']}.T_hot"],
            "h_network": sol[f"{g['id']}.htc"],
        }
        for g in segments(run)
    ]


def summary(rows: list[dict[str, Any]]) -> dict[str, float]:
    d = np.array([r["Tw_network"] - r["Tw_measured"] for r in rows])
    return {
        "n": len(d),
        "bias": float(d.mean()),
        "mae": float(np.abs(d).mean()),
        "max": float(np.abs(d).max()),
    }


def tune(run: int = TUNE_RUN, tol: float = 0.01) -> tuple[dict[str, float], dict[int, float]]:
    """Fit one gas-side Nu multiplier per segment to the measured Tw (secant
    per segment; the segments couple only through their shared holes)."""
    mult: dict[str, float] = {}
    history: dict[str, list[tuple[float, float]]] = {}
    T_in: dict[int, float] | None = None
    for _ in range(20):
        sol, T_in = solve(run, mult, T_in)
        worst = 0.0
        for r in wall_temperatures(run, sol):
            e = r["id"]
            err = r["Tw_network"] - r["Tw_measured"]
            history.setdefault(e, []).append((mult.get(e, 1.0), err))
            worst = max(worst, abs(err))
            h = history[e]
            if len(h) >= 2 and h[-1][1] != h[-2][1]:
                (m0, e0), (m1, e1) = h[-2], h[-1]
                new = m1 - e1 * (m1 - m0) / (e1 - e0)
            else:
                new = h[-1][0] * (1.5 if err < 0.0 else 0.7)
            mult[e] = float(min(max(new, 0.05), 20.0))
        if worst < tol:
            return {k: mult[k] for k in sorted(mult)}, T_in
    raise RuntimeError(f"tuning run {run} did not converge (worst {worst:.3f} K)")


def _print_rows(rows: list[dict[str, Any]], multipliers: dict[str, float]) -> None:
    print("  segment   s/arc       hole  t [mm]  Tw meas  Tw net    dTw   h net  h report  Nu mult")
    for r in rows:
        print(
            f"  {r['id']:9s} {r['s0']:.2f}-{r['s1']:.2f}  {r['hole']:4d}  {r['thickness'] * 1e3:5.2f}"
            f"  {r['Tw_measured']:7.1f}  {r['Tw_network']:6.1f}  {r['Tw_network'] - r['Tw_measured']:+5.1f}"
            f"  {r['h_network']:6.0f}  {r['h_report']:8.0f}  {multipliers.get(r['id'], 1.0):7.2f}"
        )


def _line(label: str, s: dict[str, float]) -> None:
    print(
        f"{label:40s} {s['n']} segments: bias {s['bias']:+5.1f} K, MAE {s['mae']:4.1f} K, "
        f"worst {s['max']:4.1f} K"
    )


# -- GUI file ------------------------------------------------------------------


def gui_document(run: int = TUNE_RUN, multipliers: dict[str, float] | None = None) -> dict:
    """The tuned network of one run as a GUI file (one row per segment)."""
    multipliers = TUNED if multipliers is None else multipliers
    _, T_in = solve(run, multipliers)
    cond = run_conditions(run)
    holes = coolant(run)
    air = {"source": "dry_air", "mode": "mole", "relative_humidity": 0}
    nodes: list[dict] = []
    edges: list[dict] = []

    def node(nid: str, kind: str, x: float, y: float, data: dict) -> None:
        nodes.append(
            {"id": nid, "type": kind, "position": {"x": x, "y": y}, "data": {"id": nid, **data}}
        )

    def flow(a: str, b: str) -> None:
        edges.append(
            {
                "id": f"e_{a}_{b}",
                "source": a,
                "sourceHandle": "flow-source",
                "target": b,
                "targetHandle": "flow-target",
                "type": "default",
                "data": {},
            }
        )

    for k, g in enumerate(segments(run)):
        e, i = g["id"], g["hole"]
        c = f"c{i}_{e}"
        hole = holes[i]
        y = 160.0 * k
        node(
            f"{e}_in",
            "pressure_boundary",
            0,
            y,
            {
                "label": f"{e} inlet (Pt1, Tt1)",
                "Pt": cond["Pt1"],
                "Tt": cond["Tt1"],
                "composition": air,
            },
        )
        node(
            e,
            "channel",
            220,
            y,
            {
                "label": f"{g['side']} s/arc {g['s0']:.2f}-{g['s1']:.2f} (hole {i}), tuned",
                "L": g["L"],
                "D": g["Dh"],
                "roughness": 0.0,
                "regime": "compressible",
                "Nu_multiplier": multipliers.get(e, 1.0),
                "f_multiplier": 1,
            },
        )
        node(
            f"{e}_out",
            "pressure_boundary",
            440,
            y,
            {
                "label": f"{e} measured Ps = {g['Ps_over_Pt1']:.3f} Pt1",
                "Pt": g["Ps_over_Pt1"] * cond["Pt1"],
                "Tt": cond["Tt1"],
                "composition": air,
            },
        )
        node(
            f"{c}_in",
            "mass_boundary",
            0,
            y + 70,
            {
                "label": f"hole {i}, share {g['f']:.2f}",
                "m_dot": g["f"] * hole["m_dot"],
                "Tt": T_in[i],
                "composition": air,
            },
        )
        node(
            c,
            "channel",
            220,
            y + 70,
            {
                "label": f"hole {i} share; Nu x {area_ratio(g, hole['D']):.3f} = hole/segment area",
                "L": SPAN,
                "D": hole["D"] * math.sqrt(g["f"]),
                "Dh": hole["D"],
                "roughness": 0.0,
                "regime": "compressible",
                "Nu_multiplier": area_ratio(g, hole["D"]),
                "f_multiplier": 1,
            },
        )
        node(
            f"{c}_out",
            "pressure_boundary",
            440,
            y + 70,
            {
                "label": f"hole {i} exit",
                "Pt": P_COOLANT_EXIT,
                "Tt": hole["T_avg"],
                "composition": air,
            },
        )
        flow(f"{e}_in", e)
        flow(e, f"{e}_out")
        flow(f"{c}_in", c)
        flow(c, f"{c}_out")
        layer = {"thickness": g["thickness"], "conductivity": 20.0, "material": MATERIAL}
        edges.append(
            {
                "id": f"w_{e}",
                "source": e,
                "sourceHandle": "thermal-source",
                "target": c,
                "targetHandle": "thermal-target",
                "type": "thermal",
                "data": {
                    "type": "thermal",
                    "thickness": g["thickness"],
                    "conductivity": 20.0,
                    "area": g["L"] * SPAN,
                    "layers": [layer],
                },
            }
        )
    rows = wall_temperatures(run, solve(run, multipliers, T_in)[0])
    s = summary(rows)
    return {
        "version": "1.2",
        "example": {
            "title": "Hylton 1983 C3X cooled vane",
            "source": (
                "Hylton, Mihelc, Turner, Nealy and York (1983), NASA CR-168015: the C3X "
                "radially cooled vane cascade, run 158 (Appendix A), Tables III, IV, IX, Fig. 7."
            ),
            "reproduced": (
                "Each row is one stretch of vane surface and the cooling hole nearest to it: hot "
                "gas from the cascade inlet total state to the measured local static pressure, "
                "a ss310 wall as thick as the mean distance to the hole, and the hole's share of "
                "its coolant flow. The gas-side Nu multipliers are fitted to the measured midspan "
                "wall temperature."
            ),
            "where_to_look": (
                "Solve, then select a wall edge: T_hot is the gas-side wall temperature, equal "
                "to the measured Tw of that segment (the tuning target). The gas channels' "
                "Nu_multiplier is the fitted tuner (0.8-1.2 mid-chord suction side, about 2.2 at "
                "the leading edge, where a channel correlation has no stagnation flow). The "
                "coolant channels' Nu_multiplier is geometry: the hole's area over the segment's."
            ),
            "label": "tuned",
            "label_note": (
                f"Tuned to run 158 ({s['n']} Nu multipliers for {s['n']} measured Tw; none left "
                "over). Untuned, the smooth-channel network misses run 158 by MAE 28 K (bias -23 K, "
                "worst -86 K at the leading edge): cross-source. Run 149 (Tu 6.5% vs 8.3%) with "
                "these multipliers: MAE 12 K, bias +11 K. The report's h is inferred by its own FE "
                "model, not measured, and is not a target."
            ),
            "script": "python/examples/papers/hylton1983_c3x_vane.py",
        },
        "nodes": nodes,
        "edges": edges,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--tune", action="store_true", help="refit the multipliers to run 158")
    parser.add_argument("--write-gui", action="store_true", help=f"rewrite {GUI_FILE.name}")
    parser.add_argument("--no-plot", action="store_true", help="accepted; there is no plot")
    args = parser.parse_args()

    if not DATA.is_dir():
        sys.exit(f"Transcribed data not found at {DATA}: run from a repo checkout.")

    mult = TUNED
    if args.tune:
        mult, _ = tune()
        print("TUNED = {")
        for k, v in mult.items():
            print(f'    "{k}": {v:.4f},')
        print("}")

    print("Hylton et al. (1983), NASA CR-168015, C3X vane -- TUNED (untuned: cross-source)")
    for run in (158, 149):
        c = run_conditions(run)
        print(
            f"Run {run}: Pt1 {c['Pt1'] / 1e3:.1f} kPa, Tt1 {c['Tt1']:.0f} K, Tu {c['Tu']}%, M2 {c['M2']}"
        )

    results = {}
    for run, m, label in (
        (158, {}, "run 158, untuned"),
        (158, mult, "run 158, tuned to run 158"),
        (149, {}, "run 149, untuned"),
        (149, mult, "run 149, run 158's multipliers"),
    ):
        rows = wall_temperatures(run, solve(run, m)[0])
        results[label] = rows
        print(f"\n{label}")
        _print_rows(rows, m)
    print()
    for label, rows in results.items():
        _line(label, summary(rows))

    if args.write_gui:
        GUI_FILE.write_text(json.dumps(gui_document(), indent=2) + "\n")
        print(f"\nwrote {GUI_FILE.relative_to(REPO)}")


if __name__ == "__main__":
    main()
