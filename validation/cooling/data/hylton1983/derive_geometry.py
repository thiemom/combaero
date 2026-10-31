"""C3X passage geometry from Hylton et al. 1983 (NASA CR-168015) Table III/IV + Fig. 7.

Run from the repository root: ``python validation/cooling/data/hylton1983/derive_geometry.py``.
Reads the transcribed tables beside it and rewrites the derived ones:
  c3x_cooling_holes_xy.csv     hole centres mapped from Fig. 7 (U, V) to Table III (x, y)
  c3x_wall_thickness.csv       hole perimeter -> nearest outer surface (overall, SS side, PS side)
  c3x_passage_geometry.csv     blade-to-blade width w, area w*H, Dh = 2wH/(w+H) along the passage
  c3x_station_check.csv        tabulated (s/arc -> x/c_ax) vs our arc-length parametrisation
"""

import csv
from pathlib import Path

import numpy as np

D = str(Path(__file__).resolve().parent) + "/"
PITCH = 11.773  # cm, Table IV
H = 7.620  # cm, vane height, Table IV
CAX = 7.816  # cm, axial chord
ARC_SS = 17.782
ARC_PS = 13.723
SRC = "# Source: Hylton et al. 1983, NASA CR-168015, Table III (coords, PDF p.19), Table IV (PDF p.21), Fig. 7 (holes, PDF p.23)."


def load(fn: str) -> list[dict]:
    with open(D + fn) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


P = np.array([[float(r["x_cm"]), float(r["y_cm"])] for r in load("c3x_coordinates.csv")])
HOLES = load("c3x_cooling_holes.csv")
UV = np.array([[float(r["U_cm"]), float(r["V_cm"])] for r in HOLES])
DH = np.array([float(r["D_cm"]) for r in HOLES])
N = len(P)


def closest_on_poly(
    p: np.ndarray, poly: np.ndarray, closed: bool = True
) -> tuple[float, np.ndarray, int]:
    a = poly if closed else poly[:-1]
    b = np.roll(poly, -1, axis=0) if closed else poly[1:]
    ab = b - a
    t = np.clip(((p - a) * ab).sum(1) / (ab * ab).sum(1), 0, 1)
    q = a + t[:, None] * ab
    dist = np.hypot(*(q - p).T)
    k = int(np.argmin(dist))
    return float(dist[k]), q[k], k


# --- Fig. 7 (U, V) -> (x, y): V along the pressure-side bitangent (points 39-71), U normal into the vane
a, b = P[70], P[38]  # point 71 (LE side) and point 39 (TE side)
v_hat = (b - a) / np.hypot(*(b - a))  # LE -> TE along the bitangent
u_hat = np.array([v_hat[1], -v_hat[0]])
if ((P - a) @ u_hat).mean() < 0:
    u_hat = -u_hat
V_le = (P @ v_hat).min()  # V = 0 at the LE tangent normal to the bitangent
U0 = a @ u_hat
XY = np.array([v_hat * (vv + V_le) + u_hat * (u + U0) for u, vv in UV])
setting_angle = np.degrees(np.arctan2(-v_hat[1], v_hat[0]))
V_extent = (P @ v_hat).max() - V_le

# --- arc-length parametrisation; origin at the minimum-x (LE) point on a dense polyline
# contour order: 1..29 SS (LE->TE), 30..39 TE circle, 40..78 PS (TE->LE)
dense = []
for k in range(N):
    p0, p1 = P[k], P[(k + 1) % N]
    for t in np.linspace(0, 1, 50, endpoint=False):
        dense.append(p0 + t * (p1 - p0))
dense = np.array(dense)
i_le = int(np.argmin(dense[:, 0]))
dense = np.roll(dense, -i_le, axis=0)  # starts at LE origin, runs PS?  check direction below
seg = np.hypot(*np.diff(np.vstack([dense, dense[:1]]), axis=0).T)
s_cum = np.concatenate([[0], np.cumsum(seg)])[:-1]
perim = seg.sum()
# the LE origin lies between points 75 and 76; increasing index -> 76, 77, 78, 1, 2 ... = suction side
i_te_ss = int(np.argmin(np.hypot(*(dense - P[29]).T)))  # point 30, x max
# TE tangency points: SS end ~ point 29/30, PS start ~ point 39/40; use the TE circle tangency guess:
i_29 = int(np.argmin(np.hypot(*(dense - P[28]).T)))
i_40 = int(np.argmin(np.hypot(*(dense - P[39]).T)))
arc_ss_to30 = s_cum[i_te_ss]
arc_ss_to29 = s_cum[i_29]
arc_ps_to40 = perim - s_cum[i_40]


def ss_point(frac: float, arc: float) -> np.ndarray:
    s = frac * arc
    k = int(np.searchsorted(s_cum, s))
    return dense[min(k, len(dense) - 1)]


def ps_point(frac: float, arc: float) -> np.ndarray:
    s = perim - frac * arc
    k = int(np.searchsorted(s_cum, s))
    return dense[min(k, len(dense) - 1)]


# --- check the report's s/arc vs x/c_ax pairs against our parametrisation
chk = []
for fn in ("c3x_run149_Tw_h.csv", "c3x_run158_Tw_h.csv", "c3x_run149_Ps.csv"):
    for r in load(fn):
        if r["surface"] == "LE":
            continue
        f = float(r["s_over_arc"])
        x_tab = float(r["x_over_cax"])
        if r["surface"] == "SS":
            x_ours = ss_point(f, ARC_SS)[0] / CAX
        else:
            x_ours = ps_point(f, ARC_PS)[0] / CAX
        chk.append((fn, r["surface"], f, x_tab, x_ours, x_ours - x_tab))
with open(D + "c3x_station_check.csv", "w") as f:
    f.write(
        SRC
        + "\n# x_over_cax from the Appendix A tables vs x/c_ax of our point at the same s/arc (s measured from the min-x LE point, arcs from Table IV).\n"
    )
    f.write("file,surface,s_over_arc,x_over_cax_table,x_over_cax_ours,diff\n")
    for c in chk:
        f.write(f"{c[0]},{c[1]},{c[2]:.4f},{c[3]:.4f},{c[4]:.4f},{c[5]:+.4f}\n")
diffs = np.array([c[5] for c in chk])

# --- wall thickness: hole perimeter to nearest outer surface, also per side
ss_poly = P[0:30]  # points 1..30
ps_poly = np.vstack([P[39:], P[:1]])  # points 40..78, 1
le_split = None
with open(D + "c3x_cooling_holes_xy.csv", "w") as f:
    f.write(
        SRC
        + "\n# Hole centres mapped from Fig. 7 (U,V) to Table III (x,y): V along the pressure-side bitangent through points 71 and 39, V=0 at the LE tangent normal to it, U normal into the vane.\n"
    )
    f.write(
        f"# Bitangent angle from +x: {np.degrees(np.arctan2(v_hat[1], v_hat[0])):.2f} deg; profile extent along it {V_extent:.3f} cm (Table IV true chord 14.493).\n"
    )
    f.write("hole,x_cm,y_cm,D_cm\n")
    for k, (c, d) in enumerate(zip(XY, DH, strict=False)):
        f.write(f"{k + 1},{c[0]:.4f},{c[1]:.4f},{d:.3f}\n")
wt_rows = []
for k, (c, d) in enumerate(zip(XY, DH, strict=False)):
    t_all, q_all, _ = closest_on_poly(c, P)
    t_ss, q_ss, _ = closest_on_poly(c, ss_poly, closed=False)
    t_ps, q_ps, _ = closest_on_poly(c, ps_poly, closed=False)
    wt_rows.append((k + 1, d, t_all - d / 2, t_ss - d / 2, q_ss, t_ps - d / 2, q_ps))
with open(D + "c3x_wall_thickness.csv", "w") as f:
    f.write(
        SRC
        + "\n# Wall thickness = distance from hole centre to the outer contour (linear segments between Table III points) minus D/2. cm.\n"
    )
    f.write(
        "# SS side = points 1-30, PS side = points 40-78 + 1 (TE circle 31-39 excluded from the sides, included in t_min).\n"
    )
    f.write("hole,D_cm,t_min_cm,t_ss_cm,x_ss_foot,y_ss_foot,t_ps_cm,x_ps_foot,y_ps_foot\n")
    for r in wt_rows:
        f.write(
            f"{r[0]},{r[1]:.3f},{r[2]:.3f},{r[3]:.3f},{r[4][0]:.3f},{r[4][1]:.3f},{r[5]:.3f},{r[6][0]:.3f},{r[6][1]:.3f}\n"
        )

# --- passage: suction side of vane 0 vs neighbour at +pitch (its pressure side faces our suction side)
NB = P + [0.0, PITCH]
nb_ps = ps_poly + [0.0, PITCH]
rows = []
for frac in np.round(np.arange(0.0, 1.0001, 0.05), 3):
    p = ss_point(frac, ARC_SS)
    w_all, q_all, k_all = closest_on_poly(p, NB)
    w_ps, q_ps, _ = closest_on_poly(p, nb_ps, closed=False)
    # nearest point classified: on neighbour PS (covered) or at its TE circle (uncovered, past the throat)
    kind = "PS" if 39 <= k_all <= 77 or k_all == 0 else ("TE" if 29 <= k_all <= 38 else "SS")
    area = w_all * H
    dh = 2 * w_all * H / (w_all + H)
    kk = int(np.argmin(np.hypot(*(dense + [0.0, PITCH] - q_all).T)))
    nb_ps_frac = (perim - s_cum[kk]) / ARC_PS if kind == "PS" else float("nan")
    rows.append((frac, frac * ARC_SS, p, p[0] / CAX, w_all, q_all, kind, area, dh, nb_ps_frac))
throat = min(closest_on_poly(q, P)[0] for q in NB)
with open(D + "c3x_passage_geometry.csv", "w") as f:
    f.write(SRC + "\n")
    f.write(
        "# Passage between the suction surface of vane 0 and the pressure surface of the neighbour at y + pitch (pitch 11.773 cm, x axial).\n"
    )
    f.write(
        "# Stations every 5% of the SS arc (17.782 cm) from the LE origin (min x). w = shortest distance from the SS station to the neighbour contour (not a true orthogonal-trajectory width;\n"
    )
    f.write(
        "# near the LE the passage is open to the inlet, and past the throat the nearest neighbour point is its TE (nb_kind=TE): the passage is then uncovered, w is a crude bound).\n"
    )
    f.write(
        "# nb_s_over_arc_ps = PS s/arc (fraction of 13.723 cm) of the neighbour foot point, pairing SS and PS stations. Past the throat a periodic exit width is pitch*cos(72.38 deg) = 3.566 cm.\n"
    )
    f.write(
        f"# H = {H} cm vane height. A = w*H, Dh = 2wH/(w+H). Computed throat (min vane-to-neighbour distance) = {throat:.3f} cm vs Table IV 3.292 cm.\n"
    )
    f.write(
        "s_over_arc_ss,s_cm,x_cm,y_cm,x_over_cax,w_cm,nb_x_cm,nb_y_cm,nb_kind,nb_s_over_arc_ps,A_cm2,Dh_cm\n"
    )
    for r in rows:
        f.write(
            f"{r[0]:.3f},{r[1]:.3f},{r[2][0]:.4f},{r[2][1]:.4f},{r[3]:.4f},{r[4]:.4f},{r[5][0]:.4f},{r[5][1]:.4f},{r[6]},{r[9]:.3f},{r[7]:.3f},{r[8]:.4f}\n"
        )

print(
    f"perimeter {perim:.3f}; SS arc LE->pt30 {arc_ss_to30:.3f}, LE->pt29 {arc_ss_to29:.3f} (Table IV {ARC_SS}); PS arc LE->pt40 {arc_ps_to40:.3f} (Table IV {ARC_PS})"
)
print(
    f"bitangent extent {V_extent:.3f}; angle of bitangent from x {np.degrees(np.arctan2(v_hat[1], v_hat[0])):.2f}"
)
print(f"station check: n={len(diffs)} mean {diffs.mean():+.4f} max|d| {np.abs(diffs).max():.4f}")
print("throat", throat)
for r in wt_rows:
    print("hole", r[0], f"t_min {r[2]:.3f} t_ss {r[3]:.3f} t_ps {r[5]:.3f}")
for r in rows:
    print(f"{r[0]:.2f} x/c {r[3]:.3f} w {r[4]:.3f} {r[6]} Dh {r[8]:.3f}")
