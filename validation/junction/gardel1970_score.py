"""Score junction closures against Gardel & Rechsteiner 1970 (cross-source).

Role (docs/VALIDATION_POLICY.md): CROSS-SOURCE ACCURACY for every model here.
None of them was built or tuned on this campaign: the ideal tee is
parameter-free, Bassett 2001 is his own analysis, and the shipped Mynard tee
was tuned on Bassett and Idelchik. The source states no uncertainty, so no
band is applied and nothing here is "within" anything.

Data: ``data/gardel1970/detail_points.csv``, the raw points of Figs 12-15
(water, so every model is evaluated at constant density). Gardel's four flow
modes on his branches alpha (upstream main), gamma (downstream main) and beta
(branch, at delta to alpha):

    mode 1  alpha -> beta + gamma   common alpha   Bassett type 3
    mode 2  alpha + beta -> gamma   common gamma   Bassett type 6
    mode 3  beta -> alpha + gamma   common beta    Bassett type 2
    mode 4  alpha + gamma -> beta   common beta    Bassett type 5

Each mode is scored on the two series that pair the common port with another
(the third follows from h_ya = h_yb + h_ba). A point where that other port is
dead is a pressure tap, not a loss, and is reported apart.

Every model is evaluated on Gardel's own flow pattern with outward port axes
alpha = pi, gamma = 0, beta = pi - delta and areas A_main = 1, A_beta = phi:
  - ideal: the production C++ kernel (``_core.ideal_junction_residuals_and_jacobian``)
    at f = 1;
  - Bassett: his published coefficient for the flow type, angle-corrected as
    Section 4.2 recommends for measured data, each on its own leg's fraction.
    Splits: K = 1 + r^2 - 2 r cos(3/4 d) for a turning receiver (K1/K3/K4/K6),
    q^2 - 1.5 q + 0.5 for the collinear one (K2/K5). Merges: K11/K12 (type 6,
    jet with the outlet flow), K7/K8 (type 4, jet against it, his angle
    mirrored), K9/K10 (type 5, Eqs 35/37). At delta = 90 a mode-2 merge is
    both type 4 and type 6, and his two corrected forms disagree there; type 6
    is used, as for his own 90 deg K11/K12 curves.
  - shipped: Mynard with the production joining alpha and recovery.

Errors are model loss minus measured loss, in Gardel's normalisation
h = dH / (V^2/2g) with V the total flow over the MAIN pipe area; positive =
the model predicts more loss.

    uv run python -m validation.junction.gardel1970_score
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from combaero import _core
from combaero.network._mynard2010 import junction_loss_coefficient
from validation.junction.models import bassett2001 as bassett

DATA = Path(__file__).parent / "data" / "gardel1970" / "detail_points.csv"

ALPHA, GAMMA, BETA = 0, 1, 2
NAMES = {ALPHA: "a", GAMMA: "y", BETA: "b"}

#: mode -> (common port, port flows out of the junction as a function of q,
#: where q is the figure's abscissa and the total flow is 1).
MODES = {
    "1": (ALPHA, lambda q: (-1.0, 1.0 - q, q)),
    "2": (GAMMA, lambda q: (-(1.0 + q), 1.0, q)),
    "3": (BETA, lambda q: (1.0 - q, q, -1.0)),
    "4": (BETA, lambda q: (-(1.0 + q), q, 1.0)),
}

#: series name -> (i, j) with h = H_i - H_j.
SERIES = {"h_ba": (BETA, ALPHA), "h_ya": (GAMMA, ALPHA), "h_yb": (GAMMA, BETA)}

DEAD = 1e-9

#: The shipped tee's joining energy-transfer alpha (joining_variants.MODELS['shipped']).
SHIPPED_JOINING_ALPHA = 0.2


def geometry(delta_deg: float, phi: float) -> tuple[np.ndarray, np.ndarray]:
    axis = np.array([math.pi, 0.0, math.pi - math.radians(delta_deg)])
    area = np.array([1.0, 1.0, phi])
    return axis, area


def scored_series(mode: str) -> dict[str, int]:
    """The two series pairing the common port with another: name -> other port."""
    c = MODES[mode][0]
    out = {}
    for name, (i, j) in SERIES.items():
        if c in (i, j):
            out[name] = j if i == c else i
    return out


def loss_from_h(h: float, name: str, m_out: tuple[float, ...]) -> float:
    """Measured loss (Pt drop along the flow, Gardel units) of a scored series."""
    i, j = SERIES[name]
    # H_i - H_j is a loss when the flow runs i -> j through the junction.
    return h if m_out[i] < 0.0 else -h


def ideal_losses(delta: float, phi: float, m_out: tuple[float, ...]) -> np.ndarray:
    """Total pressure per port from the production kernel, rho = 1, Gardel units."""
    axis, area = geometry(delta, phi)
    geom = _core.IdealJunctionGeometry()
    geom.area = [float(a) for a in area]
    geom.axis_rad = [float(a) for a in axis]
    geom.port_sign = [1.0, 1.0, 1.0]
    r = _core.ideal_junction_residuals_and_jacobian([0.0] * 3, [1.0] * 3, [0.0] * 3, list(m_out), 0.0, geom)
    p = -np.array(r.residual[:3])  # rows are P_i - p_jn + g_i(m), linear in P_i
    u = np.array(m_out) / area
    return 2.0 * (p + 0.5 * u * u)  # Pt / (rho V^2 / 2), V = 1


def bassett_split(delta: float, phi: float, m_out: tuple[float, ...], c: int) -> dict[int, float]:
    axis, area = geometry(delta, phi)
    u_c = 1.0 / area[c]
    out = {}
    for j in range(3):
        if j == c:
            continue
        d = abs(math.remainder(axis[j] - (axis[c] + math.pi), 2.0 * math.pi))
        if d < 1e-9:  # collinear receiver: Bassett's K2/K5
            q = m_out[j]
            out[j] = q * q - 1.5 * q + 0.5
        else:
            r = (m_out[j] / area[j]) / u_c
            out[j] = 1.0 + r * r - 2.0 * r * math.cos(0.75 * d)
    return out


def bassett_merge(delta: float, phi: float, mode: str, q: float) -> dict[int, float]:
    """Bassett's published corrected joining forms (Section 4.2), each on its own leg's fraction."""
    psi = 1.0 / phi
    if mode == "2":
        q_lat, q_str = -q, 1.0 + q
        d = math.radians(delta)  # the branch jet's deflection into gamma
        if d <= math.pi / 2:  # type 6: jet with the outlet flow
            return {BETA: bassett.K12_corr(q_lat, psi, d), ALPHA: bassett.K11_corr(q_str, psi, d)}
        th = math.pi - d  # type 4: jet against it, his angle mirrored
        return {BETA: bassett.K7_corr(q_lat, psi, th), ALPHA: bassett.K8_corr(q_str, psi, th)}
    if mode == "4":  # type 5: A = alpha, C = gamma, theta = angle of beta to alpha
        th = math.radians(delta)
        return {ALPHA: bassett.K9_corr(1.0 + q, psi, th), GAMMA: bassett.K10_corr(-q, psi, th)}
    raise ValueError(mode)


def model_losses(model: str, delta: float, phi: float, mode: str, q: float) -> dict[int, float]:
    """Loss common <-> other port, in Gardel units, keyed by the other port."""
    c, flows = MODES[mode]
    m_out = flows(q)
    axis, area = geometry(delta, phi)
    scale = (1.0 / area[c]) ** 2  # K is on the common port's head, h on V = 1
    if model == "ideal":
        pt = ideal_losses(delta, phi, m_out)
        return {j: (pt[c] - pt[j]) if m_out[c] < 0 else (pt[j] - pt[c]) for j in range(3) if j != c}
    U = -np.array(m_out) / area  # Mynard: positive INTO the junction
    if model == "bassett":
        K = bassett_split(delta, phi, m_out, c) if m_out[c] < 0.0 else bassett_merge(delta, phi, mode, q)
    elif model == "shipped":
        # K per non-common port in index order: collectors of a split, suppliers of a merge
        k = junction_loss_coefficient(U, area, axis, SHIPPED_JOINING_ALPHA, 0.0).K
        K = dict(zip([j for j in range(3) if j != c], (float(x) for x in np.atleast_1d(k)), strict=True))
    else:
        raise ValueError(model)
    return {j: k * scale for j, k in K.items()}


def load() -> list[dict]:
    with DATA.open() as fh:
        return list(csv.DictReader(fh))


MODELS = ("ideal", "bassett", "shipped")

Cell = tuple[str, str, str, str]  # fig, delta/phi/rho, mode, series


def score() -> tuple[dict[Cell, np.ndarray], int]:
    """Per cell, rows of (measured loss, error of each model in MODELS); dead-port taps left out."""
    cells: dict[Cell, list[tuple[float, ...]]] = defaultdict(list)
    taps = 0
    for r in load():
        mode, name = r["mode"], r["series"]
        other = scored_series(mode).get(name)
        if other is None:
            continue
        q, h = float(r["q"]), float(r["h"])
        delta, phi = float(r["delta_deg"]), float(r["phi"])
        m_out = MODES[mode][1](q)
        if abs(m_out[other]) < DEAD:
            taps += 1
            continue
        meas = loss_from_h(h, name, m_out)
        errs = [model_losses(m, delta, phi, mode, q)[other] - meas for m in MODELS]
        cells[(r["fig"], f"{delta:g}/{phi:g}/{r['rho']}", mode, name)].append((meas, *errs))
    return {k: np.array(v) for k, v in cells.items()}, taps


def pooled(cells: dict[Cell, np.ndarray], select=lambda k: True) -> dict[str, tuple[int, float, float]]:
    """Model -> (n, MAE, bias) over the selected cells."""
    out = {}
    for i, m in enumerate(MODELS):
        e = np.concatenate([a[:, i + 1] for k, a in cells.items() if select(k)])
        out[m] = (int(e.size), float(np.abs(e).mean()), float(e.mean()))
    return out


def _line(label: str, stats: dict[str, tuple[int, float, float]]) -> str:
    n = stats[MODELS[0]][0]
    return f"{label:31s} {n:4d}" + "".join(f"  {stats[m][1]:7.3f} ({stats[m][2]:+.3f})" for m in MODELS)


def main() -> None:
    cells, taps = score()
    print("Gardel & Rechsteiner 1970, Figs 12-15 -- role: CROSS-SOURCE ACCURACY (no stated band)")
    print("error = model loss - measured loss, Gardel units (V = total flow / main area); MAE (bias)\n")
    print(f"{'fig  delta/phi/rho  mode  series':31s} {'n':>4s}" + "".join(f"  {m:>16s}" for m in MODELS))
    for k in sorted(cells, key=lambda k: (int(k[0]), k[2], k[3])):
        print(_line(f"{k[0]:>3}  {k[1]:>13}  {k[2]:>4}  {k[3]:>6}", pooled(cells, lambda c, k=k: c == k)))
    print()
    for mode in "1234":
        print(_line(f"all, mode {mode}", pooled(cells, lambda c, mode=mode: c[2] == mode)))
    for fig in ("12", "13", "14", "15"):
        print(_line(f"all, fig {fig}", pooled(cells, lambda c, fig=fig: c[0] == fig)))
    print(_line("all", pooled(cells)))
    print(f"\ndead-port pressure taps left out: {taps}")


if __name__ == "__main__":
    main()
