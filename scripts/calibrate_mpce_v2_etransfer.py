"""
Calibrate joining_etransfer_alpha against Bassett + Idelchik analytical anchors.

Methodology (per agreed approach in PR #192):
  - Anchor: Bassett K11_corr/K12_corr + Idelchik tabulated values (analytical).
  - Validate: Bassett measured K11/K12 + Bassett Fig 12a K8 (held out).

The calibration solves for a single scalar alpha that minimizes the
sum-of-squared residuals between Mynard+correction and the analytical
anchors across a (theta, psi, q) grid covering psi>=1.25 (where the
correction matters).
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from validation.junction.models import bassett2001
from validation.junction.models.mynard2010 import junction_loss_coefficient

DATA_DIR = Path("/Users/thiemo/Projects/combaero/validation/junction/data/idelchik1966")


def mynard_K(q: float, psi: float, theta_deg: float, alpha: float) -> tuple[float, float]:
    """Compute Mynard (K_str, K_bra) for joining type-6 with optional alpha correction.

    Reference scales: rho=1, A_com=1, m_com=1, so K values are normalized by
    0.5*rho*u_com^2 = 0.5. Sign convention: U > 0 = supplier (into junction).
    """
    A_str = 1.0
    A_bra = 1.0 / psi
    A_com = 1.0
    m_str = (1.0 - q) * 1.0
    m_bra = q * 1.0
    m_com = 1.0
    U = np.array([m_str / A_str, m_bra / A_bra, -m_com / A_com])
    A = np.array([A_str, A_bra, A_com])
    theta = np.array([0.0, math.radians(theta_deg), math.pi])
    r = junction_loss_coefficient(U, A, theta, joining_etransfer_alpha=alpha)
    if r.K is None or len(r.K) != 2:
        return float("nan"), float("nan")
    return float(r.K[0]), float(r.K[1])  # K_str, K_bra


def bassett_K(q: float, psi: float, theta_deg: float) -> tuple[float, float]:
    """Bassett K11 (str), K12 (bra) at the joining type-6 cell."""
    theta = math.radians(theta_deg)
    return (
        bassett2001.K11_corr(q, psi, theta),
        bassett2001.K12_corr(q, psi, theta),
    )


def load_idelchik(K_id: str, theta_deg: int, psi: float) -> dict[float, float] | None:
    """Load Idelchik tabulated K (returns q->K mapping) or None if not on grid."""
    # Match the stored filename convention
    psi_label = (
        f"{int(round(psi))}"
        if abs(psi - round(psi)) < 1e-6
        else f"{psi:.3f}".rstrip("0").rstrip(".")
    )
    # Pick the diagram per (K_id, theta)
    diag = {
        ("K12", 30): "7-1",
        ("K11", 30): "7-2",
        ("K12", 45): "7-3",
        ("K11", 45): "7-4",
        ("K12", 90): "7-7",
        # K11 at theta=90 is psi-independent (7-7b); handled separately
    }.get((K_id, theta_deg))
    if diag is None:
        return None
    path = DATA_DIR / f"idelchik_diag{diag}_{K_id}_theta={theta_deg}_psi={psi_label}_tabulated.csv"
    if not path.exists():
        return None
    out = {}
    with open(path) as fh:
        reader = csv.reader(fh)
        next(reader)  # header
        for row in reader:
            out[float(row[0])] = float(row[1])
    return out


def load_idelchik_K11_theta90() -> dict[float, float]:
    """Special-case: psi-independent K11 at theta=90."""
    path = DATA_DIR / "idelchik_diag7-7b_K11_theta=90_psi-independent_tabulated.csv"
    out = {}
    with open(path) as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            out[float(row[0])] = float(row[1])
    return out


# Calibration grid: psi >= 1.25 (correction matters), q in [0.1, 0.9] (avoid
# degenerate endpoints), three angles spanning the range
PSI_GRID = [1.25, 1.667, 2.5, 3.333]  # practical range; extreme psi excluded
THETA_GRID = [30, 45, 90]
Q_GRID = [round(0.1 * i, 1) for i in range(1, 10)]  # 0.1 to 0.9

# Known Idelchik typos -- exclude from loss
EXCLUDED_CELLS = {
    # (theta, psi, q): reason
    (30, 10.0, 0.9): "Idelchik 7-1 typo (printed 58.0 vs formula 67.9)",
    (90, 10.0, 0.7): "Idelchik 7-7 typo (printed 42.9 vs formula 49.8)",
}


def collect_anchors() -> list[tuple[int, float, float, str, str, float]]:
    """Build the list of (theta, psi, q, K_id, source, K_value) anchors."""
    anchors = []
    idel_K11_th90 = load_idelchik_K11_theta90()
    for theta in THETA_GRID:
        for psi in PSI_GRID:
            for q in Q_GRID:
                if (theta, psi, q) in EXCLUDED_CELLS:
                    continue
                # Bassett analytical (always available)
                K11_bas, K12_bas = bassett_K(q, psi, theta)
                anchors.append((theta, psi, q, "K11", "bassett", K11_bas))
                anchors.append((theta, psi, q, "K12", "bassett", K12_bas))
                # Idelchik tabulated (if available at this cell)
                if theta == 90:
                    # K11 at theta=90 is psi-independent
                    K11_idel = idel_K11_th90.get(q)
                    if K11_idel is not None:
                        anchors.append((theta, psi, q, "K11", "idelchik", K11_idel))
                    idel_K12 = load_idelchik("K12", theta, psi)
                    if idel_K12 is not None and q in idel_K12:
                        anchors.append((theta, psi, q, "K12", "idelchik", idel_K12[q]))
                else:
                    idel_K11 = load_idelchik("K11", theta, psi)
                    idel_K12 = load_idelchik("K12", theta, psi)
                    if idel_K11 is not None and q in idel_K11:
                        anchors.append((theta, psi, q, "K11", "idelchik", idel_K11[q]))
                    if idel_K12 is not None and q in idel_K12:
                        anchors.append((theta, psi, q, "K12", "idelchik", idel_K12[q]))
    return anchors


def residuals(alpha: float, anchors: list) -> np.ndarray:
    """Per-anchor signed residual = K_mynard_corrected - K_anchor."""
    res = []
    for theta, psi, q, K_id, _src, K_anchor in anchors:
        K_str, K_bra = mynard_K(q, psi, theta, alpha)
        K_my = K_str if K_id == "K11" else K_bra
        res.append(K_my - K_anchor)
    return np.array(res)


def relative_residuals(alpha: float, anchors: list, eps: float = 0.5) -> np.ndarray:
    """Per-anchor relative residual = (K_my - K_anchor) / max(|K_anchor|, eps).

    Using eps=0.5 (comparable to typical K scale) prevents tiny-K cells
    (near sign-flips of K11) from blowing up the relative metric while
    still letting large-K cells contribute proportionally.
    """
    r = residuals(alpha, anchors)
    K_anchor = np.array([a[5] for a in anchors])
    denom = np.maximum(np.abs(K_anchor), eps)
    return r / denom


def loss_relative(alpha: float, anchors: list) -> float:
    r = relative_residuals(alpha, anchors)
    return float(np.sum(r * r))


def loss_absolute(alpha: float, anchors: list) -> float:
    r = residuals(alpha, anchors)
    return float(np.sum(r * r))


def main():
    anchors = collect_anchors()
    print(f"Calibration anchors: {len(anchors)}")
    print(
        f"  By source: bassett={sum(1 for a in anchors if a[4] == 'bassett')}, "
        f"idelchik={sum(1 for a in anchors if a[4] == 'idelchik')}"
    )

    from collections import defaultdict

    # First: signed-residual breakdown at alpha=0 (baseline) -- separated by K_id
    print("\n=== Signed residual at alpha=0, by K_id ===")
    print(f"{'theta':>5s} {'psi':>6s} {'K':>4s}  {'N':>3s}  {'signed mean':>12s}")
    r0 = residuals(0.0, anchors)
    by_baseline_K = defaultdict(list)
    for i, (theta, psi, _q, K_id, _src, _K_anchor) in enumerate(anchors):
        by_baseline_K[(theta, psi, K_id)].append(float(r0[i]))
    for k in sorted(by_baseline_K):
        rr = by_baseline_K[k]
        print(f"{k[0]:>5d} {k[1]:>6.2f} {k[2]:>4s}  {len(rr):>3d}  {float(np.mean(rr)):>+12.4f}")

    def summary(alpha: float, name: str) -> None:
        r = residuals(alpha, anchors)
        rrel = relative_residuals(alpha, anchors)
        print(
            f"  {name:>10s} alpha={alpha:+.4f}  "
            f"mean|r|={float(np.mean(np.abs(r))):.4f}  "
            f"mean|rel|={float(np.mean(np.abs(rrel))):.4f}  "
            f"max|r|={float(np.max(np.abs(r))):.4f}"
        )

    # Coarse sweep on both loss criteria
    print("\n=== Coarse sweep ===")
    print(f"{'alpha':>7s}  {'mean|r|':>10s}  {'mean|rel|':>10s}  {'max|r|':>10s}")
    for alpha in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]:
        r = residuals(alpha, anchors)
        rrel = relative_residuals(alpha, anchors)
        mae = float(np.mean(np.abs(r)))
        mrel = float(np.mean(np.abs(rrel)))
        mxe = float(np.max(np.abs(r)))
        print(f"{alpha:>7.2f}  {mae:>10.4f}  {mrel:>10.4f}  {mxe:>10.4f}")

    # Two optimizers
    print("\n=== Optimization (two loss criteria) ===")
    res_abs = minimize_scalar(loss_absolute, bounds=(-1.0, 3.0), args=(anchors,), method="bounded")
    res_rel = minimize_scalar(loss_relative, bounds=(-1.0, 3.0), args=(anchors,), method="bounded")
    alpha_abs = float(res_abs.x)
    alpha_rel = float(res_rel.x)
    summary(0.0, "baseline")
    summary(alpha_abs, "abs-opt")
    summary(alpha_rel, "rel-opt")

    # Stratified relative-residual breakdown at the relative-loss optimum
    print(f"\n=== Stratified residuals at alpha_rel = {alpha_rel:+.4f} ===")
    print(
        f"{'theta':>5s} {'psi':>6s}  {'N':>3s}  {'mean|r|':>8s}  {'mean|rel|':>10s}  {'signed mean':>12s}"
    )
    by_strat: dict[tuple[int, float], list[tuple[float, float]]] = defaultdict(list)
    r_rel_opt = residuals(alpha_rel, anchors)
    rrel_opt = relative_residuals(alpha_rel, anchors)
    for i, (theta, psi, _q, _K_id, _src, _K_anchor) in enumerate(anchors):
        by_strat[(theta, psi)].append((float(r_rel_opt[i]), float(rrel_opt[i])))
    for k in sorted(by_strat):
        rr_abs = [a for a, _ in by_strat[k]]
        rr_rel = [b for _, b in by_strat[k]]
        mae = float(np.mean(np.abs(rr_abs)))
        mrel = float(np.mean(np.abs(rr_rel)))
        signed = float(np.mean(rr_abs))
        print(
            f"{k[0]:>5d} {k[1]:>6.2f}  {len(rr_abs):>3d}  "
            f"{mae:>8.4f}  {mrel:>10.4f}  {signed:>+12.4f}"
        )


if __name__ == "__main__":
    main()
