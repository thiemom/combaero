"""Calibrate ``joining_etransfer_alpha`` against analytical and tabulated anchors.

This is the script that produces the value of
``MPCEv2Element.DEFAULT_JOINING_ETRANSFER_ALPHA``. It lives in the validation
tree so the tune is reproducible from a checkout; the earlier copy lived in a
gitignored scratch directory, which meant a production constant had no
in-repo provenance.

Method (the calibration methodology this repo uses: anchor on analytical,
validate on measured):

  - **Anchor** on Bassett 2001's angle-corrected joining coefficients
    (``K11_corr``, ``K12_corr``) and on Idelchik 1966's tabulated joining
    coefficients, over a (theta, psi, q) grid where the correction matters
    (psi >= 1.25).
  - **Validate** in-network against the digitised measured curves through
    ``validation.junction.network_runner`` (Bassett K11/K12 measured, Idelchik
    tables). That validation is NOT performed here -- an earlier docstring
    claimed a held-out validation that was never implemented. The scorecard is
    the validation.

q convention -- the trap this script fell into once
---------------------------------------------------
Throughout this file ``q`` is the **lateral inlet fraction**, ``m_bra / m_com``,
because that is what the Mynard call is built from. Each anchor source indexes
its coefficients on its own leg and must be converted:

  - Bassett Table 1: K12 on ``m_B/m_C`` (lateral) -- pass ``q``;
    **K11 on ``m_A/m_C`` (straight) -- pass ``1 - q``**.
  - Idelchik: every diagram on ``Q_b/Q_c`` (lateral) -- pass ``q`` for both.

The first version of this script passed ``q`` to ``K11_corr`` directly, fitting
alpha against a mirrored Bassett K11 curve. The 0.2 default was produced by
that version (issue #271).

Run:  uv run python validation/junction/calibrate_etransfer.py
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from combaero.network._mynard2010 import junction_loss_coefficient
from validation.junction.models import bassett2001

DATA_DIR = Path(__file__).resolve().parent / "data" / "idelchik1966"

# Calibration grid: psi >= 1.25 (where the correction is active), q in
# [0.1, 0.9] (away from the degenerate single-supplier endpoints), and the
# three angles the Idelchik diagrams cover.
PSI_GRID = [1.25, 1.667, 2.5, 3.333]
THETA_GRID = [30, 45, 90]
Q_GRID = [round(0.1 * i, 1) for i in range(1, 10)]

# Known Idelchik printing errors, excluded from the loss.
EXCLUDED_CELLS = {
    (30, 10.0, 0.9): "Idelchik 7-1 typo (printed 58.0 vs formula 67.9)",
    (90, 10.0, 0.7): "Idelchik 7-7 typo (printed 42.9 vs formula 49.8)",
}

Anchor = tuple[int, float, float, str, str, float]  # theta, psi, q, K_id, source, K


def mynard_K(q: float, psi: float, theta_deg: float, alpha: float) -> tuple[float, float]:
    """Mynard (K_str, K_bra) for joining type 6 at lateral fraction ``q``.

    Reference scales rho = 1, A_com = 1, m_com = 1, so K is normalised by the
    common dynamic head 0.5 * rho * u_com^2 = 0.5. Mynard's sign convention is
    U > 0 for a supplier.
    """
    A_str, A_bra, A_com = 1.0, 1.0 / psi, 1.0
    m_str, m_bra, m_com = 1.0 - q, q, 1.0
    U = np.array([m_str / A_str, m_bra / A_bra, -m_com / A_com])
    A = np.array([A_str, A_bra, A_com])
    theta = np.array([0.0, math.radians(theta_deg), math.pi])
    r = junction_loss_coefficient(U, A, theta, joining_etransfer_alpha=alpha)
    if r.K is None or len(r.K) != 2:
        return float("nan"), float("nan")
    return float(r.K[0]), float(r.K[1])


def bassett_K(q: float, psi: float, theta_deg: float) -> tuple[float, float]:
    """Bassett angle-corrected (K11, K12) at lateral fraction ``q``.

    K11 is indexed on the STRAIGHT inlet fraction (Table 1: m_A / m_C), so it
    is evaluated at ``1 - q``. K12 is indexed on the lateral fraction.
    """
    theta = math.radians(theta_deg)
    return (
        bassett2001.K11_corr(1.0 - q, psi, theta),
        bassett2001.K12_corr(q, psi, theta),
    )


def _read_q_to_K(path: Path) -> dict[float, float]:
    out: dict[float, float] = {}
    with open(path) as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            out[float(row[0])] = float(row[1])
    return out


def load_idelchik(K_id: str, theta_deg: int, psi: float) -> dict[float, float] | None:
    """Idelchik tabulated K as {q: K}, q = Q_b / Q_c for every diagram."""
    psi_label = (
        f"{int(round(psi))}"
        if abs(psi - round(psi)) < 1e-6
        else f"{psi:.3f}".rstrip("0").rstrip(".")
    )
    diag = {
        ("K12", 30): "7-1",
        ("K11", 30): "7-2",
        ("K12", 45): "7-3",
        ("K11", 45): "7-4",
        ("K12", 90): "7-7",
    }.get((K_id, theta_deg))
    if diag is None:
        return None
    path = DATA_DIR / f"idelchik_diag{diag}_{K_id}_theta={theta_deg}_psi={psi_label}_tabulated.csv"
    return _read_q_to_K(path) if path.exists() else None


def load_idelchik_K11_theta90() -> dict[float, float]:
    """K11 at theta = 90 is psi-independent in Idelchik (diagram 7-7b)."""
    return _read_q_to_K(DATA_DIR / "idelchik_diag7-7b_K11_theta=90_psi-independent_tabulated.csv")


def collect_anchors() -> list[Anchor]:
    anchors: list[Anchor] = []
    idel_K11_th90 = load_idelchik_K11_theta90()
    for theta in THETA_GRID:
        for psi in PSI_GRID:
            for q in Q_GRID:
                if (theta, psi, q) in EXCLUDED_CELLS:
                    continue
                K11_bas, K12_bas = bassett_K(q, psi, theta)
                anchors.append((theta, psi, q, "K11", "bassett", K11_bas))
                anchors.append((theta, psi, q, "K12", "bassett", K12_bas))
                if theta == 90:
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


def residuals(alpha: float, anchors: list[Anchor]) -> np.ndarray:
    """Signed residual per anchor: K_mynard(alpha) - K_anchor."""
    out = []
    for _theta, psi, q, K_id, _src, K_anchor in anchors:
        K_str, K_bra = mynard_K(q, psi, _theta, alpha)
        out.append((K_str if K_id == "K11" else K_bra) - K_anchor)
    return np.array(out)


def relative_residuals(alpha: float, anchors: list[Anchor], eps: float = 0.5) -> np.ndarray:
    """Residual over max(|K_anchor|, eps); eps keeps near-zero K11 cells finite."""
    r = residuals(alpha, anchors)
    denom = np.maximum(np.abs(np.array([a[5] for a in anchors])), eps)
    return r / denom


def loss_relative(alpha: float, anchors: list[Anchor]) -> float:
    r = relative_residuals(alpha, anchors)
    return float(np.sum(r * r))


def loss_absolute(alpha: float, anchors: list[Anchor]) -> float:
    r = residuals(alpha, anchors)
    return float(np.sum(r * r))


def fit(anchors: list[Anchor] | None = None) -> tuple[float, float]:
    """Return (alpha_abs, alpha_rel): the two loss criteria's optima."""
    anchors = anchors if anchors is not None else collect_anchors()
    res_abs = minimize_scalar(loss_absolute, bounds=(-1.0, 3.0), args=(anchors,), method="bounded")
    res_rel = minimize_scalar(loss_relative, bounds=(-1.0, 3.0), args=(anchors,), method="bounded")
    return float(res_abs.x), float(res_rel.x)


def _summary(alpha: float, name: str, anchors: list[Anchor]) -> None:
    r = residuals(alpha, anchors)
    rrel = relative_residuals(alpha, anchors)
    print(
        f"  {name:>10s} alpha={alpha:+.4f}  "
        f"mean|r|={float(np.mean(np.abs(r))):.4f}  "
        f"mean|rel|={float(np.mean(np.abs(rrel))):.4f}  "
        f"max|r|={float(np.max(np.abs(r))):.4f}"
    )


def main() -> None:
    anchors = collect_anchors()
    n_bas = sum(1 for a in anchors if a[4] == "bassett")
    n_idel = sum(1 for a in anchors if a[4] == "idelchik")
    print(f"Calibration anchors: {len(anchors)}  (bassett={n_bas}, idelchik={n_idel})")

    # Per (cell, source) at alpha=0, as mean|r| and signed mean. A signed
    # mean alone over a q-grid symmetric about 0.5 is invariant to mirroring
    # K11's axis (it averages the same set of q either way), so on its own it
    # could not have caught the axis bug this script once had.
    print("\n=== Residuals at alpha=0, by (theta, psi, K, source) ===")
    print(
        f"{'theta':>5} {'psi':>6} {'K':>4} {'source':>9}  {'N':>2}  {'mean|r|':>8}  {'signed':>8}"
    )
    r0 = residuals(0.0, anchors)
    by_cell: dict[tuple[int, float, str, str], list[float]] = defaultdict(list)
    for i, (theta, psi, _q, K_id, src, _K) in enumerate(anchors):
        by_cell[(theta, psi, K_id, src)].append(float(r0[i]))
    for k in sorted(by_cell):
        rr = np.array(by_cell[k])
        print(
            f"{k[0]:>5d} {k[1]:>6.2f} {k[2]:>4s} {k[3]:>9s}  {len(rr):>2d}  "
            f"{float(np.mean(np.abs(rr))):>8.4f}  {float(np.mean(rr)):>+8.4f}"
        )

    print("\n=== Coarse sweep ===")
    print(f"{'alpha':>7s}  {'mean|r|':>10s}  {'mean|rel|':>10s}  {'max|r|':>10s}")
    for alpha in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]:
        r = residuals(alpha, anchors)
        rrel = relative_residuals(alpha, anchors)
        print(
            f"{alpha:>7.2f}  {float(np.mean(np.abs(r))):>10.4f}  "
            f"{float(np.mean(np.abs(rrel))):>10.4f}  {float(np.max(np.abs(r))):>10.4f}"
        )

    alpha_abs, alpha_rel = fit(anchors)
    print("\n=== Optima ===")
    _summary(0.0, "baseline", anchors)
    _summary(alpha_abs, "abs-opt", anchors)
    _summary(alpha_rel, "rel-opt", anchors)


if __name__ == "__main__":
    main()
