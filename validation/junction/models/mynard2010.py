"""
Mynard & Valen-Sendstad 2010 Unified0D junction-loss model.

Faithful Python port of `docs/junction/JunctionLossCoefficient.m`
(BSD-licensed reference code accompanying:

    J P Mynard, K Valen-Sendstad,
    "A unified method for estimating pressure losses at vascular junctions",
    Int J Numer Methods Biomed Eng (2015), DOI:10.1002/cnm.2717).

The model handles any number of branches at any angle and any flow
direction (supplier vs collector determined by sign of U). It builds on
Bassett's CV analysis (Section 2.4.1 of the paper -- same machinery as
in `tee_junction.h::K_dat_j_closed`) and adds two refinements that the
current C++ Unified0D implementation omits:

  1. a 'pseudosupplier' branch combining all inflow, for smooth solutions
     across flow regime changes (joining flow with > 1 supplier),
  2. an empirical energy-transfer factor between diverging collectors,
     calibrated against CFD.

The Python port is structurally identical to the Matlab; variable names
follow the original where reasonable (`PseudoSupAngle`, `AreaRatio`,
`etransferfactor`). See the Matlab file for line-by-line correspondence.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def _wrap_to_pi(x: np.ndarray) -> np.ndarray:
    """Matlab's `wrapToPi`: angle to (-pi, pi]."""
    return (x + math.pi) % (2.0 * math.pi) - math.pi


def _wrap_to_2pi(x: np.ndarray) -> np.ndarray:
    """Matlab's `wrapTo2Pi`: angle to [0, 2*pi)."""
    return x % (2.0 * math.pi)


@dataclass
class MynardResult:
    """Output of one Mynard call.

    C is per-branch (collectors only; supplier entries are 0). K is per-collector
    and only computed for 3-branch junctions per the original (line 65 of the
    Matlab code). pseudo_area is per-collector (energy-transfer factor varies
    between collectors).
    """

    C: np.ndarray
    K: np.ndarray | None
    flow_ratio: np.ndarray  # |Q_collector| / Q_total_supplier, per collector
    pseudo_sup_angle: float
    pseudo_area: np.ndarray  # per-collector effective pseudosupplier area


def junction_loss_coefficient(
    U: np.ndarray,
    A: np.ndarray,
    theta: np.ndarray,
    joining_etransfer_alpha: float = 0.0,
) -> MynardResult:
    """Compute Mynard Unified0D loss coefficients for a junction.

    Args:
        U: per-branch velocities. Positive = supplier (flow INTO junction),
           negative = collector (flow OUT of junction).
        A: per-branch cross-sectional areas (same length as U).
        theta: per-branch angles relative to an arbitrary reference (radians).
        joining_etransfer_alpha: extension to Mynard for joining (converging)
            flow with asymmetric supplier areas. Mynard's original etransfer
            collapses to zero for joining (the (1 - flow_ratio) factor) since
            the single collector takes 100% of Qtot. When the supplier areas
            are asymmetric (psi != 1), the audit shows Mynard underpredicts
            K_avg vs both Bassett analytical and Idelchik tabulated; this
            term adds a joining-side correction:

                etransfer_join = alpha * (A_max - A_min) / (A_max + A_min)

            where the area asymmetry is taken over supplier ports, vanishes
            at psi=1 (equal areas), and grows monotonically toward 1 at
            extreme asymmetry. Default 0.0 preserves the faithful-port
            Mynard behavior. The single calibrated alpha is fit against
            Bassett K11_corr/K12_corr + Idelchik tables (analytical anchors)
            and validated against the held-out measured points.

    Returns:
        MynardResult with per-branch C and per-collector K (only for n<=3).
    """
    U = np.asarray(U, dtype=float).reshape(-1)
    A = np.asarray(A, dtype=float).reshape(-1)
    theta = _wrap_to_pi(np.asarray(theta, dtype=float).reshape(-1))

    Q = U * A  # volumetric flow
    Ci = Q < 0.0  # collector mask
    Si = ~Ci  # supplier mask
    Qtot = float(np.sum(Q[Si]))  # total flow into junction
    flow_ratio = -Q[Ci] / Qtot  # fraction per collector (positive)

    # ---- Reorient all branch angles so the pseudocollector angle is 0 ----
    pseudo_col_angle = float(np.mean(theta[Ci]))
    pseudo_sup_initial = math.atan2(
        float(np.sum(np.sin(theta[Si]) * Q[Si])),
        float(np.sum(np.cos(theta[Si]) * Q[Si])),
    )
    # Ensure pseudosupplier lies in the second quadrant relative to collector.
    if abs(pseudo_sup_initial - pseudo_col_angle) < math.pi / 2.0:
        pseudo_col_angle += math.pi
    theta = _wrap_to_pi(theta - pseudo_col_angle)

    # ---- Pseudosupplier angle (always positive after flipping) ----
    pseudo_direction = np.sign(np.mean(np.sin(theta[Si]) * Q[Si]))
    if pseudo_direction < 0:
        theta = -theta
    pseudo_sup_angle = math.atan2(
        float(np.sum(np.sin(np.abs(theta[Si])) * Q[Si])),
        float(np.sum(np.cos(np.abs(theta[Si])) * Q[Si])),
    )

    # ---- Effective pseudosupplier area (Mynard's empirical energy transfer) ----
    etransfer = (0.8 * (math.pi - pseudo_sup_angle) * np.sign(theta[Ci]) - 0.2) * (1.0 - flow_ratio)

    # ---- Joining-side asymmetry correction (combaero extension to Mynard) ----
    # Active only when joining_etransfer_alpha != 0 AND we have multiple
    # suppliers with different areas. Vanishes at A_max==A_min (psi=1) so
    # the original Mynard behavior is preserved in the equal-area baseline.
    if joining_etransfer_alpha != 0.0 and int(np.sum(Si)) >= 2:
        A_sup = A[Si]
        A_max = float(np.max(A_sup))
        A_min = float(np.min(A_sup))
        if A_max + A_min > 0.0:
            area_asym = (A_max - A_min) / (A_max + A_min)
            etransfer = etransfer + joining_etransfer_alpha * area_asym

    pseudo_velocity_avg = float(np.sum(U[Si] * Q[Si]) / Qtot)
    tot_pseudo_area = Qtot / ((1.0 - etransfer) * pseudo_velocity_avg)

    # ---- Area / angle ratios, then C ----
    area_ratio = tot_pseudo_area / A[Ci]
    phi = _wrap_to_2pi(pseudo_sup_angle - theta[Ci])

    C_all = np.zeros_like(U)
    # The (1 - exp(-FlowRatio/0.02)) factor avoids infinite C as FlowRatio -> 0.
    damping = 1.0 - np.exp(-flow_ratio / 0.02)
    C_all[Ci] = damping * (1.0 - (1.0 / (area_ratio * flow_ratio)) * np.cos(0.75 * (math.pi - phi)))

    # ---- K coefficients (3-branch only per the original) ----
    K = None
    if len(U) <= 3:
        Ucom = float(U[Ci][0]) if int(np.sum(Ci)) == 1 else float(U[Si][0])
        K = (U[Ci] ** 2 / Ucom**2) * (2.0 * C_all[Ci] + U[Si] ** 2 / U[Ci] ** 2 - 1.0)

    return MynardResult(
        C=C_all,
        K=K if isinstance(K, np.ndarray) else None,
        flow_ratio=flow_ratio,
        pseudo_sup_angle=pseudo_sup_angle,
        pseudo_area=np.asarray(tot_pseudo_area, dtype=float),
    )
