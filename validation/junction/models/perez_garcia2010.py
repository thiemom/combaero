"""
Perez-Garcia 2010 K-hat linking-between-branches correlations.

Reference: J Perez-Garcia, E Sanmiguel-Rojas, A Viedma,
"New Coefficient to Characterize Energy Losses in Compressible Flow at
T-Junctions", Applied Mathematical Modelling 34:4289-4305, 2010.

The paper introduces a "linking" coefficient K_hat that has lower uncertainty
amplification than Miller's standard total-pressure loss coefficient K for
compressible flow. Global correlations (Eq 44) of the form

    K_hat = s * (M3*)^m * (1+q)^(n-1)

are fitted per flow type. M3* is the extrapolated Mach number in the common
branch (after frictional losses are subtracted), q is the mass flow ratio.

Four flow types are tabulated (Fig 2):
    C1 = combining, branches 3 and 1 (K_hat_1)
    C2 = combining, branches 3 and 2 (K_hat_2)
    D1 = dividing, branches 3 and 1 (K_hat_1)
    D2 = dividing, branches 3 and 2 (K_hat_2)

The paper notes (Section 4.1) that K_hat_1 correlations can be obtained from
K_hat_2 by substituting q' = (1 - q) for q.

Range of applicability (Section 4.1):
    0.15 <= M3* <= 0.7
    q in {0, 0.25, 0.5, 0.75, 1}  (numerically tested)
    extrapolation valid within the linearity-of-regression-plane assumption
"""

from __future__ import annotations


# Perez-Garcia 2010 Table 1: K_hat = s * (M3*)^m * (1+q)^(n_m1), where
# n_m1 = n - 1 is the exponent on (1+q) per Eq 44.
#
# Symmetric flow types (C2, D2): only K_hat_2 is tabulated; K_hat_1 follows
# from the same correlation with q -> (1 - q) substitution (paper Section 4.1).
# Asymmetric flow types (C1, D1): K_hat_1 and K_hat_2 are tabulated separately.
#
# Source: Perez-Garcia, Sanmiguel-Rojas, Viedma, Appl Math Model 34 (2010)
# 4289-4305, Table 1.
TABLE_1: dict[tuple[str, str], dict[str, float]] = {
    # (flow_type, K_id) -> {s, m, n_m1, U95pct, R2}
    ("C1", "K_hat_1"): {"s": 0.6871, "m": 2.0263, "n_m1": 0.1296, "U95pct": 5.75, "R2": 1.0000},
    ("C1", "K_hat_2"): {"s": 0.7559, "m": 2.0378, "n_m1": -0.0392, "U95pct": 4.35, "R2": 0.9996},
    ("C2", "K_hat_2"): {"s": 0.7504, "m": 2.0222, "n_m1": -0.0057, "U95pct": 3.54, "R2": 1.0000},
    ("D1", "K_hat_1"): {"s": 0.7312, "m": 2.0374, "n_m1": 0.0938, "U95pct": 4.64, "R2": 0.9997},
    ("D1", "K_hat_2"): {"s": 0.6718, "m": 1.9543, "n_m1": 0.0242, "U95pct": 6.25, "R2": 0.9998},
    ("D2", "K_hat_2"): {"s": 0.7307, "m": 1.9027, "n_m1": -0.1338, "U95pct": 3.99, "R2": 1.0000},
}


_SYMMETRIC_TYPES = {"C2", "D2"}


def K_hat(flow_type: str, K_id: str, q: float, M3_star: float) -> float:
    """Evaluate K_hat per Eq 44.

    flow_type in {C1, C2, D1, D2}; K_id in {K_hat_1, K_hat_2}.
    For C2/D2, K_hat_1 reuses the K_hat_2 correlation with q -> (1-q).
    """
    key = (flow_type, K_id)
    if key in TABLE_1:
        p = TABLE_1[key]
        return p["s"] * (M3_star ** p["m"]) * ((1.0 + q) ** p["n_m1"])
    if flow_type in _SYMMETRIC_TYPES and K_id == "K_hat_1":
        p = TABLE_1[(flow_type, "K_hat_2")]
        q_eff = 1.0 - q
        return p["s"] * (M3_star ** p["m"]) * ((1.0 + q_eff) ** p["n_m1"])
    raise ValueError(f"unknown Perez-Garcia (flow_type, K_id): {key}")
