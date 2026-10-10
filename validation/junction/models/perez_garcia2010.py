"""
Perez-Garcia 2010 K-hat linking-between-branches correlations.

Reference: J Perez-Garcia, E Sanmiguel-Rojas, A Viedma,
"New Coefficient to Characterize Energy Losses in Compressible Flow at
T-Junctions", Applied Mathematical Modelling 34:4289-4305, 2010. Read from the
accepted manuscript in docs/junction/ (gitignored).

DEFINITION (Eqs 41/42). With f(M) = (1 + (gamma-1)/2 M^2)^(gamma/(gamma-1)),
which is p0/p:

    K_hat_j = (p0_3*/p_3* - 1) / (p0_j*/p_j*) = (f(M3*) - 1) / f(Mj*)

There is NO "-1" in the denominator. Branch 3 is the common branch; the same
expression serves combining and dividing flow (Eq 41). K_hat is a function of
the two branch Mach numbers ONLY, which is why the paper says it "lacks
physical significance": a total-pressure loss reaches it only second-order,
through the branch density that sets Mj at a given mass flow. See
`miller_K_per_unit_K_hat` and validation/junction/README.md ("Tier 2") for
what that means for validation -- a K_hat comparison cannot test junction
losses.

CONVENTIONS (Fig 2, nomenclature):
    q = G2 / G3 (q' = 1 - q = G1 / G3) for EVERY flow type.
    C1  combining: 3 common outlet on the main run, 1 straight inlet,
        2 lateral inlet.                     q = lateral fraction
    C2  combining: 1 and 2 opposite inlets on the main run, 3 the leg.
    D1  dividing: 3 common inlet on the main run, 2 straight outlet,
        1 lateral outlet.                    q = STRAIGHT fraction
    D2  dividing: 3 the leg inlet, 1 and 2 opposite outlets on the main run.

Table 1 fits are to the authors' NUMERICAL results (C2/D2 compared against
their experiments, C1/D1 numerical only), 90 deg, equal areas, gamma = 1.4,
0.15 <= M3* <= 0.7, q in {0, 0.25, 0.5, 0.75, 1}. The U95 column is the
correlation's expanded uncertainty (Coleman and Steele), not a measurement
band on a loss coefficient.
"""

from __future__ import annotations

import math

# Perez-Garcia 2010 Table 1: K_hat = s * (M3*)^m * (1 + q)^(n_m1), n_m1 = n - 1
# per Eq 44 (the printed exponent on (1 + q)).
#
# Re-read 2026-10-10 from the table image at 700 dpi. The previous
# transcription carried seven wrong values: C1 K_hat_1 m 2.0263, C2 m 2.0222,
# and four it shared with the (local, gitignored) docs/junction/
# tier2_reference_data.md -- C2 n_m1
# -0.0057, D1 K_hat_1 n_m1 0.0938, D2 m 1.9027 and n_m1 -0.1338.
TABLE_1: dict[tuple[str, str], dict[str, float]] = {
    # (flow_type, K_id) -> {s, m, n_m1, U95pct, R2}
    ("C1", "K_hat_1"): {"s": 0.6871, "m": 2.0283, "n_m1": 0.1296, "U95pct": 5.75, "R2": 1.0000},
    ("C1", "K_hat_2"): {"s": 0.7559, "m": 2.0378, "n_m1": -0.0392, "U95pct": 4.35, "R2": 0.9996},
    ("C2", "K_hat_2"): {"s": 0.7504, "m": 2.0223, "n_m1": -0.0957, "U95pct": 3.54, "R2": 1.0000},
    ("D1", "K_hat_1"): {"s": 0.7312, "m": 2.0374, "n_m1": 0.0398, "U95pct": 4.64, "R2": 0.9997},
    ("D1", "K_hat_2"): {"s": 0.6718, "m": 1.9543, "n_m1": 0.0242, "U95pct": 6.25, "R2": 0.9998},
    ("D2", "K_hat_2"): {"s": 0.7307, "m": 1.9927, "n_m1": -0.1538, "U95pct": 3.99, "R2": 1.0000},
}

#: The paper fixes air at gamma = 1.4.
GAMMA = 1.4

_SYMMETRIC_TYPES = {"C2", "D2"}


def K_hat(flow_type: str, K_id: str, q: float, M3_star: float) -> float:
    """Evaluate K_hat per Eq 44.

    flow_type in {C1, C2, D1, D2}; K_id in {K_hat_1, K_hat_2}.
    For C2/D2, K_hat_1 reuses the K_hat_2 correlation with q -> (1-q)
    (Section 4.2).
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


def _f(mach: float, gamma: float = GAMMA) -> float:
    return (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** (gamma / (gamma - 1.0))


def K_hat_from_mach(M3: float, Mj: float, gamma: float = GAMMA) -> float:
    """Eq 41: K_hat from the two branch Mach numbers alone."""
    return (_f(M3, gamma) - 1.0) / _f(Mj, gamma)


def miller_K_per_unit_K_hat(M3: float, qj: float, gamma: float = GAMMA) -> float:
    """How many units of Miller's K one unit of RELATIVE K_hat change is worth.

    Equal-area combining branch j carrying the fraction qj of the common mass
    flux, adiabatic: raise its total pressure by K (p0_3 - p_3) at fixed mass
    flux and stagnation temperature, re-solve its Mach from continuity, and
    see what K_hat does. Returns dK / (dK_hat / K_hat) at K = 0, so a relative
    K_hat band b corresponds to a Miller-K band of b times this.
    """
    exponent = -(gamma + 1.0) / (2.0 * (gamma - 1.0))

    def flux(mach: float, p0: float) -> float:
        # Mass flux per sqrt(gamma / (R T0)), which cancels between branches.
        return p0 * mach * (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** exponent

    p03 = 1.0
    p3 = p03 / _f(M3, gamma)
    g_j = qj * flux(M3, p03)

    def mach_at(p0j: float) -> float:
        lo, hi = 1e-12, 1.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if flux(mid, p0j) < g_j:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    dK = 1e-3
    k0 = K_hat_from_mach(M3, mach_at(p03), gamma)
    k1 = K_hat_from_mach(M3, mach_at(p03 + dK * (p03 - p3)), gamma)
    rel = (k1 - k0) / k0
    return math.inf if rel == 0.0 else dK / abs(rel)
