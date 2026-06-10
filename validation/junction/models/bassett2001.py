"""
Bassett 2001 analytical pressure-loss coefficients for T-junctions.

Reference: M D Bassett, D E Winterbone, R J Pearson,
"Calculation of steady flow pressure loss coefficients for pipe junctions",
Proc IMechE Part C, 215(8):861-881, 2001.

Notation (matches the paper):
    q     = m_dot_other / m_dot_com  (mass flow ratio; see Table 1)
    psi   = F_C / F_B                (area ratio: main / lateral)
    theta = lateral branch angle [rad]

Six flow types (Fig 2):
    1, 2, 3 = separating (common branch is inflow)
    4, 5, 6 = joining (common branch is outflow)

The twelve K coefficients (Table 2) pair up two-per-flow-type and are derived
in Section 3 using different control volumes per type. This module exposes:

  - Table 2 raw formulas (no body-text angle correction):
        K1, K2, K3, K4, K5, K6              (separating: already include
                                             Hagar (3/4) correction)
        K7_raw, K8_raw, K11_raw, K12_raw    (joining: no angle correction)
        K9_raw, K10_raw                     (joining type 5: no correction)

  - Body-text angle-corrected forms (preferred when comparing to measured
    data per Bassett Section 4.2):
        K7_corr, K8_corr  via Eq 34:  theta' = pi - (3/4)*(pi - theta)
        K11_corr, K12_corr via Eq 33: theta' = (3/4)*theta
        K9_corr, K10_corr  via Eqs 35/37: per-branch alpha = (pi-theta)/4,
                                          beta = theta/4 to capture the
                                          two-inflow-angle deviation in
                                          type 5 joining flow

All functions return float. theta is in radians.
"""

from __future__ import annotations

import math

# ---------------------------------------------------------------------------
# Separating flows (Section 3.2). Already include Hagar (3/4)*theta via the
# derivation. No "raw" / "corrected" distinction.
# ---------------------------------------------------------------------------


def K1(q: float, psi: float, theta: float) -> float:
    """Separating type 1, main-to-branch (orientation reverse of K6)."""
    return q * q * psi * psi + 1.0 - 2.0 * q * psi * math.cos(0.75 * (math.pi - theta))


def K2(q: float, psi: float = 1.0, theta: float = 0.0) -> float:
    """Separating straight loss, type 1/2. Independent of psi and theta."""
    return q * q - 1.5 * q + 0.5


def K3(q: float, psi: float, theta: float) -> float:
    """Separating type 2, branch-to-downstream (reverse of K4)."""
    return 1.0 + q * q / (psi * psi) - (2.0 * q / psi) * math.cos(0.75 * (math.pi - theta))


def K4(q: float, psi: float, theta: float) -> float:
    """Separating type 2, branch-to-downstream."""
    return 1.0 + q * q / (psi * psi) - (2.0 * q / psi) * math.cos(0.75 * theta)


def K5(q: float, psi: float = 1.0, theta: float = 0.0) -> float:
    """Separating straight loss, type 3 (Eq 15). Identical to K2."""
    return q * q - 1.5 * q + 0.5


def K6(q: float, psi: float, theta: float) -> float:
    """Separating type 3 (Eq 27), main-to-branch. Standard reference K."""
    return q * q * psi * psi + 1.0 - 2.0 * q * psi * math.cos(0.75 * theta)


# ---------------------------------------------------------------------------
# Joining flows (Section 3.1). Table 2 raw forms + body-text angle-corrected
# forms per Section 4.2.
# ---------------------------------------------------------------------------


def K7_raw(q: float, psi: float, theta: float) -> float:
    """Joining type 4, straight path (Table 2, no angle correction)."""
    c = math.cos(theta)
    return 4.0 * q - 1.0 + q * q * (psi * psi - 2.0 + 2.0 * psi * c)


def K7_corr(q: float, psi: float, theta: float) -> float:
    """Joining type 4, straight path with Eq 34 correction:
    theta' = pi - (3/4)*(pi - theta) substituted everywhere theta appears."""
    tc = math.pi - 0.75 * (math.pi - theta)
    c = math.cos(tc)
    return 4.0 * q - 1.0 + q * q * (psi * psi - 2.0 + 2.0 * psi * c)


def K8_raw(q: float, psi: float, theta: float) -> float:
    """Joining type 4, lateral path (Table 2, no angle correction).
    K8 = 1 - q^2 + 2*(1-q)*psi*cos(theta)."""
    return 1.0 - q * q + 2.0 * (1.0 - q) * psi * math.cos(theta)


def K8_corr(q: float, psi: float, theta: float) -> float:
    """Joining type 4, lateral path with Eq 34 correction."""
    tc = math.pi - 0.75 * (math.pi - theta)
    return 1.0 - q * q + 2.0 * (1.0 - q) * psi * math.cos(tc)


def K9_raw(q: float, psi: float, theta: float) -> float:
    """Joining type 5, branch A path (Table 2 simplified form)."""
    c = math.cos(theta)
    return 1.0 + (2.0 * (2.0 * q - 1.0) / psi) * c + q * q / (psi * psi)


def K9_corr(q: float, psi: float, theta: float) -> float:
    """Joining type 5, branch A path with Eq 35 two-branch correction.
    alpha = (pi - theta)/4, beta = theta/4. See Section 4.2."""
    alpha = (math.pi - theta) / 4.0
    beta = theta / 4.0
    ca = math.cos(theta + alpha)
    cb = math.cos(theta - beta)
    return (
        2.0 * q * q * ca / psi
        - 2.0 * (1.0 - q) * (1.0 - q) * cb / psi
        + q * q / (psi * psi)
        + 1.0
    )


def K10_raw(q: float, psi: float, theta: float) -> float:
    """Joining type 5, branch C path (Table 2 simplified form)."""
    c = math.cos(theta)
    return 1.0 + (2.0 * (1.0 - 2.0 * q) / psi) * c + q * q / (psi * psi)


def K10_corr(q: float, psi: float, theta: float) -> float:
    """Joining type 5, branch C path with Eq 37 two-branch correction."""
    alpha = (math.pi - theta) / 4.0
    beta = theta / 4.0
    ca = math.cos(theta + alpha)
    cb = math.cos(theta - beta)
    return (
        2.0 * (1.0 - q) * (1.0 - q) * ca / psi
        - 2.0 * q * q * cb / psi
        + q * q / (psi * psi)
        + 1.0
    )


def K11_raw(q: float, psi: float, theta: float) -> float:
    """Joining type 6, straight-to-outlet (Table 2, no angle correction)."""
    c = math.cos(theta)
    D = psi + 0.5 * c
    N = 1.0 - q * q - (1.0 - q) * (1.0 - q) * psi * c
    return (2.0 * psi / D) * N + q * q - 1.0


def K11_corr(q: float, psi: float, theta: float) -> float:
    """Joining type 6, straight-to-outlet with Eq 33 correction."""
    c = math.cos(0.75 * theta)
    D = psi + 0.5 * c
    N = 1.0 - q * q - (1.0 - q) * (1.0 - q) * psi * c
    return (2.0 * psi / D) * N + q * q - 1.0


def K12_raw(q: float, psi: float, theta: float) -> float:
    """Joining type 6, branch-to-outlet (Table 2, no angle correction).
    K12 = 2*psi/(psi + 0.5*cos(theta)) * [1 - (1-q)^2 - q^2 * psi^2 * cos(theta)]
          + q^2 * psi^2 - 1.
    Note the psi^2 inside the bracket; this is the term where tee_junction.h
    has a transcription bug (uses psi instead of psi^2)."""
    c = math.cos(theta)
    D = psi + 0.5 * c
    N = 1.0 - (1.0 - q) * (1.0 - q) - q * q * psi * psi * c
    return (2.0 * psi / D) * N + q * q * psi * psi - 1.0


def K12_corr(q: float, psi: float, theta: float) -> float:
    """Joining type 6, branch-to-outlet with Eq 33 correction."""
    c = math.cos(0.75 * theta)
    D = psi + 0.5 * c
    N = 1.0 - (1.0 - q) * (1.0 - q) - q * q * psi * psi * c
    return (2.0 * psi / D) * N + q * q * psi * psi - 1.0


# ---------------------------------------------------------------------------
# Dispatch table -- lets the runner select a K by id and form.
# ---------------------------------------------------------------------------

K_SEPARATING: dict[str, callable] = {
    "K1": K1,
    "K2": K2,
    "K3": K3,
    "K4": K4,
    "K5": K5,
    "K6": K6,
}

K_JOINING_RAW: dict[str, callable] = {
    "K7": K7_raw,
    "K8": K8_raw,
    "K9": K9_raw,
    "K10": K10_raw,
    "K11": K11_raw,
    "K12": K12_raw,
}

K_JOINING_CORR: dict[str, callable] = {
    "K7": K7_corr,
    "K8": K8_corr,
    "K9": K9_corr,
    "K10": K10_corr,
    "K11": K11_corr,
    "K12": K12_corr,
}


def evaluate(K_id: str, q: float, psi: float, theta: float, *, angle_corrected: bool = True) -> float:
    """Evaluate any Bassett K by id.

    K_id in {K1..K12}. For separating coefficients (K1..K6) the angle_corrected
    flag is ignored (the (3/4)*theta is intrinsic to the derivation).
    For joining (K7..K12), angle_corrected=True uses the Section 4.2 body-text
    forms (Eq 33/34/35/37); False uses Table 2 raw.
    """
    if K_id in K_SEPARATING:
        return K_SEPARATING[K_id](q, psi, theta)
    table = K_JOINING_CORR if angle_corrected else K_JOINING_RAW
    if K_id not in table:
        raise ValueError(f"unknown Bassett K id: {K_id!r}")
    return table[K_id](q, psi, theta)
