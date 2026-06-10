"""
Hager 1984 analytical loss coefficients for branched ducts.

Reference: W H Hager, "An approximate treatment of flow in branches and bends",
Proc IMechE Part C, 198C(4):63-69, 1984.

Foundational paper introducing the (3/4)*delta angle correction adopted by
Bassett. Considers equal-area T-junctions only (psi = 1 implicit). Both
formulas reduce to special cases of Bassett's K2/K5 (straight) and K6 (lateral)
at psi = 1.

Notation:
    q     = Delta_Q / Q  (mass flow ratio; lateral / total)
    delta = lateral branch angle [rad]

Variables:
    xi_t  = straight-through loss coefficient (= K5 / K2 of Bassett at psi=1)
    xi_l  = lateral loss coefficient (= K6 of Bassett at psi=1)
"""

from __future__ import annotations

import math


def xi_t(q: float, delta: float = math.pi / 2.0) -> float:
    """Straight-through loss coefficient (Hager Eq 8). Independent of delta.

    xi_t = q * (q - 1/2)

    Equivalent to Bassett K5 = q^2 - 1.5*q + 0.5 modulo the half-difference
    in stagnation-pressure normalization; numerically xi_t and K5 differ by
    a constant offset of 0.5. The shape (and derivative wrt q) is identical.
    """
    return q * (q - 0.5)


def xi_l(q: float, delta: float) -> float:
    """Lateral loss coefficient (Hager Eq 19). delta in radians.

    xi_l = 1 - 2*q*cos((3/4)*delta) + q^2

    Identical to Bassett K6 at psi = 1.
    """
    return 1.0 - 2.0 * q * math.cos(0.75 * delta) + q * q


def xi_l_min_q(delta: float) -> float:
    """Mass flow ratio q at which xi_l(q) attains its minimum (Hager Eq 21).

    Useful as a diagnostic / regression check on digitized lateral-loss data.
    """
    return math.cos(0.75 * delta)


def xi_l_min(delta: float) -> float:
    """Minimum value of xi_l(q) (Hager Eq 21). xi_l_min = 1 - cos^2((3/4)*delta)."""
    c = math.cos(0.75 * delta)
    return 1.0 - c * c
