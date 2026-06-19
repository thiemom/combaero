"""
Sympy-derived analytical Jacobian for the MPCE-v2 Mynard K*q_dyn term.

Covers the 3-port separating T-junction case (1 supplier at port 0,
collectors at ports 1, 2). The structure is sparse:

    dKQ_0/dm_j   = 0 for all j (supplier port has no loss contribution)
    dKQ_1/dm_2   = 0           (port 1 K does not see port 2's mdot)
    dKQ_2/dm_1   = 0           (port 2 K does not see port 1's mdot)

So only six derivatives are nontrivial:
    dKQ_1/dm_0, dKQ_1/dm_1, dKQ_2/dm_0, dKQ_2/dm_2

These are lambdified once at import time. Runtime cost ~= a single
evaluation of the polynomial/exp/cos expression -- no Mynard call per
perturbation, no FD overhead.

For configurations outside this canonical 3-port separating layout
(joining flow, more ports, asymmetric angles), the MPCE-v2 element
falls back to FD.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import sympy as sp


def _build_lambdified() -> dict[tuple[int, int], Callable]:
    """Return a dict mapping (port_i, port_j) -> dKQ_i/dm_j function.

    The function signature is (m0, m1, m2, rho0, rho1, rho2, A0, A1, A2, th2).
    """
    m0, m1, m2 = sp.symbols("m0 m1 m2", real=True)
    rho0, rho1, rho2 = sp.symbols("rho0 rho1 rho2", positive=True)
    A0, A1, A2 = sp.symbols("A0 A1 A2", positive=True)
    th2 = sp.symbols("th2", positive=True)

    U0 = -m0 / (rho0 * A0)
    U1 = -m1 / (rho1 * A1)
    U2 = -m2 / (rho2 * A2)

    Q0 = U0 * A0
    Q1 = U1 * A1
    Q2 = U2 * A2
    Qtot = Q0

    FR1 = -Q1 / Qtot
    FR2 = -Q2 / Qtot

    psi_sup = sp.pi - th2 / 2
    sign1 = -1
    sign2 = +1

    etf1 = (sp.Rational(8, 10) * (sp.pi - psi_sup) * sign1 - sp.Rational(2, 10)) * (1 - FR1)
    etf2 = (sp.Rational(8, 10) * (sp.pi - psi_sup) * sign2 - sp.Rational(2, 10)) * (1 - FR2)

    u_pseudo = U0
    tpa1 = Qtot / ((1 - etf1) * u_pseudo)
    tpa2 = Qtot / ((1 - etf2) * u_pseudo)

    AR1 = tpa1 / A1
    AR2 = tpa2 / A2

    phi1 = psi_sup - sign1 * th2 / 2
    phi2 = psi_sup - sign2 * th2 / 2

    damp1 = 1 - sp.exp(-FR1 / sp.Rational(2, 100))
    damp2 = 1 - sp.exp(-FR2 / sp.Rational(2, 100))

    C1 = damp1 * (1 - (1 / (AR1 * FR1)) * sp.cos(sp.Rational(3, 4) * (sp.pi - phi1)))
    C2 = damp2 * (1 - (1 / (AR2 * FR2)) * sp.cos(sp.Rational(3, 4) * (sp.pi - phi2)))

    Ucom = U0
    K1 = (U1**2 / Ucom**2) * (2 * C1 + Ucom**2 / U1**2 - 1)
    K2 = (U2**2 / Ucom**2) * (2 * C2 + Ucom**2 / U2**2 - 1)

    q_dyn = sp.Rational(1, 2) * rho0 * Ucom**2
    KQ1 = K1 * q_dyn
    KQ2 = K2 * q_dyn

    args = (m0, m1, m2, rho0, rho1, rho2, A0, A1, A2, th2)
    nontrivial = {
        (1, 0): sp.diff(KQ1, m0),
        (1, 1): sp.diff(KQ1, m1),
        (2, 0): sp.diff(KQ2, m0),
        (2, 2): sp.diff(KQ2, m2),
    }
    return {key: sp.lambdify(args, expr, modules=["numpy"]) for key, expr in nontrivial.items()}


# Built once at import.
_DKQ_FNS = _build_lambdified()


def dKQ_dmdot_separating_T(
    port_mdots: np.ndarray,
    rho: np.ndarray,
    A: np.ndarray,
    theta_branch_rad: float,
) -> np.ndarray:
    """Return the 3x3 dKQ_i/dm_j Jacobian for the canonical 3-port separating T.

    Assumes port 0 = supplier (canonical inlet), ports 1, 2 = collectors with
    straight at theta=0 and branch at theta_branch_rad.
    """
    m0, m1, m2 = float(port_mdots[0]), float(port_mdots[1]), float(port_mdots[2])
    rho0, rho1, rho2 = float(rho[0]), float(rho[1]), float(rho[2])
    A0, A1, A2 = float(A[0]), float(A[1]), float(A[2])
    th2 = float(theta_branch_rad)
    args = (m0, m1, m2, rho0, rho1, rho2, A0, A1, A2, th2)
    J = np.zeros((3, 3))
    for (i, j), fn in _DKQ_FNS.items():
        try:
            J[i, j] = float(fn(*args))
        except (ZeroDivisionError, FloatingPointError):
            J[i, j] = 0.0
    return J
