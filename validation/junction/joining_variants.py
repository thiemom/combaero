"""Joining-tee closures compared on Blaisdell & Manson 1963 (pre-registered,
see data/blaisdell_manson1963/PREREGISTRATION.md).

Kernel geometry: port axes point away from the junction, U is positive into
it. Joining only, so exactly one port (the collector) has U < 0. Each closure
returns {port index: K} for the supplying ports, K being the total-pressure
loss to the collector over the collector's dynamic pressure.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np

from combaero.network._mynard2010 import _wrap_to_pi, junction_loss_coefficient


def _split(U, A, th) -> tuple[np.ndarray, np.ndarray, int, list[int], dict[int, float]]:
    U, A, th = (np.asarray(x, float) for x in (U, A, th))
    c = int(np.flatnonzero(U < 0)[0])
    sup = [i for i in range(len(U)) if i != c]
    phi = {j: abs(float(_wrap_to_pi(np.array([th[j] - th[c]]))[0])) for j in sup}
    return U, A, c, sup, phi


def bassett_cv(U, A, th, correct: bool = False, wall: bool = True) -> dict[int, float]:
    """Bassett's joining CV (Eq 3-4, Table 2) written for any set of inlets.

    Benson: the inlets share one static pressure p_s. Each inlet j enters with
    momentum rho u_j^2 F_j cos d_j along the outlet, d_j its flow deflection.
    The wall facing an inclined inlet's opening carries the reaction on its
    projected area F_j cos d_j at a pressure Bassett sets by the jet's sense:
      - jet WITH the outlet flow (cos d > 0, his type 6): the mean of p_s and
        the outlet pressure (Flamang, his Eq 4), net force (p_s-p_C)/2 F_j cos d;
      - jet AGAINST it (cos d < 0, his type 4): p_s, so the two cancel -- the
        form his K7/K8 carry.
    Continuous at 90 deg, where the projection vanishes. Then
      (p_s - p_C) D = rho u_C^2 F_C - sum rho u_j^2 F_j cos d_j,
      D = F_C + sum_{inclined j, cos d_j > 0} F_j cos d_j / 2,
    with the collinear inlet (d = 0, F = F_C) closing the main duct.
    ``correct`` applies Hager's observation that the jet turns 3/4 of the way.
    ``wall=False`` drops the opening term (D = F_C): Blaisdell & Manson's
    theory, their Eqs 12 and 19.
    """
    U, A, c, sup, phi = _split(U, A, th)
    uc = abs(U[c])
    d = {j: math.pi - phi[j] for j in sup}
    if correct:
        d = {j: 0.75 * d[j] for j in sup}
    collinear = {j: phi[j] > math.pi - 1e-9 for j in sup}
    # The wall rule follows the opening, i.e. the UNcorrected geometry.
    with_flow = {j: math.cos(math.pi - phi[j]) > 1e-12 for j in sup}
    D = A[c]
    if wall:
        D += sum(0.5 * A[j] * math.cos(d[j]) for j in sup if not collinear[j] and with_flow[j])
    M = sum((U[j] / uc) ** 2 * A[j] * math.cos(d[j]) for j in sup)
    dp = 2.0 * (A[c] - M) / D  # (p_s - p_C) / (rho u_C^2 / 2)
    return {j: dp + (U[j] / uc) ** 2 - 1.0 for j in sup}


def mynard(U, A, th, alpha: float, eta: float, recovery: float | None = None) -> dict[int, float]:
    U, A, c, sup, phi = _split(U, A, th)
    kw = {} if recovery is None else {"dividing_streamline_recovery": recovery}
    K = np.atleast_1d(junction_loss_coefficient(U, A, np.asarray(th, float), alpha, eta, **kw).K)
    return dict(zip(sup, (float(k) for k in K), strict=True))


def tee(q_lat: float, psi: float, theta_deflection_deg: float) -> tuple[np.ndarray, ...]:
    """A joining tee: collector at pi, straight inlet at 0, lateral of area
    1/psi at its deflection angle (Blaisdell & Manson's theta: 0 = with the
    outlet flow, 90 = normal, > 90 against it). Index 1 = straight inlet,
    index 2 = lateral."""
    A = np.array([1.0, 1.0, 1.0 / psi])
    U = np.array([-1.0, 1.0 - q_lat, q_lat * psi])
    return U, A, np.array([math.pi, 0.0, math.radians(theta_deflection_deg)])


Closure = Callable[[np.ndarray, np.ndarray, np.ndarray], dict[int, float]]

MODELS: dict[str, Closure] = {
    "shipped": lambda U, A, t: mynard(U, A, t, 0.2, 0.0),
    "faithful": lambda U, A, t: mynard(U, A, t, 0.0, 1.0, 0.0),
    "P2 raw": lambda U, A, t: bassett_cv(U, A, t),
    "P2c": lambda U, A, t: bassett_cv(U, A, t, correct=True),
    "BM theory": lambda U, A, t: bassett_cv(U, A, t, wall=False),
}


def bm_theory(q: float, psi: float, theta_deg: float) -> tuple[float, float]:
    """Blaisdell & Manson Eqs 12 and 19 as printed: (zeta_u, zeta_b)."""
    c = math.cos(math.radians(theta_deg))
    zu = 2.0 * q - (1.0 + 2.0 * psi * c) * q * q
    zb = -1.0 + 4.0 * q - (2.0 + 2.0 * psi * c - psi * psi) * q * q
    return zu, zb
