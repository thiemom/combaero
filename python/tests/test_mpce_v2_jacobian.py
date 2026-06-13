"""
FD guardrail for the sympy-derived MPCE-v2 Jacobian.

Asserts that the analytical dKQ_i/dm_j matches a high-precision numerical
FD across the validation runner's typical (theta, psi, q) operating range.
Any divergence > 1e-3 relative trips the test loudly.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from combaero.network._mpce_v2_jacobian import dKQ_dmdot_separating_T
from validation.junction.models.mynard2010 import junction_loss_coefficient


def _fd_KQ_jacobian(
    port_mdots: np.ndarray,
    rho: np.ndarray,
    A: np.ndarray,
    theta_branch_rad: float,
) -> np.ndarray:
    """High-precision central FD reference for dKQ_i/dm_j."""
    theta_arr = np.array([math.pi, 0.0, theta_branch_rad])

    def kq(mdots):
        U = -np.array(mdots) / (rho * A)
        r = junction_loss_coefficient(U, A, theta_arr)
        K_per_port = np.array([0.0, float(r.K[0]), float(r.K[1])])
        q_dyn = 0.5 * rho[0] * (U[0] ** 2)
        return K_per_port * q_dyn

    J = np.zeros((3, 3))
    kq(port_mdots)
    for j in range(3):
        eps = max(abs(port_mdots[j]) * 1e-6, 1e-9)
        mp, mm = port_mdots.copy(), port_mdots.copy()
        mp[j] += eps
        mm[j] -= eps
        J[:, j] = (kq(mp) - kq(mm)) / (2 * eps)
    return J


CASES = [
    # (theta_deg, psi, q)
    (30, 1.0, 0.5),
    (45, 1.0, 0.5),
    (60, 1.0, 0.5),
    (90, 1.0, 0.5),
    (120, 1.0, 0.5),
    (45, 2.0, 0.5),
    (45, 3.0, 0.5),
    (45, 1.0, 0.2),
    (45, 1.0, 0.8),
    (90, 2.0, 0.3),
]


@pytest.mark.parametrize("theta_deg, psi, q", CASES)
def test_sympy_jacobian_matches_fd(theta_deg: float, psi: float, q: float) -> None:
    """Sympy analytical Jacobian must match high-precision FD within 1e-3 rel."""
    m_in = 0.1
    m_lat = q * m_in
    A_bra = 0.01 / psi
    rho = np.array([1.16, 1.16, 1.16])
    A = np.array([0.01, 0.01, A_bra])
    port_mdots = np.array([-m_in, +(m_in - m_lat), +m_lat])
    theta_branch = math.radians(theta_deg)

    J_sym = dKQ_dmdot_separating_T(port_mdots, rho, A, theta_branch)
    J_fd = _fd_KQ_jacobian(port_mdots, rho, A, theta_branch)

    # Compare nontrivial entries -- (i,j) where (i==0) is always 0 in both
    # (supplier has no loss); (1,2) and (2,1) are sparse zeros.
    for i in range(3):
        for j in range(3):
            denom = max(abs(J_fd[i, j]), 1e-6)
            rel = abs(J_sym[i, j] - J_fd[i, j]) / denom
            assert rel < 1e-3, (
                f"theta={theta_deg}, psi={psi}, q={q}: "
                f"dKQ[{i},{j}] sympy={J_sym[i, j]:.4e} fd={J_fd[i, j]:.4e} rel={rel:.2e}"
            )
