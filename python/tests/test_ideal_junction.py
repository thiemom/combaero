"""Ideal momentum junction kernel (#517) against the variant already scored.

The merge side of the ideal junction is Blaisdell & Manson's theory, which the
alpha record scores as ``BM theory`` (``joining_variants.bassett_cv(...,
wall=False)``). Those numbers are the issue's pre-registered baseline, so the
kernel must reproduce that variant, not just a formula re-derived beside it.
The C++ tests (tests/test_ideal_junction.cpp) hold the closed forms, energy,
and the Jacobian; this file holds the binding to the scored code.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from combaero import _core
from validation.junction import joining_variants as jv

RHO = 1.2
# Role-smoothing bias, relative to the collector's dynamic head (see the C++
# tests for the measured size).
SMOOTH_TOL = 5e-5


def _kernel_K(U: np.ndarray, A: np.ndarray, th: np.ndarray) -> dict[int, float]:
    """K_j = (Pt_j - Pt_c) / (rho u_c^2 / 2) from the kernel, incompressible.

    joining_variants: axes point away from the junction (the kernel's
    axis_rad as is) and U is positive INTO it, so m = -rho U A.
    """
    m = [-RHO * float(u) * float(a) for u, a in zip(U, A, strict=True)]
    geom = _core.IdealJunctionGeometry()
    geom.area = [float(a) for a in A]
    geom.axis_rad = [float(t) for t in th]
    geom.port_sign = [1.0] * len(A)
    n = len(A)
    r = _core.ideal_junction_residuals_and_jacobian([0.0] * n, [RHO] * n, [0.0] * n, m, 0.0, geom)
    P = [-x for x in r.residual[:n]]
    Pt = [p + 0.5 * RHO * float(u) ** 2 for p, u in zip(P, U, strict=True)]
    c = int(np.flatnonzero(U < 0)[0])
    head = 0.5 * RHO * float(U[c]) ** 2
    return {j: (Pt[j] - Pt[c]) / head for j in range(n) if j != c}


@pytest.mark.parametrize("psi", [1.0, 2.12, 4.0, 7.11, 16.0])
@pytest.mark.parametrize("deflection_deg", [15, 30, 45, 60, 75, 90, 105, 120, 135, 150])
def test_merge_is_the_scored_bm_theory_variant(psi: float, deflection_deg: float) -> None:
    for q in np.linspace(0.05, 0.95, 10):
        U, A, th = jv.tee(float(q), psi, deflection_deg)
        ref = jv.MODELS["BM theory"](U, A, th)
        got = _kernel_K(U, A, th)
        for j, K in ref.items():
            assert got[j] == pytest.approx(K, abs=SMOOTH_TOL), (q, j)


def test_split_branch_is_idelchik_7_6_with_A_prime_one() -> None:
    """Diverging branch: zeta = 1 + (w_b/w_c)^2 - 2 (w_b/w_c) cos a, referenced
    to the arriving head (Idelchik Section VII para 9, Eq 7-6, K'_b = 0 for a
    constant-section main, A' = 1)."""
    for psi in (1.0, 3.0):
        for a_deg in (30.0, 45.0, 60.0, 90.0):
            for qb in (0.1, 0.4, 0.7):
                a = math.radians(a_deg)
                U = np.array([1.0, -(1.0 - qb), -qb * psi])  # + into the junction
                A = np.array([1.0, 1.0, 1.0 / psi])
                th = np.array([math.pi, 0.0, a])
                m = [-RHO * float(u) * float(x) for u, x in zip(U, A, strict=True)]
                geom = _core.IdealJunctionGeometry()
                geom.area, geom.axis_rad, geom.port_sign = list(A), list(th), [1.0] * 3
                r = _core.ideal_junction_residuals_and_jacobian(
                    [0.0] * 3, [RHO] * 3, [0.0] * 3, m, 0.0, geom
                )
                P = [-x for x in r.residual[:3]]
                Pt = [p + 0.5 * RHO * float(u) ** 2 for p, u in zip(P, U, strict=True)]
                w = qb * psi  # w_b / w_c
                zeta = 1.0 + w * w - 2.0 * w * math.cos(a)
                assert (Pt[0] - Pt[2]) / (0.5 * RHO) == pytest.approx(zeta, abs=SMOOTH_TOL)


def test_binding_reports_velocities_and_column_layout() -> None:
    geom = _core.IdealJunctionGeometry()
    geom.area = [1e-3, 1e-3, 5e-4]
    geom.axis_rad = [math.pi, 0.0, math.pi / 2]
    geom.port_sign = [-1.0, 1.0, 1.0]  # inlet element flows INTO the junction
    r = _core.ideal_junction_residuals_and_jacobian(
        [1e5] * 3, [1.2] * 3, [1.2e-5] * 3, [0.01, 0.006, 0.004], 1e5, geom
    )
    assert len(r.residual) == 4
    assert len(r.jacobian) == 4 and all(len(row) == 7 for row in r.jacobian)
    assert r.velocity[0] < 0.0  # junction convention: inflow is negative
    assert r.residual[3] == pytest.approx(0.0, abs=1e-15)  # -0.01 + 0.006 + 0.004
    assert r.jacobian[3][3] == pytest.approx(-1.0)  # mass row vs outer_mdot[0]
    assert r.jacobian[0][6] == pytest.approx(-1.0)  # port row vs p_jn
