"""The structure of D(q) = K_lateral - K_straight, and what it settles.

The pressure-driven topologies constrain this difference, not the two
coefficients separately, so its shape decides whether an operating point
exists and whether it is unique (issue #271, record in
docs/archive/JUNCTION_OPERATING_POINT_271.md Finding 6).

Bassett's separating pair has closed forms, so the difference does too --
but each coefficient is indexed on the fraction in ITS OWN leg (#295), so the
straight one is evaluated at 1 - q when q is the lateral fraction:

    K5(x)         = x^2 - 1.5 x + 0.5, x the STRAIGHT fraction  (Eq 15)
    K6(q,psi,th)  = q^2 psi^2 + 1 - 2 q psi cos(0.75 th)        (Eq 27)
    D(q)          = K6(q) - K5(1-q)
                  = q^2 (psi^2 - 1) + q (0.5 - 2 psi c) + 1

CORRECTED 2026-09-05. This file first paired the two at the same q and gave
`... + q (1.5 - 2 psi c) + 0.5`. The tests were self-consistent and green
while measuring a quantity that is not the physical difference. The
conclusions below survive -- linear at psi=1, a cup with its vertex inside
(0,1) otherwise -- but every number moved.

Two things follow that are worth pinning, because a change to either
coefficient can silently break them:

* at psi = 1 the quadratic term vanishes and D is LINEAR, so the operating
  point is unique for every equal-area geometry;
* at psi != 1 D is a cup, and the second root is the mirror about its vertex,
  which is physical only for low targets.
"""

from __future__ import annotations

import math

import pytest

from validation.junction.models import bassett2001
from validation.junction.models.mpce_network import MPCENetwork

_GEOMS = [(1.0, 45), (1.0, 90), (1.0, 120), (2.0, 45), (3.0, 45), (3.333, 90)]


def _D_closed(q: float, psi: float, theta: float) -> float:
    c = math.cos(0.75 * theta)
    return q * q * (psi * psi - 1.0) + q * (0.5 - 2.0 * psi * c) + 1.0


def _D_bassett(q: float, psi: float, theta: float) -> float:
    """q is the LATERAL fraction, so the straight leg is read at 1 - q."""
    return bassett2001.K6(q, psi, theta) - bassett2001.K5(1.0 - q)


@pytest.mark.parametrize("psi, theta_deg", _GEOMS)
@pytest.mark.parametrize("q", [0.0, 0.2, 0.5, 0.8, 1.0])
def test_the_difference_is_the_quadratic_the_analysis_assumes(q, psi, theta_deg):
    theta = math.radians(theta_deg)
    assert _D_bassett(q, psi, theta) == pytest.approx(_D_closed(q, psi, theta), abs=1e-12)


@pytest.mark.parametrize("theta_deg", [30, 45, 60, 90, 120, 180])
def test_at_equal_areas_the_difference_is_linear_in_q(theta_deg):
    """psi = 1 kills the quadratic term, so a pressure-driven equal-area tee
    has exactly one operating point. 82% of the separating dataset is here."""
    theta = math.radians(theta_deg)
    slope = 0.5 - 2.0 * math.cos(0.75 * theta)

    for q in (0.0, 0.25, 0.5, 0.75, 1.0):
        assert _D_bassett(q, 1.0, theta) == pytest.approx(1.0 + slope * q, abs=1e-12)


@pytest.mark.parametrize(
    "psi, theta_deg, expected_vertex",
    [(2.0, 45, 0.471), (3.0, 45, 0.281), (4.0, 45, 0.205), (3.333, 90, 0.101)],
)
def test_at_unequal_areas_the_vertex_is_inside_the_operating_range(psi, theta_deg, expected_vertex):
    """So Bassett's own difference is NOT monotone there. Monotonicity is not a
    property of the physics and cannot be an objective for the model."""
    theta = math.radians(theta_deg)
    a = psi * psi - 1.0
    b = 0.5 - 2.0 * psi * math.cos(0.75 * theta)
    vertex = -b / (2.0 * a)

    assert vertex == pytest.approx(expected_vertex, abs=0.001)
    assert 0.0 < vertex < 1.0
    assert _D_bassett(vertex, psi, theta) < min(
        _D_bassett(0.0, psi, theta), _D_bassett(1.0, psi, theta)
    )


@pytest.mark.parametrize("q_t, has_second_root", [(0.2, True), (0.5, True), (0.7, False)])
def test_the_second_operating_point_is_the_mirror_about_the_vertex(q_t, has_second_root):
    """psi=3, theta=45: vertex 0.281, so a second physical root needs
    q_t < 0.562. The rest have a unique operating point under the true
    coefficients.
    """
    psi, theta = 3.0, math.radians(45.0)
    a = psi * psi - 1.0
    b = 0.5 - 2.0 * psi * math.cos(0.75 * theta)
    q2 = -b / a - q_t

    assert (0.0 < q2 < 1.0) is has_second_root
    if has_second_root:
        assert _D_bassett(q2, psi, theta) == pytest.approx(_D_bassett(q_t, psi, theta), abs=1e-9)


# ---------------------------------------------------------------------------
# Formerly the target for #272, now satisfied
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [45, 90, 120])
def test_model_reproduces_the_equal_area_identity(theta_deg):
    """WAS a strict xfail target; the model satisfies it since the
    dividing-streamline recovery landed (2026-09-05).

    At equal areas the quadratic term vanishes and the source gives a straight
    line whose slope and intercept are fixed by geometry alone. The model used
    to put a HUMP here, turning at q=0.30 for theta=90 and q=0.50 for
    theta=120, which manufactured a second operating point in the 82% of the
    separating dataset where the physics has exactly one.
    """
    theta = math.radians(theta_deg)
    slope = 0.5 - 2.0 * math.cos(0.75 * theta)
    model = MPCENetwork(strict=False)

    for q in (0.2, 0.5, 0.8):
        r = model.evaluate_network("bassett2001", "K6", q, 1.0, theta)
        assert r.converged, r.message
        d = r.K_lateral - r.K_straight
        assert d == pytest.approx(1.0 + slope * q, abs=0.02)
