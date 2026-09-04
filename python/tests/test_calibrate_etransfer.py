"""Pins the calibration anchors for ``joining_etransfer_alpha``.

The script that produces ``DEFAULT_JOINING_ETRANSFER_ALPHA`` once fitted
against a mirrored Bassett K11 curve: it passed the lateral fraction to
``K11_corr``, which Bassett Table 1 indexes on the straight fraction. These
tests pin the axis with an identity from the source rather than with a
stored number, so the anchors cannot silently mirror again.
"""

from __future__ import annotations

import math

import pytest

from validation.junction import calibrate_etransfer as cal
from validation.junction.models import bassett2001


@pytest.mark.parametrize("theta_deg", [30, 45, 90])
@pytest.mark.parametrize("q", [0.1, 0.3, 0.5, 0.7, 0.9])
def test_bassett_anchor_satisfies_eq32_on_the_lateral_axis(theta_deg: int, q: float):
    """Bassett Eq 32: K12 - K11 = 2q - 1 at equal area, both on q = m_B / m_C.

    Table 3 applies the same theta* = (3/4) theta to both joining coefficients,
    so the identity survives the angle correction. It holds only when K11 is
    evaluated on the straight fraction 1 - q; the mirrored anchor breaks it by
    up to an order of magnitude at 30 deg.
    """
    K11, K12 = cal.bassett_K(q, 1.0, theta_deg)

    assert pytest.approx(2.0 * q - 1.0, abs=1e-12) == K12 - K11


def test_bassett_anchor_is_not_the_mirrored_form():
    """The specific defect: K11 must not be evaluated at the lateral fraction."""
    theta = math.radians(45.0)
    q = 0.3

    K11_anchor, _ = cal.bassett_K(q, 1.0, 45)
    mirrored = bassett2001.K11_corr(q, 1.0, theta)
    correct = bassett2001.K11_corr(1.0 - q, 1.0, theta)

    assert K11_anchor == pytest.approx(correct)
    assert K11_anchor != pytest.approx(mirrored)


def test_idelchik_anchors_share_one_axis():
    """Idelchik tabulates every diagram against Q_b / Q_c, K11 included.

    Both loaders must therefore return {q_lateral: K}. Guarded by checking the
    K11 and K12 tables at the same (theta, psi) share their q keys exactly.
    """
    K11 = cal.load_idelchik("K11", 45, 2.5)
    K12 = cal.load_idelchik("K12", 45, 2.5)

    assert K11 is not None and K12 is not None
    assert set(K11) == set(K12)


def test_anchor_set_is_the_documented_size():
    """432 anchors: 3 theta x 4 psi x 9 q x 2 coefficients x 2 sources, minus
    excluded typo cells (none fall on the calibration grid)."""
    anchors = cal.collect_anchors()

    assert len(anchors) == 432
    assert sum(1 for a in anchors if a[4] == "bassett") == 216
    assert sum(1 for a in anchors if a[4] == "idelchik") == 216


def test_fit_lands_in_the_documented_band():
    """The corrected-axis optima. Pinned as a band, not a point: the loss is
    shallow in alpha (mean|r| 0.210 at 0.2 vs 0.226 at 0.3), so the exact
    optimum is criterion-dependent. What must not happen is a return to the
    mirrored-axis floor (mean|r| ~ 0.41) or an optimum outside [0.2, 0.4].
    """
    anchors = cal.collect_anchors()
    alpha_abs, alpha_rel = cal.fit(anchors)

    assert 0.2 <= alpha_abs <= 0.4
    assert 0.2 <= alpha_rel <= 0.4
    baseline = float(abs(cal.residuals(0.0, anchors)).mean())
    at_opt = float(abs(cal.residuals(alpha_rel, anchors)).mean())
    assert at_opt < baseline
    assert at_opt < 0.3, "mean|r| back near the mirrored-axis floor of ~0.41"
