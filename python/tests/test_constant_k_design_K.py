"""ConstantKTeeElement.bassett_design_K: the constant-K default.

The production copy must agree with the scored Bassett model in
validation/junction/models/bassett2001.py (which production code may not
import), with each coefficient read at its OWN leg's fraction. The golden
values are shared with gui/frontend/src/utils/teeK.test.ts, whose display
mirror must agree too.
"""

from __future__ import annotations

import math

import pytest

from combaero.network.mpce_element import ConstantKTeeElement
from validation.junction.models import bassett2001

CASES = [
    (q, psi, th) for q in (0.1, 0.5, 0.9) for psi in (1.0, 2.0, 4.0) for th in (30.0, 60.0, 90.0)
]


@pytest.mark.parametrize("q, psi, theta_deg", CASES)
def test_dividing_pair_is_bassett_k5_k6(q: float, psi: float, theta_deg: float) -> None:
    got = ConstantKTeeElement.bassett_design_K("branch", q, psi, theta_deg)
    want = bassett2001.separating_pair_at(q, psi, math.radians(theta_deg))
    assert got == pytest.approx(want, abs=1e-12)


@pytest.mark.parametrize("q, psi, theta_deg", CASES)
def test_joining_pair_is_bassett_k11_k12_angle_corrected(
    q: float, psi: float, theta_deg: float
) -> None:
    t = math.radians(theta_deg)
    got = ConstantKTeeElement.bassett_design_K("merge", q, psi, theta_deg)
    want = (bassett2001.K11_corr(1.0 - q, psi, t), bassett2001.K12_corr(q, psi, t))
    assert got == pytest.approx(want, abs=1e-12)


GOLDEN = [
    ("branch", 0.5, 4.0, 90.0, (0.0, 3.469266)),
    ("branch", 0.2, 1.0, 45.0, (-0.06, 0.707412)),
    ("merge", 0.5, 4.0, 90.0, (-0.048904, 3.701096)),
    ("merge", 0.3, 1.0, 60.0, (0.149539, -0.250461)),
]


@pytest.mark.parametrize("flow_direction, q, psi, theta_deg, want", GOLDEN)
def test_golden_values_shared_with_the_frontend(
    flow_direction: str, q: float, psi: float, theta_deg: float, want: tuple[float, float]
) -> None:
    got = ConstantKTeeElement.bassett_design_K(flow_direction, q, psi, theta_deg)
    assert got == pytest.approx(want, abs=1e-6)


@pytest.mark.parametrize("q", [0.0, 1.0, -0.2, 1.2])
def test_refuses_a_share_outside_the_open_unit_interval(q: float) -> None:
    with pytest.raises(ValueError, match="q_branch"):
        ConstantKTeeElement.bassett_design_K("branch", q, 1.0, 90.0)


def test_design_net_loss_is_the_mass_weighted_pair() -> None:
    assert ConstantKTeeElement.design_net_loss(0.3, 0.2, -0.1) == pytest.approx(0.11)


@pytest.mark.parametrize(
    "flow_direction, q, psi, theta_deg, negative",
    [
        # #779: Bassett's joining pair at equal areas and shallow angles
        # (measured at 30 deg, q = 0.5: about +0.05, Blaisdell & Manson).
        ("merge", 0.5, 1.0, 30.0, True),
        ("merge", 0.3, 1.0, 45.0, True),
        ("merge", 0.5, 1.0, 90.0, False),
        ("merge", 0.5, 2.0, 30.0, False),
        # Dividing into a narrow, shallow branch.
        ("branch", 0.25, 4.0, 15.0, True),
        ("branch", 0.5, 1.0, 90.0, False),
    ],
)
def test_bassett_pairs_that_create_total_pressure(
    flow_direction: str, q: float, psi: float, theta_deg: float, negative: bool
) -> None:
    pair = ConstantKTeeElement.bassett_design_K(flow_direction, q, psi, theta_deg)
    assert (ConstantKTeeElement.design_net_loss(q, *pair) < 0.0) is negative


def test_design_net_loss_refuses_a_share_outside_the_open_unit_interval() -> None:
    with pytest.raises(ValueError, match="q_branch"):
        ConstantKTeeElement.design_net_loss(1.0, 0.1, 0.1)
