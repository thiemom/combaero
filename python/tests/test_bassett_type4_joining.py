"""Bassett joining flow type 4, and the leg labels that were the wrong way round.

75 measured points at area ratios 1, 2 and 4 sat unscored because K7 and K8
were never mapped to a leg. Wiring them needed two things settled from the
paper rather than from our own docstrings.

**Which leg each coefficient belongs to.** Table 1, flow type 4:

    K7 = [(p_B + rho u_B^2/2) - (p_A + rho u_A^2/2)] / (rho u_A^2/2),  q = m_B/m_A
    K8 = [(p_C + rho u_C^2/2) - (p_A + rho u_A^2/2)] / (rho u_A^2/2),  q = m_C/m_A

Every ratio and the denominator are on A, so **A is the common branch**, B is
the lateral and C the other straight leg. K7 is therefore the LATERAL-to-common
coefficient on the lateral fraction and K8 the STRAIGHT one on the straight
fraction. `bassett2001.py` documented them the other way round, and nothing
caught it because neither was ever scored.

**Which network it is.** Type 6 has the common at C instead, so in type 4 the
lateral joins pointing the other way along the main duct: the same three-port
network with the lateral mirrored to `pi - theta`. Measured over the 75 points,
that mapping gives mean errors of 0.20 and 0.31 against 0.56 and 1.02
unmirrored.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

from combaero.network._mynard2010 import junction_loss_coefficient
from validation.junction.models import bassett2001
from validation.junction.models.mpce_v2_network import MPCEv2Network
from validation.junction.network_runner import (
    _LATERAL_K_IDS,
    _STRAIGHT_K_IDS,
    iter_network_records,
)
from validation.junction.runner import _read_xy_csv
from validation.junction.schema import load_dataset

_A = 0.01


@pytest.fixture(scope="module")
def model():
    return MPCEv2Network(strict=False)


@pytest.fixture(scope="module")
def measured():
    pts = []
    for f in load_dataset().files:
        kid = f.K_id or f.coefficient or ""
        if kid in {"K7", "K8"} and f.kind == "measured":
            for x, y in _read_xy_csv(f.path):
                if 0.02 < x < 0.98:
                    pts.append((kid, x, f.psi or 1.0, f.theta_deg or 45.0, y))
    assert pts, "the type-4 measured files have gone missing"
    return pts


def _closure(q_lat: float, psi: float, theta_rad: float) -> tuple[float, float]:
    """(straight, lateral) from the raw closure in joining flow."""
    U = np.array([(1.0 - q_lat) * 10.0, q_lat * 10.0 * psi, -10.0])
    A = np.array([_A, _A / psi, _A])
    ang = np.array([0.0, theta_rad, math.pi])
    r = junction_loss_coefficient(U, A, ang)
    return float(r.K[0]), float(r.K[1])


# ---------------------------------------------------------------------------
# Which leg is which
# ---------------------------------------------------------------------------


def test_K7_is_the_lateral_coefficient_and_K8_the_straight_one():
    assert "K7" in _LATERAL_K_IDS and "K7" not in _STRAIGHT_K_IDS
    assert "K8" in _STRAIGHT_K_IDS and "K8" not in _LATERAL_K_IDS


def test_the_docstrings_say_so_too():
    """They said the opposite until this was scored, so pin the words."""
    assert "LATERAL" in bassett2001.K7_raw.__doc__
    assert "STRAIGHT" in bassett2001.K8_raw.__doc__


def test_the_leg_assignment_is_what_the_measurements_prefer(measured):
    """Decided by the data, not by reading: score both assignments."""
    errs = {"as read": [], "swapped": []}
    for kid, q_file, psi, th, meas in measured:
        tm = math.pi - math.radians(th)
        q_lat = (1.0 - q_file) if kid == "K8" else q_file
        k_str, k_lat = _closure(q_lat, psi, tm)
        errs["as read"].append(abs((k_lat if kid == "K7" else k_str) - meas))
        errs["swapped"].append(abs((k_str if kid == "K7" else k_lat) - meas))

    assert np.mean(errs["as read"]) < np.mean(errs["swapped"]), (
        "the measurements prefer the other leg assignment, so Table 1 has been read wrongly again"
    )


# ---------------------------------------------------------------------------
# Which network
# ---------------------------------------------------------------------------


def test_type_4_is_type_6_with_the_lateral_mirrored(measured):
    """Both mappings scored against the measurements; the mirror must win."""
    direct, mirrored = [], []
    for kid, q_file, psi, th, meas in measured:
        t = math.radians(th)
        q_lat = (1.0 - q_file) if kid == "K8" else q_file
        for angle, bucket in ((t, direct), (math.pi - t, mirrored)):
            k_str, k_lat = _closure(q_lat, psi, angle)
            bucket.append(abs((k_lat if kid == "K7" else k_str) - meas))

    assert np.mean(mirrored) < 0.6 * np.mean(direct), (
        f"mirrored {np.mean(mirrored):.4f} vs direct {np.mean(direct):.4f}"
    )


def test_the_analytical_forms_agree_with_the_mirror_reading():
    """Bassett's own type-4 pair should sit close to his type-6 pair evaluated
    at the mirrored angle with the legs swapped. Not identical -- the two
    derivations use different control volumes -- but close enough that the
    correspondence is not a coincidence."""
    for psi in (1.0, 2.0, 4.0):
        for q in (0.2, 0.5, 0.8):
            t = math.radians(45.0)
            tm = math.pi - t
            assert bassett2001.K7_corr(q, psi, t) == pytest.approx(
                bassett2001.K12_corr(q, psi, tm), abs=0.25
            )
            assert bassett2001.K8_corr(q, psi, t) == pytest.approx(
                bassett2001.K11_corr(q, psi, tm), abs=0.25
            )


# ---------------------------------------------------------------------------
# The axis, the trap this arc keeps returning to
# ---------------------------------------------------------------------------


def test_K8_is_on_the_straight_fraction_and_K7_is_not(model):
    """K8 at q and K7 at 1-q must build the SAME network, exactly as
    Bassett's K11 pairs with K12."""
    straight = model.evaluate_network("bassett2001", "K8", 0.3, 2.0, math.radians(45.0))
    lateral = model.evaluate_network("bassett2001", "K7", 0.7, 2.0, math.radians(45.0))

    assert straight.converged and lateral.converged
    assert straight.K_straight == pytest.approx(lateral.K_straight, rel=1e-9)
    assert straight.K_lateral == pytest.approx(lateral.K_lateral, rel=1e-9)


def test_the_achieved_split_comes_back_on_the_file_axis(model):
    q = 0.3
    r = model.evaluate_network("bassett2001", "K8", q, 2.0, math.radians(45.0))

    assert r.converged, r.message
    assert r.q_converged == pytest.approx(q, abs=1e-6)


# ---------------------------------------------------------------------------
# It is actually scored
# ---------------------------------------------------------------------------


def test_type_4_reaches_the_scorecard(model):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        records = [
            r
            for r in iter_network_records(model, load_dataset(), topologies=("imposed_q",))
            if r.K_id in {"K7", "K8"}
        ]

    assert len(records) > 70, f"only {len(records)} type-4 records reach the runner"
    assert {r.psi for r in records} == {1.0, 2.0, 4.0}, "the area-ratio sweep is incomplete"
    scored = [r for r in records if r.error is not None]
    assert len(scored) > 60


def test_the_model_degrades_with_area_ratio_here_too(model):
    """The finding this data was wired to corroborate. It was previously
    resting almost entirely on Idelchik above an area ratio of about 2.6, and
    this is an independent flow type from a different figure of the paper.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        records = [
            r
            for r in iter_network_records(model, load_dataset(), topologies=("imposed_q",))
            if r.K_id in {"K7", "K8"} and r.error is not None
        ]

    by_psi = {}
    for psi in (1.0, 2.0, 4.0):
        errs = [abs(r.error) for r in records if r.psi == psi]
        assert errs, f"no scored points at psi={psi}"
        by_psi[psi] = float(np.mean(errs))

    assert by_psi[1.0] < by_psi[2.0] < by_psi[4.0], (
        f"the degradation with area ratio has changed shape: {by_psi}"
    )
