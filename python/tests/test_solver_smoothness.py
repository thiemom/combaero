"""Test the smoothness detector itself, against hazards with known answers.

A detector that quietly fails to detect is worse than no detector, because it
converts an unexamined function into one that looks examined. So every hazard
class is checked against a function whose behaviour is known analytically, and
every class is also checked NOT to fire on a smooth function.
"""

from __future__ import annotations

import math

import pytest

from validation.solver_smoothness import (
    Hazard,
    assert_smooth,
    render,
    scan,
)


def kinds(findings):
    return {g.hazard for g in findings}


# --------------------------------------------------------------------------
# It must not fire on smooth functions
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,f,lo,hi",
    [
        ("linear", lambda x: 3.0 * x - 1.0, -2.0, 2.0),
        ("quadratic", lambda x: x * x, -2.0, 2.0),
        ("sine", math.sin, -3.0, 3.0),
        ("exponential decay", lambda x: math.exp(-2.0 * x), 0.0, 4.0),
        ("steep but smooth", lambda x: math.tanh(40.0 * x), -0.25, 0.25),
        ("saturating to zero slope", lambda x: 1.0 - math.exp(-8.0 * x), 0.0, 1.0),
    ],
)
def test_no_false_positive_on_smooth_functions(name, f, lo, hi) -> None:
    assert scan(f, lo, hi) == [], f"{name}: {render(scan(f, lo, hi))}"


def test_a_steep_gradient_is_not_a_kink() -> None:
    """tanh(40x) has slope 40 at the origin. Steep is not the same as broken,
    and a detector that cannot tell them apart is useless on real physics."""
    assert scan(lambda x: math.tanh(40.0 * x), -0.25, 0.25) == []


def test_smooth_saturation_is_not_reported_as_a_floor() -> None:
    """A term whose derivative decays towards zero without reaching it is
    doing what a saturating physical effect should do."""
    findings = scan(lambda x: 1.0 - math.exp(-8.0 * x), 0.0, 1.0)
    assert Hazard.FLOOR not in kinds(findings)


def test_an_isolated_stationary_point_is_not_a_floor() -> None:
    """x**2 has a zero derivative at exactly one point and is well behaved.
    Width is what separates that from a clamp."""
    assert Hazard.FLOOR not in kinds(scan(lambda x: x * x, -2.0, 2.0))


def test_floating_point_saturation_IS_reported() -> None:
    """Pushed far enough, tanh saturates to exactly 1.0 and its derivative to
    exactly zero. That is a real flat spot -- a solver there sees nothing --
    so it is reported rather than excused as 'only floating point'."""
    findings = scan(lambda x: math.tanh(40.0 * x), -1.0, 1.0)
    assert Hazard.FLOOR in kinds(findings)


# --------------------------------------------------------------------------
# It must fire on each hazard class
# --------------------------------------------------------------------------


def test_detects_a_kink() -> None:
    findings = scan(lambda x: max(x, 0.0), -1.0, 1.0)
    assert Hazard.KINK in kinds(findings)
    kink = next(g for g in findings if g.hazard is Hazard.KINK)
    assert abs(kink.x) < 0.02


def test_detects_a_hard_clamp_as_both_floor_and_kinks() -> None:
    """min(max(x,0),1) is the shape every hard clamp has: flat outside,
    corners at both ends."""
    findings = scan(lambda x: min(max(x, 0.0), 1.0), -0.5, 1.5)
    assert Hazard.FLOOR in kinds(findings)
    assert Hazard.KINK in kinds(findings)
    corners = sorted(g.x for g in findings if g.hazard is Hazard.KINK)
    assert len(corners) == 2
    assert corners[0] == pytest.approx(0.0, abs=0.02)
    assert corners[1] == pytest.approx(1.0, abs=0.02)


def test_detects_divergence_and_distinguishes_it_from_a_kink() -> None:
    """sqrt(|x|) has unbounded slope at 0. That is a different failure from a
    corner and gets a different name, because the fix is different."""
    findings = scan(lambda x: math.sqrt(abs(x)), -1.0, 1.0)
    assert Hazard.DIVERGENCE in kinds(findings)
    assert Hazard.KINK not in kinds(findings)


def test_detects_divergence_sitting_on_a_domain_edge() -> None:
    """The hazard class that matters most and is easiest to miss.

    An interior probe cannot straddle an endpoint, so a singularity exactly
    at the edge of the domain is invisible to one -- and that is precisely
    where the crossflow term's derivative blows up (U1/Vi = 0 is both the
    boundary and the default operating point).
    """
    findings = scan(lambda x: math.sqrt(abs(x)), 0.0, 1.0)
    assert Hazard.DIVERGENCE in kinds(findings)
    assert any(g.x == 0.0 for g in findings if g.hazard is Hazard.DIVERGENCE)


def test_no_false_divergence_at_a_well_behaved_edge() -> None:
    """x**2 on [0, 1] has a perfectly finite derivative at both ends."""
    assert Hazard.DIVERGENCE not in kinds(scan(lambda x: x * x, 0.0, 1.0))
    assert Hazard.DIVERGENCE not in kinds(scan(math.sin, 0.0, 3.0))


def test_log_scan_probe_is_multiplicative() -> None:
    """On a log sweep an additive probe is enormous relative to x at the small
    end, which silently defeats the check there. Both references must behave:
    log(x) is smooth in x and must stay clean; a decade-wide flat region must
    still be caught as a floor measured in LOG space, not linear space.
    """
    assert scan(math.log, 1.0, 1.0e4, log=True) == []

    # Flat over [1, 10] of a [1, 1000] sweep: 1/3 of the log range, but only
    # 0.9% of the linear span. Measuring linearly would hide it.
    def clamped(x):
        return math.log(max(x, 10.0))

    findings = scan(clamped, 1.0, 1.0e3, log=True)
    assert Hazard.FLOOR in kinds(findings)


def test_detects_a_flat_region() -> None:
    findings = scan(lambda x: max(x, 0.5) if x > 0 else 0.5, 0.0, 2.0)
    assert Hazard.FLOOR in kinds(findings)


def test_floor_atol_lets_a_caller_demand_headroom() -> None:
    """Default reports only exactly-zero. A caller wanting 'not merely
    non-zero' passes a tolerance -- the distinction that caught eps=0.005
    in the expansion-factor work."""
    shallow = scan(lambda x: 1e-6 * x, 0.0, 1.0)
    assert Hazard.FLOOR not in kinds(shallow)

    strict = scan(lambda x: 1e-6 * x, 0.0, 1.0, floor_atol=1e-5)
    assert Hazard.FLOOR in kinds(strict)


def test_constant_function_is_wholly_floored() -> None:
    findings = scan(lambda x: 4.0, 0.0, 1.0)
    assert [g.hazard for g in findings] == [Hazard.FLOOR]


# --------------------------------------------------------------------------
# assert_smooth
# --------------------------------------------------------------------------


def test_assert_smooth_passes_and_fails_for_the_right_reason() -> None:
    assert_smooth(math.sin, -3.0, 3.0, label="sine")

    with pytest.raises(AssertionError, match="KINK"):
        assert_smooth(lambda x: max(x, 0.0), -1.0, 1.0, label="relu")


def test_allow_narrows_the_check_without_silencing_it() -> None:
    """A known, documented limitation should not force the whole check off."""
    f = lambda x: min(max(x, 0.0), 1.0)  # noqa: E731
    with pytest.raises(AssertionError):
        assert_smooth(f, -0.5, 1.5)
    # Accepting the floors still leaves the corners guarded.
    with pytest.raises(AssertionError, match="KINK"):
        assert_smooth(f, -0.5, 1.5, allow={Hazard.FLOOR})
    assert_smooth(f, -0.5, 1.5, allow={Hazard.FLOOR, Hazard.KINK})


def test_log_scan_works_and_rejects_nonpositive_lower_bound() -> None:
    assert scan(lambda x: math.log(x), 1.0, 1.0e4, log=True) == []
    with pytest.raises(ValueError):
        scan(lambda x: x, 0.0, 10.0, log=True)
