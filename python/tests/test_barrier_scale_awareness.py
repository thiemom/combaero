"""The soft barrier at every network size, with a fixed weight.

The barrier's weight multiplies a squared mass flow to produce a pressure::

    R_i = (Pt_i - Pt_jct) + alpha * max(0, -e_i * mdot_i)^2

so it carries Pa/(kg/s)^2, and its fixed point `slack* = sqrt(dP/alpha)` sits
at a sensible fraction of the flow for one network size only. Measured on one
junction scaled over five decades in 2026-09 (#272), the alpha needed to
converge followed `1/m_ref^2` exactly and the 1e11 fallback failed at a
hundredth of the size, so `NetworkSolver` derived the weight from the
network's own scales, `P_ref / (0.005 m_ref)^2`
(MPCE_CPP_PORT_DESIGN.md section 7e).

That hand-off is gone. Its control test said that when the derived weight
stops converging strictly more scaled junctions than the fixed one, it should
be deleted rather than relaxed, and the closure-consistent junction seed (#272,
2026-10) did exactly that: it starts every port in its declared basin, so the
barrier is seldom reached. Re-measured over 100 random junctions at sizes
1e-6, 1e-4, 1e-2, 1, 1e2 and 1e4 (600 pairs) and the whole junction scorecard,
derived and fixed weights no longer differ in a single outcome; on main before
the seed fix one pair in 400 still separated them.

What stays is the measurement the hand-off was built for: the same junction
must solve at every size. It guards scale robustness, not the weight -- even
the pre-#302 1e7 now passes it, and converges the same random draws and
scorecard records; it only changes how failing solves fail.
"""

from __future__ import annotations

import dataclasses
import random
import warnings

import pytest

from validation.junction import random_robustness as rr


@pytest.fixture(scope="module")
def base_case():
    rng = random.Random(20260906)
    for _ in range(400):
        case = rr.sample(rng)
        if rr.has_root(case) is not True:
            continue
        if (
            case.joining
            and case.drive == "flow_and_pressures"
            and 156.0 < case.theta_deg < 156.5
            and 2.5 < case.psi < 2.6
        ):
            return case
    pytest.skip("the traced case is no longer produced by this seed")


def _converges(case) -> bool:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return rr.classify(rr.build(case)) == "converged"


@pytest.mark.parametrize("factor", [1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0, 1.0e1, 1.0e2, 1.0e3])
def test_the_same_junction_converges_at_every_size(base_case, factor):
    """Seven decades of size. Scaling the port area at fixed pressure, Mach,
    angle, area ratio, split and loss coefficients leaves every dimensionless
    group untouched: it is the same junction, so it must still solve."""
    assert _converges(dataclasses.replace(base_case, area=base_case.area * factor))


def test_small_scaled_junctions_still_converge_with_the_fixed_weight():
    """The population the deleted hand-off's control ran on: 30 random
    junctions with a root, each at 1e-4 and 1e-2 of its size. Before the
    seed fix the derived weight converged 50 of 60 and the fixed one 48; with
    the seed both converge the same 50, pair for pair -- the 10 that fail are
    not barrier cases. A drop below 50 means the fixed weight has started to
    bite at small scale."""
    rng = random.Random(20260906)
    cases = []
    while len(cases) < 30:
        case = rr.sample(rng)
        if rr.has_root(case) is True:
            cases.append(case)

    converged = sum(
        _converges(dataclasses.replace(case, area=case.area * factor))
        for case in cases
        for factor in (1.0e-4, 1.0e-2)
    )
    assert converged >= 50
