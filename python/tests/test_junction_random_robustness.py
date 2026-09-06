"""The random-boundary-condition robustness sweep, and its own guard rails.

`validation.junction.random_robustness` measures junction convergence on
boundary conditions drawn at random inside physical ranges, with no reference
to any paper or dataset. That independence is the whole value: the scorecard's
convergence rate is measured on cases built from Bassett's analytical K, which
is circular for a robustness claim.

The sweep is only trustworthy while its sample space stays wide, so the ranges
are pinned here. **Narrowing the space would raise the convergence number
without improving anything, and doing it requires editing an assertion in this
file, which shows up in a diff.** That is deliberate.

The convergence floor is deliberately loose for the same reason. It is there to
catch a real regression, not to be optimised against.
"""

from __future__ import annotations

import math
import random

import pytest

from validation.junction import random_robustness as rr

# The sweep is a solve per draw, so keep the routine size small; the full run
# is `uv run python -m validation.junction.random_robustness 2000`.
_N = 120
_SEED = 20260906


@pytest.fixture(scope="module")
def summary():
    return rr.run(n=_N, seed=_SEED)


# ---------------------------------------------------------------------------
# The sample space stays wide
# ---------------------------------------------------------------------------


def test_the_pressure_and_temperature_span_stays_wide():
    assert rr.PT_RANGE_PA[0] <= 0.5e5 and rr.PT_RANGE_PA[1] >= 10.0e5
    assert rr.TT_RANGE_K[0] <= 250.0 and rr.TT_RANGE_K[1] >= 600.0


def test_the_mach_range_reaches_past_where_the_closure_is_trusted():
    """The closure is documented for low Mach. The sweep deliberately samples
    beyond that so the report says what happens at the edge rather than
    quietly excluding it."""
    assert rr.MACH_RANGE[0] <= 0.005
    assert rr.MACH_RANGE[1] >= 0.25


def test_the_geometry_span_covers_both_sides_of_every_ratio():
    assert rr.THETA_RANGE_DEG[0] <= 15.0 and rr.THETA_RANGE_DEG[1] >= 165.0
    # a branch both larger and smaller than the main duct
    assert rr.PSI_RANGE[0] < 1.0 < rr.PSI_RANGE[1]
    assert rr.PSI_RANGE[1] / rr.PSI_RANGE[0] >= 30.0
    # areas spanning at least two decades
    assert rr.AREA_RANGE_M2[1] / rr.AREA_RANGE_M2[0] >= 100.0
    assert rr.SPLIT_RANGE[0] <= 0.05 and rr.SPLIT_RANGE[1] >= 0.95


def test_the_imposed_pressure_differences_admit_a_negative_coefficient():
    """Both Bassett and Hager measure a mildly negative straight-leg
    coefficient in dividing flow, so the sweep must be able to ask for one."""
    assert rr.K_STRAIGHT_RANGE[0] < 0.0
    assert rr.K_BRANCH_RANGE[0] < 0.0


def test_every_drive_and_both_flow_directions_are_sampled():
    rng = random.Random(1)
    cases = [rr.sample(rng) for _ in range(300)]

    assert {c.drive for c in cases} == set(rr.DRIVES)
    assert {c.joining for c in cases} == {True, False}


def test_draws_land_inside_the_declared_ranges():
    rng = random.Random(2)
    for _ in range(200):
        c = rr.sample(rng)
        assert rr.PT_RANGE_PA[0] <= c.Pt <= rr.PT_RANGE_PA[1]
        assert rr.TT_RANGE_K[0] <= c.Tt <= rr.TT_RANGE_K[1]
        assert rr.MACH_RANGE[0] <= c.mach <= rr.MACH_RANGE[1]
        assert rr.THETA_RANGE_DEG[0] <= c.theta_deg <= rr.THETA_RANGE_DEG[1]
        assert rr.PSI_RANGE[0] <= c.psi <= rr.PSI_RANGE[1]
        assert rr.AREA_RANGE_M2[0] <= c.area <= rr.AREA_RANGE_M2[1]


# ---------------------------------------------------------------------------
# Feasibility is separated from failure
# ---------------------------------------------------------------------------


def _case(**kw) -> rr.Case:
    base = {
        "Pt": 1.0e5,
        "Tt": 300.0,
        "mach": 0.03,
        "theta_deg": 90.0,
        "psi": 1.0,
        "split": 0.5,
        "area": 0.01,
        "drive": "imposed_flows",
        "joining": False,
        "k_straight": 0.2,
        "k_branch": 1.0,
    }
    base.update(kw)
    return rr.Case(**base)


def test_imposing_both_flows_always_admits_a_root():
    """The split is a boundary condition there, so there is nothing to solve
    for and no way to ask for something unattainable."""
    assert rr.has_root(_case(drive="imposed_flows")) is True


def test_an_unattainable_pressure_difference_has_no_root():
    """Nothing in the closure produces a coefficient difference of 500."""
    assert rr.has_root(_case(drive="flow_and_pressures", k_straight=0.0, k_branch=500.0)) is False


def test_an_attainable_pressure_difference_does():
    assert rr.has_root(_case(drive="flow_and_pressures", k_straight=0.0, k_branch=1.0)) is True


def test_a_draw_with_no_root_is_not_counted_as_a_solver_failure(summary):
    """The distinction this whole harness exists to make."""
    no_root = summary.by_cut[("root", "none")]
    assert sum(no_root.values()) > 0, "the sweep must reach infeasible draws at all"
    assert no_root["no progress"] == 0, (
        "a draw with no root must be reported as such, not as the solver stalling"
    )


# ---------------------------------------------------------------------------
# What the sweep reports
# ---------------------------------------------------------------------------


def test_the_solver_never_raises(summary):
    """The one hard requirement. Any exception is a defect, whatever the
    boundary conditions."""
    assert summary.outcomes["raised"] == 0


def test_most_solvable_draws_converge(summary):
    """A LOOSE floor. Measured at 84.5% over 2000 draws and 84.7% over 300;
    the floor sits well below so it catches a regression rather than inviting
    anyone to optimise against it. If this fails, look at what changed in the
    closure or the solver, not at the ranges above.
    """
    assert summary.n_with_root >= _N // 3, "too few solvable draws to conclude anything"
    assert summary.converged_share_where_solvable > 0.70


def test_the_sweep_is_reproducible():
    """Same seed, same answer. The solver's hash-order dependence was fixed in
    #289; a failure here would mean it has come back."""
    a = rr.run(n=25, seed=7)
    b = rr.run(n=25, seed=7)

    assert a.outcomes == b.outcomes


def test_different_seeds_explore_different_cases():
    a = rr.run(n=25, seed=7)
    b = rr.run(n=25, seed=8)

    assert a.outcomes != b.outcomes or a.n_with_root != b.n_with_root


def test_the_report_names_every_outcome_it_counted(summary):
    text = rr.format_summary(summary)

    for name, count in summary.outcomes.items():
        if count:
            assert name in text, f"{name} counted but not reported"
    assert "admit a root" in text


def test_pressure_driven_cases_are_the_weak_ones(summary):
    """Documented state, not a target: with all three pressures imposed the
    flow level is free as well as the split, and that is where the solver
    struggles. Recorded so an improvement is noticed as much as a regression.
    """
    free_level = summary.by_cut[("solvable drive", "all_pressures")]
    pinned = summary.by_cut[("solvable drive", "imposed_flows")]
    if not sum(free_level.values()) or not sum(pinned.values()):
        pytest.skip("sweep too small to compare drives")

    share_free = free_level["converged"] / sum(free_level.values())
    share_pinned = pinned["converged"] / sum(pinned.values())
    assert share_free < share_pinned


def test_the_closure_curve_covers_the_split_range():
    straight, branch = rr._closure_curve(_case())

    assert straight.size == branch.size == rr._Q_GRID.size
    assert math.isfinite(float(straight[len(straight) // 2]))
