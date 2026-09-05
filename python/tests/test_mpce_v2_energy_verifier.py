"""MPCEv2's post-solve verifier must reject a junction that creates flow work.

v1 (`MultiPortChamberElement`) has had an energy check since #229. v2 had only
a flow-direction check, so a converged state where the junction manufactures
total pressure was reported as a success (issue #271, defect 10).

The condition is the model's own dissipation identity. With the element's
residual convention the dissipation across the junction is
``q_dyn_com * sum_i |m_i| K_i``, so requiring it to be non-negative is exactly

    sum_in |m_i| Pt_i - sum_out |m_i| Pt_i >= 0

which needs only Pt and m_dot from the solution dict.

Why this form and not v1's: v1 bounds every collector's Pt by the single
supplier's, which forbids ANY branch gaining total pressure. That is real
physics in dividing flow -- Bassett's K5 and Hager's xi_t both go negative --
and v1 only gets away with it because its relative tolerance is wider than the
effect at the validation dynamic head. The mass-weighted balance permits one
branch to gain at another's expense while forbidding a net gain, and it applies
to joining flow, which v1's check explicitly skips.

Audited before being imposed (tmp/audit_realizability.py):

    measured curves        worst mass-weighted mean K  +0.21   0 violations
    Bassett analytical     worst                       +0.19   0 violations
    the MODEL              worst                       -0.06   violates for
                                                               q below ~0.2
"""

from __future__ import annotations

import combaero as cb
from combaero.network.mpce_v2_element import MPCEv2Element

_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _element(flow_direction: str = "branch") -> MPCEv2Element:
    element = MPCEv2Element.__new__(MPCEv2Element)
    element.id = "jct"
    element.N = 3
    element.port_nodes = ["p0", "p1", "p2"]
    element.port_areas = [0.01, 0.01, 0.01]
    element.port_angles_deg = [180.0, 0.0, 90.0]
    element.flow_direction = flow_direction
    element._port_signs = [-1.0, 1.0, 1.0] if flow_direction == "branch" else [-1.0, -1.0, 1.0]
    element._port_element_ids = ["e0", "e1", "e2"]
    element.strict = False
    element.joining_etransfer_alpha = 0.2
    element.jacobian_method = "sympy"
    element.penalty_alpha = 0.0
    element.eta_scale = 1.0
    return element


def _sol(m: tuple[float, float, float], pt: tuple[float, float, float]) -> dict[str, float]:
    """Canonical-direction mass flows and port total pressures."""
    return {
        "e0.m_dot": m[0],
        "e1.m_dot": m[1],
        "e2.m_dot": m[2],
        "p0.Pt": pt[0],
        "p1.Pt": pt[1],
        "p2.Pt": pt[2],
    }


# ---------------------------------------------------------------------------
# What must still pass
# ---------------------------------------------------------------------------


def test_an_ordinary_dividing_junction_passes():
    """0.1 kg/s in at 100 kPa, splitting with a loss on each leg."""
    sol = _sol((0.1, 0.06, 0.04), (100_000.0, 99_980.0, 99_950.0))

    assert _element().verify_solution_consistent(sol)


def test_a_branch_may_GAIN_total_pressure_at_the_others_expense():
    """The straight leg gains 5 Pa while the lateral loses 60.

    This is the negative-K_straight physics that Bassett (K5) and Hager (xi_t)
    both measure in dividing flow. A per-branch Pt bound would reject it; the
    mass-weighted balance must not, because the junction is still net
    dissipative: 0.06*(+5) - 0.04*(60) < 0.
    """
    sol = _sol((0.1, 0.06, 0.04), (100_000.0, 100_005.0, 99_940.0))

    assert _element().verify_solution_consistent(sol)


def test_a_lossless_junction_is_admissible():
    """Exactly zero dissipation sits on the boundary and must not be rejected."""
    sol = _sol((0.1, 0.06, 0.04), (100_000.0, 100_000.0, 100_000.0))

    assert _element().verify_solution_consistent(sol)


def test_joining_flow_is_checked_too():
    """v1's energy check skips joining outright. Two inlets, one outlet, with
    a loss on each inlet path: admissible."""
    sol = _sol((0.06, 0.04, 0.1), (100_050.0, 100_030.0, 100_000.0))

    assert _element("merge").verify_solution_consistent(sol)


# ---------------------------------------------------------------------------
# What must now be rejected
# ---------------------------------------------------------------------------


def test_a_junction_that_creates_flow_work_is_rejected():
    """Both outlets leave at a higher total pressure than the inlet arrived
    with. Nothing drives that; the state is impossible."""
    sol = _sol((0.1, 0.06, 0.04), (100_000.0, 100_050.0, 100_030.0))

    assert not _element().verify_solution_consistent(sol)


def test_a_net_gain_is_rejected_even_when_one_branch_loses():
    """The straight leg gains 100 Pa on 0.06 kg/s, the lateral loses 60 Pa on
    0.04. Weighted: +6.0 - 2.4 > 0, so the junction is a net source."""
    sol = _sol((0.1, 0.06, 0.04), (100_000.0, 100_100.0, 99_940.0))

    assert not _element().verify_solution_consistent(sol)


def test_joining_flow_that_creates_work_is_rejected():
    sol = _sol((0.06, 0.04, 0.1), (100_000.0, 100_000.0, 100_050.0))

    assert not _element("merge").verify_solution_consistent(sol)


# ---------------------------------------------------------------------------
# Boundaries of the check
# ---------------------------------------------------------------------------


def test_the_direction_check_still_fires_first():
    """A reversed port is rejected regardless of the energy balance."""
    sol = _sol((0.1, -0.06, 0.04), (100_000.0, 99_980.0, 99_950.0))

    assert not _element().verify_solution_consistent(sol)


def test_an_incomplete_solution_dict_is_not_policed():
    """Same convention as v1: absent Pt keys mean the check abstains rather
    than failing a solve it cannot see."""
    sol = _sol((0.1, 0.06, 0.04), (100_000.0, 99_980.0, 99_950.0))
    del sol["p1.Pt"]

    assert _element().verify_solution_consistent(sol)


def test_rounding_level_imbalance_is_not_a_violation():
    """A converged solve leaves a residual far below the tolerance; the gate
    must not turn numerical noise into a failed solve."""
    sol = _sol((0.1, 0.06, 0.04), (100_000.0, 100_000.0, 100_000.0 + 1e-9))

    assert _element().verify_solution_consistent(sol)


def test_a_genuine_violation_is_caught_at_a_small_flow_scale():
    small = _sol((1e-4, 6e-5, 4e-5), (100_000.0, 100_100.0, 99_940.0))

    assert not _element().verify_solution_consistent(small)


def test_rounding_level_imbalance_is_still_admissible_at_a_large_flow_scale():
    """The discriminating case for the RELATIVE tolerance.

    100 kg/s through the junction with a 1e-4 Pa imbalance is round-off, and
    it must pass. A fixed absolute threshold tight enough to catch the small
    -scale violation above would reject this one, so the tolerance has to
    scale with sum|m| * max|Pt|.
    """
    big = _sol((100.0, 60.0, 40.0), (100_000.0, 100_000.0, 100_000.0 + 1e-4))

    assert _element().verify_solution_consistent(big)
