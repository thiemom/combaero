"""The soft barrier had a fixed point, and solves parked in it.

`MPCEv2Element` refuses a port flowing against its declared direction. In
non-strict mode it replaces the physics with a "soft barrier" whose docstring
promises it "pulls Newton back toward mdot_i = 0, from which a sign flip
restores the strict-physics residual on the next iteration". Measured, it did
the opposite.

The penalty is added INTO the continuity row rather than occupying one of its
own::

    R_i = (Pt_i - Pt_jct) + alpha * max(0, -e_i * mdot_i)^2

so the solver can zero that row by carrying a pressure error equal and
opposite to the penalty instead of by driving the slack to zero. The barrier
is then not a barrier: it is a fabricated pressure loss the network happily
accommodates, with a fixed point at

    slack* = sqrt(dP / alpha)

where dP is the pressure error the surrounding network can absorb. At the old
alpha of 1e7 that put slack* at 45% of the common mass flow -- nowhere near
the mdot = 0 the sign flip needs -- and a traced case sat there at a residual
floor of 1.6e4 Pa while an in-regime root existed and converged to |F| = 1e-5
when seeded at it.

The fix is the weight, not the form: alpha = 1e11 moves the fixed point to
0.45% of the common flow. Across the random boundary harness the response is
monotone in alpha and flat from 1e9 up, so this is a saturation point and not
a tuned optimum.

Issue #272. Full record in `validation/junction/MPCE_CPP_PORT_DESIGN.md`
section 7c.
"""

from __future__ import annotations

import math
import random
import warnings

import numpy as np
import pytest

from combaero.network import NetworkSolver
from combaero.network.mpce_v2_element import MPCEv2Element
from validation.junction import random_robustness as rr

_TRACED_SEED = 20260906


@pytest.fixture(scope="module")
def traced_case():
    """The plateau case the investigation took apart: joining flow, branch at
    156.3 deg, area ratio 2.565, driven by a flow and two pressures."""
    rng = random.Random(_TRACED_SEED)
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


def _solve(case, alpha: float):
    original = MPCEv2Element.soft_penalty_alpha
    MPCEv2Element.soft_penalty_alpha = alpha
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return NetworkSolver(rr.build(case)).solve(timeout=30.0)
    finally:
        MPCEv2Element.soft_penalty_alpha = original


# ---------------------------------------------------------------------------
# The mechanism
# ---------------------------------------------------------------------------


def test_the_penalty_shares_a_row_with_the_continuity_relation():
    """The structural cause. If the penalty ever gets a row of its own this
    test should be the one that fails, because the fixed point goes with it."""
    element = MPCEv2Element(
        id="jct",
        inlet_nodes=["a", "b"],
        outlet_nodes=["c"],
        inlet_angles_deg=[0.0, 90.0],
        outlet_angles_deg=[0.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="merge",
        strict=False,
    )
    assert element.N == 3
    # N ports + one mass row: no spare row for a penalty to live in.
    assert len(element.port_nodes) + 1 == 4


def test_the_barrier_fixed_point_follows_the_closed_form(traced_case):
    """slack* = sqrt(dP / alpha). Measured against the formula rather than
    against a recorded number, so it pins the mechanism and not an output."""
    sol = _solve(traced_case, 1.0e7)

    assert not sol["__success__"], "the old weight is expected to fail here"
    names = list(sol["__unknown_names__"])
    x = np.array(sol["__x_solution__"])

    # Two quantities measured independently of each other: the slack comes
    # from a mass-flow unknown, dP from two pressure unknowns. The claim is
    # that the barrier balances one against the other.
    slack = abs(float(x[names.index("lc_str.m_dot")]))
    dP = abs(float(x[names.index("port_str.Pt")]) - float(x[names.index("jct.P_jct")]))
    assert slack > 0.0 and dP > 0.0

    predicted = math.sqrt(dP / 1.0e7)
    assert predicted == pytest.approx(slack, rel=0.25), (
        f"slack {slack:.6f} kg/s against sqrt(dP/alpha) = {predicted:.6f} "
        f"for dP = {dP:.1f} Pa: the barrier is not balancing the pressure error"
    )

    # And the balance point is a large fraction of the flow, which is why the
    # sign flip the barrier is supposed to enable never happens.
    m_com = abs(float(x[names.index("lc_com.m_dot")]))
    assert slack / m_com > 0.3, (
        f"the fixed point sat at {100.0 * slack / m_com:.1f}% of the common flow"
    )


def test_an_in_regime_root_existed_all_along(traced_case):
    """The floor was not the model running out of solutions. Seeded at the
    operating point the reduced system predicts, the same network converges."""
    m_com, q = rr.operating_point(traced_case)
    solver = NetworkSolver(rr.build(traced_case))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        seed = solver.solve(timeout=30.0)
    names = list(seed["__unknown_names__"])
    x0 = np.array(seed["__x_solution__"])
    for name, value in (
        ("lc_com.m_dot", m_com),
        ("lc_str.m_dot", m_com * (1.0 - q)),
        ("lc_bra.m_dot", m_com * q),
    ):
        x0[names.index(name)] = value

    original = MPCEv2Element.soft_penalty_alpha
    MPCEv2Element.soft_penalty_alpha = 1.0e7
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solver.solve(timeout=30.0, x0=x0)
    finally:
        MPCEv2Element.soft_penalty_alpha = original

    assert sol["__success__"], "the in-regime root is not reachable even when seeded"
    assert sol["__final_norm__"] < 1e-3


# ---------------------------------------------------------------------------
# The fix
# ---------------------------------------------------------------------------


def test_the_traced_case_converges_at_the_shipped_weight(traced_case):
    sol = _solve(traced_case, MPCEv2Element.soft_penalty_alpha)

    assert sol["__success__"], f"|F| = {sol['__final_norm__']:.4e}"
    assert sol["__final_norm__"] < 1e-3


def test_it_converges_into_the_declared_regime(traced_case):
    """Not merely to something. The element is declared a merge, and the
    solution must have both inlets flowing in -- otherwise the barrier has
    been traded for a reversed root, which is a different answer."""
    sol = _solve(traced_case, MPCEv2Element.soft_penalty_alpha)
    names = list(sol["__unknown_names__"])
    x = np.array(sol["__x_solution__"])

    assert sol["__success__"]
    assert float(x[names.index("lc_str.m_dot")]) > 0.0
    assert float(x[names.index("lc_bra.m_dot")]) > 0.0


def test_the_weight_is_on_the_saturation_plateau():
    """Guards against the constant being a tuned optimum: the response must be
    flat above it, so neighbouring decades do as well."""
    assert MPCEv2Element.soft_penalty_alpha >= 1.0e9


@pytest.mark.parametrize("alpha", [1.0e9, 1.0e10, 1.0e11, 1.0e12])
def test_neighbouring_decades_all_solve_the_traced_case(traced_case, alpha):
    assert _solve(traced_case, alpha)["__success__"]


# ---------------------------------------------------------------------------
# No regression in aggregate
# ---------------------------------------------------------------------------


def test_the_random_harness_does_not_regress():
    """Banded, because it is a convergence rate over random draws and not a
    deterministic invariant. Baseline at alpha = 1e7 was 98.3%."""
    summary = rr.run(n=150, seed=_TRACED_SEED)
    text = rr.format_summary(summary)

    import re

    match = re.search(r"within Mach 0\.3\s+\d+\s+([\d.]+)%", text)
    assert match, "the harness no longer reports the in-range convergence rate"
    assert float(match.group(1)) >= 97.0, text
