"""The soft barrier had a fixed point, and solves parked in it.

`MultiPortChamberElement` refuses a port flowing against its declared direction. In
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

import combaero as cb
from combaero.network import NetworkSolver
from combaero.network.components import NetworkMixtureState
from combaero.network.mpce_element import MultiPortChamberElement
from validation.junction import random_robustness as rr

_Y = list(cb.species.dry_air_mass())

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
    original = MultiPortChamberElement.soft_penalty_alpha
    MultiPortChamberElement.soft_penalty_alpha = alpha
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return NetworkSolver(rr.build(case)).solve(timeout=30.0)
    finally:
        MultiPortChamberElement.soft_penalty_alpha = original


# ---------------------------------------------------------------------------
# The mechanism
# ---------------------------------------------------------------------------


def test_the_penalty_shares_a_row_with_the_continuity_relation():
    """The structural cause. If the penalty ever gets a row of its own this
    test should be the one that fails, because the fixed point goes with it."""
    element = MultiPortChamberElement(
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


def test_the_barrier_fixed_point_follows_the_closed_form():
    """slack* = sqrt(dP / alpha). Measured against the formula rather than
    against a recorded number, so it pins the mechanism and not an output.

    Evaluated on the element directly. This used to force the traced case with
    the old weight, assert the solve failed, and read the fixed point off the
    stalled iterate. That no longer works: with the compressible stagnation
    closure at the port momentum chambers (#357) the traced case converges at
    every alpha from 1e4 to 1e9 and never goes wrong-direction, so there is no
    stalled iterate to inspect and the barrier term is identically zero -- the
    quantity the old test read back was an ordinary port flow, not a slack.

    The barrier itself is still live (invoked 32 times across 40 sampled
    junctions), so the mechanism is worth pinning; it just has to be pinned
    where it is deterministic rather than where a solve happens to stall.
    """
    element = MultiPortChamberElement(
        id="jct",
        inlet_nodes=["a", "b"],
        outlet_nodes=["c"],
        inlet_angles_deg=[0.0, 90.0],
        outlet_angles_deg=[0.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="merge",
        strict=False,
    )
    element.soft_penalty_alpha = 1.0e7
    alpha = element.effective_penalty_alpha()
    assert alpha == pytest.approx(1.0e7)

    Pt_jct = 100_300.0
    dP = 2_500.0
    slack_star = math.sqrt(dP / alpha)

    # Port 0 is an inlet for a merge, so its canonical sign wants m_dot < 0.
    # Put it wrong-direction with exactly the predicted slack, and place its
    # Pt the predicted pressure error below the junction.
    e0 = float(element._port_signs[0])
    port_mdots = [-e0 * slack_star, -0.04, 0.10]
    states = [
        NetworkMixtureState(P=98_000.0, Pt=Pt_jct - dP, T=300.0, Tt=300.5, m_dot=0.0, Y=list(_Y)),
        NetworkMixtureState(P=98_000.0, Pt=Pt_jct, T=300.0, Tt=300.5, m_dot=0.0, Y=list(_Y)),
        NetworkMixtureState(P=98_000.0, Pt=Pt_jct, T=300.0, Tt=300.5, m_dot=0.0, Y=list(_Y)),
    ]

    residuals, _ = element._soft_barrier_residual(states, Pt_jct, port_mdots)

    # R_0 = (Pt_0 - Pt_jct) + alpha * slack^2 = -dP + alpha * (dP/alpha) = 0.
    assert residuals[0] == pytest.approx(0.0, abs=1e-6 * dP), (
        f"slack* = sqrt(dP/alpha) = {slack_star:.6f} kg/s did not zero the "
        f"barrier row for dP = {dP:.1f} Pa: got {residuals[0]:.6e}"
    )

    # And it is genuinely the penalty doing it: halve the slack and the row
    # must no longer balance.
    off, _ = element._soft_barrier_residual(states, Pt_jct, [-e0 * 0.5 * slack_star, -0.04, 0.10])
    assert abs(off[0]) > 0.5 * dP, (
        "the row balanced without the penalty carrying it; the fixed point is "
        "not the mechanism this test claims"
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

    original = MultiPortChamberElement.soft_penalty_alpha
    MultiPortChamberElement.soft_penalty_alpha = 1.0e7
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solver.solve(timeout=30.0, x0=x0)
    finally:
        MultiPortChamberElement.soft_penalty_alpha = original

    assert sol["__success__"], "the in-regime root is not reachable even when seeded"
    assert sol["__final_norm__"] < 1e-3


# ---------------------------------------------------------------------------
# The fix
# ---------------------------------------------------------------------------


def test_the_traced_case_converges_at_the_shipped_weight(traced_case):
    sol = _solve(traced_case, MultiPortChamberElement.soft_penalty_alpha)

    assert sol["__success__"], f"|F| = {sol['__final_norm__']:.4e}"
    assert sol["__final_norm__"] < 1e-3


def test_it_converges_into_the_declared_regime(traced_case):
    """Not merely to something. The element is declared a merge, and the
    solution must have both inlets flowing in -- otherwise the barrier has
    been traded for a reversed root, which is a different answer."""
    sol = _solve(traced_case, MultiPortChamberElement.soft_penalty_alpha)
    names = list(sol["__unknown_names__"])
    x = np.array(sol["__x_solution__"])

    assert sol["__success__"]
    assert float(x[names.index("lc_str.m_dot")]) > 0.0
    assert float(x[names.index("lc_bra.m_dot")]) > 0.0


def test_the_weight_is_on_the_saturation_plateau():
    """Guards against the constant being a tuned optimum: the response must be
    flat above it, so neighbouring decades do as well."""
    assert MultiPortChamberElement.soft_penalty_alpha >= 1.0e9


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
