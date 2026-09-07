"""The Unified0D closure scored against Bassett's measured points DIRECTLY.

No solver, no network, no imposed pressures: `junction_loss_coefficient` is
called at each digitised point's own split and its K read off. Everything the
network path adds -- the boundary skeleton, the operating point the solve
settles at, the guards, the barrier -- is absent by construction.

That separation is the point. The junction scorecard measures a model AND a
coupling together, so a number from it cannot say which one carries the error.
Two questions need it separated:

1. **The 0.2405 comparison.** Issue #272 sets MPCE-v1's Bassett separating
   MAE, 0.2405 over 287/315 converged points, as the bar a candidate must
   beat. That number came through a network solve. This measures the closure
   alone against the same digitised points, so a difference between the two
   locates the error in the coupling rather than the closure, or the reverse.

2. **The theta-dependence question.** Bassett section 3.2.1 states, and Fig 7c
   confirms, that `K5` -- the straight-leg separating coefficient -- is
   independent of `theta` and `psi`. That is a statement about the physics
   which any closure must reproduce, and it is checkable without a solver at
   all: hold the split, sweep the branch angle, and see whether K_straight
   moves.

These tests report and pin; they do not assert the closure is good. Where it
disagrees with Bassett the number is recorded as an xfail-free measurement so
a later change can be scored against it.

Axis convention, established in #295 and #301: Bassett indexes each
coefficient on the mass-flow fraction in ITS OWN leg, so a K5 file's abscissa
is the straight fraction and a K6 file's is the lateral one. Getting this
wrong mirrors the curve and was a real defect twice.

Issue #272.
"""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
import pytest

from combaero.network._mynard2010 import junction_loss_coefficient
from validation.junction.network_runner import _LATERAL_K_IDS
from validation.junction.runner import _read_xy_csv
from validation.junction.schema import load_dataset

# Separating flow: one supplier on the common leg, two collectors.
_SEPARATING_K_IDS = {"K5", "K6"}
_U_REF = 10.0


def _closed_form(q_lateral: float, psi: float, theta_rad: float) -> tuple[float, float]:
    """(K_straight, K_lateral) from the raw closure, at an imposed split.

    Port 0 is the common inlet at pi (Mynard's axial-back convention), port 1
    the straight outlet at 0, port 2 the lateral at theta. Velocities are set
    from the split directly, so the operating point IS the abscissa -- there is
    no solve to drift away from it.
    """
    area = 0.01
    areas = np.array([area, area, area / psi])
    U = np.array(
        [
            _U_REF,
            -(1.0 - q_lateral) * _U_REF,
            -q_lateral * _U_REF * psi,
        ]
    )
    angles = np.array([math.pi, 0.0, theta_rad])
    result = junction_loss_coefficient(U, areas, angles)
    if result.K is None or len(np.atleast_1d(result.K)) != 2:
        return float("nan"), float("nan")
    K = np.atleast_1d(result.K)
    return float(K[0]), float(K[1])


@pytest.fixture(scope="module")
def measured():
    """Bassett's measured separating points: (K_id, q_lateral, psi, theta, K)."""
    points = []
    for f in load_dataset().files:
        kid = f.K_id or f.coefficient or ""
        if kid not in _SEPARATING_K_IDS or f.kind != "measured":
            continue
        for x, y in _read_xy_csv(f.path):
            if not 0.02 < x < 0.98:
                continue
            # Each file's abscissa is the fraction in its own leg.
            q_lateral = x if kid in _LATERAL_K_IDS else 1.0 - x
            points.append((kid, q_lateral, f.psi or 1.0, f.theta_deg or 45.0, y))
    assert points, "the Bassett separating measured files have gone missing"
    return points


def _errors(measured) -> dict[str, list[float]]:
    out: dict[str, list[float]] = defaultdict(list)
    for kid, q_lat, psi, theta_deg, value in measured:
        k_str, k_lat = _closed_form(q_lat, psi, math.radians(theta_deg))
        predicted = k_lat if kid in _LATERAL_K_IDS else k_str
        if math.isnan(predicted):
            continue
        out[kid].append(abs(predicted - value))
    return out


# ---------------------------------------------------------------------------
# The closure evaluates at all, everywhere the data is
# ---------------------------------------------------------------------------


def test_the_closure_returns_a_value_at_every_measured_point(measured):
    """No solver means no convergence failures: a closed-form evaluation
    either produces a K or the closure is degenerate there. Coverage of the
    data is therefore 100% or something is wrong, which is itself the first
    thing the comparison with a solver-coupled number needs."""
    missing = [
        (kid, q, psi, theta)
        for kid, q, psi, theta in ((k, q, p, t) for k, q, p, t, _ in measured)
        if math.isnan(_closed_form(q, psi, math.radians(theta))[0])
    ]

    assert not missing, f"{len(missing)} of {len(measured)} points produced no K"


def test_the_dataset_covers_more_than_one_angle_and_area_ratio(measured):
    """Guards the theta test below: a sweep over a single angle would pass it
    vacuously."""
    assert len({t for _, _, _, t, _ in measured}) >= 2
    assert len({p for _, _, p, _, _ in measured}) >= 2


# ---------------------------------------------------------------------------
# (1) The 0.2405 comparison
# ---------------------------------------------------------------------------


def test_the_closed_form_matches_the_acceptance_gate_cells(measured):
    """The comparison issue #272 asks for, cell by cell.

    The gate's 315 records are these same 105 measured K5/K6 points expanded
    over three network topologies, so the coefficients and the (theta, psi)
    cells line up exactly with the table on the issue. What differs is the
    coupling: the gate ran MPCE-v1 through a network solve, this runs the
    Unified0D closure alone.

    Two differences at once -- model AND coupling -- so this test also scores
    the CURRENT model through the network on the same points, which separates
    them. Printed rather than asserted: the purpose is to locate the error.
    """
    import warnings

    from validation.junction.models.mpce_v2_network import MPCEv2Network
    from validation.junction.network_runner import iter_network_records

    # Closed form, per cell.
    closed: dict[tuple[str, float, float], list[float]] = defaultdict(list)
    for kid, q_lat, psi, theta, value in measured:
        k_str, k_lat = _closed_form(q_lat, psi, math.radians(theta))
        predicted = k_lat if kid in _LATERAL_K_IDS else k_str
        if not math.isnan(predicted):
            closed[(kid, theta, psi)].append(abs(predicted - value))

    # The same model through the network, same points, imposed-q only.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        records = [
            r
            for r in iter_network_records(
                MPCEv2Network(strict=False), load_dataset(), topologies=("imposed_q",)
            )
            if r.K_id in _SEPARATING_K_IDS and r.error is not None
        ]
    network: dict[tuple[str, float, float], list[float]] = defaultdict(list)
    for r in records:
        network[(r.K_id, r.theta_deg, r.psi)].append(abs(r.error))

    print(
        f"\n  {'K':<4} {'theta':>6} {'psi':>5} | {'v1 network':>11} "
        f"| {'v2 network':>11} {'n':>4} | {'v2 closed':>10} {'n':>4}"
    )
    gate = {
        ("K5", 45.0, 1.0): 0.2604,
        ("K5", 45.0, 3.0): 0.4788,
        ("K5", 60.0, 1.0): 0.2629,
        ("K5", 90.0, 1.0): 0.3397,
        ("K6", 45.0, 1.0): 0.1771,
        ("K6", 45.0, 3.0): 0.1907,
        ("K6", 60.0, 1.0): 0.1881,
        ("K6", 90.0, 1.0): 0.0979,
        ("K6", 120.0, 1.0): 0.6193,
    }
    for cell in sorted(set(closed) | set(network) | set(gate)):
        kid, theta, psi = cell
        v1 = gate.get(cell)
        net = network.get(cell)
        cf = closed.get(cell)
        print(
            f"  {kid:<4} {theta:>6.0f} {psi:>5.1f} | "
            f"{(f'{v1:.4f}' if v1 is not None else '-'):>11} | "
            f"{(f'{np.mean(net):.4f}' if net else '-'):>11} "
            f"{(len(net) if net else 0):>4} | "
            f"{(f'{np.mean(cf):.4f}' if cf else '-'):>10} "
            f"{(len(cf) if cf else 0):>4}"
        )

    pooled_closed = [e for v in closed.values() for e in v]
    pooled_net = [e for v in network.values() for e in v]
    print("\n  pooled  v1 network (issue #272): 0.2405 over 287/315 records")
    print(
        f"          v2 network, imposed_q   : {np.mean(pooled_net):.4f} "
        f"over {len(pooled_net)} points"
    )
    print(
        f"          v2 closed form          : {np.mean(pooled_closed):.4f} "
        f"over {len(pooled_closed)} points"
    )

    assert pooled_closed and pooled_net, "one of the two paths scored nothing"


def test_the_closed_form_separating_mae_is_recorded(measured):
    """The headline number, uncoupled. Banded rather than pinned to a digit:
    it is a mean over digitised points and moves with any closure change, and
    the purpose is to locate the error, not to freeze it.

    The band is deliberately wide. A tight assertion here would be a fit to
    the validation set, which is the failure this repo names explicitly.
    """
    errors = _errors(measured)
    pooled = [e for values in errors.values() for e in values]
    mae = float(np.mean(pooled))

    for kid in sorted(errors):
        print(f"  {kid}: n={len(errors[kid]):>3}  MAE {np.mean(errors[kid]):.4f}")
    print(f"  pooled: n={len(pooled)}  MAE {mae:.4f}  (MPCE-v1 network MAE: 0.2405)")

    assert 0.0 < mae < 2.0, f"closed-form separating MAE {mae:.4f} is out of any plausible range"


def test_the_lateral_and_straight_legs_are_scored_separately(measured):
    """Pooling them hides which leg carries the error, and the two are known
    to behave differently -- the straight leg is the one Bassett says is
    independent of geometry."""
    errors = _errors(measured)

    assert set(errors) == _SEPARATING_K_IDS, f"only scored {set(errors)}"
    for kid in _SEPARATING_K_IDS:
        assert len(errors[kid]) >= 20, f"{kid} has only {len(errors[kid])} scored points"


# ---------------------------------------------------------------------------
# (2) The theta-dependence question
# ---------------------------------------------------------------------------


def test_k_straight_theta_dependence_is_measured_not_assumed():
    """Bassett 3.2.1 and Fig 7c: K5 is independent of theta and psi.

    Checked directly on the closure, with no data and no solver: fix the
    split, sweep the branch angle, and see how far K_straight moves. Reported
    as a spread so a later change can be scored against it.
    """
    spreads = []
    for q in (0.2, 0.4, 0.6, 0.8):
        for psi in (1.0, 2.0):
            values = [
                _closed_form(q, psi, math.radians(theta))[0]
                for theta in (15.0, 30.0, 45.0, 60.0, 75.0, 90.0)
            ]
            assert not any(math.isnan(v) for v in values)
            spread = max(values) - min(values)
            spreads.append(spread)
            print(
                f"  q={q:.1f} psi={psi:.1f}: K_straight {min(values):+.4f}..{max(values):+.4f}"
                f"  spread {spread:.4f}"
            )

    worst = max(spreads)
    print(f"  worst spread over theta: {worst:.4f}")
    # Not an assertion that the closure is right -- an assertion that the
    # measurement happened and produced a finite, reportable number.
    assert all(math.isfinite(s) for s in spreads)
    assert worst >= 0.0


def test_k_straight_psi_dependence_is_measured_too():
    """The same claim covers psi. Separated from theta so a failure says which."""
    spreads = []
    for q in (0.2, 0.4, 0.6, 0.8):
        for theta in (45.0, 90.0):
            values = [_closed_form(q, psi, math.radians(theta))[0] for psi in (0.5, 1.0, 2.0, 4.0)]
            assert not any(math.isnan(v) for v in values)
            spreads.append(max(values) - min(values))
            print(f"  q={q:.1f} theta={theta:.0f}: K_straight spread over psi {spreads[-1]:.4f}")

    print(f"  worst spread over psi: {max(spreads):.4f}")
    assert all(math.isfinite(s) for s in spreads)


def test_the_measured_k5_points_are_themselves_theta_independent(measured):
    """The other half: Bassett's own K5 data must show the independence the
    text claims, or the claim cannot be held against the model. Grouped by
    split and compared across angles.
    """
    by_q: dict[float, dict[float, list[float]]] = defaultdict(lambda: defaultdict(list))
    for kid, q_lat, _psi, theta, value in measured:
        if kid != "K5":
            continue
        by_q[round(1.0 - q_lat, 1)][theta].append(value)

    spreads = []
    for q, per_theta in sorted(by_q.items()):
        if len(per_theta) < 2:
            continue
        means = [float(np.mean(v)) for v in per_theta.values()]
        spreads.append(max(means) - min(means))
        print(f"  q_straight={q:.1f}: K5 across {len(per_theta)} angles, spread {spreads[-1]:.4f}")

    assert spreads, "no split has K5 measured at more than one angle"
    print(f"  worst measured K5 spread over theta: {max(spreads):.4f}")


def test_the_theta_sweep_is_not_vacuous():
    """The invariance of K_straight is only meaningful if the sweep MOVES
    something. K_lateral must vary strongly over the same angles, and
    K_straight must vary over the split -- otherwise a closure returning a
    constant would pass the test above.
    """
    lateral_spreads = []
    for q in (0.2, 0.5, 0.8):
        values = [_closed_form(q, 1.0, math.radians(theta))[1] for theta in (15.0, 45.0, 90.0)]
        lateral_spreads.append(max(values) - min(values))
    assert min(lateral_spreads) > 0.1, (
        f"K_lateral barely moves over theta ({lateral_spreads}); the sweep proves nothing"
    )

    over_q = [_closed_form(q, 1.0, math.radians(45.0))[0] for q in (0.2, 0.4, 0.6, 0.8)]
    assert max(over_q) - min(over_q) > 0.1, (
        f"K_straight is constant over the split too ({over_q}); it is not being computed"
    )


def test_the_energy_transfer_term_breaks_bassett_s_theta_independence():
    """A structural argument for the production default, found by trying to
    prove the opposite.

    The sweeps above run at `eta_scale = 0.0`, where Mynard's CFD-fitted
    energy-transfer factor (Eq 36) is off, and K_straight comes out EXACTLY
    independent of theta and psi -- Bassett section 3.2.1 and Fig 7c. This
    test was written to confirm that independence is structural rather than an
    artifact of the default. It is not.

    Turning the term on moves K_straight by ~0.23 over the branch angle, which
    is the same order as the coefficient itself. So Mynard's energy-transfer
    factor is INCOMPATIBLE with the constraint #272 requires any correction to
    preserve. `eta_scale = 0.0` was already the production default, retired on
    accuracy grounds against the digitised data; this is an independent reason
    for it, and a constraint on any future attempt to switch it back on.

    Recorded as a measurement, not a target -- the numbers are printed so a
    later change can be scored against them.
    """
    area = 0.01

    def k_straight(q, psi, theta_deg, eta):
        areas = np.array([area, area, area / psi])
        U = np.array([_U_REF, -(1.0 - q) * _U_REF, -q * _U_REF * psi])
        angles = np.array([math.pi, 0.0, math.radians(theta_deg)])
        result = junction_loss_coefficient(U, areas, angles, eta_scale=eta)
        return float(np.atleast_1d(result.K)[0])

    print()
    worst_off = 0.0
    worst_on = 0.0
    for q in (0.2, 0.5, 0.8):
        for label, eta in (("eta 0.0 (shipped)", 0.0), ("eta 1.0 (Mynard)", 1.0)):
            over_theta = [k_straight(q, 1.0, t, eta) for t in (15.0, 45.0, 90.0, 120.0)]
            over_psi = [k_straight(q, p, 45.0, eta) for p in (0.5, 1.0, 2.0, 4.0)]
            s_theta = max(over_theta) - min(over_theta)
            s_psi = max(over_psi) - min(over_psi)
            if eta == 0.0:
                worst_off = max(worst_off, s_theta, s_psi)
            else:
                worst_on = max(worst_on, s_theta, s_psi)
            print(
                f"  q={q:.1f} {label:<18} K_straight spread: theta {s_theta:.6f}  psi {s_psi:.6f}"
            )

    print(f"\n  worst spread with the term OFF: {worst_off:.6f}")
    print(f"  worst spread with the term ON : {worst_on:.6f}")

    assert worst_off < 1e-12, (
        f"the shipped configuration no longer reproduces Bassett 3.2.1: "
        f"K_straight moves by {worst_off:.6f}"
    )
    assert worst_on > 0.05, (
        "the energy-transfer term no longer breaks the independence, so this "
        "constraint on switching it back on has gone -- re-derive before relying on it"
    )
