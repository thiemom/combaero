"""Junction convergence on RANDOM physical boundary conditions.

The network scorecard measures convergence on cases whose boundary conditions
are built from Bassett's analytical K at a target split. For an ACCURACY
question that is exactly right. For a ROBUSTNESS question it is circular: it
asks how often the solver reaches an operating point the paper's own
correlation predicts, on the same 2073 points the closure has been measured
and re-measured against.

This module asks the production question instead. It samples a junction and a
set of boundary conditions uniformly inside physically sensible ranges, with no
reference to any paper or dataset, and asks whether the solver returns an
admissible answer.

**The width of the sample space is the point.** A narrow space could be
flattered by tuning; this one cannot, because nothing in the closure is fitted
against a distribution the closure never sees. `test_junction_random_robustness`
pins the ranges for that reason: narrowing them to improve the number requires
editing an assertion, which shows up in a diff.

Run the full sweep with::

    uv run python -m validation.junction.random_robustness 2000

Outcome classes are kept apart on purpose. Only two of them are the solver's
fault:

``converged``       an admissible root
``no root exists``  the drawn boundary conditions constrain a function of the
                    split to a value the closure cannot produce, so there is
                    nothing to converge to. Counting these as failures is the
                    mistake that started the operating-point investigation
                    (issue #271, docs/archive/JUNCTION_OPERATING_POINT_271.md)
``rejected``        converged, then demoted by the junction's own physics
                    checks: a wrong-direction port flow, or a state where the
                    closure is not net dissipative
``no progress``     the solver could not move although a root exists
``raised``          an exception, which should never happen
"""

from __future__ import annotations

import math
import random
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

import combaero as cb
from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network._mynard2010 import junction_loss_coefficient
from combaero.network.mpce_v2_element import MPCEv2Element

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_X = list(cb.mass_to_mole(_Y))

Drive = Literal["imposed_flows", "flow_and_pressures", "all_pressures"]
DRIVES: tuple[Drive, ...] = ("imposed_flows", "flow_and_pressures", "all_pressures")

# ---------------------------------------------------------------------------
# The sample space. Widen it freely; narrowing it needs a reason in the diff.
# ---------------------------------------------------------------------------
PT_RANGE_PA = (0.5e5, 10.0e5)
TT_RANGE_K = (250.0, 600.0)
#: The closure is documented for low Mach; 0.25 is past where it is trusted,
#: on purpose, so the sweep reports what happens at the edge.
MACH_RANGE = (0.005, 0.25)
THETA_RANGE_DEG = (15.0, 165.0)
#: Area ratio A_common / A_branch, sampled log-uniformly: a branch both much
#: larger and much smaller than the main duct.
PSI_RANGE = (0.3, 10.0)
#: Common-port area, 1 cm2 to 500 cm2, log-uniform.
AREA_RANGE_M2 = (1.0e-4, 0.05)
SPLIT_RANGE = (0.05, 0.95)
#: Imposed pressure differences, as multiples of the common dynamic head. The
#: lower bound is negative because both Bassett and Hager measure a mildly
#: negative straight-leg coefficient in dividing flow.
K_STRAIGHT_RANGE = (-0.5, 3.0)
K_BRANCH_RANGE = (-0.5, 5.0)

_Q_GRID = np.linspace(0.02, 0.98, 97)


@dataclass(frozen=True)
class Case:
    """One randomly drawn junction and its boundary conditions."""

    Pt: float
    Tt: float
    mach: float
    theta_deg: float
    psi: float
    split: float
    area: float
    drive: Drive
    joining: bool
    k_straight: float
    k_branch: float


def _log_uniform(rng: random.Random, lo: float, hi: float) -> float:
    return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def sample(rng: random.Random) -> Case:
    """Draw one case uniformly inside the ranges above."""
    return Case(
        Pt=rng.uniform(*PT_RANGE_PA),
        Tt=rng.uniform(*TT_RANGE_K),
        mach=rng.uniform(*MACH_RANGE),
        theta_deg=rng.uniform(*THETA_RANGE_DEG),
        psi=_log_uniform(rng, *PSI_RANGE),
        split=rng.uniform(*SPLIT_RANGE),
        area=_log_uniform(rng, *AREA_RANGE_M2),
        drive=rng.choice(DRIVES),
        joining=rng.random() < 0.5,
        k_straight=rng.uniform(*K_STRAIGHT_RANGE),
        k_branch=rng.uniform(*K_BRANCH_RANGE),
    )


def _scales(case: Case) -> tuple[float, float]:
    """Common mass flow and dynamic head for the drawn state."""
    rho = float(cb.density(case.Tt, case.Pt, _X))
    a = float(cb.speed_of_sound(case.Tt, _X))
    m_com = rho * a * case.mach * case.area
    q_dyn = 0.5 * rho * (m_com / (rho * case.area)) ** 2
    return m_com, q_dyn


def build(case: Case) -> FlowNetwork:
    """The three-port network for this case, wired as production wires one."""
    m_com, q_dyn = _scales(case)
    a_bra = case.area / case.psi
    # Joining flow reverses which side of the junction the loss sits on.
    sign = -1.0 if case.joining else 1.0
    Pt_str = case.Pt - sign * case.k_straight * q_dyn
    Pt_bra = case.Pt - sign * case.k_branch * q_dyn

    net = FlowNetwork()
    for pid, area in (
        ("port_com", case.area),
        ("port_str", case.area),
        ("port_bra", a_bra),
    ):
        net.add_node(MomentumChamberNode(pid, area=area))

    if case.joining:
        if case.drive == "imposed_flows":
            net.add_node(
                MassFlowBoundary(
                    "b_str", m_dot=(1.0 - case.split) * m_com, Tt=case.Tt, Y=_Y
                )
            )
            net.add_node(
                MassFlowBoundary("b_bra", m_dot=case.split * m_com, Tt=case.Tt, Y=_Y)
            )
            net.add_node(PressureBoundary("b_com", Pt=case.Pt, Tt=case.Tt, Y=_Y))
        else:
            net.add_node(PressureBoundary("b_str", Pt=Pt_str, Tt=case.Tt, Y=_Y))
            net.add_node(PressureBoundary("b_bra", Pt=Pt_bra, Tt=case.Tt, Y=_Y))
            if case.drive == "flow_and_pressures":
                net.add_node(MassFlowBoundary("b_com", m_dot=m_com, Tt=case.Tt, Y=_Y))
            else:
                net.add_node(PressureBoundary("b_com", Pt=case.Pt, Tt=case.Tt, Y=_Y))
        net.add_element(LosslessConnectionElement("lc_str", "b_str", "port_str"))
        net.add_element(LosslessConnectionElement("lc_bra", "b_bra", "port_bra"))
        net.add_element(LosslessConnectionElement("lc_com", "port_com", "b_com"))
        junction = MPCEv2Element(
            id="jct",
            inlet_nodes=["port_str", "port_bra"],
            outlet_nodes=["port_com"],
            inlet_angles_deg=[0.0, case.theta_deg],
            outlet_angles_deg=[0.0],
            port_areas=[case.area, a_bra, case.area],
            flow_direction="merge",
            strict=False,
        )
    else:
        if case.drive == "imposed_flows":
            net.add_node(MassFlowBoundary("b_com", m_dot=m_com, Tt=case.Tt, Y=_Y))
            net.add_node(
                MassFlowBoundary("b_bra", m_dot=case.split * m_com, Tt=case.Tt, Y=_Y)
            )
            net.add_node(PressureBoundary("b_str", Pt=Pt_str, Tt=case.Tt, Y=_Y))
        else:
            net.add_node(PressureBoundary("b_str", Pt=Pt_str, Tt=case.Tt, Y=_Y))
            net.add_node(PressureBoundary("b_bra", Pt=Pt_bra, Tt=case.Tt, Y=_Y))
            if case.drive == "flow_and_pressures":
                net.add_node(MassFlowBoundary("b_com", m_dot=m_com, Tt=case.Tt, Y=_Y))
            else:
                net.add_node(PressureBoundary("b_com", Pt=case.Pt, Tt=case.Tt, Y=_Y))
        net.add_element(LosslessConnectionElement("lc_com", "b_com", "port_com"))
        net.add_element(LosslessConnectionElement("lc_str", "port_str", "b_str"))
        net.add_element(LosslessConnectionElement("lc_bra", "port_bra", "b_bra"))
        junction = MPCEv2Element(
            id="jct",
            inlet_nodes=["port_com"],
            outlet_nodes=["port_str", "port_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, case.theta_deg],
            port_areas=[case.area, case.area, a_bra],
            flow_direction="branch",
            strict=False,
        )
    net.add_element(junction)
    return net


def _closure_curve(case: Case) -> tuple[np.ndarray, np.ndarray]:
    """The closure's own straight and branch coefficients over the split."""
    a_bra = case.area / case.psi
    areas = np.array([case.area, case.area, a_bra])
    straight, branch = [], []
    for q in _Q_GRID:
        if case.joining:
            U = np.array([(1.0 - q) * 10.0, q * 10.0 * case.psi, -10.0])
            angles = np.array([0.0, math.radians(case.theta_deg), math.pi])
        else:
            U = np.array([10.0, -(1.0 - q) * 10.0, -q * 10.0 * case.psi])
            angles = np.array([math.pi, 0.0, math.radians(case.theta_deg)])
        result = junction_loss_coefficient(U, areas, angles)
        if result.K is None or len(result.K) != 2:
            straight.append(np.nan)
            branch.append(np.nan)
        else:
            straight.append(float(result.K[0]))
            branch.append(float(result.K[1]))
    return np.array(straight), np.array(branch)


def has_root(case: Case) -> bool | None:
    """Whether the drawn boundary conditions admit a solution at all.

    With both flows imposed the split is a boundary condition and a root always
    exists. The other two drives constrain a FUNCTION of the split: the
    DIFFERENCE of the two coefficients when the flow level is fixed by a mass
    boundary, their RATIO when the level is free too. A draw whose target lies
    outside what the closure can produce has no root, and blaming the solver
    for it would be a category error.

    None means undetermined, which happens when the closure returns nothing
    across the whole grid.
    """
    if case.drive == "imposed_flows":
        return True
    straight, branch = _closure_curve(case)
    ok = ~(np.isnan(straight) | np.isnan(branch))
    if not ok.any():
        return None
    straight, branch = straight[ok], branch[ok]
    if case.drive == "flow_and_pressures":
        target = (case.k_branch - case.k_straight) * (-1.0 if case.joining else 1.0)
        spread = straight - branch if case.joining else branch - straight
        return bool(spread.min() <= target <= spread.max())
    if abs(case.k_branch) < 1e-9:
        return None
    target = case.k_straight / case.k_branch
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = straight / branch
    crossings = np.where(np.diff(np.sign(ratio - target)) != 0)[0]
    if crossings.size == 0:
        return False
    # Matching the ratio is necessary but NOT sufficient. The level that
    # crossing implies is q_dyn = dP_branch / K_lat, and a negative q_dyn has
    # no real mass flow behind it, so such a crossing is not a root. Leaving
    # this out made the test optimistic: it called 11 of 60 draws solvable
    # that have no solution, and the solver was then blamed for failing on
    # them. The same condition is what `predicted_root` applies.
    return bool(any(case.k_branch * branch[i] > 0.0 for i in crossings))


def classify(net: FlowNetwork, timeout: float = 20.0) -> str:
    """Solve and bucket the outcome. Never raises."""
    solver = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            sol = solver.solve(timeout=timeout)
        except Exception:  # noqa: BLE001 -- the bucket IS the report
            return "raised"
    if sol.get("__success__"):
        return "converged"
    message = (sol.get("__message__") or "").lower()
    if "unphysical" in message or "dissipative" in message:
        return "rejected"
    if "not making good progress" in message:
        return "no progress"
    if "exceeds" in message:
        return "residual too large"
    return "other"


@dataclass
class Summary:
    n: int
    seed: int
    outcomes: Counter = field(default_factory=Counter)
    by_cut: dict[tuple[str, str], Counter] = field(
        default_factory=lambda: defaultdict(Counter)
    )

    @property
    def n_with_root(self) -> int:
        return sum(self.by_cut[("root", "exists")].values())

    @property
    def converged_share_where_solvable(self) -> float:
        """The headline: of the draws that admit a root, how many are found."""
        n = self.n_with_root
        return self.by_cut[("root", "exists")]["converged"] / n if n else float("nan")


def _bucket(name: str, value: float, edges: tuple[float, float]) -> str:
    lo, hi = edges
    return f"{name}<{lo:g}" if value < lo else (f"{name}{lo:g}-{hi:g}" if value < hi else f"{name}>{hi:g}")


def run(n: int = 2000, seed: int = 20260906) -> Summary:
    """Draw and solve ``n`` random cases. Deterministic for a given seed."""
    rng = random.Random(seed)
    summary = Summary(n=n, seed=seed)
    for _ in range(n):
        case = sample(rng)
        root = has_root(case)
        outcome = classify(build(case))
        if root is False and outcome in ("no progress", "residual too large", "other"):
            outcome = "no root exists"
        summary.outcomes[outcome] += 1
        label = {True: "exists", False: "none", None: "undetermined"}[root]
        summary.by_cut[("root", label)][outcome] += 1
        if root is True:
            summary.by_cut[("solvable drive", case.drive)][outcome] += 1
        summary.by_cut[("flow", "joining" if case.joining else "dividing")][outcome] += 1
        summary.by_cut[("mach", _bucket("", case.mach, (0.05, 0.15)))][outcome] += 1
        summary.by_cut[("area ratio", _bucket("", case.psi, (1.0, 3.0)))][outcome] += 1
        summary.by_cut[("angle", _bucket("", case.theta_deg, (60.0, 120.0)))][outcome] += 1
    return summary


_CLASSES = (
    "converged",
    "no root exists",
    "rejected",
    "no progress",
    "residual too large",
    "other",
    "raised",
)


def format_summary(summary: Summary) -> str:
    lines = [
        f"Random physical boundary conditions, n = {summary.n}, seed {summary.seed}",
        "",
        f"  {'outcome':<20} {'n':>6} {'share':>7}",
    ]
    for name in _CLASSES:
        count = summary.outcomes[name]
        if count:
            lines.append(f"  {name:<20} {count:>6} {100 * count / summary.n:>6.1f}%")
    lines.append(f"  {'TOTAL':<20} {summary.n:>6}")
    lines.append("")
    lines.append(
        f"  Of the {summary.n_with_root} draws that admit a root, "
        f"{100 * summary.converged_share_where_solvable:.1f}% converge to an admissible one."
    )
    lines.append("")
    lines.append(f"  {'cut':<16} {'bucket':<22} {'n':>5} {'converged':>10} {'no prog':>9}")
    for (cut, bucket), counts in sorted(summary.by_cut.items()):
        n = sum(counts.values())
        if not n:
            continue
        lines.append(
            f"  {cut:<16} {bucket:<22} {n:>5} {100 * counts['converged'] / n:>9.1f}% "
            f"{100 * counts['no progress'] / n:>8.1f}%"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    import sys

    args: list[Any] = list(sys.argv[1:] if argv is None else argv)
    n = int(args[0]) if args else 2000
    seed = int(args[1]) if len(args) > 1 else 20260906
    print(format_summary(run(n=n, seed=seed)))


if __name__ == "__main__":
    main()
