"""
Network-mode validation runner.

Unlike the correlation-mode runner (which evaluates the raw K math at each
digitized data point), the network runner builds a full FlowNetwork for
each (paper, K, psi, theta, q) case, solves it through NetworkSolver, and
records the extracted K plus convergence diagnostics.

This is the runner that actually compares the K-closure (TeeJunctionElement)
against MPCE-v1 head-to-head in the solver context they would be used in
production.
"""

from __future__ import annotations

import math
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Protocol

from validation.junction.equivalences import canonical_K
from validation.junction.models._network_builder import (
    ALL_TOPOLOGIES,
    NetworkResult,
    Topology,
)
from validation.junction.schema import Dataset, load_dataset


class NetworkModel(Protocol):
    """Junction element that supports network-mode evaluation."""

    name: str
    SUPPORTED_TOPOLOGIES: tuple[Topology, ...]

    def evaluate_network(
        self,
        paper: str,
        K_id: str,
        q: float,
        psi: float | None,
        theta_rad: float | None,
        topology: Topology = "imposed_q",
        **kwargs: float,
    ) -> NetworkResult: ...


@dataclass
class NetworkRecord:
    """One imposed-q case evaluated through a network solver."""

    paper: str
    K_id: str
    canonical_K: str
    psi: float | None
    theta_deg: float | None
    q: float
    K_measured: float
    K_extracted: float | None
    converged: bool
    residual_norm: float
    wall_time_s: float
    which: str  # "lateral" or "straight" -- which extracted K matches K_id
    topology: Topology = "imposed_q"
    message: str = ""
    # The operating point the solve actually reached, on the FILE's q axis.
    # Only `imposed_q` imposes it; elsewhere it is an outcome and the record's
    # K was extracted there, not at `q` (issue #271).
    q_converged: float | None = None
    # The measured curve's K at `q_converged`, which is the value `K_extracted`
    # must be compared against. None when the solve did not converge, or when
    # the operating point it reached lies outside the curve's digitised range.
    K_measured_at_q_converged: float | None = None

    @property
    def q_drift(self) -> float | None:
        """|target q - achieved q|. None when the solve did not converge."""
        if self.q_converged is None:
            return None
        return abs(self.q - self.q_converged)

    @property
    def off_point(self) -> bool:
        """Converged somewhere other than the operating point the dataset asked for.

        Such a record is still scored, at the point it reached. It is simply
        not evidence about `q`, and counting it as such is what item 9b fixes.
        """
        d = self.q_drift
        return d is not None and d > _Q_OFF_POINT

    @property
    def off_curve(self) -> bool:
        """Converged, but at an operating point the digitised curve does not cover."""
        return (
            self.converged
            and self.K_extracted is not None
            and self.K_measured_at_q_converged is None
        )

    @property
    def error(self) -> float | None:
        """Model minus measurement, AT THE OPERATING POINT THE SOLVE REACHED.

        Scoring against `K_measured` (the value at the TARGET q) is what this
        replaces: on 88 of 618 `mfb_two_pb` records the solve sat more than
        0.05 away in q, so that comparison was at a point never visited.
        """
        if self.K_extracted is None or self.K_measured_at_q_converged is None:
            return None
        return self.K_extracted - self.K_measured_at_q_converged


# Hager 1984 names its coefficients xi_t (straight-through) and xi_l (lateral);
# the schema stores them in `coefficient` with K_id=None, and _which_K is fed
# `K_id or coefficient`, so both spellings must be here.
_LATERAL_K_IDS = {"K1", "K6", "K12", "xi_l"}
_STRAIGHT_K_IDS = {"K2", "K5", "K11", "xi_t"}

# Ground truth for network scoring: digitised measurements, and handbook
# tabulations (Idelchik) which are reference data in their own right. Not a
# paper's own analytical curve ("calc") and not an envelope.
_GROUND_TRUTH_KINDS = {"measured", "tabulated"}

# Topologies whose extracted K is a measurement of the MODEL rather than of the
# fixture. `three_pb` is excluded: it imposes every total pressure through
# lossless connections and `_extract_K` normalises by the fixed reference flow,
# so the extracted K equals the imposed target BY CONSTRUCTION. Falsified with
# a deliberately wrong model (eta_scale=3.0, whose imposed_q K at q=0.8 moves
# 2.88 -> 3.07): it still reported 2.7677 against a target of 2.7689. See
# docs/archive/JUNCTION_OPERATING_POINT_271.md and
# python/tests/test_junction_score_at_achieved_q.py, which pins the
# falsification so this exclusion cannot be quietly reverted.
#
# `three_pb` still runs, and its convergence count is scored. Only its K is not.
_K_SCORING_TOPOLOGIES: frozenset[str] = frozenset({"imposed_q", "mfb_two_pb"})

# An achieved q within this of the target did not drift; it is float noise in
# a topology that imposes the split. Matches the tolerance the instrument is
# calibrated to in test_junction_score_at_achieved_q.py.
_Q_SNAP = 1e-6

# Beyond this, the solve probed a different operating point than the dataset
# asked for. Half the median digitised spacing of the two dominant sources
# (Bassett 0.0997, Idelchik 0.1000; Hager is finer at 0.0150), so a record
# past it has moved onto a neighbouring point's territory. Such a record can
# still be SCORED correctly -- on the curve, where it landed -- but it is no
# longer EVIDENCE about the q the dataset asked for, which is why coverage is
# reported separately from accuracy.
_Q_OFF_POINT = 0.05


def _which_K(K_id: str) -> str | None:
    """Which extracted side does this K_id correspond to."""
    if K_id in _LATERAL_K_IDS:
        return "lateral"
    if K_id in _STRAIGHT_K_IDS:
        return "straight"
    return None


def _curve_value_at(rows: list[tuple[float, float]], q: float) -> float | None:
    """The measured curve's K at ``q``, linearly interpolated between points.

    Returns None when ``q`` lies outside the curve's digitised range: the
    record then has no ground truth at the operating point it reached, and
    extrapolating a digitised curve is not measured data (policy sec 0 item 1).

    Interpolating BETWEEN two digitised points of the same curve is still that
    curve. Its cost was measured leave-one-out over all 64 curves (drop an
    interior point, rebuild it from its neighbours): median |error| 0.030,
    p90 0.175, worst 10.75 on a steep Idelchik K12 branch. That spans two gaps
    and so overstates a single-gap interpolation, but it is the floor below
    which a rescored K cannot resolve anything, and it is why
    ``NetworkCell`` reports ``n_off_curve`` rather than silently extrapolating.
    """
    if not rows:
        return None
    xs = [r[0] for r in rows]
    if q < xs[0] or q > xs[-1]:
        return None
    for i in range(len(rows) - 1):
        q_lo, k_lo = rows[i]
        q_hi, k_hi = rows[i + 1]
        if q_lo <= q <= q_hi:
            if q_hi == q_lo:
                return k_lo
            return k_lo + (k_hi - k_lo) * (q - q_lo) / (q_hi - q_lo)
    return rows[-1][1]


def iter_network_records(
    model: NetworkModel,
    dataset: Dataset,
    topologies: tuple[Topology, ...] = ALL_TOPOLOGIES,
) -> Iterator[NetworkRecord]:
    """For each (file, point, topology) triple where the model declares support,
    run the model network and yield a record. Topologies not in
    `model.SUPPORTED_TOPOLOGIES` are skipped silently.
    """
    from validation.junction.runner import _read_xy_csv

    supported = set(getattr(model, "SUPPORTED_TOPOLOGIES", ALL_TOPOLOGIES))
    for file in dataset.files:
        if file.kind not in _GROUND_TRUTH_KINDS or file.x_axis != "q":
            continue
        K_id = file.K_id or file.coefficient or ""
        which = _which_K(K_id)
        if which is None:
            continue
        rows = _read_xy_csv(file.path)
        curve = sorted(rows)
        theta_rad = math.radians(file.theta_deg) if file.theta_deg is not None else None
        paper = file.path.parent.name
        for q_val, K_m in rows:
            if not (-0.05 <= q_val <= 1.05):
                continue
            for topology in topologies:
                if topology not in supported:
                    continue
                result = model.evaluate_network(
                    paper, K_id, q_val, file.psi, theta_rad, topology=topology
                )
                K_ext = (
                    (result.K_lateral if which == "lateral" else result.K_straight)
                    if result.converged
                    else None
                )
                yield NetworkRecord(
                    paper=paper,
                    K_id=K_id,
                    canonical_K=canonical_K(paper, K_id),
                    psi=file.psi,
                    theta_deg=file.theta_deg,
                    q=q_val,
                    K_measured=K_m,
                    K_extracted=K_ext,
                    converged=result.converged,
                    residual_norm=result.residual_norm,
                    wall_time_s=result.wall_time_s,
                    which=which,
                    topology=topology,
                    message=result.message,
                    q_converged=result.q_converged if result.converged else None,
                    K_measured_at_q_converged=(
                        K_m
                        if result.q_converged is None
                        or abs(result.q_converged - q_val) <= _Q_SNAP
                        else _curve_value_at(curve, result.q_converged)
                    )
                    if result.converged
                    else None,
                )


def run_network(
    model: NetworkModel, topologies: tuple[Topology, ...] = ALL_TOPOLOGIES
) -> list[NetworkRecord]:
    """Convenience: load default dataset, return all network records."""
    return list(iter_network_records(model, load_dataset(), topologies=topologies))


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


@dataclass
class NetworkCell:
    canonical_K: str
    psi_bin: str
    theta_bin: str
    topology: Topology
    N: int = 0
    n_converged: int = 0
    pct_converged: float = 0.0
    rmse_meas: float = math.nan  # over converged records only
    bias_meas: float = math.nan
    median_wall_time_ms: float = 0.0
    # Median |target q - achieved q| over converged records. Large values mean
    # the RMSE above compares the model against the paper at an operating point
    # the solve never visited (issue #271).
    median_q_drift: float = math.nan
    # Converged records whose achieved operating point lies outside the
    # digitised curve, so they have no ground truth to be scored against.
    n_off_curve: int = 0
    # Converged records that landed further than `_Q_OFF_POINT` from the q the
    # dataset asked for. They are scored where they landed; this is the
    # COVERAGE loss, which RMSE alone cannot show.
    n_off_point: int = 0


def _median(xs: list[float]) -> float:
    if not xs:
        return 0.0
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def build_network_cells(records: list[NetworkRecord]) -> list[NetworkCell]:
    """Group by (canonical_K, psi, theta, topology) -> NetworkCell."""
    from collections import defaultdict

    groups: dict[tuple[str, float | None, float | None, Topology], list[NetworkRecord]] = (
        defaultdict(list)
    )
    for r in records:
        groups[(r.canonical_K, r.psi, r.theta_deg, r.topology)].append(r)

    topo_order = {t: i for i, t in enumerate(ALL_TOPOLOGIES)}

    def _sort_key(item):
        (cK, psi, theta, topo), _ = item
        return (
            cK,
            psi if psi is not None else -1.0,
            theta if theta is not None else -1.0,
            topo_order.get(topo, 99),
        )

    cells: list[NetworkCell] = []
    for (cK, psi, theta, topo), recs in sorted(groups.items(), key=_sort_key):
        N = len(recs)
        conv = [r for r in recs if r.converged and r.K_extracted is not None]
        n_conv = len(conv)
        if topo in _K_SCORING_TOPOLOGIES:
            errs = [e for r in conv if (e := r.error) is not None]
        else:
            errs = []
        n_off = sum(1 for r in conv if r.off_curve)
        n_off_pt = sum(1 for r in conv if r.off_point)
        rmse = math.sqrt(sum(e * e for e in errs) / len(errs)) if errs else math.nan
        bias = sum(errs) / len(errs) if errs else math.nan
        wt_ms = 1000.0 * _median([r.wall_time_s for r in recs])
        drifts = [d for r in conv if (d := r.q_drift) is not None]
        drift = _median(drifts) if drifts else math.nan
        cells.append(
            NetworkCell(
                canonical_K=cK,
                psi_bin=f"{psi:g}" if psi is not None else "n/a",
                theta_bin=f"{theta:g}" if theta is not None else "n/a",
                topology=topo,
                N=N,
                n_converged=n_conv,
                pct_converged=n_conv / N if N else 0.0,
                rmse_meas=rmse,
                bias_meas=bias,
                median_wall_time_ms=wt_ms,
                median_q_drift=drift,
                n_off_curve=n_off,
                n_off_point=n_off_pt,
            )
        )
    return cells


def format_network_scorecard(model_name: str, cells: list[NetworkCell]) -> str:
    """Render a human-readable network-mode scorecard with per-topology rows."""
    lines = []
    lines.append(f"=== Network scorecard: {model_name} ===")
    lines.append(
        f"{'canonical_K':<18} {'psi':>5} {'theta':>5} {'topology':<12} "
        f"{'N':>4} {'%conv':>6} {'RMSE':>8} {'bias':>8} {'dq':>6} {'off':>6} {'t_ms':>6}"
    )
    lines.append("-" * 95)
    n_total = 0
    n_conv = 0
    for c in cells:
        n_total += c.N
        n_conv += c.n_converged
        _taut = c.topology not in _K_SCORING_TOPOLOGIES
        rmse_str = ("taut" if _taut else "-") if math.isnan(c.rmse_meas) else f"{c.rmse_meas:.3f}"
        bias_str = ("taut" if _taut else "-") if math.isnan(c.bias_meas) else f"{c.bias_meas:+.3f}"
        drift_str = "-" if math.isnan(c.median_q_drift) else f"{c.median_q_drift:.3f}"
        off_str = (
            "-" if not (c.n_off_point or c.n_off_curve) else f"{c.n_off_point}/{c.n_off_curve}"
        )
        lines.append(
            f"{c.canonical_K:<18} {c.psi_bin:>5} {c.theta_bin:>5} {c.topology:<12} "
            f"{c.N:>4} {c.pct_converged * 100:>5.0f}% {rmse_str:>8} {bias_str:>8} "
            f"{drift_str:>6} {off_str:>6} {c.median_wall_time_ms:>5.1f}"
        )
    lines.append("-" * 95)
    overall_conv = n_conv / n_total if n_total else 0.0
    # Per-topology summary
    from collections import defaultdict

    by_topo: dict[Topology, tuple[int, int, int, int]] = defaultdict(lambda: (0, 0, 0, 0))
    for c in cells:
        n, conv, opt, ocv = by_topo[c.topology]
        by_topo[c.topology] = (
            n + c.N,
            conv + c.n_converged,
            opt + c.n_off_point,
            ocv + c.n_off_curve,
        )
    lines.append(f"Overall: {n_conv}/{n_total} converged ({100 * overall_conv:.0f}%)")
    for topo, (n, conv, opt, ocv) in by_topo.items():
        pct = 100 * conv / n if n else 0.0
        extra = ""
        if conv:
            extra = f", {opt} off-point ({100 * opt / conv:.0f}% of converged)"
            if ocv:
                extra += f", {ocv} off-curve"
        lines.append(f"  {topo}: {conv}/{n} ({pct:.0f}%){extra}")
    lines.append("")
    lines.append(
        "RMSE = model minus measurement AT THE OPERATING POINT THE SOLVE "
        "REACHED, taking the"
    )
    lines.append(
        "       measured curve there by linear interpolation between its "
        "digitised points."
    )
    lines.append(
        "dq   = median |target q - achieved q| over converged records. Only "
        "imposed_q imposes"
    )
    lines.append("       the split; elsewhere the operating point is an outcome.")
    lines.append(
        "off  = off-point / off-curve. Off-point: landed further than 0.05 "
        "from the q the"
    )
    lines.append(
        "       dataset asked for (half the median digitised spacing), so it "
        "is scored where"
    )
    lines.append(
        "       it landed but is NOT evidence about that q -- a COVERAGE "
        "loss RMSE cannot"
    )
    lines.append(
        "       show. Off-curve: landed outside the digitised range "
        "entirely, so it has no"
    )
    lines.append("       ground truth and is excluded from RMSE.")
    lines.append(
        "taut = K not scored: three_pb imposes every Pt, so the extracted K is "
        "the target by"
    )
    lines.append(
        "       construction. Its convergence column is scored; its K column "
        "would measure the"
    )
    lines.append("       fixture. See docs/archive/JUNCTION_OPERATING_POINT_271.md.")
    return "\n".join(lines)
