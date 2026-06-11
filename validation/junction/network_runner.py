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
from dataclasses import dataclass
from typing import Iterator, Protocol

from validation.junction.equivalences import canonical_K
from validation.junction.models._network_builder import (
    ALL_TOPOLOGIES,
    NetworkResult,
    Topology,
)
from validation.junction.schema import Dataset, FileMetadata, load_dataset


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


_LATERAL_K_IDS = {"K1", "K6", "K12"}
_STRAIGHT_K_IDS = {"K2", "K5", "K11"}


def _which_K(K_id: str) -> str | None:
    """Which extracted side does this K_id correspond to."""
    if K_id in _LATERAL_K_IDS:
        return "lateral"
    if K_id in _STRAIGHT_K_IDS:
        return "straight"
    return None


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
        if file.kind != "measured" or file.x_axis != "q":
            continue
        K_id = file.K_id or file.coefficient or ""
        which = _which_K(K_id)
        if which is None:
            continue
        rows = _read_xy_csv(file.path)
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
        errs = [r.K_extracted - r.K_measured for r in conv]
        rmse = math.sqrt(sum(e * e for e in errs) / len(errs)) if errs else math.nan
        bias = sum(errs) / len(errs) if errs else math.nan
        wt_ms = 1000.0 * _median([r.wall_time_s for r in recs])
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
            )
        )
    return cells


def format_network_scorecard(model_name: str, cells: list[NetworkCell]) -> str:
    """Render a human-readable network-mode scorecard with per-topology rows."""
    lines = []
    lines.append(f"=== Network scorecard: {model_name} ===")
    lines.append(
        f"{'canonical_K':<18} {'psi':>5} {'theta':>5} {'topology':<12} "
        f"{'N':>4} {'%conv':>6} {'RMSE':>8} {'bias':>8} {'t_ms':>6}"
    )
    lines.append("-" * 88)
    n_total = 0
    n_conv = 0
    for c in cells:
        n_total += c.N
        n_conv += c.n_converged
        rmse_str = "-" if math.isnan(c.rmse_meas) else f"{c.rmse_meas:.3f}"
        bias_str = "-" if math.isnan(c.bias_meas) else f"{c.bias_meas:+.3f}"
        lines.append(
            f"{c.canonical_K:<18} {c.psi_bin:>5} {c.theta_bin:>5} {c.topology:<12} "
            f"{c.N:>4} {c.pct_converged * 100:>5.0f}% {rmse_str:>8} {bias_str:>8} "
            f"{c.median_wall_time_ms:>5.1f}"
        )
    lines.append("-" * 88)
    overall_conv = n_conv / n_total if n_total else 0.0
    # Per-topology summary
    from collections import defaultdict

    by_topo: dict[Topology, tuple[int, int]] = defaultdict(lambda: (0, 0))
    for c in cells:
        n, conv = by_topo[c.topology]
        by_topo[c.topology] = (n + c.N, conv + c.n_converged)
    lines.append(f"Overall: {n_conv}/{n_total} converged ({100 * overall_conv:.0f}%)")
    for topo, (n, conv) in by_topo.items():
        pct = 100 * conv / n if n else 0.0
        lines.append(f"  {topo}: {conv}/{n} ({pct:.0f}%)")
    return "\n".join(lines)
