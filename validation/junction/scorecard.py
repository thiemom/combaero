"""
Aggregate runner records into a per-(model, K, regime) scorecard.

Metrics per cell:
    N                       count of data points
    RMSE_meas               vs measured truth
    MAE_meas                vs measured truth
    bias_meas               signed mean error (model - measured)
    pct_within_uncertainty  fraction within paper's stated uncertainty
    Delta_vs_ceiling        RMSE_meas(model) - RMSE_meas(ceiling); negative = beating literature

Headline (per model):
    core_RMSE_meas    mean over {K2, K5, K6, K12}
    extended_RMSE_meas  mean over all canonical K
    avg_Delta_ceiling  mean over all canonical K (negative = beating literature)
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field

from validation.junction.runner import Record


CORE_K_IDS = {"K_straight_sep", "K_lateral_sep", "K_lateral_join", "K_straight_join"}


@dataclass
class Cell:
    """One row of the scorecard."""

    canonical_K: str
    psi_bin: str  # "1.0" or "all" if mixed
    theta_bin: str  # "45" or "all"
    N: int = 0
    rmse_meas: float = 0.0
    mae_meas: float = 0.0
    bias_meas: float = 0.0
    max_err_meas: float = 0.0
    pct_within_uncertainty: float = 0.0
    delta_vs_ceiling: float = 0.0
    median_wall_time_s: float = 0.0
    n_unsupported: int = 0  # model returned None


def _rmse(errs: list[float]) -> float:
    if not errs:
        return 0.0
    return math.sqrt(sum(e * e for e in errs) / len(errs))


def _median(xs: list[float]) -> float:
    if not xs:
        return 0.0
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def _bin_label(values: set[float | None]) -> str:
    real = {v for v in values if v is not None}
    if not real:
        return "n/a"
    if len(real) == 1:
        v = next(iter(real))
        return f"{v:g}"
    return "all"


def build_cells(records: list[Record], *, uncertainty_K: float = 0.05) -> list[Cell]:
    """Group records by (canonical_K, psi, theta) and compute metrics."""
    groups: dict[tuple[str, float | None, float | None], list[Record]] = defaultdict(list)
    for r in records:
        groups[(r.canonical_K, r.psi, r.theta_deg)].append(r)

    def _sort_key(
        item: tuple[tuple[str, float | None, float | None], list[Record]],
    ) -> tuple[str, float, float]:
        (cK, psi, theta), _ = item
        return (cK, psi if psi is not None else -1.0, theta if theta is not None else -1.0)

    cells: list[Cell] = []
    for (cK, psi, theta), recs in sorted(groups.items(), key=_sort_key):
        N = len(recs)
        model_errs = [r.K_model - r.K_measured for r in recs if r.K_model is not None]
        ceil_errs = [r.K_ceiling - r.K_measured for r in recs if r.K_ceiling is not None]
        unsupported = sum(1 for r in recs if r.K_model is None)
        rmse_m = _rmse(model_errs)
        rmse_c = _rmse(ceil_errs)
        mae = sum(abs(e) for e in model_errs) / len(model_errs) if model_errs else 0.0
        bias = sum(model_errs) / len(model_errs) if model_errs else 0.0
        max_err = max((abs(e) for e in model_errs), default=0.0)
        within = sum(1 for e in model_errs if abs(e) <= uncertainty_K)
        pct = within / len(model_errs) if model_errs else 0.0
        wt = _median([r.wall_time_s for r in recs])
        cells.append(
            Cell(
                canonical_K=cK,
                psi_bin=f"{psi:g}" if psi is not None else "n/a",
                theta_bin=f"{theta:g}" if theta is not None else "n/a",
                N=N,
                rmse_meas=rmse_m,
                mae_meas=mae,
                bias_meas=bias,
                max_err_meas=max_err,
                pct_within_uncertainty=pct,
                delta_vs_ceiling=rmse_m - rmse_c,
                median_wall_time_s=wt,
                n_unsupported=unsupported,
            )
        )
    return cells


@dataclass
class Headline:
    model_name: str
    core_rmse: float
    extended_rmse: float
    avg_delta_ceiling: float
    n_supported: int
    n_unsupported: int


def headline(model_name: str, cells: list[Cell]) -> Headline:
    core = [c.rmse_meas for c in cells if c.canonical_K in CORE_K_IDS and c.N > 0]
    extended = [c.rmse_meas for c in cells if c.N > 0]
    deltas = [c.delta_vs_ceiling for c in cells if c.N > 0]
    n_sup = sum(c.N - c.n_unsupported for c in cells)
    n_unsup = sum(c.n_unsupported for c in cells)
    return Headline(
        model_name=model_name,
        core_rmse=sum(core) / len(core) if core else float("nan"),
        extended_rmse=sum(extended) / len(extended) if extended else float("nan"),
        avg_delta_ceiling=sum(deltas) / len(deltas) if deltas else float("nan"),
        n_supported=n_sup,
        n_unsupported=n_unsup,
    )


def format_scorecard(model_name: str, cells: list[Cell], h: Headline) -> str:
    """Render a human-readable scorecard."""
    lines = []
    lines.append(f"=== Scorecard: {model_name} ===")
    lines.append(
        f"{'canonical_K':<20} {'psi':>6} {'theta':>6} {'N':>4} "
        f"{'RMSE':>7} {'MAE':>7} {'bias':>8} {'max':>7} "
        f"{'%within':>8} {'D_ceil':>8} {'unsup':>6}"
    )
    lines.append("-" * 110)
    for c in sorted(cells, key=lambda c: (c.canonical_K, c.psi_bin, c.theta_bin)):
        lines.append(
            f"{c.canonical_K:<20} {c.psi_bin:>6} {c.theta_bin:>6} {c.N:>4} "
            f"{c.rmse_meas:>7.3f} {c.mae_meas:>7.3f} {c.bias_meas:>+8.3f} "
            f"{c.max_err_meas:>7.3f} {c.pct_within_uncertainty * 100:>7.0f}% "
            f"{c.delta_vs_ceiling:>+8.3f} {c.n_unsupported:>6}"
        )
    lines.append("-" * 110)
    lines.append(
        f"Headline:  core_RMSE={h.core_rmse:.3f}  extended_RMSE={h.extended_rmse:.3f}  "
        f"avg_Delta_ceiling={h.avg_delta_ceiling:+.3f}  "
        f"({h.n_supported} pts supported, {h.n_unsupported} unsupported)"
    )
    return "\n".join(lines)
