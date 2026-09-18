"""Aggregate runner records into a per-series scorecard.

Metrics per series:
    N             digitised points
    N_scored      points the set had a binding for
    MAE           mean absolute relative error
    RMSE          root mean square relative error
    bias          signed mean relative error (model - source)
    within        fraction inside the source's stated uncertainty
    extrap        points outside the set's advisory validity

A series the set has no binding for reports '-' rather than 0. A zero
there would read as a perfect score for a model that answered nothing,
which is the shape of mistake this harness exists to catch.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass

from validation.cooling.runner import Record


@dataclass
class Cell:
    label: str
    kind: str
    scored_by: str | None
    n: int = 0
    n_scored: int = 0
    mae: float = float("nan")
    rmse: float = float("nan")
    bias: float = float("nan")
    within: float = float("nan")
    n_extrapolated: int = 0
    reason: str | None = None  # why nothing was scored

    @property
    def unsupported(self) -> bool:
        return self.n_scored == 0


def build(records: list[Record]) -> list[Cell]:
    grouped: dict[str, list[Record]] = defaultdict(list)
    for r in records:
        grouped[r.series.label].append(r)

    cells: list[Cell] = []
    for label in sorted(grouped):
        rs = grouped[label]
        series = rs[0].series
        errs = [r.rel_error for r in rs if r.rel_error is not None]
        cell = Cell(
            label=label,
            kind=series.kind,
            scored_by=series.scores,
            n=len(rs),
            n_scored=len(errs),
            n_extrapolated=sum(1 for r in rs if r.extrapolated),
            reason=next((r.reason for r in rs if r.reason), None),
        )
        if errs:
            cell.mae = sum(abs(e) for e in errs) / len(errs)
            cell.rmse = math.sqrt(sum(e * e for e in errs) / len(errs))
            cell.bias = sum(errs) / len(errs)
            cell.within = sum(1 for r in rs if r.within_uncertainty) / len(errs)
        cells.append(cell)
    return cells


def pool(records: list[Record], prefix: str) -> tuple:
    """Aggregate every series whose label starts with prefix.

    Where class attribution is provisional -- overlapping symbols, marks
    that appear in only one panel -- the pooled cloud is what the scoring
    should rest on. For a set with no geometry exponents every class must
    land on one curve anyway, so pooling costs nothing and removes a
    dependence on labels that cannot be fully trusted.
    """
    errs = [r.rel_error for r in records
            if r.series.label.startswith(prefix) and r.rel_error is not None]
    if not errs:
        return (0, float("nan"), float("nan"), float("nan"), float("nan"))
    within = [r for r in records
              if r.series.label.startswith(prefix) and r.rel_error is not None
              and r.within_uncertainty]
    mean = sum(errs) / len(errs)
    return (
        len(errs),
        sum(abs(e) for e in errs) / len(errs),
        math.sqrt(sum(e * e for e in errs) / len(errs)),
        mean,
        len(within) / len(errs),
    )


def _pct(v: float) -> str:
    return "-" if math.isnan(v) else f"{v * 100:6.1f}%"


def render(cells: list[Cell], pools: dict | None = None) -> str:
    pools = pools or {}
    head = (
        f"{'series':<44} {'kind':<12} {'N':>3} {'scored':>6} "
        f"{'MAE':>7} {'RMSE':>7} {'bias':>7} {'within':>7} {'extrap':>6}"
    )
    lines = [head, "-" * len(head)]
    for c in cells:
        lines.append(
            f"{c.label:<44} {c.kind:<12} {c.n:>3} {c.n_scored:>6} "
            f"{_pct(c.mae)} {_pct(c.rmse)} {_pct(c.bias)} {_pct(c.within)} "
            f"{c.n_extrapolated:>6}"
        )
    unsupported = [c for c in cells if c.unsupported]
    if pools:
        lines.append("")
        lines.append("Pooled (the scored path where class attribution is provisional):")
        for name, (n, mae, rmse, bias, within) in sorted(pools.items()):
            lines.append(
                f"  {name:<44} {'':<12} {n:>3} {n:>6} "
                f"{_pct(mae)} {_pct(rmse)} {_pct(bias)} {_pct(within)}"
            )
    if unsupported:
        lines.append("")
        lines.append("Not scored:")
        for c in unsupported:
            lines.append(f"  {c.label:<44} {c.reason}")
    return "\n".join(lines)


def main() -> None:
    from validation.cooling.runner import run_all
    from validation.cooling.schema import load_dataset

    dataset = load_dataset()
    records = run_all(dataset)
    pools = {
        "han2012/fig4.46_R_*  (lower panel)": pool(records, "han2012/fig4.46_R_eD"),
        "han2012/fig4.46_G_*  (upper panel)": pool(records, "han2012/fig4.46_G_eD"),
    }
    print(render(build(records), pools))


if __name__ == "__main__":
    main()
