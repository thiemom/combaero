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


def _pct(v: float) -> str:
    return "-" if math.isnan(v) else f"{v * 100:6.1f}%"


def render(cells: list[Cell]) -> str:
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
    print(render(build(run_all(dataset))))


if __name__ == "__main__":
    main()
