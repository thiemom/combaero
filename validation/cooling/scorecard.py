"""Aggregate every cooling runner into one scorecard.

Three runners exist because three physics families score differently -- ribs
need a Reynolds-number bisection, Nu/Nu1 cancels it, and the orifice chain
works in the source's own pressure coordinates. Each owns its series and
declines the rest by returning nothing, so ``run_dataset`` can give every
series to exactly one of them. Without that, a series declared in metadata but
scored by a specialist runner appeared here as "not scored by any set", and
there was no single place to see whether cooling as a whole was healthy.

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


def run_dataset(dataset) -> list[Record]:
    """Give every series to exactly one runner.

    Ownership is an explicit predicate, not "the runner returned something".
    Every runner returns reason-carrying records for series it cannot score --
    that is the point of them, so a silent zero never masquerades as a good
    score -- which makes truthiness useless as a dispatch test and was the
    first version's bug: the jet-array runner claimed the orifice series and
    reported them unscored.
    """
    from validation.cooling import jet_array_runner, orifice_runner
    from validation.cooling.runner import run_series as rib_run

    specialists = (
        (jet_array_runner.owns, jet_array_runner.run_series),
        (orifice_runner.owns, orifice_runner.run_series),
    )

    out: list[Record] = []
    for series in dataset:
        for owns, run in specialists:
            if owns(series):
                out.extend(run(series))
                break
        else:
            # The rib runner is the fallback and reports why, rather than
            # dropping a series no specialist wanted.
            out.extend(rib_run(series))
    return out


def build(records: list[Record]) -> list[Cell]:
    grouped: dict[str, list[Record]] = defaultdict(list)
    for r in records:
        # A runner may split one series into reported groups -- Rohde's are
        # per velocity-head-ratio band, because pooling them would mix a
        # near-exact comparison with one dominated by a conversion factor.
        group = getattr(r, "group", None)
        grouped[r.series.label + (f"  [{group}]" if group else "")].append(r)

    cells: list[Cell] = []
    for label in sorted(grouped):
        rs = grouped[label]
        series = rs[0].series
        errs = [r.rel_error for r in rs if r.rel_error is not None]
        cell = Cell(
            label=label,
            kind=series.kind,
            scored_by=getattr(rs[0], "scored_by", None) or series.scores,
            n=len(rs),
            n_scored=len(errs),
            n_extrapolated=sum(1 for r in rs if getattr(r, "extrapolated", False)),
            reason=next(
                (r.reason for r in rs if getattr(r, "reason", None)), None
            ),
        )
        if errs:
            cell.mae = sum(abs(e) for e in errs) / len(errs)
            cell.rmse = math.sqrt(sum(e * e for e in errs) / len(errs))
            cell.bias = sum(errs) / len(errs)
            # A series with no stated band gets '-', not 0%. Zero would read
            # as "nothing agreed" when the truth is "the source states no
            # band to agree within".
            cell.within = (
                sum(1 for r in rs if r.within_uncertainty) / len(errs)
                if series.uncertainty is not None
                else float("nan")
            )
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


def rollup(cells: list[Cell]) -> list[Cell]:
    """One row per correlation set -- the whole-of-cooling health view.

    Grouped by (set, reporting group) rather than by set alone. Where a runner
    split a series into groups it did so because pooling them would be
    misleading, and a rollup that ignored that would put the misleading number
    back at the top of the page.

    Unscored series contribute their point count but no error, so a set that
    answered nothing cannot look perfect.
    """
    buckets: dict[tuple[str, str], list[Cell]] = defaultdict(list)
    for c in cells:
        if c.scored_by is None:
            continue
        group = c.label.split("  [")[1].rstrip("]") if "  [" in c.label else ""
        buckets[(c.scored_by, group)].append(c)

    out: list[Cell] = []
    for (set_name, group), cs in sorted(buckets.items()):
        n = sum(c.n for c in cs)
        n_scored = sum(c.n_scored for c in cs)
        agg = Cell(
            label=set_name + (f"  [{group}]" if group else ""),
            kind=f"{len(cs)} series",
            scored_by=set_name,
            n=n,
            n_scored=n_scored,
            n_extrapolated=sum(c.n_extrapolated for c in cs),
        )
        if n_scored:
            # Weight each series by the points it actually scored, so a
            # nine-point series does not count the same as a one-point one.
            w = [(c, c.n_scored) for c in cs if c.n_scored]
            agg.mae = sum(c.mae * k for c, k in w) / n_scored
            agg.bias = sum(c.bias * k for c, k in w) / n_scored
            agg.rmse = math.sqrt(sum(c.rmse**2 * k for c, k in w) / n_scored)
            agg.within = sum(c.within * k for c, k in w) / n_scored
        out.append(agg)
    return out


def _pct(v: float) -> str:
    return f"{'-':>7}" if math.isnan(v) else f"{v * 100:6.1f}%"


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
    summary = rollup(cells)
    if summary:
        lines.append("")
        lines.append("By correlation set:")
        lines.append("-" * len(head))
        for c in summary:
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
    from validation.cooling.schema import load_dataset

    dataset = load_dataset()
    records = run_dataset(dataset)
    pools = {
        "han2012/fig4.46_R_*  (lower panel)": pool(records, "han2012/fig4.46_R_eD"),
        "han2012/fig4.46_G_*  (upper panel)": pool(records, "han2012/fig4.46_G_eD"),
    }
    print(render(build(records), pools))


if __name__ == "__main__":
    main()
