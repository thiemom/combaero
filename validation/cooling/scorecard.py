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

    sampling      how completely the marks could be picked
    bnd           runs recovered as interval observations
    held          fraction of those the prediction lands inside

A series the set has no binding for reports '-' rather than 0. A zero
there would read as a perfect score for a model that answered nothing,
which is the shape of mistake this harness exists to catch.

**Rows are segregated by sampling completeness and never pooled across
it.** A figure that overplots several symbol classes yields only its
spatially isolated marks to a digitisation, and those are the ones
furthest from the cluster centre -- so a `partial` row's MAE/RMSE/within
are an UPPER BOUND on the model's error, not an estimate of it. The
difference is not academic: pooled, `han_1988_orthogonal` reported 73.0%
within against Han's own printed claim of 95%, which reads as a model
deficiency. Segregated, its completely-sampled series reach 91.1% and
its partially-sampled ones 65.3%. The gap was the pick, not the
correlation. See #393 and `recovery.py`.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass

from validation.cooling.runner import Record
from validation.cooling.schema import Point, load_points


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
    # Sampling completeness -- see `sampling_of`. "partial" means the
    # digitisation could not pick every mark, so MAE/RMSE/within are an
    # UPPER BOUND on the model's error rather than an estimate of it: the
    # marks that could be picked are the ones furthest from the cluster
    # centre (#393).
    sampling: str = "unknown"
    # "fidelity" (the correlation's own paper) or "accuracy" (cross-source).
    # Never pooled: they answer different questions and fail differently.
    basis: str = "unknown"
    n_bounded: int = 0  # runs recovered as interval observations
    held: float = float("nan")  # fraction of those the prediction lands inside

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


# Which paper each correlation set IS. Fidelity can only be judged on data
# from that paper; everything else is accuracy, and must be labelled so --
# see docs/VALIDATION_POLICY.md.
#
# This cannot be read off the C++ set: `RibCorrelationSet.provenance` is an
# Extracted/Fitted/User enum, not a citation. It is matched against a
# series' `after` field (the original author of a reprinted figure) and
# falling back to its source name, so a textbook reprint of the set's own
# paper still counts as the author's own data.
SET_ORIGIN: dict[str, str] = {
    "han_1988_orthogonal": "ASME J. Heat Transfer 110, 321",
    "han_park_1988_angled": "IJHMT 31(1), 183",
    "rallabandi_2009_high_re": "rallabandi",
    "florschuetz_1981_inline": "florschuetz1981",
    "mcgreehan_schotsch_1988_cd": "mcgreehan_schotsch1988",
    "mcgreehan_schotsch_1988_crossflow_cd": "mcgreehan_schotsch1988",
}


def basis_of(series, set_name: str | None) -> str:
    """"fidelity" when the data is the correlation's own paper, else "accuracy".

    A set with no recorded origin returns "unknown" rather than defaulting
    to fidelity, so adding a set without declaring its paper cannot quietly
    claim the stronger of the two.
    """
    origin = SET_ORIGIN.get(set_name or "")
    if origin is None:
        return "unknown"
    haystack = f"{series.after or ''} {series.source.name}".lower()
    return "fidelity" if origin.lower() in haystack else "accuracy"


def sampling_of(series, dataset) -> str:
    """How completely this series' marks could be picked.

    "complete"  a table, a drawn curve, or a panel it does not share --
                nothing could have occluded a mark.
    "partial"   MEASURED shortfall: the figure's other panel proves runs
                exist that this panel has no mark for.
    "unknown"   shares a panel with other classes, but nothing independent
                says how many runs it should have.

    Derived, never declared, so it cannot go stale against the data.
    """
    from validation.cooling.recovery import PAIR_TOL, _partner

    if series.extraction == "tabulated" or series.kind in ("correlation", "frame"):
        return "complete"
    partner = _partner(series, dataset)
    if partner is not None:
        mine = [p.x for p in load_points(series)]
        missing = sum(
            1
            for q in load_points(partner)
            if not any(abs(q.x / x - 1.0) <= PAIR_TOL for x in mine)
        )
        return "partial" if missing else "complete"
    shares = sum(
        1
        for other in dataset
        if other.kind == "measured"
        and other.source.name == series.source.name
        and other.figure == series.figure
        and other.panel == series.panel
        and other.path != series.path
    )
    return "unknown" if shares else "complete"


def score_recovered(series, dataset) -> tuple[int, float]:
    """Check the set against runs the digitisation could not pick.

    Returns (count, fraction the prediction lands inside the bound).

    Deliberately NOT folded into MAE/RMSE. A point observation gives a
    signed error; an interval gives containment. Averaging the two into
    one number would be a category error, and would let a loose bound
    flatter a set that a tight point disagrees with.

    Read `held` with the envelope in mind: it is the spread of OTHER
    classes near that abscissa, so it measures whether the prediction
    lands in the local cloud -- not whether it matches this class.
    """
    from validation.cooling.recovery import recover
    from validation.cooling.runner import run_series

    bounds = recover(series, dataset)
    if not bounds:
        return 0, float("nan")
    probes = [Point(x=b.x, y=(b.lo + b.hi) / 2.0) for b in bounds]
    records = run_series(series, probes)
    checked = [
        (b, r.predicted)
        for b, r in zip(bounds, records, strict=True)
        if r.predicted is not None
    ]
    if not checked:
        return 0, float("nan")
    inside = sum(1 for b, pred in checked if b.contains(pred))
    return len(checked), inside / len(checked)


def build(records: list[Record], dataset=None) -> list[Cell]:
    """One row per series.

    Pass `dataset` to fill in sampling completeness and the recovered
    interval observations. Without it those columns stay 'unknown' and
    empty -- the metrics are unchanged either way.
    """
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
        if dataset is not None:
            cell.basis = basis_of(series, cell.scored_by)
            cell.sampling = sampling_of(series, dataset)
            cell.n_bounded, cell.held = score_recovered(series, dataset)
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
    buckets: dict[tuple[str, str, str, str], list[Cell]] = defaultdict(list)
    for c in cells:
        if c.scored_by is None:
            continue
        group = c.label.split("  [")[1].rstrip("]") if "  [" in c.label else ""
        # Segregated by sampling completeness as well as by set and group.
        # A tabulated series and an overplotted one do not measure the same
        # thing: the first gives the model's error, the second an upper
        # bound on it, because only the marks furthest from the cluster
        # centre could be picked. Averaging them produces a number that is
        # neither (#393).
        buckets[(c.scored_by, group, c.sampling, c.basis)].append(c)

    out: list[Cell] = []
    for (set_name, group, sampling, basis), cs in sorted(buckets.items()):
        n = sum(c.n for c in cs)
        n_scored = sum(c.n_scored for c in cs)
        tag = f"  [{basis[:3]}]" + ("" if sampling == "complete" else f" <{sampling[:4]}>")
        agg = Cell(
            label=set_name + (f"  [{group}]" if group else "") + tag,
            kind=f"{len(cs)} series",
            scored_by=set_name,
            n=n,
            n_scored=n_scored,
            n_extrapolated=sum(c.n_extrapolated for c in cs),
            sampling=sampling,
            basis=basis,
            n_bounded=sum(c.n_bounded for c in cs),
        )
        held = [(c.held, c.n_bounded) for c in cs if c.n_bounded]
        if held:
            agg.held = sum(h * k for h, k in held) / sum(k for _, k in held)
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
        f"{'MAE':>7} {'RMSE':>7} {'bias':>7} {'within':>7} {'extrap':>6} "
        f"{'basis':<9} {'sampling':<9} {'bnd':>3} {'held':>7}"
    )
    lines = [head, "-" * len(head)]
    for c in cells:
        lines.append(
            f"{c.label:<44} {c.kind:<12} {c.n:>3} {c.n_scored:>6} "
            f"{_pct(c.mae)} {_pct(c.rmse)} {_pct(c.bias)} {_pct(c.within)} "
            f"{c.n_extrapolated:>6} {c.basis:<9} {c.sampling:<9} "
            f"{c.n_bounded if c.n_bounded else '':>3} {_pct(c.held)}"
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
                f"{c.n_extrapolated:>6} {c.basis:<9} {c.sampling:<9} "
                f"{c.n_bounded if c.n_bounded else '':>3} {_pct(c.held)}"
            )
        if any(c.sampling == "partial" for c in summary):
            lines.append("")
            lines.append(
                "  <partial>: only the marks a digitisation could separate are in "
                "these rows, and those are the ones"
            )
            lines.append(
                "  furthest from the cluster centre, so MAE/RMSE/within are an "
                "UPPER BOUND on the error, not an estimate."
            )
            lines.append(
                "  <unknown>: shares a panel with other classes, with nothing "
                "independent to say how many runs it should have."
            )
        lines.append(
            "  [fidelity]: the correlation's OWN paper -- a miss is our "
            "transcription. [accuracy]: cross-source, so a miss is the"
        )
        lines.append(
            "  model's limitation, not ours. Never pooled. See "
            "docs/VALIDATION_POLICY.md."
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
    print(render(build(records, dataset), pools))


if __name__ == "__main__":
    main()
