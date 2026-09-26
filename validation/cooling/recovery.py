"""Recover the marks a digitisation could not pick.

A paper prints graphs, not tables. Where a figure overplots several symbol
classes, only the spatially isolated marks can be digitised, so a committed
series is a BIASED SUBSET of the source's runs -- biased towards the tails,
because a mark far from the cluster centre is the one that can be told apart
(see #393).

The marks that were not picked are not lost, though. Two facts recover them:

1. **The run is known to exist.** Where a figure stacks two panels over one
   abscissa, each run is plotted once in each panel. A mark visible in either
   panel therefore proves the run exists in both -- so the other panel's
   silence locates a run rather than denying one.

2. **Its value is bounded.** The mark was not picked because it sits inside
   the local cluster. So its ordinate lies within the envelope of the marks
   that ARE visible near that abscissa.

Together those turn a missing mark into an INTERVAL observation at a known
abscissa. That is weaker than a point, and it is not nothing.

**Why the envelope is drawn from other classes and not from a correlation.**
The tight, tempting bound is the fitted curve plus the source's stated
scatter band -- Han prints 8% on `G` and 6% on `R`. That would bound the
observation with the very model being scored, which is the #332/#333 defect
class in a new place. The envelope here uses only positive co-location of
other digitised marks: no model, no card, no legend.

**The cost is that the envelope is only sometimes tight enough to say
anything**, and that is a property of the figure rather than a choice. On
figure 4.46, where every class hugs one correlation, the median envelope is
10% wide and 25 of 34 recoveries land under 15%. On figure 4.51, where the
classes genuinely separate, the median is 35% and only one does. Recoveries
wider than `MAX_WIDTH` are discarded rather than counted as weak evidence.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass

from validation.cooling.schema import SeriesMetadata, load_points

# A run is "the same run" across panels if the abscissae agree this closely.
# Digitisation of the same mark in two panels differs by well under a per
# cent; 5% is loose enough to survive a poor scan and tight enough that
# adjacent runs (typically 1.5x apart in Re) never merge.
PAIR_TOL = 0.05

# Marks this close in abscissa define the local envelope. Wide enough to
# catch neighbours on a sparse panel, narrow enough that a rising series
# does not contribute its own slope to the width.
ENVELOPE_TOL = 0.08

# Discard a recovery whose envelope exceeds this. At 15% the bound is
# comparable to the measurement scatter it sits in; beyond that it admits
# nearly any prediction and would inflate the count without constraining.
MAX_WIDTH = 0.15

# An envelope needs at least this many marks to mean anything.
MIN_MARKS = 2


@dataclass(frozen=True)
class Bound:
    """One recovered observation: a run known to exist, value bracketed."""

    x: float
    lo: float
    hi: float
    n_marks: int  # how many marks defined the envelope

    @property
    def width(self) -> float:
        return (self.hi - self.lo) / statistics.mean((self.lo, self.hi))

    def contains(self, y: float) -> bool:
        return self.lo <= y <= self.hi


def _partner(series: SeriesMetadata, dataset: list[SeriesMetadata]) -> SeriesMetadata | None:
    """The same physical class plotted in the figure's other panel."""
    if series.panel not in ("upper", "lower"):
        return None
    other = "lower" if series.panel == "upper" else "upper"
    for cand in dataset:
        if (
            cand.source.name == series.source.name
            and cand.figure == series.figure
            and cand.panel == other
            and cand.geometry == series.geometry
            and cand.alpha_deg == series.alpha_deg
            and cand.series == series.series
        ):
            return cand
    return None


def recover(series: SeriesMetadata, dataset: list[SeriesMetadata]) -> list[Bound]:
    """Bounded observations for runs this series is missing.

    Empty for anything without a cross-panel partner -- the recovery needs
    independent evidence that the run exists, and this is the only source
    of it the dataset currently carries.
    """
    partner = _partner(series, dataset)
    if partner is None or series.kind != "measured":
        return []

    mine = [p.x for p in load_points(series)]
    # Every mark sharing this panel, from OTHER classes, defines the envelope.
    neighbours = [
        (p.x, p.y)
        for other in dataset
        if other.kind == "measured"
        and other.source.name == series.source.name
        and other.figure == series.figure
        and other.panel == series.panel
        and other.path != series.path
        for p in load_points(other)
    ]

    out: list[Bound] = []
    for point in load_points(partner):
        x = point.x
        if any(abs(x / xm - 1.0) <= PAIR_TOL for xm in mine):
            continue  # this run was picked here too
        ys = [y for xn, y in neighbours if abs(xn / x - 1.0) <= ENVELOPE_TOL]
        if len(ys) < MIN_MARKS:
            continue
        bound = Bound(x=x, lo=min(ys), hi=max(ys), n_marks=len(ys))
        if bound.width > MAX_WIDTH:
            continue
        out.append(bound)
    return out
