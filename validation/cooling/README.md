# Cooling correlation validation

A scored harness for the cooling correlations, on the pattern
`validation/junction/` uses: digitised source points as data, a model
adapter, and a scorecard. Tracked by #333, under #339.

```bash
uv run python -m validation.cooling.scorecard
```

## Why it scores through the model, not around it

The runner never computes `G` from a correlation set's coefficients. It
calls `cb.evaluate_rib` and bisects the Reynolds number until the chain's
own `e+` reaches the digitised abscissa, so friction factor, `e+` and `G`
are all exercised. Recomputing the formula here would score the model
against a second copy of itself -- the failure this harness exists to
catch, and the one that let a 4-5x error sit inside `1.0 < x < 10.0`.

## Refusing is a result

Two kinds of "outside" are kept apart, because collapsing them produces a
number that looks like model error and is not:

- **Outside the numeric validity** is still an answer. The set extrapolates,
  the scorecard reports the point and counts it under `extrap`.
- **Outside the configuration** is not an answer. `han_1988_orthogonal` has
  `valid_alpha = [90, 90]` and no alpha term in `G`; asked about 45 degree
  ribs it returns a number that scored 26.7% bias -- which reads as model
  error when it is really the wrong correlation entirely. The runner
  refuses, and the scorecard prints the reason.

A series nothing can score shows `-`, never `0`.

## The decade trap

Han plots `e+` two ways: figures 4.46-4.48 use a plain `e+` axis, while
4.51 and 4.54 use `e+ x 10^-2`. A digitiser calibrated on the tick labels
alone drops the multiplier and shifts the data two decades -- which is what
happened to the figure 4.54 series, recovered as `x_scale: 1.0e4` in
metadata and pinned by a test that fails at every neighbouring decade.

**Check the axis multiplier before digitising**, and record it as
`x_scale` rather than folding it into the CSV. The committed coordinates
stay a faithful record of the measurement; the interpretation is metadata,
where a reviewer can disagree with it.

## Cross-verifying a digitisation

Digitised points are only as good as the axis calibration behind them, and
a plausible-looking CSV can be two decades out. Every series therefore
carries a **figure card** in its metadata: what the printed axes and
equations say, read off the page WITHOUT reference to where the digitiser
put anything.

```bash
uv run python -m validation.cooling.verify
```

The card and the points are two independent channels, and the check is
their disagreement. Which one is wrong is not decided by the tool -- it is
raised for a human.

**What the card catches on its own**, with no second digitisation:

| bug | caught by |
|---|---|
| dropped axis multiplier (the real figure 4.54 defect) | `x-span` |
| mis-calibrated axis origin or span | `x-span` / `y-span` |
| two curves on one figure swapped | `printed-exponent` |
| a drawn line that does not reproduce its own printed equation | `printed-curve` |
| double-picked or missed marks | `count`, `distinct` |
| a series read off the wrong panel | `y-span`, `monotonic` |
| an ordinate calibration carried over from another panel | `y-span` + `x-span` PASSING |

Each of those was verified by injecting the bug and confirming the right
check goes red -- except the last, which was caught on live data.

### The carried-over calibration, a worked case

Figure 4.46 stacks two panels sharing one abscissa: `G` above with ticks
8 to 40, `R/(P/e/10)^0.35` below with ticks 2 to 5. The lower panel was
digitised first. On the upper panel the ordinate calibration was left at
the lower panel's `y0 = 2, y1 = 5` where it needed `y0 = 8, y1 = 40`, and
the four points came back as 2.34 to 3.76 where the printed line runs
10.6 to 25.1.

The card reported `y-span` failing on all four and `printed-curve` at
**469.71%**. The diagnosis came from what PASSED: `x-span` was clean.
A shared abscissa cannot be disturbed by an ordinate mis-calibration, so
a clean x with a broken y points at the panel, not at the reading.

Applying the anchors `2 -> 8` and `5 -> 40` on the log axis -- a recovery
with no free parameters -- brought the RMS from 469.71% to 2.11% and
returned `3.787 (e+)^0.2722` against the printed `3.7 (e+)^0.28`, which
confirmed the cause.

**The recovered values were not used.** Their residual drifts
systematically with `e+` (-0.72% to -3.31%), so the assumed anchors are
not quite the true axis ends. A recovery good enough to identify a bug is
not automatically good enough to be data: a 2-3% systematic tilt would
consume half of Han's stated 6% band before the model is involved. The
points were re-picked against a correctly calibrated ordinate instead.

Because the panels share an abscissa, this trap recurs on every panel
switch. Set the ordinate first, on the panel actually being read.

**What still needs a human read**, and is why the cards carry
`NEEDS HUMAN READ` markers rather than guesses:

- which symbol belongs to which geometry in the legend
- whether a mark was assigned to the right series
- the point count, where marks overlap too densely to count from a scan
- a figure whose axes are unlabelled in the first place

**A bound that cannot be read is left null, not invented.** Figure 4.193c's
ordinate is log with minor ticks continuing below the lowest labelled one,
so its lower limit is `null` and unchecked. Writing a plausible number
there would have made the card agree with the data by construction, which
is the whole failure this guards against. The first run of the verifier
caught exactly that mistake in the card itself.

## Filling a card

Before or independently of digitising, read off the page:

```yaml
verification:
  x_ticks: [1, 2, 4, 6, 8, 10]   # the printed tick LABELS
  x_multiplier: 1.0e2            # the `x 10^-2` on the axis label
  y_ticks: [10, 20, 30, 40]
  y_limits: [null, 120]          # only where ticks understate the axis
  printed_curve:                 # where the figure prints the equation
    power_law: {C: 3.7, n: 0.28}
  printed_exponent: 0.35         # where only the exponent is checkable
  monotonic: increasing          # only where the figure or text demands it
  expected_points: 7             # null if the marks are too dense to count
```

Declare `monotonic` only where a violation would be a real defect. Figure
4.53's five marks sit in a tight overlapping cluster and one adjacent pair
inverts, so it is left null: a check that fails on correct data trains
people to ignore it.

## Adding a source

Create `data/<source>/` with the CSVs and a `metadata.yaml`. Every series
carries the provenance record from #333 -- `after`, `page`, `item`,
`geometry`, `extraction`, `confidence`, and a `cross_check` that could have
failed. A series without one is not loaded.

Source scans stay gitignored under `docs/heat_transfer/`. Only the
digitised coordinates are committed; they are measurements, the scans are
copyrighted.
