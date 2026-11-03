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

## A drawn line is not its printed equation

Three times now, a line Han draws has had a measurably different slope
from the equation printed on it:

| figure | drawn line measures | its printed label |
|---|---|---|
| 4.46 lower | slope `+0.00874` | `R/(P/e/10)^0.35 = 3.2`, a constant |
| 4.46 upper | `3.774 (e+)^0.2725` | `G = 3.7 (e+)^0.28` |
| 4.51 lower | `3.427 (e+)^0.3220` | `G_bar = 4.5 (e+)^0.28` |

Each was established against controls that rule out the alternatives. The
4.46 R line tilts at 23x the distortion of the panel it sits in, measured
on that panel's own frame, and in the opposite direction. The 4.51 line
drifts monotonically while the DATA in the same panel, under the same
calibration, scatters without trend -- which no calibration fault and no
picking error can produce.

In all three the data tracks the printed equation better than the drawn
line does. So:

- **Score against printed equations and against data. Never against a
  digitised drawn line.** Series with `kind: correlation` carry
  `scores: null` for this reason.
- A drawn line is still worth digitising: it verifies the axis
  calibration, and `G_bar` on figure 4.46 did that to 0.25%. It is a
  calibration standard, not a source of truth.
- Give drawn lines a wider `curve_tolerance` than data would need, and
  record the deviation rather than treating it as a defect.

## A short span is not a defect, but its slope is not a measurement

Marks hide behind each other. When only three of a class's points are
legible, three is the right answer and re-picking will not find more --
figure 4.51's 60 deg crossed series is exactly that.

What must not happen is a short span passing unnoticed and its fitted
exponent being read as physics. That series spans 0.30 decades and its
slope inverts the ordering every other class shows, which is noise, not a
finding.

So `slope-span` fails a narrow series until the metadata carries
`slope_unreliable: true`. The series still pools, and its points still
score -- only the acknowledgement is compulsory, which puts the limitation
where a later reader will find it rather than in a comment.

## G_bar/G varies by configuration

Figure 4.51 plots nine rib configurations twice: `G` (ribbed side) above,
`G_bar` (ribbed and smooth averaged) below. Digitising both measures
`G_bar/G` directly -- the first time extraction item 10 has been checked
against data rather than against printed coefficients and a figure label.

| class | ratio | | class | ratio |
|---|---|---|---|---|
| 90 deg | **1.2193** | | 60 deg V | 1.1655 |
| 60 deg // | 1.1772 | | 45 deg V | 1.1613 |
| 60 deg x | 1.1528 | | 45 deg ^ | 1.1572 |
| 60 deg ^ | 1.1763 | | | |

The 90 deg value lands on the printed `4.5/3.7 = 1.2162` to **0.3%**, which
is item 10 confirmed. The angled configurations average `1.155`, **4.5%
below** it, against a 2.6% measurement scatter. So the ratio is
configuration-dependent by a small but real margin, and applying the
90 deg value to an angled rib is an error -- which is what `_gbar_reason`
in the runner refuses.

### Two classes are unresolved and must be treated with care

`45 deg //` and `45 deg x` return `G_bar/G` of `0.988` and `0.992`. Their
two panels hold what appear to be the same marks: 5/5 and 4/5 points
agreeing within 2%, abscissae within 0.3%.

Three readings survive the evidence, and **the data cannot separate them**:

1. **The figure duplicates their `G` marks into the `G_bar` panel.** Then
   the `G` panel is right and `G_bar` for these two is unknown.
2. **The figure duplicates the other way.** Then `G_bar` is right and `G`
   is unknown.
3. **These two configurations genuinely have `G_bar ~ G`** -- their smooth
   wall performing almost as well as their ribbed wall. Striking, but not
   impossible for a rib angle that drives strong secondary flow.

A ratio strictly below 1 is unphysical, which leans against (3) -- but
only by 1%, inside the digitisation scatter, so it decides nothing.
Reconstructing the missing panel from the seven clean ratios makes (1)
consistent with Han's text and (2) inconsistent, but **that is an
inference from the text, not evidence about the plot**, and plausibility
does not promote a reading here.

Both series are therefore marked `class_confidence: disputed`. They pool
and their points score; they are excluded from any `G_bar/G` conclusion,
and from the question of whether the two panels rank the same -- which
stays **open**, since it turns on exactly these two classes.

Resolving it needs the primary paper, Han, J.C. et al. (1991), ASME JHT
113, 590, not the reprint.

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

## Redrawing a figure

```bash
uv run python -m validation.cooling.plot --list
uv run python -m validation.cooling.plot --figure 4.46
uv run python -m validation.cooling.plot --figure han2012/4.46
uv run python -m validation.cooling.plot --all
```

Puts the committed data back on axes so it can be held next to the scan.
Numbers tell you a series is self-consistent; they do not tell you the
cloud has the shape the page shows, that a class sits where your eye says
it should, or that two series were swapped in a way the arithmetic
tolerates.

It reads the same metadata the verifier and scorecard use, so a plot
cannot drift from what the checks believe. Scatter is drawn as markers
with one per class, drawn correlation lines as dashes, and the equation
each figure PRINTS is overlaid as a solid line -- so the gap between what
a figure draws and what it claims is visible directly. A disputed class
label is marked in the legend.

Panel frames are off by default (`--frames` to include them): a frame sits
outside the labelled ticks, so drawing it forces the axes open and
squashes the data into a strip, and its job is already done numerically.

Output goes to `validation/cooling/plots/`, and the images **are
committed**. A redrawn figure is derived, but a diff that shows the plot
is worth more than the few hundred kB: a reviewer can see the cloud's
shape without regenerating anything, and a change that moves points shows
up as a changed image.

Every image carries its source citation in a footer -- the book and the
primary paper behind the figure. A reproduction built from our own
measurements is ours to publish as long as the source stays named on it,
which is why the footer is not optional and is rendered before the file
is written.

**Regenerate after changing any data or card**, or the committed image
goes stale against the metadata it was drawn from:

```bash
uv run python -m validation.cooling.plot --all
```

Needs matplotlib: `uv pip install -e ".[examples]"`.

## Everything is namespaced by source

A figure number is not a unique identifier -- two books can each have a
figure 4.46 -- so `data/`, `cards/` and `plots/` are all keyed by source:

```
data/han2012/fig4.46_R_eD0.047_pe10_wh1.csv
cards/han2012/fig4.46_R_scatter.yaml
plots/han2012/fig4.46_redrawn.png
```

`--figure` accepts a bare number while it is unambiguous and the
qualified `source/figure` form always. A bare number matching more than
one source is an **error naming the alternatives**, never a silent pick.

## Adding a source

Create `data/<source>/` with the CSVs and a `metadata.yaml`. Every series
carries the provenance record from #333 -- `after`, `page`, `item`,
`geometry`, `extraction`, `confidence`, and a `cross_check` that could have
failed. A series without one is not loaded.

Source scans stay gitignored under `sources/heat_transfer/`. Only the
digitised coordinates are committed; they are measurements, the scans are
copyrighted.

## A network case: NASA TN D-7156 (#764)

Marek and Juhasz (1973) put a plate into a combustor's exhaust and cooled it
in one of three ways: a slot film on its hot face, a parallel-flow convection
passage under it, or both. They measured 62 runs, each with 10 hot-face
thermocouples. `marek_juhasz_network_runner.py` builds each run as a combaero
network:

- a hot duct carrying vitiated air;
- a `SlotFilmElement` with `film_htc="seban1960"`;
- ten `ThermalWall` segments, one per thermocouple, each carrying the film;
- a ten-segment backside `ChannelElement` chain.

The runner scores the wall temperatures per stage. The extraction is
`extractions/marek_juhasz_1973_d7156.md`, and the pins are in
`python/tests/test_marek_juhasz_d7156_network.py`.

```bash
uv run python -m validation.cooling.marek_juhasz_network_runner [--per-run]
```

Every number is **cross-source accuracy**: the same lab as D-7958's kernel,
but a different study. Nothing is tuned to these data:

- the film runs at the element's default Tu of 0.14, not the authors' fitted
  C_m of 0.03;
- the hot-side h is combaero's duct h, not their rig-tuned 0.0365 Re^0.8;
- the module docstring lists every labelled assumption.

Wall T is predicted minus measured; the band is +-11 K (p. 13).

| stage | runs | bias | RMS | TC 1-3 bias | TC 4-10 bias |
|---|---|---|---|---|---|
| convection only | 15 | -49.2 K | 59.1 K | -28.7 K | -58.0 K |
| film only (adiabatic backside) | 28 | +24.8 K | 34.8 K | +13.6 K | +29.6 K |
| combined | 17 | +2.4 K | 43.9 K | +23.9 K | -6.9 K |
| combined, m >= 1 (29-32, 37, 54, 62) | 7 | +36.6 K | 38.8 K | +28.4 K | +40.1 K |
| combined, m < 1 | 10 | -21.5 K | 47.1 K | +20.8 K | -39.7 K |

How to read the table:

- **Convection only:** the wall reads cold, so the duct h is too low: -26 K
  at 589 K and -61 K at 1033 K. This plate's boundary layer starts at its own
  leading edge, 6 cm upstream, which is why the authors tuned their hot-side
  h to this rig. The convection air's exit temperature (an energy-balance
  check) reads 12 K low.
- **Film only:** the backside is adiabatic, so this stage scores the film
  effectiveness alone.
  - The 589 K runs read +12 K, which is about the band.
  - The 1033 K runs read +36 K, and up to +104 K downstream. That fits a
    real wall that lost heat downwards, which the adiabatic backside leaves
    out. In these runs the stagnant passage's inlet thermocouple reads
    351-455 K, against coolant at 300-316 K, so heat did reach the passage.
  - At Tu 0.12, the authors' fitted C_m (a diagnostic, never adopted), the
    bias is +11.5 K.
- **TC 1-3** sit where conduction into the slot lip acts. The authors saw it
  to 3.8 cm. combaero models neither axial conduction nor the lip, so these
  thermocouples are reported separately.
- **Combined runs:** the near-zero pooled bias is two opposite errors
  cancelling, as the m split shows.

**Method comparison** (combined stage, film-side h):

| film_htc | bias | RMS | m >= 1 bias | m >= 1 RMS |
|---|---|---|---|---|
| seban1960 | +2.4 K | 43.9 K | +36.6 K | 38.8 K |
| none | -7.8 K | 37.3 K | +11.8 K | 14.4 K |
| ballal_lefebvre1973 | +0.7 K | 44.2 K | +32.5 K | 39.6 K |

All three agree at m < 1, where the film-side h hands back to the duct's.
`none` scores better at m >= 1, but that is not evidence for it. It pairs
the duct h, which the convection stage shows is too low, with a film whose
T_aw the film stage shows is too hot.

Diagnostics that moved little:

- a constant wall k of 13.8 or 22.5 W/m K instead of ss316 k(T): 0.3 K;
- the hot gas as plain air instead of vitiated: 3 K;
- the slot's flow joining after the wall, so the gas over it is undiluted:
  0.4 K.
