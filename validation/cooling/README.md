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

## Adding a source

Create `data/<source>/` with the CSVs and a `metadata.yaml`. Every series
carries the provenance record from #333 -- `after`, `page`, `item`,
`geometry`, `extraction`, `confidence`, and a `cross_check` that could have
failed. A series without one is not loaded.

Source scans stay gitignored under `docs/heat_transfer/`. Only the
digitised coordinates are committed; they are measurements, the scans are
copyrighted.
