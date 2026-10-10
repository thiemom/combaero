# Junction model validation

Consolidated dataset + scoring runner for steady-flow T-junction loss
coefficients. Used to compare candidate junction models against measured
data from the canonical references, with the published analytical
correlation as a "best-1D-physics-can-do" ceiling.

## Layout

```
validation/junction/
  models/                      # analytical correlations (= ceiling references)
    bassett2001.py             # all 12 K from Table 2 + Eq 33/34/35/37 forms
    hager1984.py               # xi_t (Eq 8), xi_l (Eq 19)
    wang2014.py                # Tables 1 + 2 exact-match lookup
    perez_garcia2010.py        # Table 1 (s, m, n-1) regression params
  data/<paper>/                # measured + paper-published calc curves
    README.md                  # paper-specific conventions, fit-tier ratings
    <paper>_figXX_Kid_<params>_<measured|calc>.csv
    metadata.yaml              # K_id, theta, psi, q axis, uncertainty per file
  runner.py                    # iterate (model, dataset) -> records
  scorecard.py                 # records -> metrics + scorecards
  random_robustness.py         # convergence on RANDOM physical BCs, no dataset
```

## Excluded sources

The directory `docs/junction/` includes several PDFs not used in this
dataset:
- **Stigler 2010**: theoretical paper, no measurements; paper itself states
  measurements "are going to be done soon" (i.e., never within the paper).
- **Torregrosa 2017**: unsteady wave propagation (transmission /
  reflection coefficients), not steady K. Could be a separate
  acoustic-tier dataset later if needed.

## What each source can measure

None of the sources is the closure's own paper (Mynard 2015), so nothing here
scores fidelity -- and no fidelity check against Mynard's own results exists
yet: `tests/test_mynard_junction.cpp` holds the C++ port to the Python one, not
the Python one to the paper. What separates them is whether a source set a constant.
`network_runner.SOURCE_ROLES` declares it and the network scorecard prints it
per row:

| source | role | why |
|---|---|---|
| Bassett 2001 | selected | chose `eta_scale = 0` and the dividing-streamline recovery (with Hager); its K11/K12 also weighed in alpha |
| Hager 1984 | selected | the same decision, with Bassett |
| Idelchik 1966 | tuned | `joining_etransfer_alpha = 0.2` is its in-network optimum. Its theta = 30 tables (diagrams 7-1/7-2) are calculating formulas with no experiment behind them (the handbook's own diagram list); 45 and 90 are formulas refined by Kinne's and Vogel's experiments |
| Wang 2014 | x-source | set nothing -- the one CROSS-SOURCE accuracy check |

## Endpoints

At q = 0 or 1 one port carries no flow. Its residual has a kink where the flow
sign flips, so an imposed_q solve there either stalls or is refused as an
unphysical root; until #272 every such record went unscored (Idelchik's
q = 0 / 1 columns, Wang's q = 0 / 1 curves, Bassett's digitised end points --
173 records). The measured value at an endpoint is the one-sided limit of a
curve whose flow has a direction, so `network_runner` solves it
`_ENDPOINT_EPS` = 1e-3 inside the range, scores it against the measured
endpoint value, and marks the record `endpoint_limit`. The limit is stable:
every endpoint cell's bias moves <= 0.03 between 1e-3 and 1e-4.

## Tier 2 (finite Mach, #272)

- **Wang 2014 is the Tier-2 source.** Cross-source: 188 of 200 points, MAE
  0.131 (interior 0.108, the q = 0 / 1 curves 0.170), flat with Mach (0.124
  below M 0.15, 0.128 above 0.45), worst at area ratio 2.44 (0.195) and at
  K_23 with a dead straight inlet (q = 1, bias -0.33). Wang states no
  uncertainty, so there is no band to be within;
  `python/tests/test_junction_tier2.py` pins regression floors and a
  falsification. The 12 unscored points are all at a = 2.44 above M 0.49.
- **Perez-Garcia 2010 is not a validation source for the closure.** Its K_hat
  (Eq 41) depends only on the two branch Mach numbers; a junction with NO loss
  reproduces its Table 1 within U95 in 57 of 72 cells, every one at
  M3* <= 0.3. Its band is worth 22 units of Miller K at M3* = 0.3 and only
  approaches the model's own error (0.1-0.4) at M3* = 0.7 with a branch
  carrying >= 3/4 of the flow. The same test file pins this so nobody wires a
  vacuous K_hat check. Table 1 was re-transcribed (seven values were wrong).

## Joining flow types

Bassett defines six flow types. Types 4 and 6 are both joining tees with one
lateral inlet, one straight inlet and one straight outlet; they differ in WHICH
straight leg is the outlet, so the lateral joins pointing the other way along
the main duct. Type 4 is therefore scored on the same network with the lateral
mirrored to `pi - theta`.

Each coefficient is indexed on the fraction in its own leg, and which leg that
is comes from Table 1, not from the coefficient's name: for type 4 the common
branch is A, so K7 is the lateral coefficient and K8 the straight one.

## Mach-indexed sources

Most files are a coefficient against the flow split. Wang 2014 is a Mach sweep
at a FIXED split, and it is the only measured compressible data in the set: 200
points from Mach 0.09 to 0.60. The network runner reads those files with the
roles swapped -- abscissa is the Mach, split comes from the metadata -- places
the network at that Mach, and extracts K with Wang's own normalisation
(the common port's total minus static, not `1/2 rho u^2`, which differ by about
9% at Mach 0.6).

The scorecard keeps a Mach band as its own axis, so an incompressible closure
is never averaged across one.

## Convergence is measured somewhere else

The scorecard's convergence column is measured on cases whose boundary
conditions are built from Bassett's analytical K at a target split. That is
right for accuracy and circular for robustness: it asks how often the solver
reaches an operating point the paper's own correlation predicts, on the same
points the closure has been measured against.

`random_robustness.py` asks the production question instead. It draws geometry
and boundary conditions uniformly inside physical ranges, with no reference to
any paper, and reports whether the solver returns an admissible answer:

    uv run python -m validation.junction.random_robustness 2000

The width of the sample space is the point -- a narrow space could be
flattered by tuning. `python/tests/test_junction_random_robustness.py` pins the
ranges so narrowing them shows up in a diff, and keeps the convergence floor
loose so it catches regressions rather than inviting anyone to optimise
against it.

It separates draws with **no root** from solver failures. A third of random
draws constrain a function of the split to a value the closure cannot produce;
counting those as failures is a category error.

## How candidate models are scored

For each `ValidationCase` (one row of a measured CSV):
- Evaluate the candidate model at the case's (q, psi, theta, M).
- Evaluate the paper's own analytical correlation (the "ceiling").
- Compute `MAPE_meas` (vs measured truth) and `MAPE_ceil` (vs the ceiling
  correlation). Headline metric: `Delta_vs_ceiling = MAPE_meas - MAPE_ceil`.
  Negative = candidate beats published 1D physics.

The scorecard reports per (model, K, regime, psi-bin, theta-bin):
- `N`, `RMSE_meas`, `MAE_meas`, `bias_meas`, `pct_within_uncertainty`,
  `Delta_vs_ceiling`. `pct_within_uncertainty` uses each paper's STATED
  measurement uncertainty and shows `-` where there is none -- which today is
  every junction source: Bassett and Hager compile other experimenters' data
  without a band, Idelchik's handbook disclaims one, Wang states none.
- For network-mode evaluation: `pct_converged`, `median_wall_time_ms`,
  `median_residual_norm`.

A candidate model can have an `excellent` score on K2/K5/K6/K12 (the
well-fit Ks) and a `poor` score on K9/K10 and still be correct — because
the underlying 1D physics itself struggles in type 5. The fit-tier ratings
in `data/bassett2001/README.md` capture this.
