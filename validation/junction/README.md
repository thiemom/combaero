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
```

## Excluded sources

The directory `docs/junction/` includes several PDFs not used in this
dataset:
- **Stigler 2010**: theoretical paper, no measurements; paper itself states
  measurements "are going to be done soon" (i.e., never within the paper).
- **Torregrosa 2017**: unsteady wave propagation (transmission /
  reflection coefficients), not steady K. Could be a separate
  acoustic-tier dataset later if needed.

## How candidate models are scored

For each `ValidationCase` (one row of a measured CSV):
- Evaluate the candidate model at the case's (q, psi, theta, M).
- Evaluate the paper's own analytical correlation (the "ceiling").
- Compute `MAPE_meas` (vs measured truth) and `MAPE_ceil` (vs the ceiling
  correlation). Headline metric: `Delta_vs_ceiling = MAPE_meas - MAPE_ceil`.
  Negative = candidate beats published 1D physics.

The scorecard reports per (model, K, regime, psi-bin, theta-bin):
- `N`, `RMSE_meas`, `MAE_meas`, `bias_meas`, `pct_within_uncertainty`,
  `Delta_vs_ceiling`.
- For network-mode evaluation: `pct_converged`, `median_wall_time_ms`,
  `median_residual_norm`.

A candidate model can have an `excellent` score on K2/K5/K6/K12 (the
well-fit Ks) and a `poor` score on K9/K10 and still be correct — because
the underlying 1D physics itself struggles in type 5. The fit-tier ratings
in `data/bassett2001/README.md` capture this.
