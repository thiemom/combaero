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
  mynard_fidelity.py           # the kernel against Mynard 2015's own figures
  gardel1970_score.py          # closures against Gardel 1970's raw points (cross-source)
  digitise_mynard2015.py       # colour extraction of those figures
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

None of the measured sources is the closure's own paper (Mynard 2015), so none
of them scores fidelity; that is done against Mynard's figures (next section).
What separates them is whether a source set a constant.
`network_runner.SOURCE_ROLES` declares it and the network scorecard prints it
per row:

| source | role | why |
|---|---|---|
| Bassett 2001 | selected | chose `eta_scale = 0` and the dividing-streamline recovery (with Hager); its K11/K12 also weighed in alpha |
| Hager 1984 | selected | the same decision, with Bassett |
| Idelchik 1966 | tuned | `joining_etransfer_alpha = 0.2` is its in-network optimum. Its theta = 30 tables (diagrams 7-1/7-2) are calculating formulas with no experiment behind them (the handbook's own diagram list); 45 and 90 are formulas refined by Kinne's and Vogel's experiments |
| Wang 2014 | x-source | set nothing -- the one CROSS-SOURCE accuracy check in the network runner |
| Gardel & Rechsteiner 1970 | x-source | set nothing; water, so scored at kernel level by `gardel1970_score.py`, not by the network runner |

## Fidelity: Mynard 2015 against his own figures

`mynard_fidelity.py` holds the Python kernel (`_mynard2010`) to the paper.
Production runs the C++ port, which is held to the Python two ways:
- **Golden headers:** `tests/test_mynard_junction.cpp` and
  `tests/test_mpce_junction.cpp` compare against values generated from the
  Python.
- **Live:** `python/tests/test_junction_cpp_python_equivalence.py` compares
  through the binding at every Mynard sample (632 cases, 1.1e-14), and fails
  when either header no longer regenerates to what is committed.

The paper has
no tables: Figs 4 and 6-11 plot his model ("Unified0D", red) and his CFD
("Ref3D", black markers) over 31 panels. `digitise_mynard2015.py` extracts both
by colour from the 150-dpi scans into `data/mynard2015/` (see its README).

- **Fidelity.** The kernel in his configuration (eta on, no joining alpha, no
  dividing-streamline recovery) is evaluated at his CFD samples and splined
  between them, as his curves are, and compared with every visible red point:
  29 panels within 2.7 px (<= 0.032 in K), Fig 4 and Fig 6d within 4.2 px.
  Fig 4's dashed eta_j = 0 curve is within 2.2 px. Each production departure
  moves the panels it acts on by 15-80 px, so the check sees it.
- **Ref3D** is reported per flow type and branch, faithful against production
  (`python -m validation.junction.mynard_fidelity`). It is laminar CFD at
  Re 350-2400 and our production closure departs from his model on purpose;
  nothing is tuned on it. All 158 samples: MAE 0.212 faithful, 0.219
  production. Fig 4, where his eta was fitted: 0.022 / 0.126. The six psi != 1
  joining samples, the only ones the joining alpha acts on: 0.198 / 0.202.
- **Conventions resolved against his curves**, stated in `_junction`:
  - the angle is measured between the side branch and the downstream main
    branch, for both flow directions;
  - Fig 9's area sweep holds both inlet velocities equal;
  - **Fig 9's rows are swapped against its caption.** His own Unified0D, and
    the dead-branch limit K -> -1, both put the side inlet on the bottom row.

## Record: what the joining alpha stands in for

`alpha_record.py` scores the frozen closure set (`joining_variants.MODELS`)
at kernel level on every joining source, at the q each network solve reached
(kernel and network agree to 1.4% relative). Nothing is tuned; alpha stays
0.2. MAE, joining only:

| source (role) | psi > 1: shipped | psi > 1: P2 raw | psi = 1: shipped | psi = 1: P2 raw |
|---|---|---|---|---|
| Wang 2014 interior (x-source) | 0.123 | 0.086 | 0.058 | 0.138 |
| Idelchik 1966 (tuned) | 0.439 | 0.159 | 0.138 | 0.178 |
| Bassett 2001 (selected) | 0.288 | 0.179 | 0.117 | 0.187 |
| Mynard Ref3D, Fig 9 area sweep (laminar CFD) | 0.202 | 0.234 | | |

P2 raw is Bassett's joining control volume with raw angles and no fitted
constant. Mynard's closure is that control volume with three changes, shown
algebraically: Hager's 3/4 angle law, one flow-weighted merge direction, and
no pressure term on the lateral opening.

- **At psi > 1 the parameter-free control volume beats the tuned alpha on
  every measured source**, including Idelchik, which alpha was fitted to. At
  psi = 1 the closure wins. The area effect is physics the closure lacks, and
  alpha patches it.
- **No single alpha fits.** Within one Bassett cell the alpha a point implies
  rises from about 0.05 to above 0.5 with the flow split; the pooled optimum
  (0.2-0.25) averages contradictory needs.
- **Blaisdell & Manson 1963** (`bm1963_score.py`, pre-registered) says the
  same per branch: for psi > 1 the zeta_u alpha grows with deflection
  (about 0.1 at 15-30 deg to 1 at 75-105) while zeta_b wants a small one, and
  at psi = 1 alpha is unidentifiable. See
  [its README](data/blaisdell_manson1963/README.md).
- **Reading:** the jet's angle law depends on jet strength. A weak jet is bent
  by the main flow (Hager's 3/4); a fast jet from a narrow lateral holds its
  line. Raw angles win at psi > 1, the 3/4 law at psi = 1, and jets opposing
  the main flow (Bassett's type 4) need the 3/4 law at any psi.
- **Wang is spent for this question.** It was unsealed once to score these
  variants. A closure form chosen from these numbers, e.g. a jet-strength
  angle law, needs new joining data at psi > 1 (the Vogel/Kinne originals,
  Ito & Imai, Oka & Ito 2005, Serre 1994) to be validated.

## Gardel & Rechsteiner 1970 (cross-source)

`gardel1970_score.py` scores three closures on the raw points of Gardel's
Figs 12-15 (`data/gardel1970/`): 1426 points, four shapes, all four flow modes
including the two the other sources lack, where the branch is the common leg.
Every model is evaluated on Gardel's own flow pattern at constant density
(water): the ideal tee through the production C++ kernel, Bassett through his
published corrected coefficient for each flow type, the shipped tee through
Mynard with the production constants. Error = model loss - measured loss, on
Gardel's head (total flow over the MAIN pipe area). No band is stated.

| | n | ideal | Bassett | shipped |
|---|---|---|---|---|
| mode 1, dividing | 282 | 0.358 (+0.358) | 0.094 (+0.090) | 0.102 (+0.084) |
| mode 2, joining | 303 | 0.398 (+0.366) | 0.155 (-0.049) | 0.164 (-0.029) |
| mode 3, branch splits both ways | 523 | 0.600 (+0.600) | 0.256 (+0.178) | 0.260 (+0.170) |
| mode 4, both mains into branch | 318 | 0.781 (+0.748) | 0.308 (+0.190) | 0.307 (+0.182) |
| all | 1426 | 0.550 (+0.535) | 0.214 (+0.115) | 0.219 (+0.113) |

MAE (bias), 2026-10-12.

- **The ideal tee over-predicts everywhere**, as a parameter-free vector
  Borda/momentum loss must: +0.36 on the classic tees, about twice that
  where the branch is the common leg. This is the baseline's measured cost,
  not a defect.
- **Bassett and the shipped tee are level** (0.21 vs 0.22) and both lean
  high. Their worst cells are the branch-common modes 3 and 4 at 90 deg
  (bias +0.09 to +0.46), and Fig 13's rounded edge, rho = 0.2, which no model
  represents (bias +0.26 and +0.28 pooled).
- **Joining at psi = 1.45** (Fig 14, mode 2) is the psi > 1 data the alpha
  question wanted. On the straight inlet (`h_ya`) Bassett and the shipped tee
  both under-predict, by -0.27, and the ideal tee over-predicts, by +0.16.
  This is one shape: it does not settle alpha by itself.
- **At delta = 90 a mode-2 merge is both Bassett type 4 and type 6**, and his
  two corrected forms disagree there; type 6 is used, as in his own 90 deg
  K11/K12 curves. The repo's general control volume (`joining_variants.
  bassett_cv`) is not used: at 90 deg it drops the opening wall term that his
  corrected K11/K12 keep, which moves K by 0.15-0.25.

`python/tests/test_junction_gardel1970.py` pins the plumbing (mode 1 equals
his K5/K6, the kernel equals the vector Borda loss, the h-to-loss signs at the
table ends), the pooled numbers, and a falsification of the metric.

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
