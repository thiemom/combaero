# Perez-Garcia 2010 — K-hat linking-between-branches coefficient

**Source paper:** J Perez-Garcia, E Sanmiguel-Rojas, A Viedma,
*"New Coefficient to Characterize Energy Losses in Compressible Flow at
T-Junctions"*, Applied Mathematical Modelling 34:4289-4305, 2010.

PDF in `docs/junction/Perez-García-2010.pdf` (gitignored — copyright).

## Conventions

- Compressible flow, 90-deg T-junctions only.
- `K_hat = (p0_03 / p3 - 1) / (p0_0j / pj - 1)`, where index 3 = common,
  j in {1, 2} = inlet/outlet branches (Eq 42).
- `M3*` = extrapolated Mach number in the common branch after frictional
  losses are subtracted (Section 3).
- `q = G_1 / G_3` = inlet 1 (or outlet 1) / common.
- Four flow types per Fig 2:
  - **C1**: combining, K_hat_1 (branch 3 <- branch 1) and K_hat_2 (3 <- 2)
  - **C2**: combining, both branches symmetric -> only K_hat_2 tabulated;
    K_hat_1 = K_hat_2(q -> 1-q)
  - **D1**: dividing, K_hat_1 (3 -> 1) and K_hat_2 (3 -> 2)
  - **D2**: dividing, symmetric -> only K_hat_2; K_hat_1 = K_hat_2(q -> 1-q)

## No digitized data

This paper is **analytical-only** in our dataset. Table 1 (6 power-law
regression correlations) is transcribed verbatim into
`validation/junction/models/perez_garcia2010.py`. Range of applicability
(Section 4.1): `0.15 <= M3* <= 0.7`; `q in {0, 0.25, 0.5, 0.75, 1}`
(numerically tested).

The paper publishes its own measured-vs-numerical comparison points (Figs 3,
4, 5 — 3D regression planes) but they are awkward to digitize from 3D
isometric plots and provide low marginal value over the closed-form Eq 44.

## Use in the scorecard

Perez-Garcia K_hat acts as a **compressible-flow ceiling reference** for
candidate junction models. Range: 0.15 < M3* < 0.7. Outside this range the
correlation is documented as extrapolated; the runner flags such cases.
