# Perez-Garcia 2010 — K-hat linking-between-branches coefficient

**Source paper:** J Perez-Garcia, E Sanmiguel-Rojas, A Viedma,
*"New Coefficient to Characterize Energy Losses in Compressible Flow at
T-Junctions"*, Applied Mathematical Modelling 34:4289-4305, 2010.

PDF in `docs/junction/Perez-García-2010.pdf` (gitignored — copyright).

## Conventions

- Compressible flow, 90-deg T-junctions only.
- `K_hat = (p0_3 / p3 - 1) / (p0_j / pj)` -- no `- 1` in the denominator --
  where index 3 = common, j in {1, 2} = the other branches (Eqs 41/42).
- `M3*` = extrapolated Mach number in the common branch after frictional
  losses are subtracted (Section 3).
- `q = G_2 / G_3` (nomenclature; `q' = 1 - q = G_1 / G_3`).
- Four flow types per Fig 2 (3 is always the common branch):
  - **C1**: combining; 1 straight inlet, 2 lateral inlet -> q is the lateral fraction
  - **C2**: combining, 1 and 2 opposite inlets on the main run, 3 the leg;
    only K_hat_2 tabulated, K_hat_1 = K_hat_2(q -> 1-q)
  - **D1**: dividing; 2 straight outlet, 1 lateral outlet -> q is the STRAIGHT fraction
  - **D2**: dividing, 3 the leg, 1 and 2 opposite outlets; K_hat_1 = K_hat_2(q -> 1-q)

## No digitized data

This paper is **analytical-only** in our dataset. Table 1 (6 power-law
regression correlations) is transcribed verbatim into
`validation/junction/models/perez_garcia2010.py`. Range of applicability
(Section 4.1): `0.15 <= M3* <= 0.7`; `q in {0, 0.25, 0.5, 0.75, 1}`
(numerically tested).

The paper publishes its own measured-vs-numerical comparison points (Figs 3,
4, 5 — 3D regression planes) but they are awkward to digitize from 3D
isometric plots and provide low marginal value over the closed-form Eq 44.

## Not a validation source for the junction closure

K_hat cannot test a loss model: a junction with no loss at all reproduces
Table 1 within its U95 in 57 of 72 cells (every cell at M3* <= 0.3). See
`validation/junction/README.md` ("Tier 2") and
`python/tests/test_junction_tier2.py`. The Table 1 transcription was corrected
on 2026-10-10 (seven values).
