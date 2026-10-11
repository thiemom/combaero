# Blaisdell & Manson 1963 -- pre-registration

Written 2026-10-10, before any Table 12 value was extracted. The commit that
adds this file is the timestamp. Exposure before writing it: the printed page
131 (pdf page 135), rows theta 15/30/45 at Ad/Ab = 1, seen once while checking
legibility. No model was evaluated on it.

## Source and convention

Blaisdell, F. W. and Manson, P. W. (1963), *Loss of energy at sharp-edged pipe
junctions in water conveyance systems*, USDA Technical Bulletin 1283. Table 12
(pdf pp. 135-149) is the run log: 2,217 runs on 64 junctions (55 centre, 9 top
laterals), 2-inch main.

- theta: the angle between the lateral and the UPSTREAM main (p. 8), i.e. the
  jet's deflection into the outlet direction. theta < 90 runs with the main
  flow (Bassett type 6), theta > 90 against it (type 4 mirrored).
- psi = Ad/Ab in {1, 2.12, 4, 7.11, 16}; q = Qb/Qd, the lateral fraction.
- zeta_u, zeta_b: total-head loss from the upstream main / the lateral to the
  downstream main, over the DOWNSTREAM velocity head (Eqs 5-6, 16). This is
  the kernel's joining K at the straight and the lateral inlet.
- Their theory (Eqs 12 and 19, last term of 18 dropped): the momentum CV with
  raw cos(theta) and NO wall-pressure term on the opening.

## Role

B&M is Bassett 2001's reference [3]: his measured joining points come from
it, via Miller. It is therefore not independent of Bassett. It sets nothing in
the closure (alpha was tuned on Idelchik, whose data are Vogel/Kinne).

Cells Bassett showed are **development** (already looked at, in Bassett's
figures):

| psi | theta (B&M) | Bassett figure |
|---|---|---|
| 1 | 30, 45, 60, 90 | 10a/b, 11b, 12a (type 6) |
| 1 | 90, 120, 135, 150 | 12a, 13b (type 4: 180 - Bassett theta) |
| 2.12, 4 | 45 | 10c, 11b |
| 2.12, 4 | 135 | 12b, 13b |

Every other centre-lateral cell is **held out**: theta 15, 75, 105, 165 at
psi 1; theta not in {45, 135} at psi 2.12 and 4; every theta at psi 7.11 and
16. Held-out results are reported as *same-rig held-out*, not cross-source.
Oka & Ito 2005 and Serre 1994 remain the cross-source checks.

Top laterals are excluded: the closure has no offset geometry.

## Models (frozen in validation/junction/joining_variants.py at this commit)

`BM theory` through the CV equals Eqs 12/19 as printed (`bm_theory`) to
3e-14 over q x psi x theta, so "their theory is the CV without the wall term"
is checked, not assumed.


- **shipped**: Mynard kernel, production config (alpha 0.2, eta 0).
- **faithful**: Mynard kernel, his config (alpha 0, eta 1, recovery 0).
- **P2 raw**: the generalised Bassett CV, raw angles, Flamang wall term for
  jets with the flow.
- **P2c**: the same with Hager's 3/4 angle law.
- **BM theory**: B&M Eqs 12/19 (raw CV, no wall term).

All evaluated on the incompressible kernel at the run's (q, psi, theta); water
at under 5 m/s, so no Mach correction.

## Metric

Per run |K_model - zeta_measured|, for zeta_u and zeta_b separately. MAE per
(theta, psi) cell; the headline is the mean of cell MAEs, so cells with many
runs do not dominate. Runs at q = 0 or 1 are scored like any other. No band:
B&M state test statistics, not a measurement uncertainty, so nothing is scored
as "within".

## Predictions

- **H1 (angle law, the question):** on held-out cells with psi >= 4 and
  theta <= 90, P2 raw has lower headline MAE than shipped, for both zeta_u
  and zeta_b. This is the prior from Idelchik, Bassett and Wang at psi > 1.
- **H2:** at psi = 1, shipped is no worse than P2 raw (the prior from all
  three sources).
- **H3 (jet strength):** the effective deflection implied by the data,
  d_eff / theta from zeta_u under the P2 form, rises towards 1 as psi grows
  at fixed theta <= 90. It is reported, not fitted into anything.
- **H4 (sanity of the extraction):** BM theory is biased high against the
  data (their own conclusion, p. 58). A failure here points at the
  extraction before it points at anything else.
- **Exploratory, no prediction:** theta > 90 (jet against the flow) and the
  implied alpha per cell. Implied alpha is reported to describe the closure,
  never to retune it.

## Fidelity check of our reading of their theory

Table 13 tabulates their theoretical coefficients. Our implementation of
Eqs 12/19 must reproduce it to the table's rounding before any scoring runs.
