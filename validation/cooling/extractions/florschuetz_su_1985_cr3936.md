# NASA CR-3936 (Florschuetz and Su, 1985) -- crossflow temperature influence factor eta_r

**Status: IMPLEMENTED 2026-10-23** as `combaero::cooling::impingement_eta_r`
and the reference temperature of `ImpingementPlateElement` with a bypass
(#467 PR 2). CONFIRMED 2026-10-23: signed off by the user after a
page-by-page check against the scans, every patched cell plus random rows,
columns and single values, about 90% of the table. No errors found. Source:
`docs/heat_transfer/flohrshuetz/NASA_CR-3936_Florschuetz_Su_1985_crossflow_temperature.pdf`
(gitignored), Appendix A, "Tabular Results in Terms of Individual Spanwise
Row Parameters", printed pp. 133-148 (PDF pp. 145-160).
Data: `data/florschuetz_su1985/appendix_a_row.csv`.

## Review questions and answers

1. **The transcription.** The user checked it against the page scans, side
   by side with the CSV: every patched cell, then random rows, columns and
   single values, about 90% of the values. No errors. The source is good
   quality, and even the cells the text layer garbled are legible on the
   page. Every block also passes the identity check below.
2. **Two blocks kept as printed despite the identity.** p150 block 1
   (B(10,4,2)I, mc/mj 0.20) closes only with a closure scale of 1.25, on
   every row, and its row 1 is still 5% off. p155 block 3 (B(5,4,3)S, mc/mj 1.00) has row 8 Gc/Gj 0.86
   printed where the identity gives 0.92. Both were checked on the image.
   *What turns on it:* 20 of the 414 fitted points. Confirmed as printed by
   the user.
3. **The measurement band.** eta_r's 95% uncertainty is given only as error
   bars on Figs. 6.9-6.20 (PDF pp. 72-83), and the bars are **not** of
   equal height. They are about +-0.01-0.02 at row 1 and low Gc/Gj, and up to
   +-0.15-0.3 at some downstream rows, where T_m - T_j has shrunk (the
   authors' own note). So a single band cannot be stated. *Open:* digitise
   the bars (per point, overplotted in places) if the closure's error is to
   be set against the measurement band. Open; it does not gate the closure.

## What the table prints

Per geometry B(xn/d, yn/d, z/d) with I = inline or S = staggered, one page
each, three tests per page. Per test: the mean jet Reynolds number Re_j
(thousands) and mc/mj (initial crossflow over total jet flow). Per row at
x/L = 0.05 ... 0.95: Nu_r, eta_r, the row's Re_j and Gc/Gj.

- **Definition (Eq. 2.11, recovery factor 1):**
  `eta_r = (T_ref,r - T_j)/(T_m,n - T_j)`. T_ref,r is the surface temperature
  for zero regional heat flux. T_m,n is the crossflow's mixed-mean total
  temperature half a hole spacing upstream of row n. The tabulated values
  keep the recovery term epsilon_r in the reduction (Sec. 6, Sec. 8.1).
- **Zero initial crossflow is absent:** there, T_m = T_j and eta_r is
  undefined. So each geometry has mc/mj about 0.2, 0.5 and 1.0.
- **12 geometries:** 10 inline (xn/d 5, 10; yn/d 4, 8; z/d 1-3) and 2
  staggered, B(5,4,3)S and B(10,8,3)S.
- **Seven rows on pp. 156-160.** The segment heater opposite row 9 burned
  out before the last runs, so row 8 became the guard heater and data stop
  at row 7 (Sec. 6.1). This covers B(10,8,3)S (p156) and the Re_j-variation
  runs on pp. 157-160 (Re_j about 6k, 10k, 20k at a fixed mc/mj).
  Three of them repeat a base test digit for digit: p158 = p148, p159 = p149
  and p160 = p153 (column `repeat_of_page`). They are dropped from any fit,
  because one test counted twice leaks across cross-validation folds.
- 48 blocks, 435 rows, of which 414 are unique: 45 tests over 12 geometries.

## Checks that could have failed

**Gc/Gj follows from the jets upstream of it:**
`Gc_i/Gj_i = (pi/4)/((yn/d)(z/d)) (S N (mc/mj) Re_mean + sum_{k<i} Re_k)/Re_i`.
S is one closure scale per test. CR-3630 prints it as MC (0.91-1.02), but
CR-3936 does not. Here it is the median over the test's rows (column
`closure`), and a row must sit within 5% of it
(`test_the_transcription_closes_on_itself`).

The check found the missing B(5,8,2) header. The OCR gave "B( 5,B,2)", so
the parser had filed p148's blocks under B(5,8,1), where the identity missed
by exactly a factor of 2 (z/d 2 vs 1). Geometry is now keyed by page.

**Glyph patches**, each read off the page image:

| PDF page | text layer | read |
|---|---|---|
| 147 | `0.':J7` (mc/mj, block 3) | 0.97 |
| 147 | `0.Y8` (eta_r row 1, block 3) | 0.98 |
| 148 | `0.Y9` (eta_r row 1, block 3) | 0.99 |
| 149 | `36.u` (Nu_r row 6, block 3) | 36.0 |
| 150 | `26.~` (Nu_r row 3, block 2) | 26.5 |
| 151 | `1 * 11` (eta_r row 10, block 3) | 1.11 |
| 158 | `26.Z` (Nu_r row 3, block 1) | 26.2 |
| 160 | `30.j` (Nu_r row 7, block 2) | 30.3 |
| 160 | `0.7';)` (eta_r row 7, block 3) | 0.75 |

Systematic quirks: `17 .9` (a space in a number), `12. /` (a 7 printed as a
slash) and `ETAA` for ETAR.

## The closure (fitted, agreed with the user 2026-10-23)

CR-3936 publishes **no** correlation for eta_r. Its Figs. 6.9-6.20 plot
eta_r against Gc/Gj per geometry, and the authors state the limits: eta_r
-> 0 as Gc/Gj -> 0, and eta_r approaches 1 when the crossflow dominates,
at Gc/Gj of about 0.5 (z/d 3) to 2 (z/d 1). It is relatively insensitive
to Re_j. Rows 1 and 2 deviate (flow history). So the closure is **ours,
fitted to this table, and labelled as fitted**:

    eta_r = X^n / (1 + X^n),   X = a (Gc/Gj) (z/d)^r
    a = 2.04, r = 1.07, n = 1.37   (least squares, all 39 inline tests, 363 points)

**Selection** (2026-10-23, throwaway venv): 13 physics-based forms were
compared against gradient boosting, a random forest and a Gaussian process,
with leave-one-geometry-out cross-validation (LOGO) over the 10 inline
geometries. RMSE in eta_r:

| form | params | inline LOGO RMSE |
|---|---|---|
| ML gradient boosting / random forest / GP (ceiling) | - | 0.069 / 0.071 / 0.073 |
| X = a G xn^p yn^q z^r | 5 | 0.082 |
| **X = a G z^r (chosen)** | **3** | **0.083** |
| X = a G | 2 | 0.128 |
| crossflow-to-row-jet mass ratio M = G yn z/(pi/4) | 2-4 | 0.13-0.14 |
| mixing cup M/(1+M) | 0 | 0.19 |

- **Stable:** across the LOGO folds a stays in 1.93-2.12, r in 1.04-1.13
  and n in 1.30-1.42. The worst held-out geometry, B(5,4,3)I, has RMSE
  0.116.
- **Matches the authors' statement:** eta_r = 0.9 is reached at Gc/Gj 0.75
  (z/d 3) and 2.44 (z/d 1).
- **No gain** from Re_j, row-1 or row-2 terms.
- **The mass-ratio picture loses:** the free yn exponent is -0.3, not +1.
- **Selection optimism:** the form was picked after seeing all the data.
  With 3 parameters it is small, but no split was held back untouched.

**Envelope**, flagged outside rather than modelled:
- **Inline only.** Every model trained on inline data, ML included, misses
  the 2 staggered geometries by a bias of about -0.12 (RMSE 0.21 for this
  form). Two geometries cannot support a pattern correction.
- **eta_r > 1** (12% of points, at Gc/Gj of order 1 and above). The authors
  attribute it to the crossflow controlling the surface once the jets no
  longer impinge, "normally not required" in design. The bounded form
  misses these points by RMS about 0.15. A free ceiling form A*hill would
  catch them only by fitting A = 1.56 at Gc/Gj -> infinity, which breaks
  the eta_r -> 1 limit, so it was rejected.
- **Ranges:** xn/d 5-10, yn/d 4-8, z/d 1-3, Re_j 2.3k-22.5k,
  Gc/Gj 0.07-4.4.

Accuracy on another rig is unknown: no second source exists.

## Review log

- 2026-10-23: transcribed from the text layer per page; glyphs patched from
  the page images; identity check with a per-test closure scale; closure
  form and labelling agreed with the user. Bars on Figs. 6.9-6.20 checked:
  not of equal height.
- 2026-10-23: signed off by the user: every patched cell checked, then random
  rows, columns and values, about 90% of the table. No errors.
- 2026-10-23: implemented, with the user's choice of how the reference
  applies: to the bypass's share of the crossflow only (T_aw = T_j + eta_r
  phi_0 (T_0 - T_j)), since the wall's self-heating of the spent air is
  already inside Florschuetz's (1981) supply-referenced h. Constants rounded
  to 4 figures: a 2.042, r 1.073, n 1.368 (in-sample RMSE 0.079, bias
  -0.009).
