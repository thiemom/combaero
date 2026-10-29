# Blaisdell & Manson 1963, Table 12 -- centre laterals

Blaisdell, F. W. and Manson, P. W. (1963), *Loss of energy at sharp-edged pipe
junctions in water conveyance systems*, USDA Technical Bulletin 1283, Table 12
(printed pp. 131-145, pdf pp. 135-149). Hypotheses, model set and metric were
frozen beforehand in [PREREGISTRATION.md](PREREGISTRATION.md); the scorer is
[`validation/junction/bm1963_score.py`](../../bm1963_score.py).

## `table12_centre.csv`

1880 rows, every centre-lateral run of the 55 (theta, Ad/Ab) junctions. Top
laterals (T) are not transcribed (excluded by the pre-registration).

| column | meaning |
|---|---|
| `printed_page`, `side` | printed page number and left/right half of the spread |
| `theta_deg` | nominal angle with the UPSTREAM main, as tabulated |
| `theta_actual` | built angle from footnote 2 (degrees and minutes, d + m/60) |
| `Ad_over_Ab` | actual area ratio (footnote 3): the printed keys 2 and 7 are 2.12 and 7.11 |
| `run`, `R` | run number and Reynolds number as printed (`R` empty where printed `--------`) |
| `Qb_over_Qd`, `zeta_u`, `zeta_b` | as printed; trailing minus signs converted to signed values; empty = printed dashes |
| `*_source` | how each value was decided (below) |
| `off_fit` | residual against a robust cubic in q for the row's (theta, psi) series, where large; informational |

Run numbering has printed gaps (e.g. 23 -> 25, 48 -> 61); all were checked on
the page, no rows are missing.

### How values were decided

Sources: the OBJ scan (`OBJ.pdf`, 600 dpi, cleanest), its own text layer, the
HathiTrust page TIFFs, and a ScanTailor dewarp of them. The scans are
copyrighted and gitignored (`sources/junction/tb1283/`); the extraction scripts
were throw-away and are not in the repo.

Each cell was cropped from the deskewed, rule-gridded pages and read by Apple
Vision OCR; the text layer was aligned by row. A value is

- `maj3` / `maj2`: three / two of {OBJ Vision, OBJ text layer, HT+ScanTailor
  (one vote, only when the two agree)} read it identically, with the bent-scan
  vote unable to win alone;
- `obj=img`: OBJ Vision agreed with a bent scan, the text layer did not;
- `reviewed`: read by eye from the OBJ page and confirmed by a second reader;
- `reviewed, print damaged (digit inferred)`: the tens digit of zeta_b is
  missing in the print (p. 142, theta 150, psi 7.11, runs 9-22) and in every
  scan; inferred from the Q-neighbours of the same series. 10 cells;
- `reviewed, low confidence`: one cell (p. 133, theta 105, psi 1, run 27,
  zeta_b .690);
- `printed blank`: dashes in the print, 2 cells (zeta_u at p. 135 R run 21 and
  p. 137 R run 8). These runs are not scored.

Audit. OBJ Vision and the bent-scan reads share one OCR engine, so their
errors correlate: two of them agreeing is weak evidence. Every value resting
on Vision reads alone (`obj=img`, and `maj2` without the text layer or against
it), every row off its series curve, and every row with conflicting reads was
proofread against the page (about 330 rows; 17 values corrected, 6 of them
damaged tens digits, plus 3 run numbers). What remains
unread by eye has the text layer agreeing with at least one Vision read, i.e.
two independent engines. Two random samples of auto-accepted rows (48 rows,
144 values) found one error, a dropped leading digit that the off-curve
proofread also caught.

`off_fit` rows are B&M's printed values, kept as printed. One is almost
certainly a sign typo in the original: p. 139 R, theta 15, psi 7.11, run 21,
zeta_b printed 24.290- where the series runs at about +23.

## Fidelity: Table 13

Table 13 (pdf pp. 150-167) tabulates B&M's theory. It is computed at the
NOMINAL angles and at Ad/Ab = (Dd/Db)^2 = 1, 2.25, 4, 9, 16 -- not at Table
12's 2.12 and 7.11. From its text layer, 4336 cells were parsed and compared
with `joining_variants.bm_theory` (Eqs 12/19 as printed):

- 4288 agree to the printed digit (counting the text layer's 5 -> 6 misreads
  and last-digit rounding);
- 30 are whole rows whose q the text layer garbled, 10 are other OCR garbles;
- 8 are a slip in the original: the q = 1.00, theta = 75 row was computed with
  cos 75 = .2598 (true .2588), consistently at all five area ratios.

With the built angles instead the agreement collapses (about half the cells
miss), which fixes the basis as nominal. Our reading of their theory is
therefore faithful.

## Results (pre-registered)

`uv run python validation/junction/bm1963_score.py` (built angles; `--nominal`
changes no verdict). 1878 scored runs, 55 cells: 11 development, 44 same-rig
held-out. Headline = mean of per-cell MAE, zeta_u / zeta_b. No band.

| model | all cells | development | held-out |
|---|---|---|---|
| shipped (Mynard, alpha 0.2, eta 0) | 1.173 / 1.084 | 0.354 / 0.239 | 1.378 / 1.295 |
| faithful (Mynard, alpha 0, eta 1) | 1.310 / 1.027 | 0.382 / 0.267 | 1.541 / 1.217 |
| P2 raw (Bassett CV, raw angles) | 0.239 / 1.904 | 0.235 / 0.605 | 0.239 / 2.229 |
| P2c (Hager 3/4 angle law) | 0.776 / 1.313 | 0.256 / 0.301 | 0.906 / 1.566 |
| BM theory (Eqs 12/19) | 0.257 / 1.905 | 0.272 / 0.644 | 0.253 / 2.221 |

- **H1 not supported.** Held-out, psi >= 4, theta <= 90 (17 cells): P2 raw
  beats shipped on zeta_u (0.262 vs 0.925) but not on zeta_b (1.369 vs
  1.252). The prediction required both.
- **H2 supported.** psi = 1: shipped 0.159 / 0.127 vs P2 raw 0.193 / 0.303
  (all 11 cells; the 4 held-out cells agree).
- **H3 partly supported.** Implied d_eff/theta from zeta_u rises with psi at
  theta 45-90 (e.g. 90 deg: 0.78, 0.95, 0.98, 1.01, 1.07) but continues past 1
  instead of approaching it. At theta 15-30 the fit is ill-conditioned (cos
  is flat) and hits its search bounds; those cells say nothing.
- **H4 supported through zeta_b only.** BM theory mean cell bias +0.013
  (zeta_u, essentially none) and +1.786 (zeta_b).

### Exploratory (no prediction)

Implied Mynard alpha (eta 0, searched on [0, 1]) per cell, reported to
describe the closure and never applied:

- at psi = 1 alpha has no effect (all cells return the search floor);
- for zeta_u it rises with theta and is nearly psi-independent for psi >= 2:
  about 0.05-0.2 at 15-30 deg, 0.4-0.5 at 45-60, at the bound 1.0 at 75-105;
- for zeta_b it is mostly 0-0.5 for theta <= 135, near 1 at 150-165 for
  psi <= 7.11;
- so no single alpha serves both branches at psi > 1: the straight branch
  wants a deflection-dependent alpha, the lateral a small one.
