# NASA CR-3630 (Florschuetz, Metzger and Su, 1982) -- initial crossflow, heat transfer

**Status: UNCONFIRMED (transcribed 2026-10-23, #467).** Source:
`docs/heat_transfer/flohrshuetz/NASA_CR-3630_Florschuetz_Metzger_Su_1982_initial_crossflow.pdf`
(gitignored), Appendix D, Table D.2, "Heat transfer data", printed pp. 150-161
(PDF pp. 162-173). Data: `data/florschuetz_metzger_su1982/table_d2_heat_transfer.csv`.
The flow side of the same report (Table D.1) is a separate, signed-off record
(`florschuetz_metzger_su_1982_cr3630.md`).

## For the reviewer

1. **Spot-check the pages against the CSV.** Every test's printed NU-bar
   matched the mean of its ten in-array NU (the self-check below), so a
   misread tens or units digit in Nu would have shown. ETA has a weaker
   check: 1.00 upstream of the array. One cell was read off the image
   (p172, mc/mj 0.51, x/L 0.45: text layer `0041`, printed 0.41), and one
   header (p165, B(5,8,1)I; the text layer lost it). *What turns on it:* a
   cross-source check of the eta_r reference temperature (#713) and of
   Florschuetz's h with initial crossflow.
2. **Matching D.2 tests to D.1 flows.** D.2 prints no Gc/Gj, so a score
   reads the row flows from the D.1 test of the same geometry and nearest
   mc/mj (within 0.05). Are these the same runs? The mc/mj and Re_j agree to
   two figures; I did not find it stated.

## What the table prints

Per geometry B(xn/d, yn/d, z/d), one page each, four tests (mc/mj 0, about
0.2, 0.5, 1.0): mean Re_j (thousands), NU-bar (mean over 0 < x/L < 1), and per
segment x/L = -0.25 ... 1.05 the Nusselt number NU and ETA.

- **Definitions (Sec. 5, Eqs. 5.1, 5.3):** q = h (T_s - T_aw), and
  ETA = (T_aw - T_j)/(T_c - T_j), with T_c the INITIAL crossflow temperature
  at the entrance. Jets are the primary flow, the crossflow the secondary.
- x/L < 0 are the segments upstream of the array, under the initial
  crossflow alone (ETA = 1.00 there). The zero-crossflow test prints
  `******` for them and for every ETA.
- 12 pages: 11 inline geometries (no B(10,4,1), which D.1 has) and
  B(5,4,3)S. 48 tests, 672 rows.

## Checks that could have failed

- **NU-bar = mean of the ten in-array NU,** every test, to 0.6% (rounding).
  All 48 pass. (The two image reads were caught by the value count and the
  header, not by this.)
- **ETA = 1.00 upstream** in every test with initial crossflow. In the
  array ETA runs from 0.05 to 1.02: one printed overshoot past 1 (p163,
  mc/mj 0.56, row 1), the kind CR-3936 reports where the crossflow controls
  the surface.
- **Review page** (scans beside the CSV, the patched cell highlighted):
  https://claude.ai/artifact/7zbTbGYpnD4LCHwazwit9F

## Previews, not gates (2026-10-23, exploratory, not committed as tests)

These are cross-source: same lab, but a different study from CR-3936, which
the closure was fitted to. They are accuracy, not fidelity.

- **The bypass reference temperature (#713, option B):** the model's
  T_aw - T_j = eta_r(Gc/Gj, z/d) phi_0 (T_0 - T_j) predicts exactly this
  table's ETA as eta_r phi_0. Per row, from the matched D.1 flows, over 28
  inline tests and 279 points with initial crossflow:

  | prediction | bias | RMSE |
  |---|---|---|
  | **eta_r phi_0 (as built)** | **+0.024** | **0.077** |
  | phi_0 alone (eta_r = 1) | +0.177 | 0.219 |
  | none (eta = 0) | -0.439 | 0.507 |

  By row, the bias runs from -0.07 at row 1 to +0.07 at row 10.

- **Florschuetz (1981) h with initial crossflow,** at each row's Re_j and
  D.1 Gc/Gj, against NU:

  | subset | n | bias | RMS |
  |---|---|---|---|
  | no initial crossflow | 100 | +4.7% | 7.7% |
  | with initial crossflow | 279 | -10.0% | 25.5% |
  | Gc/Gj <= 0.8 | 282 | -1.4% | 17.0% |
  | Gc/Gj > 0.8 (flagged out of range already) | 97 | -19.7% | 32.9% |

  This is the under-prediction CR-3936 expects where the crossflow
  dominates (Nu rising toward channel flow, its conclusion 6).

## Review log

- 2026-10-23: transcribed from the text layer per page; one cell and one
  header read off the page images; the NU-bar identity holds for all 48
  tests.
