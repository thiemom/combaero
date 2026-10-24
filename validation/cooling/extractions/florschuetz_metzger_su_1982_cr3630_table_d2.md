# NASA CR-3630 (Florschuetz, Metzger and Su, 1982) -- initial crossflow, heat transfer

**Status: REVIEWED 2026-10-24 (#467).** Signed off by the user after
checking the pages against the scans: no errors found. Source:
`docs/heat_transfer/flohrshuetz/NASA_CR-3630_Florschuetz_Metzger_Su_1982_initial_crossflow.pdf`
(gitignored), Appendix D, Table D.2, "Heat transfer data", printed pp. 150-161
(PDF pp. 162-173). Data: `data/florschuetz_metzger_su1982/table_d2_heat_transfer.csv`.
The flow side of the same report (Table D.1) is a separate, signed-off record
(`florschuetz_metzger_su_1982_cr3630.md`).

## Review questions and answers

1. **The transcription.** Checked by the user against the scans, no
   errors (2026-10-24). What the check rested on: every test's printed NU-bar
   matched the mean of its ten in-array NU (the self-check below), so a
   misread tens or units digit in Nu would have shown. ETA has a weaker
   check: 1.00 upstream of the array. One cell was read off the image
   (p172, mc/mj 0.51, x/L 0.45: text layer `0041`, printed 0.41), and one
   header (p165, B(5,8,1)I; the text layer lost it). *What turns on it:* a
   cross-source check of the eta_r reference temperature (#713) and of
   Florschuetz's h with initial crossflow.
2. **Matching D.2 tests to D.1 flows (an assumption, not confirmed by
   the report).** D.2 prints no Gc/Gj, so a score
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

## Scored (2026-10-24, cross-source)

`python/tests/test_impingement_d2_validation.py`. Cross-source: same lab,
different studies from the ones the closures were fitted to (eta_r to
CR-3936, h to Florschuetz 1981). Accuracy, not fidelity.

The rig is solved end to end at each test's own Re_j and mc/mj, with the
bypass 20 K above the jets. ETA is read from the plate's T_aw. B(10,4,2) has
no D.1 test (no P0) and is left out: 40 tests, 400 rows.

- **The bypass reference temperature (#713, option B):** T_aw - T_j =
  eta_r(Gc/Gj, z/d) phi_0 (T_0 - T_j), over 300 rows with initial crossflow.

  | prediction | bias | RMSE |
  |---|---|---|
  | **eta_r phi_0, end to end** | **+0.013** | **0.078** |
  | eta_r phi_0 at D.1's measured flows (279 matched rows) | +0.024 | 0.077 |
  | phi_0 alone (eta_r = 1) | | 0.208 |
  | none (eta = 0) | | 0.517 |

  By row, the bias runs from -0.09 at row 1 to +0.06 at row 10.

  End to end, it first scored 0.163. B(5,4,1) at mc/mj 0.5 and 1.0 starves
  its first rows, and the network ran them slightly reversed, where eta_r
  then dropped to 0. That was fixed in #728: eta_r is taken at the crossflow
  over |m_jet|.

- **Florschuetz (1981) Nu,** at the network's own Re_j and Gc/Gj
  (Nu_model / Nu - 1):

  | subset | n | bias | RMS |
  |---|---|---|---|
  | no initial crossflow | 100 | +5.8% | 8.4% |
  | with it, inside the correlation's range | 195 | -2.7% | 21.5% |
  | with it, flagged outside (`surface_extrapolated`) | 105 | -16.9% | 34.9% |
  | of which starved rows (jet reversed or Gc/Gj > 4.44) | 8 | -93.8% | 93.9% |

  The miss outside the range is the under-prediction CR-3936 expects where
  the crossflow dominates (Nu rising toward channel flow, its conclusion 6).
  The jet correlation has no channel-flow floor; that is the open #467 item.

## Review log

- 2026-10-23: transcribed from the text layer per page; one cell and one
  header read off the page images; the NU-bar identity holds for all 48
  tests.
- 2026-10-24: signed off by the user; no errors found.
- 2026-10-24: scored end to end (the previews replaced); found and fixed the
  reversed-jet eta_r (#728).
