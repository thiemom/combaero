# Hollworth & Dagan (1980) -- jet arrays venting through the target, average heat transfer

**Status: REVIEWED 2026-10-24 (#468).** Signed off by the user against the
overlays. Source: B. R. Hollworth and L. Dagan, "Arrays
of Impinging Jets with Spent Fluid Removal through Vent Holes on the Target
Surface, Part 1: Average Heat Transfer", ASME J. Eng. Power 102 (1980) 994-999
(80-GT-42). Part 2 (local heat transfer) is screened but not extracted;
it is by Hollworth, Lehmann and Rosiczkowski, ASME J. Eng. Power 105 (1983)
393-402, not Hollworth and Dagan (corrected 2026-10-25, from the PDF; the
local file is `393_1_Hollworth_1980_part_II.pdf`). Both PDFs are local only
(gitignored), in `docs/heat_transfer/hollworth/`; see `validation/sources.yaml`.
Data: `data/hollworth_dagan1980/`.

## Review questions and answers

1. **The overlays** on the review page
   (https://claude.ai/artifact/2h5JRhCfCN9vtxaSRvYGL1): the calibrated frame and ticks,
   one cross per standoff station at the cluster centre, and a bar for the
   range the cluster's symbol centres span.
   *What turns on it:* the vented-target coolant-side closure for #468,
   C(X, Z) = 100 Nu/(Re^0.8 Pr^0.33).
   *Answer (user, 2026-10-24):* reviewed, looks good.
2. **X = 15, Z = 1 (Fig. 8)** is the one scattered station. The
   inverted triangle (Re 20000) sits well below the triangle and square, so the
   centre (1.238) carries a wide range (1.12-1.37). Should it be one station,
   or the clean pair plus a separate point?
   *Answer (user, 2026-10-24):* one station, as it is. A split would cost a
   little fidelity resolution, but its points would average back to this
   centre, and a single station is unlikely ever to matter.
3. **Per-symbol values were not read.** The station centre is what the paper's
   own collapse (Eq. 7, Re^0.8) makes the closure value. A per-Re re-digitise
   would be needed to test the exponent; see "What it rests on".

## What the paper gives

- **Geometry:** square arrays of jets (d = 2.54 and 5.08 mm, no size effect),
  impinging on a target with vent holes of the SAME diameter at HALF the jet
  pitch (four vents per jet).
  - Staggered: no vent opposite any jet.
  - Inline: a vent opposite each jet.
  - Isothermal target, air near STP, Pr 0.72.
- **Ranges:** X = x_n/d 5, 10, 15, 20; Z = z/d 1 to X; Re_j 3000-35000.
  Re_j is defined on the total flow over the 16 jets (Eq. 5).
- **Form (Eq. 7):** Nu = f(X, Z, arrangement) Re^0.8 Pr^1/3.
  - Fitted exponents run 0.73-0.87, with best fits of 0.799 staggered and
    0.816 inline.
  - The exponent falls with X.
- **Uncertainty:** h is +/-4% on average and +/-9% at worst (Holman's method).
  That is a measurement band.
- **Printed (Table 2):** Nu_inline/Nu_staggered for X x Z = {5, 10, 15, 20} x
  {1, 2.5, 5, 10}, transcribed to `table2_inline_over_staggered.csv`.
- **Stated in the text:** staggered vents beat edge venting by 35/32/24% at
  X = 10/15/20, averaged over 1 <= Z <= 10. The optimum standoff is
  Z ~ 2.5 for every X.

## The files

| file | from | what |
|---|---|---|
| `fig8_staggered_stations.csv` | Fig. 8 (log axis) | staggered C per (X, Z), 18 stations |
| `fig6_x10_stations.csv` | Fig. 6 (linear axis) | staggered and inline at X = 10, 4 + 4 stations |
| `fig6_x10_fair_curves.csv` | Fig. 6 | the authors' drawn curves, every 0.25 in Z |
| `fig8_ref4_edge_vented_curves.csv` | Fig. 8, "ref. (4)" | Hollworth & Berry (1978) edge-vented curves; ANOTHER study, never data here |
| `table2_inline_over_staggered.csv` | Table 2 | printed ratios |
| `fig8_axis_ticks.csv`, `fig6_axis_ticks.csv` | the figures | tick pixel -> value |

- **Columns:** `C100` = 100 Nu/(Re^0.8 Pr^0.33), as plotted.
  `centre_lo`/`centre_hi` give the range the cluster's symbol centres span
  (cluster extent less one symbol, 27 px). It is a spread, NOT an
  uncertainty.
- **Pixels:** in the 400 dpi crop of PDF page 5 (Fig. 8: 1267 x 1489; Fig. 6:
  1244 x 822).

## Checks that could have failed

Each check uses a channel the values were not read from
(`python/tests/test_hollworth_dagan_data.py`).

- **Table 2:** Fig. 6's inline/staggered ratio reproduces the printed values
  at X = 10 within 0.009, at all four Z. The tolerance is 0.015; the
  printed values are rounded to 0.005.
- **Fig. 6 against Fig. 8 at X = 10** (the same staggered data, linear and
  log axes): within 1.3% at every Z. The tolerance is 2%.
- **The stated gains over edge venting** at X = 10/15/20, from the stations
  against the ref. (4) curves:

  | X | read | stated |
  |---|---|---|
  | 10 | +33.4% | +35% |
  | 15 | +31.4% | +32% |
  | 20 | +24.1% | +24% |

  At X = 15 the ref. (4) curve does not reach Z = 1, so that average uses
  three stations.
- **The optimum standoff** at Z = 2.5 shows in every X with a Z = 1 station.
- **Calibration:** holding out each tick in turn, Fig. 8's log axis lands
  within 0.7% and its Z axis within 0.03; Fig. 6 within 0.014 (C) and 0.008
  (Z).

## What it rests on

- **Overplotted (policy):** the five Re per station overlap, so a fidelity
  score on these stations is an UPPER bound.
  - It rests on 18 + 8 clean station centres, not on per-symbol points.
  - The spread across Re at a station (half-spread 1.7-7.9%, median 4.2%; 10.1% at X = 15, Z = 1) is what the
    exponent's 0.73-0.87 range leaves after the Re^0.8 collapse.
- **One vent geometry:** vent d = jet d, vent pitch = X/2. A closure built on
  it is valid for that pattern only. A second study with another vent
  pattern (Cho & Rhee 2001) is needed before any vent-geometry dependence.

## Review log

- 2026-10-24: automated pass at 400 dpi.
  - Axes come from tick runs; clusters are connected components assigned to
    stations by rule.
  - Fig. 6's drawn curves were fitted (degree 5, outlier-rejected) and masked
    out before reading its clusters.
  - Three independent channels agree (Table 2, the cross-figure check, the
    stated gains).
- 2026-10-24: signed off by the user; X = 15, Z = 1 kept as one station.
