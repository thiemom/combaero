# NASA CR-3630 (Florschuetz, Metzger and Su, 1982) -- initial crossflow, flow distribution

**Status: IMPLEMENTED and REVIEWED 2026-10-23** as the CR-3630 flow check of
`ImpingementArray(bypass=...)` (#467, #705). Signed off by the user after
checking every Table D.1 page against the scans: no errors found. Source:
`docs/heat_transfer/flohrshuetz/NASA_CR-3630_Florschuetz_Metzger_Su_1982_initial_crossflow.pdf`
(gitignored), Appendix D, Table D.1, printed pp. 137-148 (PDF pp. 149-160).
Data: `data/florschuetz_metzger_su1982/table_d1_flow.csv`.

## Review questions and answers

1. **The transcription.** All 11 pages (42 tests) were checked by the user
   against the page scans, side by side with the CSV: no errors. Before that
   I had checked B(10,4,1) (PDF p. 156) and B(5,8,1) (p. 153) digit by digit,
   and every page rests on the two identities below.
2. **Channel width 18.29 cm.** Inferred, not printed as such: Table 3.1 gives
   18 holes across the channel at yn/d 4, d 0.254 cm, and every plate's
   count gives the same width. Only Re_j -> m_j depends on it, not the split.
   Stands as documented.
3. **Initial crossflow temperature.** Table D.1 gives only T0 (jet plenum).
   The runs are taken as isothermal (bypass at T0); I did not find the
   flow tests' bypass temperature stated. *What turns on it:* the bypass
   density, so its momentum, in proportion to any temperature difference.
   The heat-transfer tests (Table D.2) are not used here. Stands as
   documented.

## What the table prints

Per geometry B(xn/d, yn/d, z/d)I, inline, 10 spanwise rows over L = 12.7 cm:
per test, mc/mj (initial crossflow over total jet flow), mean Re_j, T0, P0
(kPa), mass closure, and per row at x/L = 0.05 ... 0.95: P/P0 (channel static
over jet plenum), Gj/Gj_mean and Gc/Gj.

- **42 tests, 11 geometries.** B(10,4,1) has two tests (mc/mj 0, 0.49) --
  confirmed on the page, not a missing page. Others have four.
- **Gj is inferred, not metered:** from the measured channel P through
  Eq. 4.1 (isentropic orifice, plenum P0/T0) -- with the CROSSFLOW-DEPENDENT
  Cd of Table 4.1 wherever a row's Gc/Gj exceeds its plate's break point,
  solved row by row downstream (Sec. 4, printed pp. 25-27, PDF pp. 37-39).
  Corrected 2026-10-28; the earlier reading here (constant mean Cd) was
  incomplete. See "Table 4.1: the Cd behind the reduction" below.
- **Table 3.1:** d = b = 0.254 cm at xn/d 5, 0.127 cm at xn/d 10; mean Cd
  0.85 B(5,4), 0.80 B(5,8), 0.76 B(10,4) and B(10,8).
- **Rig (Sec. 3):** initial crossflow from its own plenum through two
  spanwise rows of holes, then 24.1 cm (16-95 hydraulic diameters) of
  channel before the array; L starts xn/2 ahead of row 1.

## Checks that could have failed

Both are identities of the table's own definitions, so a misread digit
breaks them (`test_the_transcription_closes_on_itself`):

1. **The jet flows average to one**, every test, to 0.003 (two decimals).
2. **Each Gc/Gj follows from the jets upstream of it:**
   `Gc_i/Gj_i = (pi/4)/((yn/d)(z/d)) (N (mc/mj)/MC + sum_{k<i} g_k)/g_i`,
   worst 4.4% (rounding) over every row with Gc/Gj >= 0.3.

The `/MC` (mass closure, 0.91-1.02) was found by this check failing, not
assumed: without it B(5,8,1) at mc/mj 0.52 (MC 0.91) row 1 misses by 11%.
The printed 1.18 is confirmed on the page. The authors' Gc counts the metered
initial crossflow against the inferred jets.

Text-layer quirks handled: `1 .96` (space in a number), `******` (Gc/Gj
overflow at a starved row, kept blank), garbled row labels (`GJ/?fJ`,
`GJkiJ`; the first parse dropped B(5,4,3) at mc/mj 0.19 on that, caught by
the case count).

## Modelling decisions

- **The development run is not modelled.** At a fixed bypass flow it only
  sets the feed pressure, which floats in the test; the split depends on
  the bypass flow and the momentum it arrives with, which the lead-in
  carries (`test_the_feed_node_does_not_matter`).
- **Constant Cd, as the authors' closed form (Eq. 4.10).** The network
  still uses each plate's constant mean Cd (Table 3.1). This is NOT the
  authors' reduction of Table D.1 at large crossflow (see the next section);
  at z/d = 1 with large crossflow the network is reported with a looser
  bound. Whether to carry Table 4.1 is #800.

## Table 4.1: the Cd behind the reduction (added 2026-10-28, #799/#800)

- **What the paper says** (Sec. 4, printed pp. 25-27, PDF pp. 37-39):
  - Special Cd tests: a single open jet row under a uniform crossflow, measured on the B(5,4) plate at z/d 1, 2 and 3 and on the B(5,8,1) plate.
  - Cd is constant (Table 3.1) up to a break point in Gc/Gj, and falls above it:
    zeta = 1.06 / (Gc/Gj + 0.806)^0.602 (Eq. 4.2, a least-squares fit to the B(5,4,1) points).

  | Plate | Cd up to the break point | Above it |
  |---|---|---|
  | B(5,4,1) | 0.85 to 0.63 | zeta |
  | B(5,4,2) | 0.85 to 0.83 | 1.08 zeta |
  | B(5,4,3) | 0.85 to 0.90 | 1.11 zeta |
  | B(5,8,1) | 0.80 to 0.54 | -0.169 x + 0.893 to 1.5, then zeta |
  | B(10,4,1), B(10,8,1) | 0.76 to 0.54 | -0.128 x + 0.825 to 1.8, then zeta |

  - The (10,x) plates were NOT tested; their law is inferred from the (5,4,1) and (5,8,1) results.
  - The D.1 reduction iterates Eq. 4.1 with this law, row by row downstream (Gc known from the upstream sum).
  - Table 4.2: with a constant Cd, the mass closures at large Gc/Gj were "much larger than unity", and the variable Cd "improved [them] significantly". The one exception the authors name is B(5,8,1) at mc/mj 0.52; they attribute it to the local crossflow at a jet being below the span average.
- **Checked against the printed table** (`tmp/799/selfcheck_var.py`, 2026-10-28): reducing the printed P/P0 with Table 4.1 reproduces the printed Gj/Gj_mean much better on the high-crossflow tests.
  - Worst row: B(5,4,2) at 1.01 goes 0.311 -> 0.079, B(5,8,1) at 0.99 0.255 -> 0.060, B(10,8,1) at 1.01 0.131 -> 0.018.
  - 12 of 42 tests close better and 3 worse.
  - What remains is mostly P/P0 printed to two decimals, which a near-stagnant row amplifies through sqrt(1 - P/P0).
- **Consequence for scoring.** A constant-Cd network is not the authors' model where Gc/Gj exceeds the break point. Scoring its Gj against D.1 there is not like for like. The measured quantity is P/P0, printed to two decimals.
- **Probe results** (2026-10-28, outer fixed point, not shipped):
  - With Table 4.1 in the network, every D.1 test solves.
  - Tests with reversed jets fall from 3 to 0, removing the #799 alternation.
  - P/P0 max error improves: mean 0.0189 -> 0.0166, worst 0.111 -> 0.090, and 14 of 16 moved tests are better.
  - The worst-row Gj error worsens: mean 0.065 -> 0.085, worst 0.27 -> 0.38. It sits at the near-stagnant rows, where the printed P/P0 cannot resolve it, and on the inferred (10,x) laws.
  - With Rogers & Hersh's (NASA CR-2681) grazing-flow Cd instead, 3 tests with reversed jets fall to 1 and P/P0 improves slightly.
  - Florschuetz 1981's chain is untouched by CR-2681's law (its max Gc/Gj, 0.77, is below where that law bites). It would NOT be untouched by Table 4.1, whose B(5,4,1) break point is 0.63.

## Review log

- 2026-10-22: transcribed from the text layer; two pages checked visually;
  closure identity added after the check failed on B(5,8,1).
- 2026-10-23: signed off by the user; all 11 pages checked against the scans,
  no errors.
