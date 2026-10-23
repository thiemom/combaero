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
- **Gj is inferred, not metered:** from the measured channel P with the
  plate's constant mean Cd (Eq. 4.1, compressible). The network computes the
  same quantity, so the comparison is like for like.
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
- **Constant Cd, as the authors' Eq. 4.1 and their closed form.** Their
  Table 4.1 Cd(Gc/Gj) is not used: it matters at z/d = 1 with large
  crossflow, where the network is reported with a looser bound, not tuned.

## Review log

- 2026-10-22: transcribed from the text layer; two pages checked visually;
  closure identity added after the check failed on B(5,8,1).
- 2026-10-23: signed off by the user; all 11 pages checked against the scans,
  no errors.
