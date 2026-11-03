# Seban (1960) and Ballal-Lefebvre SME4 -- film-side h near a slot (#962)

**Status: CONFIRMED.** The user signed off the digitised points and the
curve-B finding on 2026-11-03, in two overlay passes.

Sources (local only, gitignored):

- R. A. Seban, *Heat transfer and effectiveness for a turbulent boundary
  layer with tangential fluid injection*, ASME J. Heat Transfer 82(4),
  303-313 (1960),
  `sources/heat_transfer/film/Seban_1960_JHT_82_303_tangential_slot_film_h_and_eta.pdf`.
- D. R. Ballal and A. H. Lefebvre, *A proposed method for calculating
  film-cooled wall temperatures in gas turbine combustion chambers*,
  Cranfield SME Report 4 (1973),
  `sources/heat_transfer/film/Ballal_Lefebvre_1973_Cranfield_SME4_film_cooled_wall_T.pdf`.

Data:

- `data/seban1960/fig12_h_m_gt_1.csv`: Fig. 12 points at x/s < 40.
- `data/seban1960/near_slot_constant_robustness.csv`: the constant's
  robustness table, for both figures.
- `data/ballal_lefebvre1973/sme4_fig1_near_slot_nu.csv`: SME4 Fig. 1 points.
- `data/ballal_lefebvre1973/sme4_fig1_eq3_line.csv`: SME4's printed Eq. 3
  line, as drawn.

Code: `include/slot_film_htc.h`, `src/slot_film_htc.cpp`. Scored in
`python/tests/test_slot_film_htc_validation_962.py`.

## For the reviewer

Nothing is open. The items below were judged and signed off on 2026-11-03.

1. **The curve-B constant is 0.135 as drawn, not 0.10 as printed.** Three
   independent checks are given below. *What turns on it:* the film-side h
   is 35 % higher than the printed Eq. 13 gives.
2. **Curve C is 0.40 as drawn, against the printed 0.41.** The trace reads
   0.400. *What turns on it:* 2.5 % on h past x/s 37.
3. **The digitised points:** 38 on Seban Fig. 12 and 61 on SME4 Fig. 1, from
   two overlay review passes.

## What the sources give

- **Seban, Fig. 12** ("Heat transfer for mass-velocity ratios greater than
  unity"). The ordinate is `(u_s s / nu_s)^0.3 (h / (rho_s u_s cp)) x 10^2`
  and the abscissa x/s, so Re_s drops out of the plotted group.
  - Curve B is Eq. 13, printed `St_s = 0.10 Re_s^-0.3 (x/s)^-0.3` for x/s
    < 40.
  - Curve C is Eq. 14, printed `0.41 Re_s^-0.3 (x/s)^-0.6` beyond.
  - h is on the adiabatic-wall temperature of the same flow (Eq. 3).
  - Heat-transfer runs (Table 2): m 1.06-8.91, x/s about 4-300 (x = 1.2 in
    station pitch, s 0.250 / 0.125 / 0.063 in), Re_s about 1.4e3-8e3, air
    into air, injection 0-70 F above the free stream.
- **Seban, Fig. 8.** For m < 1, h is the plain flat-plate value from x/s of
  about 10 to 70.
- **SME4 Fig. 1** ("Heat transfer relationship for the near-slot region,
  m < 1.3").
  - The axes are Nu against Re_s x/s, with Eq. 3, `Nu = 0.090 [Re_s
    x/s]^0.7`, drawn through them.
  - The filled circles are Hartnett, Birkebak & Eckert (1961), its ref. 3.
  - The open circles are Pai & Whitelaw (1971), its ref. 7.
- **SME4 Eq. 6 and Fig. 2:** the m > 1.3 wall jet, `Nu = 0.10 Re_s^0.8
  (x/s)^0.44`. Lefebvre & Ballal reprint it as Eqs. 8.43/8.44. Its own
  Cranfield rig points lie 20-67 % below it at low Re (read in the
  2026-11-02 screening; not digitised).

## Calibration and reading

**Seban Fig. 12** (PDF page 9, journal p. 311, top plot). The scan is
200 dpi, rendered at 400 dpi.

- The axes are log-log on the printed grid:
  - x: 4, 6, 8, 10, 20, 30, 40, 60, 80, 100, 200 and 300; 808 px/decade,
    RMS 0.002 dec;
  - y: 1.25, 2, 3, 4, 6, 8 and 10; 785 px/decade, RMS 0.001 dec.
- Read uncertainty on y is about 0.004 dec for a separable symbol and
  about 0.010 dec inside a stack (overlap about half a symbol). x is the
  station column, +-0.005 dec.
- The run behind each symbol is not resolved. The 14 heat-transfer runs use
  near-identical flagged circle, square and triangle glyphs of about 10 px
  at scan resolution. The one exception is the theta at x/s 38.4, which is
  run 128 (s = 1/16 in, m = 1.10).
- 38 points at the 8 stations x/s 4.8 to 38.4 (stations 4, 6, 2, 9, 3, 5, 3
  and 6 points):
  - **5 clean** (isolated symbols);
  - **26 stacked** (touching or overlapping, centre read from the enclosed
    hole);
  - **7 ambiguous** (a possible symbol inferred from a dark overlap).

  The figure is OVERPLOTTED. Fidelity scored on it is an upper bound, and
  only 5 points are cleanly sampled.
- Nothing at x/s > 40 was digitised, so curve C is unscored.
- Curve B was also traced as drawn: 468 column samples, x/s 4.3-37.

**SME4 Fig. 1** (PDF page 17, top plot, rendered at 400 dpi).

- The axes are log-log:
  - x: ticks 1e4 (frame), 3e4-9e4, 1e5, 2e5-9e5 and 1e6 (frame),
    least-squares fit, RMS 0.004 dec, max 0.007;
  - y: ticks 60 (frame), 70, 80, 90, 100, 200-900 and 1000 (frame), RMS
    0.003 dec, max 0.005.
- Symbols were found by a disc and ring template. Read uncertainty is about
  0.005 dec (37 px symbols at about 940 px/decade), plus 0.004 dec of
  calibration: about 0.007 dec per coordinate.
- 61 points:
  - 20 are Hartnett (ref. 3) and 41 Pai & Whitelaw (ref. 7);
  - 59 are clean and 2 touch a neighbour (`overlap`);
  - 21 sit on the printed line (`on_line`, within 22 px).
- One possible point is hidden under an arrow tip (Re_s x/s 2.14e5, Nu 532).
  It is not in the CSV; the robustness table carries it both ways.
- The Eq. 3 line was traced as drawn. Its drawn slope is 0.73, and its
  local constant runs from 0.084 at 1.2e4 to 0.094 at 5e5. With the slope
  fixed at 0.7, the trace's constant is 0.0891, against the printed 0.090.

## The finding: curve B is 0.135, not 0.10

Seban's printed Eq. 13 constant, 0.10, contradicts his own figure. Three
independent checks agree on about 0.135:

1. **The B/C crossing.** The drawn curves cross at x/s 37
   (`(0.40/0.135)^(1/0.3)` = 37.4), where the text switches from Eq. 13 to
   Eq. 14 ("x/s = 40"). The printed constants would cross at x/s 110 and
   leave a 26 % step at x/s 40.
2. **Fig. 9 back-solve.** Fig. 9 back-solves to 0.140-0.144 (extraction
   pass, 2026-11-03).
3. **The data fit.** A fixed -0.3 slope through the 38 Fig. 12 points gives
   0.1336 (pooled per point).

Every data-based variant of the fit lies in 0.1296-0.1339 (robustness CSV):

| variant | n | C (Seban, St form) | x Pr 0.69 |
|---|---|---|---|
| pooled per point | 38 | 0.1336 | 0.0922 |
| excluding all 7 ambiguous | 31 | 0.1333 | 0.0920 |
| stacked blobs counted 0-2 times (9 variants) | 36-40 | 0.1333-0.1339 | 0.0919-0.0924 |
| equal weight per station | 38 | 0.1316 | 0.0908 |
| leave one station out (8) | 29-36 | 0.1296-0.1330 | 0.0894-0.0918 |
| bootstrap stations (2000), per-point pooling | 38 | 0.1334, 90 % CI 0.1273-0.1392 | 0.0879-0.0960 |
| clean points only | 5 | 0.1398 | 0.0965 |
| run 128 alone (m = 1.10) | 1 | 0.1012 | 0.0698 |
| curve B as drawn (trace) | 468 | 0.1350 | 0.0932 |
| printed Eq. 13 | - | 0.1000 | 0.0690 |

On SME4 Fig. 1 the data-based Nu constant is 0.087-0.100:

- 0.0909 pooled per point;
- 0.0930 with equal weight per source;
- 0.0869 for Pai & Whitelaw alone and 0.0996 for Hartnett alone;
- bootstrap 90 % CI 0.0881-0.0938.

This is consistent with Seban's curve and data times Pr 0.69 (0.090-0.093).

What follows:

- **SME4's 0.090** is Seban's curve and data times Pr 0.69. It is not a
  misprint of 0.069.
- **Lefebvre & Ballal Eq. 8.41's 0.069** is the printed 0.10 times 0.69. It
  inherits the misprint.
- **Run 128** is the only identified run near m = 1, and it sits on 0.069,
  25 % below curve B. That is consistent with Seban's finding that h falls
  to flat-plate values as m approaches 1.

## Modelling decisions (user, 2026-11-03)

- **Ship Seban's Stanton form as drawn:** C_B = 0.135 and C_C = 0.40. The
  printed values are recorded as `C_near_printed` / `C_far_printed`.
  - *Why:* three independent checks.
  - *Symptom that would reopen it:* a clean print of Seban's ASME paper
    59-A-177 with a different curve B.
- **C1 blend at the crossing:** a smoothstep in ln(x/s) over +-0.1 decade.
  It stays within 1.5 % of min(B, C).
  - *Why:* it keeps the Jacobian continuous. Seban calls the two-curve
    split "a convenience".
- **Hand-back to the duct's h** past x/s 300 and below m = 1.0, with a C1
  step over 0.1 in m. For the wall-jet option the limits are x/s 50 and
  m 1.3.
  - *Why:* beyond the far end Seban's data stop and Ballal-Lefebvre switch
    to Colburn. Below m = 1 Seban found flat-plate h (Fig. 8). The band
    width is numerical, not physics.
- **m > 1.3: Seban for all m > 1.** The book's wall jet ships only as the
  named option `ballal_lefebvre1973`.
  - *Why:* Seban's Fig. 12 is the high-injection data (m up to 8.91). Eq. 6
    over-predicts its own rig by 20-67 %. On a radiation-heated liner,
    where the film is cooler than the wall, a high film-side h is not
    conservative.

## Scored

See `python/tests/test_slot_film_htc_validation_962.py`.

- **Fidelity, Seban's own data** (upper bound: overplotted, 5 clean
  points). Against the implemented curve:
  - all 38 points: bias -1.0 % (geometric), RMS 0.046 dec (11 %), worst
    33 %;
  - without the ambiguous reads: -1.2 %;
  - the 5 clean points: +3.5 %.

  The RMS is the runs' scatter about one curve (several m and s), not the
  reading (1-2 %). On the same points the printed 0.10 is 25 % low.
- **Accuracy, cross-source** (SME4 Fig. 1; m < 1.3, mostly below Seban's
  box). Curve B at Pr 0.69, all 61 points:
  - bias -2.5 %, RMS 18 %, 32 of 61 within 10 %;
  - Hartnett +6.9 % and Pai & Whitelaw -6.7 %.

  These data lie 32 % above the book's 0.069.

## Review log

- 2026-11-03: Fig. 12 and SME4 Fig. 1 digitised, with the calibration,
  symbols and traced lines checked on overlays (scans and overlays kept
  locally, never committed). The digitising scripts read the scans and are
  not committed; this record and the CSV headers carry the method.
- 2026-11-03: the user reviewed the overlays in two passes and signed off
  the points, the curve-B 0.135 finding and the curve-C trace value 0.40.
- 2026-11-03: implemented (#962): `slot_film_htc.h`, scored as above.
