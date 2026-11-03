# Marek & Juhasz (1973) -- simultaneous film and convection cooling

**Status: CONFIRMED (transcribed, and signed off by the user 2026-11-03,
#764).** Source: C. J. Marek and A. J. Juhasz, *Simultaneous Film and
Convection Cooling of a Plate Inserted in the Exhaust Stream of a Gas
Turbine Combustor*, NASA TN D-7156 (1973),
`sources/heat_transfer/film/oct23/19730008203.pdf` (local only). Data:
`data/marek_juhasz1973/`. Network case:
`marek_juhasz_network_runner.py`, scored in
`python/tests/test_marek_juhasz_d7156_network.py`.

## For the reviewer

Nothing is open. The items the sign-off settled:

1. **Table II, both unit systems** (printed pp. 15-16, PDF pp. 18-19): 62 runs
   by 20 columns, in SI and in US customary units. The two reads agree in
   every cell (see "How it was read").
2. **The glyphs that needed a decision:**
   - run 59, US inlet convection temperature: printed `8b`, read 88 F (the
     SI cell is 304 K = 87.5 F);
   - run 62, US hot-gas flow: under a smudge, read 2.50 lb/s;
   - leading-dot values (`.018`) are written `0.018`;
   - the print's vertical "same as above" arrows are filled down. The CSV
     headers list every arrow cell.
3. **The wall k print error.** On p. 13 the SI conductivity of type 304 reads
   "0.0081 to 0.013 J/sec-cm-K" (and 0.01 used). The US values beside it,
   8 to 13 Btu/hr-ft-F, are 0.138 to 0.225 J/sec-cm-K. So the SI figure is
   low by a factor of 17. The US values are right for 304 stainless. The
   SI value is not used anywhere.

## What the report gives

- **Rig (pp. 8-11, Fig. 4):**
  - the test section sits in a combustor's exhaust duct, 38.1 x 10.7 cm;
  - the test surface is type 304, 0.159 cm thick, 10.2 cm wide and 12.7 cm
    long;
  - the slot is 0.38 cm high across the full width, its lip 6.0 cm (x_o)
    behind the section's leading edge;
  - a 0.3 cm convection passage runs under the plate in parallel flow;
  - under the passage is a stagnant air gap, then a Hastelloy underside
    cooled by slave air. The slave air cooled the surface by at most 11 K
    at 1033 K;
  - the test section is 3.3 cm high (Fig. 4).
- **Instrumentation (p. 11):**
  - ten Chromel-Alumel thermocouples are brazed into grooves on the hot
    surface's centreline, at a 1.22 cm pitch, x = 1.02 to 11.99 cm behind the
    slot (x/s 2.7-31.6);
  - the film-inlet thermocouple failed early. The film temperature T_s was
    inferred by extrapolating the wall temperatures to the slot exit;
  - T_h is the average of 10 Pt/Pt-13Rh rake thermocouples, pattern factor
    0.33 at 1033 K.
- **Flow system (p. 8, Fig. 2):**
  - ambient air is metered, then heated to 589 K in a direct-fired preheater
    burning ASTM A-1;
  - the test combustor (fuel ASTM A-1 in Fig. 2) fires further for 1033 K.
- **Table I (p. 13):** four conditions at 1 atm:
  - 1.59 kg/s at 589 K;
  - 1.59, 2.27 and 1.13 kg/s at 1033 K;
  - hot-gas velocities 72, 125, 180 and 90 m/s.
- **Accuracy (p. 13):** the plate temperatures are reproducible to +-5 K. The
  combustor drifts, so the test-surface temperatures are accurate to
  **+-11 K** overall. This is the measurement band the network case uses.
- **The authors' model:**
  - Eq. 2 is a wall energy balance;
  - Eq. 3 is the turbulent-mixing film eta. Its C_m 0.03 was fitted to these
    data (p. 14);
  - Eq. 22 is the hot-side h, 0.0365 Re^0.8 Pr^0.33, tuned to this rig;
  - Eq. 5 is Dittus-Boelter for the coolant side.
  - It under-predicted the combined runs by up to 80 K.

## Run classes (Table II)

| class | runs |
|---|---|
| uncooled | 1, 14 |
| film only | 2-13, 20-25, 33, 38-43, 55-57 |
| convection only | 15-19, 34-36, 44-47, 58-60 |
| combined | 26-32, 37, 48-54, 61, 62 |

## How it was read

- **Read 1:** 300 dpi strips of the page images, read in the session.
- **Read 2:** an independent agent, at 400 dpi.
- **Result:** 0 disagreements over 1240 cells (62 runs x 20 columns, the SI
  and US tables each).
- **Table I:** read at 150 and at 300 dpi, no disagreement.
- **The text layer** was used only to locate pages, never as a value.

## Checks that could have failed

- **Cross-unit:** every SI cell agrees with its US counterpart to the print's
  rounding. The one rounding outlier is run 4's hot-gas flow: 1.58 kg/s
  against 3.47 lb/s = 1.574 kg/s. The check is pinned in
  `test_the_si_table_is_the_us_table`.
- **Table I's velocities:**
  - the hot-gas flow area 0.0370 m^2 (38.1 x 10.6 cm less the 3.3 x 10.2 cm
    test section), with vitiated-air density, returns 71.7, 126.0 and
    179.9 m/s;
  - the printed values are 72, 125 and 180 m/s.
- **The run classes partition the table.** A run has film flow exactly when
  it is a film or combined run, and convection flow exactly when it is a
  convection or combined run.
- **Table I's columns:** its m/s and ft/s columns differ by up to 1.3 %
  (234 ft/s is 71.3 m/s, printed 72). They are kept as printed.

## Review log

- Transcription (#764 step 1): both tables were read twice, independently,
  and compared cell by cell (0 disagreements). The glyph cases and the k
  print error were noted.
- 2026-11-03: the user reviewed the transcription against the scans
  (side-by-side review images) and signed it off, the glyph readings and the
  k print error included. The modelling decisions for the network case were
  made in the same review (see the runner's module docstring).
