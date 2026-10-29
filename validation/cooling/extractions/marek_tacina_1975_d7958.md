# Marek & Tacina (1975) -- slot film cooling vs free-stream turbulence

**Status: CONFIRMED (transcribed and signed off by the user 2026-10-25,
#714).** Source: C. J. Marek and R. R. Tacina, *Effect
of Free-Stream Turbulence on Film Cooling*, NASA TN D-7958 (1975),
`sources/heat_transfer/film/oct23/19750017021.pdf` (local only). Data:
`data/marek_tacina1975/`.

## For the reviewer

1. **Check Table I against the page** (printed p. 14, PDF p. 15): 19 runs by 9
   stations, read off the page image at 200 dpi and checked again at 400 dpi.
   Review page (scans beside the CSVs): https://claude.ai/artifact/UnHJWCpUHSVdYmvMv2Wj3A
   - Four cells carry handwritten pen marks over the printed digits. The
     printed values are used: 0%/9.97% run at 0.95 cm, 306 (pen "298"?);
     the same run at 21.28 cm, 389; the 72%/7.54% run at 8.58 cm, 398; the
     72%/10.00% run at 8.58 cm, 364.
   - The 72% plate's 16.20 cm station is printed "---" in all five runs (blank
     in the CSV).
   *What turns on it:* the fidelity score of the slot film kernel, stage 1 of
   #714 (C_m from Tu, D-7958 Eq. 5).
2. **Check the per-blockage averages** (printed p. 8: C_m, PDF p. 9) and the
   turbulence intensities (printed p. 9, PDF p. 10, and p. 6): 7, 14, 23 and
   35% axial Tu for 0, 52, 72 and 90% blockage.

## What the report gives

- **Rig (p. 2-3):**
  - 10.5 cm high by 37 cm wide duct at atmospheric pressure;
  - hot gas at about 590 K and 62 m/s, Re 189 000;
  - a 0.63 cm slot across the full width, cooling air at ambient;
  - a 0.31 cm stainless 304 test plate, insulated underneath;
  - blockage plates 41 cm upstream: 2 rows of 7 holes, mesh 5.3 cm.
- **Definitions:**
  - eta = (T_H - T_f)/(T_H - T_s) (Eq. 1), with T_f the measured wall
    temperature (the authors' assumption: an adiabatic wall);
  - M = rho_s U_s / (rho_H U_H), which for this rig is (m_s / m_H)(H / S);
  - Juhasz-Marek: eta = 1/(1 + C_m x/(M S)) (Eq. 4);
  - C_m = (0.25 +- 0.15) Tu_u/100 (Eq. 5).
- **Measurement band:** the authors estimate a maximum eta deviation of 0.04
  (p. 6: +-5 K hot-gas variation plus up to 7 K axial wall conduction). That
  is a measurement band, usable as the fidelity band.
- **Caveats (the authors'):**
  - the 21.28 cm station reads low through conduction to the plate supports
    (p. 7);
  - Tu was measured cold at the same mass flow (Re 334 000 rather than 189 000);
  - the turbulence decays along the plate, and the averages are over the test
    area.

## Checks that could have failed

- **Cooling percent:** the printed percent is 100 m_s/(m_s + m_H) to within 0.043 (the printing's rounding).
- **Wall temperature:** at every blockage and station the wall cools as the slot flow rises.
- **Eq. 4 on Table I reproduces the authors' own average C_m:** 0.026, 0.037, 0.046 and 0.079 against the printed 0.025 +-0.01, 0.04 +-0.02, 0.05 +-0.03 and 0.1 +-0.05. This confirms the transcription and the M definition together.

## Scored (2026-10-25, #714)

`slot_film_effectiveness_marek_tacina_1975` (`include/cooling_correlations.h`)
implements Eqs. (4) and (5). `python/tests/test_slot_film.py` scores every
cooled, printed Table I cell (140):

- **With the authors' per-blockage C_m (p. 8):** 136 of 140 inside 30%,
  the worst 33%. That reproduces the report's "+-30 percent" for this case.
- **With C_m from Eq. (5) at the measured Tu:**
  - eta RMS 0.0696, bias +0.006;
  - 57 of 140 inside the +-0.04 measurement band;
  - by blockage the bias is +0.068 / +0.017 / -0.038 / -0.029 (0/52/72/90%).

  Eq. (5)'s single slope under-mixes at low Tu and over-mixes at high Tu. That
  is the line's fit to four averages, not a kernel error.
- **Not scored:** the uncooled runs (eta undefined) and the user's
  `C_m_multiplier`.

## Review log

- 2026-10-25: transcribed from the page images; the four data tests pass.
- 2026-10-25: the user checked Table I, the C_m table and the Tu values
  against the scans and signed them off, the printed digits under the four
  pen marks included.
