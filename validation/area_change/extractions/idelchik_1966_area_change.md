# Idelchik (1966) -- sudden expansion and conical diffuser

**Status: CONFIRMED (transcribed 2026-10-24, signed off by the user
2026-10-25, #731).** Source: I. E. Idel'chik, *Handbook of
Hydraulic Resistance*, 1st English edition (AEC-tr-6630, 1966),
`sources/junction/Idelchik.pdf` (local only). Data: `data/idelchik1966/`.

## For the reviewer

1. **Check the three tables against the pages:** diagram 4-1 (printed p. 128,
   PDF p. 134), parts 1 and 2, and diagram 5-2 (printed p. 167, PDF p. 173).
   Every number is a printed table entry, read off the page image at 220
   dpi. The text layer is too garbled to use. Nothing is digitised from a
   curve.
   Review page (scans beside the CSVs): https://claude.ai/artifact/F4Jh3brZt96RHL2CqDo3Gi
   *What turns on it:* the fidelity score of `sharp_area_change` (expansion
   branch) and `conical_area_change` (diffuser branch), #731 items 1 and 3.
2. **One header was corrected on reading:** diagram 4-1 part 2's last column
   is Re = 3.5e3 (printed "3.5·10³"); the text layer reads "3 5-103". Its
   values are part 1's Borda-Carnot ones, which the data test checks.

## What the diagrams give

- **4-1, part 1 (Re = w0 Dh / nu >= 3.5e3, uniform inlet profile):**
  zeta = (1 - F0/F2)^2 without baffles, ~0.6 (1 - F0/F2)^2 with baffles,
  referenced to the small section's velocity head. 10 area ratios.
- **4-1, part 2 (10 < Re < 3.5e3):** zeta(Re, F0/F2), F0/F2 0.1-0.6 by 13
  Reynolds numbers. Below Re 8 the page gives zeta = 26/Re (not
  tabulated).
- **5-2 (conical diffuser in a line, uniform inlet profile):** zeta_exp =
  phi_exp (1 - F0/F1)^2, tabulated directly as zeta_exp(F0/F1, alpha), where
  alpha is the FULL included angle. 11 ratios (0-0.6) by 14 angles (3-180
  deg). The total is zeta_exp + zeta_fr; friction is NOT in the table.

Not transcribed: diagram 3-9 (sharp contraction). Its single printed
formula, zeta = 0.5 (1 - F0/F1), Re > 1e4, is exactly the code's. Diagrams
3-5/3-6 (conical contraction, zeta'(alpha, l/D)) are follow-ups.

## Checks that could have failed

- 4-1 part 1 is (1 - n)^2 to two decimals, and 0.6x that with baffles (all
  10 ratios).
- 4-1 part 2's Re 3.5e3 column equals part 1, its Re 10 column is 3.10 at
  every ratio, and at every other Re a larger ratio loses less.
- 5-2: at every angle a larger ratio loses less; every row rises with angle
  to a peak at 60 deg; the 180 deg column sits 0.000-0.024 above
  Borda-Carnot. The data test first had 0.021 as that bound, from rows
  checked by eye, and failed on F0/F1 0.075 (0.88 printed, 0.856 Borda). The
  page confirms 0.88, so the bound was widened to the measured range.

## Scored (2026-10-25, #731)

Fidelity: the code now implements these diagrams (`include/area_change.h`),
and `python/tests/test_area_change.py` scores every printed cell.

- **Sudden expansion, 4-1:** 75 of the 78 part-2 cells exact. Three move:
  - n 0.6 and 0.5 at Re 2000, lowered to 0.41 (printed 0.50) and 0.61 (0.65), so the loss rises with the flow. The printed n 0.6 column falls between Re 2000 and 3000, which would give one pressure drop two flows;
  - n 0.3 at Re 3500, printed 0.50 for Borda-Carnot's 0.49.
  Part 1 is Borda-Carnot exactly, and below Re 8 the law is 26/Re.
- **Conical diffuser, 5-2:** all 154 cells exact (bilinear through the table).
- **Before (2026-10-24, the replaced closures):**
  - expansion, 10 < Re < 3.5e3: 72 cells, bias -50%, RMS 55%; above 3.5e3, alpha - 1 added to Borda-Carnot (+0.069 at Re 1e4, even at F0 = F1);
  - diffuser: model/Idelchik median 3.0x at alpha <= 10 deg, 1.55x at 12-30 deg, 0.84x at >= 40 deg.

## Review log

- 2026-10-24: transcribed from the page images; the three data tests pass.
- 2026-10-25: the user checked the three tables against the pages and signed
  them off.
- 2026-10-25: the closures implemented and scored (#731, above).
