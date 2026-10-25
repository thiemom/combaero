# Idelchik (1966) -- conical and low-Re contraction

**Status: CONFIRMED (transcribed and signed off by the user 2026-10-25,
#731).** Source: I. E. Idel'chik, *Handbook of
Hydraulic Resistance*, 1st English edition (AEC-tr-6630, 1966),
`docs/junction/Idelchik.pdf` (local only). Data: `data/idelchik1966/`.

## For the reviewer

1. **Check the two tables against the pages:**
   - diagram 3-6: printed p. 96, PDF p. 102;
   - diagram 3-10: printed p. 99, PDF p. 105.

   Every number is a printed table entry, read off the page image at 200 dpi
   and checked again at 300 dpi. Nothing is digitised from a curve.
   Review page (scans beside the CSVs): https://claude.ai/artifact/E5Vzv5M2iVHuaWzwupSxVh
   *What turns on it:* the conical nozzle of `conical_area_change` and
   `area_change_flow`. Today that closure is the unsourced
   C_con = 0.04 + 0.46 sin(theta)^1.5. The contraction would also gain a
   low-Re law.
2. **Diagram 3-10's last two columns are both headed "10^4"** (the image is
   clear on this). The last one prints 0.5 (1 - n) at every ratio: 3-9's
   high-Re formula. So I read it as the Re >= 1e4 limit. The other differs
   from it only at n 0.1 (0.50 against 0.45). The CSV keeps the printed
   column number (`column`, 1 to 14), with Re 1e4 on both.

## What the diagrams give

- **3-7, inlet with smooth contraction, Re > 1e4 (printed p. 96):** for a
  rectilinear converging bellmouth, which is the in-line conical contraction
  from F1 to F0,
  zeta ~= zeta' (1 - F0/F1) + zeta_fr. Here zeta' is "determined as zeta from
  diagram 3-6", and zeta_fr is the cone's friction (diagrams 5-2 to 5-5). It
  is a formula on the page, not a table, so it is recorded here only. The
  curvilinear bellmouth is friction only.
- **3-6, conical converging bellmouth with end wall, Re > 1e4:** zeta'(l/D_h,
  alpha). alpha is the FULL included angle and l the cone's axial length.
  The table has 6 lengths (0.025 to 0.60) by 9 angles (0 to 180 deg). It is
  "approximately" (the page's word). Both ends are the sharp edge's 0.50.
- **3-10, sudden contraction, inlet in the end wall, Re < 1e4 (Karev):**
  zeta(Re, F0/F1), with 6 ratios (0.1 to 0.6) by 14 printed columns (Re 10 to
  1e4). Below Re 8 the page gives zeta = 27/Re, which is not tabulated. As in
  4-1, a hump sits between Re 2e3 and 5e3, and at n 0.6 zeta Re^2 falls
  there, so a closure needs the same monotonicity check before it ships.

Not transcribed:
- diagram 3-5 (conical bellmouth *without* end wall): 3-7 points to 3-6;
- the curves of 3-6 and 3-10: their tables are the data.

## Checks that could have failed

- **3-6:** the 0 and 180 deg columns are 0.50 at every length. At every
  angle a longer cone loses no more, and every value lies in (0, 0.5].
- **3-10:** the last column is 0.5 (1 - n) to two decimals, Re 10 prints
  5.00 at every ratio, and at every Re a larger ratio loses no more.

## Review log

- 2026-10-25: transcribed from the page images; the two data tests pass.
- 2026-10-25: the user checked both tables against the pages and signed them
  off, the high-Re reading of 3-10's last column included.
