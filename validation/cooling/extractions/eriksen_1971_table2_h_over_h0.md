# Eriksen (1971) Table II -- gas-side h/h0 behind a row of 35 deg holes (#476)

**Status: CONFIRMED.** The user signed off the transcription on 2026-11-03.

Source (local only, gitignored):

- V. L. Eriksen, *Film cooling effectiveness and heat transfer with
  injection through holes*, NASA CR-72991 (HTL TR 102), University of
  Minnesota (1971); NTRS 19720007295, open access,
  `sources/heat_transfer/film/Eriksen_1971_NASA_CR-72991_holes_eta_and_h.pdf`.

Data:

- `data/eriksen1971/table2_blockA_as_printed.csv`: Table II, printed p. 160
  (PDF p. 179), laid out as printed.
- `data/eriksen1971/table2_blockB_as_printed.csv`: Table II "Continued",
  printed p. 161 (PDF p. 180), laid out as printed.
- `data/eriksen1971/table2_ha_over_h0_long.csv`: both blocks in long form,
  derived mechanically from the as-printed files.

## What is tabulated

- **Quantity.** ha/h0, dimensionless.
  - ha = q / (Tw - Taw) (eqn 10), averaged laterally over one hole pitch,
    Z/D -1.5 to +1.5 (eqn 11). For unheated jets Taw = Tinf.
  - h0 is h measured at the same X and free-stream velocity with the holes
    taped over (report pp. 41-42). It is not the Table I fit, and its
    absolute value is not tabulated.
- **Geometry.** One row of 5 tubes, inclined 35 deg in the streamwise
  direction. D = 1.18 cm, pitch 3 D, L/D about 62. X is measured from the
  hole's downstream edge. Constant wall heat flux; the boundary layer is
  tripped turbulent.
- **Block A.** Re_D 0.44e5 (about 61 m/s), delta*/D 0.149. Seven runs at
  M 0.10-0.99: four unheated (density ratio 1.0) and three heated (about
  55 K, density ratio about 0.85).
- **Block B.** About 30.5 m/s, delta*/D 0.171. Five unheated runs at
  M 0.20-1.94.
- **Stations.** 12, at X/D 6.65-80.46. ha/h0 ranges from 0.940 to 1.268.

## Print error

Block B's printed header says Re_D = 0.44e5, but the block is the
report's Fig. 68, at Re_D = 0.22e5. The figure caption and legend give the
same delta*/D and M values, and the text (report pp. 59-60) agrees. The
as-printed CSV keeps the printed header; the long CSV carries both
`Re_D_printed` and `Re_D_used`.

## Read procedure

- Two independent reads from page images: 300 dpi row-wise, and 400 dpi
  column-wise.
- **0 disagreements** over 180 value cells, the run header rows and the 24
  X/D labels.
- Dust specks beside decimal points were checked on crops.
- Plausibility (not proof): block B matches Fig. 68 by eye. In block A,
  heated minus unheated runs differ by -0.4 / -0.7 / +2.0 % at M 0.2 /
  0.5 / 1.0. The text says the heated runs read about 1 % and 1.5 % lower
  at M 0.2 and 0.5.

## Measurement band (Appendix A)

Table A.III gives h uncertainties on the centreline only: h0 1.6 %,
unheated injection 2.0 %, and heated injection 11.9-15.6 % at X/D 6.65,
falling to 1.8-2.0 % from X/D 37. Combined by root-sum-square, ha/h0 is
known to about +-2.6 % for unheated runs. Heated runs are about +-12-16 %
at X/D 6.65 and about +-2.5 % downstream. These are measurement bands,
never derived from any model's error.

## Use

The table is a cross-source accuracy set for a gas-side augmentation
closure (#476).

- Score h/h0 against M and X/D on the 9 unheated runs (108 points).
- Exclude or label the heated X/D 6.65 points.
- Limits: one 35 deg row at P/D 3, not an effusion array (Caciolli's plate
  is P/D 6-8 with 29 rows); flat plate; low turbulence.
- The ratio must be scored as a ratio: absolute h0 is not available.
