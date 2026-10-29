# Bassett 2001 — T-junction loss coefficients

**Source paper:** M D Bassett, D E Winterbone, R J Pearson,
*"Calculation of steady flow pressure loss coefficients for pipe junctions"*,
Proc IMechE Part C, 215(8):861–881, 2001.

PDF in `sources/junction/bassett.pdf` (gitignored — copyright).

## Conventions

- `q = m_other / m_com` (Table 1). For separating flows (types 1–3) `m_other`
  is one of the outlet branches; for joining flows (types 4–6) it is one of
  the inlet branches.
- `psi = F_C / F_B` (area ratio: main / lateral).
- `theta` is the lateral branch angle in degrees in filenames, radians at the
  Python API boundary. `psi = 1` = equal area.
- Coefficient ids `K1`..`K12` match Bassett Table 1 / Fig 2.
- `_measured` = digitized data points from the figure.
- `_calc` = Bassett's own published curve (not always the same as the Table 2
  analytical evaluation — some figures use Eq 33/34/35/37 angle corrections).

## File naming

`bassett_figXX[subfig]_Kid_<param=val>..._<measured|calc>[_extra].csv`

Example: `bassett_fig10c_K12_theta=45_psi=3_measured.csv` is K12 measured at
theta=45 deg, psi=3, from Fig 10c.

## Fit-tier ratings (per coefficient)

Per the paper's own measured-vs-calculated comparisons in §4:

| K   | Tier      | Notes                                                |
|-----|-----------|------------------------------------------------------|
| K2  | excellent | theta- and psi-independent; tight fit                |
| K5  | excellent | same form as K2                                      |
| K6  | very good | Fig 7a/b: ~5–10% spread vs measured                  |
| K1, K3, K4 | very good | similar derivation to K6                       |
| K12 | good      | Fig 10a/b/c: good with Eq 33 angle correction        |
| K11 | moderate  | Fig 11: predictions fall below measurements          |
| K7, K8 | moderate | Fig 12/13: only equal-area data available           |
| K9, K10 | poor    | Fig 14: paper's own text calls fit "poor"            |

A candidate junction model that fails to match measured K9 is NOT necessarily
broken — it's failing where 1D physics itself fails. The scorecard's
`Delta_vs_ceiling` metric encodes this.

## Verdicts from the runner (vs Fig 10c, Fig 12b measured)

- **K12 in `include/tee_junction.h` is CORRECT.** An earlier suspicion of a
  ψ²-vs-ψ transcription bug was wrong — the correct Bassett formula has
  ψ (not ψ²) inside the bracket. Verified at psi=3, q=0.5: C++ K12 gives
  1.47, measured K12 from Fig 10c gives 1.46.
- **K8 in `include/tee_junction.h` had a transcription error, since fixed.**
  It divided by psi where Bassett Table 2 multiplies: Table 2 prints
  `K8 = 1 - q^2 + 2 (1-q)^2 psi cos(theta)`, which is what the code and
  `models/bassett2001.py` now carry (re-read from the PDF's text layer
  2026-10-10).

## Fidelity of the transcription

`python/tests/test_junction_mynard_fidelity.py` holds `models/bassett2001.py`
to Bassett's own calculated curves (the `_calc` files): max error 0.012-0.024
at psi = 1, about 0.1 at psi = 3 / 4 for K12 (Fig 10c).

K4-K8 and K11 (Figs 7, 8, 11-13, 32 curves) were traced on 2026-10-11 from
the paper's 300-dpi images. Each panel is calibrated on its axis lines and
ticks (residual <= 3.4 px), and each curve is traced by snapping to the ink
within 8 px of our formula, every 0.02 in q, skipping q where another curve
of the panel is close. The committed points are therefore within 0.031 of
our formulas by construction; the evidence that the formulas are Bassett's
is that the unsnapped ink distance is 0.5-3 px along the curves, whereas a
re-planted K8 /psi typo sits 12-18 px off and K11 raw-for-corrected 6-21 px.
K7, K8 and K11 use the Eq 33/34 angle correction, as the captions state.

Bassett's psi > 1 curves sit where our formulas put an area ratio 1-3% off
nominal (Fig 11b "3" fits psi 2.91, Fig 13b "4" fits 4.08, Fig 7b "3" fits
2.97), at most about 0.05 in K: the ratios he computed with, not a
transcription error.
