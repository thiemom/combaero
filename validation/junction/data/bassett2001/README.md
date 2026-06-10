# Bassett 2001 — T-junction loss coefficients

**Source paper:** M D Bassett, D E Winterbone, R J Pearson,
*"Calculation of steady flow pressure loss coefficients for pipe junctions"*,
Proc IMechE Part C, 215(8):861–881, 2001.

PDF in `docs/junction/bassett.pdf` (gitignored — copyright).

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

## Known issues caught here

- `include/tee_junction.h` K12 in main has a transcription bug: the bracket
  uses `q^2 * psi * cos(theta)` instead of the correct `q^2 * psi^2 * cos(theta)`
  per Bassett Table 2. Invisible at psi=1 (ψ = ψ²) — Fig 10c psi=3 measurements
  resolve it.
- K8 in tee_junction.h uses `(1-q)^2 / psi` where Table 2 has `(1-q) * psi`.
  Fig 12a equal-area measurements resolve it.
