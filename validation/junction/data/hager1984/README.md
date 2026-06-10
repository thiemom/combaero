# Hager 1984 — Branch and bend loss coefficients

**Source paper:** W H Hager, *"An approximate treatment of flow in branches and bends"*,
Proc IMechE Part C, 198C(4):63–69, 1984.

PDF in `docs/junction/hager1984.pdf` (gitignored — copyright).

Foundational paper introducing the `(3/4)*delta` angle correction adopted by
Bassett. Considers equal-area T-junctions only (`psi = 1` implicit). The
straight-through and lateral loss correlations reduce to Bassett K2/K5 and
K6 respectively at `psi = 1`, with the q-axis convention difference noted
below.

## Conventions

- `q = DeltaQ / Q` = lateral flow / total flow. **Different from Bassett K5**,
  where `q = m_straight / m_total = 1 - q_Hager`. The runner must apply
  `q_Bassett_K5 = 1 - q_Hager` when comparing Hager Fig 4 against Bassett K5.
- `delta` is the lateral branch angle (degrees in filenames, radians at API).
- `xi_t` = straight-through loss (Hager Eq 8: `xi_t = q*(q - 1/2)`).
- `xi_l` = lateral loss (Hager Eq 19: `xi_l = 1 - 2*q*cos((3/4)*delta) + q^2`).

## Files

- `hager_fig04a_xi_t_delta=90_single_measured.csv`: Hager's own measurements
  for single-lateral T-junctions at delta=90 deg.
- `hager_fig04a_xi_t_delta=90_manifold_measured.csv`: Hager's measurements
  for manifold-lateral T-junctions at delta=90 deg.
- `hager_fig04b_xi_t_envelope=mean.csv`: Fig 4b mean curve of the overall
  data range for arbitrary branch geometry.
- `hager_fig04b_xi_t_envelope=upper.csv`, `..._envelope=lower.csv`: range
  bounds. Defines an uncertainty envelope for `xi_t` independent of paper
  source.

## Why no Fig 5 (lateral)

Hager Fig 5b plots `xi_l` measured at delta=45 deg and 90 deg using sources
Mock, Gardel/Rechsteiner, Miller, and Ito et al. — the same sources Bassett
cites for his Fig 7a K6 measured points. Bassett Fig 7a covers delta in
{45, 60, 90, 120} deg, fully superseding what would be the digitized Fig 5b.
The redundancy is documented; only Fig 4 is included here.

## Use in the scorecard

Hager provides an **independent** ground truth for the straight-through loss
coefficient. Any candidate junction model that reproduces Bassett K5 must
also reproduce Hager xi_t (with q-axis flipped) within Fig 4b's envelope.
The runner reports a single `pct_within_hager_envelope` metric for the
straight-loss coefficient.
