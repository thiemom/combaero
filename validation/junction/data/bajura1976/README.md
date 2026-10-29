# Bajura & Jones (1976), Figs. 10 and 11: manifold header pressures

Bajura, R. A., Jones, E. H., Jr. (1976). Flow distribution manifolds. ASME J. Fluids Eng. 98(4), 654-665 (76-FE-7).

The figures show measured static pressures along a dividing (inlet) header and a combining (exhaust) header joined by 10-20 laterals, plotted against x/L over the last 80 % of the header. They also show the header-to-header differential. Pressures are normalised by the paper's reference pressure difference. The rig and model coefficients are in the paper's Fig. 5 and Table 4 (see `validation/sources.yaml`).

## Files

One CSV per figure and series, with columns `x_over_L` and `p` (dimensionless):

| Series | Fig. 10 points | Fig. 11 points |
|---|---|---|
| `dividing_header` (squares) | 17 | 17 |
| `combining_header` (circles) | 16 | 16 |
| `differential` (triangles) | 16 | 17 |

## Extraction

- **Method:** an automated pass plus hand readings on 2026-10-28.
  - **Hand-read in Fig. 10:** squares at x/L 0.62, 0.75 and 1.0; circles at 0, 0.19 and 0.44; triangles at 0.31, 0.44 and 0.62.
  - **Hand-read in Fig. 11:** squares at 0.125, 0.5, 0.75 and 0.81; circles at 0 and 0.68; triangles at 0.375, 0.44, 0.5, 0.625, 0.93 and 1.0.
- **Confirmation:** the user checked every point on a scan overlay on 2026-10-29.
- **Fig. 10 calibration:**
  - The y tick printed ".8" is a misprint for 0.9: it sits at 296.0 px where 0.9 predicts 295.6 px. Only the 1.2, 0.4 and 0 ticks are used.
  - With that calibration, both headers' x = 0 point comes out at 1.00, which is the normalisation.
  - The differential then equals dividing minus combining within 0.014 at all 15 stations.
  - The print has no circle at x/L 0.375. At x = 0 the differential coincides with the dividing point and is plotted once.
- **Fig. 11 calibration:** all four y ticks are linear.

## Known inconsistency in the source (Fig. 11)

Beyond x/L of about 0.6, the plotted differential falls up to 0.18 below dividing minus combining. For example, at x/L 1.0 the dividing header reads 1.24 and the combining header 0.04, but the differential reads 1.02 instead of 1.20. The user checked the overlay on 2026-10-29: nothing is misread. The values as printed do not add up, and the plot gives no explanation.

Treat this as a property of the source. If Fig. 11 is used, score the headers and the differential separately and do not derive one from the others. Fig. 10 is self-consistent.

## Status

Extracted and unscored. These figures are for an accuracy (cross-source) check of a network of N tees with N orifice laterals against the measured header pressures (#779, #812, #813). They test no single tee model's fidelity. The paper does not measure the lateral flows.
