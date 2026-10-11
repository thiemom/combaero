# Mynard & Valen-Sendstad 2015: Figs 4 and 6-11

**Source:** J P Mynard, K Valen-Sendstad, *"A unified method for estimating
pressure losses at vascular junctions"*, Int J Numer Methods Biomed Eng
31(7):e02717, 2015, DOI 10.1002/cnm.2717. This is the accepted manuscript
(University of Melbourne repository), `docs/junction/a-unified-method-...pdf`,
which is gitignored for copyright. The digitised coordinates are committed.

## Files

| File | What it is |
|---|---|
| `mynard_figNNx_unified0d.csv` | Centre points of the red Unified0D curve, Mynard's model |
| `mynard_fig04_unified0d_eta0.csv` | Fig 4's dashed curve, his model with eta_j = 0 |
| `mynard_figNNx_ref3d.csv` | Ref3D markers, his 3D CFD |
| `calibration.yaml` | Per panel: pixels per data unit, the calibration's worst tick misfit, and the point counts |

In every CSV, `x` is the panel's abscissa: Reynolds number, flow fraction,
angle in degrees, or area ratio. `K` is Mynard's Eq 8 total-pressure loss over
the common branch's dynamic head, which is Bassett's K.

All files are written by `../../digitise_mynard2015.py`; do not edit them.
Regenerate them with that script, which needs `pdfimages` and the PDF.

## What to know before using them

- **Raster source.** The figures are 150 dpi, about 200 px per panel. One
  pixel is 0.005-0.017 in K on most panels and 0.36 on Fig 6d's -20..40 axis.
- **The x calibration from tick labels is up to 2.2 px off** (Figs 4, 6d, 7).
  `mynard_fidelity.py` re-fits x from the markers, which sit on design
  abscissae: flow fractions in steps of 0.05, angles in steps of 5 deg, and
  area ratios (k/6)^2. The insets print 0.25 and 0.44.
- **Hidden curve.** Where the Ref3D line covers the red one, no red point is
  recorded, rather than a biased one. The Re panels keep 15-160 points.
- **The red curves are splines** through the model evaluated at the CFD
  samples, not dense evaluations: on Fig 6d the curve leaves the model between
  samples. Score at the samples, or spline the model the same way.
- **Ref3D is laminar CFD, Re 350-2400.** The inlet profile is flat with a
  boundary layer (gamma = 9), and the corners are smoothed. It has no stated
  uncertainty.
- **Fig 9's rows are swapped against the caption.** The top row is the
  straight inlet, the bottom row the side inlet. See
  `../../mynard_fidelity.py`.
