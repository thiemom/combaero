# Gardel & Rechsteiner 1970 -- tables E, F, G and the raw points of Figs 12-15

Gardel, A. and Rechsteiner, G. F. (1970), *Les pertes de charge dans les
branchements en Te des conduites de section circulaire*, Bulletin technique de
la Suisse romande 96(25), pp. 363-391, doi:10.5169/seals-70882 (e-periodica,
ETH Library). The scan is copyrighted and gitignored
(`docs/junction/gardel_1970.pdf`).

**Not a reprint of Gardel 1957**: a new EPFL campaign that supersedes it.
- **Geometry:** 33 shapes, all water at D = 150 mm in asbestos-cement pipe.
  - Branch angle delta 45-135 deg (angle to the UPSTREAM main alpha).
  - Branch/main area ratio phi = (d/D)^2 = 1, 0.69, 0.44, so psi = 1/phi <= 2.27.
  - Edge fillet rho = r/D = 0, 0.1, 0.2, with the fillet made to drawing this time.
- **Flow:** four modes and more than 6000 flow cases.
- **The authors discount** their own 1957 series (fillets ill-defined), Munich
  (Vogel, Petermann, Kinne: 43 mm, near-laminar at low flows) and Blaisdell &
  Manson (13-64 mm pipes, joining only) as too scattered (pp. 369, 387-389).

## `tables_EFG.csv`

Transcribed as printed: one CSV row per printed row, one column per angle.

| column | meaning |
|---|---|
| `table` | E (q_beta = 0, straight through), F (q_beta = +1, all of alpha into the branch), G (q_beta = -1, all of the branch into gamma) |
| `q_beta` | the printed left index of h |
| `quantity` | `h_ya` = h_gamma,alpha; `h_ba` = h_beta,alpha; `h_yb` = h_gamma,beta; `h_ab` = h_alpha,beta (see below) |
| `author` | GR = this study; G = Gardel 1957; VPK = Vogel/Petermann/Kinne, Munich 1926-31; BM = Blaisdell & Manson 1963/67 |
| `phi`, `rho` | as printed; `~0.1` = Gardel 1957's approximate fillet |
| `d45` ... `d135` | h at that angle; empty = printed dash |

**Cell markers, kept from the print:**
- `*` = "voir fig. 16": the 90 deg, rho = 0 case drawn as a flow sketch.
- `(...)` = printed in parentheses, which the table legend calls "aberrant". The
  conclusions (p. 391) also call the -3.71 at phi 0.44, rho 0, 120 deg
  "certainly wrong". Exclude these from any score.
- `a;b` = printed as two values (VPK at 90 deg in table G, "-0.9, -1.1" and
  "-0.31, -0.56"). Not a single measurement; exclude.

**Two rows are kept as printed although they look swapped.** In table G the
1957 rows list rho 0.03 for phi 0.28 and 0.04 for phi 0.16. Tables E and F
print the reverse.

### Definitions (pp. 367-369)

- **Branches:** alpha is the upstream main, gamma the downstream main, beta the
  branch.
- **Head:** H is static head plus V^2/2g of the mean velocity, with no
  velocity-profile correction (eps = 0).
- **Where it is evaluated:** each branch's head line is fitted where the flow
  has recovered, with the undisturbed-pipe friction gradient removed, and
  extrapolated to the intersection of the three axes.
- **Differences:** H_ij = H_i - H_j.
- **Normalisation:** `h = H / (V^2/2g)`, with V the TOTAL flow over the
  base-pipe area pi D^2/4. In all three unsplit cases that is the velocity of
  the common (full-flow) branch, so h is on Bassett's common-branch basis.

**Signs.** A negative h is a head loss in the direction of flow:

| table | loss coefficient on the common-branch head |
|---|---|
| E, `h_ya` | -h: straight-run loss at zero branch flow |
| F, `h_ba` | -h: dividing, all flow into the branch (K for q = 1) |
| G, `h_yb` | -h: joining, all flow from the branch (K for q = 1) |
| E `h_ba`, F `h_ya`, G `h_ab` | dead-branch head relative to the live one; a pressure tap, not a loss |

The last row is the head a dead port sits at, which is directly comparable to
where a closure puts a port carrying no flow. (In table G the printed
quantity is `-1 h_alpha,beta = h_alpha - h_beta`, sign flipped by the authors
"by analogy" with table F.)

**What the print states about uncertainty:** no numeric band. The authors
report the GR series' dispersion qualitatively, and in their own
fully-separated cases it is small. Under `docs/VALIDATION_POLICY.md` the band
is therefore null.

### How values were read

- **Page images:** rendered at 400 dpi with `pdftoppm` and read by eye.
- **Cross-check:** against the PDF text layer, compared per table part as
  multisets of numeric tokens.
  - Every disagreement was a text-layer misread (`-0.8S`, `-1.B9`, `-n.27`,
    split `-0.8 3`, dashes glued to the next phi), with two exceptions:
    - one dropped cell, table F `h_ya` VPK phi 1 rho 0 90 deg = -0.35, read
      on two crops;
    - one real conflict, table E `h_ba` GR phi 0.44 rho 0 120 deg, where the
      text layer has -0.86. On a 3x zoom of the page it is -0.88.
- **Scripts:** the extraction scripts were throw-away and are not in the repo.

## `detail_points.csv`

The measured points of the four shapes the paper plots point by point, digitised
from the scan. One row per point.

| fig | delta | phi | rho |
|---|---|---|---|
| 12 | 90 | 1 | 0 |
| 13 | 90 | 1 | 0.2 |
| 14 | 90 | 0.69 | 0 |
| 15 | 135 | 1 | 0 |

| column | meaning |
|---|---|
| `mode` | flow mode as numbered in the paper: 1 dividing (alpha feeds beta and gamma), 2 joining (alpha and beta feed gamma), 3 beta feeds alpha and gamma, 4 alpha and gamma both flow into beta |
| `abscissa`, `q` | the figure's own abscissa: `q_beta` in modes 1 and 2 (0..1 and -1..0), `q_gamma` in modes 3 and 4 (0..1 and 0..-1) |
| `series` | `h_ba`, `h_ya`, `h_yb`, as in the tables; same normalisation and signs |
| `h` | the value read off the plot, to 3 decimals |

**Layout.** Each figure is one strip of four panels, left to right mode 3, 2,
1, 4. Adjacent panels share their end, an unsplit case: q_gamma = 1 / q_beta =
-1 (table G), q_beta = 0 (table E), q_beta = +1 (table F). A point exactly on
a shared end is listed once, under the left panel.

**How the points were read.**
- **Renders:** the figure pages at 500 dpi.
- **Calibration:** the scan is rotated and warped by up to 14 px, so a linear
  axis fit is biased by up to 0.03 in h. Every printed grid line was traced
  across the figure, and each point is mapped by interpolating between the
  traced lines.
- **Detection:** an open-ring detector (every angular sector of the annulus
  dark, the centre white once grid lines are removed).
- **Series:** assigned by tracking each curve across a panel. Fig. 15 needed
  rough guide curves where the three series cross with sparse points.
- **Manual pass, Fig. 15:** its rings are printed fainter, and the detector
  missed about 27. They were clicked by hand, with the series and one
  relabelling decided by eye. 24 were snapped to the ring centre and read
  through the traced grid; 3 rings too faint to snap keep the click's value
  (calibration residual 0.004).
- **Review:** every point was checked on an overlay of the scan.

**Checks.**
- **Panel ends against tables E/F/G:** 36 comparisons, 33 within 0.015, the
  largest 0.029 (Fig. 14, q_beta = +1).
- **Head identity:** `h_ya = h_yb + h_ba` holds in every flow case. On cases
  with one point of each series at the same abscissa:

  | fig | cases | MAD of residual | over 0.05 |
  |---|---|---|---|
  | 12 | 84 | 0.0026 | 0 |
  | 13 | 92 | 0.0030 | 0 |
  | 14 | 80 | 0.0027 | 0 |
  | 15 | 69 | 0.0024 | 0 |

  The identity cannot see a swap of `h_ba` and `h_yb` (the sum is the same);
  those were checked by eye on the overlay.

**Known limits.** Heavily overplotted rings can merge into one, so a cluster of repeated
measurements yields fewer rows than the paper took. No uncertainty is stated,
so the band stays null.

## Not transcribed here

- **Figs 18-27:** smoothed curves for all 33 shapes. These are curves through
  the data, not formulas, but the scatter is gone.
- **Figs 28-34:** replot tables E-G and add curves in phi.
