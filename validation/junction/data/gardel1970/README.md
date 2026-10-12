# Gardel & Rechsteiner 1970 -- unsplit-flow tables E, F, G

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

## Not transcribed here

- **Figs 12-15:** the raw points, four shapes.
- **Figs 18-27:** smoothed curves for all 33 shapes. These are curves through
  the data, not formulas, but the scatter is gone.
- **Figs 28-34:** replot tables E-G and add curves in phi.
