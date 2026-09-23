# Extraction: Han jet impingement correlations

**Status: CONFIRMED. I1, I2, I3, item 13, item 20 and item 21 all resolved (2026-09-20/21) against the primary paper (Florschuetz, Truman and Metzger, 1981), checked page-image by page-image since it is a scanned copy. Cleared for implementation.**

Every equation below has been checked against the printed page image (not just
the text layer, which garbles at least one equation -- see item 3). Nothing
here is released for implementation until a reviewer signs off, per the
project's provenance discipline.

Tracked by #337, under #339.

Source, pinned:

> Han, J., Dutta, S. and Ekkad, S. (2012). *Gas Turbine Heat Transfer and
> Cooling Technology*. 2nd Edition. CRC Press. Section 4.1, pp. 329-355.

A secondary source, as throughout this project. Two primary papers do the
actual work:

> Florschuetz, L.W., Truman, C.R. and Metzger, D.E. (1981). Streamwise Flow
> and Heat Transfer Distributions for Jet Array Impingement with Crossflow.
> *ASME J. Heat Transfer*, **103**, 337-342.
>
> **Obtained and checked directly against the source (2026-09-21).** A scan
> of the reviewer's own paper copy, page images rendered at 200dpi and read
> equation-by-equation -- see `docs/heat_transfer/flohrshuetz/` (gitignored,
> copyrighted). Being a scanned copy, its OCR text layer is noticeably less
> reliable than Han's born-digital PDF -- e.g. it interleaves the two-column
> Nomenclature block into nonsense -- so every equation and table entry below
> attributed to this source was confirmed against the rendered page image,
> not the text layer alone.
>
> Goldstein, R.J., Behbahani, A.I. and Heppelmann, K.K. (1986). Streamwise
> Distribution of the Recovery Factor and the Local Heat Transfer Coefficient
> to an Impinging Circular Air Jet. *International Journal of Heat and Mass
> Transfer*, **29**(8), 1227-1235. (page/volume as cited by Han; not yet
> checked against the paper itself)

---

## For the reviewer

One item needs judgement before this can move to implementation -- the
rest are resolved facts, not choices.

| # | what to check | why it matters |
|---|---|---|
| **I1** | **RESOLVED, see below.** Table 4.1's geometry box is not Florschuetz's validity range at all -- it is a verbatim duplicate of Eq. 4.17's rib validity range from Section 4.2, printed in the wrong place. |  |
| **I2** | **RESOLVED.** Reviewer checked the reference list directly: no Martin (1977) citation exists anywhere in Han's source. The removed implementation's attribution is a mistake, or refers to a correlation this edition has superseded with Goldstein et al. (1986). Adopting Goldstein et al. (1986) as the sole single-jet source. Since resolved, Martin (1977) turned up anyway, in Florschuetz's OWN reference list ([5], "Heat and Mass Transfer Between Impinging Gas Jets and Solid Surfaces," Advances in Heat Transfer 13, 1977) -- but as the precedent for the crossflow *flow-distribution* model (a slot-nozzle array, Eq. 2-6 below), not for any single-jet Nu correlation. So Martin (1977) is real, just not the source of Eq. 4.1 -- this confirms I2's call rather than reopening it. |  |
| **I3** | **RESOLVED (2026-09-21).** Scope for #337: single-jet (Eq. 4.1, Goldstein) + multi-jet array with crossflow (Eq. 4.9/Table 4.1, Florschuetz), using item 21's closed-form `Gc/Gj` and cross-checking `C_D` against combaero's existing discharge-coefficient correlation rather than always defaulting to the paper's `0.79`. Deferred, and per reviewer possibly never needed: leading-edge/curved-surface impingement (Section 4.1.4, "it is an edge case") and the simpler forms (Eq. 4.6 Kercher-Tabakoff -- graphical, not closed-form anyway; Eq. 4.7/4.8 -- Florschuetz's own less-tight alternate, item 13) as documented-but-unimplemented. | |

### What I still need sent

Nothing. The primary paper (Florschuetz, Truman and Metzger, 1981) closed
every remaining gap -- items 13, 20 and 21 below -- and confirmed every
Han-reprinted equation and constant digit-for-digit. See the review log for
what changed.

---

## Extracted items

### Single jet (Section 4.1.2.2, p. 333)

| # | item | as extracted | state |
|---|---|---|---|
| 1 | Eq. 4.1 form | `Nu_bar / Re^0.76 = (A - \|L/D - 7.75\|) / (B + C(R/D)^n)`, two branches for constant heat flux (`n=1.285`) and constant wall temperature (`n=1.394`) | **confirmed by two independent channels**: my direct read of the page image (`han_eq4.1-4.3_pp333.png`) and the reviewer's separate OCR extraction agree character-for-character, absolute-value bars and both exponents included. The plain pdftotext layer drops the bars and reads as nonsense subtraction without either channel |
| 2 | constants | `A=24, B=533, C=44` | confirmed |
| 3 | attribution | Goldstein et al. (1986), not Martin (1977) | confirmed as printed; see I2 |
| 4 | validity anchor | `Re=25,000, R/D=5, L/D=7.75` gives `Nu_bar=60` (const. heat flux) or `56` (const. wall temp.) | confirmed -- a closed-form check once implemented |
| 5 | optimum spacing | `L/D=7.75` maximises `Nu_bar` for both branches -- structurally guaranteed by the `\|L/D - 7.75\|` term, not a separate fitted fact | confirmed, and self-consistent: the "optimum at 7.75" claim in the text is exactly what the formula's own structure produces |
| 6 | recovery factor, Eq. 4.2 | `r = (T_r - T_j) / (u_j^2 / 2 c_p)` | confirmed |
| 7 | effectiveness, Eq. 4.3 | `eta = (T_aw - T_r) / (T_j - T_inf)` | confirmed |

### Multi-jet array without crossflow correction (Section 4.1.3.2-4.1.3.3, pp. 339-343)

| # | item | as extracted | state |
|---|---|---|---|
| 8 | Eq. 4.6, Kercher and Tabakoff (1970) form | `Nu_D,crossflow = phi1 * phi2 * Re_D^m * Pr^(1/3) * (Z/D)^0.091` | confirmed. `phi1`, `phi2` are read off correlation charts (Fig. 4.9 and a companion not yet located), not closed-form -- this correlation is graphical, not directly implementable as a formula |
| 9 | Eq. 4.7 | `Nu / Nu1 = 1 - C (xn/d)^nx (yn/d)^ny (z/d)^nz (Gc/Gj)^n` | confirmed |
| 10 | Eq. 4.8, normalising Nusselt number | `Nu1 = 0.363 (xn/d)^-0.554 (yn/d)^-0.422 (z/d)^0.068 Re_j^0.727 Pr^(1/3)` | confirmed |
| 11 | Eq. 4.7/4.8 constants | `C, nx, ny, nz, n` for Inline (`0.596, -0.103, -0.38, 0.803, 0.561`) and Staggered (`1.07, -0.198, -0.406, 0.788, 0.660`) | confirmed |
| 12 | stated accuracy | "confidence level for this correlation is 95%" | **confirmed as printed, but this is the wrong number for judging Eq. 4.9/Table 4.1's fit quality -- see item 13a.** Han's sentence is lifted from the primary paper's *experimental Nu measurement uncertainty* ("+-5 percent for a confidence level of 95 percent [3]"), not the correlation's own fit accuracy |
| 13 | attribution for Eq. 4.7/4.8 | **RESOLVED against the primary paper.** Eq. 4.7/4.8 is Florschuetz's own Eq. (11a)/(11b) -- the SAME 1981 paper as Eq. 4.9/Table 4.1, not a different or later one. The paper presents Eq. (10a)/(10b) (= Han's Eq. 4.9/Table 4.1, "recommended for detailed analysis... particularly in computer programs") FIRST, then Eq. (11a)/(11b) SECOND as "an alternate correlation, more convenient for hand computation... but not as tight overall" -- Han's book reprints them in the opposite order (4.7/4.8 before 4.9), which is what made "developed earlier" misleadingly read as a different, later paper. Confirmed against the page image (`page-6.png` in the local paper scan): constants match Han's item 11 digit-for-digit | confirmed |
| 13a | Eq. (11a)/(11b)'s actual accuracy, qualitative | Per the primary paper (`page-6.png`): "essentially as good as [Eq. 10a] in terms of standard error and 95 percent confidence levels, but is not as tight overall" -- no separate number is given, only this comparison to item 23's Eq. (10a)/Table 2 numbers | confirmed, qualitative only |

### Multi-jet array WITH crossflow (Section 4.1.3.3, pp. 343-344) -- the target for #337

| # | item | as extracted | state |
|---|---|---|---|
| 14 | Eq. 4.9 | `Nu = A * Re_j^m * {1 - B * [(z/d)(Gc/Gj)]^n} * Pr^(1/3)` | confirmed against the page image (`han_eq4.7-4.9_pp343.png`) |
| 15 | attribution | Florschuetz, Truman and Metzger (1981), ASME J. Heat Transfer 103, 337 | confirmed, printed directly under the equation and as Table 4.1's source line |
| 16 | A, m, B, n are each geometry-dependent | `A, m, B, n = C (xn/d)^nx (yn/d)^ny (z/d)^nz`, EACH with its own `C, nx, ny, nz` from Table 4.1 -- four independent power-law fits, not four constants | confirmed against the page image (`han_table4.1_fig4.13_pp344.png`) |
| 17 | Table 4.1 coefficients, Inline pattern | `A: C=1.18, nx=-0.944, ny=-0.642, nz=0.169`; `m: C=0.612, nx=0.059, ny=0.032, nz=-0.022`; `B: C=0.437, nx=-0.095, ny=-0.219, nz=0.275`; `n: C=0.092, nx=-0.005, ny=0.599, nz=1.04` | confirmed against the page image |
| 18 | Table 4.1 coefficients, Staggered pattern | `A: C=1.87, nx=-0.771, ny=-0.999, nz=-0.257`; `m: C=0.571, nx=0.028, ny=0.092, nz=0.039`; `B: C=1.03, nx=-0.243, ny=-0.307, nz=0.059`; `n: C=0.442, nx=0.098, ny=-0.003, nz=0.304` | confirmed against the page image |
| 19 | validity, five "geometries" printed under Table 4.1 | **NOT Florschuetz's geometry.** The printed box (`W/H`, `E`, `e/D`, `P/e`, `alpha`, `Re`, five rows: Square channel, Rectangular I/II/IA/IIA) is a **character-for-character duplicate** of Eq. 4.17's rib validity range from p. 377 (`docs/heat_transfer/han/han_ribs.md:131`: `P/e = 10-20, e/D = 0.047-0.078, alpha = 90-30 deg, W/H = 1-4, Re = 10,000-60,000`), broken out per test channel instead of quoted as one band. None of `e/D` (rib height ratio), `P/e` (rib pitch ratio) or `alpha` (rib angle) has any meaning for a jet impingement array -- Florschuetz's own test model (Fig. 4.11) is parameterised by `xn/d`, `yn/d`, `z/d` and `Re_j`, nothing else. This is a book production error, most likely the rib chapter's own table typeset under the wrong equation in the 2nd edition. `E` is whatever the rib chapter called that column, not an impingement-specific parameter -- there is nothing to resolve about its physical meaning, because it never described this correlation |
| 20 | Florschuetz's actual validity range | **RESOLVED against the primary paper's own stated overall ranges (p. 337, confirmed against `page-1.png`):** `Re_j = 2.5e3 to 7e4` (row-mean nominal range `5e3 to 5e4`, per p. 340); `Gc/Gj = 0 to 0.8`; `xn/d = 5-15` (inline) or `5-10` (staggered); `yn/d = 4-8`; `z/d = 1-3`; aspect ratio `xn/yn = 0.625 to 3.75`. Item 22's `xn/d = 5-15, yn/d = 4-8` (from Han's reprint) was the inline-pattern subset of this, not the full picture -- the primary paper adds staggered's tighter `xn/d` bound, `z/d`, `Re_j` and `Gc/Gj` bounds that Han's book never restated at all. This is the real replacement for Table 4.1's misprinted rib-validity box (item 19/I1) | confirmed |
| 21 | `Gc/Gj` row-by-row formula | **RESOLVED against the primary paper's Eqs. (1)-(8) (pp. 338-339, confirmed against `page-2.png`/`page-3.png`).** See "The Gc/Gj closed form" below for the full derivation. In short: `Gc/Gj` at spanwise row `i` is `[1/(sqrt(2) C_D)] * sinh(beta(x/xn - 1/2)) / cosh(beta(x/xn))`, `x = xn(i - 1/2)`, `beta = C_D sqrt(2) (pi/4) / [(yn/d)(z/d)]` -- a closed-form function of row number and geometry alone, no accumulation loop needed. `C_D` (jet plate discharge coefficient) is the one input this needs that Han's book never mentions at all; the primary paper recommends `0.79` absent a measured value, and Table 1 (measured, per configuration) ranges `0.73-0.85` | confirmed |
| 22 | jet-to-jet spacing range (running text, Han's book, pp. 342-343) | `xn/d = 5-15`, `yn/d = 4-8`, for the Florschuetz staggered test model (Fig. 4.11) | confirmed as printed in Han's book, and confirmed by item 20 to be the inline-pattern subset of Florschuetz's real range, not the full validity statement -- superseded by item 20 as the citable range |
| 23 | Eq. (10a)/Table 2's own fit accuracy (`page-5.png`) | Inline: 1400 points, standard error 5.6%, 95% of points within 11% of the fit, 99% within 16%, all within 19% except one outlier at 26%. Staggered: 680 points, standard error 6.1%, 95% within 12%, 99% within 16%, all within 18% except one outlier at 26%. Distinct from item 12's experimental measurement uncertainty (+-5% at 95% confidence) -- this is how well the fit itself tracks the measured data, not how accurate the measurements were | confirmed |
| 24 | jet plate discharge coefficient, `C_D` | Not carried at all into Han's reprint. Needed to evaluate item 21's `Gc/Gj` formula (`beta` depends on it). Primary paper recommends `C_D = 0.79` "for jet plates similar to those utilized here" absent a measured value (`page-6.png`, Concluding Remarks); Table 1 lists measured values per tested configuration, ranging `0.73` to `0.85` (`page-2.png`) | confirmed |
| 24a | cross-check against combaero's existing `Cd` correlations (`include/orifice.h`) | **Not reusable -- checked and rejected, not just skipped.** `Cd_sharp_thin_plate`/`Cd_ReaderHarrisGallagher` and siblings implement ISO 5167-2 flow-metering orifice plates: a single round hole IN A PIPE RUN, `beta = d/D` with `d < D` enforced, `Re_D >= 5000`, `D >= 50mm`, upstream/downstream pipe flow with defined tap locations. A jet-impingement plate is a different physical configuration entirely -- an array of small holes discharging from a plenum into an open/confined channel, no pipe, no `D`. Item 24's `C_D = 0.79` default (and Table 1's measured `0.73-0.85` range) is Florschuetz's own value for exactly this configuration and should be used directly, not derived from or reconciled with the pipe-orifice family | confirmed rejected -- use the literature value as-is |
| 24b | plate `C_D` is a real gap for combaero, not just for this correlation | Reviewer's rule of thumb in practice is `C_D = 0.8` -- consistent with Florschuetz's `0.79` default and Table 1's measured `0.73-0.85` band, but combaero has no jet-plate-array discharge-coefficient correlation of its own (item 24a: the existing family is for pipe-run metering orifices, a different configuration). Worth its own tracked scope, separate from #337 -- `han_park_1988_angled`-style parametrised set, or a fixed default with a literature citation, is a later decision, not one this extraction needs to make to implement Eq. 4.9 | noted, not yet actioned -- see #337 vs a possible new issue |
| 25 | open area ratio, `A0*` | `A0* = (pi/4) / [(xn/d)(yn/d)]` for a uniform rectangular array (both inline and staggered use the same formula -- staggered offsets rows spanwise, it does not change hole density) | confirmed against `page-3.png` |
| 26 | additional validation-figure candidates in the primary paper, not reprinted by Han at all | Fig. 5 (`Nu1` vs `yn/d`, `Re_j=1e4`, `xn/d in {5,10,15}`, inline); Fig. 6 (`Nu/Nu1` vs `Gc/Gj`, 9-panel matrix over `(xn/d,yn/d)` combinations, `z/d in {1,2,3}`, inline); Fig. 7 (staggered/inline `Nu` ratio vs `Gc/Gj`, 4-panel); Fig. 8 (`Nu1` vs `Re_j`, this correlation vs Kercher-Tabakoff vs Chance, at `(5,5)` and `(8,8)`); Fig. 9 (`Nu/Nu1` vs `(z/d)(Gc/Gj)`, same three correlations compared, at `(5,5,1)I`/`(5,5,3)I`/`(8,8,1)I`/`(8,8,3)I`) | Fig. 6 digitised and scored 2026-09-23 (see "Validation targets" below); Figs. 5, 7, 8, 9 seen, not digitised |
| 27 | Martin (1977) reference, as it actually appears | Martin, H., "Heat and Mass Transfer Between Impinging Gas Jets and Solid Surfaces," *Advances in Heat Transfer*, Vol. 13, Academic Press, New York, 1977, pp. 1-60 -- Florschuetz's reference [5], cited only as precedent for the one-dimensional crossflow *flow-distribution* model (Eqs. 2-6 below, for an array of slot nozzles), not for any single-jet heat-transfer correlation | confirmed against `page-3.png` and `page-6.png`; see I2 |

### The Gc/Gj closed form (Florschuetz Eqs. 1-8, item 21)

Confirmed against `page-2.png` (Nomenclature) and `page-3.png` (derivation and
Eqs. 1-8) directly, not the OCR text layer, which interleaves the two-column
Nomenclature block into nonsense.

Nomenclature:
- `Gc` = channel crossflow mass velocity, based on channel cross-sectional area
- `Gj` = jet mass velocity, based on jet hole area (this is the `Gj` that
  appears in `Re_j = Gj d / mu` and in the `Gc/Gj` ratio)
- `Gj*` = a distinct, *superficial* jet mass velocity, based on jet-plate (or
  opposing heat-transfer-surface) area rather than hole area -- an
  intermediate quantity in the derivation, related to `Gj` by the open area
  ratio: `Gj* = Gj * A0*`
- `A0*` = open area ratio, jet hole area / opposing heat-transfer surface
  area -- item 25
- `C_D` = jet plate discharge coefficient -- item 24
- `beta = C_D * sqrt(2) * (pi/4) / [(yn/d)(z/d)]`
- `M = sqrt(2) * A0* * C_D / z`
- `x = xn * (i - 1/2)`, `i = 1, 2, ..., Nc` -- streamwise location of spanwise
  row `i`, counting from upstream (`Nc` = number of spanwise rows, 10 in every
  tested configuration)

Derivation (a 1-D momentum/mass balance on a continuously-distributed
injection model, Fig. 4(a)/4(b) -- this is the model Martin (1977), item 27,
precedes for slot nozzles):

- Eq. (1): `Gj* = A0* C_D [2 rho (P0 - P)]^(1/2)`
- Eq. (2): `dP = -2 Gc dGc / rho`
- Eq. (3): `Gj* = z dGc/dx`
- Eq. (4): `d^2(Gc)/dx^2 - M^2 Gc = 0`, BCs `Gc=0` at `x=0`,
  `Gc = Gj_bar* * L/z` at `x=L` (`Gj_bar*` = mean `Gj*` over the array,
  `L = xn * Nc`)
- Eq. (5): `Gc / Gj_bar* = (L/z) * sinh(Mx) / sinh(ML)`
- Eq. (6): `Gj* / Gj_bar* = ML cosh(Mx) / sinh(ML)`
- Eq. (7), discrete-hole jet velocity distribution at row `i`:
  `Gj / Gj_bar = beta * Nc * cosh(beta * x/xn) / sinh(beta * Nc)`
- **Eq. (8), the crossflow-to-jet ratio at row `i` -- this is item 21's
  answer:**
  `Gc/Gj = [1 / (sqrt(2) * C_D)] * sinh(beta*(x/xn - 1/2)) / cosh(beta*(x/xn))`

Eq. (8) is evaluated one-half a hole spacing upstream of row `i` (i.e. at
`x/xn - 1/2`, while the denominator's jet velocity is evaluated at the row
itself, `x/xn`) -- this is what "resolved to one streamwise hole spacing"
means physically: `Gc` seen by a row is the crossflow accumulated from every
row upstream of it, not that row's own contribution. Both curves (Eqs. 7, 8)
were verified by the paper against directly measured pressure-traverse data
(Figs. 2, 3) and matched closely across the tested geometry range.

### Table 4.1's validity box, transcribed in full

| geometry | `W/H` | `E` | `e/D` | `P/e` | `alpha` | `Re x 1e-3` |
|---|---|---|---|---|---|---|
| Square channel | 1 | 0.24 | 0.047 | 10 | 90, 60, 45, 30 | 10, 30, 60 |
| Rectangular channel I | 2 | 0.32 | 0.047 | 10 | 90, 60, 45, 30 | 10, 30, 60 |
| Rectangular channel II | 2 | 0.32 | 0.078 | 10 | 90, 60, 45, 30 | 10, 30, 60 |
| Rectangular channel IA | 2/4 | 0.32 | 0.047 | 10 | 90, 60, 45, 30 | 10, 30, 60 |
| Rectangular channel IIA | 1/4 | 0.32 | 0.078 | 10 | 90, 60, 45 | 10, 30, 60 |

**RESOLVED: this whole box is misplaced, not Section 4.1's geometry at all.**
`alpha` here is rib angle from Section 4.2, appearing under a jet-impingement
table because the wrong table was typeset -- not a jet inclination angle, not
a coincidence, not a template reused deliberately. It is transcribed above
for the record (in case anyone else hits the same page and wonders the same
thing) but plays NO PART in implementing Eq. 4.9. See item 20 for what the
correlation's real validity range needs instead.

---

## Validation targets

For single-jet (Eq. 4.1, fully closed already): the closed-form check itself
(item 4, `Re=25000, R/D=5, L/D=7.75` -> `Nu=60`/`56`) is a sufficient
regression target; no figure digitisation is needed beyond it.

For multi-jet array with crossflow (Eq. 4.9/Table 4.1, the actual target of
#337), **2026-09-22 correction**: an earlier pass here filed Figs. 5/6
alongside 2/3/7 as "supporting results, not the Nu correlation itself" and
proposed Figs. 8/9 as the targets instead. That was backwards. Rendered at
400dpi and inspected directly:

- **Fig. 8** and **Fig. 9** plot the paper's own correlation curve alongside
  Kercher-Tabakoff's and Chance's -- THREE FITTED CORRELATIONS COMPARED, no
  raw measured points at all. Digitising the "Present Work" curve off either
  would only confirm this implementation reproduces Table 4.1's own formula,
  which is already confirmed digit-for-digit (items 17/18) and by the
  algebraic identity below -- not independent validation.
- **Fig. 6** ("Effect of crossflow and geometric parameters on streamwise
  resolved Nusselt numbers. Inline hole pattern", p. 340) is genuine measured
  scatter, and unambiguous: a 3x3 panel matrix, each panel printed with its
  exact geometry label -- `B(5,4)I, C(5,6)I, B(5,8)I / B(10,4)I, C(10,6)I,
  B(10,8)I / D(15,4)I, D(15,6)I, D(15,8)I` -- x-axis `Gc/Gj`, y-axis `Nu/Nu1`,
  three symbol classes for `z/d in {1,2,3}` (circle/square/triangle) printed
  in-panel. No leader-line tracing or attribution ambiguity: this is the
  **primary validation target**, ~15-20 points per panel per symbol.
- **Fig. 5** ("Effect of geometric parameters on Nusselt number for initial
  upstream row of array", p. 339) is also genuine measured scatter (`Nu1` vs
  `yn/d`, fixed `Re_j=1e4`, `z/d in {1,2,3}` by symbol, `xn/d in {5,10,15}` by
  point cluster) but geometry-to-cluster attribution runs through curved
  leader lines from a legend rather than an in-panel label, unlike Fig. 6 --
  usable, but needs careful tracing at digitisation time, flagged rather than
  guessed at.

**Why Fig. 6 alone is enough to validate the crossflow term, algebraically:**
`Nu = A*Re_j^m*{1-B[(z/d)(Gc/Gj)]^n}*Pr^(1/3)` and `Nu1` is the same
expression at `Gc/Gj=0` (bracket=1), so `Nu/Nu1 = 1 - B[(z/d)(Gc/Gj)]^n`
EXACTLY -- every `Re_j` and `Pr` dependence cancels in the ratio. Fig. 6's
`Nu/Nu1` vs `Gc/Gj` scatter therefore scores `jet_array_impingement_nu`'s
bracket term directly, with no `Re_j` or `Pr` assumption needed at all --
simpler to score than any rib figure, which all needed a Re bisection.

**Table 1** (measured `C_D` per configuration, item 24) remains relevant to
#375 (jet-plate discharge coefficient), not to Eq. 4.9/Table 4.1's own `Nu`
fit.

Not proposed as targets: Figs. 2/3 -- these genuinely are about the flow-
distribution model (Eqs. 5-8), supporting results rather than the Nu
correlation's own fit. Fig. 7 (staggered/inline `Nu` RATIO vs `Gc/Gj`, not
raw scatter) is discussed separately below, under staggered validation.

**Fig. 6, digitised 2026-09-23** (all nine panels, `xn_d in {5,10,15}` rows
by `yn_d in {4,6,8}` columns, `z_d in {1,2,3}` circle/square/triangle
symbols per panel): `validation/cooling/data/florschuetz1981/`, 27 series,
242 points. All nine panels share the SAME axis box, `x=[0,0.8]`,
`y=[0.4,1.0]` -- the original digitisation plan's guess that the x-range
differed per column was wrong, caught during digitisation itself (a shared
boundary label between adjacent rows read as a per-row bound at first).

Scored via `validation/cooling/jet_array_runner.py` (deliberately not
folded into `runner.py`, which is rib-specific and always bisects a
Reynolds number -- the algebraic identity above means this needs none):
pooled bias +3.7%, RMS 9.6% over all 242 points, close to
`florschuetz_1981_inline().standard_error` (5.6%, the FIT's own training
residual, not a bound on raw scatter). One panel is a genuine, understood
outlier: `B(5,4)I`'s `z/d=1` series sits at the correlation's validity
corner (`xn/d=5, yn/d=4, z/d=1`, all three simultaneously at their lower
bound) and underpredicts crossflow degradation there, error growing from
+14% to +42% with `Gc/Gj` -- axis calibration on that same panel reads
back within 0.9% of round numbers, so this is a real weak spot of the fit
at its own domain corner, not a digitisation artifact. See
`rc_00_circles_zd_1`'s `cross_check` in `metadata.yaml` for the full
figure.

Per-panel axis calibration was checked two ways rather than through
`verify.py`'s frame-slope machinery: tick positions against round numbers
(x-ticks within 0.32%, y-ticks within 0.97% across all nine panels) and
corners against the nominal box (within 3.4%). The `rc_XX_corners.csv` /
`_xaxis.csv` / `_yaxis.csv` files are committed but not listed in
`metadata.yaml` -- `verify.py`'s frame check is a log-log power-law fit
built for a physically-drawn distortion line spanning a real data range,
and throws a math-domain error on coordinates this close to the plot
origin (an exact `x=0` tick click takes `log(0)`).

**Staggered pattern (`florschuetz_1981_staggered()`) remains unvalidated
against scatter data.** Fig. 6 is Inline-only per its own caption -- there
is no staggered equivalent of it in this paper. Fig. 7 gives a
staggered/inline Nu RATIO (not raw Nu or Nu/Nu1) for a handful of
geometries, which could support an INDIRECT check later (multiply an
inline-scored Nu by the digitised ratio and compare against a staggered
measurement, if one existed) -- but Fig. 7 itself is not digitised, and no
figure in this paper gives raw staggered scatter directly. Recommendation:
leave `florschuetz_1981_staggered()` unvalidated by digitised data for now
(it is still confirmed digit-for-digit against Table 4.1, items 17/18,
the same as inline) rather than spend a digitisation pass on an indirect
ratio check; revisit if a staggered-specific source surfaces.

---

## Modelling decisions

**D1: a new, independent type family, not an extension of `RibCorrelationSet`.**
Ribs' schema extension (han_ribbed.md D5) worked because every rib set shares
one algebraic shape (a constant times power-law geometry terms). Impingement
shares nothing with that shape or between its own two regimes: single-jet
is `Re^0.76 * (A - |L/D-7.75|) / (B + C(R/D)^n)` (no geometry power-law at
all), and the jet array's four coefficients (`A,m,B,n`) are each their OWN
three-term geometry fit, correlated against a crossflow ratio ribs have no
analogue of. Forcing either into `RibCorrelationSet` would mean adding fields
that are `IGNORED` in every other shape, which is exactly what D5 was
already stretching to avoid. `impingement_correlation.h` is therefore a
sibling file with its own `SingleJetImpingementSet` and
`JetArrayCorrelationSet`, matching the same conventions (provenance-bearing
name/source fields, advisory validity via `ImpingementRange`, `validate_*`
hard-errors, `evaluate`-style functions that never throw from inside a
residual) without sharing a type.

**D2: `Gc/Gj` is a closed-form function, not a set field.** It depends only
on geometry and `C_D` (Florschuetz Eq. 8, item 21) -- it is not a fitted
coefficient like anything in `JetArrayCorrelationSet`, so it does not belong
on that struct. `crossflow_to_jet_ratio_at_x`/`_at_row` are free functions
a caller composes with `jet_array_impingement_nu` themselves, the same way a
caller would compute `Re_j` themselves before calling `evaluate_rib`.

**D3: Eq. 7 (the jet velocity distribution, mean-to-per-row `Re_j`) is
deliberately NOT implemented, in the correlation library OR the network
element.** It answers a different question -- how a total/mean array mass
flow splits across rows -- which stayed a design question even once
`ImpingementModel` was built: the element resolves it by treating
`ConvectiveSurface.area` as ONE ROW's own footprint and recovering that
row's own mass flow from it (the same "total flow through this element's
own area" convention every `ConvectiveSurface` model already uses, not a
new one), leaving row-to-row flow distribution to whoever chains multiple
row elements together. Eq. 7 would let a caller start from a single
array-mean flow instead; nobody has asked for that yet, and building it
speculatively would be the same mistake D3 originally flagged.

**D4: the crossflow bracket and Reynolds terms are floored the same way
ribs floor `e+`** (`smooth_magnitude`, `sqrt(x^2 + floor^2)`) -- a solver
probing reverse flow or a momentarily negative `Gc/Gj` gets a finite,
differentiable result instead of `NaN` under a non-integer power. This is a
numerical guard, not a physical claim: neither floor is stated by either
source, and both are small enough (`RE_FLOOR=1`, `CROSSFLOW_ARG_FLOOR=1e-3`)
to be invisible across the sets' real operating ranges.

**Implemented 2026-09-21**: `goldstein_1986_single_jet()` (single jet, Eq.
4.1) and `florschuetz_1981_inline()`/`florschuetz_1981_staggered()` (jet
array with crossflow, Eq. 4.9/Table 4.1), in `impingement_correlation.h`/
`.cpp`. Single-jet reproduces the closed-form check point (item 4) to
rounding precision. Table 4.1's coefficients match the primary paper
digit-for-digit (items 17/18/23). 21 C++ tests, 11 Python pybind-boundary
tests; one real bug caught before any test ran: the `Re` exponent (fixed
0.76) and the `(R/D)` exponent (`n`, boundary-condition-dependent) were
initially swapped, missing the check point by two orders of magnitude --
now pinned by a dedicated falsification test.

**Wired into the network solver and GUI, same day.** `ImpingementModel`
(jet array, one row) and `SingleJetImpingementModel`
(`python/combaero/network/components.py`), following D1-D3 above. Both
correlation functions gained an analytic-derivative sibling
(`single_jet_impingement`, `JetArrayImpingementResult::dNu_dRe_j`) for the
wall-coupling Jacobian, FD-verified in C++; the new mdot-to-h composition
chain (recovering a row's mass flow from `ConvectiveSurface.area`) is
independently FD-verified again at the element level, since that chain is
new code the correlation-level check does not cover. Wired into the GUI the
same way ribbed was (`ImpingementModelData`/`SingleJetImpingementModelData`
in `gui/backend/schemas.py`, mapped in `graph_builder.py`, reachable from
`SurfaceEnhancementInspector.tsx`'s dropdown), replacing the old,
unreachable, pre-0.7.0 impingement UI block that used the crossflow-less
field names. Figure-based validation-harness scoring (Figs. 8/9) still not
done -- proposed as candidates in "Validation targets" above, not yet
digitised.

---

## Review log

| date | reviewer | outcome |
|---|---|---|
| 2026-09-23 | reviewer + Claude | **Fig. 6 digitised and scored -- `florschuetz_1981_inline()` validated against its own source.** All nine panels, 27 series, 242 points, `validation/cooling/data/florschuetz1981/`; scoring via a new `jet_array_runner.py` using the `Nu/Nu1` algebraic identity (no Re bisection needed). Pooled bias +3.7%, RMS 9.6%. One genuine outlier found and kept, not hidden: `B(5,4)I`'s `z/d=1` series sits at the correlation's own validity corner and the fit underpredicts crossflow degradation there badly (+14% to +42% across the panel) -- axis calibration on that panel checked clean, so this reads as a real weak spot of the fit rather than a digitisation error. Two digitisation-time calibration bugs were caught and fixed along the way (a log10-instead-of-linear y-axis misread on one panel; a shared-boundary-label misread that made row 0/1 look like they spanned 0.6-1.0 instead of the correct, uniform 0.4-1.0 across all nine panels). Staggered validation considered and deferred: Fig. 6 is Inline-only, and Fig. 7's staggered/inline ratio is the only staggered-adjacent data in the paper -- not raw scatter, not digitised. See "Validation targets" above for the full writeup. |
| 2026-09-22 | Claude | **Corrected "Validation targets": Figs. 5/6 are the real targets, not Figs. 8/9.** Rendered Figs. 5, 6, 8, 9 at 400dpi and inspected directly, per reviewer's request to identify what needs digitising to validate Eq. 4.9/Table 4.1. Found Figs. 8/9 plot THREE FITTED CORRELATIONS against each other with no raw data at all -- digitising the "Present Work" curve would only re-confirm Table 4.1's coefficients, already confirmed digit-for-digit. Fig. 6 (p. 340) is genuine measured scatter, a 3x3 panel matrix with geometry printed in-panel (no attribution ambiguity), `Nu/Nu1` vs `Gc/Gj` -- and since `Nu/Nu1 = 1-B[(z/d)(Gc/Gj)]^n` exactly (the `Re_j^m*Pr^(1/3)` factor cancels between `Nu` and `Nu1`), it scores the crossflow term directly with no `Re_j`/`Pr` assumption at all. Fig. 5 is also real data but needs careful leader-line tracing for geometry attribution, unlike Fig. 6. Digitisation itself not yet done. |
| 2026-09-21 | Claude | **Item 24a: cross-checked `C_D` against combaero's existing discharge-coefficient correlations (`include/orifice.h`) and rejected them as inapplicable.** That family (`Cd_sharp_thin_plate`/Reader-Harris-Gallagher, Stolz, Miller, thickness/rounded-entry corrections) implements ISO 5167-2 flow-metering pipe orifices -- a single hole in a pipe run with `beta=d/D`, `Re_D>=5000`, `D>=50mm` -- not a plenum-fed multi-hole jet plate. Florschuetz's own `C_D=0.79` default (item 24) stands as-is, with Table 1's measured `0.73-0.85` available if a specific plate's value is known. Nothing left blocking I3. |
| 2026-09-21 | reviewer | **I3 resolved -- cleared for implementation.** Scope confirmed: single-jet (Goldstein) + multi-jet array with crossflow (Florschuetz), using item 21's closed-form `Gc/Gj`. `C_D` (item 24) to be cross-checked against combaero's existing discharge-coefficient correlation rather than always defaulting to the paper's `0.79` for a bare/unmeasured plate -- see item 24a for the outcome. Leading-edge/curved-surface impingement deferred as an edge case; the simpler forms (Eq. 4.6, Eq. 4.7/4.8) deferred and may never be needed. |
| 2026-09-21 | reviewer + Claude | **Items 13, 20, 21 resolved; basis review complete.** Reviewer supplied a scan of their own paper copy of Florschuetz, Truman and Metzger (1981) (`docs/heat_transfer/flohrshuetz/`, gitignored). Read page-image by page-image (it is a scanned copy, OCR text layer unreliable -- confirmed by the mangled two-column Nomenclature). Found: (a) item 13 -- Eq. 4.7/4.8 is the SAME paper's own Eq. (11a)/(11b), an alternate hand-computation form presented right after Eq. (10a)/Table 2, not a separate later paper as Han's presentation order suggested; (b) item 20 -- the paper's own stated overall ranges (`Re_j=2.5e3-7e4`, `Gc/Gj=0-0.8`, `xn/d=5-15`/`5-10`, `yn/d=4-8`, `z/d=1-3`) replace the misprinted Table 4.1 box entirely; (c) item 21 -- the paper's Eqs. (1)-(8) give a full closed-form `Gc/Gj` at any spanwise row, needing only `C_D` (recommended default `0.79`) and geometry, no accumulation loop; (d) Martin (1977) turned up as Florschuetz's own reference [5] -- real, but cited only for the crossflow flow-distribution model precedent, not for any single-jet correlation, which supports rather than reopens I2. Also found richer accuracy data (item 23) and new validation-figure candidates (item 26, Figs. 8/9) not present in Han's reprint at all. Nothing from Han's reprinted equations or Table 4.1 constants disagreed with the primary source -- every digit checked matched. |
| 2026-09-21 | reviewer | **Item 20's geometry resolved.** Confirmed that the running-text passage on p. 342-343 (`xn/d = 5-15`, `yn/d = 4-8`) is Florschuetz's own restated geometry, i.e. the geometric range the correlation covers -- not just a description of the Fig. 4.11 test rig in the abstract. Reviewer's explicit caveat: this is the geometric range, not full validity -- `Re_j` and `z/d` bounds are still not stated anywhere in the pages read and remain open. |
| 2026-09-21 | Claude | Searched Sections 4.1.3, 4.1.3.1, 4.1.3.2 in full for a `Gc`/`Gj` formula or a restated Florschuetz validity box; found neither. Item 21 (`Gc/Gj` definition, including row-by-row crossflow accumulation for a multi-row array) is genuinely absent from Han's text as extracted so far, not a re-reading gap -- flagged as needing the primary Florschuetz (1981) paper, matching the `rallabandi_2009_high_re` precedent (D4). |
| 2026-09-21 | reviewer | **I2 resolved.** Checked the reference list directly: no Martin (1977) citation exists in this source. The removed implementation's attribution is a mistake or has been superseded -- treat Goldstein et al. (1986) as the sole single-jet source and do not pursue Martin further. Also: before any implementation, the basis must be fully reviewed and validation targets defined, matching the discipline #333 already established for ribs. |
| 2026-09-20 | reviewer + Claude | **I1 resolved: Table 4.1's geometry box is a misprint, not Florschuetz's data.** Its `W/H, E, e/D, P/e, alpha, Re` values are character-for-character identical to Eq. 4.17's rib validity range (p. 377, `han_ribs.md:131`), which makes no physical sense under a jet-array correlation -- there are no ribs, no rib pitch, no rib angle in an impingement chamber. Almost certainly the rib chapter's own table typeset in the wrong place during the 2nd edition's production. `E`'s meaning needs no further guessing: it never described this correlation. Florschuetz's actual validity range (item 20) is still missing and needs the primary paper or an unchecked earlier page. |
| 2026-09-20 | reviewer | supplied an independent OCR extraction of Eq. 4.1 and its surrounding text. Matches my page-image read exactly, including the absolute-value bars and both exponents (1.285/1.394) -- item 1 promoted to confirmed by two channels. I1, I2, I3 still open. |
| 2026-09-20 | extracted by Claude | UNCONFIRMED. Every equation cross-checked against the rendered page image (not the plain-text layer, which drops the absolute-value bars in Eq. 4.1). Three items need a reviewer: the `E` column's meaning (I1), the Goldstein-vs-Martin citation question (I2), and the proposed implementation scope (I3). |
