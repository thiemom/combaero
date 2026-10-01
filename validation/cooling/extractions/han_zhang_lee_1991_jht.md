# Han, Zhang and Lee (1991) JHT 113, 590: the primary behind figure 4.51

**Status: extracted 2026-10-01, single channel (visual read of the
publisher PDF, text layer cross-read). Not committed as scored series.**

Source: `docs/heat_transfer/film/han_1991_590_vol_113.pdf` (gitignored --
copyrighted; filed under `film/` by mistake, it is a rib paper).
Han, J.C., Zhang, Y.M. and Lee, C.P. (1991), "Augmented Heat Transfer in
Square Channels With Parallel, Crossed, and V-Shaped Angled Ribs", ASME J.
Heat Transfer 113(3), 590-596. Supported by General Electric.

**Not the conference twin of 91-GT-3.** 91-GT-3
(`han_zhang_lee_91gt3_table3.md`) is the same authors' companion study on
ribbed-to-smooth heat FLUX RATIO. Their `R` agrees; their `G` does not --
see "Cross-check against 91-GT-3" below.

## Rig and conventions (stated, pp. 590-592)

| item | value | where |
|---|---|---|
| channel | square, 5.08 x 5.08 cm, L = 101.6 cm, L/D 20, two opposite ribbed walls | p. 591 |
| ribs | square brass, `e/D` 0.0625, `P/e` 10, **in-line** (directly opposite) | p. 590 |
| Re | 15,000-90,000 | abstract |
| `f` | channel with two opposite ribbed walls, isothermal, Eq. (1) | p. 592 |
| `f_0` | `0.046 Re^-0.2`, "Blasius", smooth circular tube, Eq. (2) | p. 592 |
| `f_r` | four-sided ribbed equivalent, `f_r = 2f - f_0`, Eq. (6) | p. 592 |
| `Nu_0` | `0.023 Re^0.8 Pr^0.4`, McAdams/Dittus-Boelter, Eq. (4) | p. 592 |
| area basis | **smooth-wall (projected)** area for both walls | p. 592 |
| averaging | ten isolated copper sections; channel-averaged Nu is over the whole heated length, **X/D 0-20, entrance region included** | pp. 591, 594 |
| `St_bar` | `(St_r + St_s)/2`, Eq. (10) | p. 592 |
| uncertainty | `f` < 8%, `Nu` < 8% (Kline-McClintock, Re > 10,000) | p. 592 |
| smooth check | smooth duct within 10% of Dittus-Boelter | p. 594 |

`R`, `G`, `G_bar` are built on `f_r` exactly as Han (1988), Eqs. (5)-(9),
so they are directly comparable with `han_1988_orthogonal`.

## Table 2 (p. 595)

`R = a(e+)^b`, `G = a(e+)^b`, `G_bar = a(e+)^b`. `b = 0` for `R` in every
case -- `R` is independent of `e+`, as Han states elsewhere.

| case | figure label | `R` a | `G` a | `G` b | `G_bar` a | `G_bar` b | `G_bar/G` at e+ 100 / 300 / 1000 |
|---|---|---|---|---|---|---|---|
| 1 | 90 | 3.18 | 3.97 | 0.28 | 4.86 | 0.28 | 1.224 / 1.224 / 1.224 |
| 2 | 60 parallel | 2.05 | 1.52 | 0.41 | 2.41 | 0.36 | 1.259 / 1.192 / 1.122 |
| 3 | 60 crossed | 3.18 | 3.24 | 0.32 | 4.57 | 0.28 | 1.173 / 1.123 / 1.070 |
| 4 | 60 V | 1.72 | 1.35 | 0.42 | 1.76 | 0.40 | 1.189 / 1.163 / 1.135 |
| 5 | 60 Lambda | 1.42 | 1.59 | 0.43 | 2.12 | 0.41 | 1.216 / 1.190 / 1.161 |
| 6 | 45 parallel | 3.05 | 2.07 | 0.36 | 3.01 | 0.32 | 1.209 / 1.157 / 1.103 |
| 7 | 45 crossed | 4.40 | 2.11 | 0.37 | 2.43 | 0.37 | 1.152 / 1.152 / 1.152 |
| 8 | 45 V | 2.04 | 1.36 | 0.43 | 1.93 | 0.40 | 1.236 / 1.196 / 1.154 |
| 9 | 45 Lambda | 1.70 | 1.83 | 0.41 | 2.49 | 0.38 | 1.185 / 1.147 / 1.106 |

**`G_bar/G` is above 1 for all nine, 1.07-1.26** -- inside the 1.096-1.413
the #401 pool found, and Case 1's 1.224 matches 91-GT-3 Case 1's 1.225.
Nothing in the primary supports a ratio below 1.

The 90 deg row, `G = 3.97 (e+)^0.28`, sits within 7% of Han (1988)'s
`3.7 (e+)^0.28` and the paper says so ("about the same as the previous
correlation"); `R` 3.18 against 3.2.

## #403 item 1: which figure 4.51 panel duplicated the other

Every digitised figure 4.51 series scored against both printed fits:

| class | G panel vs printed G | G_bar panel vs printed G | G_bar panel vs printed G_bar |
|---|---|---|---|
| six confirmed classes (60 par/crs/V/Lambda, 45 V/Lambda) | -0.3% to +1.5% bias, MAE <= 2.8% | -12.4% to -14.5% | -1.7% to +1.8%, MAE <= 2.5% |
| 45 parallel (was disputed) | +2.8% | **+3.8%** | +18.6% |
| 45 crossed (was disputed) | -2.9% | **-2.1%** | +12.8% |

The confirmed classes show the split a correct pair must: each panel on its
own correlation, 12-15% apart. **For the two disputed classes, BOTH panels
land on the printed `G`.** So the `G` series are right and the `G_bar`
series carry the `G` marks a second time. Reading (1) of the old
`cross_check` -- G duplicated into the G_bar panel -- is the one the primary
supports. Whether the duplication happened in the textbook redraw or in our
digitisation is not established, and the fix does not depend on it.

Applied in `data/han2012/metadata.yaml`: `fig4.51_G_45par` and
`fig4.51_G_45crs` -> `class_confidence: confirmed`;
`fig4.51_Gbar_45par` and `fig4.51_Gbar_45crs` -> `scores: null`, label kept
`disputed` so nothing pools them as `G_bar`.

## Cross-check against 91-GT-3, same lab and authors

Same `R` (3.18 at 90 deg and 60 crossed, 2.05 at 60 parallel, 1.72 at 60 V).
Different `G`: 91-GT-3 Case 1 over this paper's Table 2 is 0.77 / 0.90 /
1.07 at e+ 100 / 300 / 1000 for 90 deg, and 0.85-1.11 for the three angled
rows -- steeper exponents in 91-GT-3 (0.41-0.48 against 0.28-0.42). Two
papers from one rig, nominally the same uniform-heating case, disagree by
up to 23% in `G` at the low end. Recorded, not reconciled: it bounds how
far any `G` fit in this family can be read as a property of the rib rather
than of the data reduction.

## What this paper makes possible (not done here)

- **#434, rib shape.** Table 2 is a printed `R`/`G`/`G_bar` source for
  crossed and Lambda ribs as well as V and parallel -- the half #434 called
  unavailable. Fidelity data is figure 4.51 itself (this paper's Fig. 12).
  Validity is a single point in `e/D`, `P/e`, `W/H`, in-line ribs only.
- **Figure 4.53's continuous ribs** (#435/#440) appear to be this paper's
  data: its text gives 60 deg V at Nu ratio 2.7-3.5 for f ratio 8-11, which
  the digitised series reproduce (2.7-3.6 at 8.3-11.1). Its Figs. 8/9 plot
  both ratios against Re, which would remove the Re recovery step.
- **Figure 4.51 geometry.** Its series carry `W/H` 1 or nothing; the rig is
  now stated (`e/D` 0.0625, `P/e` 10). Adding it changes what
  `han_park_1988_angled` predicts (its square-channel `G` carries
  `(p/e/10)^0.1`, currently evaluated at the probe's `P/e` 15), so it is a
  separate, measured change.
