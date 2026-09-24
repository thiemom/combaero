# Extraction: Lau, Kukreja and McMillin (1990), V-shaped rib arrays

**Status: CONFIRMED for the 90 deg baseline (case 1a/-90), which is the only
configuration an implemented set can be compared against. The V-shaped cases
are transcribed and NOT cleared for implementation -- see "For the reviewer".**

Tracked by #385, under #333.

> Lau, S.C., Kukreja, R.T. and McMillin, R.D. (1990). "Effects of V-shaped rib
> arrays on turbulent heat transfer and friction of fully developed flow in a
> square channel." *Int. J. Heat Mass Transfer* 34(7), 1605-1616.
> `docs/heat_transfer/ribs_lau_1990.pdf` (gitignored, copyrighted).

**The first rib source outside Han**, and a primary paper rather than a
textbook reprint. #333 named this gap explicitly: agreement with the curve a
correlation was fitted to is a consistency check, not a validation.

It also arrives **tabulated**. Table 2 prints the roughness functions as
closed forms, so nothing about the correlation needed digitising -- which is
the preference `extractions/README.md` states.

---

## For the reviewer

| # | what needs judgement | what turns on it |
|---|---|---|
| **L1** | **Is a Lau V-rib `RibCorrelationSet` wanted?** Table 2 gives `R`/`G`/`Gbar` coefficients for V-45, V-60, V-120 and V-135, a shape `han_park_1988_angled` cannot express (it carries angle only, which is why Han figure 4.51's own V classes are deliberately unscored). | A code task shaped like #334, not a data task. Deferred here; this document transcribes the coefficients so the decision has something to rest on. |
| **L2** | **Is the ~13% `R` gap acceptable as a recorded disagreement?** Definitional causes were checked and ruled out (below). Two labs differ by ~12% in ribbed-wall friction for nominally the same 90 deg configuration. | Nothing breaks either way, but it bounds how much trust the friction side of `han_1988_orthogonal` should carry. A third source would settle whether Han or Lau is the outlier. |

### What I still need sent

Nothing. The scan's right margin is clipped on page 2, so the rib-geometry
sentence was read in fragments -- but the numbers are confirmed by an exact
arithmetic identity (below) and by the reviewer independently.

---

## Geometry and conditions

| quantity | value | |
|---|---|---|
| channel | square, aluminium, 1.52 m long, 7.62 x 7.62 cm | confirmed; reviewer read the same sentence independently |
| `D` | 0.0762 m | square, so `D_h` = side |
| `L/D` | 19.9 | |
| rib height `e` | 4.76 mm | read from a clipped column, confirmed by the identity below |
| **`e/D`** | **0.0625** | |
| ribbed walls | two opposite; two smooth | |
| `p/e` | 10 and 20 | |
| `alpha` | 45, 60, 90, 120, 135 deg | V-shaped; case 1a is 90 deg full ribs |
| `Re_D` | 10 000 - 60 000 | |
| measurement window | `x/D` = 7.66 - 15.16, thermally fully developed | |
| properties | at average bulk temperature | |

**A check that could have failed.** `e/D` was read off a page whose right
margin the scan cuts. It is exact by design, not by rounding: 7.62 cm is
exactly 3 inches and 4.76 mm exactly 3/16 inch, so
`e/D = (3/16)/3 = 1/16 = 0.0625`. An imperial design reproducing an exact
binary fraction is not something a misread digit produces.

**Stated experimental uncertainty** (Kline and McClintock, ref [6]):
`Re` +/-2.9%, Stanton numbers +/-5.8%, friction factor +/-10.9%. This is a
property of the measurement, so it is what populates `uncertainty` in the
committed metadata -- the distinction #389 exists to enforce.

---

## Extracted items

| # | item | content | status |
|---|---|---|---|
| 1 | Eq. (5), roughness Reynolds number | `e+ = (e/D) Re_D [(2 fbar - f_ss)/2]^0.5` | confirmed |
| 2 | Eq. (6), friction roughness function | `R = [(2 fbar - f_ss)/2]^-0.5 + 2.5 ln[2(e/D)] + 2.5` | confirmed |
| 3 | Eq. (7), heat transfer roughness function | `G = [(2 fbar - f_ss)/2]^0.5 / St_r + 2.5 ln[2(e/D)] + 2.5` | confirmed |
| 4 | Eq. (8), average heat transfer roughness function | `Gbar = [(2 fbar - f_ss)/2]^0.5 / St_avg + 2.5 ln[2(e/D)] + 2.5` | confirmed -- **and not Han's `G_bar`, see item 8** |
| 5 | `(2 fbar - f_ss)/2` | the ribbed-wall friction extracted from the channel average, because only two of four walls are ribbed. `fbar` alone is the channel average | confirmed |
| 6 | `St_r` basis | ribbed-wall heat flux divided by the **projected** area, "not including the increased rib surface area" | confirmed |
| 7 | Table 2 | seven cases, `R`/`G`/`Gbar` each as `a(e+)^b` -- reproduced below | confirmed, single channel; see the review log |
| 8 | **`Gbar` is a four-wall average** | uses `St_avg` over two ribbed and two smooth walls. Han's `G_bar` is the **Prandtl-normalised** `G Pr^-0.57`, which `han_ribbed.md` item 10 resolved and explicitly recorded as "not a four-wall average" | confirmed -- **a false-confirmation hazard, see below** |
| 9 | geometry group | Lau's `2(e/D)`; combaero's `(2 e/D)(2W/(W+H))`. For a square channel `2W/(W+H) = 1`, so they coincide | confirmed by derivation |

### Table 2, as printed

| case | `R` a | `R` b | `G` a | `G` b | `Gbar` a | `Gbar` b |
|---|---|---|---|---|---|---|
| 1a/-90 | 3.674 | -0.0033 | 4.218 | 0.257 | 5.450 | 0.250 |
| 4/-V-45 | 1.605 | 0.0015 | 1.819 | 0.355 | 2.656 | 0.335 |
| 5/-V-60 | 1.232 | 0.0454 | 1.299 | 0.399 | 2.163 | 0.360 |
| 5a/-V-60 | 3.687 | -0.0571 | 1.685 | 0.376 | 2.719 | 0.336 |
| 6/-V-120 | 1.537 | -0.0228 | 1.983 | 0.352 | 3.460 | 0.315 |
| 6a/-V-120 | 3.617 | -0.0734 | 2.739 | 0.324 | 4.351 | 0.291 |
| 7/-V-135 | 2.070 | -0.0384 | 1.992 | 0.362 | 3.685 | 0.313 |

---

## The false-confirmation hazard: `Gbar`

The most important thing in this document, because it fails *quietly and in
the flattering direction*.

Lau's `Gbar` and Han's `G_bar` are different quantities sharing a symbol, and
their ratios to `G` land 2-3% apart:

| | Lau `Gbar/G` | Han `G_bar/G` |
|---|---|---|
| `e+` = 50 | 1.2572 | 1.2162 |
| `e+` = 100 | 1.2511 | 1.2162 |
| `e+` = 400 | 1.2390 | 1.2162 |

Scoring one against the other returns a few percent and reads as an
independent confirmation. Nothing in the schema objects: same symbol, same
units, overlapping range, plausible number.

Worse, the mistake is **automated**: `runner.py` multiplies any series tagged
`y_axis: G_bar` by `G_BAR_OVER_G`, Han's Prandtl factor. Mislabelling Lau's
`Gbar` would silently apply a conversion that does not belong to it.

**So Lau's `Gbar` is deliberately not committed.**

A structural tell, free here and worth generalising (see #389): Lau's ratio
**drifts with `e+`** (exponents 0.250 vs 0.257) while Han's is **constant**
(both 0.28). A Prandtl normalisation is a constant factor; a wall-averaging
ratio is not.

---

## Results

### `G` -- an independent confirmation

`han_1988_orthogonal` prints `G = 3.7 (e+)^0.28`; Lau prints
`4.218 (e+)^0.257`. Different coefficient *and* exponent, yet the two forms
cross near `e+ = 300` and agree across the whole overlapping range:

| `Re` | `e+` | Lau `G` | combaero `G` | |
|---|---|---|---|---|
| 10 000 | 106 | 13.98 | 13.65 | +2.4% |
| 30 000 | 318 | 18.54 | 18.57 | -0.1% |
| 60 000 | 636 | 22.16 | 22.55 | -1.7% |

Scored through the harness: **N=9, MAE 1.2%, RMSE 1.4%, bias -0.1%, 100%
within** Lau's stated +/-5.8% Stanton uncertainty.

Robust to the geometry: `han_1988_orthogonal` carries no `e/D` exponent for
`G`, and the comparison is made at matched `e+`, varying under 0.9% across
`e/D = 0.047-0.078`.

### `R` -- a ~13% disagreement that is not definitional

| `Re` | `e+` | Lau `R` | combaero `R` | |
|---|---|---|---|---|
| 10 000 | 106 | 3.618 | 3.200 | +13.1% |
| 30 000 | 318 | 3.605 | 3.200 | +12.7% |
| 60 000 | 636 | 3.597 | 3.200 | +12.4% |

Ruled out: `han_1988_orthogonal` has no `e/D` exponent for `R`, and at
`p/e = 10` its `p/e` term sits exactly at reference, so the guessed geometry
cannot be the cause. The geometry group matches for a square channel (item 9),
and both `R` are built on the ribbed-wall friction (item 5).

Converted via `R = sqrt(2/f_r) + 2.5 ln(2 e/D) + 2.5`: Han's `R = 3.2` implies
`f_r = 0.0575`, Lau's `~3.6` implies `0.0504`. **Two labs about 12% apart in
ribbed-wall friction factor for nominally the same configuration.** Recorded,
not reconciled -- see reviewer item L2.

Not scored by the harness: `runner.py`'s `e+` path resolves `y_axis` to `G`,
`G_bar` or `R_normalised` and has no absolute-`R` branch. The series is
committed with `scores: null` so the gap is visible rather than absent.

### Reflection folding is falsified

`han_ribbed.md` recorded reflection folding (`alpha -> 180 - alpha`) as the
bounded alternative if Eq. 4.17 were ever guarded past 90 deg. Folding
predicts V-120 behaves as V-60 and V-135 as V-45. Lau measured both:

| `e+` | `G(V-60)` | `G(V-120)` | | `G(V-45)` | `G(V-135)` | |
|---|---|---|---|---|---|---|
| 100 | 8.16 | 10.03 | **+22.9%** | 9.33 | 10.55 | **+13.1%** |
| 400 | 14.19 | 16.34 | +15.2% | 15.26 | 17.43 | +14.2% |

Higher `G` means lower heat transfer, so the reversed arrays transfer less --
and the paper's own conclusion 4 says so independently: "Reversing the 45 and
60 V-shaped rib arrays lowers the heat transfer from the channel walls and
increases the channel pressure drop."

**The open item in `han_ribbed.md` is answered: folding is not physical.** No
guard should be written on that basis.

---

## Review log

| date | reviewer | outcome |
|---|---|---|
| 2026-09-24 | reviewer | Channel dimensions confirmed independently: "a straight, square channel made of aluminium, 1.52 m long, cross section 7.62 by 7.62 cm". Also raised the `Gbar`/`fbar` question, which is what surfaced the false-confirmation hazard above -- it had not been checked, and the `G` comparison already made was one symbol away from being invalid. |
| 2026-09-24 | extracted by Claude | Table 2 read off a 300 dpi render, **one channel only** -- it is a clean typeset table with no ambiguous glyphs, but a second read is still outstanding and is the main residual risk in this document. Definitions (Eqs. 5-8) read at 320 dpi. Geometry read from a column whose right margin the scan clips, and confirmed by the 3 inch / (3/16) inch identity. `R`/`G` findings and the reflection-folding falsification computed from the printed coefficients, not from any digitised point. |

Change **Status** at the top when reviewed, and record corrections here rather
than silently editing the tables above.
