# Extraction: Han rib correlations

**Status: UNCONFIRMED -- awaiting review against the book.**

Nothing here may be implemented until this document is reviewed and its status
changed. See the review log at the bottom. Tracked by #334, under #339.

---

## Source

> Han, J., Dutta, S. and Ekkad, S. (2012). *Gas Turbine Heat Transfer and
> Cooling Technology*. 2nd Edition. CRC Press.

Secondary source. The rib correlations are **after Han, J.C., ASME J. Heat
Transfer, 110, 321, 1988**, per the caption of Figure 4.46. Figure 4.46 also
carries data attributed to Han (1984).

Page images are kept in `docs/sources/`, which is gitignored -- they are
copyrighted and must not be committed. This document is the tracked artifact.

| item | file in `docs/sources/` |
|---|---|
| Eq. 4.15, 4.16, e+, four-sided f | `han_2012_rib_eq_4.15_4.16.png` |
| Fig. 4.46, Eq. 4.17 | `han_2012_rib_fig_4.46_eq_4.17.png` |

---

## Extracted items

Each line below is what a reviewer should compare against the book, one at a
time. State is one of **confirmed**, **suspect**, or **missing**.

### Definitions

| # | item | as extracted | state |
|---|---|---|---|
| 1 | dimensionless average velocity | `u+ = (2/f)^(1/2)` | confirmed |
| 2 | dimensionless average temperature | `T+ = (f/2)^(1/2) / St` | confirmed |
| 3 | roughness Reynolds number | `e+ = (e/D) * Re * (f/2)^(1/2)` | confirmed |
| 4 | Eq. 4.15, roughness function | `R(e+) = (2/f)^(1/2) + 2.5 ln( (2e/D) * (2W/(W+H)) ) + 2.5` | confirmed |
| 5 | Eq. 4.16, heat-transfer roughness function | `G(e+,Pr) = R(e+) + [ (f/(2*St_r)) - 1 ] / (f/2)^(1/2)` | confirmed |
| 6 | four-sided ribbed channel friction factor | `f = f_bar + (H/W) * (f_bar - f_s) * f_bar` | **suspect** |

**Item 6 is the one to check most carefully.** Two independent readings -- a
visual read and macOS Vision OCR -- both report a trailing multiplication by
`f_bar` after `(f_bar - f_s)`, and a high-resolution crop shows a centred dot
followed by an overbarred `f`. So the *transcription* is settled. What is not
settled is whether it is correct **in the book**: as printed the sidewall
correction is negligible, which is odd for something given its own equation.

| H/W | as printed | without the trailing `* f_bar` | ratio |
|---|---|---|---|
| 0.25 | 0.020075 | 0.023750 | 1.18x |
| 1.0 | 0.020300 | 0.035000 | 1.72x |
| 4.0 | 0.021200 | 0.080000 | 3.77x |

(evaluated at `f_bar = 0.020`, `f_s = 0.005`)

Resolving this needs a worked example in the book, or Han (1988) itself. It
must **not** be resolved by judging which looks more physical -- that is the
move that produced the defects this rebuild exists to undo.

### Correlations, 90 degree orthogonal ribs (Figure 4.46)

| # | item | as extracted | state |
|---|---|---|---|
| 7 | roughness function | `R / (P/e/10)^0.35 = 3.2`, independent of e+ | confirmed |
| 8 | heat-transfer roughness function, solid line | `G = 3.7 * (e+)^0.28` | confirmed |
| 9 | heat-transfer roughness function, dashed line | `G_bar = 4.5 * (e+)^0.28` | confirmed |
| 10 | which of items 8 and 9 is the ribbed wall vs the channel average | unbarred = ribbed sidewall, barred = channel average | **assumed, needs the surrounding text** |
| 11 | figure validity, this study | `alpha = 90 deg`, `10,000 <= Re <= 60,000` | confirmed |
| 12 | figure validity, Han (1984) data | `alpha = 90 deg`, `8,000 <= Re <= 80,000` | confirmed |
| 13 | plotted e+ range | roughly 40 to 1000 | confirmed |

Item 7 is the important structural fact: the lower panel is a **horizontal
line**, so `R` does not depend on `e+`. Consequences below.

Geometries carried in the figure legend, usable as harness cases:

| set | e/D | P/e | W/H |
|---|---|---|---|
| this study | 0.047 | 10 | 1 |
| this study | 0.047 | 20 | 1 |
| this study | 0.047 | 10 | 2 |
| this study | 0.047 | 20 | 2 |
| this study | 0.078 | 10 | 4 |
| this study | 0.078 | 20 | 4 |
| Han (1984) | 0.063 | 10 | 1 |
| Han (1984) | 0.063 | 20 | 1 |
| Han (1984) | 0.042 | 10 | 1 |
| Han (1984) | 0.021 | 10 | 1 |

### Correlation with rib angle (Eq. 4.17)

| # | item | as extracted | state |
|---|---|---|---|
| 14 | Eq. 4.17 | `R / [ (P/e/10)^0.35 * (W/H)^m ] = 12.3 - 27.07*(alpha/90) + 17.86*(alpha/90)^2` | confirmed |
| 15 | the exponent `m` | not given on the extracted page | **missing** |
| 16 | Figure 4.47, angled ribs | referenced in the text, not extracted | **missing** |

---

## Cross-checks performed

These are checks that could have failed and did not. They test the extraction,
not the physics.

**1. Eq. 4.15 and Eq. 4.16 combine into the analogy structure.**
The numerator of 4.16 is `T+ - u+` in disguise, and `u+` is the leading term of
`R`, so `u+` cancels exactly:

```
G = T+ + 2.5 ln( (2e/D) * (2W/(W+H)) ) + 2.5
```

which is `R` with temperature substituted for velocity -- what the
heat-momentum analogy requires. The expectation is external to the page, so
this is not circular. Verified numerically at four operating points, agreeing
to 3.55e-15.

Falsified: five plausible misreadings of 4.16 move `G` by 3.65 to 25.00 at a
reference point where the correct reading gives exactly 0.

**2. Eq. 4.17 reduces to Figure 4.46 at alpha = 90 degrees.**

```
Eq. 4.17 at alpha=90 : 12.3 - 27.07 + 17.86 = 3.0900
Fig. 4.46            :                        3.2
agreement            : 3.44%
```

`(W/H)^m` drops out at `W/H = 1`, so this holds without knowing `m` (item 15).
Two separately extracted equations from different parts of the book agreeing at
their common point.

**3. The prediction chain lands the figure's own legend geometries inside the
figure's own e+ range.** Eleven of twelve fall within 40 to 1000; the twelfth
(e/D=0.078, W/H=4, Re=60,000) reaches 1122, 12% beyond the axis, at the extreme
corner of the plotted set.

**4. Two independent extraction channels.** Visual reading and macOS Vision OCR
(`scripts/ocr_page.swift`). Both report `12.3`, `27.07`, `17.86`, `4.5`, `3.7`,
`0.28` and `3.2` identically. OCR garbles the typeset fractions -- it read
`(3)"/2` for `(2/f)^(1/2)` -- so it is a disagreement detector, not a source of
truth. Its per-line confidence self-flags: the one badly garbled region
reported 0.30 where everything else reported 1.00.

---

## Derived prediction chain

Follows from items 4, 5, 7 and 8. **Closed form -- no iteration.**

```
R    = 3.2 * (P/e/10)^0.35                                  item 7
f    = 2 / [ R - 2.5 ln( (2e/D)*(2W/(W+H)) ) - 2.5 ]^2      invert item 4
e+   = (e/D) * Re * (f/2)^(1/2)                             item 3
G    = 3.7 * (e+)^0.28                                      item 8
St_r = f / ( 2 * [ 1 + (G - R) * (f/2)^(1/2) ] )            invert item 5
```

This corrects an earlier claim on #334 that the correlation is implicit and
would need an inner iteration carrying its own analytic derivative. It does
not: `R` is independent of `e+` (item 7), so Eq. 4.15 inverts directly for `f`.
The error came from reasoning about the structure instead of reading the
figure.

Worked values from the chain, for review against the book if it carries an
example:

| e/D | P/e | W/H | Re | f | e+ | G | St_r |
|---|---|---|---|---|---|---|---|
| 0.047 | 10 | 1 | 10,000 | 0.04576 | 71.1 | 12.21 | 0.00968 |
| 0.047 | 10 | 1 | 60,000 | 0.04576 | 426.6 | 20.16 | 0.00642 |
| 0.078 | 10 | 4 | 10,000 | 0.11503 | 187.1 | 16.01 | 0.01413 |

Note `f` is independent of `Re` in this chain, which follows from item 7 and is
worth a reviewer's attention.

---

## Validation target

Figure 4.46 carries **no data table**, so the scatter would have to be
digitised. But the correlation lines are printed on the plot as equations
(items 7, 8, 9), which splits the job in two:

**The correlation line is the target for correctness, and it needs no
digitising.** An implementation of Han's correlation must reproduce
`R/(P/e/10)^0.35 = 3.2` and `G = 3.7 (e+)^0.28` exactly, because those *are*
the correlation. Reproducing them is a closed-form check.

**The scatter is the target for tolerance.** How far the experimental points
fall from the printed line is the accuracy of Han's correlation against his own
measurements -- which is what a harness tolerance should be set from, rather
than a number chosen for convenience. Digitising the scatter is therefore worth
doing, but for the band, not for the values.

Keeping these apart matters. Our code is not obliged to reproduce Han's
measurements; it is obliged to reproduce Han's correlation, and to be honest
about how well that correlation matched reality. Conflating the two is how a
tolerance ends up wide enough to admit a 4-5x error, which is what
`assert 1.0 < multiplier < 10.0` did in the removed implementation.

Not yet done: digitising the scatter of either panel. Until then the harness
has a correctness target but no defensible tolerance.

---

## Bearing on the removed implementation

Not part of the extraction; recorded because it explains what was removed
in #332 and should inform the re-add.

Han gives `R = 3.2 * (P/e/10)^0.35`. The deleted `rib_enhancement_factor` was:

```
3.5 * (e/D)^0.35 * (P/e/10)^-0.15 * f_alpha
```

The constant is close (3.5 against 3.2), the `(P/e/10)` group is exactly Han's,
and the exponent `0.35` is present -- but attached to `(e/D)` rather than to
`(P/e/10)`, and the `-0.15` has no counterpart. It reads as a garbled
recollection of `R`.

More seriously: **`R` is a roughness function, inversely related to the
friction factor.** It is not a Nusselt ratio. The removed code returned it as
`Nu_rib / Nu_smooth`, which is a category error rather than a constant error,
and is consistent with the friction side ending up 4-5x low.

---

## Open items

- **6** four-sided friction factor: transcription settled, correctness suspect
- **10** which of `G` and `G_bar` is the ribbed wall
- **15** the exponent `m` in `(W/H)^m`
- **16** Figure 4.47, angled ribs
- whether the book carries a worked example, which would resolve item 6 and
  give the chain an end-to-end check
- digitise the Figure 4.46 scatter, to set the harness tolerance from the
  source's own spread rather than from a chosen number

---

## Review log

| date | reviewer | outcome |
|---|---|---|
| 2026-09-10 | extracted by Claude | UNCONFIRMED -- submitted for review |

Change **Status** at the top of this file when reviewed, and record corrections
here rather than silently editing the tables above -- a correction is evidence
about the extraction process, not just about the number.
