# Fitting a superposition correction on Murray

Attempted 2026-09-30, after #425 established that `murray2018` is the first
family whose curves all fall the same side of the model.

**Outcome.** A one-knob correction works well and nearly reaches the
measurement floor. Gao's alpha is the *wrong knob for this dataset* -- not
because the form is wrong, but because it is parameterised along an axis
Murray does not vary. The older "third category" coupling form fits better
here and is the right family for data like this.

## TWO CORRECTIONS to earlier readings, both found by re-reading the paper

**1. Gao DOES publish the coefficients.** Section 4.3.1:

> "The empirical coefficients in Equation (5), a and b, were determined to
> be 12 and 0.9465, respectively."

This repo previously recorded that they were never printed and that the
search had been exhaustive. Wrong on both counts -- the search looked for
`a = 12` and the paper states the value in prose. A regex over a paper's
prose does not justify the word "exhaustive". Now recorded in code as
`film_superposition::gao_a_case1` / `gao_b_case1`.

**2. Gao calibrates on ONE geometry across blowing ratios**, which is the
same axis `murray2018` varies:

> "The empirical coefficients in the correction model were calibrated using
> the cooling efficiency distributions from Case 1 at blowing ratios of 0.3
> and 1.0."

So the approach taken here -- fit `(a, b)` on one plate across blowing
ratios -- matches Gao's own procedure. An intermediate version of this
record claimed the opposite, that Gao fits across geometry at fixed blowing
ratio, inferred from their sentence that "the correction coefficients ...
increased as the blowing ratio decreased". That was an over-reading of one
sentence against an explicit statement of method elsewhere.

## The published pair cannot reach what this plate needs

Not a fitting failure -- arithmetic on the constants. Since
`a r / (a r + 1) >= 0` for `a > 0` and `r >= 0`:

    alpha >= b = 0.9465    for EVERY r

| M | alpha Murray needs | below the floor by |
|---|---|---|
| 0.19 | 0.850 | 0.097 |
| 0.48 | 0.864 | 0.083 |
| 0.96 | 0.686 | 0.260 |

All three sit below the floor, so no choice of `r` reaches them. The
published pair can damp a row by at most **5.35%** while staying at or
below 1; this plate needs 15% to 31%.

**This argument is deliberately scale-free, and that matters.** A first
version claimed the coefficients give `alpha > 1` above M = 0.20 on
Murray's rig, from an assumed `r = M x A_holes_row / A_duct`. Withdrawn:
Gao's test-section dimensions are not in the extractable text, so `r`
cannot be put on their scale. Note the asymmetry that caught me out -- the
*fit* below is legitimately invariant to `r`'s scale, because `r -> k r` is
absorbed by `a -> a/k`, but applying PUBLISHED coefficients is not.

**Checked and not the cause: the `M` versus `Me` distinction.** Gao's
equivalent blowing ratio `Me = M0 A0/Ae` (Eq. 10) normalises the single-row
efficiency lookup `eta = f(X/(M s))` onto a baseline hole configuration. It
is not alpha's `r`, which is `m_i/m_g` from the Eq. (3) energy balance, and
the quoted calibration values 0.3 and 1.0 are `M`. combaero needs no `Me`
equivalent because Baldauf takes lateral spacing as an explicit argument
where Gao's single-row database is fixed at one spacing.

**The likely cause is structural.** Eq. (5) has no streamwise-spacing term:

| rig | spanwise | streamwise |
|---|---|---|
| Gao Case 1 (calibration) | 3.5 d | 10.5 d |
| Gao Case 4 (tightest) | 2.8 d | 8.4 d |
| Murray 5.75D staggered | 5.75 D | **2.875 D** |

Murray's rows are 3.7x closer than Gao's calibration plate, and Murray
proves streamwise spacing is what drives superposition error -- tripling it
to 17.25 D restores plain superposition to under 10%. Two plates with the
same coolant fraction and different row spacing get the same alpha from
this form.

Gao's own limits are consistent: "in Case 4, ... when the blowing ratio
exceeded 0.64, the prediction accuracy decreased compared to the
traditional Sellers model" -- their tightest plate, above M = 0.64. Murray
at M = 0.96 and 2.875 D is further into that corner than anything Gao
tested.

## A one-knob correction works

For one geometry `r` is fixed per series, so alpha collapses to ONE number
per blowing ratio. The fit is then a scan over a single parameter -- no
optimiser, no local minima, and the objective's shape is visible.

| M | Sellers | best alpha | MAE | best C | MAE |
|---|---|---|---|---|---|
| 0.19 | 33.1% | 0.850 | 15.8% | **2.15** | **13.7%** |
| 0.48 | 43.6% | 0.864 | 21.4% | **2.16** | **17.1%** |
| 0.96 | 112.6% | 0.686 | 23.7% | **4.29** | **16.2%** |

At M = 0.96 a single knob takes the error from 113% to 16%, which is close
to the paper's own stated 15% experimental uncertainty.

The minima are genuine, not plateaus -- checked, because a flat objective
would mean the "best" value is not a measurement:

    M0p19   alpha  0.70=25.4%  0.80=17.7%  0.85=15.8%  0.90=17.5%  0.95=23.6%
    M0p48   alpha  0.70=32.7%  0.80=23.5%  0.85=21.5%  0.90=24.3%  0.95=33.0%
    M0p96   alpha  0.65=24.3%  0.70=23.8%  0.75=27.9%  0.85=53.4%  0.95=89.9%

M = 0.19 and M = 0.48 share a minimum at alpha 0.85, so their apparent
0.850-vs-0.864 difference sits inside the flat region and is not real.

## The coupling form fits better, and is sourced

Gao's own survey of the "third category" of corrections gives the form:

> "scholars established a general form: eta = eta1 + eta2 - C eta1 eta2,
> where C is adjusted to modify the predicted cooling efficiency."

Applied recursively over rows. **C = 1 is exactly Sellers**, verified to ten
digits against `film_superposition_sellers` -- the same identity property
that makes alpha = 1 the classical model.

It beats alpha at every blowing ratio here (13.7 / 17.1 / 16.2 against
15.8 / 21.4 / 23.7), and most at M = 0.96 where the correction matters.
Its fitted values 2.15, 2.16, 4.29 are monotone in blowing ratio.

**And critically it is parameterised along the axis Murray varies.** Gao
attributes the M dependence to Xu et al. [36], who "suggested that the
correction factor for turbine blades depends on both the blowing ratio, M,
and the dimensionless distance, X/D". Others in the same family: Huo et al.
[26] add dimensionless hole pitch; Zhang et al. [37] drop X/D and make C
depend on hole-row count and hole Reynolds number; Zhang et al. [38] handle
>50 rows via a spanwise non-uniformity factor and streamwise dissipation
rate.

Gao's stated reason for not using it: the third category "lacks a
consideration of streamwise temperature variations and requires extensive
databases to handle continuous changes in the hole diameter". That is a
fair objection for their design-iteration purpose. It does not apply to
scoring one fixed geometry.

## What is NOT shipped

**No fitted `(a, b)`, and no fitted `C`.** The values above are tuned to
Murray's 5.75D staggered plate and to no one else's. Matching a rig is the
user's job and the harness never applies a tuner.

**No implementation of the C form yet.** It is one line of recursion and it
is sourced, so it is a reasonable candidate -- but Xu's actual `C(M, X/D)`
is in a paper this repo does not have, and shipping the form with `C` as a
required argument would repeat exactly the `(a, b)` situation: a shape with
no coefficients. Worth doing only alongside a source for `C`, or as an
explicit user knob.

**No invented replacement.** A monotone-decreasing alpha or a
blowing-indexed table would both fit better and both be made up.

## The confound that no fit on this family can escape

Every point shares one plate, so `r` and `M` move together. Nothing fitted
here can tell which of the two a correction actually responds to -- and
more blowing ratios on the same plate would not help. That is precisely why
Gao used four plates.

What would settle it:
- data at a **different geometry** at the same blowing ratios -- Murray's
  own 3.0D plate would do it, but Figure 4 is contours with no
  superposition reference;
- Gao's `(a, b)`, if they ever surface;
- Xu et al. [36] for the `C(M, X/D)` form;
- Murray's Figure 9a at 17.25D streamwise, where superposition nearly
  works: any correction must approach its identity value there (alpha -> 1,
  C -> 1), and one that cannot is wrong regardless of its fit at 5.75D.

Related: `film_superposition.md`,
`murray_ireland_2018_effusion_superposition.md`, and #420's finding that
`andrei2014` cannot support this fit at all.
