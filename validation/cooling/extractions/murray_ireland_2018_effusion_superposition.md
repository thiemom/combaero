# Murray & Ireland (2018) -- what it settles about superposition

Murray, A.V., Ireland, P.T., Wong, T.H., Tang, S.W. and Rawlinson, A.J. (2018).
"High Resolution Experimental and Computational Methods for Modelling Multiple
Row Effusion Cooling Performance." *Int. J. Turbomachinery, Propulsion and
Power* **3**(1), 4. Open access (MDPI, CC BY).
`docs/heat_transfer/film/c32eeeb88f403f5f28556753359153c1b85b.pdf`

**Read 2026-09-30, before digitising anything.** No data is extracted yet.
This record exists because the paper already answers the question #420 left
open, and the answer changes what the next extraction is *for*.

## Why this source was sought

#420 scored Baldauf + Sellers against `andrei2014` and it failed: MAE 45.5%,
with the error changing sign from +61% at BR 1 to -81% at BR 3. Two
explanations were live -- the superposition, or the per-row closure -- and
Andrei alone cannot separate them.

## What the paper does

Two flat plates, primary hole pitches **3.0D and 5.75D**, PSP via the
heat-mass transfer analogy, blowing ratios **0.1-1.2**, turbulence intensity
**4.43-4.72% (stated)**, density ratio ~1. Plate depth 2.5D. Geometry 1 is a
7x8 primary array (112 holes with staggered), Geometry 2 a 4x5 array (40).
Stated experimental error **does not exceed 15%**, worst at the lowest
blowing ratios where the mass-flow controllers lose precision.

Crucially it also applies **Sellers' superposition to a single-hole CFD
result** and compares against both multi-hole CFD and the PSP experiment.
That is the same chain combaero implements, with a CFD single-hole
distribution in place of Baldauf -- so it isolates the superposition from
the closure.

## The finding

1. **Sellers OVER-PREDICTS at effusion pitches, by about a factor of two.**
   Verbatim: "by a blowing ratio of approximately one, spanwise averaged film
   effectiveness values via the superposition method were around twice those
   displayed by the multi-hole CFD simulation and PSP experiment."
2. **The error grows with blowing ratio** -- "the deviation in effectiveness
   becoming more prominent as the blowing ratio was elevated."
3. **The cause is streamwise jet interaction, and they proved it by
   experiment.** Holding spanwise pitch at 5.75D and tripling the STREAMWISE
   pitch to 17.25D brings superposition back to "less than 10% discrepancy by
   the third row of holes". Velocity and turbulent-kinetic-energy contours
   (Figs. 10, 11) show the flow field re-establishing between rows at the
   wider pitch and not at the tighter one.
4. Their conclusion: "without alteration to the superposition method to
   account for variations in the flow field, its applicability to effusion
   cooling with closely pitched streamwise holes is somewhat limited."

## What this settles for combaero

**The two error sources in #420 are now separable, and both are identified.**

| error | direction | grows with | cause | evidence |
|---|---|---|---|---|
| superposition | OVER-predicts | blowing, row count, tighter streamwise pitch | streamwise jet interaction | Murray, CFD single-hole -- no correlation involved |
| per-row closure | UNDER-predicts above M ~ 1 | blowing | Baldauf's lift-off collapse | #420, `test_film_effectiveness.cpp` |

`andrei2014`'s numbers fit this exactly. Its streamwise pitch `s_x/d` = 9.15
sits between Murray's failing 5.75D and working 17.25D, so a moderate
superposition over-prediction is expected -- and +61% / +32% at BR 1 is what
was measured. At BR 3 the closure collapse (-81%) takes over and masks it.

**This qualifies #420's conclusion that Gao's `(a, b)` cannot be fitted.**
That remains true *on Andrei*, but the reason is narrower than stated: a
closure artefact contaminates the superposition measurement at BR 2-3. Gao's
alpha is exactly the right SHAPE for the superposition error -- bounded
[0, 1], only reducing, per-row, growing with accumulation -- and Murray's
independent CFD confirms the error it corrects is real and is an
over-prediction. Murray's blowing range 0.1-1.2 is entirely below Baldauf's
lift-off peak, so it is the range where alpha can be measured without the
closure contaminating it.

## What to extract, and why

**CORRECTION, 2026-09-30.** An earlier draft of this record said Figure 6 was
the 3.0D pitch and therefore inside Baldauf's `s/D` envelope. Both halves
were wrong, found by opening the figure rather than trusting the caption --
the caption states no pitch, and the panel titles read **S = 5.75D**. What
the paper actually offers, per pitch:

| pitch | figure | what is on it | superposition shown? |
|---|---|---|---|
| **3.00D** | Fig. 4 | contours, M = 0.19 / 0.47 / 0.93 | **NO** -- experiment and multi-hole CFD only |
| **5.75D** | Fig. 6 | spanwise-averaged eta vs x/D, M = 0.19 / 0.48 / 0.96 | **yes**, all three curves |
| **5.75D** | Fig. 5 | contours, three blowing ratios | no |
| 17.25D streamwise | Fig. 9a | spanwise-averaged eta vs x/D, M = 0.48 | yes |
| both + Ling 10D/16D | Fig. 7 | average eta vs normalised mass-flow | no |

**So there is no 3.0D spanwise-averaged plot, and no 3.0D superposition
reference at all.** The clean in-envelope test hoped for does not exist in
this paper.

**Figure 6 is still the target**, but on a narrower claim:

- it carries the PSP measurement, a CFD reference AND the superposition
  prediction on one axis, so the superposition error is read directly rather
  than inferred -- no other source does this;
- at **M 0.19-0.96 it is entirely below Baldauf's lift-off peak** (M ~ 0.7-1.0
  at these spacings), which was the confound that wrecked `andrei2014` at
  BR 2-3;
- `s/D` = 5.75 is **outside** Baldauf's envelope of 2-5, but by 15% against
  `andrei2014`'s 47% at 7.37. Extrapolation is reduced, not eliminated, and
  the series must still be flagged.

Figure 9a (17.25D streamwise) is the control: superposition should come good
there, and a fitted alpha should approach 1.

Figure 4's 3.0D contours could in principle be spanwise-averaged by image
analysis to recover a line plot, but that is a far harder extraction than
digitising a curve and it would still lack a superposition reference. Not
recommended as a first step.

**`uncertainty` for these series is 0.15** -- the paper's own stated
experimental error -- which is a MEASUREMENT band as #389 requires, not a
model error. Note it is an upper bound stated for the worst case (lowest
blowing), so using it uniformly is conservative; if a per-blowing-ratio
figure is readable from the paper, prefer it and record which.

## Extraction will be harder than andrei2014

Every figure page is a **raster image** -- checked, zero vector path
operators on pages 7-13. So this needs pixel digitisation with the
automated-pass + manual re-digitisation + cross-check workflow, not the exact
vector extraction `andrei2014` allowed. Budget accordingly, and do not claim
`extraction: pdf-vector-exact` for it.

## Not yet read

- Equation (6), the paper's own statement of Sellers -- worth checking
  against `film_superposition_sellers` for convention, though the form is
  standard.
- Ling et al. [27], the 10D and 16D pitch comparison in Fig. 7. A third and
  fourth pitch, if that data is reachable, would make the pitch dependence a
  curve rather than three points.
- Whether the blowing ratio quoted per figure is the plate average or
  per-hole. The paper is explicit that M "was not constant" across the
  multi-hole plate and that superposition "inherently disregards" this --
  which is itself one of their hypothesised causes, so it matters for how a
  scored comparison is labelled.
