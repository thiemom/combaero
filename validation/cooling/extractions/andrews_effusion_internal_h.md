# Effusion plate internal heat transfer -- Andrews 86-GT-225

Andrews, G.E., Alikhanizadeh, M., Asere, A.A., Hussain, C.I., Khoshkbar
Azari, M.S. and Mkpadi, M.C. (1986). "Small Diameter Film Cooling Holes:
Wall Convective Heat Transfer." ASME 86-GT-225. Department of Fuel and
Energy, University of Leeds.
`docs/heat_transfer/film/V004T09A022-86-GT-225.pdf`

This closes #387's internal side. The external side is
`baldauf_2002_film_effectiveness.md` and `film_superposition.md`.

## Why this paper and not the two it rests on

Andrews 88-GT-290 computes its wall heat transfer from two correlations it
does not derive: **Sparrow** for the coolant-side hole approach and
**Mills** for the short-hole throat. 86-GT-225 curve-fits both, so neither
original paper is needed -- which is what made this implementable from what
the repo already has.

## The physical picture

A coolant hole cools its wall in two places:

1. the **approach** flow over the coolant-side plate surface, converging
   into the hole;
2. the **throat**, a short tube with a sharp-edged entry.

Andrews' headline is that the first dominates: "the hole approach flow heat
transfer is much larger than the internal hole heat transfer". He sums
them -- "the authors have treated the wall heat transfer as the summation
of Equations 13 or 14 and 15" -- after Eq. (18) puts both on the same
Nusselt definition.

## The equations, as transcribed

| | |
|---|---|
| Eq. (12) | `Nu = 0.023 Re^0.8 Pr^(1/3) R_Nu` -- Mills, throat |
| Eq. (13) | `R_Nu = 0.13 z^3 - 0.75 z^2 + 1.04 z + 2.24`, `z = L/D <= 2` |
| Eq. (14) | `R_Nu = 1 - 49.3 w^4 + 58.6 w^3 - 26.5 w^2 + 7.48 w`, `w = D/L`, `L/D > 2` |
| Eq. (15) | `Nu_l = 0.881 Re^0.476 Pr^(1/3)` -- Sparrow, on the surface dimension `l` |
| Eq. (16) | `l = plate area / hole pitch` |
| Eq. (18) | `Nu = 0.881 Re^0.476 Pr^(1/3) X/(pi L)` -- Eq. 15 rebased onto the hole internal area |
| Eq. (19) | `Nu = (0.27 Re^0.476 + 0.023 Re^0.8 R_Nu) Pr^(1/3)` -- summed, for one worked geometry |

`Re` is on the hole diameter throughout.

## THE SCAN IS POOR, SO EVERY EQUATION WAS CHECKED THREE WAYS

Text extraction returns things like `RNu 0.13(n) - 0*75(r) + 1.04(r) +
2.24` and `Nu = 0.023 Re03Prel'33RNu`. The page was rendered at 400 dpi and
read visually, and then each reading was tested against something the
transcription could not fake:

1. **The two R_Nu branches must meet at L/D = 2.** They are independent
   polynomials in different variables -- `L/D` and `D/L` -- and the
   candidate reading gives 2.3600 and 2.3588, agreeing to **1.3e-3**. A
   misread coefficient would not reproduce that.
2. **R_Nu must approach 1 as L/D grows.** A long hole is fully developed
   and has no entrance enhancement. Eq. (14) gives 1.54 at L/D = 10, 1.14
   at 50, 1.015 at 500.
3. **Eq. (19)'s printed 0.27 must equal 0.881 X/(pi L).** For Andrews' own
   comparison geometry, X = 6.11 mm and L = 6.35 mm, that is **0.26983**
   against the printed 0.27. This is the only place the paper evaluates
   Eq. (18)'s geometry factor numerically, so it is the one arithmetic
   check available on the approach term -- and it passes.

Eq. (16) needed judgement: the typesetting shows `l = plate area/hole pitch
= X^2 - pi D^2/4`, whose right-hand side is an AREA. The paper then says
`l` "is approximately the same as the hole pitch, X" for small holes, which
only holds for `(X^2 - pi D^2/4)/X`. For plate C that is 14.69 mm against
X = 15.24 mm. Read as the division.

## Scoring against 88-GT-290's Figure 8

`andrews1988/fig8_h_effusionC`, 10 points, plate C: D = 3.27 mm,
X = 15.24 mm, thickness 6.3 mm, so L/D = 1.93 -- the Eq. (13) branch.

Two conversions stand between Eq. (19) and that figure's axes, and both are
pure geometry:

    Re = 4 (G X^2) / (pi D mu)        coolant per hole from G
    h  = Nu k / D * A_h / A           hole-area Nu to plate-area h

**`A/A_h` is 3.46 for this plate**, matching the paper's tabulated 3.4, so
omitting the second conversion would overstate `h` by that factor. That is
a definitional error of the #389 class, not a modelling one, and it is
pinned by its own test.

**`G` is per GROSS plate area.** The nomenclature separates "plate area"
(for G) from "hole approach surface area A" (for h), and the hole count
settles it: N = 4306 per square metre against 1/X^2 = 4305.8.

### Result

| | |
|---|---|
| bias | **-9.9%** |
| RMSE | 11.6% |

(re-scored 2026-10-08 after #485: air's viscosity was 11.8% low before). It was -13.5% bias and 14.7% RMSE with the low viscosity.
| basis | **accuracy** (cross-source) |

**This is ACCURACY, not fidelity, and the distinction was nearly missed.**
Same author and same group invites "fidelity", but the correlations are
86-GT-225 and the data is 88-GT-290 -- same lab, different study, which the
validation policy counts as cross-source. The first `SET_ORIGIN` entry used
the bare author name, which matched the `andrews1988` source and reported
fidelity; it now names `86-GT-225` so the label is right.

The miss is **one-sided**: every point under-predicts. That is a systematic
gap, not scatter, and it should not be described as "within scatter" and
left there -- though it does sit inside the source's own +/-16% spread
about its power-law fit, and the two worst points (-28% and -20%) are the
two the `andrews1988` metadata already flags as off-trend in the source.

### What was ruled out as the cause

- **Coolant temperature.** Not stated in the paper. Across 280-320 K the
  bias moves -11.8% to -8.1%, so the assumption is worth about 4 points
  and cannot explain 9.9%. Reported by `temperature_sensitivity`, never
  tuned -- picking the temperature that scores best would be fitting an
  unmeasured input to the metric judging it.
- **Property values.** A first spike used hand-typed air properties and got
  -8%; combaero's own properties then gave -13.5%. This note used to blame
  the hand-typed values (closer to 400 K air). It had it backwards: the
  library's air viscosity was 11.8% low, from flagged transport data
  (#485). With that fixed combaero gives -9.9%, close to the hand-typed
  -8%. The lesson stands, inverted: a property swing this size is a reason
  to check the library's properties against literature, not to trust them.
- **The `G` definition.** Per net area rather than gross changes the flow
  by 3.6% and moves the bias the wrong way.

### Not yet ruled out

The remaining candidates need work this record does not do:

- Sparrow's multi-hole correlation is, in Andrews' own words, "based on a
  relatively small number of geometries", and plate C may sit outside them.
- Andrews notes the two Sparrow papers -- single-hole (21) and multi-hole
  (22) -- "are incompatible with each other", and only the multi-hole one
  is used here.
- 88-GT-290 Fig. 8 carries two more plates (A, A/C) that are digitised but
  not committed. The 1986 plates below answer the geometry question for
  the 1986 rig; they do not answer it for the 1988 one.

## Scoring against 86-GT-225's own Figure 8 -- the FIDELITY basis

`andrews1986/fig8_hm_plate_{a,b,c,d}`, 41 points over the four Table 5
plates. Same paper as the correlations, so this is fidelity where the 1988
comparison above is accuracy, and the scorecard reports them as two rows.

### THE TWO PAPERS DO NOT PLOT THE SAME h

This is the error worth the most here, and nothing in the numbers says
which convention is meant. The nomenclatures:

| paper | symbol | printed definition | basis |
|---|---|---|---|
| 86-GT-225 | `h_m` | "Average heat transfer coefficient, W/m2K, over the hole length" | hole internal, `pi D L` |
| 88-GT-290 | `h` | "Convective heat transfer coefficient based on the surface area, A" | plate, `A = X^2 - pi D^2/4` |

`A/A_h` runs 1.4 (plate d) to 9.8 (plate a). Scoring the 1986 points on
the plate-area convention gives **-60%**; on the hole-length convention,
**-10.4%**. The two are carried as distinct `y_axis` values -- `h_internal`
and `h_hole_length` -- so the runner cannot pick the wrong one silently,
and `test_the_two_papers_h_conventions_are_different_quantities` asserts
both the ratio and the collapse.

86-GT-225's `h_m` needs NO area conversion at all: Eq. (19)'s Nusselt
number is already on the hole internal area, so `h_m = Nu k / D`.

### Geometry: Table 5, with D recovered

The scan loses Table 5's `n` and `D` columns. `D = X/(X/D)` recovers them,
and two independent things then check the recovery:

| plate | D mm | L mm | L/D | X mm | X/D | array | technique |
|---|---|---|---|---|---|---|---|
| a | 1.178 | 6.35 | 5.38 | 15.20 | 12.9 | 10 x 10 | Drilled |
| b | 0.642 | 6.35 | 9.92 | 6.10 | 9.5 | 25 x 25 | Spark Eroded |
| c | 0.897 | 6.35 | 7.05 | 6.10 | 6.8 | 25 x 25 | Drilled |
| d | 1.298 | 6.35 | 4.85 | 6.10 | 4.7 | 25 x 25 | Drilled |

1. `6.35/D` reproduces each printed `L/D` to 1%.
2. The four `L/D` land on Figure 10's four abscissae (4.81, 5.35, 7.10,
   9.90) to 1% -- a different figure, digitised in a different session.

Plate b is the `X = 6.11, D = 0.64, L = 6.35` geometry Eq. (19)'s printed
`0.27` is evaluated for.

### THE AXES ARE LOG AND WERE DIGITISED ON A LINEAR CALIBRATION

Caught, as the two before it, by a calibration file failing to reproduce
round numbers. `fig8_calibration_as_read.csv` is committed UNCONVERTED
because it is the evidence: the x ticks at true 0.5 and 1.0 read 1.1299 and
1.5689, and a linear pick of a log axis spanning 0.1 to 2 predicts 1.1218
and 1.5619. The committed data CSVs are converted through

    true = t0 * (t1/t0) ** ((read - l0) / (l1 - l0))

anchored on the frame, which recovers the two ticks it does NOT use to
+1.3%/+1.1% (x) and -0.8%/-1.5% (y). A third channel confirms it: the paper
states G ran "from 0.1 to 1.7 kg/sm2" and the converted abscissae span
0.089 to 1.79.

Page skew between the two top corners is 0.21% of value and moves a
mid-range point by 0.39% -- recorded, not corrected, being an order of
magnitude below the residual.

### Result

| | |
|---|---|
| bias | **-7.1%** |
| MAE | 9.2% |
| n | 41 |

(re-scored 2026-10-08 after #485: air's viscosity was 11.8% low before); it was -10.4% bias, 11.1% MAE.
| basis | **fidelity** (the correlations' own paper) |

Per plate: a -2.4%, b -16.5%, c -2.0%, d -9.3% (before #485: a -6.2%, b -19.6%, c -5.5%, d -12.5%).

### THE FIDELITY MISS IS NOT A TRANSCRIPTION ERROR

A 10% miss on the author's own data would normally be our bug -- that is
what the fidelity label means. Here it demonstrably is not, and Figure 10
is why: it plots Andrews' own evaluation of his own Eq. (19), and combaero
reproduces it to 1.5%. Fig. 10's Re = 2200 sits inside Fig. 8's range (18
points below, 23 at or above), so the transcription check covers the
regime the miss is measured in.

So the two figures of one paper separate the two failure modes cleanly:
**our transcription is right to 1.5%, and Andrews' correlation
under-predicts Andrews' own measurements by about 7%.** That is a stronger
statement than either figure alone supports, and it is what made this
digitisation worth doing.

### What the residual is NOT explained by

- **Reynolds number** (band figures from before #485; not re-run). Pooled bias by band: -20.4% below Re 500 (n = 4),
  -6.9% for 500-1500, -9.1% for 1500-3000, -11.6% above 3000. No monotone
  trend, and the low-Re band is four points from three different plates
  (b 405, c 288, d 199, d 404) -- too few to read as a regime.
- **L/D.** Ordered by L/D the biases run d 4.85 (-9.3%), a 5.38 (-2.4%),
  c 7.05 (-2.0%), b 9.92 (-16.5%) -- not monotone either.
- **Pitch.** Plate a is the only plate at X = 15.2 mm and is the second
  best fit, so the `X/(pi L)` term is not carrying the error.

### An observation, NOT a demonstrated cause

**Plate b is the worst fit and is the only SPARK ERODED plate in Table 5.**
Spark erosion leaves a different bore finish and entry geometry than
drilling, and the throat term assumes a sharp-edged entry -- a plausible
mechanism. It is recorded and nothing more: one plate is not evidence, no
correction is applied, and the fit is not improved by excluding it (pooled
bias without plate b is -8.2% over 33 points, which is a smaller number
for a smaller dataset, not a better model).

### Unlike the 1988 comparison, this is not one-sided

The worst positive error is +3.4%, so this is a strong negative bias with
scatter rather than a pure offset. Worth distinguishing, because the 1988
test asserts one-sidedness and this one must not.

## Deliberately not implemented

**The temperature chain.** An earlier note in this repo said to "chain them
as Andrews does: approach outlet temperature becomes throat inlet
temperature". The paper does no such thing -- Eq. (19) is a plain sum of
Nusselt numbers, and the text says so. The chain idea appears to have been
an inference; it is not in 86-GT-225.

**A plate-area `h` convenience function.** The conversion needs only
geometry and no correlation, so it belongs with the caller that knows the
geometry rather than in the correlation library.
