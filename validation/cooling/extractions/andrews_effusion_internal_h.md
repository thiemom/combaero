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

## Scoring against Andrews' own Figure 8

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
| bias | **-13.5%** |
| RMSE | 14.7% |
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
  bias moves -15.3% to -11.7%, so the assumption is worth about 4 points
  and cannot explain 13.5%. Reported by `temperature_sensitivity`, never
  tuned -- picking the temperature that scores best would be fitting an
  unmeasured input to the metric judging it.
- **Property values.** A first spike used hand-typed air properties and got
  -8%; combaero's own properties at the same temperature give -13.5%. The
  hand-typed values were closer to 400 K air. Recorded because the 6-point
  swing came entirely from not using the library's own properties.
- **The `G` definition.** Per net area rather than gross changes the flow
  by 3.6% and moves the bias the wrong way.

### Not yet ruled out

The remaining candidates need work this record does not do:

- Sparrow's multi-hole correlation is, in Andrews' own words, "based on a
  relatively small number of geometries", and plate C may sit outside them.
- Andrews notes the two Sparrow papers -- single-hole (21) and multi-hole
  (22) -- "are incompatible with each other", and only the multi-hole one
  is used here.
- Fig. 8's other two plates (A, A/C) are digitised in the same figure but
  not committed; scoring them would show whether the -13.5% is geometry
  specific or uniform.

## Deliberately not implemented

**The temperature chain.** An earlier note in this repo said to "chain them
as Andrews does: approach outlet temperature becomes throat inlet
temperature". The paper does no such thing -- Eq. (19) is a plain sum of
Nusselt numbers, and the text says so. The chain idea appears to have been
an inference; it is not in 86-GT-225.

**A plate-area `h` convenience function.** The conversion needs only
geometry and no correlation, so it belongs with the caller that knows the
geometry rather than in the correlation library.
