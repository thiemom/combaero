# Extraction: Colban, Thole and Bogard (2011), shaped-hole film effectiveness

**Status: IMPLEMENTED and REVIEWED 2026-10-21** as
`combaero::cooling::film_effectiveness_colban_2011` (#386). One closed-form
equation and three coefficients; the implementation lands on the line the
authors draw within about a pixel. Signed off by the user on the items
below, with their answers.

## Review questions and answers

1. **Is `D` the metering diameter?** The nomenclature says "D = hole throat
   diameter" and Section 4 substitutes `A_in = pi/4 D^2`, so D is the inlet
   (metering) bore and M is referenced to the metering velocity (Eq. 7).
   The implementation and docstrings say so. *Turns on it:* a caller using
   the exit width would shift xi by the area ratio, a factor of 2.5-4.7.
   **Reviewer: yes** -- Fig. 2 sketches D as the metering diameter too.
2. **Should alpha and density ratio be inputs at all?** They do not enter
   Eq. (19). They are taken only to flag extrapolation, because every data
   set was at 30 deg and DR 1.7-2.0 (Sections 3.4, 3.8). The alternative --
   leaving them out -- would let a DR = 1 lab case read as in-domain.
   **Reviewer:** neither is absent from the physics, only from the form.
   - alpha is carried by the hole geometry it produces: the breakout length
     and so t and AR.
   - Density ratio enters through the blowing ratio, M = DR * VR.
   - The authors present Eq. (19) as general on that basis. Flagging
     outside the fitted 30 deg and DR 1.7-2.0 is fine, but read the flag as
     "outside the fitting data", not "invalid" -- the authors appear to
     accept use beyond it.

   **Decision:** both stay range-check inputs; the flag is information and
   changes no value.
3. **Fidelity basis.** The scored points are Gritsch et al. (2005), one of
   the four studies Eq. (19) was fitted on, with Colban a co-author. The
   scorecard labels them fidelity, which is right by the policy (the
   correlation's own data) but optimistic: it is fit quality, not prediction.
   **Reviewer: agreed.**

## What I still need sent

Nothing to close this extraction. Figs. 16-19 (pp. 10-11) plot Eq. (19)
against further fitting data -- TTCRL at M 0.6/1.5/2.5, and t/P, AR and P/D
sweeps -- with several symbol sets per panel; digitising them by hand would
widen the fidelity basis from one geometry to five.

## Source, pinned

> Colban, W.F., Thole, K.A. and Bogard, D. (2011). "A Film-Cooling
> Correlation for Shaped Holes on a Flat-Plate Surface." *ASME J.
> Turbomachinery* **133**(1), 011002. DOI 10.1115/1.4002064.
> `sources/heat_transfer/film/2011-Colban-JT-1fygyl9.pdf` (gitignored).

One row of fan-shaped (FS) or laidback fan-shaped (LFS) holes on a flat
plate, laterally averaged adiabatic effectiveness, low freestream turbulence,
negligible pressure gradient.

## Extracted

| Item | Value | Where | State |
|---|---|---|---|
| Eq. (19) | `eta = 1 / (P/t + C1 M^C2 xi^C3)` | p. 5 and Table 2 | confirmed |
| Eq. (18) | `xi = (4/pi) (X/D)(P/D) / (M AR)` | p. 5 | confirmed |
| Eq. (12) | `eta(X/D = 0) = t/P` | p. 4 | confirmed |
| C1, C2, C3 | 0.1721, -0.2664, 0.8749 | Table 3 | confirmed |
| M limits | 0.2-2.5 | Table 4 | confirmed |
| t/P limits | 0.31-0.65 | Table 4 | confirmed |
| AR/(M P/D) limits | 0.17-1.17 | Table 4 | confirmed |
| Fitted spans | AR 2.5-4.70, P/D 4-8, t/P 0.31-0.75, M 0.2-4.0 | Section 5 | confirmed |
| Fitted conditions | alpha 30 deg, DR 1.7-2.0, Tu 0.2-4%, Ma_ext <= 0.3 | 3.4, 3.8, Table 1 | confirmed |
| eta definition | `(T_ad - T_jet)/(T_inf - T_jet)`, lateral mean | Eqs. (1)-(2) | see note |

**eta convention.** Colban writes `(T_ad - T_jet)/(T_inf - T_jet)`, which is
`1 - (T_inf - T_ad)/(T_inf - T_jet)` -- the complement of the usual form as
printed. Read against the data it cannot be: eta is t/P at the hole and
falls downstream, so it is the cooled fraction, `(T_inf - T_ad)/(T_inf -
T_jet)`, the same quantity as Baldauf's and `adiabatic_wall_temperature`'s.
Recorded as a printing slip in Eq. (1), not a different convention -- the
printed form, confirmed on the page at 300 dpi, would give eta = 0 at the
breakout where Section 3.5 and Eq. (12) put "eta = 1 across the hole
breakout".

The paper's own accuracy statement is Table 5: 46 of 57 data sets fall
inside the Table 4 limits and 93.5% of those have R^2 > 0.7. No RMS is given.

## The check that could have failed

The paper draws Eq. (19) for Gritsch configuration B (AR 3.5, P/D 6,
t/P 0.49) at M = 0.5 (Fig. 6) and M = 2.5 (Fig. 7). Its line was read off
the 300 ppi raster, raw, at 28 stations, independent of the implementation:

| | stations | worst | pixel |
|---|---|---|---|
| Fig. 6, x/D 10-47.5 | 12 | < 0.003 | 0.0008 |
| Fig. 7, x/D 3-48 | 19 | < 0.003 | 0.0013 |

Median offset over every clean column is +0.0002 (Fig. 6) and -0.0001
(Fig. 7). Perturbing C3 by 1.1% fails the test, so it discriminates.
Fig. 6 below x/D 10 is crossed by the Eq. 20/23 lines and is not read.

## On its own fitting data

The circles in the same figures are Gritsch et al. (2005) data, inside the
Table 4 limits:

| | points | bias | RMS | pattern |
|---|---|---|---|---|
| M = 0.5 | 15 | +1.6% | 3.6% | uniform |
| M = 2.5 | 17 | +4.8% | 12.1% | low by 15% near the hole, high by 10-20% past x/D 20 |

Errors are predicted/measured - 1, as the scorecard reports them. The M = 2.5
shape is what the paper's Fig. 7 shows and is the correlation's, not the
digitisation's: the line and the circles both reproduce. The x/D = 2 point
at M = 2.5 reads 0.514, above t/P = 0.49, which the paper attributes to
near-hole conduction error (Section 7).

## Modelling decisions

- **Shaped holes only.** Not merged with Baldauf: the paper's own Fig. 7
  shows the cylindrical correlation an order of magnitude low at M = 2.5.
  A caller chooses the hole type; nothing switches between them.
- **No DR, alpha or Tu dependence invented.** The fit never varied them, so
  they flag extrapolation and change nothing. alpha acts through the
  geometry (t, AR) and DR through M, which is why the authors call the form
  general. Revisit if a shaped-hole source with DR or angle variation is
  added.
- **Only M carries a derivative.** Geometry is fixed per element and DR is
  absent from the equation.

## Review log

- 2026-10-21: extracted from the PDF text layer and the page at 400 dpi;
  Table 3 and Table 4 read in both; Eq. (1) checked on the page image. Implemented with the figure check above.
- 2026-10-21: signed off by the user. Answers are under "Review questions
  and answers"; the Fig. 6/7 redraws were judged well drawn.
