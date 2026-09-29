# Extraction: Lichtarowicz, Duggins and Markland (1965), long orifices

**Status: IMPLEMENTED 2026-09-29** as
`DischargeCdCorrelation::Lichtarowicz1965`, replacing the member that had
declared-and-refused since #409.

Equations read off the page rendered at 600 dpi with `pdftoppm -r 600`, not
the PDF text layer, which garbles both (Eq. 7 OCRs as `0+327-0.00851/d`).

## Source, pinned

> Lichtarowicz, A., Duggins, R.K. and Markland, E. (1965). "Discharge
> coefficients for incompressible non-cavitating flow through long
> orifices." *J. Mech. Engng Sci.* **7**(2), 210-219.
> `docs/orifices/lichtarowicz-et-al-1965-...pdf` (gitignored, copyrighted).

## Why it earns a place -- the coverage criterion

A source is worth adding only if it plugs a hole no other source covers:

| correlation | Re range | l/d range |
|---|---|---|
| McGreehan-Schotsch (1988) | >= 1e4 | long holes, WITH crossflow |
| Idelchik (1966) 4-18a | 25 to 1e6 | l/Dh up to 4 |
| **Lichtarowicz (1965)** | **10 to 2e4** | **2 to 10** |

It extends l/d from 4 to 10 and, more importantly, it is the only one of the
three that is *about* long holes at low Reynolds number. **That is where a
cooling hole actually runs**: Andrews' own effusion plate C spans Re 432 to
8718 across its measured range, entirely inside McGreehan's floored region.

## The equations

**Eq. (7)** -- the ultimate (high-Re) discharge coefficient, stated to 1.5%:

    C_du = 0.827 - 0.0085 l/d           for 2 <= l/d <= 10
    C_du = 0.810                        for 1.5 <= l/d < 2

**Eq. (12)** -- the low-Re form, "fits all but a few points to better than
0.02 in the range of l/d from 2 to 10 and of Re from 10 to 2 x 10^4":

    1/Cd = 1/C_du + (20/Re)(1 + 2.25 l/d)
           - (0.005 l/d) / (1 + 7.5 (log10(0.00015 Re))^2)

Re is based on the orifice diameter. The paper's `Cd = Q(1-m^2)^(1/2) /
A(2gh)^(1/2)` reduces to our plenum definition as the area ratio m -> 0, so
the conventions are compatible and this is not a #389 case.

## Decisions

**D1. Refuses below l/d = 1.5, rather than clamping.** This is the source's
own design recommendation (1): "Avoid l/d less than 1.5, since the discharge
coefficient varies rapidly with l/d below this value, and **there is the
possibility of hysteresis in operation**." A single-valued correlation cannot
represent hysteresis, so returning a number there would be asserting
something the source explicitly denies. Callers with a short hole want
`Idelchik1966Thick`, which is built for one.

**D2. l/d is HELD at 10 above the range, not extrapolated.** Eq. (7) is
linear in l/d and, extrapolated, walks Cd to zero near l/d = 97.

**D3. Re is floored at 1, where the value holds and the derivative is
reported as ZERO.** This differs deliberately from the Idelchik floor, where
the derivative is *continued*: there the floor is the edge of a table and the
curve is still moving, so a Newton step below it must still see sensitivity.
Here the floor exists only to stop `20/Re` diverging, and the value genuinely
stops changing, so a zero derivative is the truthful report. The source plots
Eq. (12) down to Re ~ 1 in Fig. 10, so the floor sits below the drawn curve.

**D4. No crossflow term, and `dCd/d(U1_over_Vi)` is EXACTLY zero.** The
experiments are plenum-fed. Same treatment and same reason as Idelchik.

**D5. Analytic derivative, no dual numbers needed.** Unlike the Idelchik
members this is closed form, so `dCd/dRe` is differentiated by hand:

    d(1/Cd)/dRe = -b/Re^2 + 15 c u / (den^2 Re ln10)
    dCd/dRe     = -Cd^2 d(1/Cd)/dRe

with `b = 20(1 + 2.25 l/d)`, `c = 0.005 l/d`, `u = log10(0.00015 Re)`,
`den = 1 + 7.5 u^2`. Matches central differences to ~1e-9 relative.

## Verification against the source's own figure

Eq. (12) at l/d = 2, against the curve Fig. 10 plots:

| Re | Cd |
|---|---|
| 10 | 0.0817 |
| 100 | 0.4284 |
| 1000 | 0.7446 |
| 1e4 | 0.8081 |
| plateau | 0.810 = Eq. (7) |

The figure shows the curve crossing 0.5 near Re ~ 160; the implementation
crosses between 100 and 300. Across l/d = 2 to 10 the high-Re limit lands on
Eq. (7) to within 0.003.

## An artifact recorded rather than smoothed away

**Eq. (12) is not monotone forever.** Its two Reynolds terms pull opposite
ways: `20(1 + 2.25 l/d)/Re` falls without limit, while the log-squared term
peaks at Re = 1/0.00015 = 6667 and decays on both sides. Above Re ~ 8.7e5 the
viscous term is spent while the log term is still decaying, so Cd overshoots
Eq. (7) by **2.2e-4** and settles back.

This was found by a monotonicity test written over the whole sweep, and the
test was then narrowed to the validated range rather than loosened -- the
overshoot is 43x beyond the source's validated ceiling of 2e4 and roughly
100x SMALLER than its own stated accuracy of +/-0.02. It is a property of
extrapolating the fit, not a defect to correct.

## Cross-source: the reason this matters

Andrews plate C geometry (l/d = 1.93), over the Re range its own data spans:

| Re | Lichtarowicz | Idelchik 4-18a | McGreehan-S | M-S error |
|---|---|---|---|---|
| 432 | 0.675 | 0.799 | 0.822 | **+21.7%** |
| 1381 | 0.764 | 0.869 | 0.822 | +7.6% |
| 4298 | 0.799 | 0.868 | 0.822 | +2.9% |
| 8718 | 0.808 | 0.850 | 0.822 | +1.7% |
| 2e4 | 0.809 | 0.828 | 0.813 | +0.5% |

Cd varies 20% over this range. McGreehan-Schotsch is floored at Re = 1e4 and
returns a near-constant, reading +21.7% high at the bottom.

**Idelchik and Lichtarowicz disagree by 10-15% at low Re, and both are
primary.** Recorded, not hidden. The geometries differ -- Idelchik 4-18a is a
hole between two infinite plena, Lichtarowicz a long orifice in a hydraulic
line -- but the Cd conventions are compatible, so this is a genuine
source-to-source disagreement of the same kind as the 13% Han/Lau ribbed-wall
friction gap. Neither is adopted as "the" answer; the caller picks the
geometry their hole matches.

## A guard added alongside

`McGreehanSchotsch1988` now WARNS when handed Re below its `re_min = 1e4`
floor, naming Lichtarowicz and Idelchik as the valid alternatives. It warns
rather than refuses because a network can transit low Re during Newton
iteration and refusing would break convergence; the point is that a frozen Cd
should not be silent.

## Falsification

Eight perturbations, each applied alone and reverted, all RED:

| # | perturbation | result |
|---|---|---|
| P1 | Eq. (7) slope 0.0085 -> 0.0100 | RED |
| P2 | viscous coefficient 20 -> 25 | RED |
| P3 | transition term omitted | RED |
| P4 | l/d not held at 10 | RED |
| P5 | flat branch for 1.5 <= l/d < 2 ignored | RED |
| P6 | refusal below l/d = 1.5 removed | RED |
| P7 | derivative sign flipped | RED |
| P8 | Re floor removed | RED |

## What is not here

- **Cavitation.** The title says non-cavitating; the source's Section on
  cavitation inception is not extracted, and nothing here warns about it.
- **Compressibility.** Incompressible throughout.
- **The approach-pipe area ratio `m`.** Our use is plenum-fed, m -> 0. A hole
  fed by a duct of comparable area would need the full definition.
