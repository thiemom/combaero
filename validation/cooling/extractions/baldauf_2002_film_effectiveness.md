# Extraction: Baldauf et al. (2002), film-cooling effectiveness

**Status: IMPLEMENTED 2026-09-29** as
`combaero::cooling::film_effectiveness_baldauf_2002`. Forty equations, all
transcribed exactly; one of them contradicts the paper's own worked example
and is implemented as printed. See "The one discrepancy" below.

## Source, pinned

> Baldauf, S., Scheurlen, M., Schulz, A. and Wittig, S. (2002). "Correlation
> of Film-Cooling Effectiveness From Thermographic Measurements at Enginelike
> Conditions." *ASME J. Turbomachinery* **124**(4), 686-698.
> DOI 10.1115/1.1504443. Conference version ASME GT-2002-30180.
> `docs/heat_transfer/film/686_1_Baldauf_Film_2002.pdf` and the publisher HTML
> alongside it (both gitignored, copyrighted).

A row of CYLINDRICAL, streamwise-inclined holes on a flat plate, measured by
infrared thermography and laterally averaged.

## Why this source

The review by Xia, Chen and Ellis (2024, *Energies* **17**, 4480) is blunt
about the alternatives: single-hole models "are completely wrong in the near
wake region", and the older row correlations "give an unrealistic maximum at
the ejection position and enforce a monotonic decay", so "the vicinity of the
ejection is usually excluded".

Baldauf is built to avoid both. It is valid "from the point of the ejection
to far downstream", and it carries the ADJACENT JET INTERACTION -- the
lateral hole-spacing effect that drives jet lift-off -- as a correlated
parameter. Colban, Thole and Bogard (2011) benchmark it and find it "works
well for cylindrical hole injection".

## How it was extracted -- and why that matters here

The PDF text layer is unusable for these equations (Eq. 8 OCRs as
`h¯ ~sin a!0.06 s D P0.9/ s D`). **The publisher HTML carries MathML**, so all
40 equations came out exactly rather than by reading images.

Cross-checked four ways, which is what makes the one discrepancy credible:

1. **HTML MathML** -- the publisher's own markup.
2. **The PDF at 600 dpi**, read visually. Eqs. 8-15 agree with the MathML
   character for character.
3. **The PDF text layer**, garbled but structurally consistent.
4. **A second reader's independent transcription** of Eq. 31.

## Angle units -- measured, not assumed

Table 3 reports `alpha [deg]`, but every trigonometric function consumes
RADIANS. Established from the worked example rather than assumed: seven
coefficients spanning `cos(a)`, `cos(1.5a)`, `cos(2.3a)`, `cos(2.5a)`,
`sin(2a)` and `cos^0.65(a)` reproduce Table 4 on radians and fail on degrees
by 10-30%.

**One of them would NOT have caught it.** `b*_T` comes out 0.70138 on radians
against 0.70178 on degrees -- 0.06% -- because its `cos 2.5a` sits inside an
`exp[...]` heavily damped at `Tu = 1.5%`. A units test written against `b*_T`
alone would have passed a units bug. The test asserts on `xi_c`, `k` and
`a_1`, which move 10-30%.

The API therefore takes DEGREES (matching how geometry is specified, and
`EffusionPlateElement`) and converts internally.

## Structure

    x/D    -> xi'      Eq. (40), inverted in CLOSED FORM -- no iteration
    xi'    -> eta*'    Eq. (36), the turbulence-dependent base curve
    eta*'  -> eta*     Eq. (37), undoing the adjacent-jet downstream branch
    eta*   -> eta_c    Eq. (38), blowing up to this case's peak
    eta_c  -> eta_bar  Eq. (39), backscaling out the geometry normalisation

Table 2's base-curve constants (`xi_0 = 9`, `eta_0 = 5.8`, `a* = 4`,
`b* = 0.7`, `c* = 0.24`) are valid for all geometry and density ratios --
that is the point of collapsing every measurement onto one base curve.

**A free consistency check on Eq. (17):** its velocity-ratio exponent
`((s/D)/3)^-0.75` reproduces the paper's own stated values of 1.0 at
`s/D = 3`, 1.37 at 2 and 0.68 at 5, to within 0.015.

## The one discrepancy

Running the transcription on Table 3's inputs reproduces **17 of 19** Table 4
coefficients to better than 2e-6 relative. Eq. (31) does not, and Eq. (32)
inherits it:

    b_0   Eq. (31) as printed -> 0.83612467    Table 4 -> 0.61626073
    b_1   Eq. (32) follows    -> 0.74322193    Table 4 -> 0.54778731

Everything ruled out before concluding it is the paper's:

| hypothesis | result |
|---|---|
| misread equation | four independent channels agree on the printed form |
| radians vs degrees, conversion inside or outside | bit-identical, 0.83612467 |
| whole argument taken as radians, `sin(28 rad)` | 0.767 |
| nine angle conventions (sign, from-normal, supplement, ...) | 0.388 to 1.017 |
| ANY alpha in +/-360 deg | only -10.31, -182.52, +203.15 deg |
| single-token typo in each of the six constants | none lands plausibly |
| a stray `s/D` | would need 5.87, but 17 other values used 3 correctly |

Table 4 is self-consistent with itself: `b_1/b_0 = 0.88889 = 1/(1+M^-3)`
exactly, per Eq. (32). So its pair came from some correct implementation --
just not the printed Eq. (31). Not enough information survives to recover
what the spreadsheet actually did.

**Decision: implement Eq. (31) AS PRINTED.** Fidelity means implementing what
the paper states, not what it might have meant. `BaldaufTable4.
EquationThirtyOneContradictsThePapersOwnTable` pins both numbers so a future
reader meets the contradiction immediately, and so anyone "fixing" the
constant must delete an assertion explaining why not to.

### How much it matters

`b_1` enters `(1 + (xi'/xi_1)^(b_1 c_1))^(1/c_1)`, so it bites once `xi'`
passes `xi_1 = 65/(M/2.5)^a_1` -- at HIGH blowing rate:

| s/D | M | alpha | max relative difference in eta_bar |
|---|---|---|---|
| 2 | 0.5 | 30 | 1.0% |
| 2 | 0.5 | 90 | 5.2% |
| 2 | 2.5 | 30 | 28.0% |
| 5 | 2.5 | 90 | **49.6%** |

Under 5% below `M ~ 0.5`. Andrei's effusion data is BR 1-3, so this sits
inside the range we care about, and the Andrei/Murray scoring will be run on
both branches and REPORTED -- not used to silently pick one.

## Validation

**End to end against Fig. 14**, the paper's own measurement-versus-correlation
plot for exactly the Table 3 case. The implementation rises from ~0 at the
ejection point, peaks near `x/D = 30`, and settles at `eta = 0.117-0.120`
against the figure's plotted plateau of roughly 0.11-0.12. This validates the
full 40-equation chain and the assembly order, and it holds for either `b_0`
branch (the two give 0.118 and 0.115 here, which is why Fig. 14 cannot
arbitrate).

**Magnitude and ordering against Fig. 8**, which plots peak effectiveness
against `M` per hole spacing with bands of roughly 0.35-0.6 at `s/D = 2`,
0.15-0.45 at 3 and 0.05-0.3 at 5. Peaks land inside, and closer spacing
always cools better.

**Analytic derivatives** `(eta, d eta/dM, d eta/dP)` via forward-mode dual
numbers over the SAME templated equation chain the value uses, so the two
cannot drift apart. Agreement with central differences: 9e-9 relative worst
over 18 operating points. `d eta/dM` changes sign through the lift-off peak,
which is physical and is what a solver needs.

Paper's own accuracy: RMS deviation 5.5%, 5% at the apex, 3% on the
descending branch.

## Falsification

Nine perturbations, each applied alone and reverted.

| # | perturbation | result |
|---|---|---|
| P1 | Table 2 `a*` 4 -> 3 | RED |
| P2 | Eq. (9) cos term dropped | RED |
| P3 | Eq. (15) 0.048 -> 0.060 | RED (after tightening, see below) |
| P4 | alpha fed as degrees | RED |
| P5 | Eq. (40) exponent forced to 1 | RED (after tightening) |
| P6 | Eq. (39) density backscale dropped | RED (after tightening) |
| P7 | Eq. (34) turbulence term neutralised | RED |
| P8 | Eq. (32) collapsed to `b_1 = b_0` | RED (after tightening) |
| P9 | Eq. (31) replaced by Table 4's value | RED |
| P10 | Eq. (13) 336 -> 300 | RED (after tightening) |
| P11 | Eq. (18) `xi_c` forced to 1 | RED (after tightening) |

**The first run found three real gaps.** P3, P6 and P8 originally changed no
test result: the figure-based checks are necessarily loose (Fig. 14 is
readable to about +/-0.005, Fig. 8's bands span all nine alpha/P curves), so
a single coefficient could move `eta` several per cent unnoticed. Closed by
adding `PinnedValuesAcrossTheEnvelope`, eight pinned values at 1e-9 across
the parameter space.

Those pins are REGRESSION pins and are labelled as such in the test -- they
are this implementation's own output, not numbers the paper prints. Their
literature anchor is the Fig. 14 and Fig. 8 tests: those say the curve is in
the right place, the pins say it has not moved.

## What is not here

- **Only one row.** Multi-row superposition (Sellers, and Gao et al. 2025's
  mainstream-temperature correction for the error that accumulates with row
  count) is the next piece, and it is what #386 and #387 both need.
- **Cylindrical holes only.** Shaped holes are Colban, Thole and Bogard
  (2011), in `docs/heat_transfer/film/`, not extracted.
- **No heat transfer coefficient.** The companion paper (ref [25]) gives the
  matching `h` correlation and is not in hand.
- **`delta_1/D` and `L/D`** appear in Eq. (2)'s parameter list but are not
  correlated parameters in the final form.
- **The envelope is recorded, not enforced.** A network solve transits odd
  states during Newton iteration; refusing there would break convergence
  rather than protect anyone. `baldauf2002::M_min` and friends are there for
  a caller who wants to check.
