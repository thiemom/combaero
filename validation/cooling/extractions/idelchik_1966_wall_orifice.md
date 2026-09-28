# Extraction: Idelchik (1966), orifice in a large wall (Section IV)

**Status: IMPLEMENTED 2026-09-28.** All four edge types ship as
`DischargeCdCorrelation::Idelchik1966*`. No digitisation was needed: every
value below is tabulated in the source and was read off the page rendered at
400 dpi with `pdftoppm -r 400`, not the PDF text layer, which garbles these
tables badly (diagram 4-12's rounded row OCRs as `0.I12.oil06 .2`).

Tracked under #339. Feeds #387 (`EffusionPlateElement`), whose geometry this
is exactly.

## Source, pinned

> Idelchik, I.E. *Handbook of Hydraulic Resistance*, 1st English edition,
> AEC-tr-6630 (1966). Israel Program for Scientific Translations.
> `docs/junction/Idelchik.pdf` (gitignored, copyrighted) -- the same copy the
> junction work digitised diagrams 7-1..7-7 from.

Section IV, Figure 4-8 enumerates the four edge types this implements:
(a) sharp-edged, (b) thick-walled, (c) beveled, (d) rounded.

| diagram | page (PDF) | covers | form |
|---|---|---|---|
| 4-17 | 149 | sharp-edged hole in a large wall | `zeta = 2.85` (Re >= 1e5); `zeta = zeta_phi0(Re) + eps_re(Re)` below |
| 4-18a | 150 | thick-walled ("deep") hole | `zeta = zeta'(l/Dh) + lam l/Dh` |
| 4-18b | 150 | beveled edges, 40-60 deg | `zeta = f(l/Dh)` |
| 4-18c | 150 | rounded edges | `zeta = f(r/Dh)` |

## Why THIS geometry and not the in-a-pipe diagrams

Idelchik gives three distinct orifice geometries and they are not
interchangeable:

- **4-9..4-12** -- passage from a conduit of one size to another (`F1 != F2`);
- **4-13..4-16** -- orifice in a straight conduit (`F1 = F2`, the metering
  plate case);
- **4-17/4-18** -- orifice in a large wall (`F1 = F2 = infinity`).

Only the last has no pipe at all, which is why it maps onto
`DischargeHoleGeometry` (which carries no `D` and forms no `beta`) and why it
is the effusion-plate case.

## The definitional point, which is the trap here

`zeta = DH / (rho w0^2 / 2)`, referenced to the **hole** velocity `w0`, and
`DH` is the **full permanent loss** -- there is no downstream recovery
between two infinite volumes. Therefore

    Cd = 1 / sqrt(zeta)

is **exact**, not a convention-dependent conversion.

This is NOT true of the in-a-pipe diagrams (4-13..4-16): there `zeta` is
referenced to the pipe velocity `w1`, and `DH` is still the permanent loss,
whereas ISO 5167's `Cd` is referenced to the **tapping differential**, which
is larger because of downstream pressure recovery. Converting one to the other
needs the recovery relation and is not algebra. This is the #389
"definitional comparability" class -- the same trap as Rohde's total- vs
static-referenced `Cd` and its factor of 1.33.

It is also the second reason the metering and discharge families cannot share
a selector (the first being that a wall hole has no `beta`).

## What this replaced, and why

`Cd_thick_plate`, `Cd_rounded_entry`, `Cd_orifice`, `orifice::thickness_correction`
and `orifice::Cd_rounded` were removed. They computed an ISO 5167
Reader-Harris/Gallagher `Cd` and multiplied it by a correction factor. Two
things were wrong with that:

1. **The base does not apply.** A thick-edged or rounded orifice is not the
   normed device ISO 5167 describes, so its `Cd` correlation is not a valid
   starting point. Idelchik gives the whole `zeta` directly and needs no base.
2. **The correction was an unlabelled fit that had lost its tail.** The
   rounded `K` was `0.5 (1 - (r/d)/0.15)^2` with a plateau at `0.04`. Measured
   against Idelchik diagram 4-12b, which it was evidently fitted to:

   | r/Dh | Idelchik | removed code | error |
   |---|---|---|---|
   | 0.00 | 0.50 | 0.500 | +0.0% |
   | 0.02 | 0.37 | 0.376 | +1.5% |
   | 0.06 | 0.19 | 0.180 | -5.3% |
   | 0.08 | 0.15 | 0.109 | **-27.4%** |
   | 0.12 | 0.09 | 0.020 | **-77.8%** |
   | 0.16 | 0.06 | 0.040 | -33.3% |
   | 0.20 | 0.03 | 0.040 | +33.3% |

   Good to `r/Dh ~ 0.06`, then it falls apart. The thick-plate factor had the
   same shape of problem: `0.35 (1 - e^{-8 t/d})(1 - beta^2) - 5.67 f t/d`,
   clamped to `[0.5, 1.3]`, where `5.67` was commented "friction loss
   calibration factor". Idelchik's own form is `tau`-term plus `lam l/Dh` --
   structurally the same two components, with a real curve and a real friction
   factor in place of both fudges.

Two defects went with them, independent of provenance:

- **An 11.1% jump discontinuity in `Cd` at `Re_D = 1e5`** (0.881391 ->
  0.979326), from an invented `1 - 0.1 (1e5/Re)^0.2` low-Re factor applied
  only below the threshold. A hard C0 break in a solver input. Idelchik has a
  real low-Re branch there instead, and 1e5 is his own regime boundary --
  almost certainly where the threshold came from.
- **`Cd_rounded_entry` at `r/d = 0` silently returned `Cd_Stolz`.** The same
  silent-substitution class as the `Bohl -> Idelchik` alias.

None of the 16 tests over that code would have caught either: all were
directional (`Cd_thick > Cd_thin`, `approx(1.0)` for a thin plate). The
replacements assert the source's tabulated values at every knot.

## Decisions

**D1. One enum member per diagram; the geometry does not auto-select.** The
removed `Cd_orifice` picked a correlation from `r` and `t` behind the
caller's back. A caller who asks for rounded and supplies `r = 0` now gets the
rounded correlation evaluated at `r/d = 0`, which is continuous with sharp
(below), rather than a different correlation.

**D2. All four tables anchor on `zeta = 2.85`.** 4-18a and 4-18b tabulate it
at `l/Dh = 0` directly. 4-18c's tabulated row starts at `r/Dh = 0.01`, but its
graph c starts at 2.85 and the sharp value forces it. So every edge type
degrades continuously to sharp, and the selector is safe to sweep. Pinned by
`EveryEdgeTypeDegradesToSharpAtZero`.

**D3. Monotone cubic (Fritsch-Carlson) interpolation, not linear, not a
natural spline.** Linear is C0: `dCd/dRe` would jump at each of 14 knots --
the hazard class we spent #383 regularising out of the McGreehan chain, and 14
of them is worse than the one removed. A natural cubic overshoots on the flat
tail (`zeta'` is 1.58, 1.55, 1.55 over its last three knots) and would invent
a `Cd` above the source's. Fritsch-Carlson cannot overshoot, so the
interpolant stays inside the tabulated envelope by construction. Pinned by
`MonotoneInterpolantNeverLeavesTheTabulatedEnvelope`.

The derivative is the Hermite form differentiated analytically -- exact, not a
difference, per the `(f, J)` rule. Pinned against central differences by
`AnalyticDerivativeAgreesWithFiniteDifferences`.

**D4. `lam` comes from `friction.h`, and that is reproduction, not
substitution.** Idelchik takes `lam` from his diagrams 2-2..2-5 as
`f(Re, Delta/Dh)`. Those diagrams ARE the Colebrook/Nikuradse family, so
Haaland's explicit form computes the same quantity. A drilled or laser-cut
cooling hole is hydraulically smooth at these Reynolds numbers
(`default_roughness_over_d = 0`). The term is worth about 2% of `zeta` at
`l/Dh = 2`.

Below `Re = 2300` the laminar branch `lam = 64/Re` is used -- Hagen-Poiseuille,
exact, not a correlation -- because at `Re = 25` with `l/Dh = 2` bore friction
dominates `zeta'` entirely and omitting it would be wrong by more than a
factor of three. A smoothstep in `log(Re)` carries one branch to the other
over `2300 < Re < 4000`. **That blend is a numerical device, not physics**:
neither branch is the source's inside the window, and Idelchik has no
transition correlation either. It moves `zeta` by at most 0.9% at `l/Dh = 2`
and nothing outside the window.

**D5. No crossflow term, and `dCd/d(U1_over_Vi)` is returned as EXACTLY
zero.** Idelchik's geometry is plenum to plenum; there is no approach velocity
to form `U1/Vi` from. The absence is a property of the source, not a
truncation, and the test says so. A hole that does see inlet crossflow wants
`McGreehanSchotsch1988`, whose Eq. (17) is precisely the term Idelchik lacks.

**D6. Held at the table edges.** Below `Re = 25` the value is held and the
derivative continued -- the same treatment, and the same reason, as
`mcgreehan_schotsch::re_min`: a Newton step below the floor must still see a
Re sensitivity or it stalls. Above the last knot every table is already flat,
so holding is the source's own behaviour rather than an extrapolation.

**D7. The printed `0.342` is a rounded `1/2.85`, and we use the exact value.**
Diagram 4-18a's low-Re form is
`zeta = zeta_phi0 + k eps_re zeta' + lam l/Dh`, printed with `k = 0.342`. It
must reduce to item 1's `zeta' + lam l/Dh` as `eps_re -> 2.85` and
`zeta_phi0 -> 0`, which forces `k = 1/zeta_sharp = 0.350877` exactly. The
printed three-figure rounding IS the entire 2.5% by which the source's own
two formulas disagree at high Re.

Using the printed value would leave a **2.5% step in `Cd`** at the top of the
table -- a C0 break in a solver input, and the same defect class as the 11%
jump this correlation replaced. `thick_low_re_coef_as_printed` keeps the
printed figure in the header so the discrepancy is on the record rather than
silently corrected.

**D7b. Neither branch is a separate regime; the table runs everywhere.**
Diagram 4-17 item 1 reads "Re >= 1e5: zeta = 2.85". That is the coarse
statement of the same curve item 2 gives finely: the table reaches 2.85 at
its LAST knot, `Re = 1e6`, and says **2.60** at `Re = 1e5`. An early
implementation treated the two as separate branches and put a 4.5% step at
1e5 -- caught by `test_cd_is_continuous_across_the_1e5_boundary`, which was
written to guard against exactly the defect being removed and instead caught
a fresh instance of it. Both correlations now evaluate one expression over
the whole range, with PCHIP holding the endpoint above the last knot.

One consequence worth stating: **`Cd` is not monotone in `Re`.** `zeta`
bottoms out near `Re = 400` (1.91) and rises in both directions, so `Cd`
peaks there at 0.72 and falls to 0.58 at `Re = 25` and 0.59 at `Re = 1e6`.
That is Idelchik's own table, not an artifact.

**D8. Bohl was deleted rather than kept as a refusing member.** The criterion
is whether a source plugs a coverage hole. Bohl/Elmendorf, *Technische
Stroemungslehre*, covers DIN 1952 normed orifices (the ISO 5167 lineage) and
standard loss coefficients (the Idelchik class). Both rows are already held by
primary sources, so it plugs nothing. Before 2026-09-28 `BohlThick` and
`BohlRounded` fell through to the Idelchik implementations, silently returning
a different correlation than the one asked for.

## Cross-source accuracy (LABELLED cross-source, per the validation policy)

Idelchik (1966, Russian handbook, `zeta` measurements) against McGreehan &
Schotsch (1988, gas-turbine cooling paper, `Cd` measurements). Independent
sources, independent decades, independent apparatus. At `Re = 1e5`, `r/d = 0`,
no crossflow:

| l/Dh | Cd Idelchik | Cd M-S | diff |
|---|---|---|---|
| 0.0 | 0.5923 | 0.5926 | **+0.0%** |
| 0.2 | 0.6063 | 0.6031 | -0.5% |
| 0.4 | 0.6202 | 0.6379 | +2.9% |
| 0.6 | 0.6537 | 0.6849 | +4.8% |
| 0.8 | 0.7161 | 0.7305 | +2.0% |
| 1.0 | 0.7538 | 0.7660 | +1.6% |
| 2.0 | 0.8032 | 0.8055 | +0.3% |
| 4.0 | 0.8032 | 0.7888 | -1.8% |

And on radius, at `l/Dh = 0`: within 1.1% to 1.9% over `r/Dh = 0.01..0.12`,
worst +3.4% at 0.20.

**This is an accuracy result, not a fidelity one.** It says the two sources
agree; it does not say either implementation mirrors its own paper. Fidelity
for Idelchik is the knot-equality tests against his tabulated `zeta`; fidelity
for McGreehan-Schotsch is its own extraction record.

The agreement at `l/Dh = 0` is the striking one -- 0.5923 against 0.5926, when
McGreehan's Eq. (8) asymptote `0.5885` and Idelchik's `1/sqrt(2.85)` were
derived 22 years and one iron curtain apart. Pinned by
`AgreesWithMcGreehanSchotschAcrossTheSharedRange`, which would notice a
regression on either side.

## What is NOT covered

- **The in-a-pipe family (4-13..4-16).** Diagram 4-16 carries a full 8x21
  `zeta(r/Dh, F0/F1)` table, but 4-14's thick-edged branch needs `tau(l/Dh)`
  off diagram 4-11 (p.143), which is a **graph, not a table** and is the one
  piece here that needs digitising. Not started.
- **Low-Re for beveled and rounded.** Diagram 4-18 gives the low-Re branch
  only for sharp (4-17) and thick (4-18a). The beveled and rounded curves are
  stated for `Re >= 1e5`; below that the value is held rather than borrowed
  from the sharp branch.
- **Lichtarowicz (1965)**, long orifices `L/d = 2..10` including cavitation
  limits. PDF is in `docs/orifices/`; declared as
  `DischargeCdCorrelation::Lichtarowicz1965`, which **refuses explicitly**
  rather than substituting.
- **Compressibility.** All of Section IV is incompressible. Diagram 4-19
  carries a Mach correction that has not been extracted.

## Falsification

Seven perturbations, each built and run against the suite, then reverted.

| # | perturbation | result |
|---|---|---|
| P1 | last rounded table value `1.37 -> 1.30` | RED in Python only -- see below |
| P2 | revert `k` to the printed `0.342` | RED: cross-source + degrades-to-sharp |
| P3 | PCHIP tangent scaled `4x` | RED: envelope test |
| P4 | `dCd/dRe` off by 2% | RED: analytic-vs-FD |
| P5 | crossflow derivative `0.0 -> 1e-30` | RED: exact-zero test |
| P6 | laminar friction branch removed | **NOTHING WENT RED** -- gap, now closed |
| P7 | restore the hard `Re >= 1e5` switch | RED: continuity + low-Re knots |

**P1 is the division-of-labour finding.** The C++ knot tests read the table
from the header they are testing, so they pin *interpolation exactness* and
are tautological with respect to the DATA. The Python tests transcribe the
numbers independently from the page, so they are what guards the data. Both
are needed and neither substitutes for the other. P1 is red in Python, green
in C++, and that is correct rather than a weakness -- but it is only correct
because the Python copy exists.

**P3 needed `4x`, not `2.5x`.** The first attempt scaled the tangent by 2.5
and nothing went red, which looked like a missing test. It is not: the
classic monotone-Hermite bound is `m/delta <= 3`, so `2.5x` is still inside
the no-overshoot region. `4x` crosses it and the envelope test fires. Worth
recording so the next reader does not re-derive it.

**P6 found a real gap.** Swapping the laminar branch (`64/Re`) for the
turbulent one changed no test result. It matters: at `Re = 25` with
`l/d = 2`, `lam = 2.56` and the friction term `lam l/Dh = 5.12` is larger
than `zeta'` itself (1.55), so dropping the laminar branch would overstate
`Cd` by more than a factor of two. A deep hole at creeping flow is a
Poiseuille pipe, not an orifice. Closed by
`LaminarBoreFrictionDominatesADeepHoleAtLowRe` and
`FrictionBlendIsContinuousAcrossTheTransitionWindow`; P6 now goes red on
both.

A separate bug surfaced while probing P3: `DischargeHoleGeometry::bevel` was
never bound to pybind11, so the beveled correlation was unreachable from
Python and `OrificeElement`'s `'IdelchikBeveled'` arm would have raised
`AttributeError` on first use. Bound, with a `units_data.h` entry.
