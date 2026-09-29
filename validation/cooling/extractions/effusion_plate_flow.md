# Extraction: effusion plate, flow side (#387)

**Status: IMPLEMENTED 2026-09-29, FLOW SIDE ONLY.** `EffusionPlateElement`
ships the geometry and discharge path. The thermal model is deliberately not
in it -- see "What is not here" below.

Values read off pages rendered at 400 dpi with `pdftoppm -r 400`, not the PDF
text layer.

## Sources, pinned

> Andrews, G.E., Asere, A.A., Hussain, C.I., Mkpadi, M.C. and Nazari, A.
> (1988). Impingement/Effusion Cooling: Overall Wall Heat Transfer. ASME
> 88-GT-290, Department of Fuel and Energy, University of Leeds.
> `docs/heat_transfer/film/Impingement_Effusion_Cooling_Overall_Wal.pdf`
> (gitignored, copyrighted).

> van de Noort, M. and Ireland, P. (2022). A Low Order Flow Network Model for
> Double-Wall Effusion Cooling Systems. *Int. J. Turbomach. Propuls. Power*
> **7**, 5. DOI 10.3390/ijtpp7010005. Open access (CC BY-NC-ND).
> `docs/heat_transfer/film/A_Low_Order_Flow_Network_Model_for_Doubl.pdf`.

> Andrei, L., Andreini, A., Bianchini, C., Caciolli, G., Facchini, B.,
> Mazzei, L., Picchi, A. and Turrini, F. (2014). Effusion cooling plates for
> combustor liners: experimental and numerical investigations on the effect
> of density ratio. *Energy Procedia* **45**, 1402-1411.
> `docs/heat_transfer/film/1-s2.0-S1876610214001489-main.pdf`. **Not used
> yet** -- it is the PR 2 film closure.

## Geometry, validated against Andrews Table 1

| Plate | | N [m^-2] | D [mm] | X [mm] | X/D | A/A_h |
|---|---|---|---|---|---|---|
| Impingement | A | 4306 | 1.38 | 15.2 | 11.0 | 8.38 |
| Effusion I | B | 4306 | 2.16 | 15.2 | 7.1 | 5.28 |
| Effusion II | C | 4306 | 3.27 | 15.2 | 4.7 | 3.4 |

Square array, plate thickness 6.3 mm, 8 mm impingement gap.

**X is 0.6 inch, not 15.2 mm.** `1/(15.2e-3)^2 = 4328`, which is not the
tabulated 4306; `1/(15.24e-3)^2 = 4305.8` is. The table's `15.2` is a rounded
print of an imperial pitch. This matters because `N = 1/X^2` is the relation
the element uses to derive its hole count, and checking it against `15.2`
would have shown a spurious 0.5% error.

Our derivation reproduces the table: `N = 4300 m^-2` on a 0.1 x 0.1 m panel
(43 holes, exact 43.056) and `A/A_h = 5.35` against the tabulated 5.28, a
1.3% gap consistent with `D` and `t` being nominal values.

`A` is the hole approach area `X^2 - pi D^2/4` and `A_h` the hole internal
area `pi D L`, per the paper's own nomenclature.

## Decisions

**D1. Geometry is given as it is designed, not as an area.** Pitch, hole
diameter, wall thickness and inclination; the hole count follows. The count
is ROUNDED to a whole number, and `hole_count_exact` plus `pitch_actual`
record what the rounding did -- the element flows whole holes, and the
implied pitch is not quite the one asked for.

**D2. An inclined hole is longer than the wall is thick**, `L = t/sin(alpha)`.
That is the point of inclining them: more internal surface, and a longer
throat, for the same wall. At 30 deg a 1 mm wall gives a 2 mm hole.

**D3. The correlation is fed ONE hole, the flow equation the summed area.**
Conflating them would hand a 2.16 mm hole's correlation a ~100 mm equivalent
bore. Pinned by `test_correlation_sees_one_hole_not_the_equivalent_bore`.

**D4. Default `IdelchikThick`; ISO 5167 correlations REFUSED.** Diagram 4-18a
is a thick-walled hole in a large wall between two plena -- exactly a
plenum-fed effusion plate, and valid to Re = 25. A panel has no pipe and no
`beta`, so a normed metering correlation is undefined here, not merely
inaccurate. See [[project_idelchik_orifice_coverage]].

**D5. No velocity-of-approach `beta`, and the reason is physical.** van de
Noort's Equation (8) applies `1/sqrt(1 - beta^4)` with `beta` built on a
half-pitch square inlet area. We do not, because a PLENUM has no approach
velocity to correct for. A CHANNEL-fed panel does have an approach flow, but
that is a **crossflow** -- McGreehan-Schotsch's `U1/Vi` -- not a
velocity-of-approach term. This is why both correlations are carried, and it
is the same definitional-comparability class as #389. At van de Noort's own
pitch the factor is 1.0032; at `X/d = 3` it would be 11.6%, so the distinction
is not academic at tight pitch.

**D6. One panel homogenises, and the answer is more panels.** van de Noort
and Ireland show the effect is large: with uniform inlet AND outlet pressure,
holes at the end of an array pass ~75% of what a central hole passes; under a
spanwise gradient one row took ~10% of the total against ~20% for others.
Their LOM is structurally a node-pressure network of discharge links, which
is what combaero's solver already is, so the decomposition is theirs and is
validated (10.9% on total mass flow, ~20% on individual channel shares).

**D7. The effusion CHANNEL is a composition, not a new element.** A panel
hung off each segment of a coolant channel gives
`m_channel_in = m_channel_out + m_effusion` at every node -- coolant flow as
`f(x)`. That is van de Noort's Level 4 with lateral links. No channel-specific
element was written and none is needed.

**D8. Ingestion is reported, not averaged away.** A reversed drive means hot
gas entering (van de Noort's `CMF > 0.5`). A homogenised panel would fold it
into a healthy net outflow, so `is_ingesting` is a diagnostic.

## Measured behaviour worth recording

**A uniform gas side does NOT give a uniform wall.** Channel friction plus
the mass bled off drops the coolant pressure along the wall. Four 20 mm
segments, feed at 1.20 bar, gas at 1.00 bar throughout:

| panel | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| node Pt [kPa] | 114.2 | 109.8 | 106.5 | 104.0 |
| drive [kPa] | 14.2 | 9.8 | 6.5 | 4.0 |
| share of bleed | 0.335 | 0.274 | 0.220 | 0.172 |

The first panel passes ~1.9x the last. This is the coolant-side counterpart
of van de Noort's migration and is invisible to a single panel.

**A falling external pressure EVENS IT OUT.** With gas at 1.03, 1.01, 0.99,
0.97 bar the shares become 0.294, 0.256, 0.231, 0.219 -- spread 0.075 against
0.163. The two gradients partly cancel, which is the opposite of the
intuition that an external gradient always worsens maldistribution, and it is
a design lever: matching the external gradient to the channel loss is what
makes the bleed uniform. Both pinned by tests.

## A bug this work surfaced

`OrificeElement._effective_Cd` passed the PIPE Reynolds number `Re_D` to the
discharge-hole correlations, which are functions of the HOLE Reynolds number.
A regression from #409, one day old. Wrong by 1-8% for an orifice in a pipe
and by up to **38% for a plenum-fed hole**, where `D_up = 0` froze `Re_D` at
a `1e5` fallback so `Cd` stopped responding to flow at all -- which is
precisely the effusion case, so it would have silently poisoned this element.
Fixed via an overridable `_hole_reynolds` hook.

## Falsification

Seven perturbations, each applied alone and reverted:

| # | perturbation | result |
|---|---|---|
| P1 | hole count truncated, not rounded | RED |
| P2 | hole length ignores the angle | RED |
| P3 | correlation fed the equivalent bore | RED |
| P4 | hole Re not divided by the hole count | RED |
| P5 | ingestion never flagged | RED |
| P6 | normed correlations allowed | RED |
| P7 | flow area is the panel, not the holes | RED |

P1 initially reddened only the "rounds to none" test: the 43.056-hole case
does not distinguish rounding from truncation, so on its own it would have
let a truncating implementation pass. A 43.7-hole case was added and P1 now
reds the rounding test directly.

**The first falsification run was invalid and was redone.** The restore path
`/tmp/comp.bak` was not writable under the sandbox, so `cp` failed silently
and the seven perturbations ACCUMULATED -- every result after the first was
measured against a compound perturbation. Redone one at a time with the
backup inside the project tree and `__pycache__` cleared between runs. Same
trap as [[feedback_extraction_tooling_traps]], different path.

## What is not here

**No thermal model.** Effusion cools in three places -- the coolant-side hole
approach, the throat, and the external film -- and Andrews shows they are not
additive in effectiveness:

> "The film cooling reduces the mean gas temperature adjacent to the wall.
> This leads to the temperature difference between the wall and the coolant
> being reduced which in turn reduces the heat flux removed by the internal
> wall cooling."

So an OVERALL effectiveness correlation cannot be the closure: it already
contains the internal convection a resistance network would compute again.
PR 2 takes the split explicitly -- internal convection from the resistance
network, external film from Andrei et al. (2014), which is ADIABATIC by
construction (PSP mass-transfer analogy) and therefore cannot contain
internal convection. Overall `eta` becomes an OUTPUT, scored against Andrews
Fig. 10's effusion-only curves.

Also absent: Andrews Fig. 10 is not digitised; the low-Re branch for beveled
and rounded holes; compressibility (all of this is incompressible).
