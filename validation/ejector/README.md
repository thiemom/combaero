# Ejector validation (Huang 1999 + Kracik & Dvorak 2016)

Groundwork and reference model for a supersonic-ejector network element, plus a
language-neutral golden dataset that the eventual C++ implementation (analytic
Jacobian + ctests), the pybind layer, and the GUI all validate against **without
recreating a CoolProp environment**.

## Source

- `../../docs/ejector/Huang_1d_analysis_ejector.pdf` — B.J. Huang, J.M. Chang,
  C.P. Wang, V.A. Petrenko, "A 1-D analysis of ejector performance",
  *Int. J. Refrigeration* 22 (1999) 354-364. Source of `entrainment_ratio`
  (Eqs. 1-8) and of the digitized Tables 3-4 validation data.
- `../../docs/ejector/Development_of_an_Analytical_Method_for.pdf` — J. Kracik,
  V. Dvorak, "Development of an Analytical Method for Predicting Flow in a
  Supersonic Air Ejector", EPJ Web of Conferences 114, 02059 (2016),
  DOI: 10.1051/epjconf/201611402059. Source of the mixing closure
  `critical_back_pressure` implements (their Eqs. 7-13) -- see Accuracy below
  for why this replaced Huang's original P_c* chain.
- `../../docs/ejector/art_10.1007_s40430-025-05536-7.pdf` — S. Akbarnejad,
  M. Ziabasharhagh, "A novel 1D model for the analysis of double-choked
  ejectors validated by CFD simulations", J. Braz. Soc. Mech. Sci. Eng.
  47:253 (2025), DOI: 10.1007/s40430-025-05536-7. CFD confirmation that a
  naive momentum balance across the mixing region is unreliable (evaluated,
  not adopted -- see Accuracy below).
- `../../docs/ejector/Lienhard_One dimensional.pdf` — R.K. McGovern,
  K.V. Bulusu, M.A. Antar, J.H. Lienhard V, "One-dimensional Model of an
  Optimal Ejector and Parametric Study of Ejector Efficiency" (MIT/GWU/KFUPM).
  Independent corroboration of the same momentum-retaining, shock-free mixing
  pattern as Kracik & Dvorak (evaluated, not adopted -- see Accuracy below).
- `../../docs/ejector/Henry_2007_HEFAT_HF1_supersonic_ejector.pdf` — Henry, Leclaire, Hemidi,
  Seynhaeve, Bartosiewicz, HEFAT2007: independent air-ejector CFD + experiment
  (validation cross-check target; not a model source). Its measured Fig. 5 is
  `data/henry2007_fig5.py`, with the rig's geometry.
- `../../docs/ejector/Theoretical Model for Predicting Performance of a Gas
  Ejector.pdf` — S. Akbarnejad, M. Ziabasharhagh, Thermal Engineering 72(1)
  (2025). Steam-ejector CFD (Tables 1-3, `data/akbarnejad2025_tables.py`), a
  cross-check target. Its own model is a scaling law from a reference curve
  of the same ejector, so it gives no geometry-only anchor to adopt.
- I.W. Eames, A. Milazzo, D. Paganini, M. Livi, "The design, manufacture and
  testing of a jet-pump chiller for air conditioning and industrial
  application", *Applied Thermal Engineering* 58 (2013) 234-240,
  DOI: 10.1016/j.applthermaleng.2013.04.028. R245fa chiller, measured (not
  committed). Fig. 8 and Table 1 are `data/eames2013_fig8.py`, a cross-check
  target.

## Model scope

`models/huang1999.py` implements the **critical (double-choked) mode** 1-D
model: primary flow choked at the nozzle throat, entrained flow choked at the
hypothetical throat y-y, constant-pressure mixing (Huang's Eqs. 1-8,
`entrainment_ratio`). In this regime the entrainment ratio omega is fixed and
independent of back pressure -- the flat plateau of the ejector characteristic.
It also carries the off-plateau operating-regime closures (see "Operating
regimes" below): `cd_nozzle_mass_flow` / `cd_nozzle_exit_static` (choked or
unchoked primary), `dead_head_back_pressure` + `subcritical_entrainment_ratio`
+ `blended_entrainment_ratio` (the subcritical droop), and
`jet_pump_entrainment_ratio` (the unchoked subsonic jet pump).

It also implements a critical back pressure P_c* (`critical_back_pressure`)
using Kracik & Dvorak's mixing closure. Since the element's default became
`anchor_closure="huang1999"`, that is the selectable alternative; the default
P_c* is Huang's printed chain (`huang_critical_back_pressure`, see "Cross-source
accuracy of the anchors" below). Kracik & Dvorak's closure (their Eqs. 7-13)
works differently from Huang's original Eqs. 9-18: mass, momentum and energy are solved together directly
from the y-y state to the fully-mixed, subsonic state -- no separate loss
coefficient on momentum, no explicit shock model; the conjugate subsonic root
of the same conservation laws plays the shock's role. A `recovery_efficiency`
parameter (default 1.0, no artificial loss) scales the resulting lossless
mixed stagnation pressure down to P_c*, exposed for a user to tune against
their own data rather than fitted to this project's validation set. Unlike
omega, this needs the specific gas constant R (each stream's own
critical/sonic velocity is referenced to its own stagnation temperature).

`EjectorGeometry` is pure dimensionless area ratios (A_p1/A_t, A_3/A_t), so
Huang's single-centered-nozzle geometry and a peripheral multi-nozzle geometry
(e.g. Kracik & Dvorak 2016's 12-nozzle wind-tunnel ejector) are the same model
-- only the *inputs* differ. `EjectorGeometry.from_count(count, ...)` builds
the aggregate ratios for `count` identical, independently-choking nozzles
sharing one mixing chamber (their own stated assumption): A_p1/A_t comes out
invariant to count, while A_3/A_t shrinks as count grows, since the aggregate
throat area scales with count but the shared mixing chamber does not.
`count=1` reduces exactly to the direct constructor.

## Operating regimes

The device passes through these regimes as the forced primary flow rises: (1)
`mp < mdot_choked(Pb)` -> primary **unchoked**, a subsonic jet pump (primary Pt
floats just above the back pressure); (2) primary **choked** -> supersonic
ejector, whose entrained flow is then **critical** (`Pb <= P_c*`, the flat
plateau), **subcritical** (`P_c* < Pb < P_b0`, drooping omega), or **backflow**
(`Pb >= P_b0`). The critical plateau was the original model; the subcritical
droop and the unchoked jet pump are now implemented too (Phase 1A/1B reference
closures listed under "Model scope"; the network `EjectorElement` composes all
of them into one C1 residual system that resolves every regime in a single
solve -- Phase 2). Full diagnosis (of the original non-converging
low-primary-flow case, which used to converge to a self-contradictory
choked-branch root with primary Pt *below* the back pressure) and the
provenance of the fix are in `OPERATING_REGIMES_DESIGN.md`. Backflow
(`Pb >= P_b0`) remains out of scope.

## gamma provenance (the only CoolProp step)

Huang's model uses ideal-gas constant-gamma relations. gamma is evaluated per
suction condition at the **entrained-flow choking plane** y-y (the plane that
sets the entrainment), from the CoolProp R141b EOS:

    T_sy = T_e / ((gamma + 1) / 2),   gamma = cp0(T_sy) / (cp0(T_sy) - R)

giving ~1.111 (Table 3: 8 C -> 1.11210; Table 4: 12 C -> 1.11103). Evaluating at
the primary *stagnation* temperature instead fits far worse (mean 2.2% vs 0.8%),
because entrainment is governed by the cold choking plane, not the hot inlet.

These values are **baked into** `data/huang1999_tables.py`; CoolProp is only
needed to regenerate them, which is rare.

P_c*'s validation target uses CoolProp the same way: each row's T_c* (the
saturated-vapour temperature the paper prints alongside its theory A_3/A_t and
omega) is converted to a saturation pressure P_c* via the CoolProp EOS and
baked into `PC_STAR_PA_BY_TC_STAR_C`.

## Files

| File | What |
|------|------|
| `models/huang1999.py` | Reference model: `entrainment_ratio` + `critical_back_pressure`. |
| `data/huang1999_tables.py` | Digitized Tables 3 & 4 (theory column) + baked gamma and P_c* targets. |
| `data/generate_gamma.py` | Regenerate gamma from CoolProp (**needs CoolProp**). |
| `data/generate_reference.py` | Emit the golden JSON + C++ header (no CoolProp). |
| `data/huang_pc_real_gas.py` | Huang's original P_c* chain, ideal gas vs the R141b EOS (**needs CoolProp**; not a test). |
| `data/huang1999_reference.json` | Canonical golden fixture (inputs, outputs, Jacobian, theory). |
| `data/huang1999_reference_data.h` | Same data as a `constexpr` array for ctests. |
| `data/henry2007_fig5.py` | Henry 2007 Fig. 5, measured air characteristic, and the rig's geometry. |
| `data/akbarnejad2025_tables.py` | Akbarnejad & Ziabasharhagh 2025 Tables 1-3, steam CFD rows. |
| `data/chen2013_figs.py` | Chen et al. 2013 Figs. 5a and 7a: Chen's model curves and the measured points (Hemidi/Henry rig; Chen's rectangular rig, with error bars). |
| `data/eames2013_fig8.py` | Eames et al. 2013 Fig. 8 and Table 1: R245fa chiller, two ducts on one nozzle, measured T_c* and COP, the authors' design point; R245fa properties baked. |
| `data/generate_eames2013.py` | Regenerate those baked R245fa values (**needs CoolProp**). |
| `eames2013_check.py` | Cross-source accuracy of ω and both P_c* closures on that rig. |
| `data/mazzelli2015_tables.py` | Mazzelli & Milazzo 2015 Tabs. 3-4: the same rig's third duct (A_3/A_t 9.72), measured mass flows with 95 % bands; the 2014 original-nozzle point; R245fa properties baked. |
| `data/generate_mazzelli2015.py` | Regenerate those baked R245fa values (**needs CoolProp**). |
| `mazzelli2015_check.py` | Cross-source accuracy of ω against measured flows, the primary choke, and both P_c* closures against the measured bracket. |
| `anchor_scorecard.py` | Cross-source accuracy of omega, P_c*, P_b0 on all three rigs (run it to print the table). |
| `printed_chain_check.py` | Huang's printed P_c* chain and single-change hypotheses against Huang's and Chen's own numbers; Chen's printed subcritical branch against his droop. |
| `bouhanguel2011_regime_check.py` | Bouhanguel et al. 2011: which anchor pair agrees with the observed flow regime, with P2 bracketed by the measured U. |
| `data/lucas2012_choke_check.py` | Lucas & Koehler 2012 transcritical CO2: ideal-gas vs two-phase choke against the measured driving flow (**needs CoolProp**; not a test). |

The Python regression test lives at `python/tests/test_ejector_huang.py` and
consumes the JSON fixture directly.

## Golden fixture contract

Each case in `huang1999_reference.json` (and the `HuangCase` struct in the
header) carries:

- `inputs` — everything the model needs, including the baked `gamma`.
- `reference` — this model's `omega` and intermediates; the C++ **forward solve
  must match these tightly** (same math).
- `jacobian` — central-difference `d(omega)/d(input)`; the target the C++
  **analytic Jacobian must match to finite-difference tolerance**.
- `paper_theory` — Huang's published theory `omega` (physics check, `omega_rel_tol`
  = 3%) and each row's `T_c*`-derived `p_c_pa` target (physics check,
  `pc_rel_tol` = 15%; see Accuracy below).

## Regeneration workflow

```bash
# Rare: only when the fluid/gamma recipe changes (needs CoolProp).
uv run --with CoolProp python validation/ejector/data/generate_gamma.py
#   -> paste values into GAMMA_BY_SUCTION_C in huang1999_tables.py

# Routine: after any change to the model or table (no CoolProp).
uv run python validation/ejector/data/generate_reference.py
#   -> rewrites huang1999_reference.json and huang1999_reference_data.h

uv run pytest python/tests/test_ejector_huang.py
```

## Accuracy

The reference model reproduces Huang's theory column across all 31 rows (both
primary pressures 0.40-0.60 MPa, both suction conditions, both nozzles, omega
0.16-0.86) to a **mean 0.8%, max ~2%**. The residual is the hard floor of a
constant-gamma ideal-gas model (the true ideal-gas gamma drifts 1.091 -> 1.123
across the expansion) plus table round-off -- not an implementation error. That
floor only comes down by dropping the ideal-gas assumption, which is the job of
the production element on combaero's real-fluid thermo.

### P_c*: why the model changed, and what was tried

Huang's original P_c* chain (Eqs. 9-18: a `phi_m`-weighted momentum balance,
then a separate Rankine-Hugoniot shock, then isentropic diffuser recovery)
was systematically **25% low (max 35%)** against each row's `T_c*`-derived
target. Investigation ruled out gamma by direct test (P_c* is monotonic in
gamma over its entire valid range; the gamma -> 1 limit is still ~20-30%
short; re-anchoring gamma's reference temperature moves P_c* under 1%) and
ruled out an equation error (hand-verified against the printed paper; Eqs.
16-17 reproduce classical normal-shock-table values exactly). That left the
momentum treatment itself as the leading suspect: reproducing the target
needed `phi_m` ~ 0.90-0.95, above Huang's own stated calibrated range
(0.80-0.84), even for rows where `phi_m` sat well inside that range.

**Real-fluid properties are ruled out too** (`data/huang_pc_real_gas.py`).
Huang's working fluid is R141b, and his primary is far from ideal (saturated
vapour at 0.40-0.60 MPa, Z = 0.87-0.90). The script runs his chain on the
CoolProp R141b EOS: real-gas velocities, mixing enthalpy, a Rankine-Hugoniot
shock solved on the EOS, and an isentropic diffuser on the EOS. P_c* against
the measured-`T_c*` target, all 31 rows:

| Chain | mean | range |
|---|---|---|
| ideal gas, constant gamma | -22.7% | -29.3..-16.1% |
| EOS mixing, shock and diffuser (y-y pressure and omega from the ideal chain) | -24.1% | -31.0..-17.0% |
| EOS throughout (also the choke pressure and both mass fluxes) | -22.9% | -29.7..-16.2% |

- **Real gas moves P_c* by at most 1.5 points, and not toward the target.**
  Mixing, the shock and the diffuser all sit at 20-150 kPa, where R141b is
  close to ideal. The non-ideal primary inlet enters mainly through the
  primary mass flux, which shifts omega rather than the pressure recovery.
- **The fully real-gas omega lands 6% below Huang's theory omega**
  (-2..-14%), while the ideal chain reproduces it to 0.8%. So Huang computed
  with ideal-gas relations, as his assumption 1 states, and the miss is in
  the chain as printed rather than in its properties.
- **The script is not flat:** with `phi_m` forced to 0.93 the EOS chain
  lands at -7..+4%, the same response as the ideal chain.
- **The R141b vapour line is slightly wet here.** On the EOS, the
  saturated-vapour line has a negative slope at 8-12 C, against the paper's
  statement that R141b is a dry fluid. An equilibrium isentrope reaches
  quality 0.985-0.994 by the y-y plane, so the script holds the vapour
  metastable. A 1% liquid fraction cannot move P_c* by 23%.

Four independent models/reformulations were then implemented against the
same y-y starting state and tested against all 31 rows:

| Approach | mean \|err\| | max \|err\| | sign |
|---|---|---|---|
| Huang (original: `phi_m` momentum + shock + isentropic diffuser) | 25% | 35% | low |
| Lienhard/McGovern (zero-loss momentum, static-T energy, no shock) | 9.0% | 16.0% | high |
| **Kracik & Dvorak (zero-loss, stagnation-referenced, no shock) -- adopted** | **6.2%** | **13.2%** | high |
| Akbarnejad (drops momentum, imposes exact-sonic condition) | 32-60% | -- | low |

Lienhard/McGovern and Kracik & Dvorak both retain momentum in full (mass +
momentum + energy solved together) but drop `phi_m` and the separate shock
step -- both beat Huang's chain by 3-4x and, tellingly, **flip sign**
(over-predict where Huang under-predicts): a fully lossless treatment is a
theoretical upper bound, so landing above target while Huang's damped
treatment lands below is exactly the expected bracketing. Kracik & Dvorak's
`lambda`/`q(lambda)`/`z(lambda)` formalism edges out Lienhard/McGovern
because it preserves each stream's stagnation enthalpy along its own
isentropic path before combining, while Lienhard/McGovern's energy equation
explicitly drops the kinetic-energy term (a mass-weighted *static*
temperature average only, confirmed by re-reading their printed equation
directly).

Akbarnejad & Ziabasharhagh's CFD independently measured the same defect from
a different angle: their `Psi_5` metric (the relative error of assuming
momentum is simply conserved across the mixing region) reaches **up to 30%**
at critical back pressure in their own water-vapor ejector -- close to
Huang's ~25-35% gap here, and real, CFD-quantified confirmation that the
momentum treatment specifically is the weak link. Their own proposed fix
(drop momentum entirely; close on an *externally imposed* exact-sonic
condition) was implemented and tested too, but performed worse than Huang's
own chain (32-60% depending on the reading of their ambiguous mixing-pressure
convention) -- their model is explicitly a geometry-free **design** tool (it
never consumes a given area ratio, only inlet conditions), not a
given-geometry forward analysis like this one, and the paper says so
directly. Dropping momentum is not the right direction for this use case.

**Adopted then: Kracik & Dvorak's closure.** It is no longer the element's
default: on the independent rigs below, Huang's printed chain is the better
P_c* despite its fidelity gap on Huang's own rows, and it took over as
`anchor_closure="huang1999"`. Kracik & Dvorak's closure was adopted with
`recovery_efficiency` defaulting
to **1.0** (no artificial loss) rather than a value fitted to this dataset --
fitting `phi_m` or `recovery_efficiency` against Huang's own 31 rows would be
exactly the "never fit against measurements" trap this project avoids
elsewhere; a principled default, if one is ever adopted, should come from
further validated literature or real test/CFD data, not a curve fit to this
one dataset. `test_ejector_huang.py` checks P_c* by regression, by trend
(monotonicity in `A_3/A_t`, and `P_e < P_c* < P_g`), and now also by a real
(if loose) physics tolerance against `paper_theory.p_c_pa` (`pc_rel_tol` =
15%, comfortably above the observed 13.2% max).

### Cross-source accuracy of the anchors

Everything above is the critical closure against Huang's own theory columns,
for R141b. `anchor_scorecard.py` runs the closures, as the network element
takes them (γ and R at the choke plane, `recovery_efficiency` = 1), on three
rigs that are not their sources. It is **accuracy, labelled cross-source**:
a miss is the model's limitation.

| Rig | Kind | A_3/A_t | ω | P_c* | **P_b0** | P_b0/P_c*: model / rig |
|---|---|---|---|---|---|---|
| Henry 2007, air, P1 absolute | measured | 5.30 | +11 to +22% | +20 to +31% | **+78 to +83%** | 1.80-1.82 / 1.21-1.32 |
| Henry 2007, air, P1 gauge | measured | 5.30 | −8 to −11% | +34 to +45% | **+93 to +108%** | 1.76-1.82 / 1.21-1.32 |
| Chen 2013 rectangular, air, P_p absolute | measured | 4.46 | +13 to +16% | +14 to +21% | **+62 to +68%** | 1.66-1.70 / 1.15-1.24 |
| Akbarnejad 2025 E1, steam | CFD | 5.44 | −4 to +43% | +38 to +63% | **+86 to +135%** | 1.44-1.74 / 1.06-1.21 |
| Akbarnejad 2025 E2, steam | CFD | 33.3 | −18 to −22% | +18 to +21% | **+257 to +274%** | 3.17-3.37 / 1.06-1.07 |

- **The dead-head anchor P_b0 is the weak point.**
  - It is Kracik & Dvorak's mixing at ω = 0: lossless, and independent of A_3/A_t by construction.
  - So the linear droop's far end runs about twice above every rig's breakdown, and 3.7 times above it at A_3/A_t = 33.
  - Over the window the straight line therefore keeps too much entrainment. The roots the second law refuses (#598) are this: stress draws at A_3/A_t of 15 to 33, deep in the droop, entraining more than the primary jet can drive.
- **P_c* lands 14-63% high.** That is above the reference, as a lossless closure must be, and about where Huang's own rows put it (+6% mean, 13% max).
- **Henry's droop shape is still linear** (R² ≥ 0.94, `test_ejector_subcritical.py`). The chord's anchor is the error, not its shape.
- **Caveats:**
  - Henry's ω depends on whether P1 is gauge or absolute, so both are shown. The efficiency-free test (ratios of plateau ω across P1, which cancel every efficiency coefficient) reads it as **gauge**: ratio RMS 0.018 against 0.048 for absolute. The same test reads Chen's rectangular rig as absolute (0.010 against 0.038). Both are pinned in `test_ejector_chen2013.py`.
  - Chen's rectangular rig: the knee is the midpoint between the last plateau point and the first droop point. The curves stop at ω 0.13-0.45, so its breakdown is the last two points extrapolated to ω = 0.
  - The CFD's stream temperatures are assumed saturated, because the paper does not print them. ω follows that assumption, and P_b0/P_c* hardly does.

**A candidate dead-head, not adopted yet.** `candidate_dead_head` in
`anchor_scorecard.py` treats ω = 0 as a sudden expansion. The primary leaves
the nozzle (M_p1, P_p1 on A_p1) into the mixing tube, and the dead-headed
suction annulus A_3 − A_p1 sits at the secondary pressure P_e. Mass, momentum
and energy over A_3 then give a subsonic exit. This is Keenan & Neumann's
zero-secondary-flow constant-area mixer: no fitted constant, and A_3/A_t
enters.
- **P_b0 accuracy:** 3-39% high on every row against today's 62-274%. Henry absolute is +3-4%, and E2 at A_3/A_t = 33 is +12-14%.
- **On a rig added after its form was fixed** (Chen's rectangular rig): +3 to +5% on all five curves.
- **The catch:** against today's P_c*, which is 18-63% high, its window P_b0/P_c* is only 0.98-1.06, where the rigs show 1.06-1.32. Adopted alone, it would make the droop vanish, and invert it on some rows.
- **Paired with a P_c* that is not lossless** (Huang's printed chain, next paragraph), the window comes out positive on every rig.
- **Variant tried and set aside:** the primary core at the y-y state, with the annulus at P_py. It fell 12-24% below the breakdown and below P_c*.

**The candidate pair: Huang's printed P_c* chain with the sudden-expansion
P_b0.** The scorecard's "cand" columns, also not in the model. The printed
chain (`printed_chain_check.printed_chain`) adds Huang's φ_m mixing loss and
normal shock, so it lands below the lossless closure:
- **P_c\*:** −18 to +1% on the measured rigs and −19 to +17% on the CFD, against +14 to +63% today.
- **Window P_b0/P_c\*:** 1.37-1.65 everywhere, so the droop never collapses or inverts. That includes A_3/A_t = 33, where today's window is 3.2-3.4.
- **Whole measured curve** (the element's plateau and straight chord against every measured point, RMS over the plateau): 12-23% on all 11 measured curves, against 29-59% today.
- **The catch** is fidelity: as printed, the chain does not reproduce its authors' own numbers (next section).

**The element's default: `EjectorElement(anchor_closure="huang1999")`**
(`huang_critical_back_pressure` and `sudden_expansion_dead_head`, with analytic
Jacobians in C++). `"kracik_dvorak"`, the only pair up to 0.6.0, stays
selectable:
- **Network level.** Every measured back pressure of both measured air rigs solves with either pair (91 of 91). The curve RMS over the plateau drops from 29-48% to 12-23% on all eight curves (`test_ejector_anchor_closure_element.py`).
- **Regime, a fourth lab** (Bouhanguel 2011, below): only Huang's P_c* agrees with the regimes the authors saw without an edge-case reading.
- **Stress family `ejector`, seeds 0-400.** Run A (`kracik_dvorak`) against run B (`huang1999`), same seeds and jobs:
  - **Second law.** Run A converges 52 roots that the second law refuses (#598: σ < 0 at the ejector, an entrainment the primary cannot drive). Run B has none.
  - **Successful solves:** 347 → 277. Run B loses 122 draws:
    - 98 do not converge;
    - 24 converge to a root with the suction flowing backwards (ω < 0), which #619 refuses.
    - The 98 and 23 of the 24 have their outlet boundary 2-107% above the new dead head, where a real ejector backflows and the model does not reach. The 24th sits on the dead head (ω = −0.004), its outlet pipe's loss putting the port just above it.
    - On the measured rigs that dead head sits 3-5% above the breakdown. The solver's message names it (`EjectorElement.no_root_hint`).
  - **Gained draws:** none.
  - **Time.** The median is unchanged (10 ms) and p90 improves (450 → 300 ms).
- **Incompressible ejectors** use no anchors, but their compressible warm-start proxy (`NetworkSolver._compressible_seed`) does.
  - Under Huang's narrower window that proxy converged on 1 of the 61 incompressible draws that call it, against 24 under K&D, and the family took 51 s instead of 40 s.
  - So the proxy runs with K&D's window (`EjectorElement.compressible_proxy_overrides`): 24 proxy solves converge again, the family takes 40 s again, and the same 190 of 191 solve.
  - **Set aside:** skipping the proxy and leaning on the Newton-homotopy retry also solves 190 of 191, in 11 s. But that retry is opt-in (the GUI and the stress harness turn it on). With it off, the proxy is what rescues them: 154 of 191 solve with it, 130 without.

### Why the printed chain misses its authors' own numbers

Huang's original chain (his Eqs. 9-18 with φ_m, Eq. 19) is the best P_c* on
every independent rig (the candidate pair, above), but as printed it does not reproduce the
numbers its authors computed with it. Two authors, the same chain:

- **Huang 1999.** His theory A_3/A_t is what his code solved to hit each row's
  measured T_c*. Run forward at that area, the printed chain lands **−22.7%**
  (−29..−16%). Run in his design direction (P_c* in, A_3/A_t out) it misses
  his A_3/A_t by −38%.
- **Chen 2013** (`data/chen2013_figs.py`). Chen runs the same chain, same
  coefficients, on two air rigs. His model's ω we reproduce to +1-2%
  (Eqs. 1-8, `test_ejector_chen2013.py`), but its P_c* (the plateau end of his
  curves) sits **8-18% above** the printed chain on all 8 curves.
- **Effective φ_m.** Both authors behave as if φ_m were about 0.88-0.91
  (Chen, all 8 curves) or 0.90-0.95 (Huang), against the printed 0.80-0.84.

These are model-against-model comparisons: each author's code took the same
inputs we do, so no rig bias (gauge readings, thermocouple error) can enter.
`printed_chain_check.py` tests single changes their codes could have made:

| One change to the printed chain | Huang, 31 rows (row std) | Chen, 8 curves (in his bracket) |
|---|---|---|
| none (as printed) | −22.7% (3.0%) | −12.7% (0/8) |
| φ_m as a kinetic-energy efficiency (√φ_m) | −5.2% (3.4%) | +3.9% (2/8) |
| φ_m on the primary stream only | −17.9% (4.6%) | −2.6% (2/8) |
| φ_m on the slip velocity only | +8.5% (2.6%) | +17.7% (0/8) |
| mixed T from static enthalpies | −15.7% (2.6%) | −2.1% (3/8) |
| shock at the primary-core Mach | +71% (22%) | +109% (0/8) |
| no shock, isentropic recovery from the mix | −9.1% (5.7%) | −11.9% (0/8) |
| incompressible diffuser | −24.2% (2.6%) | −18.0% (0/8) |
| diffuser exponent inverted | −37.9% (1.7%) | −43.3% (0/8) |
| M_3 read as M_3² (Eq. 17) | −32.1% (2.4%) | −22.7% (0/8) |
| a_py at the stagnation T_g (static/total slip) | +1.3% (3.0%) | +44% (0/8) |
| a_sy at the stagnation T_e | −22.1% (3.1%) | −7.8% (0/8) |
| V_py from the nozzle-exit state | −31.2% (5.3%) | −17.3% (0/8) |
| a_m at the mixed stagnation T | −30.9% (3.9%) | −26.4% (0/8) |

- **The test for a real explanation:** it should bring the mean to zero AND collapse the row-to-row scatter to table rounding (about 0.5%). Two authors with their own codes could have made different slips, so a hypothesis need not fit both.
  - The closest for Huang alone is the sound speed at the stagnation T_g: mean +1.3%, but the scatter is unchanged at 3.0%.
  - The closest for Chen alone is the mixed temperature from static enthalpies: −2.1%, 3 of 8 inside his bracket.
  - Neither collapses the scatter, so each is a coincidental shift of the mean, not the explanation.
- **Unit slips inside their codes are excluded too** (`printed_chain_check.py`):
  - **A pressure offset** (gauge read as absolute, a "+1 bar" approximation, 0.1 MPa against 1 atm) would shift P_c* by a constant. The authors sit +0.25 ± 0.07 bar (Huang) and +0.19 ± 0.06 bar (Chen) above the printed chain: not constant, and nowhere near 1 bar.
  - **A gauge P_p used only in the nozzle-exit relation** moves P_c* the wrong way (−26.5% Huang, −18% Chen).
  - **Temperature units.** The chain gives the identical P_c* in K and in °R, since only Mach numbers and T_e/T_g enter it. °C or °F would scale ω by about 3× through T_g/T_e, and Huang's ω matches ours to 0.8%.
- **Already ruled out:**
  - γ (above);
  - an equation-transcription error (each equation checked against the PDF);
  - real-fluid R141b (`data/huang_pc_real_gas.py`, above);
  - an inconsistent property set. Neither paper states the state at which C_p and γ were taken (Huang: "constant properties"; Chen: "computed by the NIST database"). Taking both at the generator, evaporator or condenser saturated-vapour state moves Huang by at most 3 points (−19.7 to −22.0%; `huang_pc_real_gas.py`). Air is ideal at these conditions, so properties cannot explain Chen's gap at all.
- **Not ruled out:** an unprinted step in their codes. Only the authors could
  confirm one.
- **Huang's gap is about twice Chen's,** but no γ-dependent change found so
  far explains the difference.

**Chen's subcritical branch fails its fidelity check too.** As printed (the
entrained flow unchokes, P_sy rises until the chain's back pressure equals
P_c), it gives a droop much shallower than Chen's own plotted one: ω +0.1 to
+0.2 high one or two points into the droop (`printed_chain_check.py`, second
table). So it is not a candidate for the droop.

**Measurement-side checks (model against rig).** These matter for accuracy, not for the fidelity gap above:
- **Pressure reading.** Huang's pressures are absolute: each printed P_g and P_e matches its printed saturation temperature to 0.1 K. Henry's are gauge and Chen's rectangular rig is absolute (ratio test, above).
- **Suction loss.** On Chen's rectangular rig, one suction pressure of about 0.91 bar (−10%) reproduces all five measured plateaus. That would fit a constant loss in the secondary's metering run (an ISO 5167 orifice) at the roughly constant ṁ_s. A lower η_s fits the plateaus equally well; the two separate only on P_c*, which a lower P_s lowers and η_s barely moves.

### A refrigerant rig: Eames 2013 (R245fa)

`eames2013_check.py` scores the kernel on a fourth lab: Eames, Milazzo,
Paganini & Livi 2013 (`data/eames2013_fig8.py`, Fig. 8 digitised and checked by
the user). It is **accuracy, labelled cross-source**, and it is the only
refrigerant rig scored.
- **The rig:** a 40 kW R245fa chiller with two mixing ducts, A_3/A_t 8.6 (A) and 13.2 (B), fitted to ONE primary nozzle.
- **The data:** ten critical points, giving the measured T_c* (a pressure) and the plateau COP.
- **The paper has no nozzle exit diameter. That does not matter:** ω and both P_c* closures are exactly independent of A_p1/A_t.

| Ejector | A_3/A_t | P_c\* Huang chain (default) | P_c\* K&D | ω / ω(COP) |
|---|---|---|---|---|
| A | 8.6 | −5 to −12% (bias −8.9%) | +30 to +34% (+31.1%) | 1.4-1.9 |
| B | 13.2 | −14 to −24% (bias −19.5%) | +3 to +15% (+8.9%) | 1.8-3.3 |
| all | | bias −14.2%, RMS 15.4% | bias +20.0%, RMS 23.1% | 2.0 |

- **The two closures bracket every measured P_c\*.** Huang's chain is low on all ten points and closer on average. That is the same pattern as on the air rigs.
- **Geometry, not the closure.** Going from duct A to duct B with the same nozzle, BOTH closures drop by 10-22 points, and the ω overprediction grows.
  - The closures share only Huang's Eqs. 1-8 (the y-y state and ω), so the effect sits in that front end.
  - The perturbation confirms the duct drives it: with B scored at A's ratio, the shift reverses.
  - **Candidates:**
    - the entrained flow filling all of A_3 − A_py at y-y;
    - the nozzle's position and exit area, which the kernel cannot see (the paper's Fig. 5 moves COP by 30-45% over ±4-6 mm of NXP);
    - mixing length (B's diffuser throat sits 32 mm further downstream).
  - The other rigs show no clean trend with A_3/A_t. Their errors (table above) are confounded by fluid and conditions, and only this rig changes the duct alone.
- **ω against the COP is scored, with the authors' doubt attached.**
  - The model is 1.4-3.3× the ω the plotted COP implies.
  - The authors doubt that COP: chilled-water fluctuations and a small evaporator ΔT. Their own 1-D code is also about 2× their measurement (0.31 against 0.15).
  - So this column carries their stated uncertainty, not a measurement band.
  - The COP-to-ω conversion itself checks out: their code's COP 0.31 maps to ω 0.45 through the CoolProp enthalpies, against the 0.44 they print.
- **Model against model:** at their design point (106/39/4 °C) the kernel gives ω 0.11 on ejector A, the duct their code designed, against the code's 0.44. Their code is a different front end, not Huang's.

### The same rig's third duct: Mazzelli & Milazzo 2015 (#697)

`mazzelli2015_check.py` scores the third ejector of the same chiller
(`data/mazzelli2015_tables.py`). The source is Mazzelli & Milazzo, *Int. J.
Refrigeration* 49 (2015) 79-92, Tables 3-4, read from the accepted manuscript.
The geometry also appears in Table 5.3 of Grazzini, Milazzo & Mazzelli,
*Ejectors for Efficient Refrigeration* (2018). It is **accuracy, labelled
cross-source**, and nothing is fitted.
- **The rig:** a 31.8 mm diffuser throat on the same 10.2 mm primary throat, so A_3/A_t = 9.72, between ducts A (8.65) and B (13.2).
- **The data are tabulated, not digitised:** 8 points in two back-pressure sweeps, T_g 89-90 °C and T_e 5 and 10 °C.
  - **The bands:** measured primary and secondary flows carry 95 % bands from the authors' uncertainty analysis. So ω = m_s/m_p has a **measurement** band of ±9-14 %, the first on this rig.
  - **What is read from the data alone:** the plateau, i.e. points within the band of the sweep's largest ω. Also the T_c* bracket: [28.3, 28.8] °C and [29.4, 30.2] °C.

| Quantity | Kernel | Measured | |
|---|---|---|---|
| ω (5 plateau points) | 0.38-0.51 | 0.51-0.72 ± 9-12 % | **−25 to −30 %** (bias −27.8 %), 0/5 in band |
| primary flow | 0.326-0.333 kg/s (Eq. 1, γ 1.084) | 0.35 ± 0.03 kg/s | −5 to −8 %; the real-gas choke (0.355-0.363) is inside |
| P_c\* Huang chain (default) | | bracket | −18.5 %, −14.7 % outside; −17.6 % to the midpoint |
| P_c\* K&D | | bracket | +10.6 %, +12.7 % outside; +13.0 % to the midpoint |

- **ω flips sign against #687.**
  - On ducts A and B, the kernel was about 2× the ω derived from the authors' doubted COP.
  - Against this rig's measured flows, it is 28 % LOW.
  - The 2014 original-nozzle point (T_g 100 °C, COP 0.23, which they did not flag) reads 25 % low too.
  - So the overprediction in #687 most likely reflects that COP, not the kernel. The two are different ducts, so this is not proven.
- **Part of the ω gap is the primary flow.** The kernel's ideal-gas choke at the ω-convention γ is 5-8 % below both the measurement and the real-gas choke. That leaves the entrained flow about 30 % low.
- **The two closures again bracket the measured P_c\*.** Huang is below the bracket and K&D above it, as on Eames 2013 and the air rigs.
- **PRE-REGISTERED and MISSED.** Interpolating #687's A → B shift linearly in A_3/A_t predicted Huang −11.5 % and K&D +25.7 % at 9.72.
  - **What happened:** both land 6-13 points lower. The third duct reads closer to B than its area ratio suggests.
  - **Confounders:** the conditions differ (T_g 89-90 °C here, 100-110 °C in Eames 2013), and the third duct is longer (950 mm, a 272 mm mixing region).
  - **What follows:** a linear trend in A_3/A_t alone does not carry over between ducts. Reported, not refitted.
- **The kernel cannot tell the two nozzles apart.** ω and P_c\* are independent of A_p1/A_t (pinned in `test_ejector_eames2013.py`). The authors' CFD gives +16 % ω and +2 K T_c* for the smaller exit (book Fig. 5.15). The 2014 point and the tables also differ in T_g, so they are not a clean A/B.
- **Sensitivity (falsification):**
  - Scoring the same points at A's ratio gives ω −42 % and P_c\* Huang −13 %; at B's ratio, ω +18 % and Huang −28 %. The metric follows the geometry.
  - PHI_P from 0.88 to 0.80 moves ω to −22.6 % and Huang to −18.9 %.
- **Not a reference:** the paper's CFD matches these data only with a fitted 20 µm sand-grain roughness, against 4-6 µm measured Ka. It is tuned, so it is not scored against.

### Other sources reviewed

None of these has a back-pressure sweep at fixed motive and suction
pressures, the only kind of data that scores the anchors directly. What each
does give:

- **Bouhanguel, Desevaux & Gavignet 2011** (Int J Refrig 34:1633; air; d* 8, d 12, D 24 or 16 mm, discharging to atmosphere). Laser-sheet visualisation, so the data are regime statements plus two entrainment ratios. `bouhanguel2011_regime_check.py`, accuracy, cross-source:
  - **What it states:** A_3/A_t = 9 runs with the secondary *not* choked at P1 = 4-6 bar (subcritical at the atmospheric outlet), while A_3/A_t = 4 chokes at 5.5 bar. It also gives U = 0.47 and 0.50 at P1 = 5 and 4.5 bar.
  - **How P2 is bracketed:** P2 is not printed (a valve throttles the ambient intake, so P2 ≤ Pa). A measured U bounds it from below, because a point that is not choked lies on or below the plateau at its own P2.
  - **Huang:** fits both A_3/A_t = 9 points over wide P2 windows, with P1 read as absolute or as gauge.
  - **K&D, gauge reading:** contradicts "not choked" at both points.
  - **K&D, absolute reading:** fits P1 = 4.5 bar only if U is at least 73% of the plateau, and P1 = 5 bar only if U is at least 99.9% of it (a 0.04% window in P2).
  - **A = 4:** both pairs call it choked.
  - **A model limit, not scored:** the measured unmixed length reaches 9.2 D at high P1, longer than the 8.3 D tube. Mixing is not complete inside it, while the 1-D model assumes it is.
- **Lucas & Koehler 2012** (Int J Refrig 35:1595; transcritical CO₂ ejector cycle). It is **outside the element's scope**:
  - the motive flow is supercritical CO₂ at 73-103 bar and 500-700 kg/m³, flashing to two phases in a converging nozzle;
  - the ideal-gas choke is 34-51% below the measured driving flow, against +2 to +8% for an equilibrium two-phase choke on the CoolProp EOS (`data/lucas2012_choke_check.py`; Fig. 5 read by eye);
  - the outlet pressure is set by the cycle, so no plateau end or breakdown is measured.
- **Liao 2008** (PhD dissertation, Texas A&M). **No usable data:**
  - its one measured set holds suction and outlet at atmosphere and varies only the supply pressure, through a converging orifice primary;
  - its entrainment ratio is defined as (ṁ_p + ṁ_s)/ṁ_p, and the data are plotted only;
  - its 1-D models are never compared with measurements.

### Validated envelope

Where each closure has been scored, and against what. Outside these ranges a
result is extrapolation, and nothing in the element flags it.

| Closure | Scored against | Range |
|---|---|---|
| critical ω, P_c* (fidelity) | Huang 1999 theory columns (R141b, γ ≈ 1.11) | A_p1/A_t 2.9-3.3, A_3/A_t 6.6-10.9, P_g/P_e 8.5-15.1, ω 0.16-0.86 |
| critical ω, P_c*, P_b0 (accuracy) | Henry 2007 and Chen 2013 rectangular (air, measured), Akbarnejad 2025 (steam, CFD) | A_3/A_t 4.5-5.4 and 33; P_p/P_s 3-35; ω 0.19-1.4; table above |
| P_c*, P_b0 of `anchor_closure="huang1999"` (accuracy) | Henry 2007 and Chen 2013 rectangular (air, measured), Akbarnejad 2025 (steam, CFD) | as the row above; P_c* −19 to +17%, P_b0 +3 to +39%, curve RMS 12-23% |
| critical ω, P_c* (second author's own model) | Chen 2013 model curves (air) | A_3/A_t 4.5 and 5.3; ω reproduced, P_c* not (above) |
| critical ω, P_c* of both pairs (accuracy) | Eames 2013 (R245fa, measured T_c*; COP doubted by the authors) | A_3/A_t 8.6 and 13.2, P_g/P_e 15-19, T_c* 32-37 °C; P_c* Huang −24 to −5%, K&D +3 to +34% |
| critical ω (measured band), primary choke, P_c* bracket, both pairs (accuracy) | Mazzelli & Milazzo 2015 (R245fa, same rig, third duct; measured flows) | A_3/A_t 9.72, P_g/P_e 12-15, T_c* 28-30 °C; ω −25 to −30 %, P_c* Huang −18 %, K&D +13 % |
| regime at a fixed outlet, both pairs (accuracy) | Bouhanguel 2011 (air, visualisation) | A_p1/A_t 2.25, A_3/A_t 4 and 9, P1 4.5-5.5 bar; only Huang's P_c* fits without an edge-case reading |
| droop shape (linear) | Henry 2007 (air, measured) | P_p/P_s about 4-7, one geometry |
| jet pump (incompressible) | Sanger, NASA TN D-4759 / D-4445 (water) | A_p1/A_3 0.066 and 0.197 |
| unchoked jet pump (compressible) | analytic low-Mach limit only | no measured data |

Back-flow (outlet above P_b0) is outside the model. The element refuses a
converged root with a reversed port (#619) and one that violates the second
law (#598).

Also outside it, and not flagged by the element:
- **A dense or flashing motive fluid**, such as transcritical CO₂. The ideal-gas choke is 34-51% low there (Lucas 2012, above).
- **A mixing tube too short to finish mixing.** The 1-D chain assumes a fully mixed exit, but Bouhanguel 2011 measured the unmixed core reaching 9.2 D.

## Roadmap

1. ~~Add the P_c* / critical-back-pressure chain to the reference model.~~
   Done -- see Accuracy above. Kracik & Dvorak's closure first; Huang's
   original Eqs. 9-18 (`anchor_closure="huang1999"`) are the default since
   the cross-source scoring.
2. ~~Draft the production 3-port `EjectorElement` (primary-in, secondary-in,
   outlet) in `python/combaero/network/`, real-fluid thermo, analytic `(f, J)`.~~
   Done (PR #252).
3. ~~Port to C++ with analytic Jacobians (derived symbolically with sympy where
   useful) and ctests consuming `huang1999_reference_data.h`'s central-
   difference Jacobian columns to tight tolerance.~~ Done (PR #253).
4. ~~pybind bindings, in a separate PR.~~ Done (PR #254).
5. ~~GUI work and validation.~~ Done (PR #259).
6. ~~**Operating-regime extension:** choked/unchoked primary (R0), subcritical
   entrainment droop (R1), and the unchoked subsonic jet-pump mode, in one C1
   residual system so critical <-> subcritical <-> jet-pump solve in a single
   Newton pass.~~ Done: reference closures (Phase 1A/1B) and the whole-element
   C++ `(f, J)` composed into `EjectorElement` (Phase 2). Full findings,
   design, and provenance in `OPERATING_REGIMES_DESIGN.md`.
