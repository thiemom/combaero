# Ejector operating-regime extension: choked/unchoked primary, subcritical entrainment, and unchoked jet-pump mode

Design + findings note for extending `EjectorElement` beyond the critical
(double-choked) plateau. Written 2026-08-26. Audience: developers of the
network module (not end users). Captures the diagnosis of a real
non-converging low-primary-flow case and proposes the model extension that
resolves it. Nothing here is implemented yet -- this is the spec that the
implementation PRs will follow.

Companion docs:
- `README.md` (this directory) -- the reference model (critical plateau +
  subcritical/jet-pump closures), accuracy comparison, golden-fixture contract.
- `../../python/combaero/network/_ejector_huang1999.py` -- the closed-form
  physics across regimes (Huang 1999 + Kracik & Dvorak 2016).

---

## 1. Motivation: a physically posed network that does not converge

A GUI network (`gui/tmp/ejector_test_low_flow_not_converged.json`) drives an
ejector with a **forced primary mass flow of 0.0051 kg/s** (a `mass_boundary`),
a secondary suction reservoir at 101 325 Pa, and an outlet reservoir at
101 325 Pa. The solver reaches a root (|F| ~ 2.5e-10) but `EjectorElement`'s
`verify_solution_consistent` demotes it to non-converged (outlet Pt exceeds the
computed critical back pressure P_c*). The user's intuition -- that a physical
solution exists -- is correct. The failure is a **modeling-regime** limitation,
not a solver bug, and the root cause is deeper than "subcritical back
pressure": the current residual set has no branch for an **unchoked primary
nozzle**, and picks a self-contradictory choked root.

## 2. Findings

### 2.1 R0 hard-codes the choked branch and lands on an unphysical root

The element's primary-flow residual is (see `ejector_element.py`):

    R0 = mp - ejector_choked_mass_flow(Pt_primary, Tt_primary, A_t, gamma, R, eta_p)

`ejector_choked_mass_flow` is the **always-choked** (sonic-throat) relation:
`mdot = const * Pt * A_t / sqrt(Tt)`, monotone in `Pt`. For the forced
`mp = 0.0051 kg/s` it back-solves to `Pt_primary = 71 552 Pa`. But a nozzle can
only be choked when `Pt/P_back >= ((gamma+1)/2)^(gamma/(gamma-1)) ~= 1.89`. At
`Pt = 71 552 Pa` against a ~101 kPa mixing/back pressure the ratio is 0.71 --
the root **asserts choked flow at a pressure that cannot choke**, and puts the
"motive" stream *below* the pressure it is supposed to discharge against.

The physically consistent root is the **unchoked / subsonic** one:

| | Primary Pt | Pt/Pb | State | dp = Pt - Pb |
|---|---|---|---|---|
| R0 choked root (solver artifact) | 71 552 Pa | 0.71 | "choked" (impossible) | **-29 773 Pa** (backward) |
| Physical subsonic root | ~102 435 Pa | 1.01 | unchoked, M_throat ~ 0.13 | **+1 110 Pa** |

For forward flow the pressure drop is `>= 0` and the primary supply sits
*slightly above* both the secondary and outlet pressures -- exactly what a weak
forced flow through a subsonic nozzle requires. R0 has no unchoked branch, so
Newton falls into the bogus low-pressure choked root.

The choke threshold is sharp: at `Pb = 101 325 Pa` the choked mass flow is
**0.00722 kg/s**. Any forced `mp` below that cannot choke and the primary Pt
must float up to just above `Pb`; any `mp` above it chokes and Pt rises well
above `Pb`. The observed "convergence returns" at `mp ~ 0.00722 kg/s` in a
primary-flow sweep is this choke threshold, not a coincidence of P_c*.

### 2.2 The regime hierarchy (ordered by forced primary mass flow)

At fixed geometry, secondary condition, and back pressure `Pb`, the device
passes through these regimes as the forced primary flow (equivalently, primary
supply pressure) increases:

1. **`mp < mdot_choked(Pb)` -> primary UNCHOKED.** No supersonic core; the
   device is a **subsonic jet pump / mixing junction**. Primary Pt floats to
   `~Pb + dynamic head`; entrainment is small. The double-choked ejector
   construction (Huang Eqs. 1-8) does not apply.
2. **`mp` large enough to choke -> supersonic ejector.** The primary chokes at
   the throat and expands supersonically; the entrained flow's regime is then
   set by the back pressure:
   - **`Pb <= P_c*`** -> **critical** (double-choked): `omega = omega_crit`,
     independent of `Pb`. *(modeled today)*
   - **`P_c* < Pb < P_b0`** -> **subcritical**: `omega` droops from
     `omega_crit` toward 0 as `Pb` rises. *(not modeled)*
   - **`Pb >= P_b0`** -> **backflow/breakdown**: `omega <= 0`, secondary
     reverses out the suction port. *(not modeled)*

`P_b0` is the **zero-entrainment (dead-head) back pressure** -- the highest back
pressure the ejector can hold with no secondary flow. It is computable from the
*existing* Kracik & Dvorak mixing closure evaluated at `omega = 0` (their
`gam = 0`): `z3 = z1`, `lambda3` = subsonic conjugate of `lambda1`, and
`P_b0 = recovery * p_g * q(lambda1)/q(lambda3)`. No new calibration constant.

So there are **two coupled branch selections**: R0 (choked <-> unchoked
primary) *gates* R1 (critical <-> subcritical <-> backflow entrainment). The
entrainment sub-branches only exist once the primary is choked.

### 2.3 The `M_py < 1` degeneracy is a symptom of the mis-forced choked root

When the secondary stagnation pressure approaches or exceeds the primary
(`P_e >~ P_g`), the primary core is compressed back to **subsonic** at the
mixing plane y-y (`M_py < 1`). In that state Huang's constant-pressure-mixing
entrainment formula still returns a number, but a large and spurious one: the
entrained area `A_sy = A_3/A_t - A_py/A_t` becomes almost the whole mixing
chamber, so `omega` blows up. For the case in Section 2.4 this yields the
nonsensical `omega_crit = 32.8` (a 0.0051 kg/s jet nominally entraining 33x its
own mass). `M_py >= 1` (primary core still supersonic where the secondary
chokes) is the true validity condition for the double-choked construction. In a
healthy motive ejector (`P_g >> P_e`) `M_py > 1` and `P_b0 > P_c*` with a
comfortable margin; the degenerate corner is exactly the weak/reversed-pressure
case, and it coincides with the primary being unchoked -- i.e. it is the same
regime-1 failure seen from the entrainment side.

### 2.4 Worked example (the reported case)

Geometry: `A_t = 3.14e-5`, `A_p1 = 1.0e-4`, `A_3 = 8.0e-4` m^2, so
`A_p1/A_t = 3.185`, `A_3/A_t = 25.48`. Air, `gamma ~ 1.40`, `R ~ 288.1`.
Forced `mp = 0.0051 kg/s`; `Pb = 101 325 Pa`.

- R0 choked artifact: `Pt_primary = 71 552 Pa` (< Pb; cannot choke).
- Physical subsonic branch: `Pt_primary ~ 102 435 Pa` (dp = +1 110 Pa).
- Critical-mode outputs at the artifact state: `M_py = 0.66` (< 1, degenerate),
  `omega_crit = 32.8`, `P_c* = 100 024 Pa`; outlet `Pb = 101 325 Pa > P_c*`.
- Frozen-y-y dead-head collapses to `P_b0 = p_g = 71 552 Pa` (degenerate).

Both the primary-flow branch (R0) and the entrainment branch (R1) are outside
their modeled regimes. The correct model here is the **unchoked subsonic
jet-pump**: primary Pt just above Pb, small forward entrainment, outlet =
primary + secondary. Fixing R0 alone is necessary but not sufficient -- with Pt
floated up to ~102 kPa the primary core is still subsonic at y-y, so the
entrainment closure still needs its own regime handling.

## 3. Proposed design

Keep the residual **layout** (`P_jct` + N+1 rows on the `MultiPortChamberElement`
topology) and the C++ `(f, J)` discipline. Only the **content** of R0 and R1
changes, plus `verify_solution_consistent`. Recovery of a single self-consistent
root across regimes is achieved with C1 blends, not discrete re-solves, matching
the house soft-barrier idiom already used by `MPCEv2Element`
(`alpha * max(0, -e_i * mdot_i)^2`, C1 at the kink).

### 3.1 R0: choked/unchoked primary via the compressible-nozzle closure

Replace `ejector_choked_mass_flow` with the **existing** compressible-orifice
nozzle relation, which already implements isentropic nozzle flow with a *smooth
choked transition* and an analytic Jacobian:

    OrificeResult orifice_compressible_residuals_and_jacobian(
        m_dot, P_total_up, T_up, Y_up, P_static_down, Cd, area, beta);

(`include/solver_interface.h`; used by `OrificeElement(regime="compressible")`).
Map: `P_total_up = Pt_primary`, `T_up = Tt_primary`, `area = A_t`,
`P_static_down = P_py` (the common mixing-plane static pressure the entrainment
construction already computes as `P_py = P_sy = P_e/((gamma+1)/2)^(gamma/(gamma-1))`),
`Cd` from `eta_p` (or 1.0), `beta` the throat/upstream area ratio. When choked
this reduces to the current relation (independent of `P_static_down`); when
unchoked, `mp` depends on `P_py` and the root floats `Pt_primary` up to just
above the back pressure -- the physical branch of Section 2.1. **Decision: use
`P_py` (internal mixing-chamber pressure), not the outlet `Pb`, as the
downstream pressure** (user-confirmed); `Pb` is a cheaper proxy retained only as
a fallback if the `R0 <-> P_py` coupling hurts convergence.

New Jacobian couplings introduced by R0: `d R0 / d P_py` (and `P_py`'s own
dependence on `Pt_secondary`, `Tt_secondary` through the entrainment
construction), on top of the existing `d R0 / d Pt_primary`, `d R0 / d Tt_primary`.

### 3.2 R1: critical <-> subcritical entrainment via a C1 smooth-min

Model the realized entrainment as the smaller of the choke-limited and the
back-pressure-limited value:

    omega_eff = smoothmin(omega_crit, omega_sub(Pb); eps)
              = 0.5 * (omega_crit + omega_sub
                       - sqrt((omega_crit - omega_sub)^2 + eps^2))

    R1 = ms - omega_eff * mp

- Below `P_c*`, `omega_sub > omega_crit` -> choke limits -> `omega_eff = omega_crit`
  (the flat plateau, unchanged).
- Above `P_c*`, `omega_sub < omega_crit` -> back pressure limits ->
  `omega_eff = omega_sub(Pb)`, drooping to 0 at `P_b0`.
- `eps` is a small fraction of `omega_crit` (e.g. `1e-3 * omega_crit`) so the
  corner is rounded over a physically negligible width but the Jacobian is
  smooth -- the same C-infinity `sqrt`-smoothing flavor as the MPCE penalty.

**`omega_sub(Pb)` fidelity -- proposed Tier 1 (linear droop):**

    omega_sub(Pb) = omega_crit * (P_b0 - Pb) / (P_b0 - P_c*)

Uses only quantities the critical closure already produces (`omega_crit`,
`P_c*`, and `P_b0` = the same `p03` formula at `omega = 0`). Exact at both
physical anchors `(P_c*, omega_crit)` and `(P_b0, 0)`; monotone; clean analytic
Jacobian; and the near-linear subcritical droop is well documented empirically
(Huang's own characteristic, ESDU 86030, and multiple CFD studies). **Tier 2**
(a later refinement) inverts the actual mixing curve `p03(omega) = Pb/recovery`
for `omega` via a 1-D inner Newton (frozen y-y), giving the true curve rather
than the chord, with the Jacobian via the implicit function theorem. Tier 1
first; Tier 2 if validation shows the chord is too coarse.

New Jacobian coupling introduced by R1: **`d R1 / d(outlet.Pt)`** (currently R1
has no back-pressure dependence). This is the structural change that lets the
subcritical branch exist as a network unknown.

### 3.3 Unchoked regime = subsonic jet-pump closure (chosen behavior)

When the primary is unchoked (regime 1), model the device as a **subsonic jet
pump** rather than rejecting the point (user-confirmed choice). The primary Pt
floats to `~Pb + dynamic head` via the R0 compressible-nozzle branch (3.1); the
entrainment comes from a **subsonic mixing/momentum balance** (mass + momentum +
energy across the constant-area mixing section with both streams subsonic and
pressure-matched), giving a small physical `omega` instead of the spurious
double-choked value. This closure blends into the choked-ejector closure at the
choke threshold so the overall `omega(mp)` and `omega(Pb)` surfaces stay C1.
(Selection between the jet-pump and ejector entrainment closures follows the
same primary-choke smooth indicator that R0 already exposes -- to be finalized
during implementation; the guiding constraint is a single C1 residual, one
solve, no discrete regime detection.)

### 3.4 `verify_solution_consistent` changes

Today it rejects `outlet.Pt > P_c*` outright (treats every subcritical point as
unphysical). Under this design:
- Subcritical points `P_c* < Pb < P_b0` become **valid** (the model now covers
  them) -- no longer rejected.
- Retain a guard only for genuinely unmodeled states: `Pb >= P_b0` (backflow,
  if that phase is deferred) and any residual non-physical corner, with a clear
  diagnostic naming the regime rather than a bare "artifact root".
- The unchoked jet-pump regime becomes a **valid converged state**, not a
  demotion.

### 3.5 Interop: one unified residual set, one solve

The two smooth blends (R0 choke transition; R1 critical/subcritical/jet-pump
selection) make the whole element a single C1 residual system that Newton solves
in one pass across all regimes -- no outer regime detection, no two-phase
re-solve, no chatter across the boundary. This is the "critical <-> subcritical
interop" requirement, and it mirrors how `OrificeElement` already spans
choked/unchoked and how `MPCEv2Element` spans merge/branch within one C1
residual.

## 4. Validation strategy

Off-design validation data is **available in air for the subcritical branch**
(regime 2's droop) and **absent for the unchoked jet-pump branch** (regime 1),
which sets two different bars for the two branches.

### 4.1 Subcritical droop -- validate against air experiment + CFD

Two in-repo, air-fluid sources cover the plateau -> knee (`P_c*`) -> droop ->
breakdown (`P_b0`) structure directly, so no refrigerant EOS / CoolProp is
needed (unlike the R141b critical fixture):

- **Henry et al. (2007), HEFAT2007** (`../../docs/ejector/Bartosiewicz_ejector_2004.pdf`;
  filename predates the confirmed citation). Fig. 5 is the measured air-ejector
  characteristic `omega(Pb/P2^0)` at three primary pressures (4/5/6 bar), each
  showing the full plateau -> droop -> near-breakdown shape, with a reported
  CFD-vs-experiment deviation < 10%. **Digitized and committed** as
  `validation/ejector/data/henry2007_fig5.py` (cross-checked against the source
  figure; see that module's docstring for the provenance + cross-check verdict).
  This is the primary droop-*shape* target.
- **Akbarnejad & Ziabasharhagh (2025)**, *Theoretical Model for Predicting
  Performance of a Gas Ejector*, Thermal Engineering 72(1). Air runs (case E2;
  E1 temperature sweep) with numeric **tables of both critical and breakdown
  back pressures** across ~7 primary/secondary combinations -- directly
  digitizable anchor values for `P_c*` and `P_b0` without curve tracing. The
  anchor-pressure cross-check. *(Not yet digitized -- add when Phase 1 lands.)*

Use these as **validation, not calibration**: check the model's droop shape and
`P_c*`/`P_b0` anchors against them (trend + loose relative tolerance, in the
spirit of the existing critical-mode `pc_rel_tol`); do not fit
`recovery_efficiency` or any droop parameter to them. If the Tier-1 linear chord
deviates materially from the Henry Fig. 5 droop, that deviation is the concrete
trigger to implement Tier-2 (mixing-curve inversion), per sec 3.2.

*(Galindo et al. (2020) has the richest fitted subcritical + breakdown surface
but uses R1234yf -- outside combaero's air/combustion EOS, so it is a
methodology reference only, not a usable fixture.)*

### 4.2 Unchoked jet-pump branch -- validate analytically, no dataset needed

None of the in-repo ejector references (nor the surveyed literature) report a
measured `omega` characteristic for the unchoked-primary / subsonic jet-pump
regime; only Lienhard/McGovern even name it, qualitatively. This branch does
**not** require a bespoke dataset, because it rests on exact analytic relations
rather than an empirical droop:

- R0 is isentropic compressible-nozzle flow (choked/unchoked), already validated
  in the suite via the compressible `OrificeElement` path R0 reuses (sec 3.1) --
  the jet-pump branch **inherits** that validation.
- The entrainment is fixed by mass/momentum/energy conservation across the
  subsonic mixing section (sec 3.3), checkable by conservation + consistency
  assertions rather than a fitted curve.

So validate regime 1 by **analytical exactness + consistency checks**:
`dp >= 0` for forward flow, mass/energy balance, and C1 continuity of the
`omega(mp)` / `omega(Pb)` surfaces with the ejector branch at the choke
threshold. (An external eductor/jet-pump reference such as ESDU 86030 would be a
nice-to-have independent number, but is not required for the branch to be
validated.)

### 4.3 Golden-fixture + regression mechanics

- Extend the golden fixture (`validation/ejector/data/huang1999_reference.json`
  + generated `.h`) with **subcritical rows** (a `Pb`-sweep at fixed inlet
  conditions across `P_c* -> P_b0`) and **unchoked jet-pump rows**, each with
  central-difference Jacobian targets, so the C++ analytic `(f, J)` is checked
  the same way the critical rows are today.
- Regression + trend checks (monotonicity of `omega` in `Pb`; continuity at
  `P_c*` and at the choke threshold; `dp >= 0` for forward flow).

## 5. Open decisions (flagged, not yet locked)

1. **Backflow coverage.** Implement subcritical droop + unchoked jet-pump now,
   defer true backflow (`omega < 0`, sign-aware suction port) to a later phase
   behind a guard? (Recommended: yes, defer backflow.)
2. **Subcritical fidelity.** Ship Tier 1 (linear droop) first, add Tier 2
   (mixing-curve inversion) only if validation demands it? (Recommended: yes.)
3. **R0 downstream pressure resolution.** `P_py` (chosen) vs `Pb` proxy if the
   `R0 <-> P_py` coupling degrades Newton robustness -- to be confirmed
   empirically during implementation.

## 6. Phased implementation plan

Each phase carries the mandatory syncs from `CLAUDE.md` (units_data.h,
API_PYTHON.md / API_CPP.md, CHANGELOG `[Unreleased]`, `test_units_sync.py` for
new exported symbols, and the doc-hygiene rules).

1. **Reference model** (`validation/ejector/` + `_ejector_huang1999.py`): add
   `P_b0` (dead-head), `omega_sub` (Tier 1), the subsonic jet-pump entrainment,
   and the smooth-min blend as pure closed forms. Extend the golden fixture.
   - **1A (DONE):** `dead_head_back_pressure` (P_b0), `subcritical_entrainment_ratio`
     (Tier-1), `blended_entrainment_ratio` (smooth-min); golden fixture extended
     with per-case P_b0 + a subcritical `Pb`-sweep; see the provenance record in
     sec 8.1.
   - **1B (DONE):** `jet_pump_entrainment_ratio` (unchoked primary), a Keenan
     constant-area mixing solve reusing the Kracik & Dvorak closure with a
     subsonic primary; golden fixture extended with an air jet-pump sweep; see
     the provenance record in sec 8.2.
2. **Production `EjectorElement`**: new R0 (compressible-nozzle, `P_py`
   downstream) and R1 (smooth-min); updated `verify_solution_consistent`;
   analytic `(f, J)` including the new `d R0/d P_py` and `d R1/d(outlet.Pt)`
   couplings. Regression against the reported case: it must converge to the
   subsonic jet-pump root (`Pt_primary ~ 102 kPa`, `dp >= 0`, small `omega`).
3. **C++ port** (`include/ejector.h` / `src/ejector.cpp`): analytic Jacobians
   (sympy cross-check where useful), ctests consuming the extended fixture. Reuse
   the compressible-orifice C++ path for R0 where possible.
4. **pybind bindings** (separate PR), same fixture.
5. **GUI**: surface the operating regime (critical / subcritical / jet-pump /
   backflow) as an ejector diagnostic; validate the reported network converges.

## 6b. Phase 2 C++ port -- progress and the element-assembly blueprint

Phases 2 (element) and 3 (C++ port) were merged into one combined effort: the
element is C++-`(f, J)`-only, so a Python element path would be pure churn. Done
and committed on `feat/ejector-phase2-element` (WIP, not yet a PR):

- **Increment 1-2 -- C-D nozzle R0** (`ejector_cd_nozzle_mass_flow[_and_jacobian]`):
  subsonic exit-`A_p1` flux smooth-min'd with the sonic-throat cap; value parity
  vs the Python `cd_nozzle_mass_flow` (1e-9) and analytic `d/d(p0,t0,P_py)` vs
  FD. Introduced a generic `DualN<N>` (the validated 4-seed `Dual4` untouched).
- **Increment 3 -- jet-pump mixed discharge R3**
  (`ejector_jetpump_discharge_and_jacobian`): `recovery * p03`, the Kracik &
  Dvorak mixing reused at the subsonic-primary state, `DualN<6>` Jacobian
  (seeds p_g,t_g,p_e,t_e,P_py,omega), value parity + FD-checked all six seeds.
  The secondary entrained flux reuses `ejector_cd_nozzle_mass_flow`
  (area_throat = area_exit = A_s), so no new closure.

**Increment 4a (DONE):** pybind for `ejector_cd_nozzle_mass_flow_and_jacobian`
and `ejector_jetpump_discharge_and_jacobian` (Python values match the ctest to
machine precision) so the element can call them.

**Increment 4b -- element structure settled by a Jacobian-rank spike.** The
first blueprint (five rows, `P_py` exposed as a SECOND owned unknown with an
`R_py` row) was RULED OUT: a throwaway Jacobian-rank spike showed the 5-row
block is **rank-deficient in the jet-pump regime** (smallest singular value
~1e-8), and this is structural, not tunable. Any `R_py` that pins `P_py` is, at
the jet-pump solution, linearly parallel to an existing row -- pin via the
nozzle and `grad(R_py) = c * grad(R0)`; pin via the discharge and it is parallel
to `R3`; and at the jet-pump solution nothing zero-valued is independent of
R0-R3. So **`P_py` cannot be a Newton unknown**. (The earlier "A' reduces to
critical" claim was never validated -- that spike case was buggy.)

Settled structure -- **four rows, `P_py` DERIVED** (not a solver unknown); the
element keeps the base `MultiPortChamberElement` unknown/row count, so **no
solver-dispatch change is needed**:

    P_py = s_choke*P_sy + (1 - s_choke)*nozzle_inverse(mp, Pt_p)   derived
    s_choke = smoothstep(mp / choked_mass_flow(Pt_p))    keyed on mp, NOT P_py
    R0  = mp - cd_nozzle(Pt_p, Tt_p, P_py, A_t, A_p1)          primary C-D nozzle
    R1  = ms - [s_choke*omega_crit*mp + (1 - s_choke)*ms_sec]  blended entrainment
          ms_sec = cd_nozzle(Pt_s, Tt_s, P_py, A_s, A_s)        (secondary flux)
    R2  = mdot_out - mp - ms                                   mass conservation
    R3  = outlet.Pt - recovery*[s_choke*p03_crit + (1 - s_choke)*p03_jetpump(P_py)]

- **Full rank in both regimes** (spike-confirmed: 4/4 in jet-pump, 4/4 in
  critical). `R0` transitions smoothly from a real constraint (`mp = choked` in
  critical) to nearly-redundant-but-nonzero (jet-pump); the jet-pump block is
  ill-conditioned in deep-unchoked (cond ~1e4) but well within the solver's LM
  fallback.
- **`s_choke` is keyed on `mp / choked_mass_flow(Pt_p)`, NOT on `P_py`** -- this
  is essential: a `P_py`-based indicator is circular (`s` depends on `P_py`
  depends on `s`) and mislabels the choked critical state. It must be a
  SATURATING C1 form (clamped smootherstep on `mp/cap` with `lo~0.90`,
  `hi~0.999`) so `s = 1` EXACTLY at the choke (`mp -> cap`) -- a `tanh` that
  only asymptotes to 1 leaves an ~10 Pa R3 residual and a ~0.3% omega error in
  critical mode. Verified: with the smootherstep, both physical roots evaluate
  to all-four-residuals ~= 0 (jet-pump: R0=-6e-17, R3=1e-10, omega~6.6;
  critical: s=1, P_py=P_sy, omega_eff=omega_crit, R3=0 all exact).
- **`nozzle_inverse(mp, Pt_p)`** is the new closure needed: the unchoked C-D
  nozzle exit static `P_py` such that `cd_nozzle(Pt_p, P_py) = mp` -- a 1-D
  inversion with an implicit-function derivative (`dP_py/dmp`, `dP_py/dPt_p`).
  Python reference -> C++ -> FD-checked, same recipe as the C-D nozzle.
- **Spurious-root exclusion still A' + B:** R3 pins `outlet.Pt = recovery*p03`
  (blended discharge), so the choked artifact (`discharge ~= P_c* != outlet`) is
  a non-root; B is the physical warm-start (Pt_primary >= outlet) keeping Newton
  out of the sticky choked basin.
- **Jacobian VALIDATED** (prototype matches FD at non-degenerate points, both
  regimes): assembled via a Python dual (value + grad-dict over the 8 primitives
  mp, ms, mdot_out, p_g, t_g, p_e, t_e, p_out), chaining the committed closures'
  partials -- `P_py` is injected as a dual carrying its implicit partials from
  `cd_nozzle_exit_static` + `P_sy` + `s_choke`. Two physics details this pinned
  down: the discharge mixing uses `omega_eff = ms / mp` (the ACTUAL ratio; R1
  separately pins `ms = ms_model`), and the symmetric reported root has
  genuinely ~0 R3/d(mp,ms) partials (FD is unreliable there -- validate the
  Jacobian at asymmetric points).
- **`P_jct` bookkeeping:** R3 must reference `outlet.Pt` (the A' pin), so it
  cannot also define the base's `P_jct`. The ejector is not an impulse-CV
  junction, so drop `P_jct` from `unknowns()` (the 4 rows R0-R3 constrain the 3
  port pressures + mass flows directly; `flow_at_node` uses the port outer-mdot
  indices, not `P_jct`, so it is unaffected) -- verify the DOF balances via the
  end-to-end network solve.
- Then: assemble the 4-row analytic `(f, J)` (blend of the committed closures +
  the derived `P_py` chain, whose Jacobian folds in `dP_py` via the chain rule);
  `verify_solution_consistent` accepts subcritical + jet-pump roots (drop the
  `outlet <= P_c*` demotion); rewrite the element unit tests for the new
  structure; extend the golden fixture with coupled rows; and the end-to-end
  regression: `gui/tmp/ejector_test_low_flow_not_converged.json` converges to
  `Pt_primary ~= Pb`, `dp >= 0`, `omega ~= 7`.

## 6c. Phase 2 -- final implemented structure (supersedes 4b's "`P_py` DERIVED")

**What actually shipped, and why 4b's derived-`P_py` plan was wrong.** During
element assembly the 4b structure (`P_py` DERIVED from the primary nozzle
inverse) converged but ALWAYS to the spurious double-choked artifact
(`s_choke -> 1`, `omega ~ 33`), never the jet pump. Root cause: deriving
`P_py = nozzle_inverse(mp, Pt_p)` makes `R0 = mp - cd_nozzle(Pt_p, P_py)`
**VACUOUS in the unchoked branch** -- `P_py` is defined AS the exit static that
makes the nozzle flow `mp`, so `R0 == 0` identically for any `Pt_p`. R0 then
fails to pin `Pt_p`, and the solver slides to the choked basin. (Deriving `P_py`
from the SECONDARY nozzle instead just moves the vacuity to R1.) The two nozzle
relations genuinely share one `P_py`; it must be a free unknown for both to stay
live.

**Resolution -- `P_py` is the single owned unknown (NOT a second one).** The 4b
rank spike ruled out `P_py` as a SECOND unknown carrying its own `R_py` row (that
row is parallel to R0/R3 in the jet pump). But the base already owns exactly one
scalar (`{id}.P_jct`) and emits N+1 = 4 rows -- the DOF the junction contract
requires (4 rows - 3 skipped port-MCN mass rows - 1 owned = 0). So **repurpose
that single owned unknown as `P_py`** (the mixing-plane static), add NO new row,
and let R0/R1/R3 close it through the network coupling. Full rank, DOF-balanced,
no solver-dispatch change. Shipped rows (`s = s_choke`):

    R0 = mp - cd_nozzle(Pt_p, P_py, A_t, A_e)                  primary C-D nozzle
    R1 = ms - [s*omega_crit*mp + (1 - s)*ms_sec(P_py)]         blended entrainment
    R2 = mdot_out - mp - ms                                    mass conservation
    R3 = [w_pin*outlet.Pt + (1 - w_pin)*P_py] - recovery*discharge
         discharge = s*P_c* + (1 - s)*p03_jetpump(P_py, omega_eff = ms/mp)
         w_pin     = 1 - s*(1 - s_sub)
         s_sub     = smootherstep(outlet.Pt / (recovery*P_c*), 0.98, 1.02)

**`R3` is regime-dependent by `w_pin`** -- this is the key correction over 4b's
unconditional `outlet.Pt = recovery*discharge` pin, which BROKE critical mode
(it over-constrains a genuine double-choked point whose outlet is legitimately
below `P_c*`, forcing `outlet.Pt = recovery*P_c*` against the fixed downstream
boundary -> `|F|` stall). The three regimes:

  * jet pump (`s=0`) and subcritical droop (`s=1, outlet > P_c*, s_sub=1`):
    `w_pin = 1` -> pin `outlet.Pt`. The pin is what excludes the artifact
    (its recovered discharge cannot match the imposed back pressure), and the
    downstream boundary forces `P_py` through the pin, which R0/R1 turn into
    `Pt_p` and entrainment.
  * critical (`s=1, outlet < P_c*, s_sub=0`): `w_pin = 0` -> `P_py = recovery*P_c*`
    (diagnostic) and the outlet floats free. This is exactly 4b's diagnostic
    form, recovered as the `w_pin -> 0` limit, so all the existing critical-mode
    element tests pass unchanged.

**`0 * NaN` guard.** The critical closures (`entrainment_ratio`,
`critical_back_pressure`) return NaN for an unchoked primary; since smootherstep
is flat at `s in {0, 1}` the inactive branch's weight AND gradient vanish, but
`0 * NaN = NaN` would still poison the value/Jacobian. Each branch is therefore
evaluated ONLY where its weight is nonzero (`crit_active = s.v > 0`,
`jet_active = s.v < 1`); the skipped branch is a hard-zero dual, C1-consistent
with the flat endpoint. This was the concrete source of the `|F| = nan` stall
before the guard.

**Warm-start (jet-pump regime).** The GUI ejector warm-start seeds a mass-flow
primary at its choked-inverse Pt; when that is BELOW the outlet back pressure
the primary cannot choke, so it now seeds `Pt_p = outlet.Pt` (jet-pump branch)
and `P_jct = 0.99*outlet.Pt` instead of the choked-inverse / `P_c*` seed. Cold
Newton lands in the jet-pump basin directly (`graph_builder._seed_ejector_warmstart`).

**Verified:** element Jacobian matches FD at the jet-pump root (the earlier
`d(R3)/d(mc_p.Pt) = NaN` gone; the lone residual `e4.m_dot R3` "mismatch" is FD
round-off on a genuinely ~0 partial). End-to-end: the reported low-flow network
converges to `omega = 6.62`, `s_choke = 0`, `P_py = 100.1 kPa`, `Pt_p = Pb`
(`dp = 0`); the 700/22.85/40 kPa double-choked reference still converges in
critical mode (`critical_mode = 1`). Regression tests:
`test_ejector_runner_solves_unchoked_jet_pump_regime` and the retained cold
critical solve.

**Assembly moved to C++ (whole-element (f, J)).** The 4-row blend was first
assembled in Python via a forward-mode dual (`_D`) chaining the C++ scalar
closures' partials -- correct and FD-validated, but a pattern used by no other
element. A repo survey established that the house practice for a multi-row
COUPLED junction is whole-element `(f, J)` in C++ (the base
`MultiPortChamberElement` and `TeeJunctionElement` both do this; the ejector is
a `MultiPortChamberElement` subclass), so the assembly was ported to a single
C++ function `ejector_element_residuals_and_jacobian` (`src/ejector.cpp`) that
seeds a `DualN<9>` over the nine unknowns and runs the identical smootherstep /
branch-skip / `w_pin` blend, returning the 4 residuals + a 4x9 Jacobian. The
Python `residuals()` is now a thin relabeling shim (physical-flow signs +
seed->column names), mirroring the base class. The C++ output is **bit-identical**
to the retired Python dual (residuals to 0.0, Jacobian to ~1e-11 round-off,
across jet-pump / critical / non-converged points) since it is the same
chain-rule over the same closures; and it is FD-checked directly in
`tests/test_ejector_jacobians.cpp::ElementResidualsMatchCentralDifference`. The
`_and_jacobian` scalar closures and this element function are solver-internal
(`combaero._solver_tools`), a tier below the documented public API, so no
`units_data.h` / `API_*.md` entry is required (consistent with the other
ejector closures).

## 7. References

- Huang, Chang, Wang, Petrenko (1999), *A 1-D analysis of ejector performance*,
  Int. J. Refrigeration 22:354-364.
  `../../docs/ejector/Huang_1d_analysis_ejector.pdf`.
- Kracik, Dvorak (2016), *Development of an Analytical Method for Predicting
  Flow in a Supersonic Air Ejector*, EPJ Web Conf. 114:02059.
  `../../docs/ejector/Development_of_an_Analytical_Method_for.pdf`.
- Henry, Leclaire, Hemidi, Seynhaeve, Bartosiewicz (2007), *Analysis of
  Supersonic Ejector Operation: Experimental Validation and Two-Phase Aspects*,
  HEFAT2007, Sun City, paper HF1 -- air-ejector experiment + CFD; Fig. 5 is the
  subcritical-droop validation target. `../../docs/ejector/Bartosiewicz_ejector_2004.pdf`
  (filename predates the confirmed citation). Digitized:
  `validation/ejector/data/henry2007_fig5.py`.
- Akbarnejad, Ziabasharhagh (2025), *Theoretical Model for Predicting
  Performance of a Gas Ejector*, Thermal Engineering 72(1) -- air runs with
  numeric critical + breakdown back-pressure tables (`P_c*`/`P_b0` anchors).
- Galindo et al. (2020), R1234yf ejector (fitted subcritical + breakdown
  surface) -- methodology reference only; refrigerant fluid, not an air fixture.
- Sanger, N.L. (1968), *Noncavitating Performance of Two Low-Area-Ratio Water
  Jet Pumps Having Throat Lengths of 7.25 Diameters*, NASA TN D-4445 -- the
  lossless incompressible jet-pump relation N(M,R) (Appendix B) that the Phase
  1B compressible closure reduces to as Mach -> 0, and the measured water data
  that anchors the ~1.4x lossless overprediction (recovery ~= 0.7).
  `../../docs/ejector/19680008095.pdf`.
- ESDU 86030, *Ejectors and jet pumps* (subcritical/breakdown characteristic
  shape reference; optional independent jet-pump cross-check).

## 8. Implementation provenance records

Short WHY-only records closing each implemented phase (per the
`model-provenance` skill): formulation chosen + alternatives, literature basis,
the invariant it must preserve, dead ends, and the test that pins it.

### 8.1 Phase 1A -- subcritical droop (P_b0, Tier-1 omega_sub, smooth-min)

**Formulation & alternatives.** P_b0 is Kracik & Dvorak's *existing* mixing
closure evaluated at `omega = 0` (their `gam = 0`), not a new correlation -- so
the subcritical window carries no new calibration constant, only the `omega`
anchor the critical closure already produced. `critical_back_pressure` was
refactored to share a single `_yy_lambda_state` + `_mixed_flow_stagnation`
helper with `dead_head_back_pressure` so the two anchors cannot drift apart;
the refactor is behaviour-preserving (the 31 golden values regenerate
byte-identically). omega_sub is the Tier-1 linear chord between the two anchors;
Tier-2 (mixing-curve inversion) was deliberately deferred, not implemented, per
the "Tier 1 first" open decision (sec 5.2). The critical<->subcritical join is a
C1 smooth-min (0.5*(a+b-sqrt((a-b)^2+eps^2))), the same sqrt-smoothing idiom as
`MPCEv2Element`'s soft barrier; a bare `min()` was rejected because its
derivative jump would break the single-Newton-solve requirement (sec 3.5).

**Literature basis.** Mixing closure: Kracik & Dvorak 2016 (Eqs. 7-13). The
near-linear-droop assumption is cross-checked against measured AIR data, Henry
et al. 2007 HEFAT2007 Fig. 5 (`henry2007_fig5.py`): a least-squares line through
each curve's droop fits with R^2 0.992 / 0.989 / 0.939 (4/5/6 bar). Absolute
anchor placement (P_c*, P_b0) is NOT validated here -- that needs Akbarnejad &
Ziabasharhagh 2025's numeric tables (not yet digitized); Henry validates the
*shape* only.

**Invariant preserved.** P_e < P_c* < P_b0 < P_g on all 31 rows; the omega -> 0
limit is exact (at A_3/A_t = A_py/A_t, `critical_back_pressure` == P_b0); the
smooth-min stays within eps/2 of `min` and is C1 through P_c*; the pre-existing
critical golden fixture is unchanged.

**Dead ends.** First measured linearity via deviation from a chord anchored at an
*eyeballed* knee -- discarded: it conflated "is the droop linear" (yes, R^2 >=
0.94) with "is my knee estimate right", and flagged a false 26% error on the 6
bar curve. The least-squares metric isolates linearity. The residual concavity
it does show (6 bar, up to 8.5% of plateau) is the documented trigger for Tier 2.

**Pinned by.** `python/tests/test_ejector_subcritical.py` (anchors, omega -> 0
continuity, smooth-min bound + C1, Henry linearity, fixture regression) and the
extended `huang1999_reference.json` / `_data.h` golden data.

### 8.2 Phase 1B -- unchoked-primary subsonic jet-pump

**Formulation & alternatives.** `jet_pump_entrainment_ratio` is a Keenan-Neumann-
Lustwerk constant-area mixing solve: both streams expand isentropically to a
common mixing static pressure p1 (primary through A_p1, entrained through
A_3 - A_p1) so omega = ms/mp = f(p1), and the physical p1 is the one whose mixed
flow recovers to the back pressure (a 1-D root solve). Two grounded candidates
were tried in the spike and REJECTED with numbers: (a) geometric area-filling
(as in the choked closure) reproduces the omega ~ 33 degeneracy -- a weak primary
cannot really entrain a whole annulus; (b) a momentum balance that leaves omega
free is under-constrained. The mass+momentum+energy mixing is NOT a new closure:
it was verified bit-identical to Kracik & Dvorak's `_mixed_flow_stagnation`
(0.00 Pa difference), which is reused here with a subsonic lambda1 -- so the
jet-pump and ejector branches share one mixing model and meet at the
primary-choke threshold.

**Literature basis.** Ground truth is the lossless incompressible jet-pump
relation N(M, R) (Sanger, NASA TN D-4445, Appendix B), the same conservation
laws in the constant-density limit. The compressible closure is verified to
converge to it as Mach -> 0 (relative error 1.1% -> 0.034% as the driving
pressure falls), the exact-limit signature.

**Constant audit (recovery_efficiency).** Kept at 1.0, no fitted constant, per
the calibration policy -- but AUDITED against data: the lossless form
overpredicts Sanger's measured water-jet-pump head by ~1.4x (real pump ~= 70% of
ideal head, consistent across his R = 0.066 and 0.197 pumps), so recovery ~= 0.7
is the data-anchored value a caller would set. Documented, not baked in. (A
first audit pass reported ~2x / recovery ~0.5; that came from a mis-transcribed
Eq. 14 off the 1968 scan -- the Mach -> 0 reduction check caught the unphysical
"N rising with M", and the equation was re-derived from first principles. The
provenance discipline caught the error before it reached code or tests.)

**Invariant preserved.** Bounded, forward omega O(few) on the reported case
(NOT 33); omega monotone increasing as the back pressure falls; named boundary
handoffs (`back_pressure` -> omega 0; `primary_choke` -> the ejector branch);
refusal when the primary is too strong to be unchoked.

**Dead ends.** Area-filling (omega ~ 33) and free-omega momentum balance (both
above). Also an inverted interior-bisection direction, caught by the low-Mach
reduction test landing 50% off before the fix.

**Pinned by.** `python/tests/test_ejector_jetpump.py` (low-Mach reduction +
monotone convergence, ~1.4x Sanger overprediction, bounded reported-case omega,
monotonicity, boundary flags, fixture regression) and the `jetpump_cases`
section of the golden fixture.

### 8.3 Phase 2 -- element integration (all three regimes in one solve)

**Formulation & the two rejected structures.** The `EjectorElement` residuals
were rebuilt to resolve critical / subcritical-droop / unchoked-jet-pump in one
C1 Newton system (full structure in sec 6c). Two grounded structures were built
and REJECTED with evidence: (a) the 4b "`P_py` DERIVED from the primary nozzle
inverse" plan converged but ALWAYS to the double-choked artifact (`omega ~ 33`),
because deriving `P_py` makes `R0` algebraically vacuous in the unchoked branch,
so it stops pinning `Pt_p`; (b) an unconditional outlet-pin `R3` (right for the
jet pump) BROKE genuine critical mode -- it forces `outlet.Pt = recovery*P_c*`
against a fixed downstream boundary whose outlet is legitimately below `P_c*`,
stalling at `|F| ~ 1.8e4`. The shipped structure owns `P_py` as the SINGLE
repurposed junction unknown (no new row, DOF unchanged) and makes `R3`
regime-dependent via `w_pin = 1 - s_choke*(1 - s_sub)`, recovering 4b's
diagnostic form as the `w_pin -> 0` (critical) limit.

**Why the artifact is a genuine math root and how it is excluded.** With a
back-pressure-independent closure the choked artifact (`Pt_p < Pb`, primary
sonic, `omega ~ 33`) satisfies R0-R2 and is a real root; only the OUTLET PIN
excludes it -- its recovered discharge cannot match the imposed back pressure,
so it is a non-root of the FULL network. The pin is active exactly where it must
be (jet pump + droop) and off where it must be (critical).

**Constant audit.** No new fitted constants. The only new numeric knobs are
smootherstep windows: `s_choke in [0.90, 0.999]` on `mp/cap` (must SATURATE to
exactly 1 at choke -- a non-saturating form leaves a ~10 Pa R3 residual and
~0.3% omega error in critical, per 4b), and `s_sub in [0.98, 1.02]` on
`outlet.Pt/(recovery*P_c*)` (centred on the physical critical/subcritical
boundary `outlet = P_c*`). `recovery_efficiency` stays 1.0 (audited in 8.2).

**Invariant preserved.** The pre-existing critical-mode element tests
(`residuals` cross-check vs the reference physics, all-residuals-zero at the
choked root, analytic-Jacobian-vs-FD, full double-choked network) pass
UNCHANGED -- the critical regime is the `s_choke=1, s_sub=0` limit of the new
blend. Mass conservation (`R2`) and `omega = ms/mp` (actual, regime-independent)
hold at every root.

**Dead ends.** Derived-`P_py` (vacuous R0) and unconditional outlet-pin (breaks
critical), both above. Also a `0 * NaN = NaN` Jacobian poison from eagerly
evaluating the critical closures (which return NaN for an unchoked primary) even
at `s_choke = 0` -- fixed by evaluating each branch only where its weight is
nonzero, C1-safe because smootherstep is flat at the endpoints. Also a warm-start
that seeded the primary at its choked-inverse Pt (the artifact basin) for a
mass-flow primary below the choke threshold -- fixed by detecting that case and
seeding `Pt_p = outlet.Pt`.

**Pinned by.** `python/tests/test_ejector_element.py` (new
`test_diagnostics_reports_actual_omega_and_regime`,
`test_diagnostics_reports_jet_pump_regime`,
`test_verify_solution_consistent_always_true_all_regimes_modeled`, plus the
retained critical residual/Jacobian/network tests) and
`python/tests/test_gui_ejector.py::test_ejector_runner_solves_unchoked_jet_pump_regime`
(the reported low-flow case -> `omega ~ 6.6`, `s_choke = 0`, `dp >= 0`)
alongside the retained cold critical solve. The C++ whole-element assembly is
FD-guarded by
`tests/test_ejector_jacobians.cpp::ElementResidualsMatchCentralDifference`.
