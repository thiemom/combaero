# Ejector operating-regime extension: choked/unchoked primary, subcritical entrainment, and unchoked jet-pump mode

Design + findings note for extending `EjectorElement` beyond the critical
(double-choked) plateau. Written 2026-08-26. Audience: developers of the
network module (not end users). Captures the diagnosis of a real
non-converging low-primary-flow case and proposes the model extension that
resolves it. Nothing here is implemented yet -- this is the spec that the
implementation PRs will follow.

Companion docs:
- `README.md` (this directory) -- current (critical-mode) model, accuracy
  comparison, golden-fixture contract.
- `../../python/combaero/network/_ejector_huang1999.py` -- the closed-form
  critical-mode physics (Huang 1999 + Kracik & Dvorak 2016).
- `../../docs/ejector/handover.md` -- the original (2026-07-13) groundwork
  handover (untracked, alongside the source PDFs; now largely superseded, kept
  for history).

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

- Extend the golden fixture (`validation/ejector/data/huang1999_reference.json`
  + generated `.h`) with **subcritical rows** (a `Pb`-sweep at fixed inlet
  conditions across `P_c* -> P_b0`) and **unchoked jet-pump rows**, each with
  central-difference Jacobian targets, so the C++ analytic `(f, J)` is checked
  the same way the critical rows are today.
- Cross-check the subcritical droop shape and `P_b0` against the air-ejector
  CFD + experiment in `../../docs/ejector/Bartosiewicz_ejector_2004.pdf` and the
  characteristic curves in Huang (1999).
- Regression + trend checks (monotonicity of `omega` in `Pb`; continuity at
  `P_c*` and at the choke threshold; `dp >= 0` for forward flow).
- Keep the calibration policy: anchor closures on analytical/CFD references, do
  not fit `recovery_efficiency` or any droop parameter to a measurement set.

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

## 7. References

- Huang, Chang, Wang, Petrenko (1999), *A 1-D analysis of ejector performance*,
  Int. J. Refrigeration 22:354-364.
  `../../docs/ejector/Huang_1d_analysis_ejector.pdf`.
- Kracik, Dvorak (2016), *Development of an Analytical Method for Predicting
  Flow in a Supersonic Air Ejector*, EPJ Web Conf. 114:02059.
  `../../docs/ejector/Development_of_an_Analytical_Method_for.pdf`.
- Bartosiewicz et al. (2004), air-ejector CFD + experiment (validation
  cross-check). `../../docs/ejector/Bartosiewicz_ejector_2004.pdf`.
- ESDU 86030, *Ejectors and jet pumps* (subcritical/breakdown characteristic
  shape reference).
