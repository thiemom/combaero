# Inner Iterative Solves in the C++ Kernels: Audit Record (issue #868)

Follow-up to #861, which found two Python inner solves returning a non-root
on their iteration cap (`_flow_consistent_residuals`, fixed in #863;
`_self_consistent_residuals`, fixed in #867). This audits every iterative
inner solve in `src/` and the header-only code in `include/` (none of the
headers holds one: their loops run over ports, species or quadrature nodes).
The code owns the solves; this owns what was measured, what was changed and
what was left for a decision.

## Method

1. **Inventory.** Every loop with an iteration cap, a convergence test or a
   bracket search in `src/*.cpp` and `include/*.h` (grep for `max_iter`,
   `for (...; it < N; ...)`, `while`, Newton, bisection, Brent, secant,
   Muller, fixed point). Loops over species, ports, panels, layers or a fixed
   RK4 grid are not solves and are not listed.
2. **Instrumentation** (a throwaway patch, never committed): each cap exit,
   flat-derivative exit and unbracketed Brent call appended its solve name,
   last residual, last step or bracket width and the inputs to a log, keyed
   by the test that was running. Run on the full Python suite (5926 tests),
   the Cantera validation (678) and the network F and J of all 14 stress
   families x seeds 0-39 at x0 plus three perturbed points (2240 evaluations).
3. **Kinks and clamps:** every cap hit was traced to its inputs and the
   failure reproduced directly (the falsifying tests below).
4. **Fixes only for clear defects** (a silent non-root on the cap). Every fix
   leaves converged results bit-identical: the same 2240 network F and J
   evaluations and 13 kernel sweeps hash identically on main and on the branch.

Cap hits, Python suite on main (Cantera validation and stress x0 in brackets):

| Solve | Hits | Where | What the exit returned |
|---|---|---|---|
| `solve_T_isentropic` | 148 849 (0) | `test_solver_balance_552` (2 tests), `test_energy_boundary_target_792` (2) | all at `T_FLOOR` = 10 K: the root is below the floor (P/Pt ~ 1e-5 and less) and every step is clamped back to it, 50 times per call |
| `newton_damped` (equilibrium) | 787 (+42 Cantera) | WGS and equilibrium tests, `test_fuel_boundary_892` | see row E1: 236 at an empty interval, ~546 stopped at ln(Q/K) = 13-14, 5 adiabatic T off by 1729 J/kg |
| `fanno_inlet_mach` | 995 (0) | `test_paper_examples::test_c3x_panel_numbers` (992), `test_out_of_range_481` (2), `test_ejector_reverse_flow_619` (1) | all at G < 1.4e-11 G_max, relative flux residual up to 0.69 |
| `fanno_mach` `brent`, unbracketed | 55 (0) | `test_c3x_panel_numbers` | a point of an interval with no sign change, as the choked flux |
| `chamber_static_state` | 48 (0) | `test_tee_network` cascaded manifold (2 xfail tests) | probes at Tt = -1.4e19 K: T = 1e11 K or -1.4e19 K |
| `nozzle_quasi1d` area-Mach | 57 (0) | `test_compressible::test_de_laval_nozzle_subsonic` | A/A*_eff 0.81-0.99999 < 1: no subsonic root; held at the 0.999 clamp |
| `jobson_expansibility` | 2 (0) | Jobson tests | last step 3e-14 / 4e-15 relative: a rounding cycle, a root |
| `solve_P_back_from_mdot`, `solve_Pt_from_mdot` | 1 + 1 (0) | | bracket at 1 ulp (5.8e-11 Pa at 4e5 Pa): a root |

Every other solve in the table below never reached its cap.

## The table

Verdicts: **safe** (the cap exit is a root, or the exit is reported),
**latent** (silent non-root on the cap, not reached in the suite),
**bug** (silent non-root reached). Actions: **fixed** in this PR,
**rec** = recommendation left for a decision (next section), **none**.
"Throw" means `std::runtime_error` naming the function and the last residual
(`include/inner_solve.h`). "Cap rule" means a bracketed solve's cap exit is
returned when its bracket or last step has closed to
`inner_solve::kCapAcceptRel` = 1e-10 of the iterate (a root to rounding, which
a 1e-13..1e-15 tolerance can cycle on), else throws.

Line numbers are on the branch.

| # | file:line | function | method | cap | tolerance | exit on cap (main) | bracket / safeguard | verdict | action |
|---|---|---|---|---|---|---|---|---|---|
| T1 | thermo.cpp:388 | `solve_increasing_T` (`calc_T_from_h/s/u/sv`) | rtsafe Newton | 50 (caller) | 1e-6 K abs step, or bracket < tol | throws | factor-4 search, then bisection when Newton leaves or fails to halve | safe | none |
| T2 | thermo.cpp:751 | `calc_T_from_cp` | Newton, step <= 0.1 T, clamped [200, 6000] K | 50 (caller) | 1e-6 abs on step or cp | silent: returns the clamp for an unreachable target (6000 K for cp = 1000 J/(mol K)) | none; cp(T) need not be monotone | bug (reachable by the caller) | fixed: throws. rec R6 fixed (#939, PR #954): the root nearest T_guess in [200, 6000] K, outward scan + bracketed Newton, no root in range throws; range extended to [1, 6000] K (PR #961, below) |
| S1 | state.cpp:118 | `State::set_VH_mass` | Newton on u + R T = h | 20 | 1e-4 K abs step | silent | none (monotone) | latent | fixed: throws. rec R5 fixed (#938, PR #954): `calc_T_from_h_mass` |
| G1 | stagnation.cpp:77 | `Tt_from_static_v` (and `Tt_from_static`, `Pt_from_static`, `stagnation_from_static`) | Newton, step <= 0.5 T, floor `T_FLOOR` | 50 (caller) | 1e-8 (1 + T) step | silent (also on a flat derivative) | none (h monotone) | latent (reachable by the caller) | fixed: throws; a root below the floor still returns the floor. R1 fixed (#935, PR #949): no floor, the root is returned |
| G2 | stagnation.cpp:132 | `T_from_stagnation` (and `P_from_stagnation`) | Newton with approximate dF/dT (gamma' dropped), floor | 50 (caller) | 1e-8 (1 + T) step | silent | none | latent | fixed: throws, floor kept. rec R3. R1 fixed (#935, PR #949): no floor. R3 fixed (#936, PR #954): exact slope, quadratic |
| G3 | stagnation.cpp:324 | `chamber_static_state` | Newton from Tt | 100 | 1e-13 Tt step | silent; implicit J then evaluated at the non-root | none (f falls in T, one root) | bug (48 probe hits at Tt = -1.4e19 K) | fixed: cap rule; R4 fixed (#937): Tt <= 0 or non-finite refused |
| G4 | stagnation.cpp:408 | `sonic_static_state` | Newton | 100 | 1e-13 Tt step | silent, implicit J at non-root | none | latent | fixed: cap rule |
| G5 | stagnation.cpp:489 | `duct_static_state` | Newton in ln P on [ln P*, ln Pt] | 100 | 1e-15 abs step in ln P | silent | bisection when Newton leaves the bracket | safe | R9 fixed (#941): cap rule |
| C1 | compressible.cpp:41 | `solve_T_isentropic` (`nozzle_flow`, `mass_flux_isentropic`, `mach_from_pressure_ratio`, `critical_pressure_ratio`, the solves below) | Newton, step <= 0.5 T, floor | 50 (caller) | 1e-8 (1 + T) step | silent; 148 849 hits, all at the floor | none (s monotone) | latent (floor hits are the documented clamp) | fixed: a root below the floor returns the floor at once (same value; the 49 clamped repeats are skipped), any other cap exit throws. rec R1. R1 fixed (#935, PR #949): no floor; below 10 K the step is in ln T |
| C2 | compressible.cpp:120 | `find_critical_pressure` | golden section on [0.3, 0.99] Pt | >= 200 | width < tol (1 + a + b) | silent, bracket 0.618^200 | fixed bracket (P*/Pt 0.49-0.61 for gamma 1.67-1) | safe | none |
| C3 | compressible.cpp:292 | `solve_P_back_from_mdot` | bisection on [P*, Pt] | >= 200 | G within tol, width tol (Pt - P) | silent; 1 hit at 1 ulp | bracket | safe | none |
| C4 | compressible.cpp:341 | `solve_Pt_from_mdot` bracket | x10 expansion | 12 | | throws if not bracketed | | safe | none |
| C5 | compressible.cpp:358 | `solve_Pt_from_mdot` | bisection | >= 200 | as C3 | silent; 1 hit at 1 ulp | bracket | safe | none |
| C6 | compressible.cpp:462 | `solve_T_from_energy` (Fanno marches) | Newton, step <= 0.3 T, floor | 50 | 1e-8 T step | silent | none (F rises in T: u grows with T at fixed P) | latent | fixed: throws, floor kept. R1 fixed (#935, PR #949): no floor (no root for h0 < h(0 K): throws) |
| C7 | compressible.cpp:647 | `fanno_channel` march | adaptive RK4 | `kFannoMaxSteps` = 2e5 | | silent: outlet at x < L reported unchoked | step halving near M = 1 | latent | fixed: throws |
| C8 | compressible.cpp:894 | `fanno_channel_rough` march | adaptive RK4 | 2e5 | | silent, as C7 | | latent | fixed: throws |
| C9 | compressible.cpp:1094 | `solve_station_state` | Newton | caller | tol T step | returns `false` | | safe (reports); `[[maybe_unused]]`, no caller | R10 fixed (#942): deleted |
| C10 | compressible.cpp:1320 | `nozzle_quasi1d` area-Mach | Newton, M clamped to [0.001, 0.999] subsonic, <= 10 | 50 | abs 1e-10 in A/A* | silent; 57 hits where A < A*_eff (no subsonic root; at a choked throat the root is M = 1, held at 0.999) | clamps | best effort by design | rec R7 |
| C11 | compressible.cpp:1372 | `nozzle_quasi1d` (T, P) per station | sequential Newton, approximate da/dT, clamps | 50 (caller) | tol relative on s and h | silent | clamps T <= Tt, P in [100, Pt] | latent | rec R7 |
| F1 | fanno_mach.cpp:64 | `static_T` | Newton in a bracket | 100 | 1e-14 T step | silent | bisection when Newton leaves | safe | fixed: cap rule |
| F2 | fanno_mach.cpp:135 | `brent` | Brent | 200 | caller rel_tol (1e-13, 1e-14, `kFannoMachRelTol`) | silent at cap (unreachable with a bracket); NO check for a sign change | Brent's bisection fallback | bug: called unbracketed 55 times (F7) | fixed: throws without a sign change |
| F3 | fanno_mach.cpp:323 | `fanno_inlet_mach` | Newton, central FD slope h = 1e-7 max(M, 1e-3), in [0, 1] | 60 | 1e-15 M step | silent; 995 hits at G < 1.4e-11 G_max, M up to 73% high | bisection when Newton leaves | bug: below M = h the stencil crossed M = 0, where the flux is held at 0, read half the slope, and Newton 2-cycled | fixed: one-sided stencil when M - h <= 0 (only those iterations change), cap rule. rec R8 fixed (#940, PR #954): exact slope from `isentropic_flux_slopes` |
| F4 | fanno_mach.cpp:384 | `fanno_duct` exit Mach | `brent` on [M_in, 1] | | `kFannoMachRelTol` | | f(M_in) = -L < 0 < L* - L by construction | safe | none |
| F5 | fanno_mach.cpp:410 | `fanno_choked_mass_flux` | 60 e-fold search, then `brent` | 30 | 1e-14 | unbracketed `brent` possible | L* -> inf as G -> 0 | safe (and F2 now throws) | none |
| F6 | fanno_mach.cpp:610, 636 | `raw_flow` exit Mach | widening from the guess, halving to 1e-6, `brent` | | 1e-13 | unbracketed `brent` possible | psi > 0 at low Me | safe (F2 guards) | none |
| F7 | fanno_mach.cpp:994 | `choked_flux` (Rayleigh) | 60 e-fold search, then `brent` | 30 | 1e-14 | `brent` ran without a sign change at unphysical exit Tt probes (2.1e5 K, 1.8e12 K) and returned G = 6.3e-5, choked | | bug | fixed via F2 |
| F8 | fanno_mach.cpp:1057, 1099, 1108 | `rayleigh_raw_flow` warm and cold searches | widening / scan, `brent` only when a sign change was found | | 1e-13 | | checked | safe | none |
| A1 | area_change.cpp:581, 646 | `solve_param` / `raw_flow` bracket searches | x0.1 down, x4 up | 2000 / 200 | | silent if no bracket | loss -> 0 and -> inf | safe | none |
| A2 | area_change.cpp:589 | `solve_param` | Newton in ln p in a bracket | 200 | 1e-15 step or width | silent | bisection | safe | fixed: cap rule |
| A3 | area_change.cpp:772 | `area_change_loss` Mach | Newton in [0, 1] | 200 | 1e-15 M step | silent | bisection | safe | fixed: cap rule |
| O1 | orifice.cpp:1307 | `pipe_equivalent_c` | Newton in [1e-6, 1e3] | 100 | 1e-14 C step | silent | bisection | safe | fixed: cap rule |
| O2 | orifice.cpp:1462 | `solve_orifice_mdot` | fixed point on Cd(Re) | caller | rel change < tol | throws | | safe | none |
| O3 | orifice.cpp:1778 | `jobson_expansibility` | fixed point (approach term O(beta^4)) | 200 | 1e-15 relative | silent; 2 hits at 3e-14 / 4e-15 (rounding cycle) | | safe at the hits | fixed: cap rule (the hits are kept, bit-identical). rec R11 fixed (#943, PR #954): 1e-13 |
| H1 | humidair.cpp:205 | `dewpoint` | Newton, FD slope, clamped [T_lo, T] | 100 | 1e-6 K | throws; RH > 1 would put the root above T (Newton on the clamp, "converged" to T), but `vapor_pressure` refuses RH > 1 first | clamp | safe | R12 (#944): documented |
| H2 | humidair.cpp:317 | `wet_bulb_temperature` | Newton in a bracket | 200 | 1e-6 K | throws | bisection | safe | none |
| R1 | friction.cpp:92 | `friction_colebrook` | Newton from Serghides | 20 (caller) | 1e-10 relative | silent (also on a flat slope) | none (monotone) | latent (reachable by the caller) | fixed: throws |
| E1 | equilibrium.cpp:80 | `newton_damped` (WGS extent at 165, 276, 441, 845; adiabatic T at 217, 972) | damped Newton, FD slope eps = 1e-6 max(1, abs x), clamped | 40 / 50 | abs 1e-12 on f, or xtol_rel | silent on cap and on a flat slope | clamps only | bug (787 + 42 hits) | fixed #950 (R2): rtsafe on ln Q with the analytic slope; adiabatic T bracketed + Brent; cap rule |
| E2 | equilibrium.cpp:386 | `solve_reforming_isothermal` | Newton in a bracket | 200 | abs 1e-12 on F, width 1e-15 span | silent; a non-finite F breaks to a silent clamp | bisection | latent | fixed #950 (R9): +-inf at a bracket end bisects, NaN and cap throw |
| E3 | equilibrium.cpp:605, 672 | `joint_newton` + line search | scaled Newton, curved step, backtracking | `REFORMING_NEWTON_*_IT` / 60 | `REFORMING_NEWTON_FTOL` | returns false -> `EquilibriumResult.converged` | fraction to the boundary | safe (reports) | rec R2 (one caller drops the flag) |
| E4 | equilibrium.cpp:~750 | `element_potential_start` | Gordon-McBride Newton | `REFORMING_NEWTON_MAX_IT` | 1e-9 | returns false | RP-1311 damping | safe | none |
| B1 | combustion.cpp:1104 | `bisection_solve` (fuel/oxidizer stream inverses) | Brent | 100 (caller) | tol_f, 1e-10 rel width | silent (returns b); unreachable with its checked bracket | throws if not bracketed | safe | none |
| M1 | acoustics.cpp:1041 | `refine_root_muller` | Muller | `kMullerMaxIter` = 50 | 1e-6 rad/s | returns last iterate by design | | safe: both callers (can_annular_solvers.cpp:188, 316) accept a mode only if abs(D) < `kDispersionTolerance` y_scale | none |
| M2 | can_annular_solvers.cpp:432 | `refine_root_muller_annular` | Muller | 50 | 1e-6 | as M1 | | safe: caller (549) checks abs(D) | none |
| X1 | crossflow_cd.cpp:144 | `florschuetz_metzger_su_1982_cd` falling branch | Newton, clamped (0, 4 cd0] | 100 | 1e-15 cd0 step | silent; implicit dCd/dq then at a non-root | none (one root, dF/dCd > 0) | latent | fixed: cap rule |
| J1 | ejector.cpp:25 | `ejector_mach_from_area_ratio_supersonic` | bisection on [1 + 1e-12, 50] | 200 | | silent end point for an area ratio outside [1, A/A*(50)] | | latent | R13 fixed (#945): flag and clamp |
| J2 | ejector.cpp:168 | `ejector_cd_nozzle_exit_static_and_jacobian` | bisection on [P*, Pt (1 - 1e-12)] | 200 | | silent end point (and implicit J there) for m_dot outside the nozzle's range | | latent | R13 fixed (#945): flag and clamp, J of the clamp |
| M3 | side_stream_momentum.cpp:172 | `chamber_merge_face_state` | fixed point (quadratic in P_f, then static T) | 100 | 1e-12 Tt step | silent; implicit (P_f, T_f) derivatives at a non-root | discriminant floor | latent | fixed: cap rule |

Not solves, checked and excluded: the fixed-grid RK4 `march_zeta`
(`kRayleighSteps`), quadrature (`composite_gl`, `gl8_split`), the pole
enumeration in `can_annular_solvers.cpp:231`, the reforming start passes
(`REFORMING_START_PASSES`, a warm start whose result Newton then solves), the
table interpolation searches in `area_change.cpp` and `orifice.cpp`, and every
loop in `include/` (ports, species, dual-number components).

## Per-question findings

**1. Exit on the cap.** 47 rows (52 loops). 20 solves returned their last
iterate unchecked on the cap or on a flat derivative, or (F2) ran without a
sign change: T2, S1, G1-G4, C1, C6-C8, F1-F3, A2, A3, O1, O3, R1, X1, M3.
Five of them are bracketed with bisection (F1, A2, A3, O1, O3), so a cap exit
there sits in a bracket that is small in practice; the other 15 had no such
guarantee. Reached in the suite with a non-root: T2 (by its caller's target),
G3, F2/F7, F3. All 20 now throw or apply the cap rule. Left silent by
decision: C10/C11, E1/E2, J1/J2, H1 and the `T_FLOOR` return (R1, since
fixed: #935).

**2. Brackets and safeguards.** The bracketed solves (T1, G5, F1, F3, A2, A3,
O1, H2, E2) fall back to bisection when Newton leaves the bracket; only T1
also bisects when a step fails to halve (rtsafe), which is what breaks a
2-cycle that stays inside the bracket. F3's 2-cycle was of that kind: the
fix is the stencil, not the safeguard. The unbracketed Newtons (T2, S1, G1-G4,
C1, C6, R1, X1) are on monotone functions except T2 (cp(T)).

**3. Tolerance against the outer Jacobian.** The kernels the network
differentiates or calls per evaluation converge to 1e-13..1e-15 (G3-G5, F1,
F3, A2, A3, X1, M3) or exactly (bisection to 1 ulp), far below the 1e-6
relative FD steps the network uses (e.g. `fanno_mach.cpp` choked-flux
derivatives, the `d_ddir` 1e-5 composition step). The default 1e-8 relative
step test of G1, G2 and C1 is fine for quadratic Newton (the error after a
1e-8 step is ~1e-16), but G2 converges only linearly (its dF/dT drops the
dgamma/dT term), so its error is ~1e-8 T times the derivative ratio: 1e-2
relative noise in a derivative differenced at a 1e-6 step. G2 is called from
`solver_interface.cpp:666` (an entropy diagnostic, integrated, not
differenced). See R3.

**4. Derivative consistency at a non-converged exit.** G3, G4, F3, F7, X1
and M3 feed implicit-function derivatives; at their former cap exits those
derivatives described a non-root. With the throws they never do. J2 still
returns implicit J at a clamped end point (R13). E1's composition feeds the
equilibrium derivatives at a non-root (R2).

**5. Multiple roots and branch selection.** G3, G4, C1, C6, R1, X1, M3, H1,
H2 have one root. F1 and F3 select the subsonic root by their bracket; J1 the
supersonic one by its bracket [1, 50]; C10 by its start (0.3 or 2.0) and its
clamps. T2 is the one solve whose root depends on the start: cp(T) of a
mixture is not monotone everywhere, so `T_guess` picks the branch (R6).

## Fix evidence

- **Falsifying tests.** `tests/test_inner_solve_cap_868.cpp` (ctest) and
  `python/tests/test_inner_solve_cap_868.py`: built against main, 5 of the 6
  C++ cases and all 5 Python cases fail; on the branch all pass. The sixth C++
  case pins that the `T_FLOOR` return is unchanged (passes on both). Cases:
  `chamber_static_state` at Tt = -1.393e19 K (main: T = 1e11 K);
  `calc_T_from_cp(1000)` (main: 6000 K); `fanno_inlet_mach` at 1e-12, 1e-20,
  1e-25 of the sonic flux (main: 36-73 % high); `fanno_rayleigh_channel_flow`
  with Tt_exit = 2.1e5 K, the inputs of a C3X probe (main: G = 6.3e-5,
  choked); `friction_colebrook`, `Tt_from_static`, `T_from_stagnation` with
  `max_iter = 1` (main: the half-converged iterate).
- **Bit identity.** Network F and dense J (SHA-256) at x0 and three seeded
  2 % perturbations for every stress family x seeds 0-39, 2240 evaluations,
  0 exceptions, all identical to main; sweeps of `Tt_from_static`,
  `T_from_stagnation`, `Pt_from_static`, `P_from_stagnation`,
  `Tt_from_static_v`, `stagnation_from_static`, `mach_from_pressure_ratio`,
  `mass_flux_isentropic`, `critical_pressure_ratio`, `nozzle_flow`,
  `friction_colebrook`, `calc_T_from_cp`, `calc_T_from_h` (T 150-2500 K, M
  0.01-2.5) identical.
- **Where results move.** Only where a solve used to end at a non-root: the
  C3X example's tuned run (149) lands on the same root within 1e-11 relative
  (wall temperatures within 1.4e-11 K, |F| 1.95e-13 vs 1.78e-13), its probes
  at unphysical states now rejected as steps; run 158 is bit-identical. The
  suite, ctest and the Cantera validation pass unchanged.

## Recommendations left for a decision

- **R1 (C1, G1, G2, C6) The `T_FLOOR` return is not a root.** FIXED (#935,
  PR #949, decision 2026-11-02 "can go to zero"; section below). An expansion
  to P/Pt below ~1e-5 returns the 10 K floor as its static T, silently; the
  mass flux and Mach built on it are then the floor's. Kept (it is the
  documented clamp, and 148 849 calls in the suite end there, all from
  probes). Options: a flag on the results that carry it, or let T go to the
  root (the gas is ideal to 0 K since #626).
- **R2 (E1, E3) `newton_damped` does not converge on trace species.** Two
  classes: (a) when CO or H2O is a trace the feasible WGS interval is
  [-0.11, 3e-9..2e-7] but the FD step is 1e-6, so the stencil reads past the
  bound (where the 1e-300 floor makes f jump) and Newton stops at xi ~ -8e-7
  with ln(Q/K) = +13..14; the composition error is bounded by the trace
  amount (~1e-6 mole). (b) the 1500 K rich equilibrium combustor
  (`test_fuel_boundary_892`) adiabatic T stops at 1390 K with 1729 J/kg left
  (~1.5 K). A fix (a bracketed Newton on the monotone ln Q, the stencil kept
  inside the bounds or an analytic slope) moves results and needs the Cantera
  comparison redone, so it is not in this PR. Throwing instead would fail the
  829 calls that end there today. 236 more cap exits are on an empty
  interval (xi_min = xi_max = 0), where 0 is the answer: an early return
  would only save time. `ReformingAdiabaticFun` also drops
  `solve_reforming_wgs_isothermal`'s converged flag.
  **Fixed in #950 (#927).** Both classes were larger than stated here: (a)
  where the equilibrium lies away from the trace the whole shift is missed
  (CO 7e-7 instead of 0.015 from CO 1e-8 at 800 K); (b) the 1729 is J/g
  (`h/mwmix`), and the root is 924.18 K, 466 K below the 1390 K returned
  (Cantera HP on the same species: 923.85 K). The WGS extent is now rtsafe on
  ln Q - ln Kp with the analytic slope (empty interval returns 0), the
  adiabatic T a bracket from the guess plus Brent on [T_FLOOR, 4500] K, the
  inner converged flag is propagated, and every non-root exit throws. Cap
  exits: suite 787 -> 0, Cantera validation 42 -> 0; every other suite and
  Cantera-validation equilibrium result unchanged to 5e-12 K and 3e-15 in
  mole fraction, Cantera deviations per test identical.
- **R3 (G2) `T_from_stagnation`'s derivative.** dF/dT drops the
  T dgamma/dT term, so Newton converges linearly and stops ~1e-8 T short;
  the exact slope (as `fanno_mach.cpp`'s `sound_speed_sq_dT`) makes it
  quadratic and moves results at ~1e-9 relative.
  **Fixed (#936, PR #954):** see "R3, R5, R6, R8, R11 fixed" below.
- **R4 (G3) Input validation for `chamber_static_state`.** It accepts a
  negative Tt (garbage from an upstream state) and, when the Newton happens
  to converge, returns a negative T. `duct_static_state` floors its inputs;
  this one could refuse Tt <= 0. Behaviour change for probes.
  **Fixed (#937, #948):** Tt <= 0 or non-finite raises `std::runtime_error` (the
  cap exit's type). Main returned T = -1.00001 K at Tt = -1 K and 0 at 0. The
  cascaded branching manifold (`test_tee_network`, xfail) now records its
  garbage-Tt probes (-4.3e18 K) as rejected node evaluations instead of
  Newton cap exits; it ends at the same |F| = 1.910e4 as main.
- **R5 (S1) `set_VH_mass` is `calc_T_from_h_mass`.** For an ideal gas
  u + R T = h, so the hand Newton (20 steps, 1e-4 K) is the safeguarded
  enthalpy inverse at a coarser tolerance; switching moves results by up to
  1e-4 K.
  **Fixed (#938, PR #954):** measured 2.2e-8 K at most (below).
- **R6 (T2) `calc_T_from_cp` branch.** cp(T) is not monotone for every
  mixture, so the root depends on `T_guess`, and the [200, 6000] K clamp
  predates the #606 tails. A bracketed solve over a stated range would make
  the branch explicit.
  **Fixed (#939, PR #954):** the root nearest `T_guess` in [200, 6000] K
  (below). **Range extended to [1, 6000] K** (decision 2026-11-03; "R6 range
  extended below 200 K" below).
- **R7 (C10, C11) `nozzle_quasi1d` is a best-effort profile.** Where
  A < A*_eff the area-Mach Newton has no subsonic root and holds the 0.999
  clamp (57 stations in one test); at a choked throat the root is M = 1
  exactly. A per-station converged flag, or solving M = 1 at A/A* = 1,
  would report it; the (T, P) loop uses an approximate da/dT.
- **R8 (F3) Exact slope in `fanno_inlet_mach`.** `isentropic_flux_slopes`
  already gives d flux/dM exactly; using it instead of the FD stencil moves
  converged results at rounding level (not bit-identical, so not done here).
  **Fixed (#940, PR #954)** (below).
- **R9 (G5, E2) Cap rule for the remaining bracketed solves.** Not reached
  and bracketed with bisection; the same rule could be applied for
  uniformity. E2's non-finite break returns a clamped extent silently.
  **G5 fixed (#941, #948):** `duct_static_state` applies the cap rule on its last
  ln P step (`inner_solve::kDuctStaticStateMaxIter` = 100, unchanged). Not
  reached: network F and J bit-identical (below). **E2 fixed (#950):** an
  amount rounded to zero at a bracket end is F's +-inf limit (bisect); a NaN
  and the cap throw.
- **R10 (C9) Delete `solve_station_state`.** Dead code (`[[maybe_unused]]`).
  **Fixed (#942, #948):** deleted; nothing else used only it.
- **R11 (O3) Jobson tolerance.** 1e-15 relative cycles in the last bits
  (2 cap exits at 3e-14); 1e-13 would end those calls ~190 iterations sooner
  with no change beyond 3e-14.
  **Fixed (#943, PR #954)** (below).
- **R12 (H1) `dewpoint` with RH > 1.** Supersaturated input returns T (the
  clamp) as the dewpoint; refuse RH > 1, or document dewpoint = T.
  **Done (#944, #948), and the finding was narrower than stated:** the public
  `dewpoint` never reaches the Newton with RH > 1, because `vapor_pressure`
  refuses RH outside [0, 1] (on main too: `dewpoint(300, 101325, 1.0000001)`
  raises; `std::invalid_argument` since #928). Only RH = 1 lands on the clamp,
  where T is the right answer. Documented as such (header, binding,
  API_CPP.md, API_PYTHON.md); no behaviour change.
- **R13 (J1, J2) Ejector bisections outside their range.** An area ratio
  below 1 (or above A/A*(50)) and an m_dot above the C-D nozzle's choked flow
  (or below its near-zero end) return a bracket end silently, J2 with its
  implicit Jacobian there. Clamp or refuse is a modelling decision for the
  ejector element.
  **Fixed (#945, #948), flag and clamp (user decision):** both keep the bisection's
  value. J2 (`ejector_cd_nozzle_exit_static_and_jacobian`) sets
  `out_of_range` (+1 above the flow at P*, -1 below the flow at
  Pt (1 - 1e-12)) and returns the clamped function's partials
  (dP/dm_dot, dP/dPt, dP/dTt) = (0, k, 0), k = r_crit or 1 - 1e-12,
  FD-verified. J1 is reported by the new internal `ejector_nozzle_exit_mach`
  and `EjectorElementResidualJacobian.nozzle_exit_mach_out_of_range`;
  `EjectorElement.diagnostics` reports `<id>.nozzle_exit_mach_out_of_range`.
  M_p1 is a constant in every Newton unknown, so its clamp's zero derivative
  is what the chain already used; the element rows FD-check at a clamped
  point (A_p1/A_t = 1e4, gamma 5/3) to 6e-7. Which paths reach them: J2 has
  no caller inside the library (the element solves P_py as an unknown; the
  Python `cd_nozzle_exit_static` is a value-only twin), so only a direct
  call can be out of range. J1's area ratio is a fixed parameter (> 1,
  checked at the element and, since #923, at the public entry); with gamma
  in (1, 5/3] it is clamped only for A_p1/A_t above ~7.8e3 (gamma 5/3) to
  ~1.5e6 (gamma 1.4), which no stress draw reaches (#948: 1786 ejector element evaluations over the 40 stress seeds, 0 flagged).

## R1 fixed: the static T goes to its root (#935, PR #949)

Decision (2026-11-02): "R1 -> can go to zero". The gas is ideal to 0 K
(#626, the C2 tails keep cp > 0 at every T), so a static-T solve returns its
root however small. T stays positive: an isentrope at P > 0 has a positive
root, and so do the energy solves when h0 > h(0 K).

**What changed.**

- C1 `solve_T_isentropic`: no floor. The step limit is still half of T. A
  step that would land below `thermo::T_LOG_STEP` (10 K, the former floor) is
  taken in ln T, as T exp(-F/cp). s is nearly cp ln T there, so P/Pt = 1e-300
  converges within the 50-step budget. A plain floor removal would need one
  halving per iteration and threw beyond P/Pt ~ 1e-20. Above 10 K the iterates
  are the T Newton's, so every solve that never reached the floor is
  unchanged. The step test is relative below 1 K (min(1 + T, 2 T)), which is
  bit-identical above 1 K.
- G1 `Tt_from_static_v`, G2 `T_from_stagnation`, C6 `solve_T_from_energy`:
  no floor. Their step limits (0.5 T, 0.3 T) keep the iterate positive, and
  their starts are close to the root. The tests stay absolute below 1 K. h
  carries the heat of formation, so its rounding resolves T only to ~1e-12 K,
  and a relative test cycled there. C6 has no root for h0 < h(0 K). It now
  throws on its cap, where it returned 10 K. That case cannot be reached from
  `fanno_channel`, whose h0 is an inlet state's.
- `s_R` held its value below 10 K, so s was flat there while ds/dT = cp/T was
  not. That hold is what made the floor the only answer. s_R is now exact for
  T > 0 and held only below `thermo::S_R_T_MIN` = 1e-300 K, so a non-positive
  probe T stays finite. That is far below any root a double can carry: T/Tt =
  (P/Pt)^(R/cp) >= ~1e-123. The low tail's `tail_S_R` took its end-value
  branch (-inf at T_end = 0 K) once t rounded to 1, which returned NaN below
  ~1e-14 K. It now keeps the fade branch there.
- `Pt_from_static_and_jacobian_M` zeroed dP0/dT0 for T0 < 1 K. The guard is
  now T0 > 0, which is unchanged above 1 K.

**Every `T_FLOOR` (and literal 10 K) use, decided.**

| Use | What it is | Decision |
|---|---|---|
| `solve_T_isentropic` (C1), `Tt_from_static_v` (G1), `T_from_stagnation` (G2), `solve_T_from_energy` (C6) | root clamp | fixed (above) |
| `s_R` hold (thermo.cpp) | property hold that made the root unreachable | fixed: hold moved to `S_R_T_MIN` |
| `h_RT` hold, `g_over_RT` (thermo.cpp) | property floor of g/RT for the equilibrium constants (1/T); every equilibrium caller floors its T at 10 K | kept; `g_over_RT` now holds both terms at `T_FLOOR` explicitly, so it is bit-identical to before |
| `complete_adiabatic`, `adiabatic_T_equilibrium_and_jacobians`, the combustor's equilibrium input (solver_interface.cpp) | input floor on a probe's inlet T, derivative zeroed when clamped | kept: combustion at < 10 K is a probe, not a state |
| `duct_static_state`, `duct_htc_static_state(_dY)` (stagnation.cpp) | input floor on a probe's Tt (with `P_FLOOR`), partials zeroed | kept: an input guard, not a root |
| `nozzle_flow` Tt, `fanno_channel(_rough)` T_in (literal `10.0`, now `thermo::T_FLOOR`) | input floor on a probe's stagnation or inlet T | kept (named, bit-identical) |
| `nozzle_quasi1d` station T clamp [T_FLOOR, Tt] (C11) | clamp in a best-effort profile | kept: R7, left as it is by decision |
| `solve_station_state` (C9) | dead code | kept: R10 deletes it |
| `T_adiabatic_wall` at v = 0 (literal 10 K), the leaky 10 K property floors in solver_interface.cpp (`T_eff = 10 + 0.1 (T - 10)`) | property and probe guards, no solve; the v = 0 value is pinned by `test_nan_guards_769` | kept, not R1 |
| side_stream_momentum.cpp, impingement_correlation.cpp, rib_correlation.cpp | no temperature floor: `MERGE_DISCRIMINANT_FLOOR`, `JET_ROW_MDOT_FLOOR`, `RE_FLOOR`, `BRACKET_FLOOR`, `EPLUS_FLOOR` | not R1 |

**Evidence.**

- *Kernel sweep, bit identity* (56 700 calls). The kernels are
  `mass_flux_isentropic`, `mach_from_pressure_ratio`, `nozzle_flow`,
  `critical_pressure_ratio`, `orifice_compressible_mdot_and_jacobian`,
  `T_from_stagnation`, `P_from_stagnation`, `Tt_from_static`,
  `Tt_from_static_v`, `Pt_from_static`, `stagnation_from_static`, both
  `*_and_jacobian_M`, `s`, and `fanno_channel`. The ranges are P/Pt 1e-9 to 1
  (115 ratios), Tt 100-2500 K, M 0.01-1000, static T 1-2500 K, and air, burnt
  gas, CH4, H2 and Ar. A throwaway counter marked every call whose old path
  touched the floor: an iterate below 10 K, or s at a T below 10 K. Of the
  calls that never touched it, 0 changed (44 683 identical). Calls that touched
  it: 11 908 changed and 109 came out identical.
- *The new values are roots.* Over 1800 isentropic points (P/Pt 1e-9 to 1,
  four Tt, five gases), the mass flux and Mach agree with an independent
  inverse (`calc_T_from_s_mass`, rtsafe) to 7e-14. The lowest static T is
  3e-4 K. An example: 300 K air at P/Pt = 1e-8 gives ~1.5 K, where it gave
  10 K. The stagnation residuals are below 1.4e-8 cp T. That is G2's linear
  convergence (R3), the same as above 10 K. Tt_from_static is below 1e-12.
- *Derivatives consistent* (FD). The isentropic identity dG/dP = (M^2 - 1)/v
  holds against a central difference of `mass_flux_isentropic` to 2.1e-8 for
  P/Pt from 0.99 down to 1e-9. `Tt_from_static_and_jacobian_M` and
  `Pt_from_static_and_jacobian_M` agree to 6.7e-5 for static T down to 1 K
  and M 0.01-1000. The worst case is at M = 0.01, where dT0 is ~1e-9 K per
  step and the inner tolerance limits the difference.
- *Network F and J* (SHA-256 at x0 + 3 seeded perturbations, all 14 families
  x seeds 0-39, 2240 evaluations): all identical to main.
- *Stress outcome A/B* (all 14 families x seeds 0-39, `--jobs 4`, serial,
  main then branch). 0 seeds lost convergence. Four seeds gained it:
  bleed_station 2, impingement 1, random 3 and random 34 (residual_too_large,
  timeout, timeout, no_progress on main). Rerun alone (`--jobs 1`), all four
  converge on main as well, with the same final norm to the bit. All four are
  wall-clock flips, not the fix. 15 seeds had the same outcome but a different
  final norm. Rerun alone, main's own final norms move as much between runs
  (the stall detector and the timeouts are wall-clock). Converged counts per
  family (branch run): bleed_station 37, effusion 40, ejector 30, hotcold 40,
  hotcold_mild 40, ideal_nport 40, impingement 40, impingement_bypass 40,
  jet_chain 38, legacy_tee 32, pin_fin 40, random 32, series 40, tee3 32.
- *Tests.* `tests/test_inner_solve_cap_868.cpp` drops the old pin of the
  floor return and adds two cases. The first checks that the static T goes to
  its root (T_from_stagnation at M = 50, Tt_from_static from 5 K, the
  isentrope to P/Pt = 1e-8 against `calc_T_from_s_mass`, and P/Pt = 1e-300).
  The second checks that s is exact and finite to 1e-20 K and at a
  non-positive T. Both fail on main. The second also fails with the tail
  branch fix reverted (NaN). The full suite (6808 passed), ctest (53) and the
  paper examples pass.

## R3, R5, R6, R8, R11 fixed (#936, #938, #939, #940, #943, PR #954)

Decision (2026-11-02): the R-items are fixed. These five move results
slightly. None is bit-identical.

**What changed.**

- R3 (G2) `T_from_stagnation`: dF/dT = cp + (M^2/2) R (gamma + T dgamma/dT),
  with dgamma/dT = -R_u cp'/cv^2, the slope `sound_speed_sq_dT` uses.
- R5 (S1) `State::set_VH_mass`: `calc_T_from_h_mass(h, X, T)`, as
  `set_HP_mass`. `inner_solve::kSetVHMaxIter` and `kSetVHTolK` are deleted.
- R6 (T2) `calc_T_from_cp`: the range is [`thermo::CP_INVERSE_T_MIN`,
  `CP_INVERSE_T_MAX`] = [200, 6000] K, where every species has NASA-9 data.
  From 6000 K (H2O, NH3, CH4, the alkanes) the #606/#607 tails level cp off
  toward a constant, so cp no longer identifies T. Below 200 K the low tail
  fades cp to the rigid-rotor value; extending the range below 200 K was not
  done here (it was later: "R6 range extended below 200 K"). The branch rule: scan outward from `T_guess`, clamped into the range,
  10 K at a time on both sides, then run a bracketed Newton on the first sign
  change, so the root returned is the one nearest the guess to within a step.
  N2 is the case that needs it: cp has a minimum near 214 K, so cp(205 K) is
  reached again at 226.74 K. No sign change in the range throws
  `std::runtime_error` ("no root in [200, 6000] K", with cp at both ends).
  Before, #924's cap exit threw only after 50 steps. A non-finite or
  non-positive `cp_target`, `T_guess` or `tol` is `std::invalid_argument`.
  On main a negative guess returned 530 K for air.
- R8 (F3) `fanno_inlet_mach`: the Newton slope is
  `isentropic_flux_slopes(...).dM`, not the finite-difference stencil.
- R11 (O3) `jobson_expansibility`: the tolerance is
  `inner_solve::kJobsonExpansibilityRelTol` = 1e-13.

**Evidence.**

- *Kernel sweep, main vs branch* (9162 calls). The gases are air, a burnt
  gas, CH4, H2, Ar and steam/N2.
  - R3, `T_from_stagnation`: Tt from 100 to 2500 K and M from 0.01 to 10,
    2250 calls. Static T moves by up to 1.1e-9 relative (4.8e-7 K), with a
    median of 5e-15; 388 values are identical. `P_from_stagnation` moves by
    up to 5.7e-9 relative. The energy residual |F|/(cp T) was up to 2.4e-8
    (median 3.7e-14) and is now up to 9.0e-14 (median 2e-16): quadratic
    convergence. FD of dT/dTt agrees with the analytic slope to 3.5e-10,
    against 1.7e-9 on main.
  - R5, `set_VH_mass`: T from 150 to 3000 K and two start temperatures, 360
    calls. T moves by up to 2.2e-8 K (median 5.8e-9). |h(T) - h|/cp was up
    to 1.2e-8 K and is now up to 7e-13 K. The 1e-4 K the issue quotes is the
    old step tolerance; the actual error was smaller.
  - R6, `calc_T_from_cp`: targets cp(T) for T from 200 to 6000 K and guesses
    of 300, 1000 and 3000 K. Excluding Ar, 885 calls: every one is the same
    root as on main, moved by up to 8.0e-7 K (median 2e-12). Ar's cp is
    constant to within data noise (dcp/dT changes sign at 950-1642 K), so
    every T is a root to rounding. On main 44 of its 177 calls threw and the
    rest returned arbitrary T. On the branch all return a point where cp
    matches, near the guess.
  - R8, `fanno_inlet_mach`: G/G_max from 1e-14 to 1 and four Tt, 1440
    calls. M_in moves by up to 3.4e-14 relative; 772 values are identical.
  - R11, `jobson_expansibility`: r from 0.05 to 0.999, beta up to 0.8, Ci
    from 0.6 to 1 and n from 1.2 to 1.667, 1800 calls. Results move by up to
    9.8e-14; 960 are identical. 108 calls throw on both main and the branch:
    at beta 0.7-0.8 with Ci 0.8-1, and at beta 0.5 with Ci 1, for r below
    ~0.56, the fixed point does not converge (step -2.3 at the cap). That
    failure predates this PR and is not changed here; see "Not verified".
- *Network F and J* (SHA-256 at x0 + 3 seeded perturbations, 14 families x
  seeds 0-39). 344 of 560 differ. These families are unchanged: ejector,
  legacy_tee and series. The largest change relative to the largest entry,
  per family: impingement 5.8e-7, random 2.7e-7, jet_chain 1.3e-8, bleed_station
  7.2e-9, effusion 6.7e-9, impingement_bypass 2.6e-9, tee3 1.9e-13, and
  1e-16 or below elsewhere. A build with R8 reverted made random 4,
  impingement 15 and effusion 0 identical to main, so R8 accounts for those
  two largest entries. Both are d(row)/d(m_dot) of a low-flow ChannelElement
  (m_dot ~ 2e-3 kg/s). There the entry moves by 8e-5 (random 4) and 1.4e-5
  (impingement 15) relative. A central difference disagrees with the
  analytic value by 1.5e-3 (random 4) and ~1e-6 (impingement 15) on main and
  on the branch alike, so the change lies within the existing FD scatter at
  that point.
- *Stress outcome A/B* (all 14 families x seeds 0-39, `--jobs 4`, serial,
  main then branch): 0 seeds changed outcome. Converged per family, the same
  on both: bleed_station 38, effusion 40, ejector 30, hotcold 40,
  hotcold_mild 40, ideal_nport 40, impingement 40, impingement_bypass 40,
  jet_chain 39, legacy_tee 32, pin_fin 40, random 33, series 40, tee3 32.
  343 final norms differ (the norms are wall-clock noisy).
- *Tests.* Python (`test_inner_solve_cap_868.py`) and C++
  (`test_inner_solve_cap_868.cpp`) cover the following:
  - R3: the energy residual is below 1e-13 cp T at 2500 K and M 1.8-5.6. On
    main it was 2e-10 to 8e-10.
  - R6: the nearest root (N2 at 213 K and 218 K guesses, guesses outside the
    range), no root in range (1000 and 20 J/(mol K)), and ValueError /
    `std::invalid_argument` for invalid input.
  - R5: `set_VH_mass` equals `calc_T_from_h_mass`.

  All fail on main. Falsified on the branch: the old slope fails the R3
  tests, and an upward-only scan fails the nearest-root tests, in Python and
  C++. R8 and R11 have no new test, because their results move at rounding
  level. The existing vanishing-flux test still passes with the slope scaled
  by 1.5, since the bracket keeps Newton on the root. The full suite (6901
  passed), ctest (54) and the paper examples pass.

## R6 range extended below 200 K (2026-11-03, PR #961)

Decision (2026-11-03): `calc_T_from_cp` searches below 200 K. No issue; it
follows #939 (closed by PR #954).

**The tail below 200 K.** Below the NASA-9 data (200 K for every species)
the #606 low tail is a quintic in t = (200 K - T)/200 K: C2 with the data at
200 K, cp/R = the rigid-rotor value (5/2 atom and H2, 7/2 linear, 4
nonlinear) at 0 K, with zero slope and curvature there. So the gas is ideal
down to 0 K (no lower validity limit; #626 retired the composition floors),
cp is continuous, and cp - cp(0 K) ~ (T/200 K)^3 near 0 K. cp/R on
[1e-3, 200] K, 20001 log points:

| species | cp/R at 0 K | cp/R - cp/R(0) at 1 K | at 10 K | extrema below 200 K |
|---|---|---|---|---|
| N2 | 3.5 | +1.7e-8 | +1.6e-5 | max 3.50568 at 127 K |
| O2 | 3.5 | -4.8e-8 | -4.4e-5 | min 3.48491 at 124 K |
| CO | 3.5 | -5.7e-10 | -5.2e-7 | min 3.49988 at 101 K |
| CO2 | 3.5 | -4.0e-9 | +9.5e-7 | min 3.5 - 2.7e-7 at 6.2 K |
| H2O | 4 | -5.1e-8 | -4.7e-5 | min 3.98602 at 117 K |
| CH4 | 4 | -1.5e-7 | -1.4e-4 | min 3.95224 at 123 K |
| NH3 | 4 | -8.7e-8 | -7.9e-5 | min 3.97829 at 111 K |
| H2 | 2.5 | +6.5e-7 | +6.1e-4 | none (increasing) |
| C2H6 .. NC7H16 | 4 | +4.8e-7 .. +1.1e-5 | +4.6e-4 .. +1.0e-2 | none (increasing) |
| AR | 2.5 | 0 | 0 | constant |

**The lower bound: `thermo::CP_INVERSE_T_MIN` = 1 K.** The tail has no
validity limit of its own, so the bound is set by conditioning: a root's
relative error from cp rounded to the last bit is ~eps cp / (3 |cp - cp(0)|).
For CO, the flattest species, that is ~5e-7 at 1 K and ~5e-4 at 0.1 K; below
1 K every species' cp is within 3e-6 relative of its 0 K value. A mixture can
cancel the cubic further (air: +3.5e-9 at 1 K, N2 and O2 of opposite sign).

**The scan.** Nodes are 10 K apart (`CP_INVERSE_SCAN_STEP`) down to
`thermo::CP_INVERSE_LOG_SCAN_T` = 100 K, and a factor 0.9 apart below it (the
same step relative to T, continuous at the switch; 44 nodes from 100 K to
1 K), so that each step spans a like share of the T^3 approach and roots a
few K apart near 0 K are separated (CO2's cp(4 K) is also reached at 7.6 K,
within one 10 K step). The downward scan steps onto `thermo::CP_INVERSE_T_DATA`
= 200 K. Each step advances the side whose next node is nearer `T_guess` (a
tie, within 1e-12 relative, goes up), so the root returned is nearest in K
to within the local step. With equal 10 K steps on both sides this is the
former up-then-down alternation; one step per side per round would have made
the 10 % steps below 100 K lose to the 10 K steps above (O2's cp(1.5 K) is
also reached at 189.9 K: from a 20 K guess the alternation returned 189.9 K).

**Evidence.**

- *Sweep, branch vs main* (45825 calls): 24 compositions (14 pure species
  without Ar, air, 10 random mixtures), targets cp(T) for 60 T log-spaced in
  [1, 199] K and 81 T in [200, 5990] K, 13 guesses from 1 to 9000 K. An
  independent reference finds every root in [1, 6000] K (sign changes on
  40000 log-spaced points, brentq to 1e-14 K).
  - Targets with roots only below 200 K (19240 calls): main throws on all.
    The branch returns a root on 19203. Every result is a reference root to
    1.0e-7 relative (CO near 1 K, the conditioning above) and has
    |cp(T) - target| / target <= 3.0e-15. It is the nearest root on 19190;
    on 13 it is not, by at most 0.97 of the local step. 37 throw "no root":
    the target sits within a step of a tail extremum, so both roots lie
    inside one step (O2 at 127.07 K and 120.76 K, H2O near 116 K, CO2 near
    6.0 K): the documented miss.
  - Targets with roots only at or above 200 K (26143 calls): for a guess in
    [200, 9000] K all 14077 are bit-identical to main. For a guess below
    200 K, 10303 are bit-identical (the scan reaches 200 K on main's nodes)
    and 1763, all at a 199 K guess, differ: main clamped the guess to 200 K,
    so the nodes are offset by 1 K; same roots.
  - Targets with roots on both sides (442 calls): 248 differ from main (a tail
    root nearer the guess, or a guess below 200 K); 427 are the nearest root and
    15 are not (N2 near its 127 K tail maximum, two roots 0.76 K apart, and
    a near-tie at 210/220 K, within the step).
- *Network F and J* (SHA-256 at x0 + 3 seeded perturbations, 14 families x
  seeds 0-39): 560 of 560 identical to main. No network code calls
  `calc_T_from_cp`.
- *Stress outcome A/B* (all 14 families x seeds 0-39, `--jobs 4`, serial,
  main then branch): 6 seeds changed outcome (bleed_station 2, 18; jet_chain
  1, 19; random 3 lost, random 34 gained). Rerun one at a time (`--jobs 1`,
  load average ~30 from other agents) they are wall-clock flips: each ends at
  ~20 s (the stall and timeout budgets), and bleed_station 2 converged on the
  branch on the second rerun, with the same final norm as main
  (8.13e-12). Since F and J are identical, the solve is the same.
- *Tests* (`test_inner_solve_cap_868.py` / `.cpp`): the tail is reached (H2
  at 50 K, C3H8 at 1.5 K, NC7H16 at 10 K, guess 300 K); CO2's two roots near
  its 6.2 K minimum are separated; the nearer-node rule (O2, 1.5 K from a
  20 K guess); N2's tail root (37.4 K) from a 50 K guess; no root in
  [1, 6000] K. Falsified in Python and C++: `CP_INVERSE_T_MIN` = 199 K fails
  the tail tests, 10 K steps below 100 K fail the CO2 test, and one step per
  side per round fails the O2 test.

## Not verified

- The stress families were evaluated at x0 and perturbations (F, J hashes),
  not solved: no stress A/B was run, so whether a solve that used to accept a
  non-root now takes a different path in a full stress solve is not measured.
- MSVC and MinGW builds (CI covers them).
- The Fanno marches' step budget (C7, C8) and most "latent" cap exits could
  not be reached with constructed inputs; their throws are untested beyond
  compilation and the unchanged suite.
- R3/R5/R6/R8/R11 (PR #954): the `calc_T_from_cp` scan can miss two roots
  closer than a step (a cp extremum within it of the target): 10 K, or
  10 % of T below 100 K since the range reaches 1 K; then it throws "no root".
  Below 1 K and for a mixture whose tail cubic cancels, the root is
  ill-conditioned (not measured beyond the sweep above).
  `jobson_expansibility` fails to converge on valid input at high beta Ci and
  low r: 108 of the 1800 sweep calls, the same on main. It is not fixed here.
