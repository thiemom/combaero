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
| T2 | thermo.cpp:751 | `calc_T_from_cp` | Newton, step <= 0.1 T, clamped [200, 6000] K | 50 (caller) | 1e-6 abs on step or cp | silent: returns the clamp for an unreachable target (6000 K for cp = 1000 J/(mol K)) | none; cp(T) need not be monotone | bug (reachable by the caller) | fixed: throws. rec R6 |
| S1 | state.cpp:118 | `State::set_VH_mass` | Newton on u + R T = h | 20 | 1e-4 K abs step | silent | none (monotone) | latent | fixed: throws. rec R5 |
| G1 | stagnation.cpp:77 | `Tt_from_static_v` (and `Tt_from_static`, `Pt_from_static`, `stagnation_from_static`) | Newton, step <= 0.5 T, floor `T_FLOOR` | 50 (caller) | 1e-8 (1 + T) step | silent (also on a flat derivative) | none (h monotone) | latent (reachable by the caller) | fixed: throws; a root below the floor still returns the floor |
| G2 | stagnation.cpp:132 | `T_from_stagnation` (and `P_from_stagnation`) | Newton with approximate dF/dT (gamma' dropped), floor | 50 (caller) | 1e-8 (1 + T) step | silent | none | latent | fixed: throws, floor kept. rec R3 |
| G3 | stagnation.cpp:324 | `chamber_static_state` | Newton from Tt | 100 | 1e-13 Tt step | silent; implicit J then evaluated at the non-root | none (f falls in T, one root) | bug (48 probe hits at Tt = -1.4e19 K) | fixed: cap rule. rec R4 |
| G4 | stagnation.cpp:408 | `sonic_static_state` | Newton | 100 | 1e-13 Tt step | silent, implicit J at non-root | none | latent | fixed: cap rule |
| G5 | stagnation.cpp:489 | `duct_static_state` | Newton in ln P on [ln P*, ln Pt] | 100 | 1e-15 abs step in ln P | silent | bisection when Newton leaves the bracket | safe | none. rec R9 |
| C1 | compressible.cpp:41 | `solve_T_isentropic` (`nozzle_flow`, `mass_flux_isentropic`, `mach_from_pressure_ratio`, `critical_pressure_ratio`, the solves below) | Newton, step <= 0.5 T, floor | 50 (caller) | 1e-8 (1 + T) step | silent; 148 849 hits, all at the floor | none (s monotone) | latent (floor hits are the documented clamp) | fixed: a root below the floor returns the floor at once (same value; the 49 clamped repeats are skipped), any other cap exit throws. rec R1 |
| C2 | compressible.cpp:120 | `find_critical_pressure` | golden section on [0.3, 0.99] Pt | >= 200 | width < tol (1 + a + b) | silent, bracket 0.618^200 | fixed bracket (P*/Pt 0.49-0.61 for gamma 1.67-1) | safe | none |
| C3 | compressible.cpp:292 | `solve_P_back_from_mdot` | bisection on [P*, Pt] | >= 200 | G within tol, width tol (Pt - P) | silent; 1 hit at 1 ulp | bracket | safe | none |
| C4 | compressible.cpp:341 | `solve_Pt_from_mdot` bracket | x10 expansion | 12 | | throws if not bracketed | | safe | none |
| C5 | compressible.cpp:358 | `solve_Pt_from_mdot` | bisection | >= 200 | as C3 | silent; 1 hit at 1 ulp | bracket | safe | none |
| C6 | compressible.cpp:462 | `solve_T_from_energy` (Fanno marches) | Newton, step <= 0.3 T, floor | 50 | 1e-8 T step | silent | none (F rises in T: u grows with T at fixed P) | latent | fixed: throws, floor kept |
| C7 | compressible.cpp:647 | `fanno_channel` march | adaptive RK4 | `kFannoMaxSteps` = 2e5 | | silent: outlet at x < L reported unchoked | step halving near M = 1 | latent | fixed: throws |
| C8 | compressible.cpp:894 | `fanno_channel_rough` march | adaptive RK4 | 2e5 | | silent, as C7 | | latent | fixed: throws |
| C9 | compressible.cpp:1094 | `solve_station_state` | Newton | caller | tol T step | returns `false` | | safe (reports); `[[maybe_unused]]`, no caller | rec R10 |
| C10 | compressible.cpp:1320 | `nozzle_quasi1d` area-Mach | Newton, M clamped to [0.001, 0.999] subsonic, <= 10 | 50 | abs 1e-10 in A/A* | silent; 57 hits where A < A*_eff (no subsonic root; at a choked throat the root is M = 1, held at 0.999) | clamps | best effort by design | rec R7 |
| C11 | compressible.cpp:1372 | `nozzle_quasi1d` (T, P) per station | sequential Newton, approximate da/dT, clamps | 50 (caller) | tol relative on s and h | silent | clamps T <= Tt, P in [100, Pt] | latent | rec R7 |
| F1 | fanno_mach.cpp:64 | `static_T` | Newton in a bracket | 100 | 1e-14 T step | silent | bisection when Newton leaves | safe | fixed: cap rule |
| F2 | fanno_mach.cpp:135 | `brent` | Brent | 200 | caller rel_tol (1e-13, 1e-14, `kFannoMachRelTol`) | silent at cap (unreachable with a bracket); NO check for a sign change | Brent's bisection fallback | bug: called unbracketed 55 times (F7) | fixed: throws without a sign change |
| F3 | fanno_mach.cpp:323 | `fanno_inlet_mach` | Newton, central FD slope h = 1e-7 max(M, 1e-3), in [0, 1] | 60 | 1e-15 M step | silent; 995 hits at G < 1.4e-11 G_max, M up to 73% high | bisection when Newton leaves | bug: below M = h the stencil crossed M = 0, where the flux is held at 0, read half the slope, and Newton 2-cycled | fixed: one-sided stencil when M - h <= 0 (only those iterations change), cap rule. rec R8 |
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
| O3 | orifice.cpp:1778 | `jobson_expansibility` | fixed point (approach term O(beta^4)) | 200 | 1e-15 relative | silent; 2 hits at 3e-14 / 4e-15 (rounding cycle) | | safe at the hits | fixed: cap rule (the hits are kept, bit-identical). rec R11 |
| H1 | humidair.cpp:205 | `dewpoint` | Newton, FD slope, clamped [T_lo, T] | 100 | 1e-6 K | throws; but RH > 1 puts the root above T, Newton sits on the clamp and "converges" to T | clamp | latent | rec R12 |
| H2 | humidair.cpp:317 | `wet_bulb_temperature` | Newton in a bracket | 200 | 1e-6 K | throws | bisection | safe | none |
| R1 | friction.cpp:92 | `friction_colebrook` | Newton from Serghides | 20 (caller) | 1e-10 relative | silent (also on a flat slope) | none (monotone) | latent (reachable by the caller) | fixed: throws |
| E1 | equilibrium.cpp:80 | `newton_damped` (WGS extent at 165, 276, 441, 845; adiabatic T at 217, 972) | damped Newton, FD slope eps = 1e-6 max(1, abs x), clamped | 40 / 50 | abs 1e-12 on f, or xtol_rel | silent on cap and on a flat slope | clamps only | bug (787 + 42 hits) | rec R2 |
| E2 | equilibrium.cpp:386 | `solve_reforming_isothermal` | Newton in a bracket | 200 | abs 1e-12 on F, width 1e-15 span | silent; a non-finite F breaks to a silent clamp | bisection | latent | rec R9 |
| E3 | equilibrium.cpp:605, 672 | `joint_newton` + line search | scaled Newton, curved step, backtracking | `REFORMING_NEWTON_*_IT` / 60 | `REFORMING_NEWTON_FTOL` | returns false -> `EquilibriumResult.converged` | fraction to the boundary | safe (reports) | rec R2 (one caller drops the flag) |
| E4 | equilibrium.cpp:~750 | `element_potential_start` | Gordon-McBride Newton | `REFORMING_NEWTON_MAX_IT` | 1e-9 | returns false | RP-1311 damping | safe | none |
| B1 | combustion.cpp:1104 | `bisection_solve` (fuel/oxidizer stream inverses) | Brent | 100 (caller) | tol_f, 1e-10 rel width | silent (returns b); unreachable with its checked bracket | throws if not bracketed | safe | none |
| M1 | acoustics.cpp:1041 | `refine_root_muller` | Muller | `kMullerMaxIter` = 50 | 1e-6 rad/s | returns last iterate by design | | safe: both callers (can_annular_solvers.cpp:188, 316) accept a mode only if abs(D) < `kDispersionTolerance` y_scale | none |
| M2 | can_annular_solvers.cpp:432 | `refine_root_muller_annular` | Muller | 50 | 1e-6 | as M1 | | safe: caller (549) checks abs(D) | none |
| X1 | crossflow_cd.cpp:144 | `florschuetz_metzger_su_1982_cd` falling branch | Newton, clamped (0, 4 cd0] | 100 | 1e-15 cd0 step | silent; implicit dCd/dq then at a non-root | none (one root, dF/dCd > 0) | latent | fixed: cap rule |
| J1 | ejector.cpp:25 | `ejector_mach_from_area_ratio_supersonic` | bisection on [1 + 1e-12, 50] | 200 | | silent end point for an area ratio outside [1, A/A*(50)] | | latent | rec R13 |
| J2 | ejector.cpp:168 | `ejector_cd_nozzle_exit_static_and_jacobian` | bisection on [P*, Pt (1 - 1e-12)] | 200 | | silent end point (and implicit J there) for m_dot outside the nozzle's range | | latent | rec R13 |
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
decision: C10/C11, E1/E2, J1/J2, H1 and the `T_FLOOR` return (R1).

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

- **R1 (C1, G1, G2, C6) The `T_FLOOR` return is not a root.** An expansion
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
- **R3 (G2) `T_from_stagnation`'s derivative.** dF/dT drops the
  T dgamma/dT term, so Newton converges linearly and stops ~1e-8 T short;
  the exact slope (as `fanno_mach.cpp`'s `sound_speed_sq_dT`) makes it
  quadratic and moves results at ~1e-9 relative.
- **R4 (G3) Input validation for `chamber_static_state`.** It accepts a
  negative Tt (garbage from an upstream state) and, when the Newton happens
  to converge, returns a negative T. `duct_static_state` floors its inputs;
  this one could refuse Tt <= 0. Behaviour change for probes.
- **R5 (S1) `set_VH_mass` is `calc_T_from_h_mass`.** For an ideal gas
  u + R T = h, so the hand Newton (20 steps, 1e-4 K) is the safeguarded
  enthalpy inverse at a coarser tolerance; switching moves results by up to
  1e-4 K.
- **R6 (T2) `calc_T_from_cp` branch.** cp(T) is not monotone for every
  mixture, so the root depends on `T_guess`, and the [200, 6000] K clamp
  predates the #606 tails. A bracketed solve over a stated range would make
  the branch explicit.
- **R7 (C10, C11) `nozzle_quasi1d` is a best-effort profile.** Where
  A < A*_eff the area-Mach Newton has no subsonic root and holds the 0.999
  clamp (57 stations in one test); at a choked throat the root is M = 1
  exactly. A per-station converged flag, or solving M = 1 at A/A* = 1,
  would report it; the (T, P) loop uses an approximate da/dT.
- **R8 (F3) Exact slope in `fanno_inlet_mach`.** `isentropic_flux_slopes`
  already gives d flux/dM exactly; using it instead of the FD stencil moves
  converged results at rounding level (not bit-identical, so not done here).
- **R9 (G5, E2) Cap rule for the remaining bracketed solves.** Not reached
  and bracketed with bisection; the same rule could be applied for
  uniformity. E2's non-finite break returns a clamped extent silently.
- **R10 (C9) Delete `solve_station_state`.** Dead code (`[[maybe_unused]]`).
- **R11 (O3) Jobson tolerance.** 1e-15 relative cycles in the last bits
  (2 cap exits at 3e-14); 1e-13 would end those calls ~190 iterations sooner
  with no change beyond 3e-14.
- **R12 (H1) `dewpoint` with RH > 1.** Supersaturated input returns T (the
  clamp) as the dewpoint; refuse RH > 1, or document dewpoint = T.
- **R13 (J1, J2) Ejector bisections outside their range.** An area ratio
  below 1 (or above A/A*(50)) and an m_dot above the C-D nozzle's choked flow
  (or below its near-zero end) return a bracket end silently, J2 with its
  implicit Jacobian there. Clamp or refuse is a modelling decision for the
  ejector element.

## Not verified

- The stress families were evaluated at x0 and perturbations (F, J hashes),
  not solved: no stress A/B was run, so whether a solve that used to accept a
  non-root now takes a different path in a full stress solve is not measured.
- MSVC and MinGW builds (CI covers them).
- The Fanno marches' step budget (C7, C8) and most "latent" cap exits could
  not be reached with constructed inputs; their throws are untested beyond
  compilation and the unchanged suite.
