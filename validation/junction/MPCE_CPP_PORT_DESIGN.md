# MPCE Junction: C++ Port Design and Provenance Audit (issue #271)

Deep-dive read of the three source papers on disk (Bassett 2001, Hager 1984,
Mynard & Valen-Sendstad 2015), the internal compressible spec
(`compressible_junction_model_v3.pdf`), the Matlab reference
(`JunctionLossCoefficient.m`), and the code. Written to decide *how* to port
`MPCEv2Element` to C++, what has to be fixed before that is safe, and what the
digitised data can and cannot prove.

## 0. Policy checklist -- re-read before starting ANY step

Context drift is the enemy over a long piece of work. Each step opens by
re-reading this list and stating which items bind it. Items are the user's
decisions (2026-09-04), not mine.

1. **Ground truth is digitised measured data for our regime** (engine air,
   high Re, compressible). Bassett, Hager, Idelchik, Wang, Perez-Garcia. Not
   a paper's own CFD, not the current code's output, not a 3-point fixture.
2. **Every tuned constant proves itself with an on/off (or old/new)
   before-after table** on the scorecard cells covering its range, and is
   labelled at its definition: source, fitted range, date, proving cell.
   CFD-derived corrections are tunings; the source CFD is not reproduced.
3. **Retune after every fix, physics or code.** A constant fitted with a bug
   present may have absorbed it. The intermediate state (bug fixed, tuning
   stale) can look worse -- that is incomplete follow-through, not a wrong
   fix. Follow it through: fix, re-fit, measure.
4. **Never falsify a change by removing one half of a matched pair.** If two
   terms were designed together, test them together (#272 gate + cross term).
5. **Residual-form change is stepwise and tested.** Faithful `Pt` form first
   with an exact equivalence gate; Mynard-native `C_j` later, data-gated.
   Never both in one PR.
6. **Compressible corrections meet the same bar** as any other tuning.
7. **Baseline before, score after, revert on a no.** Write the numbers down
   either way. A negative result reported honestly is a good outcome.
8. **A declared limitation is an unmeasured assumption until it has a
   number.** ("K5 not reproduced", "dead code", "gate is load-bearing" --
   all three were wrong this arc.)
9. **Reachable-but-unexercised code gets a test on its precondition; a known
   defect left unfixed is a `strict=True` xfail asserting the correct
   behaviour**, never a green test on the wrong one.

## 1. Verdict in five lines

1. **The port is an assembly job, not a rewrite.** A complete C++ Unified0D
   kernel with an O(M^2) compressibility correction already exists
   (`include/tee_junction.h`, retired with the old tee), and the v3 spec is a
   20-page port design with the analytical Jacobian tables and C++ fragments
   already written. What is missing in C++ is Mynard's pseudosupplier
   construction and energy-transfer factor -- ~60 lines of Python.
2. **MPCEv2 is Mynard, not Bassett.** It computes Mynard's static-pressure
   coefficient `C_j` and converts to `K` via his Eq 18. That conversion is what
   restricts it to 3 ports; Mynard's native `C_j` form works for any N.
3. **MPCEv2 has no compressibility correction.** The retired tee had O(M^2).
   The production junction regressed on this and nothing documents it.
4. **Three constants in the closure have no provenance in this repo**, and one
   was calibrated on plumbing that has since been fixed.
5. **Most of the digitised data is not scored.** 57+5+36+30 files exist;
   the network adapters score 34 of them, all Bassett. Hager's own K_straight
   measurements and 36 Idelchik files sit unused.

## 2. What already exists

| artefact | what it holds | state |
|---|---|---|
| `include/tee_junction.h` | Bassett K1-K12 raw + blended forms; `K_dat_j_closed<T>` = Unified0D K with O(M^2) kappa; `BranchState<T>` templated on `double`/`complex`; complex-step `cstep_deriv` | retired with `TeeJunctionElement` (#209), still compiled |
| `compressible_junction_model_v3.pdf` sec 8-9 | n+1 residual set, full analytical Jacobian tables for pseudodatum and branch primitives, C++ pseudodatum + kernel + outer-chain fragments, choking guard | design only; never wired to MPCE |
| `src/ejector.cpp` `DualN<N>` | forward-mode dual with `+ - * /`, `dsqrt`, `dpow` | used by the ejector port (#270) |
| `_mynard2010.py` | faithful port of the Matlab + combaero's `joining_etransfer_alpha` extension | production, 100% Python |
| `_mpce_v2_jacobian.py` | sympy `dKQ/dmdot` for the 3-port separating T only | production; joining has none |

So the C++ that the port needs to *write* is: the pseudosupplier angle/area
construction (Mynard Eqs 31-35), the energy-transfer factor (Eq 36), the
Matlab's damping regulariser, combaero's alpha term, and the C -> K conversion
(Eq 18). Everything else is reuse.

**The v3 spec already hit the port's hardest problem and solved it.** Its
"Changes in v2" section records that templating the pseudodatum on `complex`
broke complex-step derivatives because `std::atan2` and `std::abs` are
non-analytic -- so it keeps the pseudodatum in `double` and differentiates only
the kernel, chaining through the outer partials analytically (sec 9.8-9.9).
Mynard's reorientation logic (the +pi flip and the theta -> -theta mirror in
`_mynard2010.py`) is the same class of problem, plus two more branch points.

## 3. What MPCEv2 actually computes

The residual is `Pt_i - Pt_jct + sign * K_i * q_dyn_com = 0`, with `K_i` from
`_mynard2010.py`:

    C_j = damping * (1 - (1/(AR_j * FR_j)) * cos(3/4 (pi - phi_j)))     Mynard Eq 30
    K_j = (U_j^2/U_com^2) * (2 C_j + U_S^2/U_j^2 - 1)                  Mynard Eq 18, rearranged

`AR_j` is the *collector-specific* pseudosupplier area ratio (Eq 35), which is
where Mynard's energy-transfer factor `eta_j` (Eq 36) enters. So:

- It is **not** Bassett's Table 2. Bassett's K5/K6/K11/K12 are the *validation
  targets*; the closure is Mynard's generalisation of Bassett's CV with a
  pseudosupplier and a CFD-fitted energy exchange.
- `K` is only defined by the Matlab for `n <= 3` (its line 65; the Python
  `len(U) <= 3`). Mynard's `C_j` (Eq 16, static-pressure form) is defined for
  any N -- the paper's whole point is "any number of branches". The N > 3
  defect in #271 is an artefact of converting to `K` when no single common
  branch exists.
- **No Mach dependence anywhere.** `q_dyn_com` is `0.5 rho u^2` with rho from
  the port static state; `K` is incompressible Mynard. The retired tee's
  `K_dat_j_closed` carried `K_inc (1 + kappa M_dat^2)`; the header states the
  incompressible assumption needs `Ma < 0.2`. Nothing in `MPCEv2Element`
  says this.

## 4. Provenance audit of the closure

Every constant and structural choice, traced to its source or flagged.

| item | where | provenance | status |
|---|---|---|---|
| Hager 3/4 inflow-deviation angle | `cos(0.75 (pi - phi))` | Hager 1984 sec 2.2, observed "epsilon = delta/4" at mid-lateral; Bassett Eq 22, 33-34; Mynard Eq 19-20 | **grounded** -- three independent papers |
| Contraction coefficient `1/xi = 1 + sqrt(1 + 1/(lambda psi)^2 - 2 cos(phi)/(lambda psi))` | inside `C_j` via `AR * FR` | Hager Eq 14 = Bassett Eq 22 = Mynard Eq 28 | **grounded** |
| Dividing-streamline pressure `p* = (p_1 + p_s)/2` | inside Bassett K5 (Eq 14), Hager Eq 3, Mynard Eq 26 | all three papers, 2D plane-duct origin (Hager: "appropriate for plane duct flow") | **grounded, with a stated 2D limitation** |
| `eta_j = [0.8 (pi - theta_dat) sign(theta_j) - 0.2] (1 - lambda_j)` | `_mynard2010.py` etransfer | Mynard Eq 36, "determined only from CFD data in Figure 4", Re 350-2400 blood flow, 3-branch, T and Y only | **empirical fit outside combaero's regime**; Mynard sec 4.5.1 argues high-Re validity via Gardel/Levin agreement, not by test |
| `damping = 1 - exp(-FR/0.02)` | `_mynard2010.py` | Matlab line 60-63 comment only: "avoids infinite C when FlowRatio approaches zero". **Not in the paper.** | **numerical regulariser presented as physics**; the 0.02 is a knob |
| `joining_etransfer_alpha = 0.2` | `mpce_v2_element.py` `DEFAULT_JOINING_ETRANSFER_ALPHA` | combaero invention (docstring says so); fitted by `tmp/calibrate_etransfer_join.py` to Bassett K11_corr/K12_corr + Idelchik at psi 1.25-3.33, "validated" on held-out measured points; memory note: "calibrated on pre-#212 degenerate plumbing -- revisit" | **stale calibration**; never validated in-network |
| Pseudosupplier angle via `atan2` of flow-weighted sin/cos | Mynard Eq 34 | grounded | fine, but non-analytic for complex-step |
| `+pi` flip and `theta -> -theta` mirror | `_mynard2010.py` lines 91-102 | Matlab angle normalisation (Fig 5 flowchart "reorient branch angles") | **implementation, not physics**; creates derivative discontinuities |
| `K` conversion restricted to `n <= 3` | `_mynard2010.py` | Matlab line 65 | design limit of the *conversion*, not the model |
| `Pt`-based residual instead of Mynard's static `C_j` form | `MPCEv2Element.residuals` | combaero choice ("ITERATION-2" comment) | the choice that imports the n <= 3 limit |

Three items are fudge-shaped in the sense the model-provenance discipline
means: the damping knob, the alpha term, and -- less severely -- `eta`, which
is real physics but fitted in a regime combaero does not operate in.

## 4a. Policy for tuned corrections, and where the data does not reach

Decision (user, 2026-09-04): **any empirical or CFD-derived correction must
prove itself against validation data, be labelled and documented as tuned,
and the validation data must cover the range it was tuned over.** Compressible
corrections are held to the same bar. This turns sec 4's audit into a
coverage check:

| constant | tuned over | in-network data covering that range | scored today | gap |
|---|---|---|---|---|
| `joining_etransfer_alpha = 0.2` | psi 1.25-3.33, theta {30, 45, 90}, joining type 6 (Bassett K11_corr/K12_corr + Idelchik anchors) | **Idelchik K11/K12: theta 30/45/90 x psi 1, 1.25, 1.67, 2.5, 3.33, 5, 10** (36 files) -- covers the range and beyond; Bassett Fig 10c/11b psi 1-4 at 45 deg | **yes, since step 0** (#282) | none in coverage. **Third provenance bug found at step 0:** `tmp/calibrate_etransfer_join.py` feeds `bassett2001.K11_corr(q, ...)` the lateral fraction, but Bassett K11 takes the straight one -- alpha was fitted against a mirrored Bassett K11 anchor. Retune (checklist #3) fixes the script's axis first, then re-fits, then the alpha 0 vs 0.2 table |
| Mynard `eta`: a0=0.8, a1=-0.2 | Mynard Fig 4: dividing Y-junction, Re 1363-1817 blood flow -- **a different scope**. A CFD-derived correction is a tuning, and its proof is not reproducing the CFD it came from but improving agreement on the data for **our** regime | Bassett K5/K6 (dividing, engine Re) and Hager xi_t | **yes, since step 0** (#282): Hager xi_t 45/45, MAE 0.286 | none in coverage. Proof = before/after table, `eta` on vs off (a0=a1=0), on Bassett K5/K6 + Hager. If it does not improve the dividing cells it is a tuning that does not earn its place here |
| Matlab damping `1 - exp(-FR/0.02)` | not tuned; a regulariser for FR -> 0 | Bassett curves reach q ~ 0 and q ~ 1, i.e. FR -> 0 on one collector | yes | none in coverage. Proof = sensitivity of the endpoint K to the 0.02 knob; document as regulariser, not physics |
| Hager 3/4 deviation, Bassett CV, dividing-streamline `p*` | derived, not tuned (three independent papers) | Bassett K5/K6/K11/K12, Hager xi_t | partially (Hager unscored) | none once Hager is scored |
| compressible `kappa M^2` (`K_dat_j_closed`, v3 spec sec 5) | derived analytically, **not yet proven against any data** | Wang 2014: joining, 45 deg, M 0.1-0.6, area ratio 1/1.56/2.44 (30 files); Perez-Garcia: 90 deg, C1/C2/D1/D2, M* 0.15-0.7 (**not digitised**) | **no** | **coverage gap**: no compressible dividing data at 45 deg, no compressible psi != 1 at 90 deg; Perez-Garcia must be digitised before dividing-flow compressibility can be claimed at all |

One of the five has a coverage gap that scoring alone cannot close: the
compressible dividing regime, which needs Perez-Garcia digitised from the
tables in `tier2_reference_data.md` before `kappa` can be claimed for
dividing flow. **Mynard's CFD is deliberately not reproduced** -- his scope
(blood flow, Re 350-2400) is not ours, and a CFD-fitted term is a tuning, so
the bar it has to clear is our validation data, not his figures.

**Retuning after a fix is expected, not optional.** A tuned constant was fitted
*with* whatever bugs existed at the time, so it may have absorbed part of one.
Fixing a bug -- physics or code -- therefore changes what the constant should
be, and the fix is not finished until the constant is re-fit and the impact
measured. Multiple bugs can cancel or dampen each other, which means the
intermediate state (one bug fixed, compensating tuning still in place) can
look *worse* than before. That is not evidence the fix was wrong; it is
evidence the follow-through is incomplete. Two live instances:

- `alpha = 0.2` was fitted on pre-#212 plumbing. Its retune on current
  plumbing is the measurement of what #212 changed, and belongs to step 1.
- The #272 gate removal was falsified by K6 collapsing -- but the kernel
  comment states the `sin^2` gate and the cross-coupling term were designed
  *as a pair* to hit Bassett K6 at 90 deg. Removing one half of a matched
  pair and observing failure shows the two are **coupled**, not that the gate
  is correct. The follow-through -- gate off *and* the cross term re-derived
  from Bassett Eq 27 without the gate assumption -- was never run. #272's
  "do not remove the gate" stands as "do not remove it *alone*".

**Labelling.** Each tuned constant carries, at its definition: the source
figure/table it was fitted to, the range, the date, and the scorecard cell
that proves it. `alpha` and `damping` currently carry none of this;
`DEFAULT_JOINING_ETRANSFER_ALPHA` reads as a physics constant.

## 4b. Step 1.3 result: the alpha x eta table (2026-09-05)

Every (alpha, eta_scale) combination scored on every digitised cell, imposed_q
topology, per source. Prediction written down first and confirmed: the two
terms are structurally separable in this dataset -- `eta` moves only dividing
cells (its `(1 - lambda)` factor is zero for a single collector), `alpha`
moves only psi > 1 joining cells -- to 4 decimals. The baseline `(0.2, 1)`
reproduces step 0 to the digit.

**`eta` -- keep at (0.8, -0.2).** On/off at fixed alpha:

| cell | eta off | eta on |
|---|---|---|
| Hager xi_t (45 pts) | MAE 0.3385, bias +0.338 | **0.2859, +0.055** |
| Bassett K5 | 0.3812, +0.235 | **0.3440, +0.033** |
| Bassett K6 | 0.0700, -0.047 | **0.0526, -0.002** |

Improves every dividing cell of both independent sources, chiefly by
removing a large positive bias -- without it the closure over-predicts
straight-port loss by a third of a dynamic head. A CFD-fitted term from
another regime that nonetheless earns its place on ours.

**`alpha` -- keep at 0.2.** Swept at eta on, joining sources:

| alpha | Bassett K11 | Bassett K12 | Idelchik K11 | Idelchik K12 | Idelchik all |
|---|---|---|---|---|---|
| 0.00 | 0.2268 (-0.159) | 0.1128 (+0.020) | 0.3095 (-0.309) | 0.6561 (-0.543) | 0.5116 |
| 0.10 | 0.1980 | 0.1025 | 0.1851 | 0.5599 | 0.4037 |
| 0.15 | 0.1887 | **0.1018** | **0.1611** | 0.5366 | 0.3801 |
| **0.20** | 0.1803 (-0.027) | 0.1032 (+0.048) | 0.1646 (-0.025) | **0.5287** (-0.274) | **0.3770** |
| 0.25 | 0.1738 | 0.1060 | 0.1744 | 0.5346 | 0.3845 |
| 0.30 | **0.1708** (+0.039) | 0.1094 | 0.2107 (+0.116) | 0.5539 | 0.4109 |

In-network optimum 0.2, shallow over [0.15, 0.25]; the sources pull slightly
apart (Bassett K11 prefers ~0.25, Bassett K12 and Idelchik K11 ~0.15). A
single scalar times area asymmetry cannot satisfy both angles at once --
that is the form's limit, not a tuning question.

**The anchor refit lost to the data.** #283's corrected-axis fit proposed
0.28-0.31; in-network 0.3 is *worse* than 0.2 (Idelchik +0.034, and it
over-corrects K11 at 30 deg badly). The anchors are Bassett's analytical
curves -- including his own K11@90 anomaly (sec 5) -- while the table is
measured and tabulated data. Checklist #1: the data decides. The calibration
script remains the "anchor on analytical" half of the method; it does not
set the value.

**What no alpha fixes:** Idelchik K12 keeps a bias of -0.27 at the MAE
minimum. That is the K12@90 under-prediction from step 0, growing with psi,
and it is a closure-shape limit (and partly Idelchik-vs-Bassett source
disagreement), not something the joining correction can reach.

Pinned executably by `python/tests/test_junction_tuned_constants_proof.py`:
each term improves its own block by a margin and does not touch the other.

## 4c. Step 1.4 result: the damping regulariser (2026-09-05)

`damping = 1 - exp(-FlowRatio / tau)` on Mynard's collector `C`, from the
Matlab reference only ("avoids infinite C when FlowRatio approaches zero").
Swept tau over {off, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2} on all 150 dividing
points (Bassett K5/K6, Hager xi_t) at the proven eta and alpha. A regulariser
is not tuned to data; the task was to show it does not matter to it.

| tau | interior MAE (101 pts) | endpoint MAE (49 pts) | endpoint converged |
|---|---|---|---|
| off | 0.1437 | **0.3809** | 38 / 49 |
| 0.005 | 0.1437 | 0.3822 | 38 / 49 |
| 0.01 | 0.1437 | 0.3837 | 38 / 49 |
| **0.02** | 0.1437 | 0.3859 | 38 / 49 |
| 0.05 | 0.1445 | 0.3986 | 38 / 49 |
| 0.1 | 0.1491 | 0.4196 | 38 / 49 |
| 0.2 | 0.1595 | 0.4425 | 38 / 49 |

- **(a) inert in the interior** over (0, 0.02], to 4 decimals.
- **(c) value-insensitive** over that band; degrades monotonically above it.
- **(b) necessary for convergence -- no.** Convergence is identical with the
  knob off at every endpoint (the same 11 Bassett failures at every tau,
  which are the artifact-root guard at true q -> 0/1, not a C blow-up), and
  the endpoint MAE is slightly *better* off. The "infinite C" case occurs
  nowhere in the data.

**Why it is cancelled.** Probing the FlowRatio -> 0 limit directly: `C`
diverges as 1/FlowRatio without damping (-21 -> -2.2e5 from 1e-2 to 1e-6)
and damping bounds it (~ -10.9). But MPCEv2 uses `K`, and Mynard Eq 18
multiplies `C` by FlowRatio^2, so `K_lat -> 1` -- the correct dead-branch
limit -- identically with or without the knob. The regulariser protects a
quantity the production path never exposes. It *would* matter for Mynard's
native `C_j` residual (sec 8 step 3b), which is the honest reason to keep it
at the reference value rather than delete it.

**Verdict:** band (0, 0.02], retained at 0.02 for fidelity to the reference
and for the C-form option. Labelled as such; pinned by
`python/tests/test_junction_damping_band.py` with exact tolerances (a limit,
not a tuning).

**Edge found on the way:** at FlowRatio exactly 0, `K_lat = -1` with or
without damping -- a sign-flipped dead branch the knob never guarded. Pinned
by test so the degenerate-iterate work (step 2) removes the test when it
fixes the edge, rather than rediscovering it.

## 4d. The strict-mode artefact, and Bassett's own case (2026-09-05)

Found while characterising one of the solves the step-2 fallback could not
rescue. The validation adapter defaulted to `strict=True`; production
(`gui/backend/graph_builder.py`) builds `MPCEv2Element` with `strict=False`.
`strict=True` raises on any transient wrong-sign Newton iterate instead of
steering it back through the soft barrier, **and** disables the post-solve
direction verifier. Measured on the full scorecard:

| topology | strict=True | strict=False |
|---|---|---|
| imposed_q | 602 / 691 | 602 / 691 |
| mfb_two_pb | 587 | **629** |
| three_pb | 389 | **480** |
| all | 1578 | **1711** |

- The `imposed_q` gate tables (sec 4b, 4c) are unaffected and stand.
- Every all-topology convergence comparison made earlier was a strict-mode
  artefact, including "v2 converges far less than v1" -- v1 has no strict
  raise. Withdrawn; the provenance record is corrected.
- **The "mirror root" reading of Bassett Fig 7b `mfb_two_pb` q=0.8 is
  withdrawn (2026-09-05, same day).** It converges to q=0.857, not 0.8, and
  its extracted K of 3.4370 equals the model's own `imposed_q` K at q=0.857.
  It is the model's correct root at a shifted operating point; the 24%
  "error" was a comparison against Bassett at the target q. See
  [JUNCTION_OPERATING_POINT_271.md] for the measurement and for what the
  direction-only verifier does and does not need to catch.
- **The run-to-run variation is real, and its cause is now known.** The same
  case converges in 5 of 10 identical processes. It is `PYTHONHASHSEED`:
  `_propagate_pressure_guess` seeds its BFS from `list(set(...))` over node-ID
  strings, so the propagated x0 pressures depend on the process hash seed.
  Fixed hash seed gives perfectly reproducible outcomes, and different seeds
  give different roots. **Fixed** (defect 12 below): the scorecard is now
  1700/2073 under every seed tested, against 1700 or 1711 before. Every
  convergence count recorded in this document before 2026-09-05 carries that
  uncertainty.
- The adapter now defaults to `strict=False`. Never compare v1 and v2
  convergence without matching it.

**Bassett's own Fig 7b case** (theta=45, psi=3, M ~ 0.03 -- inside his
measured conditions) converges in `three_pb` at q=0.4/0.6/0.8 and in
`mfb_two_pb` at q=0.6/0.8. The low-lateral-fraction points (q=0.2/0.4) do
not, and the follow-up investigation showed **why: those boundary-value
problems have no solution.** The model's `K_lat - K_str` has a minimum of
0.551 for this geometry and the targets sit at 0.122 and 0.385. No seed
reaches a root that does not exist. They stay pinned as xfail targets in
`python/tests/test_bassett_fig7b_case.py`, now labelled infeasible rather
than solver-path failures, and they come off when #272 closes the
`K_straight` gap. **The `three_pb` assertions in that test were tautological
and have been removed** -- see the next section.

## 4e. The pressure-driven topologies score less than they appear to (2026-09-05)

Full record: [JUNCTION_OPERATING_POINT_271.md]. Three results that change what
the port can be gated on:

1. **`three_pb`'s K extraction is a tautology.** Every `Pt` is imposed by a
   boundary through a lossless connection and `_extract_K` normalises by the
   fixed `m_dot_ref`, so the extracted K is the target by construction. A
   deliberately wrong model (`eta_scale=3.0`, `imposed_q` K 2.88 -> 3.07)
   still reproduces the q=0.8 target to 0.05%. `three_pb` measures
   convergence only.
2. **The converged operating point drifts from the target q**, and the score
   compares against the paper at the target anyway: `three_pb` median drift
   0.134, `mfb_two_pb` 12.6% of converged rows more than 0.25 away. Scoring K
   at the **achieved** q is the missing instrument.
3. **Where the BC set is feasible it often has two roots**, both exactly
   reproducing the model at their own q, with the initial guess deciding.
   Uniqueness is a property of `K_lat - K_str` being monotonic, i.e. of
   #272, not of the solver.

Consequences for sec 8: the step-4 gate cannot include "reproduce the
`three_pb` K table", and the `mfb_two_pb` table is only meaningful once K is
scored at the achieved q.

## 5. What the papers settle about the K_straight question (#272)

**Second motivation, added 2026-09-05:** `K_straight`'s size at low q also
decides whether the pressure-driven problem has a unique solution. The
model's `K_lat - K_str` is U-shaped, which makes low-q targets infeasible
and mid-range targets doubly-rooted. A monotonic difference leaves one
operating point. See [JUNCTION_OPERATING_POINT_271.md].

Three independent sources say **negative `K_straight` in dividing flow is real
physics**, not a solver artefact:

- Hager Eq 8: `xi_t = q (q - 1/2)`; Fig 4a shows the measured data negative for
  q < 0.5; sec 3: "It arises mainly due to acceleration forces acting on
  particles in the lateral outflow domain. The effect might be compared with a
  water jet pump".
- Mynard Fig 6e-h (Type 1, inlet-to-straight): CFD negative at low lambda;
  sec 4.3: "slightly negative loss coefficients for the straight branch in
  diverging T junction flow".
- Bassett Fig 7c: K5 measured minimum -0.06 at q ~ 0.8 (straight fraction).

`MPCEv2` uses Mynard's `K`, which carries this. The v1
`MultiPortChamberElement` uses the gated impulse residual, which cannot. The
scored data agrees:

| K5 (straight, separating) | v1 bias | v2 bias |
|---|---|---|
| theta=45 psi=1 | +0.166 | -0.013 |
| theta=45 psi=3 | -0.017 | +0.174 |
| theta=60 psi=1 | +0.153 | +0.005 |
| theta=90 psi=1 | +0.182 | -0.022 |

**v2 reproduces the sign structure of K_straight; v1 does not.** The
"redesign the junction pressure representation" item parked on #272 is
therefore a v1-only problem, and the cheaper resolution is to retire v1 as a
solver element (keep it as the M -> 0 regression target it was always meant to
be) rather than build a second pressure representation into it.

Also from the papers, for calibrating expectations: **Bassett's own joining
model has a documented anomaly at 90 deg.** Fig 11a shows K11 measured flat at
~+0.5 while the corrected prediction runs 0.1-0.4; the text concedes the angle
correction "is having too strong an effect". Any closure in the Bassett lineage
inherits this, so K11 at 90 deg is not a clean acceptance target.

## 6. Digitised data: what exists vs what is proven against

This is the provenance step that decides whether a port is *verified* or just
*transcribed*. Inventory of `validation/junction/data/`:

| source | files | coefficients | regime | scored by a network adapter? |
|---|---|---|---|---|
| Bassett 2001 | 57 | K4, K5, K6, K7, K8, K9, K11, K12, xi | separating + joining, high-Re engine air | **K5, K6, K11, K12 only** (34 files) |
| Hager 1984 | 5 | xi_t measured (single + manifold laterals), envelope | separating straight, plane duct | **yes since #282** (2 measured files, 45 points; envelopes excluded) |
| Idelchik 1966 | 36 | K11, K12 tabulated at theta 30/45/90, psi 1..10 | joining, psi sweep | **yes since #282** (396 points). NB: Idelchik indexes *every* diagram on Q_b/Q_c, K11 included -- unlike Bassett. And Idelchik disagrees with Bassett at theta=90, psi=1 (q=0.8: 1.56 vs ~0.9); K12@90 misses are partly source-vs-source, so show both, never average |
| Wang 2014 | 30 | K_13, K_23 | joining, **compressible** M 0.1-0.6, 45 deg | **no** (Tier 2, no tests) |
| Perez-Garcia 2010 | 0 csv | correlations in `tier2_reference_data.md` | **compressible** 90 deg, M* 0.15-0.7 | **no** -- not even digitised |
| Mynard 2015 Figs 6-11 | 0 | CFD Ref3D, six flow types, Re 350-2400 | different scope; `eta` treated as a tuning | **by decision, not reproduced** |

Consequences:

- **Hager's data is the independent check on K_straight.** It is the one
  dataset that is not Bassett, is specifically the coefficient #272 is about,
  and it has been digitised and ignored. `equivalences.py` already knows its
  `q` is `1 - q_Bassett`.
- **Idelchik is the joining psi-sweep** -- 36 files across psi = 1 to 10 -- and
  it is exactly what `alpha = 0.2` was fitted to, then never checked
  in-network. Whether alpha helps or hurts at the network level is unknown.
- **Bassett K4, K7, K8, K9 (16 files) are digitised and unscored.** K7/K8 are
  joining type 4 -- the second joining topology -- and the adapters' joining
  path only builds type 6.
- **Nothing compressible is scored**, and MPCEv2 has no compressibility to
  score. Wang is digitised; Perez-Garcia needs digitising from the tables in
  `tier2_reference_data.md`.
- **Mynard's own CFD is not digitised, and by decision will not be.** His
  scope is blood flow at Re 350-2400; the `eta` term it produced is a tuning
  whose bar is our engine-Re data, not a reproduction of his figures.

## 7. Bugs found, fixable now, independent of the port

| # | defect | evidence | fix | cost |
|---|---|---|---|---|
| 1 | ~~`N > 3` returns finite residuals with an all-zero `dKQ` block~~ | **fixed, step 2.** `MPCEv2Element.__init__` raises for N > 3 naming the Mynard 3-branch limit; ConstantKTee inherits it. The strict xfail is retired (the state is unconstructible); the closure-level precondition test stays | done |
| 2 | ~~Residual-level silent physics switch: `mpce_v2_element.py` catches *any* exception from Mynard and returns a **lossless** continuity residual with an **empty** Jacobian `{}`~~ | **fixed, step 2.1.** Measured first: 0 hits across 2073 scorecard records + the full suite (the pre-check guard ahead of the call fired 145 times -- that is the legitimate degenerate path). All three wide `except`s (residual, FD loop, diagnostics) narrowed to `(IndexError, ValueError)`; residual and FD paths raise a named `RuntimeError` with the cause chained, diagnostics annotates `closure_error`. Falsified 8/9 against pre-fix; scorecard identical to the digit | done |
| 3 | Hager and Idelchik unscored | sec 6 | extend `MPCEv2Network.evaluate_network` to `hager1984` (separating, `q -> 1 - q`) and `idelchik1966` (joining, `q` = Bassett's) | small; data + `q_transform` exist |
| 4 | `alpha = 0.2` calibrated pre-#212 | memory + `tmp/calibrate_etransfer_join.py` | re-run the calibration on current plumbing **after** #3 lands, so its validation is in-network; if it no longer earns its place, set the default to 0 | script exists |
| 5 | `damping` 0.02 undocumented as a regulariser | Matlab comment | name it (`FLOW_RATIO_DAMPING = 0.02`) with the Matlab citation and the sentence "numerical regulariser, not physics" | trivial |
| 6 | No compressibility, no stated limit | sec 3 | short-term: document `Ma < 0.2` on the element as `tee_junction.h` does; long-term: sec 8 step 5 | trivial now |
| 7 | Bassett K4/K7/K8/K9 unscored | sec 6 | add a type-4 joining builder and a K4 separating case to the adapter | medium |
| 8 | `three_pb`'s K score is a tautology (sec 4e) | `eta_scale=3.0` still reproduces the q=0.8 target to 0.05% | score K at the achieved q, or demote `three_pb` to a convergence-only row in the scorecard | small |
| 9 | K is scored at the TARGET q while the solve sits at the achieved q | drift: `three_pb` median 0.134; `mfb_two_pb` 12.6% of rows > 0.25 | return the converged q from `solve_and_extract` and score against the paper there | small |
| 10 | `verify_solution_consistent` passes near-degenerate roots | a root at q=0.0102 (lateral leg at 1% of the flow) is reported as a success; its only criterion is `m_dot > 1e-6` | add a minimum-flow-fraction criterion | small |
| 11 | `analytical_pt_prop` seeds no mass flow through `LosslessConnectionElement` | x0 violates continuity at every junction (0.1 in, 0.2 out); a junction-aware mass-conserving seed gains ~70 solves | extend the seed, but land it WITH item 9 -- it moves 22 already-converged roots between the two legitimate branches | medium |
| 12 | ~~Solver results depend on `PYTHONHASHSEED`~~ **fixed** | `_propagate_pressure_guess` seeds its BFS from `list(set(p_guess.keys()))`; the same case converges in 5 of 10 identical processes, and is deterministic per fixed hash seed | **done.** `queue = list(p_guess.keys())`; scorecard identical (1700/2073) under seeds 0/1/7/13, against 1700 or 1711 before | done |

Items 1, 2, 5, 6 are an afternoon. Item 3 is the one that changes what the
port can be held to. Item 12 is done; it blocked reproducible measurement of everything else.
Items 8-11 come from
[JUNCTION_OPERATING_POINT_271.md] and change what the port can be *gated*
on: 9 before 11, and 8 before either.

## 7a. Step 2.1 note: the wide-except audit (2026-09-05)

Every `except` in the junction model, adapters and Jacobian helper was
listed. Only three were wide, all in `mpce_v2_element.py`, all around
`junction_loss_coefficient`. `_network_builder.py`'s `except Exception` is
the correct pattern -- it reports `converged=False` *with the message*.
`_mpce_v2_jacobian.py:118` is narrow in type (`ZeroDivisionError`,
`FloatingPointError`) but wide in effect (a zero Jacobian entry); it belongs
with the degenerate-iterate work below.

**Measured, not assumed.** Instrumented all three and the pre-check guard,
then ran the full scorecard (3 sources x 3 topologies, 2073 records) and the
full suite: the three `except`s fired **0** times; the pre-check fallback
fired 145 times on the scorecard and once in the suite. The structural
argument -- the guard ahead of the call already routes every degenerate
split to the continuity fallback, so the closure cannot raise its two known
exceptions inside the `try` -- is confirmed. The handlers only ever caught
programming errors, and disguised them.

**What is now the open degenerate-iterate item.** The pre-check fallback
itself returns the continuity residual with `jac = {}` -- 145 times per
scorecard, on legitimate transient iterates. It works (the solves converge),
but an empty Jacobian for a live element is the same shape of thing at a
smaller scale, and it sits next to the `FlowRatio = 0 -> K = -1` edge from
step 1.4 and the `N > 3` rejection. Those three are the remainder of step 2.

## 7b. Step 2 result: the degenerate-iterate question (2026-09-05)

Measured on the full scorecard (2073 solves) before changing anything.

**(a) The pre-check fallback was where solves went to die.** Its 145 hits
came from **26 solves, median 5 hits each, 0 of 26 converged**. It returned
the continuity residual with `jac = {}` -- so even the trivial Jacobian it
owns (+-1 on the `Pt` rows, the port signs on the mass row) was withheld,
handing the solver N+1 zero rows. The iterates were all-same-sign with real
magnitude (e.g. `[0.0021, 0.0011, 0.0048]`, every port flowing out), which
violates mass conservation and is a wrong-direction state for at least one
port -- exactly what the existing soft barrier handles, except the pre-check
returned first. **Fix:** the all-ports-zero case (a cold start; never seen on
the scorecard) returns the soft barrier at zero slack, which *is* the
continuity residual with its Jacobian; the same-sign case now falls through
to the wrong-direction check. No new code path.

**(b) The mask mismatch was the more frequent silent lossless junction.**
384 hits. Example `[-0.1, -0.0, 0.125]`: MPCEv2 excludes the `-0.0` port
(`|U| <= 1e-9`), Mynard's `Q < 0` classifies it as a supplier, the two
disagree, `K` comes back mis-shaped, the shape guard trips, and `K_per_port`
stays all-zero. At what looks like initialization. **Fix:** an excluded port
is snapped to a tiny flow in its *declared* direction and the masks are
rebuilt from the snapped values, so the classifications agree by
construction and the `FlowRatio -> 0` limit proven in step 1.4 applies
(`K -> 1` for a dead outlet). After the fix the shape guard trips **0**
times on the scorecard.

**(c) N > 3** is refused at construction. Nothing outside the strict xfail
ever constructed one.

**Verified.** 9 of 10 new tests red against the pre-fix element (the 10th
pins "three ports still accepted"); full suite green; the scorecard
**identical to the digit** -- 1578/2073 converged, every imposed_q cell
unchanged. That last number is the honest part: the 26 doomed solves are
still doomed. The empty Jacobian was wrong, but it was not what killed them
-- they are the true q-endpoint / artifact-root cases. What the fix buys is
correctness of what the solver is handed on every iterate, which the port
must preserve; it does not buy convergence here.

**Dead end recorded:** re-instrumenting the element for the after-measurement
and then reverting with `git checkout` reverted the fix as well. Commit
first, then instrument.

## 8. The port, sequenced by provenance

Each step has a gate that is a table against digitised data, not a green suite.

**Step 0 -- widen the gate (items 3, 7 above).** Score Hager and Idelchik
through the current Python MPCEv2 and record the tables on #271 next to the
Bassett ones. Until this is done the port can only be shown equivalent to a
model that was never checked against half its available evidence.

**Step 1 -- resolve the three unprovenanced constants (items 4, 5).**
Re-calibrate or drop alpha with the wider gate; name the damping knob. Do not
transcribe a constant into C++ that cannot be explained.

**Step 2 -- fix the silent failures (items 1, 2).** The port must not inherit a
`continue` or an `except Exception` that turns "no value" into zero.

**Step 3 -- decide the residual form.** Two options:

- (a) Keep the `Pt`-based `K` form. Faithful to current behaviour; equivalence
  gate is exact; keeps the `n <= 3` limit (with item 1's rejection).
- (b) Move to Mynard's static `C_j` form, `p_i - p_j = C_j rho u_j^2`. Lifts
  `n <= 3` natively, matches the paper's residual, and matches the v3 spec's
  residual set. Changes results, so the gate is the data tables, not
  equivalence.

Recommend (a) for the port itself and (b) as a separate, later change --
porting and changing the physics in one step makes the equivalence gate
useless.

**Step 4 -- the port.** Reuse in this order:

1. `DualN<N>`: add `dexp`, `dsin`, `dcos`, `datan2` (the closure needs all
   four; the helper has only `dsqrt`/`dpow`). `d atan2(y, x) = (x dy - y dx)
   / (x^2 + y^2)`.
2. Pseudosupplier: port `_mynard2010.py` lines 80-118 on `DualN`, **evaluating
   the two reorientation decisions on the primal value** (branch-on-primal).
   Within a branch the derivative is exact; at a flip boundary it is one-sided,
   which is correct and is better than a zero column. Document this in the
   header. (The v3 spec's alternative -- pseudodatum in `double`, analytic
   outer chain -- is equivalent and already written; pick one, do not mix.)
3. `C_j` and `K` conversion: ~10 lines on `DualN`.
4. Whole-element `(f, J)`: `mpce_v2_residuals_and_jacobian` seeded over
   `(P_i, Pt_i, mdot_i)` per port plus `P_jct`, returning the `N+1` rows and
   the full Jacobian. Follow `ejector_element_residuals_and_jacobian`'s
   shape. The Python element becomes a relabelling shim, exactly as the
   ejector did.
5. Delete `_mpce_v2_jacobian.py` and the FD fallback. The sympy derivation
   stays as an offline cross-check of the canonical case (the ejector recipe),
   not as a runtime path.

**Step 4 gate:** the `imposed_q` Bassett, Hager and Idelchik tables
reproduced to solver tolerance with no worse convergence. **Not** the
`three_pb` K table, which any model reproduces (sec 4e), and the
`mfb_two_pb` table only once K is scored at the achieved q (item 9);
until then `mfb_two_pb` contributes convergence counts, not accuracy.
Plus whole-row FD tests
(`test_mpce_v2_jacobian_rows.py`, extended to joining and to `N = 4` once
item 1 is lifted) at `< 1e-6`; a ctest comparing the C++ kernel against golden
Python values (the ejector recipe's `.h` golden-data pattern).

**Step 5 -- compressibility.** Port `K_dat_j_closed`'s `kappa M_dat^2`
correction into the closure (it is already in `tee_junction.h`, templated, with
the v3 spec's derivation). Gate: Wang 2014 (digitised) and Perez-Garcia
(digitise first). This is the step that restores what the retired tee had.

**Step 5a -- operating-point instrumentation (items 8-11).** Independent of
the port and a prerequisite for gating it: return the converged q, score
against the paper there, demote the `three_pb` K column, tighten the
verifier, then land the junction-aware seed. Doing the seed first would
silently move 22 roots with nothing able to tell which branch is wanted.

**Step 6 -- retire v1 as a solver element** (sec 5). Keep
`MultiPortChamberElement` as the M -> 0 regression target; close the #272
pressure-representation item as superseded.

## 9. Decisions

Taken (user, 2026-09-04):

1. **Residual form: stepwise and tested.** Port the faithful `Pt`-based form
   first with an exact equivalence gate; move to Mynard-native `C_j` as a
   separate, later, data-gated change. Never both in one step.
2. **Tuned corrections** (`alpha`, `eta`, damping, and any future one) prove
   themselves against validation data for our regime, and are labelled as
   tuned. CFD-derived corrections count as tuning; the CFD they came from is
   **not** reproduced. Sec 4a is the coverage check; its one remaining gap
   (compressible dividing flow) is a digitisation job that precedes `kappa`.
3. **Compressible corrections** are appealing and held to the same bar --
   Wang plus a digitised Perez-Garcia before `kappa` is claimed.

Still open:

4. **What to do with `alpha`** if, once Idelchik is scored, it does not earn
   its place: drop (provenance-clean default) or re-fit on current plumbing.
5. **Retiring v1** as a solver element (sec 5).

## References on disk (`docs/junction/`, gitignored, present)

- `bassett.pdf` -- Table 1 (definitions, q conventions), Eqs 13-15 (K5),
  22-27 (K6, xi), 28-34 (joining, angle corrections), Figs 7, 9-13.
- `hager1984.pdf` -- Eqs 1-8 (xi_t from momentum + dividing streamline),
  14 (contraction), 19-20 (xi_l), Fig 4 (xi_t data), sec 3 (jet-pump
  explanation of negative xi_t).
- `a-unified-method-...pdf` (Mynard & Valen-Sendstad 2015) -- Eqs 13-18
  (K vs C), 19-30 (CV analysis), 31-36 (pseudodatum, eta), Figs 4, 6-11
  (CFD reference), sec 4.5.1 (Re range argument).
- `compressible_junction_model_v3.pdf` -- sec 5 (closed O(M^2) K), 8
  (residuals + Jacobian tables), 9 (C++ fragments), 11 (limitations).
- `JunctionLossCoefficient.m` -- line 60-63 (damping), 65 (n <= 3 K).

[JUNCTION_OPERATING_POINT_271.md]: ../../docs/archive/JUNCTION_OPERATING_POINT_271.md
