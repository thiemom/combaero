# MPCE Junction: C++ Port Design and Provenance Audit (issue #271)

Deep-dive read of the three source papers on disk (Bassett 2001, Hager 1984,
Mynard & Valen-Sendstad 2015), the internal compressible spec
(`compressible_junction_model_v3.pdf`), the Matlab reference
(`JunctionLossCoefficient.m`), and the code. Written to decide *how* to port
`MPCEv2Element` to C++, what has to be fixed before that is safe, and what the
digitised data can and cannot prove.

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
| `joining_etransfer_alpha = 0.2` | psi 1.25-3.33, theta {30, 45, 90}, joining type 6 (Bassett K11_corr/K12_corr + Idelchik anchors) | **Idelchik K11/K12: theta 30/45/90 x psi 1, 1.25, 1.67, 2.5, 3.33, 5, 10** (36 files) -- covers the range and beyond; Bassett Fig 10c/11b psi 1-4 at 45 deg | **no** | none in coverage; entirely in scoring. Proof = before/after table, alpha=0 vs 0.2, on Idelchik + Bassett psi sweeps |
| Mynard `eta`: a0=0.8, a1=-0.2 | Mynard Fig 4: dividing Y-junction, Re 1363-1817 blood flow -- **a different scope**. A CFD-derived correction is a tuning, and its proof is not reproducing the CFD it came from but improving agreement on the data for **our** regime | Bassett K5/K6 (dividing, engine Re) and Hager xi_t | partially (Hager unscored) | none in coverage. Proof = before/after table, `eta` on vs off (a0=a1=0), on Bassett K5/K6 + Hager. If it does not improve the dividing cells it is a tuning that does not earn its place here |
| Matlab damping `1 - exp(-FR/0.02)` | not tuned; a regulariser for FR -> 0 | Bassett curves reach q ~ 0 and q ~ 1, i.e. FR -> 0 on one collector | yes | none in coverage. Proof = sensitivity of the endpoint K to the 0.02 knob; document as regulariser, not physics |
| Hager 3/4 deviation, Bassett CV, dividing-streamline `p*` | derived, not tuned (three independent papers) | Bassett K5/K6/K11/K12, Hager xi_t | partially (Hager unscored) | none once Hager is scored |
| compressible `kappa M^2` (`K_dat_j_closed`, v3 spec sec 5) | derived analytically, **not yet proven against any data** | Wang 2014: joining, 45 deg, M 0.1-0.6, area ratio 1/1.56/2.44 (30 files); Perez-Garcia: 90 deg, C1/C2/D1/D2, M* 0.15-0.7 (**not digitised**) | **no** | **coverage gap**: no compressible dividing data at 45 deg, no compressible psi != 1 at 90 deg; Perez-Garcia must be digitised before dividing-flow compressibility can be claimed at all |

One of the five has a coverage gap that scoring alone cannot close: the
compressible dividing regime, which needs Perez-Garcia digitised from the
tables in `tier2_reference_data.md` before `kappa` can be claimed for
dividing flow. **Mynard's CFD is deliberately not reproduced** -- his scope
(blood flow, Re 350-2400) is not ours, and a CFD-fitted term is a tuning, so
the bar it has to clear is our validation data, not his figures.

**Labelling.** Each tuned constant carries, at its definition: the source
figure/table it was fitted to, the range, the date, and the scorecard cell
that proves it. `alpha` and `damping` currently carry none of this;
`DEFAULT_JOINING_ETRANSFER_ALPHA` reads as a physics constant.

## 5. What the papers settle about the K_straight question (#272)

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
| Hager 1984 | 5 | xi_t measured (single + manifold laterals), envelope | separating straight, plane duct | **no** -- `paper != "bassett2001"` returns unsupported |
| Idelchik 1966 | 36 | K11, K12 tabulated at theta 30/45/90, psi 1..10 | joining, psi sweep | **no** -- used only by the alpha calibration script |
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
| 1 | `N > 3` returns finite residuals with an all-zero `dKQ` block | pinned by strict xfail in `test_mpce_v2_fd_fallback_guards.py` | reject `N > 3` at construction with a message naming the Mynard 3-branch K limit; the C-form is the eventual lift | 1 line + test |
| 2 | Residual-level silent physics switch: `mpce_v2_element.py:402` catches *any* exception from Mynard and returns a **lossless** continuity residual with an **empty** Jacobian `{}` | same pattern as the FD guards, one level up: an unrepresentable state becomes a plausible-looking lossless junction mid-solve | narrow the `except` to the two exception types Mynard actually raises (`IndexError`, `ValueError` -- pinned in the guard tests), and either raise or return the soft-barrier residual, never `{}` | small |
| 3 | Hager and Idelchik unscored | sec 6 | extend `MPCEv2Network.evaluate_network` to `hager1984` (separating, `q -> 1 - q`) and `idelchik1966` (joining, `q` = Bassett's) | small; data + `q_transform` exist |
| 4 | `alpha = 0.2` calibrated pre-#212 | memory + `tmp/calibrate_etransfer_join.py` | re-run the calibration on current plumbing **after** #3 lands, so its validation is in-network; if it no longer earns its place, set the default to 0 | script exists |
| 5 | `damping` 0.02 undocumented as a regulariser | Matlab comment | name it (`FLOW_RATIO_DAMPING = 0.02`) with the Matlab citation and the sentence "numerical regulariser, not physics" | trivial |
| 6 | No compressibility, no stated limit | sec 3 | short-term: document `Ma < 0.2` on the element as `tee_junction.h` does; long-term: sec 8 step 5 | trivial now |
| 7 | Bassett K4/K7/K8/K9 unscored | sec 6 | add a type-4 joining builder and a K4 separating case to the adapter | medium |

Items 1, 2, 5, 6 are an afternoon. Item 3 is the one that changes what the
port can be held to.

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

**Step 4 gate:** all Bassett, Hager and Idelchik tables reproduced to solver
tolerance with no worse convergence; whole-row FD tests
(`test_mpce_v2_jacobian_rows.py`, extended to joining and to `N = 4` once
item 1 is lifted) at `< 1e-6`; a ctest comparing the C++ kernel against golden
Python values (the ejector recipe's `.h` golden-data pattern).

**Step 5 -- compressibility.** Port `K_dat_j_closed`'s `kappa M_dat^2`
correction into the closure (it is already in `tee_junction.h`, templated, with
the v3 spec's derivation). Gate: Wang 2014 (digitised) and Perez-Garcia
(digitise first). This is the step that restores what the retired tee had.

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
