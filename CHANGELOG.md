# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- **MPCE networks now default to ``analytical_pt_prop`` initialization.**
  ``NetworkSolver.solve(init_strategy="default")`` auto-upgrades to
  ``analytical_pt_prop`` when the network contains a
  ``MultiPortChamberElement``; passing any other strategy or an
  explicit ``x0`` opts out, and the strategy actually used is recorded
  in ``solver_settings_used``. Evidence: on the certified inverse-design
  audit (``tmp/mpce_audit_v2_generator.py`` + ``_runner.py``, 32
  feasible fixtures with stored ground-truth roots), analytical
  converged 32/32 to the certified root and was the fastest strategy
  overall (41 s total), vs 28/32 for the plain cold start and
  homotopy. The GUI ``SolverSettings`` schema now accepts
  ``analytical_pt_prop`` explicitly. NOTE: the earlier LHS-32 audit
  numbers (complementary default/analytical partition) are void -- that
  fixture family was later shown to be infeasible for any physical tee
  closure (equal-Pt sinks + constant-area legs demand static recovery
  the model set cannot represent).

### Deprecated
- **``init_strategy="incompressible_warmstart"``** now emits a
  ``DeprecationWarning``. On the certified audit it converged 10/32 at
  roughly 10x the wall time of ``analytical_pt_prop`` (409 s vs 41 s
  total); on a real GUI tee network it needed 445 evaluations vs 31
  for the cold start. The incompressible proxy seeds a low-residual but
  poorly conditioned region for the compressible restart. Will be
  removed once the GUI drops the option.

### Added
- **GUI: all supported init strategies selectable in the sidebar.**
  ``analytical_pt_prop`` added to the strategy dropdown (the backend
  accepted it since the auto-upgrade PR but the frontend never offered
  it); ``incompressible_warmstart`` is labelled "(deprecated)".
- **GUI: the solver status badge shows the init strategy that actually
  ran** (from ``solver_settings_used``), flagged "(auto)" when the
  solver upgraded the requested "default" -- previously the MPCE
  auto-upgrade was invisible from the GUI.

### Changed
- **GUI: "Reset All Initial Guesses" is now "Reset Guesses & Results"
  and truly resets.** It previously cleared only
  ``node.data.initial_guess`` -- but the graph builder warm-starts from
  ``node.data.result`` whenever ``initial_guess`` is empty, so the old
  button silently RE-ENABLED result-based seeding instead of giving a
  cold start. It now clears stored solve results as well; the next
  solve behaves like a freshly built network.
- **GUI: ``mpce_tee`` is displayed as "Tee".** Sidebar drag entry, node
  card default label, and Inspector header now say "Tee"; the Inspector
  section carries a short explainer naming the underlying model
  (momentum-based multi-port chamber junction, MPCE / Mynard closure).
  The internal node type and save-file format are unchanged.

### Fixed
- **The MPCE ``analytical_pt_prop`` auto-upgrade now requires >= 2
  ``PressureBoundary`` nodes.** The strategy's Pt propagation needs a
  boundary pressure gradient to be meaningful. On MFB-driven networks
  with a single PB (the common GUI pattern: mass inlet -> tee tree ->
  one pressure sink), propagated per-node Pt differences are
  path-length artifacts and the Bernoulli sign logic seeded REVERSED
  channel m_dots -- parking the first iterate inside the junctions'
  soft-barrier penalty region (|F(x0)| ~ 1.15e6 of pure penalty on a
  real 5-tee GUI network, hybr immediately "not making good progress")
  while the plain cold start converges in 16 evaluations. Such
  networks now keep the plain default initialization. The certified
  audit result behind the auto-upgrade is unaffected (all 32 fixtures
  are PB-driven).
- **GUI ``mpce_tee`` junctions are now built with ``strict=False``.**
  Strict mode raises whenever a Newton ITERATE (not the converged
  answer) explores a wrong-sign port m_dot, killing GUI solves that
  would otherwise converge with "Unexpected error during residual
  evaluation: MPCEv2Element ... declared flow_direction=... but
  observed wrong flow direction". Soft mode lets the solver explore
  wrong-direction states transiently; the post-solve
  ``verify_solution_consistent`` guard (see entry below) rejects any
  wrong-direction artifact root at convergence, so this switch does
  not admit silently wrong answers.
- **Soft-barrier artifact roots are demoted to honest failures.**
  ``MPCEv2Element`` with ``strict=False`` replaces a wrong-direction
  port's physics row with Pt continuity plus a one-sided penalty
  ``alpha * mdot^2`` -- and that penalty can exactly cancel the
  continuity mismatch, so the solver converges (|F| ~ 1e-10) to a state
  with a slightly reversed feed that the strict residual rejects
  (found via constant-K warm starts on the certified audit: a merge
  feed reversed to -0.0035 kg/s while the certified all-forward root
  has it at +0.213). ``NetworkSolver.solve`` now runs the element's
  ``verify_solution_consistent`` post-solve and reports such solutions
  as non-converged with an explanatory message, instead of returning a
  silently wrong flow split. This is the prerequisite for building GUI
  junctions with ``strict=False`` (Newton may explore wrong-direction
  states transiently without dying, and artifact roots can no longer
  masquerade as answers).

### Removed
- **``tee_junction`` (Bassett K-closure) retired from the GUI** (#177).
  ``MPCE Tee`` (``mpce_tee``, momentum-CV closure) is now the sole
  user-drawable tee type. Removed the Sidebar drag-create entry,
  ``TeeJunctionNode.tsx`` node card, ``TeeJunctionData`` schema,
  ``graph_builder`` dispatch, ``runner`` state projection, and all
  Inspector conditionals + display fields. The ``TeeJunctionElement``
  Python class + Python tests + validation adapters
  (``validation/junction/models/tee_junction_element_network.py``,
  ``tee_junction_raw.py``) are retained as the M -> 0 regression
  baseline for the junction validation suite -- the K-closure is no
  longer a competing user-facing model but remains authoritative for
  regression against Bassett 2001. GUI networks saved with a
  ``tee_junction`` node no longer load; they need to be recreated with
  ``mpce_tee``.

### Added
- **``init_strategy="analytical_pt_prop"`` on ``NetworkSolver.solve``**.
  Opt-in initialization strategy targeting MPCE compressible cold-start
  cases where the solver's default per-element uniform-dP propagation
  lands Newton in the wrong basin. Seeds two element unknowns the
  default path leaves at arbitrary constants:
  ``ChannelElement.m_dot`` from a Bernoulli estimate on the propagated
  Pt drop (default fallback is ``ref["m_dot"] = 0.1``), and
  ``MultiPortChamberElement.P_jct`` from the propagated Pt at the
  junction's "common" port (default fallback is ``min(boundary Pt)``,
  almost never near the true junction static). User-provided
  ``initial_guess`` entries take precedence; single-shot semantics
  (overrides cleared after ``_build_x0``). LHS-32 audit
  (``tmp/mpce_cold_start_audit.py``) partitions cases where the new
  strategy is complementary to ``default``: high-theta (> 50 deg)
  topologies favour analytical, low-theta favour default. Combined
  success rate 21/32 = 66% vs 14/32 = 44% for either alone; ship as an
  alternative attractor, not a default replacement. Follow-up work:
  fallback/auto strategy that tries multiple starts.
- **``__convergence_history__`` field on NetworkSolver solution dict**.
  Previously the per-iteration residual-norm trace lived only on
  ``solver._diagnostic_data`` (consumed by the GUI API path) and was
  invisible to plain ``solver.solve()`` callers. Now also emitted as
  ``sol["__convergence_history__"]`` so validation scripts / Python
  notebooks can plot the |F| vs iteration curve without reaching into
  solver internals. Format: list of ``{eval, t, norm}`` dicts, throttled
  at 0.1s intervals (or whenever the iterate improves), capped at 1000
  points. Both sources read from the same ``_history`` accumulator;
  this just exposes it via the public sol dict too.
- **CSV download button on ``ConvergenceChart`` GUI modal**. The chart
  already showed live the residual-norm trace + worst-residuals table
  after each solve; users can now click a download icon to save the
  raw trace as ``convergence-{timestamp}.csv`` (columns ``eval,t_s,norm``,
  plus the worst-residuals list as a comment block at the end). Lets
  you grab a slow-convergence trace, close the modal, and analyse
  later without re-solving.

### Changed
- **MPCE-v2 default Jacobian method flipped from ``"fd"`` to ``"sympy"``**
  (``MPCEv2Element.jacobian_method``). The sympy analytical Jacobian for
  the canonical 3-port separating T has been wired since the prototype
  PR but was kept disabled because the iteration-4 commit message
  noted FD's tiny truncation acting as regularization on some hard
  low-q mfb_two_pb cases. After the soft-barrier iteration and PR #197
  joining-side calibration shipped, that empirical concern no longer
  shows up: a focused audit comparing both methods across the full
  separating measured grid (Bassett K6 + K5, 105 points x 3 topologies)
  shows IDENTICAL convergence counts (94/32/11 of 105 across imposed_q
  / three_pb / mfb_two_pb) and IDENTICAL mean/max K errors, with
  ~1.14x wall-time speedup from cutting 4 Mynard evaluations to 1 per
  Newton step on the canonical case. Joining flow and non-canonical
  geometries (theta != 0 on straight, N != 3, port 0 not supplier)
  auto-fall-back to FD via the existing canonical-check, so the flip
  only affects separating canonical cases. Set ``MPCEv2Element.jacobian_method = "fd"``
  on a per-instance basis to revert.

### Added
- **K diagnostics for MPCE-v2 tee element**. ``MPCEv2Element.diagnostics``
  re-evaluates Mynard at the converged state and emits per-port K values
  plus topology-aware named aliases:
    - separating (``flow_direction='branch'``): ``K_straight``, ``K_branch``,
      ``mass_flow_ratio = m_dot_branch / m_dot_com``
    - joining (``flow_direction='merge'``): ``K11``, ``K12``,
      ``mass_flow_ratio`` (same convention as Bassett joining-T q)
  Plus convenience ``Pt_jct`` field (renames the legacy ``P_jct`` slot
  which already stores Pt in MPCE-v2). GUI Inspector / save files can now
  show K coefficients in the same shape that ``TeeJunctionElement``
  emits, enabling side-by-side comparison and validation against Bassett
  / Hager / Idelchik reference data without back-calculating from port
  pressures. Falls back to the parent fields only when ``port_mdots`` is
  not provided (legacy callers).
- Per-port ``port_{i}_m_dot`` field on ``MultiPortChamberElement.diagnostics``
  when the solver supplies ``port_mdots``.

### Fixed
- **``L_choke`` is interpolated within the breaking march step instead of
  snapping to the RK4 grid.** ``fanno_channel`` / ``fanno_channel_rough``
  reported the length-to-choking quantized to ``L / n_steps`` (the
  position of the last completed step). The in-channel choke barrier
  (see the barrier entry below) builds its ramp on ``L_choke``, so the
  channel's dP(m_dot) closure was a staircase near choke onset: jumps of
  ``kappa * P_in / n_steps`` (~1 kPa at 1 bar with the default 100-step
  march) every ~0.008 kg/s. Any pressure-driven network whose residual
  zero fell inside a staircase jump had NO root -- a single air channel
  (L=1 m, D=0.05 m) between two pressure boundaries failed to converge
  at EVERY Pt ratio from 1.2 to 10 under every method/init combination,
  stalling at a tread edge (best |F| ~ 90 at ratio 1.2). The sonic
  station is now interpolated: linearly in Mach between the last two
  steps for the ``M >= kFannoChokeMach`` exit, and linearly in the
  pressure probe for the ``P <= 0`` RK4-stage exits. The same
  single-channel sweep now converges at all 13 sampled ratios (1.2-10),
  landing on the monotone barrier ramp at the choke boundary as the
  barrier design intended. The choke threshold Mach (0.999) is hoisted
  to ``constexpr kFannoChokeMach`` in ``compressible.h``. NOTE: for BCs
  beyond the subsonic envelope the converged ``m_dot`` sits on the
  barrier ramp slightly above the physical choked flow -- check the
  channel's ``choked`` flag when interpreting such solutions. This also
  invalidates (again) the pre-fix LHS-32 audit numbers: the 2026-07-03
  re-run showed 25/32 all-strategy failures, most of which trace to this
  staircase because the audit's channel geometry leaves the subsonic
  envelope at Pt ratio ~1.19.
- **MultiPortChamberElement ports now report their real throughflow to the
  solver's state propagation.** ``flow_at_node`` returned 0 "for safety",
  so every port-MCN fed by the junction (all collector/outflow ports)
  received zero-mass upstream streams. Three consequences, all fixed:
  (1) the collector-port MCN closure ``Pt = P + 0.5 rho v^2`` silently
  degenerated to ``Pt = P`` -- the outflow dynamic head (~29 kPa for
  0.4 kg/s through a 5 cm port at 1 bar) vanished from the bookkeeping,
  and MPCE-v2's Mynard residual received a static pressure where its
  design expects the port stagnation pressure; (2) merging streams with
  unequal temperatures mixed with all-zero weights, so the outlet
  temperature snapped to the FIRST inlet stream's value (a 0.3 kg/s @
  400 K + 0.1 kg/s @ 300 K merge reported 400 K instead of the
  enthalpy-weighted ~375 K) -- previously masked because all junction
  tests are isothermal; (3) the degenerate closure made non-physical
  network states exactly satisfiable (dead-branch / reversed-supply roots
  at |F| ~ 1e-7). Ports now resolve their outer connecting elements'
  ``m_dot`` unknowns via an index map stashed by ``NetworkSolver`` at
  unknown-vector build time; supplier ports and single-supplier
  collectors carry their own outer flow, multi-supplier collectors the
  sum of the supplier feeds, with matching ``flow_jac_at_node`` entries
  so the closure Jacobian stays consistent. Additionally,
  ``MultiPortChamberElement.resolve_topology`` pushes inferred port areas
  onto auto-sized port-MCNs (collector ports have no upstream channel to
  inherit ``Dh`` from and previously kept the 0.1 m^2 fallback, hiding
  the face velocity even when flows were known). Converged results for
  MPCE networks shift by the collector dynamic head; junction loss
  calibrations (e.g. MPCE-v2 ``joining_etransfer_alpha``) were fitted on
  the degenerate plumbing and should be revisited against the validation
  suite. Known consequence: with honest collector stagnation pressures,
  the v1 impulse rows (``sin^2(theta)`` per-port form, Bassett
  correspondence derived for separating flow) are energetically
  inconsistent for JOINING flow -- a v1 merge junction can manufacture
  flow work and may lose its all-forward root entirely.
  ``test_step_4_adding_bypass`` migrated its junctions to
  ``MPCEv2Element`` (the production model, joining-capable by design);
  the v1 vs K-closure comparison test now seeds the forward basin
  explicitly since v1's sign-symmetric residuals admit the mirrored
  all-reverse root from a cold start.
- **Choked compressible channels no longer report zero resistance to the
  network solver.** ``channel_compressible_mdot_and_jacobian`` passed
  through ``fanno_channel_rough``'s truncated march: a supersonic inlet
  (M >= 1) returned ``dP = 0`` with zero ``m_dot`` sensitivity, and
  in-channel choking (``L_choke < L``) collapsed the reported drop toward
  zero as the requested velocity approached sonic. The infeasible region
  therefore looked frictionless to Newton and acted as a flat attractor
  with exact spurious network roots -- e.g. the canonical 3-port MPCE
  branch tee (pt_ratio 2.0, theta 45 deg, psi 1.0) "converged" at
  |F| ~ 3e-7 to a dead branch arm with the straight arm at M ~ 1.1 and
  zero total-pressure drop. The kernel now applies a monotone
  infeasibility barrier (``kChannelChokeBarrierKappa``,
  ``kChannelChokeBarrierBeta`` in ``solver_interface.h``): in-channel
  choking adds ``kappa * P_in * (1 - L_choke/L)`` on top of the truncated
  drop, a supersonic inlet returns ``kappa * P_in`` plus quadratic growth
  in velocity, keeping ``dP(u)`` strictly increasing across feasible,
  choked, and supersonic regions so Newton is always pushed back toward
  the feasible branch. Feasible (un-choked) results are bit-identical to
  before. The previously-spurious dead-branch state now evaluates at
  |F| > 1e4 (regression-tested); such cases honestly report
  non-convergence instead of silently returning a wrong flow split.
  Notable interplay with the port-throughflow fix (#212): the |F|~3111
  cold-start stall on audit case #4 was the flat region itself (with the
  barrier alone the default cold start converges the case in seconds),
  but the honest collector closure changes the v1 impulse merge's root
  structure so the default cold start stalls again -- while
  ``analytical_pt_prop`` now converges the case directly in ~1-2 s
  instead of the earlier ~110-eval marathon. The #210 regression test
  keeps its rescue semantics with the updated timings. The LHS-32
  audit's strategy comparison predates both fixes and should be re-run
  before designing the auto/fallback init strategy.
- **Merging/splitting MomentumChamberNode now refused at build time** (#174).
  MCN's scalar ``Pt = P + 0.5 rho v^2`` closure models a single bulk
  velocity and cannot represent merging or splitting streams within the
  chamber itself; the silent-corruption case at the time of #174 was
  ``test_step_4_adding_bypass``, where the face-convention proxy
  (channel -> Pt, orifice -> static) produced a non-physical 17 kPa
  asymmetry at a high-q merging MCN. ``FlowNetwork.validate`` now raises
  ``ValueError`` for any MCN with > 1 incoming or > 1 outgoing edges,
  pointing the user at ``MultiPortChamberElement`` (momentum-CV junction
  with one MCN per port) as the correct migration target. The
  formerly-xfailed ``test_step_4_adding_bypass`` is restored as a
  passing test on the MPCE topology and now exhibits the expected
  diameter-driven split (smaller-bore bypass < main branch). Pass-through
  MCNs (1-in / 1-out, including the GUI auto-insert pattern and MPCE
  port faces) remain valid.
- **MPCE Tee GUI shows ``m: 0.000 kg/s`` on the node card and in Live
  Telemetry** even when the junction is actually flowing. Root cause:
  ``MPCEv2Element.diagnostics`` only emitted per-port ``port_i_m_dot``
  values, so ``runner._build_result_objects`` (which derives
  ``ElementResult.m_dot`` from ``diag["m_dot_com"]``) fell through to
  ``0.0``. Now emits the topology-aware aliases ``m_dot_com``,
  ``m_dot_branch``, ``m_dot_straight`` (matching ``TeeJunctionElement``'s
  convention) for the canonical 3-port case. The Live Telemetry panel
  also annotates ``port_X_*`` keys with the arm letter (C/S/B) so the
  port index displayed matches the C/S/B labels on the node card.
- **Tee area inheritance** (both ``tee_junction`` and ``mpce_tee``): the
  schema docstring promised ``F_C / F_branch = None`` means "inherit from
  connected channel", but the graph_builder dispatch never implemented the
  walk. New ``_resolve_tee_areas`` helper inherits ``F_C`` from the
  common-arm channel and ``F_branch`` from the branch-arm channel,
  computes ``psi = F_C / F_branch``. Explicit user values still win;
  fallback to ``F_C = 0.01`` and stored ``psi`` when no channels are
  connected. Frontend ``NetworkCanvas`` drag-create defaults no longer
  bake in ``F_C=0.01`` / ``psi=1.0`` so the schema-level ``None`` default
  propagates and inheritance kicks in for new tees. Backwards-compatible:
  existing saved networks with explicit ``F_C: 0.01`` retain that value;
  users can clear the Inspector field to opt into inheritance.

### Added
- **Idelchik 1966 joining-T dataset** (`validation/junction/data/idelchik1966/`):
  396 raw tabulated K-coefficient cells from Section VII diagrams 7-1 through
  7-7 (converging wye with `Fs+Fb > Fc, Fs = Fc`, theta=30/45/90 deg, F_b/F_c
  in {0.1..1.0}, q in {0..1}). Two-pass independent transcription with
  formula cross-check; two confirmed Idelchik typos tagged in metadata.yaml.
  Mapping: Idelchik zeta_c.b -> K12, zeta_c.s -> K11; psi = 1/(F_b/F_c)
  matches Bassett convention. Reformatted into 36 per-(psi, theta) CSVs
  under the standard dataset naming. Used as second analytical cross-anchor
  for the MPCE-v2 joining-side calibration.
- **MPCE-v2 joining-side etransfer correction** (combaero extension to
  faithful Mynard 2010): `joining_etransfer_alpha: float = 0.2` on
  `MPCEv2Element`. Mynard's original etransfer collapses to zero for
  converging flow (the `(1 - flow_ratio)` factor); the asymmetric supplier
  areas in joining-T networks then leave Mynard underpredicting K_avg at
  psi > 1. The correction adds `alpha * (A_max - A_min) / (A_max + A_min)`
  to etransfer, vanishing at psi=1 (preserves the equal-area baseline) and
  growing with area asymmetry. Single parameter, calibrated against Bassett
  K11/K12 analytical + Idelchik 1966 tabulated values at psi in [1.25,
  3.33] and theta in {30, 45, 90} (analytical-only anchors, measured points
  held out as independent validation). Independent validation: imposed_q
  joining-flow audit gives mean error -10%, max error -14% vs faithful
  Mynard; three_pb and mfb_two_pb topologies essentially tied. Pass
  `joining_etransfer_alpha=0.0` for faithful Mynard, or a custom float for
  re-calibrated workflows. The default value is locked by a regression
  test (`test_mpce_v2_joining_etransfer.py`).
- **Calibration script** `scripts/calibrate_mpce_v2_etransfer.py`:
  reproduces the alpha=0.2 fit and reports stratified residuals by
  (theta, psi) for any future re-calibration.

### Changed
- `validation/junction/schema.py`: extended `Kind` Literal with
  `"tabulated"` (handbook semi-empirical values, distinct from `"calc"`
  analytical formulas and `"measured"` experimental points). Extended
  `QConvention` Literal with `"idelchik"` (equivalent to `"bassett"` for
  joining type-6).
- `validation/junction/models/mynard2010.py`: `junction_loss_coefficient`
  gains optional `joining_etransfer_alpha: float = 0.0` parameter. Default
  preserves faithful-port Mynard behavior; non-zero values enable the
  combaero joining-side correction described above.
- **Python bindings for all 12 Bassett K coefficients**: `cb._core.tee_K1`
  through `tee_K12`. Previously only K5, K6, K11, K12 were exposed; the
  rest required separate C++ tests. The new bindings let the
  `validation/junction` scorecard evaluate all 12 K against measured
  data in one pass (K1, K3, K4, K7, K8, K9, K10 now contribute real
  numbers instead of being silently unsupported). `units_data.h` and
  `docs/API_CPP.md` / `docs/API_PYTHON.md` updated to match.
- **`MultiPortChamberElement` + `BorderCarnotLossElement`** (Phase 1 of the
  momentum-CV junction). Sanctioned successor to `TeeJunctionElement` for
  n-port manifolds, ejector / merge-split-reverse regimes, and high-Mach
  operation outside the K-closure's structural `K_run ~ 0.30` ceiling
  (see `docs/junction/junction_model_v3_addendum.md` Finding 6 and the
  implementation guide `docs/junction/momentum cv implementation guide.pdf`).
  The junction owns one scalar `P_jct` and emits `N` per-port impulse-function
  residuals plus a global sum-mass residual; loss content moves to companion
  per-port `BorderCarnotLossElement` instances applying the Hager (3/4)
  effective-angle correction (reproduces Hager `xi_l` / Bassett `K_inc`
  exactly at M -> 0 on a sharp-edged lateral). Both elements coexist with
  `TeeJunctionElement` as first-class options; the user selects per junction.
  Solver change: at each port-MCN the mass-conservation row is skipped (the
  junction's sum-mass residual covers conservation; otherwise the channel and
  junction declare equal-and-opposite flows at the port, producing a zero
  row). Tier-4 cross-check against the K-closure and GUI selection dropdown
  follow in subsequent commits.
- **Solver convergence chart**: after every solve a "View convergence" link appears in
  the solver status badge (both success and failure). Clicking opens a modal with a
  log-scale `||F||` vs. time chart (no new dependency — pure SVG). The hybr phase is
  shown in orange; the LM fallback phase (if triggered) in blue with a dashed boundary
  line. A "Worst residuals" table below the chart shows the top-5 unknowns by absolute
  residual value, helping identify which node or element is hardest to converge.
- **Convergence history in API response**: `SolveResponse` now includes
  `convergence_history` (list of `{eval, t, norm}` sampled at every `||F||` improvement
  and at most every 0.1 s, capped at 1 000 points), `worst_residuals` (top 10 per-unknown
  residuals at the final iterate), `solver_settings_used`, and `lm_started_at_eval`.
  Useful for LLM-assisted diagnostics: the response is self-contained with full solver
  context.

### Changed
- **`TeeJunctionElement`**: replaced Bassett 2001 empirical K tables with the
  Unified0D compressible model (Mynard & Valen-Sendstad 2015, extended to O(M²)
  by closed-form Borda-Carnot correction). Residuals now use stagnation pressure
  (p₀) and mass-flow continuity, consistent with compressible duct elements.
  Branching tees use two K-equation residuals (loss from supplier to each collector).
  Merging tees use R_sp (equal static pressure between suppliers) + R_K,com (loss
  from mass-weighted pseudodatum to the collector), matching the LaTeX derivation
  in `docs/junction/junction_model.tex` exactly. Tee port nodes in the GUI backend
  are now `MomentumChamberNode` so static pressure P can differ from Pt at junction
  ports, allowing R_sp to be satisfied at non-zero velocity.
  Branching split closure: both collector legs now use one consistent blended
  turning-loss closure `Pt_com - Pt_leg - K_turn(theta)*q_leg - beta*x_leg*K_bc*q_ref`
  (straight leg: `theta = 0` so the turning term vanishes, leaving the Borda-Carnot
  extraction term). This removes the prior over-charge of a near-closed branch
  (`K -> 1`) and lets a low-Mach `MomentumChamberNode`-bounded branching cascade reach
  a physical all-forward, mass-conserving root. Known envelope limit: the straight run
  is a low-loss through-leg whose effective loss coefficient is capped near
  `beta*4/27 ~ 0.30` of the common dynamic head, so a large imposed straight-run
  stagnation drop (and ejector/reverse behaviour generally) is outside this closure -
  a momentum-CV junction is the sanctioned successor (see
  `docs/junction/junction_model_v3_addendum.md` Finding 6 and
  `docs/junction/momentum cv implementation guide.pdf`). The merging tee is unaffected.

### Fixed
- **Junction validation scorecard**: cells where the candidate model could
  not evaluate any record (e.g. K not yet bound, paper not supported) now
  display `-` instead of `0.000` and a misleading negative `Delta_vs_ceiling`.
  Fully-unsupported cells are also excluded from the headline aggregates.
- **`tee_junction.h` K8 transcription**: K8 (joining type 4, lateral path) had
  `c / psi` where Bassett Table 2 + Eq 34 angle correction call for `psi * c`.
  Invisible at `psi = 1` (the only point covered by the existing
  `K7-K8` complementary-identity test, since both forms reduce to the same
  expression there); diverges at `psi != 1`. Verdict surfaced by the new
  junction-validation runner (`#183`) on Bassett Fig 12b digitized data
  (`psi in {1, 2, 4}` at theta=45 deg). Corrected K8 matches Bassett's own
  published calc line in Fig 12b; the residual gap vs measured (~0.6 at
  psi=4, q=0.5) is the paper's own acknowledged moderate fit, not a code
  defect. Same fix applied to `dK8_dq`. K7-K8 identity test still passes
  unchanged. New regression test pins K8 across psi=1, 2, 4 at theta=45
  deg.
- **Channel-to-MomentumChamberNode face convention**: `ChannelElement.residuals` coupled
  the upstream node's stagnation pressure to the downstream node's *static* pressure
  (`Pt_up - P_down - dP_friction = 0`). Across an inline `MomentumChamberNode` (where
  `Pt = P + q`) this leaked the node dynamic head `q_N` as a free pressure gain - negligible
  at low Mach, fatal near choke. Since the C++ `dP_calc` is a friction-only stagnation-loss
  (the downstream static pressure is not an input), a channel now couples stagnation to
  stagnation (`Pt_up - Pt_down - dP_friction = 0`). For `PlenumNode`/boundary terminations
  (`Pt = P`) this is identical to the previous form, so only inline-MCN networks change.
  This does not address high-dynamic-head *merging* chambers, which need an axial momentum
  balance over the chamber (angle-dependent transverse loss); that closure is tracked in
  #174 and `test_step_4_adding_bypass` is `xfail` until it lands.
- **Tee residual scaling**: `_build_residual_scales` classified `TeeJunctionElement`
  rows with the mass-flow scale (`ref_mdot`, ~0.1 kg/s) although they are stagnation-
  pressure residuals (~10^5 Pa). The ~10^6 mismatch inflated the tee Jacobian rows and
  pushed the scaled condition number to ~5e8, stalling Newton. Tee rows now use the
  pressure scale (`ref_p`), dropping the scaled condition number to ~7e2.
- **Tee compressibility-correction sign**: the O(M^2) correction kappa in
  `K_dat_j_closed` (`include/tee_junction.h`) was missing the bracket factor, giving the
  wrong sign so that compressibility *raised* K. The corrected closed form
  `kappa = (s/2)[-1 - 2x/sqrt(K_inc) + s/(2 K_inc)]` is <= 0 over the physical domain
  (compressibility lowers K), matching `docs/junction/junction_model_v2.tex`.
- **LM fallback for tee networks**: the Levenberg-Marquardt fallback now also fires for
  any network containing a `TeeJunctionElement` (previously only elements flagged as
  compressible). The two supplier stagnation pressures of a tee share a near-singular
  common mode that defeats hybr; LM's regularised step resolves it, so merging tees now
  converge by default.
- **Branching-tee port areas**: a branching `TeeJunctionElement` now propagates its duct
  geometry (`F_C`, `F_branch`) to auto-area `MomentumChamberNode` ports, which cannot see
  an upstream channel during their own topology resolution and otherwise fall back to a
  sentinel area.
- **Tee phi angle convention**: the branching and merging tee models now use
  `phi_j = 0.75 * |theta_j - theta_dat|` (absolute-difference form) instead of
  the LaTeX formula `0.75*(pi - (theta_dat - theta_j))` which assumed the
  "outward-pointing" angle convention. Python passes `theta_com = theta_str = 0`
  (axis-aligned), so the old formula produced phi_str = 135° instead of 0° and
  phi_bra ≈ 168° instead of 67.5° for a 90° branch, making K values 4–8× too
  large and preventing manifold convergence. The Jacobian `dphi/d(theta_dat)` sign
  is also corrected per branch in the merging model.
- **Tee FD Jacobian eps_P scaling**: the finite-difference pressure step in
  `tee_fd_density_jacobians` is now `max(1.0, |P| × 10⁻⁴)` Pa instead of a fixed 1 Pa,
  preventing a numerically negligible step at high absolute pressures (e.g. 300 kPa inlet).
- **MomentumChamberNode Jacobian coupling for tee junctions**: `MomentumChamberNode`
  was silently dropping the `d(residual)/d(m_dot)` coupling when the upstream element
  was a `TeeJunctionElement`, because the old code emitted `{"tee.m_dot": coeff}` — a
  non-existent variable name. The solver now calls `flow_jac_at_node(nid, indices)` once
  per element in `_propagate_states()` to obtain the correct variable names
  (`tee.m_dot_com`, `tee.m_dot_branch`, or `lc.m_dot`) and propagates them via the new
  `_m_dot_jac_names` attribute on each upstream state. `MomentumChamberNode` accumulates
  these into `_upstream_m_dot_jac` and uses them in `residuals()`. This fixes a Newton
  plateau (solver stall at non-zero residual) that occurred in both incompressible and
  compressible mode whenever a `MomentumChamberNode` was connected to a tee junction.
- **Levenberg-Marquardt fallback for hybr oscillation** on compressible networks: when
  `hybr` (Powell's method) fails to converge on a compressible network — typically because
  the Jacobian is ill-conditioned near a choked operating point and Newton steps bounce
  without making progress — the solver now automatically retries from the best iterate
  using Levenberg-Marquardt.  LM's `(J^TJ + lI)` regularisation limits the step even
  when J is nearly singular, providing monotone reduction in `||F||`.  The timeout budget
  is split 65/35 between hybr and the LM fallback so both methods have working time.
  The initial hybr trust-region is also reduced (`factor=10`, down from the scipy default
  of 100) to dampen oscillating overshoots from the start.
- **Solver convergence hardening** for networks with dead-end branches (Wall nodes) and
  inherited-geometry TeeJunctions (`F_C: null`):
  - `_build_x0` now guards the TeeJunction Bernoulli estimate against `F_C=None`,
    eliminating the `TypeError: NoneType * float` crash when using *Continuation* init
    on a network whose tee junctions inherit area from adjacent channels.
  - Initial guess for dead-end channels (connected to a `WallNode`) is now seeded at
    `m_dot=0` instead of the full upstream flow, so Newton starts from the physically
    correct state.
  - `WallNode` initial pressure is now set equal to the upstream node pressure (no
    friction drop at zero flow) instead of the BFS-propagated underestimate.
  - Default `maxfev`/`maxiter` limit is now `max(500, 200*(n+1))` (matching scipy's
    hybr internal default) computed from the actual problem size, replacing the
    hard-coded 500 that caused convergence failures on larger networks or near-choking
    operating points.
  - `ChannelElement.diagnostics` applies the same Re floor (`sqrt(Re²+64²)`) used by
    the C++ friction correlations, preventing `f` from reporting astronomically large
    values (e.g. 5.7 × 10¹⁶) for dead-end channels with near-zero mass flow.

### Added
- **Copy-paste nodes** (Ctrl/Cmd+C / Ctrl/Cmd+V): selected nodes can be duplicated
  on the canvas. Physics and configuration data (including orientation) are copied;
  custom labels and stale solver results are not. Connections are not duplicated —
  the pasted nodes start unconnected. Pasted nodes land 50 px offset from the
  originals and are immediately selected so they can be repositioned.
- **Upstream geometry inheritance**: geometry fields on elements and nodes can now
  be left blank (schema value `null`); the GUI Inspector shows an inherited value
  as a greyed placeholder with a "Reset to inherited" / "Reset to circular" button.
  Affected fields:
  - `OrificeElement.diameter` — inherits from adjacent `ChannelElement`
  - `MomentumChamberNode.Dh` — inherits from upstream `ChannelElement.D`
  - `CombustorNode.Dh` — inherits from upstream `ChannelElement.D`
  - `AreaChangeElement.F0` / `F1` — inherit from upstream/downstream `ChannelElement`
  - `TeeJunctionElement.F_C` — inherits from common- or straight-arm `ChannelElement`
  - `TeeJunctionElement.F_branch` — new field; inherits from branch-arm `ChannelElement`;
    `psi = F_C / F_branch` is now computed automatically from the two areas, replacing
    the manual psi input
  - `ChannelElement.D` — inherits from upstream `AreaChangeElement.F1`,
    `ChannelElement`, or `TeeJunctionElement` (correct arm area used for branch vs. straight)
  - `ChannelElement.Dh` — new optional field; defaults to `D` (circular assumption);
    override for non-circular ducts to use correct hydraulic diameter in friction/Nu
- **TeeJunction Inspector overhaul**: three-arm area layout (C, S read-only, B) with
  per-arm inheritance and a derived ratio card showing `psi`, `D_C → D_B`, and an
  expansion/contraction label — mirrors the `AreaChange` inspector style. The `F_straight
  = F_C` Bassett (2001) model constraint is now clearly annotated.
- **AreaChange Inspector**: non-editable area ratio summary card showing `F₁/F₀`,
  expansion/contraction/straight label, and equivalent circular diameters `D₀ → D₁`.
- **Discrete loss area inheritance**: `PressureLossElement` now inherits flow
  area from an adjacent `ChannelElement` when none is explicitly set, so discrete
  losses can be placed after any element (not only combustors / momentum chambers).
  The Inspector placeholder and hint text update automatically when a channel is
  connected upstream.
- **Topology resolution robustness**: two-pass topology resolution now handles
  mutual dependencies (e.g. channel and tee both with null geometry) without
  ghost fallback values; deferred fallbacks applied only after second pass.
- `WallNode`: closed-end (zero mass-flow) boundary node for modelling dead-end
  branches and symmetry planes of ring manifolds. Drag from the GUI palette or
  instantiate directly via `combaero.network.WallNode("id")`.

## [0.3.1] - 2026-05-21

### Fixed
- `combaero-gui` dependency pin updated from `~=0.2` to `~=0.3`; the old pin
  resolved to `combaero 0.2.6` on PyPI which predates `TeeJunctionElement` and
  `VortexElement`, causing an `ImportError` at startup (#156).
- `publish-gui.yml` CI workflow upgraded from Node 20 to Node 22; `pnpm latest`
  (v11) requires Node ≥22 and crashed with `ERR_UNKNOWN_BUILTIN_MODULE: node:sqlite`
  on Node 20 (#157).

## [0.3.0] - 2026-05-20

### Changed
- **Solver hardening** against unphysical intermediate states during Newton iteration:
  - Temperature floor raised from 50 K to 200 K (both in `_propagate_states` and
    `_get_node_state`); ceiling of 5000 K added to bound combustion runaway.
  - Static pressure floor of 1000 Pa added in `_get_node_state` to prevent NaN
    propagation in thermo calls when P drifts near zero.
  - Penalty residual (triggered on unphysical state exceptions) now encodes
    `F = x_scaled − x_best_scaled` with `J = I`, so the Newton step jumps back
    exactly to the best physical iterate seen so far — previously the constant
    penalty `[10,…,10]` could drive iterates further into unphysical territory.
  - `CombustorNode.compute_derived_state` now wraps the C++ combustion call in a
    try/except and re-raises with a diagnostic message; the solver penalty path
    already handled C++ exceptions, but the new wrapper provides traceability.
  - `graph_builder` now seeds `initial_guess` from the previous solve result
    (`data.result` in node/element data) when no explicit user override exists,
    giving repeated solves a warm start without requiring the `continuation`
    strategy.

### Added
- Global ambient conditions (T, P, RH) in the sidebar "Global Settings" panel, defaulting
  to ISO 2314 gas-turbine reference conditions (288.15 K, 101325 Pa, RH 0.6). When set,
  these override the per-node `ambient_T`, `ambient_P`, and `relative_humidity` fields on
  all boundary nodes that use the `humid_air` composition source. The CompositionEditor
  shows a live `global: —` / `global: value` hint below each affected field.
- All new boundary nodes now default to `humid_air` composition (was `dry_air`). Existing
  saved networks with explicit `source: "dry_air"` are unaffected.
- `test_humid_air_rh0_equals_dry_air`: verifies `species.humid_air_mass(T, P, RH=0)` is
  identical to `species.dry_air_mass()` within floating-point tolerance.
- Global Nu/f multipliers in the sidebar "Global Settings" panel. Setting either
  scales all channel, combustor, momentum-chamber, and discrete-loss elements
  multiplicatively on top of their per-element values — useful for network-wide
  sensitivity studies without touching individual elements.
- `VortexElement` network element: models centrifugal pressure rise in rotating cavities
  (disc pumps, pre-swirl chambers, inter-stage labyrinth spaces) using the Vatistas n-vortex
  model. Parameters: `r_c` (core radius), `r_out`, `r_in` (evaluation radii), `omega_rpm`
  (shaft speed; `None` = inherit global setting), `n` (Vatistas shape parameter, default 2).
  Residual `Pt_out - Pt_in - dP_vortex = 0` with analytical Jacobians w.r.t. Pt and inlet
  density. GUI: Tornado icon (violet), inspectable from the sidebar with a global shaft-speed
  field and per-element override.
- Vatistas (1991) n-vortex model: `vatistas_v0_bar`, `vatistas_vr_bar`,
  `vatistas_pressure_integral` and their analytical derivatives, plus dimensional
  `vatistas_v_theta` / `vatistas_delta_p` with full analytical Jacobians w.r.t.
  r, Γ and r_c. Closed-form antiderivatives for n=1 and n=2; composite Simpson
  quadrature for general n. Python class `VatistasVortex` wraps the free functions
  with vectorised `V_theta`, `dV_theta_dr`, `delta_P`, and `delta_P_bar` methods.
  Validated against digitised experimental data from the original paper
  (Vatistas, Kozel, Mih, *Exp. Fluids* 11, 73–76, 1991).
- Channel element friction model selector in the GUI Inspector: Haaland (default),
  Serghides, Colebrook-White, or Petukhov. The Roughness field is greyed out when
  Petukhov is selected (smooth-pipe model, ignores roughness). Petukhov is the natural
  pairing for channels that already use the Petukhov/Gnielinski heat transfer correlation.
- `NetworkRunner` and `NetworkResult` in `gui/backend/runner.py`: programmatic
  network driver for use in Python scripts and Jupyter notebooks without
  starting the web server.
  - `NetworkRunner.from_file(path)` / `from_dict(d)`: load a GUI-saved JSON.
  - `NetworkRunner.solve(overrides, method, init_strategy, timeout)`: stateless
    single solve with optional boundary-condition overrides keyed by node label
    (``"air_inlet.m_dot": 1.2``).  Each call deep-copies the schema so repeated
    calls inside optimisation loops do not interfere.
  - `NetworkRunner.sweep(params, metrics, ...)`: parametric sweep over a
    DataFrame of BC overrides; returns a compact result DataFrame (one row per
    solve) when ``metrics`` is given, or the full per-entity detail export
    stacked with a ``_sweep_index`` column when omitted.
  - `NetworkResult.get(key)`: scalar result by raw solver key **or by
    ``<label>.<quantity>`` format** (e.g. ``result.get("combustor.T")``),
    resolving GUI node labels to internal IDs automatically.
  - `NetworkResult.node_state(label)`: full solved thermodynamic state dict for
    a node by GUI label.
  - `NetworkResult.to_dataframe()`: full-detail DataFrame identical to the GUI
    CSV export, with unit-annotated column names.
  - `NetworkResult.swap_boundary(label, new_type)`: retype a boundary node
    (mass ↔ pressure) and auto-populate its required field from the solved
    state (``Pt`` when going mass→pressure; ``m_dot`` when going
    pressure→mass).  Returns a new ``NetworkRunner`` pre-seeded with the full
    solved state — including junction-node pressures and element flows — so the
    first re-solve starts from the exact operating point and converges without
    requiring a warm-start strategy.  Enables off-design pressure sensitivity
    studies and matched boundary-condition workflows.
- `python/examples/network_runner_bc_swap.py`: example loading the compressible
  CH4/air combustor, solving the design point at fixed air mass flow, swapping
  the air inlet to a pressure boundary, and sweeping inlet total pressure
  ±10 % around the design point (11/11 points converge).
- `_schema_maps` and `_build_result_objects` helpers extracted from the HTTP
  solve/export handlers into `runner.py`; shared by `NetworkRunner` and
  `main.py` to eliminate duplication.
- `python/examples/network_runner_combustion_sweep.py`: end-to-end example
  loading a GUI-saved compressible CH4/air combustor network, injecting node
  labels, and sweeping equivalence ratio 0.3–0.9 via `NetworkRunner.sweep()`.
- `scripts/check-python-style.sh` now scans `gui/backend/` alongside the core
  Python source directories so backend style is checked locally, not only in CI.
- `scripts/check-gui-style.sh` now supports a `--fix` flag; default (no flag)
  is read-only, matching the CI `lint-gui` job.

### Fixed
- Network solver no longer crashes on combustor networks where Newton steps temporarily
  produce non-physical states (T ≤ 0 or P ≤ 0). `adiabatic_T_complete_and_jacobian_T()`
  and `adiabatic_T_equilibrium_and_jacobians()` now clamp T and P to safe floors instead of
  throwing; `_propagate_states` and `_get_node_state` in the Python solver clamp T ≥ 50 K
  so that reverse-flow mixing artefacts (negative T_mix from negative m_dot streams) never
  propagate into physics calls.
- `friction_petukhov()` and `friction_and_jacobian_petukhov()` no longer throw
  for Re < 3000; the Petukhov formula is mathematically continuous so values below
  the validated range are returned as extrapolated rather than raising an error.
  `CorrelationValidity::EXTRAPOLATED` is set when Re < 3000; `VALID` otherwise.
- `NetworkGraphSchema` now normalises the camelCase ``solverSettings`` key
  (written by the GUI's JSON export) to ``solver_settings`` via a Pydantic
  ``model_validator``, preventing solver regime, method, and timeout settings
  from being silently dropped when loading a GUI-saved file.

### Changed
- `gui/backend/main.py`: result-building and DataFrame construction replaced
  with calls to the shared helpers in `runner.py` (~135 lines removed).

- Thermal wall GUI: arrow direction, probe hot/cold labels, dropdown options, depth hint
  ("0 = Hot, L = Cold"), and "Hot Side Wall Temp" now all reflect the actual heat flow
  direction from the solver (Q sign) rather than the topology of the drawn edge. Wiring a
  thermal wall in either direction gives correct visual feedback without manual adjustment.
- Polar molecule transport model (Monchick-Mason Omega*(2,2) table, 37 T* x 8 delta* points):
  - `Transport_Props` struct gains `dipole_moment` [Debye] and `z_rot` fields.
  - `omega22()` refactored to accept `(T_star, delta_star)` and bilinearly interpolate
    the full Monchick-Mason table (same table as Cantera's `MMCollisionInt`).
  - `compute_delta_star()` helper computes reduced dipole moment in SI.
  - Thermal conductivity upgraded from simple Eucken to Mason-Monchick modified Eucken
    with Parker Z_rot temperature correction (matches Cantera's `fitProperties` formula).
  - All 15 species entries updated with `dipole_moment` and `z_rot` (GRI-Mech 3.0 values).
  - H2O dipole was previously missing; now correctly set to 1.844 D.
  - NH3 polarizability corrected from 0.0 to 2.100 A^3.
  - `generate_thermo_data.py` updated to emit both new fields; `TransportProps` dataclass
    extended with `dipole_moment` field.
- Polar species Cantera validation tests (`TestPolarSpeciesTransport`): pure NH3 and H2O
  viscosity/conductivity vs GRI-Mech at 300/1000/2000 K; NH3/air and H2O/N2 mixtures.

### Fixed
- Compressible solver convergence with mass-flow boundary inlets: `_propagate_mdot_guess`
  was applying a Mach-0.3 cap even to elements directly seeded by a `MassFlowBoundary`,
  causing the initial guess to underestimate the known mass flow and leading hybr to
  converge to a spurious local minimum. The cap is now skipped for elements whose
  immediate upstream node is a `MassFlowBoundary`.
- Edge connection paths for rotated nodes: `FlowEdge` and `ThermalEdge` now apply a
  `rotatePosition` helper that maps each handle's static `position` prop (e.g.
  `Position.Left`) through the connected node's CSS rotation in 90° steps. This gives
  `getBezierPath` the correct perpendicular exit/entry direction after rotation — for
  example, a 180°-rotated node's `flow-target` (logically Left) is correctly treated as
  Right, producing a smooth C-curve instead of an initial vertical leg routing behind the
  element bounding box.
- Rotated-node handle placement and connector triangle orientation: each `Handle` now
  receives an inline style from `handleStyle(basePos, rotation)` in `nodeUtils.ts`.
  This overrides ReactFlow's position CSS classes (which apply in the rotated coordinate
  system and would land a handle on the wrong visual edge) and counter-rotates the handle
  by `−rotation` so the `::after` CSS triangle defined in `index.css` always points
  perpendicular to the edge it sits on.  The centering uses `calc(50% − 6px)` rather than
  `transform:translate`, so `getBoundingClientRect()` returns the correct screen centre
  regardless of node rotation — keeping edge spline endpoints accurate.
- `combaero-gui` white page on fresh PyPI install: `frontend/dist/assets/` was not bundled in
  the wheel because the `**` glob in `package-data` is unreliable below setuptools 69; added an
  explicit `assets/*` pattern and a `MANIFEST.in` as belt-and-suspenders.
- Silent NH3 composition drop in `TestTransportProperties` and `TestCombustionEquilibrium`
  Cantera validation tests: hardcoded 14-item species lists were missing NH3 (index 14), so
  any NH3 mole fraction was silently discarded when building the Cantera composition dict.
- Floating-point precision failure in `test_channel_ribbed_jacobians` on Linux: the
  `dq/dT_hot == -h` identity is exact analytically but accumulated ~4 × 10⁻¹⁴ rounding
  error after NH3 extended `dry_air()` from 14 to 15 elements; replaced exact equality with
  `abs(...) < 1e-9` (matching the pattern used elsewhere in the same file).

### Added
- Ammonia (NH3) species to the thermodynamic and transport database.
- Dynamic species count support in C++ and Python test suites to prevent regressions when modifying the database.
- Cantera validation tests fully migrated to dynamic species API (`num_species()` /
  `species_name(i)` from `combaero._core`); hardcoded species lists removed from all
  `set_cantera_composition` helpers and `SPECIES_ORDER` / `CB_SPECIES` constants.
- `scripts/test-pypi-wheel.sh`: local pre-publish wheel verification — builds both wheels,
  inspects the GUI wheel for bundled assets, installs in an isolated venv, and smoke-tests that
  the server starts and serves JS with the correct MIME type.
- CI `verify` job in `publish-gui.yml` gates the PyPI publish step on the same wheel checks.
- `test_gui_export_dataframe.py`: unit tests for the CSV export pipeline's DataFrame
  construction, column-rename, and `to_csv` logic; covers pandas 3.0 compatibility.

### Changed
- Narrowed `requires-python` to `>=3.12`; dropped `from __future__ import annotations` throughout
- Self-referential `classmethod` return types migrated to `typing.Self`

### Added
- Tee junction pressure-loss model (Bassett 2001) in `include/tee_junction.h`:
  - Pure K-coefficient functions K5, K6, K11, K12 with analytical dK/dq derivatives.
  - `merging_tee_K_straight/branch` and `branching_tee_K_straight/branch`: blended,
    smooth K functions valid for any real inputs (soft-lower protection, tanh blend for
    topology reversal, no NaN or throws outside the validated range).
  - `tee_check_inputs(q, psi, theta)` returns `TeeInputStatus` with per-parameter
    validity flags and an overall `CorrelationStatus` (Valid / Extrapolated).
  - `soft_lower(x, lo)`: smooth max(x, lo) using sqrt regularisation.
- `TeeJunctionResult` struct and `merging_tee_residuals_and_jacobian` /
  `branching_tee_residuals_and_jacobian` in `solver_interface.h/.cpp`: full Jacobians
  w.r.t. mass-flow rates; FD Jacobians w.r.t. P, T, and species mass fractions.
- Python bindings in `_core.cpp` exposing `TeeJunctionResult`, both residual functions,
  and utility functions (`tee_K5`, `tee_K6`, `tee_K11`, `tee_K12`, `tee_blend_weight`,
  `tee_check_inputs`, `merging/branching_tee_K_straight/branch`).
- GTest suite `tests/test_tee_junction.cpp` (13 tests) and Python test suite
  `python/tests/test_tee_junction.py` (19 tests) covering reference values, derivative
  accuracy, zero-flow regularisation, robustness outside the validated range, and
  Jacobian finite-difference verification.
- `TeeJunctionElement` network component: three-port flow element using Bassett 2001
  coefficients, solvable via `NetworkSolver`.  Supports both merging and branching
  tee configurations.  Unknowns are `{id}.m_dot_com` (total common-arm flow) and
  `{id}.m_dot_branch`; straight flow is implicit.  Analytical Jacobian provided for
  all unknowns including pressures, temperature, and species.
- Multi-port element protocol on `NetworkElement`: `all_source_nodes()`,
  `all_sink_nodes()`, `flow_at_node()`, `flow_jac_at_node()` extension points allow
  N-port elements to integrate with the existing solver and graph topology without
  breaking 2-port elements.
- `NetworkSolver` Bernoulli-based initial guess for `TeeJunctionElement` mass flows
  (A * sqrt(2*rho*dP)) to ensure Newton converges from a physically reasonable start.
- Network-level tee tests in `python/tests/test_tee_network.py` (15 tests): graph
  connectivity, solver convergence and mass conservation for both merging and
  branching configurations, Jacobian FD verification, and `validate()` error paths.
- `TeeJunctionElement.diagnostics()` extended to expose C++-computed `K_straight`,
  `K_branch`, and `correlation_extrapolated` (1.0 when Bassett inputs are outside the
  validated range) alongside the existing `m_dot_com/straight/branch` and `q` fields.
- `test_tee_network.py`: conservation tests for mixed-composition merging streams
  (`test_merging_tee_node_mass_conservation`, `test_merging_tee_species_conservation`,
  `test_merging_tee_enthalpy_conservation`) and full central-difference Jacobian coverage
  for both merging and branching tees including mixed-composition sensitivities.
- GUI tee junction inspector, palette item, and React Flow node visual.
- GUI solver status panel now shows full error messages in a scrollable area (previously
  truncated to a single line).
- Tee junction node now shows port labels S / C / B (straight / common / branch) in
  blue (input) or amber (output) based on the current tee type, so wiring direction
  is unambiguous at a glance.
- Inspector shows a port guide card that updates with the selected tee type and
  describes which arms are inlets vs. outlets.
- `graph_builder` auto-inserts a `PlenumNode` junction when an element is connected
  directly to a tee port, removing the requirement to manually place a plenum on every
  tee arm.

### Changed
- `TeeJunctionElement.diagnostics()` key `q` renamed to `mass_flow_ratio`
  (unit kg/kg) for clarity; updated in Python, GUI inspector, CSV export, and
  `quantities.ts` catalogue.
- `TeeJunctionElement.validate()` now accepts negative branch angles: `|theta|`
  must be in `(0, pi/2]` and `abs(theta)` is passed to the C++ correlation,
  since the model uses `cos(0.75*theta)` which is symmetric. The inspector hint
  text and frontend validation both reflect the `(-90, 0) U (0, 90]` deg domain.
- GUI FLOW DISTRIBUTION section in the tee inspector now matches the LIVE
  TELEMETRY visual style (card background, `text-xs` heading, stacked label+value
  cells with `gap-y-3`).
- `TeeJunctionData` schema gains `initial_guess: dict[str, float]` field;
  `graph_builder` wires it through `_expand_initial_guess`, and the inspector
  shows an `InitialGuessEditor` for `m_dot_com` and `m_dot_branch`.
- Frontend validation (`validation.ts`) covers tee junction `theta`, `F_C`, and
  `psi` with the same yellow warning-box UX as other elements.
- CSV export mole fractions X for tee junction common-arm state (consistent with
  other element types).

### Fixed
- `NetworkSolver` Bernoulli initial guess for `TeeJunctionElement` now uses the
  propagated `p_guess` dict for non-`PressureBoundary` nodes (plenums, mass-flow
  boundaries, momentum chambers).  Previously `getattr(node, "Pt", ref["P"])` fell back
  to the reference pressure for interior nodes whose `Pt` is a solver unknown, giving
  `dP_total = 0` and a degenerate `m_dot` starting point.
- CSV export for tee junction rows now includes the common-arm thermodynamic state
  (T, P, Pt, Tt, rho, and full species composition Y[N2] ... Y[NH3]) so that the
  canonical mixed state is accessible without joining to node rows.
- `_propagate_mdot_guess` overwrote rather than accumulated flow for multi-source elements
  (merging tee): the branch-arm visit (0.1 kg/s) silently clobbered the straight-arm visit
  (1.0 kg/s), leaving the tee initial-guess at 0.1 instead of 1.1 kg/s.  The resulting 7x
  underestimate prevented Newton from converging for small `F_C` (e.g. 0.01) where the
  quadratic pressure-drop scaling amplifies initial-guess errors.

## [0.2.6] - 2026-05-03

### Fixed
- `combaero-gui` wheel did not bundle `frontend/dist/assets/`; the `**` glob in
  `package-data` is unreliable below setuptools 69 — added explicit `assets/*`
  pattern and `MANIFEST.in` as belt-and-suspenders (#124).
- `setuptools_scm` appended a local version identifier (`+g…`) to the GUI wheel
  version, causing PyPI to reject uploads; suppressed with
  `local_scheme = "no-local-version"` (#125).
- Tracked `gui/frontend/dist/` in git caused dirty working tree during CI builds,
  making the version suffix appear; removed from index (#125).
- `publish-gui.yml` skipped after re-pushing a tag because the `workflow_run`
  trigger resolved to `failure`; replaced with a direct `push: tags: v*` trigger
  (#126).

### Added
- `scripts/test-pypi-wheel.sh`: local pre-publish wheel verification — builds both
  wheels, inspects bundled assets, installs in an isolated venv, and smoke-tests
  MIME types (#124).
- CI `verify` job in `publish-gui.yml` gates the PyPI publish step on the same
  checks (#124).

## [0.2.4] - 2026-05-02

### Fixed
- GUI bundle not rebuilt before release; SPA catch-all route incorrectly served
  JSON for unmatched paths; development banner visible in production build (#118).

## [0.2.3] - 2026-05-02

### Fixed
- Relative API URLs broke on non-root deployments (#117).
- SPA catch-all route returning JSON instead of HTML for deep-linked paths (#117).
- PyPI documentation links pointed to wrong locations (#117).

## [0.2.2] - 2026-05-02

### Fixed
- `combaero-gui` frontend assets not served when installed from PyPI; broken
  documentation code examples (#116).

## [0.2.1] - 2026-05-02

### Fixed
- `gui` extra missing from `combaero-gui` package metadata (#115).
- `gui/frontend/dist/` artefacts incorrectly tracked in git (#115).

## [0.2.0] - 2026-04-26

### Added
- Network GUI: React Flow visual network builder backed by FastAPI (`combaero-gui`)
- Analytical combustor Jacobians and full network Jacobian integration
- Continuation solver with GUI method-selection sidebar
- Stagnation and transport property helpers with C++ performance optimisations
- Inverse solver API: `calc_T_from_u`, `dg_over_RT_dT` bindings
- Multi-stream axial staging combustor network example
- Diagnostic telemetry schema (Pt/Tt stagnation properties)
- Solver hard timeouts and empty-network guard
- Area-change element diagnostics and validation hardening

### Changed
- uv workspace consolidation; project-wide ruff/uv modernisation
- pybind11 requirement raised to `>=3.0.4`
- Frontend toolchain standardised on pnpm, Vitest, Node 24 LTS
- CI: separated `gui` path filter from C++ source matrix; added `lint-gui` job

### Fixed
- State property units mismatch in example runner
- Solver stability under degenerate initial conditions

[Unreleased]: https://github.com/thiemom/combaero/compare/v0.2.6...HEAD
[0.2.6]: https://github.com/thiemom/combaero/compare/v0.2.4...v0.2.6
[0.2.4]: https://github.com/thiemom/combaero/compare/v0.2.3...v0.2.4
[0.2.3]: https://github.com/thiemom/combaero/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/thiemom/combaero/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/thiemom/combaero/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/thiemom/combaero/releases/tag/v0.2.0
