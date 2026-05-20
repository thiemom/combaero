# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
