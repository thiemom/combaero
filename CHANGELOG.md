# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
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
