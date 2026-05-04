# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- `combaero-gui` white page on fresh PyPI install: `frontend/dist/assets/` was not bundled in
  the wheel because the `**` glob in `package-data` is unreliable below setuptools 69; added an
  explicit `assets/*` pattern and a `MANIFEST.in` as belt-and-suspenders.

### Added
- `scripts/test-pypi-wheel.sh`: local pre-publish wheel verification — builds both wheels,
  inspects the GUI wheel for bundled assets, installs in an isolated venv, and smoke-tests that
  the server starts and serves JS with the correct MIME type.
- CI `verify` job in `publish-gui.yml` gates the PyPI publish step on the same wheel checks.
- `test_gui_export_dataframe.py`: unit tests for the CSV export pipeline's DataFrame
  construction, column-rename, and `to_csv` logic; covers pandas 3.0 compatibility.

### Changed
- Narrowed `requires-python` to `>=3.12`; dropped `from __future__ import annotations` throughout
- Self-referential `classmethod` return types migrated to `typing.Self`

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

[Unreleased]: https://github.com/thiemom/combaero/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/thiemom/combaero/releases/tag/v0.2.0
