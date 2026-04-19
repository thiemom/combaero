# 🧹 Unified `uv` Workspace & `ruff` Consolidation

This PR completes the migration to a unified `uv` workspace and strictly enforces `ruff` as the project-wide linting and formatting tool, removing all legacy Poetry, Black, and Pip fragments.

## Key Changes

### 🔧 Workspace Completion
- **`thermo_data_generator` Migration**: Converted the thermo data generator from Poetry to a PEP 621 `pyproject.toml`. It is now a formal member of the root `uv` workspace.
- **Root Workspace Integration**: Updated the root `pyproject.toml` to include all sub-projects, enabling single-command environment synchronization.

### 🛡️ Tooling Consolidation
- **Makefile Overhaul**: Refactored the root `Makefile` to strictly use `uv run` and `uv sync`.
- **Legacy Removal**: Deleted all surviving references to `black`, `isort`, and `flake8` in Makefiles and requirements.
- **Pre-commit Modernization**: Updated `.pre-commit-config.yaml` hooks to execute via `uv run`, ensuring consistent dependency resolution during git hooks.

### 📚 Documentation & Help
- **Helper Scripts**: Updated `cantera_validation_tests/run_tests.sh` and internal Makefiles to use `uv`.
- **README Cleanups**: Replaced all `poetry` command examples with their `uv` equivalents across all sub-project documentation.

## ✅ Verification
- **Full Workspace Sync**: `uv sync` resolved and installed all dependencies for the core, validation suite, and thermo generator.
- **Unified Linting**: `make style-check` (now running `ruff` via `uv`) passes across the entire codebase.
- **Test Execution**: Verified that `cantera_validation_tests` and python units run correctly through the new `uv` shortcuts.
