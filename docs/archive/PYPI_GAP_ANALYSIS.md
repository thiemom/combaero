# PyPI Publishing Gap Analysis

**Date**: 2026-05-01
**Scope**: PyPI readiness for `combaero` and `combaero-gui`

## 1. Review of Proposed Structure (`docs/PACKAGING.md`)

The strategy outlined in `docs/PACKAGING.md` proposes a monorepo approach with a strict boundary between the core C++ physics library (`combaero`) and the React/FastAPI desktop application (`combaero-gui`).

**Assessment**: This is a robust and highly recommended architecture. It prevents users who only want the Python math library from downloading heavy Node.js assets or web frameworks, while keeping development friction low in a single repository.

**Implementation Status**: The steps proposed in `docs/PACKAGING.md` are **already completed**:
- ✅ `gui/` is excluded from the core sdist via `[tool.scikit-build.sdist.exclude]` in the root `pyproject.toml`.
- ✅ The `combaero-gui` shell entry point exists in `gui/pyproject.toml` (`[project.scripts]`).
- ✅ The `run_server` launcher is implemented in `gui/backend/main.py`.
- ✅ The CI pipeline (`.github/workflows/ci.yml`) cleanly separates `source` and `gui` paths, preventing UI tweaks from triggering heavy C++ cross-compilation matrices.

## 2. Gap Analysis: Open Items for PyPI Publication

While the internal separation is complete, several gaps remain before the package can be successfully and professionally published to the Python Package Index (PyPI).

### Critical / High Priority
1. ~~**Cross-Platform Wheel Builds (`cibuildwheel`)**~~ — **Resolved**: Added `[tool.cibuildwheel]` to `pyproject.toml` (cp312, manylinux_2_28/x86_64, macOS x86_64+arm64, Windows AMD64) and created `.github/workflows/publish.yml` triggering on `v*.*.*` tags via OIDC trusted publisher. `cmake>=3.21` and `ninja` added to `build-system.requires` so cibuildwheel environments get a self-contained toolchain.
2. ~~**License Metadata**~~ — **Resolved**: Added `license-files = ["LICENSE"]` to `[project]` in `pyproject.toml`.
3. ~~**Source Distribution (sdist) Bloat**~~ — **Resolved**: `sdist.exclude` in `pyproject.toml` now also excludes `cantera_validation_tests/`, `thermo_data_generator/`, and `.github/`.

### Medium Priority (Developer Experience)
4. ~~**`combaero-gui` Missing PyPI Metadata**~~ — **Resolved**: `gui/pyproject.toml` now includes `readme`, `license`, `authors`, `keywords`, `classifiers`, and `[project.urls]`, mirroring the core package. `gui/README.md` (pre-existing) is used as the landing-page description. Note: `license = {text = "MIT"}` (inline text) is used rather than `license-files` because the `LICENSE` file is at the repo root, not inside `gui/`; PyPI will display the license type but will not attach the full license text to the `combaero-gui` wheel. Copying `LICENSE` to `gui/LICENSE` and switching to `license-files = ["LICENSE"]` would close this minor gap if needed.
5. ~~**`combaero-gui` Version Not Dynamic**~~ — **Resolved**: `gui/pyproject.toml` switched to `dynamic = ["version"]` with `[tool.setuptools_scm] root = ".."` so it reads the same git tag as the core package.
6. ~~**Python Version Classifier Mismatch**~~ — **Resolved**: `requires-python` narrowed to `">=3.12"` in both `pyproject.toml` and `gui/pyproject.toml`, matching the single `Programming Language :: Python :: 3.12` classifier and the CI-tested version. Six files in `python/combaero/network/` use `X | Y` union syntax without `from __future__ import annotations` (requires Python 3.10+), and the codebase style mandates `|` union syntax throughout, so 3.12 is the correct floor. `from __future__ import annotations` removed from all source files; two self-referential `classmethod` return types migrated to `typing.Self`.
7. ~~**Type Stubs & `py.typed`**~~ — **Resolved**: Empty `python/combaero/py.typed` marker added; included in the wheel via `wheel.packages = ["python/combaero"]`. Generating full `.pyi` stubs with `pybind11-stubgen` remains a future improvement.
8. ~~**Changelog**~~ — **Resolved**: `CHANGELOG.md` created in Keep a Changelog format with the `v0.2.0` entry and an `[Unreleased]` section. Linked from `[project.urls]` as `Changelog`.
9. **Trusted Publisher Setup**:
   - **PyPI projects created**: `combaero` and `combaero-gui` projects exist on PyPI.
   - **Remaining**: Configure Trusted Publishers in each project's PyPI dashboard (`Manage → Publishing`):

   | PyPI project | Repository | Workflow file | Environment |
   |---|---|---|---|
   | `combaero` | `thiemom/combaero` | `publish.yml` | `release` |
   | `combaero-gui` | `thiemom/combaero` | `publish-gui.yml` | `release` |

   The GitHub Actions `release` environment must also be created in the repo settings (`Settings → Environments`).

## 3. Version Number Assessment

- **Current Version**: `v0.2.0` tag exists; HEAD is 11 commits ahead (`v0.2.0-11-g0727427`), so the installed dev version currently resolves to `0.2.0.dev11+g0727427` rather than a clean `0.2.0`. A clean release requires tagging HEAD.
- **Project Status**: The project is highly advanced, featuring a robust thermodynamic engine, a generalized network solver with Jacobian analytics, and a full frontend GUI.
- **Assessment**: The `0.2.0` version is **appropriate** for an initial public release.
  - `0.x` correctly signals to users that the API (especially C++ bindings and network JSON schemas) may still undergo breaking changes before a stable `1.0.0`.
  - The `pyproject.toml` classifier is set to `Development Status :: 4 - Beta`, which perfectly aligns with a `v0.2.x` release candidate.

## Recommendation

Items 1–8 are resolved. The remaining pre-release steps are:
1. Complete **Trusted Publishers** (item 9): create the `release` environment in GitHub repo settings, then add the two Trusted Publisher entries on PyPI as listed in the table above.
2. Tag HEAD with a new version to produce a clean release version string and trigger both publish workflows.
3. Smoke-test the core wheel locally before tagging: `uv build && pip install dist/combaero-*.whl && python -c "import combaero"`.
