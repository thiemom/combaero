# CombAero Architecture: Package Separation Strategy

This document details the architecture and rationale for the separation between the core `combaero` physics/solver library and the `combaero-gui` application.

## 1. Rationale for Separation

The core C++ engine, Python bindings (`combaero`), the network solver (`combaero.network`), and the frontend GUI (`gui/`) are decoupled to ensure optimal performance, maintenance, and user experience.

### 1.1 User Perspective
- **The API User (Data Scientist / Engineer)**: Can run `pip install combaero` for a lightweight, fast installation without web dependencies or Node.js assets.
- **The GUI User (Analyst / Designer)**: Receives a standalone application (e.g., `.exe` or `.dmg`) focused on the visual network builder experience.

### 1.2 Maintenance & CI/CD
- **Build Times**: Decoupling prevents UI changes from triggering heavy C++ cross-platform compilation.
- **Release Cadence**: The GUI can iterate rapidly on UX without forcing version bumps on the core mathematical library.
- **Issue Tracking**: Separation of concerns isolates UI state bugs from thermodynamic solver logic.

### 1.3 Coding & Architecture
- **Strict API Boundaries**: Forces a formal, stable contract between the backend (FastAPI) and the frontend.
- **Ecosystem Integrity**: Avoids fragile tool clashing between Python (`pyproject.toml`) and JavaScript (`npm`) build systems.

---

## 2. Monorepo Architecture

The project utilizes a monorepo approach (`thiemom/combaero`) to ease cross-cutting feature development while maintaining distinct distribution packages.

- **`combaero` (Core)**: C++ `src/`, `include/`, and `python/combaero`. Published as PyPI wheels.
- **`combaero-gui` (App)**: `gui/` (React frontend) and FastAPI wrapper. Distributed as a standalone desktop app and optional PyPI package.

---

## 3. Implementation Details

## Current State

### Already done
- Core `pyproject.toml` dependencies: only `numpy>=1.24`, `scipy>=1.10` — no GUI deps
- `gui/pyproject.toml` exists as a separate workspace package (`combaero-gui`) with fastapi/uvicorn/pydantic/pandas isolated there
- `wheel.packages = ["python/combaero"]` — core wheel already excludes `gui/`
- Backend code lives in `gui/backend/` (clean separation of files)

### Gaps to fix
1. **sdist includes `gui/`** — scikit-build-core bundles all git-tracked files into the source distribution by default; `gui/` is tracked but not excluded
2. **No `combaero-gui` CLI entry point** — no `[project.scripts]` in `gui/pyproject.toml`, no launcher function in the backend
3. **CI: GUI changes trigger C++ builds** — `gui/**` is in the `source` path filter in `ci.yml`, causing expensive 3-platform C++ matrix builds on every UI change

---

## Implementation

### 1. Exclude `gui/` from core sdist — `pyproject.toml`

Add under `[tool.scikit-build]`:

```toml
[tool.scikit-build]
wheel.packages = ["python/combaero"]
sdist.exclude = ["gui/*", "gui/**/*"]
build.verbose = true
```

The wheel is already clean (restricted by `wheel.packages`). This closes the sdist gap.

### 2. Add `combaero-gui` entry point — `gui/pyproject.toml`

Add a scripts section:

```toml
[project.scripts]
combaero-gui = "gui.backend.main:run_server"
```

### 3. Add `run_server()` launcher — `gui/backend/main.py`

Add at the bottom of the file (after the FastAPI `app` definition):

```python
def run_server(host: str = "127.0.0.1", port: int = 8000) -> None:
    import uvicorn
    uvicorn.run(app, host=host, port=port)
```

This gives users `combaero-gui` as a shell command after `pip install combaero-gui`.

### 4. Split GUI out of C++ CI path — `.github/workflows/ci.yml`

In the `changes` job filters, move `gui/**` from `source` to a new `gui` filter:

```yaml
filters: |
  source:
    - 'src/**'
    - 'include/**'
    - 'python/**'
    - 'tests/**'
    - 'scripts/**'
    - 'examples/**'
    - '.github/workflows/**'
    - 'CMakeLists.txt'
    - 'pyproject.toml'
    - '.pre-commit-config.yaml'
    - 'ruff.toml'
    - 'Makefile'
  gui:
    - 'gui/**'
  documentation:
    - 'docs/**'
    - 'README.md'
    - '**/*.md'
```

Then expose `gui` as an output from the `changes` job and add a lightweight `lint-gui` job that runs only on `gui` changes (biome + Python style check on `gui/`), without triggering `build-and-test`.

---

## Critical Files

| File | Change |
|------|--------|
| `/Users/thiemo/Projects/combaero/pyproject.toml` | Add `sdist.exclude` under `[tool.scikit-build]` |
| `/Users/thiemo/Projects/combaero/gui/pyproject.toml` | Add `[project.scripts]` with `combaero-gui` entry |
| `/Users/thiemo/Projects/combaero/gui/backend/main.py` | Add `run_server()` function |
| `/Users/thiemo/Projects/combaero/.github/workflows/ci.yml` | Split `gui/**` into separate filter; add `lint-gui` job |

---

## Verification

1. **sdist exclusion**: `python -m build --sdist` on the root package, inspect the `.tar.gz` — `gui/` must not appear
2. **Entry point**: after `uv pip install -e gui/`, run `combaero-gui` — uvicorn should start on `127.0.0.1:8000`
3. **CI**: open a PR that only changes a file in `gui/frontend/` — verify only the `lint-gui` job runs, not `build-and-test`
