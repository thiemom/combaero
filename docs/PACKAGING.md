# CombAero Packaging & GUI Separation Strategy

This document details the rationale and implementation plan for separating the core `combaero` physics/solver library from the `combaero-gui` application.

## 1. Rationale for Separation

Currently, the core C++ engine, Python bindings (`combaero`), the network solver (`combaero.network`), and the frontend GUI (`gui/`) live together and are conceptually packaged as a single unit. As the project matures toward PyPI publishing and end-user distribution, they must be decoupled.

### 1.1 User Perspective
- **The API User (Data Scientist / Engineer)**: Wants to run `pip install combaero`. They use the library in Jupyter notebooks, optimization loops, or larger simulation frameworks. They want a lightweight, fast installation. They **do not** want Node.js assets, FastAPI servers, React bundles, or heavy web dependencies cluttering their environment.
- **The GUI User (Analyst / Designer)**: Wants a standalone application (e.g., double-clicking a `.exe` or `.dmg`). They might not know or care about Python environments, virtual environments, or `pip`. They just want the visual network builder.

### 1.2 Maintenance & CI/CD
- **Build Times**: The core package involves compiling C++ across multiple platforms (macOS, Windows, Linux) which is inherently slow. The GUI builds React/TypeScript via Node.js. Decoupling means changing a UI button styling doesn't trigger the heavy C++ CI pipeline.
- **Release Cadence**: The GUI will likely iterate rapidly (UX improvements, new visualizations). The core solver will iterate slower, requiring strict semantic versioning to protect mathematical correctness and API contracts. Tying them together forces unnecessary version bumps for the core library just to release a UI fix.
- **Issue Tracking**: Separation of concerns allows developers to isolate UI/state bugs from core thermodynamic/solver bugs more easily.

### 1.3 Coding & Architecture
- **Strict API Boundaries**: Physical separation forces a strict boundary (API contract) between the backend (FastAPI) and the frontend. It prevents "leakage" where the GUI relies on internal, undocumented Python states instead of formal, stable endpoints.
- **Ecosystem Clash**: Python packaging (wheels, PyPI) and JavaScript packaging (npm, Electron/Tauri) use vastly different tools. Mixing them in a single `pyproject.toml` build backend (like `scikit-build-core`) becomes a fragile, highly complex setup.

---

## 2. Recommended Architecture

**Recommendation: The Monorepo Approach with Split Packaging**

We should keep both components in the same `thiemom/combaero` Git repository (a monorepo) to ease cross-cutting feature development (e.g., adding a new component requires both C++ solver logic and a GUI node). However, we must **package and distribute** them entirely separately.

- **`combaero` (Core)**: Contains C++ `src/`, `include/`, and `python/combaero` (excluding GUI/server code). Published to PyPI as wheels.
- **`combaero-gui` (App)**: Contains `gui/` (React frontend) and the FastAPI wrapper backend. Distributed primarily as a standalone bundled desktop app via GitHub Releases, and optionally as a separate PyPI package (`pip install combaero-gui`).

---

## 3. Implementation Plan

# Plan: CombAero Package Separation

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
