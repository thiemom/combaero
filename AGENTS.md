# CombAero — Agent Instructions

This file provides instructions for AI coding agents (OpenAI Codex, Cursor, etc.) working on the CombAero project. Follow these rules exactly.

## 1. Hard Rules

- **Virtual environment only.** Use `uv run` and `uv pip`. Never install into system Python or global `pip`.
- **Never skip CI.** Do not bypass pre-commit hooks (`--no-verify`) or GitHub Actions checks. All code must pass cleanly.
- **No manual edits to auto-generated files.** `docs/UNITS.md` and `include/thermo_transport_data.h` are generated — edit their sources and re-run the generator scripts.
- **API sync is mandatory.** Adding, removing, or changing any function, property, or unit requires updating `include/units_data.h` and the relevant API reference (`docs/API_CPP.md` or `docs/API_PYTHON.md`) in the same commit.
- **CHANGELOG is mandatory.** Every user-visible change (feature, bug fix, breaking change) must get an entry under `[Unreleased]` in `CHANGELOG.md` in the same commit. Use sub-headings Added / Changed / Fixed / Removed. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). At release time the maintainer promotes `[Unreleased]` to `[x.y.z] - YYYY-MM-DD`.
- **Documentation hygiene.** Follow `docs/HYGIENE.md`. Archive feature reports to `docs/archive/`; keep core guides high-signal.
- **Unit sync test.** Every exported Python symbol needs a `units_data.h` entry, or an explicit `IGNORE_LIST` entry in `python/tests/test_units_sync.py`.
- **Solver (f, J) rule.** Solver-facing calculations must expose a C++ PyBind11 API returning `std::tuple<double, double>` (value, derivative). Use analytical chain-rule derivatives; no finite differences exposed to Python.
- **Define Once.** Physics constants go as `constexpr` in the relevant public header namespace — never as magic numbers in `.cpp` files.
- **ASCII only.** No non-ASCII characters in any C++ or Python source file.

## 2. C++ Style (C++17)

- `#pragma once`, sorted minimal includes, `//` line comments only (no `/* */`).
- Smart pointers and RAII. No raw `new` / `delete`.
- Use `math_constants.h` for `M_PI` (MSVC compatibility).
- Explicit standard library includes — macOS implicit includes break Linux CI.

## 3. Python Style

- `ruff` for lint + format. Run `./scripts/check-python-style.sh --fix` before committing.
- Type annotations on all signatures; `-> None` for void functions; `|` union syntax (Python 3.12+).
- No non-ASCII characters, including docstrings and comments.

## 4. Frontend Style (`gui/frontend`)

- Package manager: `pnpm`. Use `pnpm install`, `pnpm run dev`, etc.
- Linting and formatting: `biome`. Run `./scripts/check-gui-style.sh` before committing.
- Biome config: `biome.json` at repo root.

## 5. Key Scripts

```bash
./scripts/check-source-style.sh             # C++ style
./scripts/check-python-style.sh --fix       # Python lint + format
./scripts/check-gui-style.sh                # Frontend lint + format (biome)
uv pip install -e .                         # rebuild C++ extension
uv run pytest python/tests/                 # Python tests
cd build && ctest                           # C++ tests
uv run scripts/generate_units_md.py         # regenerate docs/UNITS.md
./scripts/test-pypi-wheel.sh                # verify combaero-gui wheel before PyPI publish
```

## 6. GUI Wheel Verification

Before any `combaero-gui` PyPI release, run:

```bash
./scripts/test-pypi-wheel.sh
# or, if gui/frontend/dist/ is already fresh:
./scripts/test-pypi-wheel.sh --skip-frontend
```

This script:
1. Builds the React frontend (`pnpm build`)
2. Builds the Python wheel (`uv build --package combaero-gui`)
3. Inspects the wheel zip to confirm `frontend/dist/assets/*.js` is bundled — if this directory is missing, FastAPI skips the `/assets` StaticFiles mount at startup and every JS request gets `text/html`, causing a blank page
4. Installs the wheel in a clean isolated venv (pulling `combaero` from PyPI, exactly as an end-user would) and smoke-tests that `/health` responds and JS assets are served with the correct MIME type

The CI `publish-gui.yml` workflow runs an equivalent `verify` job and gates the publish step on it. Run the script locally before tagging a release.

## 7. Auto-Generated Files

| Target | Source | Regenerate with |
|--------|--------|-----------------|
| `docs/UNITS.md` | `include/units_data.h` | `uv run scripts/generate_units_md.py` |
| `include/thermo_transport_data.h` | `thermo_data_generator/` | see `thermo_data_generator/README.md` |

## 8. Version Sync Checklist

When changing Python version or tool versions, keep in sync:
`.python-version` · `ruff.toml` · `ci.yml` · `validation-tests.yml` · `.pre-commit-config.yaml` · `Makefile`

When changing frontend tool versions, keep in sync:
`gui/frontend/package.json` · `biome.json` (`$schema` version) · `.pre-commit-config.yaml` (biome rev + `@biomejs/biome` pin)
