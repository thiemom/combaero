# CombAero — Claude Code Instructions

## Hard Rules
- **Virtual environment only:** use `uv run` and `uv pip`. Never install into system Python or use global `pip`.
- **Never bypass pre-commit hooks** or CI checks. All code must pass cleanly before pushing.
- **No manual edits to auto-generated files** (`docs/UNITS.md`, `include/thermo_transport_data.h`). Use the scripts below.
- **API sync is mandatory:** adding/removing/changing any function or property requires updating `include/units_data.h` and the relevant API reference ([docs/API_CPP.md](docs/API_CPP.md) or [docs/API_PYTHON.md](docs/API_PYTHON.md)) in the same commit.
- **CHANGELOG is mandatory:** every user-visible change must have an entry added under `[Unreleased]` in `CHANGELOG.md` in the same commit. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) (Added / Changed / Fixed / Removed). At release time the `[Unreleased]` block becomes `[x.y.z] - YYYY-MM-DD`.
- **Documentation Hygiene:** Follow [docs/HYGIENE.md](docs/HYGIENE.md) classification. Archive feature reports to `docs/archive/`.
- **Unit sync test:** every exported Python symbol needs a `units_data.h` entry, or an `IGNORE_LIST` entry in `python/tests/test_units_sync.py`.
- **Solver (f, J) rule:** solver-facing calculations must expose a C++ PyBind11 API returning `std::tuple<double, double>` (value, derivative). Analytical derivatives via chain rule; no finite differences exposed to Python.
- **Define Once:** physics constants go as `constexpr` in the relevant public header namespace — never as magic numbers in `.cpp` files.
- **ASCII only:** no non-ASCII characters in any C++ or Python source file.

## C++ Style (C++17)
- `#pragma once`, sorted minimal includes, `//` comments only (no `/* */`).
- Smart pointers, RAII. No raw `new`/`delete`.
- Use `math_constants.h` for `M_PI` (MSVC compatibility).
- Explicit standard library includes — macOS-only implicit includes break Linux CI.

## Python Style
- `ruff` for lint + format. Run `./scripts/check-python-style.sh --fix` before committing.
- Type annotations on all signatures; `-> None` for void functions; `|` union syntax (Python 3.12).
- No non-ASCII characters, including docstrings and comments.

## Frontend Style (gui/frontend)
- Package manager: `pnpm`. Use `pnpm install`, `pnpm run dev`, etc.
- Linting and formatting: `biome`. Run `./scripts/check-gui-style.sh` before committing.
- Biome config: `biome.json` at repo root.

## Key Scripts
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

## Auto-Generated Files
| Target | Source | Regenerate with |
|--------|--------|-----------------|
| `docs/UNITS.md` | `include/units_data.h` | `uv run scripts/generate_units_md.py` |
| `include/thermo_transport_data.h` | `thermo_data_generator/` | see `thermo_data_generator/README.md` |

## Version Sync Checklist
When changing Python version or tool versions, keep in sync:
`.python-version` · `ruff.toml` · `ci.yml` · `validation-tests.yml` · `.pre-commit-config.yaml` · `Makefile`

When changing frontend tool versions, keep in sync:
`gui/frontend/package.json` · `biome.json` (`$schema` version) · `.pre-commit-config.yaml` (biome rev + `@biomejs/biome` pin)
