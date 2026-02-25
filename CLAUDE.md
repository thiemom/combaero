# CombAero — Claude Code Instructions

## Hard Rules
- **Virtual environment only:** use `.venv/bin/python` and `.venv/bin/pip`. Never install into system Python.
- **Never bypass pre-commit hooks** or CI checks. All code must pass cleanly before pushing.
- **No manual edits to auto-generated files** (`docs/UNITS.md`, `include/thermo_transport_data.h`). Use the scripts below.
- **API sync is mandatory:** adding/removing/changing any function or property requires updating `include/units_data.h` and `docs/API_REFERENCE.md` in the same commit.
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

## Key Scripts
```bash
./scripts/check-source-style.sh          # C++ style
./scripts/check-python-style.sh --fix    # Python lint + format
.venv/bin/pip install -e . --no-build-isolation  # rebuild C++ extension
.venv/bin/pytest python/tests/           # Python tests
cd build && ctest                        # C++ tests
python scripts/generate_units_md.py      # regenerate docs/UNITS.md
```

## Auto-Generated Files
| Target | Source | Regenerate with |
|--------|--------|-----------------|
| `docs/UNITS.md` | `include/units_data.h` | `python scripts/generate_units_md.py` |
| `include/thermo_transport_data.h` | `thermo_data_generator/` | see `thermo_data_generator/README.md` |

## Version Sync Checklist
When changing Python version or tool versions, keep in sync:
`.python-version` · `ruff.toml` · `ci.yml` · `validation-tests.yml` · `.pre-commit-config.yaml` · `Makefile`
