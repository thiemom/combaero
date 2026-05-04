# CombAero Scripts

Utility scripts for code quality, documentation generation, and testing.

## Code Quality Scripts

### check-source-style.sh

Check C++ source code style and quality.

```bash
./scripts/check-source-style.sh
./scripts/check-source-style.sh --verbose
./scripts/check-source-style.sh --strict
```

**Checks:** Non-ASCII characters, block comments, M_PI without `math_constants.h`.
**Scans:** `src/`, `include/`, `examples/`, `tests/`

### check-python-style.sh

Check Python code style and quality.

```bash
./scripts/check-python-style.sh        # Check
./scripts/check-python-style.sh --fix  # Auto-fix formatting
```

**Checks:** Non-ASCII characters, ruff lint, ruff format, optional mypy.
**Scans:** `python/`, `cantera_validation_tests/`, `thermo_data_generator/`

### check-gui-style.sh

Check frontend code style (Biome).

```bash
./scripts/check-gui-style.sh
```

### run-clang-tidy.sh

Run clang-tidy static analysis (requires a CMake build with `compile_commands.json`).

```bash
./scripts/run-clang-tidy.sh
```

## Documentation Scripts

### generate_units_md.py

Regenerate `docs/UNITS.md` from `include/units_data.h`.

```bash
uv run scripts/generate_units_md.py
```

## Testing Scripts

### test-pypi-wheel.sh

Build and verify the `combaero-gui` wheel end-to-end before a PyPI release. Builds
the frontend and both wheels, inspects wheel contents, installs in an isolated venv,
and smoke-tests that the server starts and serves JS assets with the correct MIME type.

```bash
./scripts/test-pypi-wheel.sh                        # Full run
./scripts/test-pypi-wheel.sh --skip-frontend        # Reuse existing gui/frontend/dist/
./scripts/test-pypi-wheel.sh --skip-combaero-build  # Use existing dist/combaero-*.whl
```

## Quick Reference

| Task | Command |
|------|---------|
| Check C++ style | `./scripts/check-source-style.sh` |
| Check Python style | `./scripts/check-python-style.sh` |
| Auto-fix Python formatting | `./scripts/check-python-style.sh --fix` |
| Check GUI style | `./scripts/check-gui-style.sh` |
| Run static analysis | `./scripts/run-clang-tidy.sh` |
| Regenerate units docs | `uv run scripts/generate_units_md.py` |
| Verify GUI wheel before publish | `./scripts/test-pypi-wheel.sh` |

All scripts return `0` on success and `1` on failure.
