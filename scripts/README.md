# CombAero Scripts

Utility scripts for code quality, documentation generation, and testing.

## Code Quality Scripts

### check-source-style.sh

Check C++ source code style and quality.

```bash
# Basic check
./scripts/check-source-style.sh

# Verbose output
./scripts/check-source-style.sh --verbose

# Strict mode (fail on warnings)
./scripts/check-source-style.sh --strict
```

**Checks:**
- Non-ASCII characters outside comments and strings
- Block comments (only line comments allowed)
- M_PI usage without math_constants.h

**Scans:** `src/`, `include/`, `examples/`, `tests/`

### check-python-style.sh

Check Python code style and quality (PEP 8 compliance).

```bash
# Basic check
./scripts/check-python-style.sh

# Auto-fix formatting issues
./scripts/check-python-style.sh --fix

# Verbose output
./scripts/check-python-style.sh --verbose

# Strict mode (fail on type hints)
./scripts/check-python-style.sh --strict
```

**Checks:**
1. Non-ASCII characters outside comments and strings
2. PEP 8 compliance (flake8)
3. Code formatting (black)
4. Import sorting (isort)
5. Type hints (mypy, optional)

**Scans:** `python/`, `cantera_validation_tests/`, `thermo_data_generator/`

**Dependencies:**
```bash
pip install black isort flake8 mypy
# or: uv pip install black isort flake8 mypy
```

See [INSTALL.md](INSTALL.md) for detailed installation options (Poetry, uv, pre-commit).

### run-clang-tidy.sh

Run clang-tidy static analysis on C++ code.

```bash
# Requires existing CMake build with compile_commands.json
./scripts/run-clang-tidy.sh
```

**Scans:** `src/`, `include/` (skips tests and Python bindings)

## Documentation Scripts

### generate_units_md.py

Generate `docs/UNITS.md` from `include/units_data.h`.

```bash
python scripts/generate_units_md.py
```

Parses unit annotations and creates formatted documentation organized by module.

## Usage in CI/CD

### Pre-commit Hooks

Add to `.pre-commit-config.yaml`:

```yaml
- repo: local
  hooks:
    - id: check-cpp-style
      name: Check C++ Style
      entry: ./scripts/check-source-style.sh
      language: system
      types: [c++]
      pass_filenames: false

    - id: check-python-style
      name: Check Python Style
      entry: ./scripts/check-python-style.sh
      language: system
      types: [python]
      pass_filenames: false
```

### GitHub Actions

```yaml
- name: Check C++ Style
  run: ./scripts/check-source-style.sh

- name: Check Python Style
  run: |
    pip install black isort flake8
    ./scripts/check-python-style.sh
```

## Quick Reference

| Task | Command |
|------|---------|
| Check C++ style | `./scripts/check-source-style.sh` |
| Check Python style | `./scripts/check-python-style.sh` |
| Auto-fix Python formatting | `./scripts/check-python-style.sh --fix` |
| Run static analysis | `./scripts/run-clang-tidy.sh` |
| Generate units docs | `python scripts/generate_units_md.py` |

## Exit Codes

All scripts return:
- `0` - All checks passed
- `1` - Violations found

This allows easy integration into CI/CD pipelines.
