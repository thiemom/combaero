# Testing Guide

Complete guide to testing CombAero at all levels.

## Quick Start

```bash
# Run all tests
make test-all

# Run style checks before building
make pre-build

# Full CI workflow
make ci
```

## Test Categories

### 1. C++ Unit Tests (GoogleTest)

**Location**: `tests/`

**Run:**
```bash
make build
make test

# Or directly with ctest
cd build
ctest --output-on-failure
```

**Tests**:
- Humid air calculations
- Orifice flow
- Thermodynamic properties
- Transport properties
- Water and ice equation accuracy

### 2. Python Unit Tests (pytest)

**Location**: `python/tests/`

**Run:**
```bash
make test-python

# Or directly
cd python/tests
python -m pytest -v
```

**Tests**:
- Compressible flow
- Equilibrium (reforming, WGS, SMR)
- Equivalence ratio
- Incompressible flow
- Orifice Cd
- Thermo properties
- Transport properties

### 3. Cantera Validation Tests (pytest + Poetry)

**Location**: `cantera_validation_tests/`

**Run:**
```bash
make test-validation

# Or with specific targets
cd cantera_validation_tests
make test-combustion    # Combustion only
make test-equilibrium   # Equilibrium only
make test-transport     # Transport only
make test-mixing        # Mixing only
```

**Tests** (38 total):
- **Combustion** (12 tests): Complete combustion, adiabatic flame temperature
- **Equilibrium** (6 tests): WGS isothermal/adiabatic, Kp validation
- **Transport** (12 tests): Viscosity, thermal conductivity, Prandtl number
- **Mixing** (8 tests): Stream mixing, enthalpy conservation

**Validation Results**:
- Combustion: < 5 K deviation (NASA-7 vs NASA-9)
- Equilibrium: < 0.002% composition, 0.0 K temperature
- Transport: 10-22% (expected, different correlations)

## Code Style Checking

### C++ Style

```bash
make style-check-cpp

# Or directly
./scripts/check-source-style.sh
./scripts/check-source-style.sh --verbose
```

**Checks**:
- Non-ASCII characters
- Block comments (only line comments allowed)
- M_PI usage without math_constants.h

### Python Style

```bash
make style-check-python

# Auto-fix issues
make style-fix

# Or directly
./scripts/check-python-style.sh
./scripts/check-python-style.sh --fix
./scripts/check-python-style.sh --verbose
```

**Checks**:
- Non-ASCII characters
- PEP 8 compliance (flake8)
- Code formatting (black)
- Import sorting (isort)
- Type hints (mypy)

### All Style Checks

```bash
make style-check
```

## CI/CD Workflows

### GitHub Actions

Three workflows in `.github/workflows/`:

1. **ci.yml** - Main CI (C++ build + tests)
   - Runs on: Ubuntu, macOS, Windows
   - Checks: C++ style, Python style
   - Builds: Release configuration
   - Tests: C++ unit tests (ctest)

2. **validation-tests.yml** - Python validation
   - Runs on: Ubuntu, macOS
   - Python: 3.11, 3.12
   - Jobs:
     - Python style check
     - Cantera validation tests
     - Python unit tests
   - Artifacts: Test results

3. **cmake.yml** - CMake build verification
   - Basic build and test
   - Accuracy tests

### Pre-commit Hooks

See `.pre-commit-config.yaml`:
- Black formatting
- Flake8 linting
- Validation tests (on relevant changes)

## Make Targets Reference

### Root Makefile

| Target | Description |
|--------|-------------|
| `make build` | Build C++ library (Release) |
| `make build-debug` | Build with debug symbols |
| `make test` | Run C++ tests |
| `make test-python` | Run Python unit tests |
| `make test-validation` | Run Cantera validation tests |
| `make test-all` | Run all tests |
| `make style-check` | Check all code style |
| `make style-check-cpp` | Check C++ style only |
| `make style-check-python` | Check Python style only |
| `make style-fix` | Auto-fix Python style |
| `make pre-build` | Style checks before building |
| `make ci` | Full CI workflow |
| `make clean` | Clean all build artifacts |
| `make install-deps` | Install Python dependencies |

### Validation Tests Makefile

```bash
cd cantera_validation_tests
```

| Target | Description |
|--------|-------------|
| `make install` | Install Poetry dependencies |
| `make test` | Run all validation tests |
| `make test-verbose` | Verbose output |
| `make test-parallel` | Parallel execution |
| `make test-combustion` | Combustion tests only |
| `make test-equilibrium` | Equilibrium tests only |
| `make test-transport` | Transport tests only |
| `make test-mixing` | Mixing tests only |
| `make clean` | Clean cache files |

## Test Dependencies

### C++ Tests
- CMake ≥ 3.14
- GoogleTest (fetched automatically)
- C++17 compiler

### Python Unit Tests
- Python ≥ 3.11
- NumPy ≥ 1.24
- pytest ≥ 8.0

### Validation Tests
- Python ≥ 3.11
- Cantera ≥ 3.0
- NumPy ≥ 1.26
- pytest ≥ 8.0
- Poetry (for dependency management)

### Style Tools
- black ≥ 24.0
- isort ≥ 5.13
- flake8 ≥ 7.0
- mypy ≥ 1.8 (optional)

**Install:**
```bash
pip install black isort flake8 mypy
# or: uv pip install black isort flake8 mypy
```

See `scripts/INSTALL.md` for detailed installation options.

## Continuous Integration

### Local CI Simulation

```bash
# Full CI workflow
make ci

# Or step by step
make style-check
make build
make test-all
```

### GitHub Actions Triggers

**ci.yml**:
- Push to main/master
- Pull requests to main/master

**validation-tests.yml**:
- Push to main/master (if src/, include/, python/, or validation tests changed)
- Pull requests (same path filters)

## Test Results and Artifacts

### Validation Test Results

After running validation tests, check:
- `cantera_validation_tests/TEST_RESULTS.md` - Complete combustion results
- `cantera_validation_tests/EQUILIBRIUM_TEST_RESULTS.md` - Equilibrium results
- `cantera_validation_tests/FINAL_TEST_REPORT.md` - Summary

### GitHub Actions Artifacts

Validation test results are uploaded as artifacts:
- `validation-results-ubuntu-py3.11`
- `validation-results-ubuntu-py3.12`
- `validation-results-macos-py3.11`
- `validation-results-macos-py3.12`

## Troubleshooting

### Validation Tests Fail

1. **CombAero not found**:
   ```bash
   python -m build --wheel
   pip install dist/*.whl
   ```

2. **Cantera not installed**:
   ```bash
   cd cantera_validation_tests
   poetry install
   ```

3. **Import errors**:
   - Check that CombAero wheel is installed in the same environment
   - Verify Poetry virtual environment is activated

### Style Checks Fail

1. **Python style issues**:
   ```bash
   make style-fix  # Auto-fix with black/isort
   ```

2. **C++ style issues**:
   - Fix manually (no auto-fix for C++)
   - See error output for specific violations

### Build Fails

1. **CMake version**:
   - Requires CMake ≥ 3.14
   - Update: `pip install --upgrade cmake`

2. **Compiler issues**:
   - Requires C++17 support
   - GCC ≥ 7, Clang ≥ 5, MSVC ≥ 2017

## Best Practices

1. **Before committing**:
   ```bash
   make pre-build  # Style checks
   make test-all   # All tests
   ```

2. **During development**:
   ```bash
   make build-debug  # Debug build
   make test         # Quick C++ tests
   ```

3. **Before PR**:
   ```bash
   make ci  # Full CI workflow
   ```

4. **Validation changes**:
   ```bash
   cd cantera_validation_tests
   make test-verbose  # See detailed output
   ```

## Documentation

- **API Reference**: `docs/API_REFERENCE.md`
- **Units**: `docs/UNITS.md`
- **Validation**: `cantera_validation_tests/README.md`
- **Scripts**: `scripts/README.md`
- **Installation**: `scripts/INSTALL.md`
