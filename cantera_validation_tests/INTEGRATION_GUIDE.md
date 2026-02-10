# Cantera Validation Tests - Integration Guide

## Overview

This test suite validates CombAero's combustion, mixing, and transport property calculations against Cantera as a reference implementation. The tests are designed to be part of the pre-commit workflow to ensure all changes maintain consistency with established thermodynamic models.

## Architecture

### Test Structure

```
cantera_validation_tests/
├── pyproject.toml              # Poetry dependencies (Cantera, pytest)
├── conftest.py                 # Shared fixtures and configuration
├── test_combustion_validation.py   # Combustion tests
├── test_mixing_validation.py       # Stream mixing tests
├── test_transport_validation.py    # Transport property tests
├── run_tests.sh                # Convenience script
├── Makefile                    # Make targets
└── README.md                   # Documentation
```

### Complete Combustion Validation Approach

**Critical**: CombAero's `complete_combustion()` produces only CO2 and H2O (no equilibrium species like CO, H2, OH, etc.). To validate this correctly, we use a **restricted species gas phase** in Cantera:

```python
# Create gas with only complete combustion species
species = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
complete_species = [species[S] for S in ("N2", "O2", "AR", "CO2", "H2O", "CH4", "C3H8", "H2")]
gas = ct.Solution(thermo="ideal-gas", species=complete_species)

# Now equilibrate - restricted to complete combustion products only
gas.TP = T_in, P_in
gas.equilibrate("HP")
```

This ensures Cantera produces the same products as CombAero's complete combustion model.

### Test Categories

1. **Combustion Validation** (`test_combustion_validation.py`)
   - Complete combustion calculations
   - Adiabatic flame temperature (Tad)
   - Product composition (CO2, H2O, O2)
   - Oxygen requirement calculations
   - Equivalence ratio calculations
   - Temperature and pressure variations

2. **Mixing Validation** (`test_mixing_validation.py`)
   - Two-stream mixing (equal and unequal mass flows)
   - Three-stream mixing
   - Enthalpy conservation
   - Density calculations at various conditions

3. **Transport Validation** (`test_transport_validation.py`)
   - Viscosity calculations
   - Thermal conductivity calculations
   - Prandtl number
   - Temperature and pressure variations
   - High-temperature transport properties
   - Thermodynamic properties (Cp, enthalpy)

## Setup Instructions

### 1. Install Dependencies

```bash
cd cantera_validation_tests
poetry install
```

This will install:
- Cantera (≥3.0.0)
- NumPy (≥1.26.0)
- pytest (≥8.0.0)
- pytest-xdist (for parallel execution)

### 2. Verify Installation

```bash
poetry run python -c "import cantera; print(cantera.__version__)"
```

### 3. Build CombAero Python Bindings

Before running tests, ensure CombAero Python bindings are built:

```bash
cd ..
python -m build --wheel
pip install dist/combaero-*.whl
```

## Running Tests

### Basic Usage

```bash
# Run all tests
cd cantera_validation_tests
poetry run pytest

# Run with verbose output
poetry run pytest -v

# Run specific test file
poetry run pytest test_combustion_validation.py -v

# Run specific test class
poetry run pytest test_combustion_validation.py::TestCompleteCombustion -v

# Run specific test
poetry run pytest test_combustion_validation.py::TestCompleteCombustion::test_methane_air_stoichiometric -v
```

### Using Make

```bash
make test              # Run all tests
make test-verbose      # Run with verbose output
make test-parallel     # Run tests in parallel
make clean             # Clean cache files
```

### Using Shell Script

```bash
./run_tests.sh
```

## Pre-commit Integration

The tests are integrated into the pre-commit workflow. To set up pre-commit hooks:

```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install

# Run manually
pre-commit run --all-files
```

The Cantera validation tests will run automatically on commits that modify:
- `src/` (C++ implementation)
- `include/` (C++ headers)
- `python/` (Python bindings)
- `cantera_validation_tests/` (test files)

## Tolerance Configuration

Tests use the following tolerances (defined in `conftest.py`):

| Property | Tolerance | Notes |
|----------|-----------|-------|
| Temperature | ±5 K | Accounts for polynomial fit differences |
| Mole fractions | ±0.01 | 1% absolute difference |
| Enthalpy | ±1% | Relative difference |
| Transport properties | ±5% | Typical for correlation differences |
| Density | ±1% | Relative difference |

These tolerances are conservative and account for:
- Different polynomial fits (NASA-7 vs NASA-9)
- Different transport property correlations
- Numerical precision differences

## Understanding Test Failures

### Temperature Mismatches

If adiabatic flame temperature tests fail:
1. Check NASA polynomial coefficients in `thermo_transport_data.h`
2. Verify species composition mapping
3. Ensure consistent reference states (298.15 K)

### Composition Mismatches

If product composition tests fail:
1. Verify complete combustion logic in `combustion.cpp`
2. Check species balance equations
3. Ensure proper normalization of mole fractions

### Transport Property Mismatches

If transport property tests fail:
1. Check Lennard-Jones parameters in `thermo_transport_data.h`
2. Verify mixing rules (Wilke's method for viscosity, etc.)
3. Ensure temperature-dependent correlations are correct

## Debugging Tips

### Enable Verbose Output

```bash
poetry run pytest -v -s  # -s shows print statements
```

### Run Single Test with Full Traceback

```bash
poetry run pytest test_combustion_validation.py::TestCompleteCombustion::test_methane_air_stoichiometric -vv --tb=long
```

### Compare Values Directly

Add print statements in tests to see exact values:

```python
print(f"CombAero: T={burned_cb.T:.2f} K")
print(f"Cantera:  T={gri30_gas.T:.2f} K")
print(f"Diff:     {abs(burned_cb.T - gri30_gas.T):.2f} K")
```

### Check Species Mapping

Verify species names match between CombAero and Cantera:

```python
print("CombAero species:", cb.species_names())
print("Cantera species:", gri30_gas.species_names)
```

## Continuous Integration

For CI/CD pipelines, add the following step:

```yaml
- name: Run Cantera Validation Tests
  run: |
    cd cantera_validation_tests
    poetry install
    poetry run pytest -v --tb=short
```

## Extending the Test Suite

### Adding New Test Cases

1. Add test method to appropriate test class:

```python
def test_new_case(self, combaero, cantera, gri30_gas, species_mapping, tolerance_config):
    """Test description."""
    cb = combaero

    # Setup
    X_mix = ...

    # CombAero calculation
    result_cb = cb.some_function(...)

    # Cantera calculation
    self.set_cantera_composition(gri30_gas, X_mix, species_mapping)
    gri30_gas.TP = T, P
    result_ct = ...

    # Assertion
    assert abs(result_cb - result_ct) < tolerance_config["..."]
```

2. Run the new test:

```bash
poetry run pytest -v -k test_new_case
```

### Adding New Test Files

1. Create `test_new_module_validation.py`
2. Import fixtures from `conftest.py`
3. Follow existing test structure
4. Add to pytest discovery (automatic if named `test_*.py`)

## Performance Considerations

### Parallel Execution

Run tests in parallel to reduce execution time:

```bash
poetry run pytest -n auto  # Uses all CPU cores
poetry run pytest -n 4     # Uses 4 cores
```

### Selective Test Execution

Run only fast tests during development:

```bash
poetry run pytest -m "not slow"
```

Mark slow tests with decorator:

```python
@pytest.mark.slow
def test_expensive_calculation(...):
    ...
```

## Troubleshooting

### Cantera Not Found

```bash
# Verify Cantera installation
poetry run python -c "import cantera"

# If fails, reinstall
poetry install --no-cache
```

### CombAero Not Found

```bash
# Rebuild and install CombAero
cd ..
python -m build --wheel
pip install --force-reinstall dist/combaero-*.whl
```

### GRI-Mech 3.0 Not Found

Cantera should include GRI-Mech 3.0 by default. If not:

```bash
# Check available mechanisms
poetry run python -c "import cantera; print(cantera.Solution('gri30.yaml'))"
```

### Test Hangs or Times Out

Some equilibrium calculations can be slow. Increase timeout:

```bash
poetry run pytest --timeout=300  # 5 minutes per test
```

## Best Practices

1. **Run tests before committing**: Ensure all tests pass locally
2. **Update tolerances carefully**: Document reasons for tolerance changes
3. **Add tests for new features**: Maintain validation coverage
4. **Keep tests independent**: Each test should be self-contained
5. **Use descriptive names**: Test names should explain what is being validated
6. **Document expected behavior**: Add docstrings to test methods

## Support

For issues or questions:
1. Check this guide and README.md
2. Review test output and error messages
3. Compare with Cantera documentation
4. Open an issue with detailed error information
