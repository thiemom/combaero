# Cantera Validation Tests

This test suite validates CombAero's combustion, mixing, and transport property calculations against Cantera as a reference implementation.

## Purpose

These tests ensure that CombAero's complete combustion, stream mixing, and transport property calculations match Cantera's results within acceptable tolerances. This is critical for:

- Verifying thermodynamic consistency
- Validating adiabatic flame temperature calculations
- Ensuring accurate species composition predictions
- Confirming transport property correlations

## Setup

Install dependencies using Poetry:

```bash
cd cantera_validation_tests
poetry install
```

## Running Tests

```bash
# Run all tests
poetry run pytest

# Run with verbose output
poetry run pytest -v

# Run specific test file
poetry run pytest test_combustion_validation.py -v

# Run in parallel
poetry run pytest -n auto

# Run with coverage
poetry run pytest --cov=. --cov-report=html
```

## Validation Methodology

**Complete Combustion**: Uses enthalpy balance (NOT equilibrium) to find adiabatic flame temperature. This correctly validates CombAero's complete combustion model (CO2 + H2O only) without equilibrium effects (no CO, H2, dissociation).

**WGS Equilibrium**: Uses Cantera's equilibrium solver with restricted species (CO, H2O, CO2, H2) to match CombAero's partial equilibrium model. Validates both isothermal and adiabatic equilibrium, plus equilibrium constant calculations.

## Test Structure

- `test_combustion_validation.py`: Complete combustion, adiabatic flame temperature, product composition (12 tests)
- `test_equilibrium_validation.py`: WGS equilibrium (isothermal, adiabatic), equilibrium constants (6 tests)
- `test_mixing_validation.py`: Stream mixing, enthalpy balance, mixture properties (8 tests)
- `test_transport_validation.py`: Viscosity, thermal conductivity, Prandtl number (12 tests)

**Total**: 38 validation tests

## Tolerances

Based on measured deviations (NASA-7 vs NASA-9 polynomials):

- **Temperature**: ±5 K (measured max: 4.6 K for C3H8 combustion)
- **Mole fractions**: ±0.01 (1% absolute, stoichiometry exact)
- **Enthalpy**: ±1.5% relative (measured max: 1.02%)
- **Transport properties**: ±25% relative (measured max: 21.5% at high T)
- **Density**: ±1% relative (ideal gas law, very accurate)
- **Equilibrium composition**: ±0.0001 (0.01% absolute, measured max: 0.002%)
- **Equilibrium temperature**: ±1 K (measured max: 0.0 K)
- **Equilibrium constant (Kp)**: ±0.2% relative (measured max: 0.12%)

See `FINAL_TEST_REPORT.md` and `EQUILIBRIUM_TEST_RESULTS.md` for detailed validation results.

## Pre-commit Integration

These tests are integrated into the pre-commit workflow to ensure all changes maintain consistency with Cantera.
