# Quick Start Guide - Cantera Validation Tests

## 5-Minute Setup

### 1. Install Dependencies

```bash
cd cantera_validation_tests
poetry install
```

This installs Cantera, pytest, and all required dependencies in an isolated virtual environment.

### 2. Build CombAero (if not already done)

```bash
cd ..
python -m build --wheel
pip install dist/combaero-*.whl
```

### 3. Run Tests

```bash
cd cantera_validation_tests
poetry run pytest -v
```

## What Gets Tested?

### Combustion Tests
- ✓ Adiabatic flame temperature for CH4, C3H8, H2 + air
- ✓ Product composition (CO2, H2O, O2)
- ✓ Lean and rich combustion
- ✓ Temperature and pressure variations
- ✓ Oxygen requirement calculations
- ✓ Equivalence ratio calculations

### Mixing Tests
- ✓ Two-stream mixing (equal and unequal mass flows)
- ✓ Three-stream mixing
- ✓ Enthalpy conservation
- ✓ Density at various T and P

### Transport Tests
- ✓ Viscosity (air, methane, mixtures)
- ✓ Thermal conductivity
- ✓ Prandtl number
- ✓ Temperature variations (300-2000 K)
- ✓ High-temperature post-combustion properties

## Expected Results

All tests should pass with tolerances:
- Temperature: ±5 K
- Mole fractions: ±0.01
- Transport properties: ±5%

## Common Issues

### "Cantera not found"
```bash
poetry install --no-cache
```

### "CombAero not found"
```bash
cd .. && python -m build --wheel && pip install --force-reinstall dist/combaero-*.whl
```

### Tests fail
Check that you're using the latest CombAero build and that NASA polynomial data is up to date.

## Next Steps

- Read [README.md](README.md) for detailed documentation
- See [INTEGRATION_GUIDE.md](INTEGRATION_GUIDE.md) for CI/CD integration
- Run `make help` for available commands
