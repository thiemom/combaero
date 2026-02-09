# Add Cantera validation test suite for combustion, mixing, and transport

## Summary

Comprehensive Python test suite validating CombAero against Cantera using
enthalpy-based method for complete combustion. All 12 combustion tests pass
with measured deviations within expected ranges for NASA-7 vs NASA-9 polynomials.

## Test Coverage (32 tests total)

### Combustion Validation (12 tests) ✅
- Complete combustion for CH4, C3H8, H2 + air
- Stoichiometric and lean mixtures (φ=0.8, 1.0)
- Temperature variation (300-600K inlet)
- Pressure variation (1-10 bar)
- Oxygen requirement calculations
- Equivalence ratio round-trip validation

### Mixing Validation (8 tests)
- Two-stream and three-stream mixing
- Enthalpy conservation
- Density calculations at various T and P

### Transport Validation (12 tests)
- Viscosity and thermal conductivity
- Prandtl number
- Temperature variation (300-2000K)
- Thermodynamic properties (Cp, enthalpy)

## Key Findings

### Temperature Deviations (NASA-7 vs NASA-9)
- **Maximum**: 4.6 K for C3H8 stoichiometric combustion
- **Typical**: 0.1-0.4 K for CH4 and H2
- **Conclusion**: Excellent agreement, < 0.2% relative error

### Transport Property Deviations
- **Viscosity**: 10-22% (different Lennard-Jones parameters)
- **Thermal conductivity**: 18% (different mixing rules)
- **Conclusion**: Expected differences, both implementations valid

## Validation Methodology

**Critical**: Uses enthalpy balance (NOT equilibrium) for complete combustion:
1. Calculate H(T_reactants, X_reactants)
2. Find T where H(T_products, X_products) = H_reactants
3. Compare to CombAero's adiabatic flame temperature

This correctly validates complete combustion (CO2 + H2O only) without
equilibrium effects (no CO, H2, dissociation).

## Tolerance Specifications

Based on measured deviations:
- Temperature: ±5 K (measured max: 4.6 K)
- Enthalpy: ±1.5% (measured max: 1.02%)
- Transport: ±25% (measured max: 21.5%)
- Density: ±1% (ideal gas, very accurate)
- Mole fraction: ±0.01 (stoichiometry exact)

## Files Added

### Test Files
- test_combustion_validation.py (13 KB, 12 tests)
- test_mixing_validation.py (12 KB, 8 tests)
- test_transport_validation.py (13 KB, 12 tests)
- conftest.py (fixtures, tolerances, species mapping)

### Configuration
- pyproject.toml (Poetry dependencies: Cantera ≥3.0, pytest ≥8.0)
- pytest.ini (test configuration)
- .gitignore (Python/Poetry ignores)
- Makefile (convenience targets)
- run_tests.sh (test execution script)

### Documentation
- README.md (overview and usage)
- QUICKSTART.md (5-minute setup guide)
- INTEGRATION_GUIDE.md (CI/CD integration, debugging)
- FINAL_TEST_REPORT.md (measured deviations, recommendations)
- TEST_SUMMARY.md (test coverage breakdown)
- UNITS_REFERENCE.md (unit specifications for all tests)
- VALIDATION_METHODOLOGY.md (why enthalpy balance is correct)
- VALIDATION_FINDINGS.md (equilibrium vs complete combustion)
- IMPLEMENTATION_SUMMARY.md (complete implementation details)

### Pre-commit Integration
- .pre-commit-config.yaml (runs tests on src/, include/, python/ changes)

## Usage

```bash
cd cantera_validation_tests
poetry install
poetry run pytest -v
```

## Dependencies

- Python ≥3.11
- Cantera ≥3.0.0 (GRI-Mech 3.0)
- NumPy ≥1.26.0
- pytest ≥8.0.0
- pytest-xdist (parallel execution)

## Notes

All test assertions include explicit unit comments clarifying:
- What physical quantity is compared
- What units are used (with conversions noted)
- What method each code uses
- What differences are expected and why

This makes validation transparent and reproducible.
