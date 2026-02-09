# Cantera Validation Test Suite - Summary

## Overview

Comprehensive validation test suite comparing CombAero calculations against Cantera as a reference implementation. Tests use GRI-Mech 3.0 mechanism for consistency.

**Important**: Complete combustion tests use a **restricted species gas phase** (N2, O2, AR, CO2, H2O, fuels) to match CombAero's complete combustion model, which produces only CO2 and H2O without equilibrium species (CO, H2, OH, etc.).

## Test Coverage

### Combustion Validation (test_combustion_validation.py)

#### TestCompleteCombustion
- `test_methane_air_stoichiometric`: CH4 + air at φ=1.0, validates Tad and products
- `test_methane_air_lean`: CH4 + air at φ=0.8, validates excess O2
- `test_propane_air_stoichiometric`: C3H8 + air at φ=1.0
- `test_hydrogen_air_stoichiometric`: H2 + air at φ=1.0, validates H2O production
- `test_temperature_variation`: Combustion at T_in = 300, 400, 500, 600 K
- `test_pressure_variation`: Combustion at P_in = 1, 2, 5, 10 bar

#### TestOxygenRequirement
- `test_methane_oxygen_requirement`: Validates CH4 + 2 O2 → CO2 + 2 H2O
- `test_propane_oxygen_requirement`: Validates C3H8 + 5 O2 → 3 CO2 + 4 H2O
- `test_hydrogen_oxygen_requirement`: Validates 2 H2 + O2 → 2 H2O

#### TestEquivalenceRatio
- `test_stoichiometric_mixture`: φ=1.0 round-trip validation
- `test_lean_mixture`: φ=0.8 validation
- `test_rich_mixture`: φ=1.2 validation

**Total: 12 combustion tests**

### Mixing Validation (test_mixing_validation.py)

#### TestStreamMixing
- `test_two_stream_mixing_equal_mass`: Equal mass flow mixing, validates T and mdot
- `test_two_stream_mixing_unequal_mass`: 10:0.5 mass ratio mixing
- `test_three_stream_mixing`: Air + fuel + steam mixing
- `test_enthalpy_conservation`: Validates H conservation in mixing

#### TestDensityCalculation
- `test_air_density`: Air at standard conditions
- `test_methane_density`: Pure CH4 density
- `test_density_temperature_variation`: ρ(T) at 300, 500, 1000, 1500, 2000 K
- `test_density_pressure_variation`: ρ(P) at 1, 2, 5, 10 bar

**Total: 8 mixing tests**

### Transport Validation (test_transport_validation.py)

#### TestTransportProperties
- `test_air_viscosity`: Air μ at standard conditions
- `test_viscosity_temperature_variation`: μ(T) at 300, 500, 1000, 1500, 2000 K
- `test_methane_viscosity`: Pure CH4 viscosity
- `test_air_thermal_conductivity`: Air k at standard conditions
- `test_thermal_conductivity_temperature_variation`: k(T) at 300-2000 K
- `test_methane_thermal_conductivity`: Pure CH4 thermal conductivity
- `test_prandtl_number`: Pr calculation validation
- `test_mixture_transport_properties`: CH4-air mixture at φ=1.0
- `test_high_temperature_transport`: Post-combustion transport properties

#### TestThermodynamicProperties
- `test_air_cp`: Specific heat capacity
- `test_cp_temperature_variation`: Cp(T) at 300-2000 K
- `test_enthalpy`: Enthalpy calculation

**Total: 12 transport tests**

## Grand Total: 32 Validation Tests

## Test Execution Time

- Sequential: ~30-60 seconds
- Parallel (`-n auto`): ~10-20 seconds

## Tolerances

| Property | Tolerance | Rationale |
|----------|-----------|-----------|
| Temperature | ±5 K | NASA-7 vs NASA-9 polynomial differences |
| Mole fractions | ±0.01 | 1% absolute, accounts for equilibrium solver differences |
| Enthalpy | ±1% | Relative, polynomial integration differences |
| Transport | ±5% | Correlation and mixing rule differences |
| Density | ±1% | Relative, ideal gas law consistency |

## Species Coverage

Tests validate calculations for:
- **Fuels**: CH4, C3H8, H2
- **Oxidizer**: Air (N2, O2, AR, CO2)
- **Products**: CO2, H2O, N2, O2 (lean)
- **Inerts**: AR

## Temperature Range

- Inlet: 300-600 K
- Post-combustion: 1500-2300 K (depending on fuel and φ)
- Transport validation: 300-2000 K

## Pressure Range

- 1-10 bar (101325-1000000 Pa)

## Test Quality Metrics

- **Coverage**: All major combustion, mixing, and transport functions
- **Parametric**: Multiple T, P, φ, fuel types
- **Reference**: Industry-standard Cantera with GRI-Mech 3.0
- **Robustness**: Tolerances account for legitimate model differences
- **Speed**: Fast enough for pre-commit hooks (<1 minute)

## Continuous Integration

Tests are designed for:
- Pre-commit hooks (local development)
- CI/CD pipelines (GitHub Actions, GitLab CI, etc.)
- Regression testing after code changes
- Validation after data updates (NASA polynomials, transport data)

## Maintenance

Update tests when:
- Adding new species to CombAero
- Changing thermodynamic data (NASA polynomials)
- Modifying combustion algorithms
- Updating transport property correlations
- Changing mixing logic

## Known Limitations

1. **Equilibrium vs Complete Combustion**: Tests use complete combustion (no WGS), not full equilibrium
2. **Mechanism Differences**: GRI-Mech 3.0 may have different species than CombAero's 14-species set
3. **Transport Correlations**: Different mixing rules may cause small differences within tolerance
4. **Polynomial Fits**: NASA-7 (Cantera) vs NASA-9 (CombAero) can differ by a few Kelvin

These limitations are acceptable and within engineering tolerances for the intended use cases.
