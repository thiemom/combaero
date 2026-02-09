# Test Units Reference

## Overview

This document clarifies the units used in all validation tests and what is being compared between CombAero and Cantera.

## Combustion Tests

### Temperature Comparisons

**Units**: Kelvin [K]

**What's compared**:
- **CombAero**: Adiabatic flame temperature from NASA-9 polynomials, complete combustion model (CO2 + H2O only)
- **Cantera**: Temperature from enthalpy balance using NASA-7 polynomials with same product composition

**Method**: Find T where H(T_products, X_products) = H(T_reactants, X_reactants)

**Expected differences**: < 5 K due to NASA-7 vs NASA-9 polynomial fits

### Oxygen Requirement

**Units**: mol O2 / mol fuel [dimensionless ratio]

**What's compared**:
- Stoichiometric oxygen requirement from molecular formula
- Example: CH4 + 2 O2 → CO2 + 2 H2O gives 2.0 mol O2/mol CH4

**Expected differences**: Exact (0.0) - pure stoichiometry

### Equivalence Ratio

**Units**: Dimensionless [φ]

**What's compared**:
- Round-trip test: Set φ → calculate mixture → calculate φ from mixture
- Tests consistency of equivalence ratio calculations

**Expected differences**: < 1e-6 (numerical precision only)

## Transport Property Tests

### Dynamic Viscosity

**Units**: Pascal-second [Pa·s] or [kg/(m·s)]

**What's compared**:
- **CombAero**: Sutherland or polynomial correlations with specific Lennard-Jones parameters
- **Cantera**: GRI-Mech 3.0 transport data with different Lennard-Jones parameters

**Expected differences**: 10-22% due to:
- Different Lennard-Jones σ and ε/k parameters
- Different temperature-dependent correlations
- Different mixing rules for mixtures

### Thermal Conductivity

**Units**: Watt per meter-Kelvin [W/(m·K)]

**What's compared**:
- **CombAero**: Eucken or modified Eucken formula
- **Cantera**: GRI-Mech 3.0 transport correlations

**Expected differences**: 15-20% due to:
- Different correlations (Eucken vs others)
- Different mixing rules
- Different transport data sources

### Prandtl Number

**Units**: Dimensionless [Pr]

**Formula**: Pr = (Cp · μ) / k

**What's compared**:
- Calculated from each code's own Cp, μ, and k values
- Differences propagate from transport property differences

**Expected differences**: 10-15% (combined effect of Cp, μ, k differences)

## Thermodynamic Property Tests

### Molar Heat Capacity (Cp)

**Units**: Joule per mole-Kelvin [J/(mol·K)]

**What's compared**:
- **CombAero**: NASA-9 polynomials (7 coefficients, single range 200-6000 K)
- **Cantera**: NASA-7 polynomials (7 coefficients per range, with temperature break point)

**Note**: Cantera returns J/(kmol·K), converted to J/(mol·K) by dividing by 1000

**Expected differences**: < 1% for well-fitted polynomials

### Molar Enthalpy

**Units**: Joule per mole [J/mol]

**Formula**: H(T) = H°(298.15K) + ∫[298.15 to T] Cp(T) dT

**What's compared**:
- Integration of Cp polynomials from reference temperature
- **CombAero**: NASA-9 integration
- **Cantera**: NASA-7 integration

**Note**: Cantera returns J/kmol, converted to J/mol by dividing by 1000

**Expected differences**: < 1.5% due to polynomial integration differences

## Mixing Tests

### Mixed Stream Temperature

**Units**: Kelvin [K]

**What's compared**:
- Temperature after adiabatic mixing at constant pressure
- Both codes use enthalpy balance: H_mixed = Σ(mdot_i · H_i) / Σ(mdot_i)

**Expected differences**: < 5 K (NASA-7 vs NASA-9 in enthalpy calculations)

### Mass Flow Rate

**Units**: Kilogram per second [kg/s]

**What's compared**:
- Total mass flow rate after mixing
- Should be exact: mdot_total = Σ(mdot_i)

**Expected differences**: < 1e-10 (numerical precision only)

### Density

**Units**: Kilogram per cubic meter [kg/m³]

**Formula**: ρ = P·MW / (R·T) (ideal gas law)

**What's compared**:
- **CombAero**: Uses mixture molecular weight from mole fractions
- **Cantera**: Same ideal gas calculation

**Expected differences**: < 0.1% (only from MW calculation precision)

## Unit Conversion Notes

### Cantera to CombAero Conversions

| Property | Cantera Units | CombAero Units | Conversion |
|----------|---------------|----------------|------------|
| Cp | J/(kmol·K) | J/(mol·K) | Divide by 1000 |
| Enthalpy | J/kmol | J/mol | Divide by 1000 |
| Viscosity | Pa·s | Pa·s | No conversion |
| Thermal conductivity | W/(m·K) | W/(m·K) | No conversion |
| Temperature | K | K | No conversion |
| Pressure | Pa | Pa | No conversion |
| Density | kg/m³ | kg/m³ | No conversion |

### Molar vs Mass-Specific

**CombAero convention**: Molar basis for thermodynamic properties
- Cp: J/(mol·K)
- H: J/mol
- S: J/(mol·K)

**Cantera provides both**:
- `cp_mole`: J/(kmol·K) → convert to J/(mol·K)
- `cp_mass`: J/(kg·K) → for Prandtl number calculation
- `enthalpy_mole`: J/kmol → convert to J/mol
- `enthalpy_mass`: J/kg → for enthalpy balance

## Tolerance Specifications

Based on measured deviations:

| Property | Tolerance | Measured Max | Reason |
|----------|-----------|--------------|--------|
| Temperature | 5.0 K | 4.6 K | NASA-7 vs NASA-9 polynomials |
| Mole fraction | 0.01 | 0.0 | Stoichiometry exact |
| Enthalpy | 1.5% | 1.02% | Polynomial integration differences |
| Transport | 25% | 21.5% | Different correlations and data |
| Density | 1% | < 0.1% | Ideal gas law, very accurate |

## Physical Interpretation

### Why Transport Properties Differ More

Transport properties (μ, k) depend on:
1. **Molecular parameters**: Lennard-Jones σ and ε/k
2. **Correlations**: Sutherland, power law, polynomial
3. **Mixing rules**: Wilke, Herning-Zipperer, others
4. **Data sources**: Different experimental databases

These are **legitimate differences** between implementations, not errors.

### Why Thermodynamic Properties Agree Better

Thermodynamic properties (Cp, H, S) depend on:
1. **Polynomial fits**: NASA-7 vs NASA-9 to same JANAF data
2. **Integration**: Numerical integration of polynomials
3. **Reference states**: Both use 298.15 K

Differences are **only from polynomial fitting**, much smaller.

### Why Stoichiometry is Exact

Element balance is exact:
- C, H, O, N atoms conserved
- No polynomial fitting involved
- Pure integer arithmetic

## Summary

All tests include explicit unit comments in assertions to clarify:
1. What physical quantity is being compared
2. What units are used
3. What method each code uses
4. What differences are expected and why

This makes the validation transparent and reproducible.
