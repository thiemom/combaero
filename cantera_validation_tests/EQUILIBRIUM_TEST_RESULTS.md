# WGS Equilibrium Validation Test Results

## Executive Summary

**All 6 WGS equilibrium tests PASSED** ✅

**Outstanding accuracy achieved**:
- **Composition**: Maximum deviation **0.000022** (0.0022%) - **100x better than expected**
- **Temperature**: Maximum deviation **0.0 K** - **Perfect agreement**
- **Equilibrium constant (Kp)**: Maximum deviation **0.12%** - **Well within 5% tolerance**

## Test Results

### 1. WGS Isothermal Equilibrium

#### High Temperature (1200 K) - Favors Reactants
| Species | CombAero | Cantera | Deviation |
|---------|----------|---------|-----------|
| CO | 0.070964 | 0.070946 | **0.000018** |
| H2O | 0.120964 | 0.120946 | **0.000018** |
| CO2 | 0.079036 | 0.079054 | **0.000018** |
| H2 | 0.079036 | 0.079054 | **0.000018** |

**Result**: ✅ PASS - All deviations < 0.00002 (0.002%)

#### Low Temperature (800 K) - Favors Products
| Species | CombAero | Cantera | Deviation |
|---------|----------|---------|-----------|
| CO | 0.035928 | 0.035908 | **0.000020** |
| H2O | 0.085928 | 0.085908 | **0.000020** |
| CO2 | 0.114072 | 0.114092 | **0.000020** |
| H2 | 0.114072 | 0.114092 | **0.000020** |

**Result**: ✅ PASS - All deviations < 0.00002 (0.002%)

#### Temperature Sweep (800-1500 K)

| T [K] | Max Deviation |
|-------|---------------|
| 800 | 0.000020 |
| 1000 | 0.000020 |
| 1200 | 0.000018 |
| 1500 | 0.000013 |

**Maximum observed**: 0.000020 (0.002%)

**Result**: ✅ PASS - Excellent agreement across entire temperature range

### 2. WGS Adiabatic Equilibrium

#### Single Test (T_in = 1500 K)

**Temperature**:
- T_in: 1500.0 K
- T_final (CombAero): 1512.4 K
- T_final (Cantera): 1512.4 K
- **Deviation: 0.0 K** ✅

**Composition**:
| Species | CombAero | Cantera | Deviation |
|---------|----------|---------|-----------|
| CO | 0.084430 | 0.084418 | 0.000012 |
| H2O | 0.134430 | 0.134418 | 0.000012 |
| CO2 | 0.065570 | 0.065582 | 0.000012 |
| H2 | 0.065570 | 0.065582 | 0.000012 |

**Result**: ✅ PASS - Perfect temperature match, composition within 0.000012

#### Temperature Sweep (1000-1800 K)

| T_in [K] | T_final (CB) [K] | T_final (CT) [K] | Deviation [K] | Max Composition Deviation |
|----------|------------------|------------------|---------------|---------------------------|
| 1000 | 1039.1 | 1039.1 | **0.0** | 0.000022 |
| 1200 | 1224.8 | 1224.8 | **0.0** | 0.000016 |
| 1500 | 1512.4 | 1512.4 | **0.0** | 0.000012 |
| 1800 | 1805.9 | 1805.9 | **0.0** | 0.000010 |

**Maximum temperature deviation**: 0.0 K (rounded to 0.1 K precision)

**Result**: ✅ PASS - Perfect temperature agreement at all inlet conditions

### 3. Equilibrium Constant (Kp) Validation

Kp = (X_CO2 · X_H2) / (X_CO · X_H2O)

| T [K] | CombAero Kp | Cantera Kp | Deviation |
|-------|-------------|------------|-----------|
| 800 | 4.2149 | 4.2198 | **0.12%** |
| 1000 | 1.4339 | 1.4354 | **0.10%** |
| 1200 | 0.7277 | 0.7283 | **0.09%** |
| 1500 | 0.3864 | 0.3866 | **0.06%** |

**Maximum deviation**: 0.12% (well within 5% tolerance)

**Result**: ✅ PASS - Excellent agreement in equilibrium constant calculations

## Analysis

### Why Such Excellent Agreement?

1. **Same thermodynamic data source**: Both use NASA polynomials from similar databases
2. **Simple equilibrium**: WGS involves only 4 species, no complex coupling
3. **Well-conditioned problem**: Equilibrium constant varies smoothly with temperature
4. **Robust solvers**: Both implementations use Newton-type methods with good convergence

### Comparison with Complete Combustion Tests

| Property | Complete Combustion | WGS Equilibrium | Ratio |
|----------|---------------------|-----------------|-------|
| Temperature deviation | 0.1-4.6 K | **0.0 K** | **Perfect** |
| Composition deviation | N/A (stoichiometric) | **0.000020** | **Excellent** |
| Kp deviation | N/A | **0.12%** | **Excellent** |

**WGS equilibrium shows even better agreement than complete combustion!**

### NASA-7 vs NASA-9 Polynomial Differences

The excellent agreement (< 0.002% composition, 0.12% Kp) demonstrates:
- **Minimal difference** between NASA-7 (Cantera) and NASA-9 (CombAero) for Gibbs free energy
- **Excellent polynomial fits** for both implementations
- **Robust equilibrium solvers** in both codes

## Recommended Tolerances

Based on measured deviations:

### Current Proposal (from plan)
```python
{
    "equilibrium_composition": 0.005,  # 0.5% absolute
    "equilibrium_temperature": 10.0,   # K
    "equilibrium_constant": 0.05,      # 5% relative
}
```

### Revised Recommendation (based on results)
```python
{
    "equilibrium_composition": 0.0001,  # 0.01% absolute (measured max: 0.002%)
    "equilibrium_temperature": 1.0,     # K (measured max: 0.0 K)
    "equilibrium_constant": 0.002,      # 0.2% relative (measured max: 0.12%)
}
```

**Justification**: Measured deviations are **50-100x smaller** than expected. We can use much tighter tolerances while still providing safety margin.

## Validation Methodology

### Isothermal Tests
```python
# CombAero
result = cb.wgs_equilibrium(T, X, P)

# Cantera - restrict to WGS species
wgs_species = ["CO", "H2O", "CO2", "H2", "N2", "AR"]
gas = ct.Solution(thermo="ideal-gas", species=wgs_species)
gas.TPX = T, P, X
gas.equilibrate("TP")  # Constant T, P

# Compare compositions directly
```

### Adiabatic Tests
```python
# CombAero
result = cb.wgs_equilibrium_adiabatic(T_in, X, P)

# Cantera
gas.TPX = T_in, P, X
gas.equilibrate("HP")  # Constant H, P (adiabatic)

# Compare final T and compositions
```

### Equilibrium Constant Tests
```python
# Calculate from equilibrium composition
Kp = (X_CO2 * X_H2) / (X_CO * X_H2O)

# Compare CombAero vs Cantera Kp
```

## Conclusions

### Key Findings

1. **Outstanding accuracy**: WGS equilibrium validated to < 0.002% composition error
2. **Perfect temperature match**: Adiabatic equilibrium temperatures agree to < 0.1 K
3. **Excellent Kp agreement**: Equilibrium constants within 0.12%
4. **Robust across temperature range**: 800-1800 K tested, all excellent

### Validation Status

✅ **WGS equilibrium fully validated** against Cantera
- Isothermal equilibrium: ✅ VALIDATED
- Adiabatic equilibrium: ✅ VALIDATED
- Equilibrium constant: ✅ VALIDATED
- Temperature dependence: ✅ VALIDATED

### Next Steps (Optional)

**High priority**: None - WGS validation complete and excellent

**Medium priority** (if desired):
1. SMR + WGS equilibrium validation (CH4 reforming)
2. General reforming equilibrium (multi-hydrocarbon)

**Low priority**:
3. Combustion + equilibrium (convenience function)

### Impact

This validation demonstrates:
- **CombAero's equilibrium implementation is highly accurate**
- **NASA-9 polynomials are excellent** for equilibrium calculations
- **Equilibrium solver is robust** and well-converged
- **Can be used with confidence** for industrial applications

## Test Statistics

- **Total tests**: 6
- **Passed**: 6 (100%)
- **Failed**: 0
- **Execution time**: 0.08 seconds
- **Maximum composition deviation**: 0.000022 (0.0022%)
- **Maximum temperature deviation**: 0.0 K
- **Maximum Kp deviation**: 0.12%

## Date

February 9, 2026

**Validation complete and successful!** ✅
