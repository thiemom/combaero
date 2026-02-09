# Equilibrium Validation Plan

## Current Status

### Equilibrium Functions in CombAero (`src/equilibrium.cpp`)

CombAero implements **partial equilibrium** models (not full chemical equilibrium):

1. **WGS Equilibrium** (Water-Gas Shift)
   - `wgs_equilibrium(State)` - isothermal
   - `wgs_equilibrium_adiabatic(State)` - adiabatic
   - Reaction: CO + H2O ⇌ CO2 + H2

2. **SMR + WGS Equilibrium** (Steam Methane Reforming + WGS)
   - `smr_wgs_equilibrium(State)` - isothermal
   - `smr_wgs_equilibrium_adiabatic(State)` - adiabatic
   - Reactions:
     - SMR: CH4 + H2O ⇌ CO + 3H2 (endothermic)
     - WGS: CO + H2O ⇌ CO2 + H2 (exothermic)

3. **General Reforming + WGS Equilibrium** (All hydrocarbons)
   - `reforming_equilibrium(State)` - isothermal
   - `reforming_equilibrium_adiabatic(State)` - adiabatic
   - Handles: CH4, C2H6, C3H8, iC4H10, nC5H12, nC6H14, nC7H16
   - Reactions:
     - CnHm + n·H2O ⇌ n·CO + (n + m/2)·H2
     - CO + H2O ⇌ CO2 + H2

4. **Combustion + Equilibrium** (Convenience function)
   - `combustion_equilibrium(State)` - adiabatic
   - Workflow:
     1. Complete combustion (fuel + O2 → CO2 + H2O)
     2. Reforming + WGS equilibrium on products

### Existing Tests (`python/tests/test_reforming_equilibrium.py`)

**Current coverage**: 11 tests, all passing

#### TestReformingEquilibrium (6 tests)
- ✅ Isothermal reforming with CH4
- ✅ Adiabatic reforming (temperature drop verification)
- ✅ Multiple hydrocarbons (CH4 + C2H6 + C3H8)
- ✅ C2+ hydrocarbons only (no CH4)
- ✅ Element conservation (C, H, O ratios)

#### TestSmrWgsEquilibrium (3 tests)
- ✅ Isothermal SMR+WGS
- ✅ Adiabatic SMR+WGS
- ✅ Fallback to WGS when no CH4

**Test characteristics**:
- Sanity checks (temperature trends, species presence)
- Element conservation verification
- No quantitative validation against reference implementation

## Gap Analysis

### What's Missing: Cantera Validation

Current tests verify:
- ✅ Code runs without errors
- ✅ Qualitative behavior (T decreases for endothermic, species produced)
- ✅ Element conservation

**NOT validated**:
- ❌ Equilibrium constant calculations (Kp from Gibbs free energy)
- ❌ Equilibrium composition accuracy
- ❌ Temperature convergence for adiabatic cases
- ❌ Comparison with reference implementation (Cantera)

## Recommended Additions to Cantera Validation Suite

### Phase 1: WGS Equilibrium Validation (Simplest)

**Test file**: `test_equilibrium_validation.py`

#### Tests to add:
1. **WGS isothermal equilibrium**
   - Compare equilibrium composition (CO, H2O, CO2, H2) vs Cantera
   - Test at multiple temperatures (800-1500 K)
   - Verify Kp calculation

2. **WGS adiabatic equilibrium**
   - Compare final temperature vs Cantera
   - Compare equilibrium composition
   - Verify enthalpy conservation

**Validation approach**:
```python
# CombAero
result_cb = cb.wgs_equilibrium(T, X)

# Cantera - restrict to WGS species only
species = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
wgs_species = [species[S] for S in ("CO", "H2O", "CO2", "H2", "N2", "AR")]
gas = ct.Solution(thermo="ideal-gas", species=wgs_species)
gas.TPX = T, P, X
gas.equilibrate("TP")  # Isothermal

# Compare compositions
```

### Phase 2: SMR + WGS Equilibrium Validation

#### Tests to add:
1. **SMR+WGS isothermal** (CH4 + H2O ⇌ products)
   - Compare equilibrium composition vs Cantera
   - Test at multiple temperatures (800-1200 K)
   - Verify coupled equilibrium constants

2. **SMR+WGS adiabatic**
   - Compare final temperature vs Cantera
   - Compare equilibrium composition
   - Verify enthalpy conservation

**Validation approach**:
```python
# Restrict to SMR+WGS species
smr_species = [species[S] for S in ("CH4", "H2O", "CO", "CO2", "H2", "N2", "AR")]
gas = ct.Solution(thermo="ideal-gas", species=smr_species)
gas.TPX = T, P, X
gas.equilibrate("HP")  # Adiabatic
```

### Phase 3: General Reforming Equilibrium Validation

#### Tests to add:
1. **Multi-hydrocarbon reforming**
   - CH4 + C2H6 + C3H8 mixtures
   - Compare equilibrium composition vs Cantera
   - Verify all hydrocarbons are handled correctly

2. **Combustion + equilibrium**
   - Start from unburned fuel+air
   - Compare to two-step Cantera calculation:
     1. Complete combustion (enthalpy balance)
     2. Equilibrium on products

## Expected Tolerances

Based on equilibrium validation literature:

### Composition Tolerances
- **Major species** (X > 0.01): ±0.005 absolute (0.5%)
- **Minor species** (X < 0.01): ±0.001 absolute
- **Trace species** (X < 0.001): Not validated (too sensitive)

### Temperature Tolerances
- **Adiabatic equilibrium**: ±10 K
  - Larger than complete combustion (±5 K) due to:
    - Equilibrium solver convergence differences
    - Iterative temperature solution
    - Coupled reactions

### Equilibrium Constant Tolerances
- **Kp**: ±5% relative
  - Due to Gibbs free energy polynomial differences

## Implementation Priority

### High Priority (Should add)
1. ✅ **WGS isothermal** - Simplest, validates Kp calculation
2. ✅ **WGS adiabatic** - Validates temperature solver

### Medium Priority (Nice to have)
3. ⚠️ **SMR+WGS isothermal** - More complex, coupled reactions
4. ⚠️ **SMR+WGS adiabatic** - Full validation

### Low Priority (Optional)
5. ⏸️ **General reforming** - Complex, many hydrocarbons
6. ⏸️ **Combustion + equilibrium** - Convenience function, already tested components

## Key Differences from Complete Combustion Tests

### Complete Combustion (Already validated)
- **Model**: Stoichiometric, no equilibrium
- **Method**: Enthalpy balance to find Tad
- **Products**: CO2 + H2O only
- **Validation**: Temperature within ±5 K

### Equilibrium (To be validated)
- **Model**: Chemical equilibrium (Gibbs minimization)
- **Method**: Equilibrium constants from G(T)
- **Products**: CO, H2, CO2, H2O (equilibrium mixture)
- **Validation**: Composition + temperature

## Cantera Setup for Equilibrium

```python
import cantera as ct

# Create restricted mechanism with equilibrium species
species = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}

# WGS only
wgs_species = [species[S] for S in ("CO", "H2O", "CO2", "H2", "N2", "AR")]
gas_wgs = ct.Solution(thermo="ideal-gas", species=wgs_species)

# SMR + WGS
smr_species = [species[S] for S in ("CH4", "H2O", "CO", "CO2", "H2", "N2", "AR")]
gas_smr = ct.Solution(thermo="ideal-gas", species=smr_species)

# Equilibrate
gas.TPX = T, P, X
gas.equilibrate("TP")  # Isothermal (constant T, P)
gas.equilibrate("HP")  # Adiabatic (constant H, P)
```

## Conclusion

**Recommendation**: Add WGS equilibrium validation tests (isothermal + adiabatic) to the Cantera validation suite.

**Rationale**:
1. WGS is the simplest equilibrium model
2. Validates fundamental equilibrium constant calculations
3. Provides confidence in more complex equilibrium models
4. Relatively easy to implement (similar to complete combustion tests)

**Estimated effort**: 2-3 hours
- Create `test_equilibrium_validation.py`
- Implement WGS isothermal tests (4-5 tests)
- Implement WGS adiabatic tests (2-3 tests)
- Document validation methodology
- Run and measure deviations

**Expected outcome**: Equilibrium compositions within ±0.5% (major species), temperatures within ±10 K.
