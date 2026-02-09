# Validation Methodology

## Complete Combustion vs Full Equilibrium

### CombAero's Complete Combustion Model

CombAero's `complete_combustion()` function implements **stoichiometric complete combustion** that produces only CO2 and H2O:

```
CxHyOz + (x + y/4 - z/2) O2 → x CO2 + (y/2) H2O
```

**Key characteristics:**
- No equilibrium species (CO, H2, OH, radicals, etc.)
- Simple stoichiometric balance
- Fast calculation (no equilibrium solver)
- Suitable for preliminary design and screening

### Why Restricted Species in Cantera?

If we use full GRI-Mech 3.0 (53 species) and call `equilibrate("HP")`, Cantera will produce:
- CO, H2 (from dissociation and WGS equilibrium)
- OH, H, O (radicals)
- Other intermediate species

This is **full chemical equilibrium**, not complete combustion.

### Solution: Restricted Species Gas Phase

We create a Cantera gas phase with **only complete combustion species**:

```python
species = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
complete_species = [
    species[S] for S in ("N2", "O2", "AR", "CO2", "H2O", "CH4", "C3H8", "H2", "CO")
]
gas = ct.Solution(thermo="ideal-gas", species=complete_species)
```

Now when we call `gas.equilibrate("HP")`, Cantera can only produce CO2 and H2O because:
- CO and H2 are not available as products (only as fuels)
- No intermediate species exist in the mechanism
- Result matches CombAero's complete combustion model

## Validation Test Types

### 1. Complete Combustion Tests
**Fixture**: `complete_combustion_gas` (restricted species)
**Validates**: Adiabatic flame temperature, product composition
**Approach**: Restricted species equilibrium = complete combustion

### 2. Mixing Tests
**Fixture**: `gri30_gas` (full mechanism) or state-based
**Validates**: Enthalpy conservation, mixture properties
**Approach**: Direct property comparison, no combustion

### 3. Transport Tests
**Fixture**: `gri30_gas` (full mechanism)
**Validates**: Viscosity, thermal conductivity, Prandtl
**Approach**: Direct transport property comparison

## When to Use Each Approach

### Use Restricted Species When:
- Validating `complete_combustion()` function
- Testing adiabatic flame temperature for complete combustion
- Comparing product composition (CO2, H2O only)

### Use Full Mechanism When:
- Validating equilibrium functions (WGS, reforming)
- Testing transport properties
- Validating mixing (no combustion involved)
- Testing density, Cp, enthalpy (state properties)

## Example: Correct vs Incorrect

### ❌ Incorrect (Full Equilibrium)
```python
gas = ct.Solution("gri30.yaml")  # 53 species
gas.TP = 300, 101325
gas.set_equivalence_ratio(1.0, "CH4", "O2:1, N2:3.76")
gas.equilibrate("HP")  # Produces CO, H2, OH, etc.
# Result won't match CombAero's complete_combustion()
```

### ✅ Correct (Restricted Species)
```python
species = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
complete_species = [species[S] for S in ("N2", "O2", "AR", "CO2", "H2O", "CH4")]
gas = ct.Solution(thermo="ideal-gas", species=complete_species)
gas.TP = 300, 101325
gas.set_equivalence_ratio(1.0, "CH4", "O2:1, N2:3.76")
gas.equilibrate("HP")  # Can only produce CO2 and H2O
# Result matches CombAero's complete_combustion()
```

## Thermodynamic Consistency

Both approaches are thermodynamically consistent:
- **Complete combustion**: Assumes instantaneous, complete reaction to CO2/H2O
- **Full equilibrium**: Accounts for dissociation and reverse reactions

The difference is in the **chemical model**, not the thermodynamics. Our validation ensures CombAero's complete combustion model is implemented correctly for its intended use case.

## Future: Equilibrium Validation

When validating CombAero's equilibrium functions (WGS, reforming), we will:
1. Use full GRI-Mech 3.0 mechanism
2. Compare equilibrium compositions
3. Use appropriate tolerances for equilibrium solver differences

This is separate from complete combustion validation.
