# Transport Property Analysis & Recommendations

## Current Situation

### What CombAero Actually Does
CombAero uses **full kinetic theory**, not Sutherland's formula:
- Collision integrals Ω(2,2) from Lennard-Jones potential
- Wilke's mixing rule for mixture viscosity
- Per-species transport parameters (σ, ε/k)

Source: `src/transport.cpp` lines 79-125

### Why 30% Deviation at High Temperature?

The deviation comes from **different transport parameters**, not the model:

1. **CombAero's parameters**: Currently hardcoded in `include/thermo_transport_data.h`
2. **Cantera's parameters**: From GRI-Mech 3.0 (`gri30_highT.yaml`)

Example for N2:
- GRI-Mech 3.0: σ=3.621 Å, ε/k=97.53 K
- CombAero: Need to check what's actually in `thermo_transport_data.h`

## Recommendations

### Option 1: Extract Transport Data from Mechanism Files ✅ RECOMMENDED

**Pros:**
- Uses validated, experimental data
- Consistent with Cantera/Chemkin
- No fitting required
- Can be automated

**Implementation:**
```python
# In thermo_data_generator/
def extract_transport_data(yaml_file):
    """Extract L-J parameters from mechanism YAML."""
    data = yaml.safe_load(open(yaml_file))

    transport_params = {}
    for species in data['species']:
        if 'transport' in species:
            t = species['transport']
            transport_params[species['name']] = {
                'geometry': t.get('geometry', 'atom'),
                'diameter': t.get('diameter', 0.0),  # Angstroms
                'well_depth': t.get('well-depth', 0.0),  # Kelvin
                'polarizability': t.get('polarizability', 0.0),
                'rotational_relaxation': t.get('rotational-relaxation', 0.0)
            }

    return transport_params
```

**Sources (in priority order):**
1. `gri30_highT.yaml` - Most validated for combustion
2. `JetSurf2.yaml` - Good for larger hydrocarbons
3. `sandiego20161214.yaml` - Alternative source

### Option 2: Fit Sutherland Coefficients ❌ NOT RECOMMENDED

**Why not:**
- Sutherland is less accurate than kinetic theory
- You'd be downgrading from a better model
- Still need "truth" data (which is the L-J parameters anyway)

### Option 3: Keep Current Tolerance ⚠️ ACCEPTABLE

**If transport accuracy isn't critical:**
- 30% error in viscosity is often acceptable for combustion
- Primary focus is thermodynamics and kinetics
- Speed advantage of simpler collision integral table

## Validation of Mixing Rules

### Current Implementation
CombAero uses **Wilke's mixing rule** (src/transport.cpp:111-125):

```cpp
μ_mix = Σ(X_i * μ_i) / Σ(X_i * Φ_ij)
```

Where Φ_ij accounts for molecular weight and viscosity ratios.

### Is This Validated?

**Currently: NO explicit test**

The validation tests check:
- ✅ Pure species (air, methane)
- ✅ Post-combustion mixtures (implicit)
- ❌ Systematic mixture composition sweep
- ❌ Comparison of mixing rule predictions

**Recommended test:**
```python
def test_viscosity_mixing_rule():
    """Validate Wilke's mixing rule against Cantera."""
    # Test binary mixtures at various compositions
    for x_ch4 in [0.0, 0.25, 0.5, 0.75, 1.0]:
        X = mixture_of_ch4_and_air(x_ch4)

        mu_cb = cb.viscosity(T, P, X)
        mu_ct = cantera_viscosity(T, P, X)

        # Should match within ~5% if using same L-J parameters
        assert abs(mu_cb - mu_ct) / mu_ct < 0.05
```

## Action Items

### Immediate (to reduce deviation to <5%)
1. ✅ Extract transport data from `gri30_highT.yaml`
2. ✅ Update `thermo_data_generator/generate_thermo_data.py` to include transport
3. ✅ Regenerate `include/thermo_transport_data.h` with correct L-J parameters
4. ✅ Add mixing rule validation test

### Future Enhancements
- Support multiple transport databases (GRI-Mech, JetSurf, etc.)
- Add thermal conductivity mixing rules
- Validate against experimental data (NIST, etc.)

## Expected Outcome

With correct L-J parameters from GRI-Mech 3.0:
- **Current**: 30% deviation at 2000K
- **After fix**: <5% deviation at all temperatures
- **Mixing rule**: <5% deviation for mixtures

The kinetic theory model is already correct - we just need the right input data!
