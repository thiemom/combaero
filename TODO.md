# CombAero TODO - Composite Dataclass States

## Overview
Next target: Composite Python dataclasses for thermodynamic, transport, and combustion states. These will provide convenient bundles of related properties with IDE autocomplete and type safety.

**Requirements:**
1. Every function needs comprehensive tests
2. Composite functions must match single-call results with machine precision
3. Commit after each phase if all tests and pre-commit hooks pass
4. **CODE QUALITY:**
   - **NO non-ASCII characters** in Python or C++ code (use ASCII equivalents: ° → deg, · → *, ² → 2, ³ → 3, etc.)
   - Run `./scripts/check-python-style.sh` before committing Python changes
   - Use `pre-commit run --all-files` to catch issues early
   - Target **100% pass** for every commit (all tests + CI checks green)
5. **DOCUMENTATION: Every phase must update:**
   - `docs/API_REFERENCE.md` - function signatures, parameters, units
   - `include/units_data.h` - add unit entries for new functions
   - Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`
6. **PHASE COMPLETION: At the end of each phase:**
   - **Update `TODO.md`** - mark all tasks as [x] complete
   - Add status: ✅ COMPLETE (Commit: hash)
   - Document test results, bonus achievements, and final status
   - This is the FINAL STEP of every phase before moving to the next
7. **COMPOSITE FUNCTIONS:**
   - **All composite Python functions must return dataclasses**
   - Follow the established pattern: C++ struct → pybind11 binding → Python dataclass-like behavior
   - Read-only attributes with `.def_readonly()` in pybind11
   - Nice `__repr__` for debugging
   - **Future consideration:** Add `.to_dict()` methods for serialization (design with this in mind)

---

## Design Philosophy

**Dataclass Pattern (Established in Phases 8, 10, 11):**
- C++ side: `struct PropertyBundle { double prop1; double prop2; ... };`
- Python side: Bound as class with readonly attributes (behaves like frozen dataclass)
- Benefits: IDE autocomplete, type safety, self-documenting, single function call
- Consistency: All composite functions return dataclasses

**Future Enhancement:**
- `.to_dict()` methods for serialization/JSON export
- Keep this in mind during design to ease future implementation
- Consider property naming that works well for dict keys

---

## Phase 12: ThermoState Dataclass

**Priority:** HIGH - Core thermodynamic state bundle
**Difficulty:** MEDIUM - Dataclass pattern + comprehensive property set

**Rationale:** Users frequently need multiple thermodynamic properties for a given state. Currently they must call individual functions. A `ThermoState` dataclass provides all properties in one call.

**Design Decision:** ✅ **Dataclass for thermodynamic state**
- **C++ side:** `struct ThermoState { double T; double P; double rho; double h; double s; ... };`
- **Python side:** Bound as class with readonly attributes
- **Consistent with:** AirProperties, OrificeFlowResult, AcousticProperties patterns

**Properties to include (16):**
- [x] `T` - Temperature [K] (input, echoed back)
- [x] `P` - Pressure [Pa] (input, echoed back)
- [x] `rho` - Density [kg/m³]
- [x] `cp` - Specific heat at constant pressure [J/(mol·K)]
- [x] `cv` - Specific heat at constant volume [J/(mol·K)]
- [x] `h` - Specific enthalpy [J/mol]
- [x] `s` - Specific entropy [J/(mol·K)]
- [x] `u` - Specific internal energy [J/mol]
- [x] `gamma` - Isentropic expansion coefficient [-]
- [x] `a` - Speed of sound [m/s]
- [x] `cp_mass` - Mass-specific cp [J/(kg·K)]
- [x] `cv_mass` - Mass-specific cv [J/(kg·K)]
- [x] `h_mass` - Mass-specific enthalpy [J/kg]
- [x] `s_mass` - Mass-specific entropy [J/(kg·K)]
- [x] `u_mass` - Mass-specific internal energy [J/kg]
- [x] `mw` - Molecular weight [g/mol]

**Function signature:**
```python
thermo_state(T: float, P: float, X: list[float], P_ref: float = 101325.0) -> ThermoState
```

**Implementation:**
- [x] Create `struct ThermoState` in `include/thermo.h`
- [x] Implement `thermo_state()` function in `src/thermo.cpp`
- [x] Add pybind11 binding in `python/combaero/_core.cpp`
- [x] Export in `python/combaero/__init__.py`
- [x] Write comprehensive tests in `python/tests/test_thermo_state.py`
- [x] Update `docs/API_REFERENCE.md`
- [x] Add unit entry to `include/units_data.h`
- [x] Regenerate `docs/UNITS.md`
- [x] **Design consideration:** Property names suitable for `.to_dict()` serialization

**Testing requirements:**
- [x] Each property must match individual function call within machine precision
- [x] Test with various gas compositions (air, methane, combustion products)
- [x] Test at various conditions (ambient, high T, high P)
- [x] Verify all properties are consistent (e.g., h = u + P/rho)

**Test Results:**
- ✅ 35 new tests (14 test classes)
- ✅ Machine precision verification for all 16 properties
- ✅ Physical relationships verified (gamma=cp/cv, h=u+PV, ideal gas law)
- ✅ Mass-molar conversion consistency
- ✅ Temperature range: 200-2000 K
- ✅ Pressure range: 1 kPa - 10 MPa
- ✅ All 357 tests pass (35 new + 322 existing)

**Status:** ✅ COMPLETE (Commit: 17ff696)

---

## Phase 13: TransportState Dataclass

**Priority:** HIGH - Core transport property bundle
**Difficulty:** MEDIUM - Dataclass pattern + transport properties

**Rationale:** Transport properties (viscosity, thermal conductivity, diffusivity) are frequently needed together for heat transfer and fluid dynamics calculations.

**Design Decision:** ✅ **Dataclass for transport state**
- **C++ side:** `struct TransportState { double mu; double k; double nu; double alpha; double Pr; };`
- **Python side:** Bound as class with readonly attributes
- **Consistent with:** Established dataclass patterns

**Properties to include (9):**
- [x] `T` - Temperature [K] (input, echoed back)
- [x] `P` - Pressure [Pa] (input, echoed back)
- [x] `rho` - Density [kg/m³]
- [x] `mu` - Dynamic viscosity [Pa·s]
- [x] `k` - Thermal conductivity [W/(m·K)]
- [x] `nu` - Kinematic viscosity [m²/s]
- [x] `alpha` - Thermal diffusivity [m²/s]
- [x] `Pr` - Prandtl number [-]
- [x] `cp` - Specific heat at constant pressure [J/(kg·K)]

**Function signature:**
```python
transport_state(T: float, P: float, X: list[float]) -> TransportState
```

**Implementation:**
- [x] Create `struct TransportState` in `include/transport.h`
- [x] Implement `transport_state()` function in `src/transport.cpp`
- [x] Add pybind11 binding in `python/combaero/_core.cpp`
- [x] Export in `python/combaero/__init__.py`
- [x] Write comprehensive tests in `python/tests/test_transport_state.py`
- [x] Update `docs/API_REFERENCE.md`
- [x] Add unit entry to `include/units_data.h`
- [x] Regenerate `docs/UNITS.md`

**Testing requirements:**
- [x] Each property must match individual function call within machine precision
- [x] Test with various gas compositions
- [x] Test temperature dependence (Sutherland's law validation)
- [x] Verify Pr = mu * cp / k

**Test Results:**
- ✅ 31 new tests (13 test classes)
- ✅ Machine precision verification using individual function calls as reference
- ✅ Physical relationships verified (nu=mu/rho, alpha=k/(rho*cp), Pr=mu*cp/k)
- ✅ Various compositions tested (air, N2, CH4, combustion products)
- ✅ Temperature range: 200-2000 K
- ✅ Pressure range: 1 kPa - 10 MPa
- ✅ All 388 tests pass (31 new + 357 existing)

**Status:** ✅ COMPLETE (Commit: 1d9769f)

---

## Phase 14: CombustionState Dataclass

**Priority:** MEDIUM - Combustion analysis bundle
**Difficulty:** MEDIUM-HIGH - Dataclass pattern + combustion calculations

**Rationale:** Combustion calculations involve many related properties (equivalence ratio, adiabatic temperature, product composition). A dataclass bundles these for convenient analysis.

**Design Decision:** ✅ **Nested dataclass with ThermoState objects**
- **C++ side:** `struct CombustionState { ThermoState reactants; ThermoState products; double phi; ... };`
- **Python side:** Bound as class with readonly attributes, nested ThermoState access
- **Consistent with:** Established dataclass patterns + composition over duplication
- **Key insights:**
  - Fuel/oxidizer are composition vectors (not strings) to support fuel mixtures
  - Reactants and products are full thermodynamic states - reuse ThermoState!
  - Avoids duplication: instead of storing T, P, X, h separately, nest ThermoState objects

**Properties (~4 core + 2×25 from nested CompleteStates = ~54 total):**
- [ ] `phi` - Equivalence ratio [-] (input, echoed back)
- [ ] `fuel_name` - Optional fuel label (string, default: "")
- [ ] `reactants` - **CompleteState** at (T_reactants, P, X_reactants)
  - Includes all 16 ThermoState properties: T, P, rho, cp, cv, h, s, u, gamma, a, cp_mass, cv_mass, h_mass, s_mass, u_mass, mw
  - Includes all 9 TransportState properties: T, P, rho, mu, k, nu, alpha, Pr, cp
- [ ] `products` - **CompleteState** at (T_ad, P, X_products)
  - Includes all 16 ThermoState properties: T, P, rho, cp, cv, h, s, u, gamma, a, cp_mass, cv_mass, h_mass, s_mass, u_mass, mw
  - Includes all 9 TransportState properties: T, P, rho, mu, k, nu, alpha, Pr, cp
- [ ] `mixture_fraction` - Bilger mixture fraction [-]
- [ ] `fuel_burn_fraction` - Fraction of fuel burned [0-1] (useful for O2-limited cases)

**Function signatures:**

**Variant 1: From equivalence ratio (for calculations)**
```python
combustion_state(
    X_fuel: list[float],      # Fuel composition (can be pure or mixture)
    X_ox: list[float],        # Oxidizer composition
    phi: float,               # Equivalence ratio [-] (INPUT)
    T_reactants: float,       # Reactant temperature [K]
    P: float,                 # Pressure [Pa]
    fuel_name: str = ""       # Optional label for fuel (e.g., "CH4", "Natural Gas")
) -> CombustionState
```

**Variant 2: From streams (for lab measurements where flows are measured)**
```python
combustion_state_from_streams(
    fuel_stream: Stream,      # Fuel stream with mdot, T, X
    oxidizer_stream: Stream,  # Oxidizer stream with mdot, T, X
    fuel_name: str = ""       # Optional label for fuel
) -> CombustionState
# In this variant, phi is COMPUTED from mass flow rates (output, not input)
```

**Usage examples:**

**Example 1: From phi (typical for calculations)**
```python
X_CH4 = [0, 0, 0, 0, 0, 1.0, 0, ...]  # Pure methane
X_air = cb.standard_dry_air_composition()
state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)

# Combustion properties (phi was INPUT)
print(state.phi)                         # 1.0 (echoed back)
print(state.mixture_fraction)            # 0.055
```

**Example 2: From streams (typical for lab measurements)**
```python
fuel = cb.Stream().set_T(300).set_X(X_CH4).set_mdot(0.01)  # 0.01 kg/s
air = cb.Stream().set_T(298).set_X(X_air).set_mdot(0.17)   # 0.17 kg/s (measured)
state = cb.combustion_state_from_streams(fuel, air, fuel_name="CH4")

# Combustion properties (phi was COMPUTED from mdot)
print(state.phi)                         # 1.0 (computed from flow rates)
print(state.mixture_fraction)            # 0.055

# Reactant thermo properties (via CompleteState)
print(state.reactants.thermo.T)          # 300 K
print(state.reactants.thermo.h)          # -103.6 J/mol
print(state.reactants.thermo.rho)        # 1.177 kg/m³

# Reactant transport properties (via CompleteState)
print(state.reactants.transport.mu)      # 1.68e-5 Pa·s
print(state.reactants.transport.Pr)      # 0.781

# Product thermo properties (via CompleteState)
print(state.products.thermo.T)           # 2328 K (T_ad)
print(state.products.thermo.h)           # -103.6 J/mol (energy balance)
print(state.products.thermo.gamma)       # 1.28 (lower at high T)

# Product transport properties (via CompleteState)
print(state.products.transport.mu)       # 6.8e-5 Pa·s (higher at high T)
print(state.products.transport.Pr)       # 0.73
```

**Implementation:**
- [ ] **Prerequisite:** Phase 15 (CompleteState) must be implemented first
- [ ] Create `struct CombustionState` in `include/combustion.h` (with nested CompleteState)
- [ ] Implement `combustion_state()` function in `src/combustion.cpp` (from phi)
  - Compute X_reactants from set_equivalence_ratio_mole()
  - Compute X_products from complete_combustion_to_CO2_H2O()
  - Call complete_state() for reactants and products
  - Compute mixture_fraction from bilger_mixture_fraction_from_moles()
- [ ] Implement `combustion_state_from_streams()` function in `src/combustion.cpp`
  - Mix streams to get X_reactants and T_reactants
  - Compute phi from equivalence_ratio_mole()
  - Compute X_products from complete_combustion_to_CO2_H2O()
  - Call complete_state() for reactants and products
  - Compute mixture_fraction from bilger_mixture_fraction_from_moles()
- [ ] Add pybind11 bindings for both functions in `python/combaero/_core.cpp`
- [ ] Export both functions in `python/combaero/__init__.py`
- [ ] Write comprehensive tests in `python/tests/test_combustion_state.py`
  - Test both phi-based and stream-based variants
  - Verify phi matches between variants for same conditions
- [ ] Update `docs/API_REFERENCE.md`
- [ ] Add unit entries to `include/units_data.h`
- [ ] Regenerate `docs/UNITS.md`

**Testing requirements:**
- Verify nested CompleteState properties match individual complete_state() calls
- Verify energy balance: state.reactants.thermo.h == state.products.thermo.h (for adiabatic)
- Test with various fuels (CH4, C3H8, fuel mixtures)
- Test at various equivalence ratios (lean, stoichiometric, rich)
- **Test both variants produce same results:**
  - combustion_state(X_fuel, X_ox, phi=1.0, ...)
  - combustion_state_from_streams(fuel, oxidizer) with mdot set for phi=1.0
- Verify phi matches equivalence_ratio_mole() calculation (both variants)
- Verify mixture_fraction matches bilger_mixture_fraction_from_moles()
- Compare T_ad with Cantera validation tests
- Verify transport properties are reasonable (mu increases with T, etc.)
- **Stream variant specific:** Verify phi is correctly computed from mass flow rates

**Status:** Not started (implement Phase 15 first to use in Phase 14)

---

## Phase 15: CompleteState Dataclass

**Priority:** HIGH - Building block for CombustionState
**Difficulty:** LOW - Simple composition of existing dataclasses

**Rationale:** Combine thermodynamic and transport properties in one dataclass. This will be used as a building block for CombustionState (reactants and products) and is also useful standalone.

**Design Decision:** ✅ **Nested dataclass combining thermo + transport**
- **C++ side:** `struct CompleteState { ThermoState thermo; TransportState transport; };`
- **Python side:** Nested structure with `state.thermo.X` and `state.transport.X` access
- **Reusable:** Will be used in Phase 14 for reactants and products

**Properties (25 total = 16 thermo + 9 transport):**
- [x] `thermo` - **ThermoState** with all 16 thermodynamic properties
  - T, P, rho, cp, cv, h, s, u, gamma, a, cp_mass, cv_mass, h_mass, s_mass, u_mass, mw
- [x] `transport` - **TransportState** with all 9 transport properties
  - T, P, rho, mu, k, nu, alpha, Pr, cp

**Function signature:**
```python
complete_state(T: float, P: float, X: list[float], P_ref: float = 101325.0) -> CompleteState
```

**Usage example:**
```python
X = cb.standard_dry_air_composition()
state = cb.complete_state(T=300, P=101325, X=X)

# Thermodynamic properties
print(state.thermo.h)        # -103.6 J/mol
print(state.thermo.gamma)    # 1.400

# Transport properties
print(state.transport.mu)    # 1.68e-5 Pa·s
print(state.transport.Pr)    # 0.781
```

**Implementation:**
- [x] Create `struct CompleteState` in `include/thermo.h`
- [x] Implement `complete_state()` function (calls thermo_state + transport_state)
- [x] Add pybind11 binding in `python/combaero/_core.cpp`
- [x] Export in `python/combaero/__init__.py`
- [x] Write comprehensive tests in `python/tests/test_complete_state.py`
- [x] Update `docs/API_REFERENCE.md`
- [x] Add unit entry to `include/units_data.h`
- [x] Regenerate `docs/UNITS.md`

**Testing requirements:**
- [x] Verify state.thermo matches thermo_state() call
- [x] Verify state.transport matches transport_state() call
- [x] Test with various compositions and conditions
- [x] Verify all nested properties are accessible

**Test Results:**
- ✅ 23 new tests (8 test classes)
- ✅ Nested states match individual function calls at machine precision
- ✅ Consistency verified between thermo and transport (T, P, rho)
- ✅ Various compositions tested (air, N2, CH4)
- ✅ Temperature range: 200-1500 K
- ✅ Pressure range: 1 kPa - 2 MPa
- ✅ All 411 tests pass (23 new + 388 existing)

**Status:** ✅ COMPLETE (Commit: cc9a10c)

---

## Notes

- Do NOT commit this TODO.md file itself
- Each phase must be fully tested before moving to next
- Composite functions are the most valuable but require careful testing
- Machine precision testing is critical for composite functions
- **All composite functions must return dataclasses** (not tuples or individual values)
- Design with future `.to_dict()` serialization in mind
- Property names should be valid Python identifiers and good dict keys
