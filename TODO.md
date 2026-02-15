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

**Design Decision:** ✅ **Dataclass for combustion state**
- **C++ side:** `struct CombustionState { double phi; double T_ad; double LHV; ... };`
- **Python side:** Bound as class with readonly attributes
- **Consistent with:** Established dataclass patterns

**Properties to include (~8-10):**
- [ ] `phi` - Equivalence ratio [-]
- [ ] `T_ad` - Adiabatic flame temperature [K]
- [ ] `LHV` - Lower heating value [J/mol or J/kg]
- [ ] `AFR_mass` - Air-fuel ratio (mass basis) [-]
- [ ] `AFR_mole` - Air-fuel ratio (mole basis) [-]
- [ ] `X_products` - Product mole fractions (vector)
- [ ] `Y_products` - Product mass fractions (vector)
- [ ] `h_reactants` - Reactant enthalpy [J/mol]
- [ ] `h_products` - Product enthalpy [J/mol]
- [ ] `fuel_name` - Fuel name (string). <-- what is "fuel_name" ?
- [ ] `mixture_fraction' - Bilger Mixture fraction [-]

**Function signature:**
```python
combustion_state(fuel: str, phi: float, T_reactants: float, P: float) -> CombustionState
```

**Implementation:**
- [ ] Create `struct CombustionState` in `include/combustion.h`
- [ ] Implement `combustion_state()` function in `src/combustion.cpp`
- [ ] Add pybind11 binding in `python/combaero/_core.cpp`
- [ ] Export in `python/combaero/__init__.py`
- [ ] Write comprehensive tests in `python/tests/test_combustion_state.py`
- [ ] Update `docs/API_REFERENCE.md`
- [ ] Add unit entry to `include/units_data.h`
- [ ] Regenerate `docs/UNITS.md`

**Testing requirements:**
- Each property must match individual function call within machine precision
- Test with various fuels (CH4, C3H8, etc.)
- Test at various equivalence ratios (lean, stoichiometric, rich)
- Verify energy balance: h_reactants = h_products (for adiabatic)
- Compare with Cantera validation tests

**Status:** Not started

---

## Phase 15: CompleteState Dataclass (Optional)

**Priority:** LOW - Ultimate convenience bundle
**Difficulty:** MEDIUM - Combines all previous dataclasses

**Rationale:** For maximum convenience, provide a single dataclass that includes thermodynamic, transport, and derived properties. This is the "kitchen sink" option for users who want everything.

**Design Decision:** ✅ **Mega-dataclass combining thermo + transport**
- **C++ side:** `struct CompleteState { ThermoState thermo; TransportState transport; };`
- **Python side:** Nested dataclass or flat structure
- **Alternative:** Just document using both `thermo_state()` and `transport_state()`

**Properties to include (~20):**
- All properties from ThermoState
- All properties from TransportState
- Possibly additional derived properties

**Function signature:**
```python
complete_state(T: float, P: float, X: list[float]) -> CompleteState
```

**Implementation:**
- [ ] Decide on nested vs flat structure
- [ ] Create `struct CompleteState` in appropriate header
- [ ] Implement `complete_state()` function
- [ ] Add pybind11 binding
- [ ] Export in `__init__.py`
- [ ] Write tests
- [ ] Update documentation

**Status:** Not started (may skip in favor of separate dataclasses)

---

## Notes

- Do NOT commit this TODO.md file itself
- Each phase must be fully tested before moving to next
- Composite functions are the most valuable but require careful testing
- Machine precision testing is critical for composite functions
- **All composite functions must return dataclasses** (not tuples or individual values)
- Design with future `.to_dict()` serialization in mind
- Property names should be valid Python identifiers and good dict keys
