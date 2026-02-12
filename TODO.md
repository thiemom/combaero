# Missing Utility Functions Implementation Plan

## Overview
Add helper functions identified from example analysis to simplify user code and reduce repetition.

**Requirements:**
1. Every function needs comprehensive tests
2. Composite functions must match single-call results with machine precision
3. Split into phases: easy → difficult, sorted by priority
4. Commit after each phase if all tests and pre-commit hooks pass
5. Do NOT commit this TODO.md file itself
6. **DOCUMENTATION: Every phase must update:**
   - `docs/API_REFERENCE.md` - function signatures, parameters, units
   - `include/units_data.h` - add unit entries for new functions
   - Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`

---

## Phase 1: Simple Geometry Helpers ✅

**Priority:** HIGH - Simple, frequently used, no dependencies
**Difficulty:** EASY - Pure geometric calculations

**Functions to implement (3):**
- [x] `pipe_area(D)` - Circular pipe cross-sectional area [m²]
- [x] `annular_area(D_outer, D_inner)` - Annular cross-sectional area [m²]
- [x] `pipe_volume(D, L)` - Cylindrical pipe volume [m³]

**Implementation:**
- [x] Add functions to `src/geometry.cpp`
- [x] Add declarations to `include/geometry.h`
- [x] Add Python bindings to `python/combaero/_core.cpp`
- [x] Write unit tests in `python/tests/test_geometry.py` (16 tests)
- [x] Verify formulas: A = π(D/2)², A_annular = π((D_o/2)² - (D_i/2)²), V = πD²L/4
- [x] Update `docs/API_REFERENCE.md`
- [x] **Add unit entries to `include/units_data.h`**
- [x] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [x] Run all tests: `pytest python/tests/` - **137/137 PASSED**
- [x] Commit: "Add simple geometry helper functions (pipe_area, annular_area, pipe_volume)"

**Testing requirements:**
- Verify against manual calculations with known values
- Test edge cases (D=0 should error, D_outer < D_inner should error)
- Verify units are correct [m²] and [m³]

**Status:** ✅ COMPLETED

---

## Phase 2: Mass-Specific Thermodynamic Properties ⏳

**Priority:** HIGH - Frequently needed, simple conversions
**Difficulty:** EASY - Simple wrappers around existing functions

**Rationale:** Maintain symmetric interface - if we have molar-basis functions (cp, h, cv, s),
we should have mass-basis equivalents for all of them to avoid user confusion.

**Functions to implement (4):**
- [ ] `cp_mass(T, X)` - Mass-specific heat capacity at constant pressure [J/(kg·K)]
- [ ] `cv_mass(T, X)` - Mass-specific heat capacity at constant volume [J/(kg·K)]
- [ ] `h_mass(T, X)` - Mass-specific enthalpy [J/kg]
- [ ] `s_mass(T, P, X, P_ref)` - Mass-specific entropy [J/(kg·K)]

**Implementation:**
- [ ] Add functions to `src/thermo.cpp`
- [ ] Add declarations to `include/thermo.h`
- [ ] Add Python bindings to `python/combaero/_core.cpp`
- [ ] Export in `python/combaero/__init__.py`
- [ ] Write unit tests in `python/tests/test_thermo_transport.py`
- [ ] **CRITICAL TESTS:**
  - Verify cp_mass(T, X) == cp(T, X) / mwmix(X) * 1000
  - Verify cv_mass(T, X) == cv(T, X) / mwmix(X) * 1000
  - Verify h_mass(T, X) == h(T, X) / mwmix(X) * 1000
  - Verify s_mass(T, P, X) == s(T, P, X) / mwmix(X) * 1000
  - All must match within machine precision
- [ ] Test with multiple compositions (air, fuel, products)
- [ ] Update `docs/API_REFERENCE.md`
- [ ] **Add unit entries to `include/units_data.h` (4 entries)**
- [ ] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add mass-specific thermodynamic property functions"

**Testing requirements:**
- Each function must match molar / mwmix * 1000 within machine precision
- Test with air, fuel mixtures, combustion products
- Verify units [J/(kg·K)] for cp/cv/s, [J/kg] for h
- Test consistency: gamma_mass = cp_mass / cv_mass should equal gamma_molar

**Status:** Not started

---

## Phase 3: Velocity from Mass Flow ⏳

**Priority:** HIGH - Common conversion in flow calculations
**Difficulty:** EASY - Simple formula with geometry

**Functions to implement (1):**
- [ ] `velocity_from_mdot(mdot, rho, D)` - Velocity from mass flow [m/s]

**Implementation:**
- [ ] Add function to `src/incompressible.cpp` (or new `src/flow_utils.cpp`)
- [ ] Add declaration to `include/incompressible.h`
- [ ] Add Python binding to `python/combaero/_core.cpp`
- [ ] Write unit tests in `python/tests/test_flow_utils.py`
- [ ] **CRITICAL TEST:** Verify v = mdot / (rho * π(D/2)²)
- [ ] Test round-trip: mdot → v → mdot should match
- [ ] Update `docs/API_REFERENCE.md`
- [ ] **Add unit entry to `include/units_data.h`**
- [ ] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add velocity from mass flow function (velocity_from_mdot)"

**Testing requirements:**
- Must match manual calculation: v = mdot / (rho * A) where A = π(D/2)²
- Test round-trip consistency
- Verify units [m/s]

**Status:** Not started

---

## Phase 4: Pressure Drop in Pipe (Composite) ⏳

**Priority:** HIGH - Very common, combines multiple steps
**Difficulty:** MEDIUM - Composite function, multiple dependencies

**Functions to implement (1):**
- [ ] `pressure_drop_pipe(T, P, X, v, D, L, roughness=0.0, correlation='haaland')`
  - Returns: tuple (dP [Pa], Re [-], f [-])

**Implementation:**
- [ ] Add function to `src/incompressible.cpp` or new `src/pipe_flow.cpp`
- [ ] Add declaration to header
- [ ] Add Python binding with tuple return
- [ ] Write comprehensive unit tests in `python/tests/test_pipe_flow.py`
- [ ] **CRITICAL TESTS:**
  - Verify dP matches manual calculation: f * (L/D) * (ρv²/2)
  - Verify Re matches: ρvD/μ
  - Verify f matches friction correlation
  - Test all correlations: 'haaland', 'serghides', 'colebrook', 'petukhov'
  - Test smooth vs rough pipes
- [ ] Update `docs/API_REFERENCE.md`
- [ ] **Add unit entry to `include/units_data.h`**
- [ ] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add composite pressure drop function (pressure_drop_pipe)"

**Testing requirements:**
- Decompose and verify each step matches individual function calls
- Test with air at different conditions
- Test with different roughness values
- Verify all correlation options work
- Machine precision match with manual multi-step calculation

**Status:** Not started

---

## Phase 5: Heat Transfer Coefficient in Pipe (Composite) ⏳

**Priority:** HIGH - Very common heat transfer calculation
**Difficulty:** MEDIUM - Composite function, multiple dependencies

**Functions to implement (1):**
- [ ] `htc_pipe(T, P, X, v, D, correlation='dittus_boelter', heating=True, mu_ratio=1.0, roughness=0.0)`
  - Returns: tuple (h [W/(m²·K)], Nu [-], Re [-])

**Implementation:**
- [ ] Add function to `src/heat_transfer.cpp`
- [ ] Add declaration to `include/heat_transfer.h`
- [ ] Add Python binding with tuple return
- [ ] Write comprehensive unit tests in `python/tests/test_heat_transfer.py`
- [ ] **CRITICAL TESTS:**
  - Verify h matches: Nu * k / D
  - Verify Nu matches correlation (Dittus-Boelter, Gnielinski, Sieder-Tate)
  - Verify Re matches: ρvD/μ
  - Test all correlations
  - Test heating vs cooling
  - Test with friction factor for Gnielinski
- [ ] Update `docs/API_REFERENCE.md`
- [ ] **Add unit entry to `include/units_data.h`**
- [ ] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add composite heat transfer coefficient function (htc_pipe)"

**Testing requirements:**
- Decompose and verify each step matches individual function calls
- Test all correlation options
- Test heating vs cooling for Dittus-Boelter
- Test viscosity correction for Sieder-Tate
- Machine precision match with manual multi-step calculation

**Status:** Not started

---

## Phase 6: Pipe Roughness Database ⏳

**Priority:** MEDIUM - Convenience for users
**Difficulty:** EASY - Just data storage

**Functions to implement (2):**
- [ ] `pipe_roughness(material)` - Returns absolute roughness ε [m]
- [ ] `standard_pipe_roughness()` - Returns dict of all materials

**Implementation:**
- [ ] Add to `src/utils.cpp` or new file
- [ ] Create static map/dict of standard roughness values
- [ ] Materials: 'smooth', 'drawn_tubing', 'commercial_steel', 'galvanized_iron', 'cast_iron', 'concrete', 'riveted_steel'
- [ ] Add Python bindings
- [ ] Write unit tests
- [ ] Update `docs/API_REFERENCE.md`
- [ ] **Add unit entries to `include/units_data.h`**
- [ ] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add pipe roughness database functions"

**Testing requirements:**
- Verify values match engineering handbooks
- Test case-insensitive lookup
- Test error handling for unknown materials

**Status:** Not started

---

## Phase 7: LMTD Convenience Wrappers ⏳

**Priority:** MEDIUM - Convenience, saves users from calculating dT
**Difficulty:** EASY - Simple wrappers around existing lmtd()

**Functions to implement (2):**
- [ ] `lmtd_counterflow(T_hot_in, T_hot_out, T_cold_in, T_cold_out)` - Returns LMTD [K]
- [ ] `lmtd_parallelflow(T_hot_in, T_hot_out, T_cold_in, T_cold_out)` - Returns LMTD [K]

**Implementation:**
- [ ] Add to `src/heat_transfer.cpp`
- [ ] Add declarations to `include/heat_transfer.h`
- [ ] Add Python bindings
- [ ] Write unit tests in `python/tests/test_heat_transfer.py`
- [ ] **CRITICAL TESTS:**
  - Verify counterflow: lmtd(T_hot_in - T_cold_out, T_hot_out - T_cold_in)
  - Verify parallelflow: lmtd(T_hot_in - T_cold_in, T_hot_out - T_cold_out)
  - Must match existing lmtd() function exactly
- [ ] Update `docs/API_REFERENCE.md`
- [ ] **Add unit entries to `include/units_data.h`**
- [ ] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add LMTD convenience wrapper functions"

**Testing requirements:**
- Must match manual lmtd() calls exactly
- Test with realistic heat exchanger temperatures
- Verify counterflow > parallelflow for same conditions

**Status:** Not started

---

## Phase 8: Air Properties Bundle (Optional) ⏳

**Priority:** LOW - Nice to have, but multiple calls aren't too bad
**Difficulty:** MEDIUM - Returns struct/dict, needs careful design

**Functions to implement (1):**
- [ ] `air_properties(T, P, humidity=0.0)` - Returns struct/dict with all properties

**Implementation:**
- [ ] Design return type (struct in C++, dict or named tuple in Python)
- [ ] Add to `src/thermo.cpp` or `src/air_utils.cpp`
- [ ] Properties: rho, mu, k, cp, cv, Pr, nu, alpha, gamma, a
- [ ] Add Python binding
- [ ] Write unit tests
- [ ] **CRITICAL TESTS:** Each property must match individual function call
- [ ] Update `docs/API_REFERENCE.md`
- [ ] **Add unit entry to `include/units_data.h`**
- [ ] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add air properties bundle function"

**Testing requirements:**
- Each property must match individual function call exactly
- Test with dry air and humid air
- Test at different temperatures and pressures

**Status:** Not started

---

## Phase 9: Update Examples to Use Convenience Functions ⏳

**Priority:** HIGH - Demonstrate value of new functions, improve code quality
**Difficulty:** EASY - Refactoring existing working code

**Goal:** Update all 4 new example files to use the new convenience functions where applicable.

**Files to update:**
- [ ] `python/examples/friction_and_pressure_drop.py`
- [ ] `python/examples/heat_transfer_pipe_flow.py`
- [ ] `python/examples/reactor_design.py`
- [ ] `python/examples/combustor_heat_transfer.py`

**Changes to make:**
- [ ] Replace manual Reynolds number calculations with existing `ca.reynolds()` or composite functions
- [ ] Replace manual area calculations with `ca.pipe_area()` and `ca.annular_area()`
- [ ] Replace manual volume calculations with `ca.pipe_volume()`
- [ ] Replace multi-step pressure drop with `ca.pressure_drop_pipe()`
- [ ] Replace multi-step heat transfer with `ca.htc_pipe()`
- [ ] Replace manual cp_mass conversions with `ca.cp_mass()`
- [ ] Replace manual velocity calculations with `ca.velocity_from_mdot()`
- [ ] Use `ca.pipe_roughness()` for standard materials
- [ ] Use `ca.lmtd_counterflow()` / `ca.lmtd_parallelflow()` where applicable

**Implementation:**
- [ ] Update each example file one at a time
- [ ] **CRITICAL:** Verify output remains identical (or very close) after refactoring
- [ ] Run each example to ensure it still works: `python python/examples/<file>.py`
- [ ] Ensure examples are cleaner and more readable
- [ ] Add comments showing the convenience of new functions
- [ ] Run all tests to ensure nothing broke: `pytest python/tests/`
- [ ] Commit: "Refactor examples to use new convenience functions"

**Testing requirements:**
- Each example must still run successfully
- Output should be identical (or negligibly different due to rounding)
- Code should be noticeably cleaner and shorter
- Examples should better demonstrate CombAero's ease of use

**Status:** Not started

**Note:** This phase should only be done AFTER all previous phases are complete and working.

---

## Phase 10: Orifice Flow Utilities ⏳

**Priority:** MEDIUM - Useful for flow metering and injection applications
**Difficulty:** EASY - Simple geometric and flow calculations

**Rationale:** The orifice module is comprehensive but lacks simple helper functions for common calculations. No examples exist demonstrating the Cd correlation system.

**Functions to implement (4):**
- [ ] `orifice_velocity(mdot, rho, d)` - Velocity through orifice throat [m/s]
- [ ] `orifice_area_from_beta(D, beta)` - Orifice area from beta ratio [m²]
- [ ] `beta_from_diameters(d, D)` - Beta ratio from diameters [-] (standalone)
- [ ] `orifice_Re_d(mdot, d, mu)` - Orifice Reynolds number from mass flow [-]

**Implementation:**
- [ ] Add functions to `src/orifice.cpp`
- [ ] Add declarations to `include/orifice.h`
- [ ] Add Python bindings to `python/combaero/_core.cpp`
- [ ] Export in `python/combaero/__init__.py`
- [ ] Write unit tests in `python/tests/test_orifice.py`
- [ ] **CRITICAL TESTS:**
  - Verify orifice_velocity = mdot / (rho * π(d/2)²)
  - Verify orifice_area_from_beta = π(D*beta/2)²
  - Verify beta_from_diameters = d/D
  - Verify orifice_Re_d = 4*mdot / (π*d*mu)
  - All must match within machine precision
- [ ] Update `docs/API_REFERENCE.md`
- [ ] **Add unit entries to `include/units_data.h` (4 entries)**
- [ ] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add orifice flow utility functions"

**Testing requirements:**
- Each function must match manual calculation within machine precision
- Test with realistic orifice geometries (beta = 0.2 to 0.7)
- Verify units: [m/s], [m²], [-], [-]
- Round-trip tests where applicable

**Example to create:**
- [ ] `python/examples/orifice_flow_metering.py` - Demonstrate Cd correlations, beta effects, Re effects, flow measurement

**Status:** Not started

---

## Phase 11: Acoustics Utilities ⏳

**Priority:** MEDIUM - Useful for combustor design and instability analysis
**Difficulty:** EASY - Simple frequency and wavelength calculations

**Rationale:** The acoustics module is comprehensive but lacks simple helper functions. No examples exist demonstrating acoustic mode analysis for combustors.

**Functions to implement (5):**
- [ ] `wavelength(f, c)` - Acoustic wavelength [m]
- [ ] `frequency_from_wavelength(lambda, c)` - Frequency from wavelength [Hz]
- [ ] `acoustic_impedance(rho, c)` - Characteristic acoustic impedance [Pa·s/m]
- [ ] `sound_pressure_level(p_rms, p_ref)` - SPL in dB [dB]
- [ ] `particle_velocity(p, rho, c)` - Particle velocity from pressure [m/s]

**Implementation:**
- [ ] Add functions to `src/acoustics.cpp`
- [ ] Add declarations to `include/acoustics.h`
- [ ] Add Python bindings to `python/combaero/_core.cpp`
- [ ] Export in `python/combaero/__init__.py`
- [ ] Write unit tests in `python/tests/test_acoustics.py`
- [ ] **CRITICAL TESTS:**
  - Verify wavelength = c / f
  - Verify frequency_from_wavelength = c / lambda
  - Verify acoustic_impedance = rho * c
  - Verify SPL = 20*log10(p_rms/p_ref)
  - Verify particle_velocity = p / (rho*c)
  - All must match within machine precision
- [ ] Update `docs/API_REFERENCE.md`
- [ ] **Add unit entries to `include/units_data.h` (5 entries)**
- [ ] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add acoustics utility functions"

**Testing requirements:**
- Each function must match manual calculation within machine precision
- Test with realistic combustor conditions (T=1500K, P=15 bar)
- Verify units: [m], [Hz], [Pa·s/m], [dB], [m/s]
- Test consistency: wavelength * frequency = c

**Example to create:**
- [ ] `python/examples/combustor_acoustics.py` - Demonstrate mode analysis, Helmholtz resonators, quarter-wave dampers, instability screening

**Status:** Not started

---

## Completion Checklist

- [x] All Phase 1 tests passing (geometry helpers) ✅
- [ ] All Phase 2 tests passing (mass-specific thermo properties)
- [ ] All Phase 3 tests passing (velocity_from_mdot)
- [ ] All Phase 4 tests passing (pressure_drop_pipe)
- [ ] All Phase 5 tests passing (htc_pipe)
- [ ] All Phase 6 tests passing (roughness database)
- [ ] All Phase 7 tests passing (LMTD wrappers)
- [ ] All Phase 8 tests passing (air_properties - optional)
- [ ] **Phase 9: Examples updated to use new convenience functions**
- [ ] All examples run successfully with new functions
- [ ] Code is cleaner and demonstrates ease of use
- [ ] **Phase 10: Orifice utilities and example**
- [ ] All orifice utility tests passing
- [ ] Orifice flow metering example created and working
- [ ] **Phase 11: Acoustics utilities and example**
- [ ] All acoustics utility tests passing
- [ ] Combustor acoustics example created and working
- [ ] Documentation updated for all new functions
- [ ] All commits pushed to main

---

## Testing Strategy

### For Simple Functions:
1. Test against manual calculations with known values
2. Test edge cases and error conditions
3. Verify units are correct

### For Composite Functions:
1. **Decompose and verify:** Each step must match individual function calls
2. **Machine precision:** Results must match within floating-point precision (~1e-15 relative error)
3. **Round-trip tests:** Where applicable (e.g., mdot → v → mdot)
4. **Multiple scenarios:** Test with different fluids, conditions, correlations

### General:
- All tests must pass before committing
- Pre-commit hooks must pass (black, ruff, etc.)
- Document all units in docstrings
- Follow existing code patterns in CombAero

---

## Notes

- Do NOT commit this TODO.md file itself
- Each phase must be fully tested before moving to next
- Composite functions are the most valuable but require careful testing
- Start with simple geometry helpers to build confidence
- Machine precision testing is critical for composite functions
- Update examples in a separate task after all functions are implemented

---

## Estimated Effort

- **Phase 1:** 1-2 hours (simple geometry) ✅ DONE
- **Phase 2:** 1-1.5 hours (4 simple wrappers, symmetric interface)
- **Phase 3:** 30 min (simple formula)
- **Phase 4:** 2-3 hours (composite, multiple correlations)
- **Phase 5:** 2-3 hours (composite, multiple correlations)
- **Phase 6:** 1 hour (data entry)
- **Phase 7:** 1 hour (simple wrappers)
- **Phase 8:** 2 hours (struct design, optional)
- **Phase 9:** 2-3 hours (refactor 4 examples)
- **Phase 10:** 1.5 hours (4 orifice utilities + example)
- **Phase 11:** 1.5 hours (5 acoustics utilities + example)

**Total Core (1-7,9):** ~11-15 hours
**Total Extended (1-7,9-11):** ~14-18 hours
**Total Complete (1-11):** ~16-20 hours

---

## Priority Order (Recommended Implementation)

### **Core Phases (High Priority):**
1. **Phase 1** - Geometry helpers ✅ DONE
2. **Phase 2** - Mass-specific thermo (symmetric interface)
3. **Phase 3** - velocity_from_mdot (simple, high value)
4. **Phase 4** - pressure_drop_pipe (high value composite)
5. **Phase 5** - htc_pipe (high value composite)
6. **Phase 6** - Roughness database (convenience)
7. **Phase 7** - LMTD wrappers (convenience)
9. **Phase 9** - Update examples (demonstrate value)

### **Extended Phases (Medium Priority):**
10. **Phase 10** - Orifice utilities + example (flow metering)
11. **Phase 11** - Acoustics utilities + example (combustor design)

### **Optional Phase:**
8. **Phase 8** - Air properties bundle (nice to have)

**Recommended scope:** Phases 1-7 + 9 for maximum ROI (core utilities + examples).
Add Phases 10-11 for domain-specific applications (orifice metering, acoustics).
Phase 8 can be added later if needed.

---

## Notes

- Do NOT commit this TODO.md file itself
- Each phase must be fully tested before moving to next
- Composite functions are the most valuable but require careful testing
- Start with simple geometry helpers to build confidence
- Machine precision testing is critical for composite functions
- Update examples in a separate task after all functions are implemented
