# CombAero TODO - Pending Work

## Overview
Remaining utility functions and improvements to implement.

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

---

## Phase 2: Mass-Specific Thermodynamic Properties ✅ COMPLETE

**Priority:** HIGH - Frequently needed, simple conversions
**Difficulty:** EASY - Simple wrappers around existing functions

**Rationale:** Maintain symmetric interface - if we have molar-basis functions (cp, h, cv, s, u),
we should have mass-basis equivalents for all of them to avoid user confusion.

**Functions implemented (5):**
- [x] `cp_mass(T, X)` - Mass-specific heat capacity at constant pressure [J/(kg·K)]
- [x] `cv_mass(T, X)` - Mass-specific heat capacity at constant volume [J/(kg·K)]
- [x] `h_mass(T, X)` - Mass-specific enthalpy [J/kg]
- [x] `s_mass(T, P, X, P_ref)` - Mass-specific entropy [J/(kg·K)]
- [x] `u_mass(T, X)` - Mass-specific internal energy [J/kg] ⭐ ADDED

**Implementation:**
- [x] Add functions to `src/thermo.cpp`
- [x] Add declarations to `include/thermo.h`
- [x] Add Python bindings to `python/combaero/_core.cpp`
- [x] Export in `python/combaero/__init__.py`
- [x] Write unit tests in `python/tests/test_mass_specific_thermo.py`
- [x] **CRITICAL TESTS:**
  - Verify cp_mass(T, X) == cp(T, X) / mwmix(X) * 1000
  - Verify cv_mass(T, X) == cv(T, X) / mwmix(X) * 1000
  - Verify h_mass(T, X) == h(T, X) / mwmix(X) * 1000
  - Verify s_mass(T, P, X) == s(T, P, X) / mwmix(X) * 1000
  - Verify u_mass(T, X) == u(T, X) / mwmix(X) * 1000
  - All match within machine precision ✅
- [x] Test with multiple compositions (air, fuel, products)
- [x] Update `docs/API_REFERENCE.md`
- [x] **Add unit entries to `include/units_data.h` (5 entries)**
- [x] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [x] Run all tests: `pytest python/tests/` - **202 tests pass** ✅
- [x] Commit: "Complete Phase 2: Add internal energy functions (u and u_mass)"

**Testing results:**
- All functions match molar / mwmix * 1000 within machine precision ✅
- Tested with air, fuel mixtures, combustion products ✅
- Verified units [J/(kg·K)] for cp/cv/s, [J/kg] for h/u ✅
- Verified consistency: gamma_mass = cp_mass / cv_mass equals gamma_molar ✅
- Verified u = h - RT relationship for mass basis ✅

**Bonus:** Also added missing `u(T, X)` molar function to Python bindings

**Status:** ✅ COMPLETE (Commit: 078aa14)

---

## Phase 5: Heat Transfer Coefficient in Pipe (Composite) ✅ COMPLETE

**Priority:** HIGH - Very common heat transfer calculation
**Difficulty:** MEDIUM - Composite function, multiple dependencies

**Functions implemented (1):**
- [x] `htc_pipe(T, P, X, v, D, correlation='gnielinski', heating=True, mu_ratio=1.0, roughness=0.0)`
  - Returns: tuple (h [W/(m²·K)], Nu [-], Re [-]) ✅
  - Correlation options: 'gnielinski' (default), 'dittus_boelter', 'sieder_tate', 'petukhov' ✅

**Implementation:**
- [x] Add function to `src/heat_transfer.cpp` (112 lines with full correlation logic)
- [x] Add declaration to `include/heat_transfer.h` (with comprehensive documentation)
- [x] Add Python binding with tuple return (`python/combaero/_core.cpp`)
- [x] Write comprehensive unit tests in `python/tests/test_htc_pipe_composite.py` (17 tests)
- [x] **CRITICAL TESTS:**
  - Verify h matches: Nu * k / D ✅
  - Verify Nu matches correlation (all 4 correlations tested) ✅
  - Verify Re matches: ρvD/μ ✅
  - Test all correlations (Gnielinski, Dittus-Boelter, Sieder-Tate, Petukhov) ✅
  - Test heating vs cooling ✅
  - Test with roughness for friction factor ✅
  - Machine precision decomposition verification ✅
- [x] Update `docs/API_REFERENCE.md`
- [x] **Add unit entry to `include/units_data.h`**
- [x] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`**
- [x] Run all tests: `pytest python/tests/` - **219 tests pass** ✅
- [x] Commit: "Complete Phase 5: Add composite htc_pipe function with correlation selection"
- [x] **Update TODO.md to mark phase complete** ✅

**Testing results:**
- All 4 correlations match manual decomposition within machine precision ✅
- Tested heating vs cooling for Dittus-Boelter ✅
- Tested viscosity correction for Sieder-Tate ✅
- Tested roughness effects with Gnielinski ✅
- Tested laminar flow handling (Re < 2300) ✅
- Tested temperature, velocity, diameter effects ✅
- Edge cases and error handling verified ✅

**Bonus achievements:**
- Fixed all non-ASCII × characters in codebase (replaced with x)
- Automatic property computation from (T, P, X)
- Handles laminar flow automatically
- Clear error messages for invalid inputs

**Status:** ✅ COMPLETE (Commit: e3b2542)

---

## Phase 8: Air Properties Bundle (Optional) ✅ COMPLETE

**Priority:** LOW - Nice to have, but multiple calls aren't too bad
**Difficulty:** MEDIUM - Returns struct/dataclass, needs careful design

**Design Decision:** ✅ **Dataclass-like struct** (C++ struct → pybind11 class binding)
- **Rationale:** IDE autocomplete, type safety, self-documenting, better UX for engineers
- **C++ side:** `struct AirProperties { double rho; double mu; ... };`
- **Python side:** Bound as class with readonly attributes (behaves like frozen dataclass)
- **Future consideration:** Add `.to_dict()` method to all exported dataclasses for serialization

**Functions implemented (1):**
- [x] `air_properties(T, P, humidity=0.0)` - Returns AirProperties struct/class ✅

**Properties implemented (10):**
- [x] `rho` - Density [kg/m³] ✅
- [x] `mu` - Dynamic viscosity [Pa·s] ✅
- [x] `k` - Thermal conductivity [W/(m·K)] ✅
- [x] `cp` - Specific heat at constant pressure [J/(kg·K)] ✅
- [x] `cv` - Specific heat at constant volume [J/(kg·K)] ✅
- [x] `Pr` - Prandtl number [-] ✅
- [x] `nu` - Kinematic viscosity [m²/s] ✅
- [x] `alpha` - Thermal diffusivity [m²/s] ✅
- [x] `gamma` - Heat capacity ratio [-] ✅
- [x] `a` - Speed of sound [m/s] ✅

**Implementation:**
- [x] Create `struct AirProperties` in `include/thermo.h` (45 lines with documentation)
- [x] Implement `air_properties()` function in `src/thermo.cpp` (45 lines with validation)
- [x] Add pybind11 class binding in `python/combaero/_core.cpp` (67 lines)
  - Used `py::class_<AirProperties>` with `.def_readonly()` for each property ✅
  - Added docstrings with units for each property ✅
  - Added `__repr__` for nice debugging output ✅
- [x] Export `AirProperties` class in `python/combaero/__init__.py` ✅
- [x] Write comprehensive unit tests in `python/tests/test_air_properties.py` (32 tests, 13 classes)
- [x] **CRITICAL TESTS:**
  - Each property matches individual function call with machine precision ✅
  - Tested with dry air (humidity=0.0) ✅
  - Tested with humid air (humidity>0.0) ✅
  - Tested at various temperatures (200K - 600K) ✅
  - Tested at various pressures (50kPa - 200kPa) ✅
  - Verified derived properties: nu = mu/rho, alpha = k/(rho*cp), gamma = cp/cv ✅
- [x] Update `docs/API_REFERENCE.md` with struct definition and usage examples ✅
- [x] **Add unit entry to `include/units_data.h`** ✅
- [x] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`** ✅
- [x] Run all tests: `pytest python/tests/` - **251 tests pass** ✅
- [x] Commit: "Add air properties bundle function with AirProperties dataclass" ✅
- [x] **Update TODO.md to mark phase complete** ✅

**Testing results:**
- All 10 properties match individual function calls within machine precision ✅
- Tested with dry air and humid air (50%, 70%, 100% RH) ✅
- Tested temperature effects (200K - 600K range) ✅
- Tested pressure effects (50kPa - 200kPa range) ✅
- Verified derived properties computed correctly ✅
- Edge cases and error handling verified ✅
- Typical engineering conditions tested (standard, hot day, cold day) ✅

**Achievements:**
- IDE autocomplete works perfectly for all properties
- Type-safe attribute access (read-only)
- Nice `__repr__` for debugging
- Comprehensive error messages for invalid inputs
- Supports both dry and humid air
- Single function call replaces 10+ individual calls

**Status:** ✅ COMPLETE (Commit: 48da87d)

---

## Phase 10: Orifice Flow Utilities with Real Gas Support ✅ COMPLETE

**Priority:** MEDIUM - Useful for flow metering and injection applications
**Difficulty:** MEDIUM - Dataclass pattern + real gas corrections

**Rationale:** The orifice module is comprehensive but lacks:
1. Simple helper functions for common calculations
2. Real gas support via compressibility factor Z (for high-pressure applications)
3. Dataclass pattern for bundled orifice flow results
4. Examples demonstrating the Cd correlation system

**Design Decision:** ✅ **Dataclass-like struct for orifice flow results**
- **Rationale:** Consistent with `AirProperties` pattern, IDE autocomplete, type safety
- **C++ side:** `struct OrificeFlowResult { double mdot; double v; double Re_d; double Cd; ... };`
- **Python side:** Bound as class with readonly attributes (behaves like frozen dataclass)
- **Real gas support:** Add optional `Z` parameter (default=1.0) to density-dependent functions

**Functions implemented (1 composite + 4 utilities):**

### Composite Function (returns dataclass):
- [x] `orifice_flow(geom, dP, T, P, mu, Z=1.0, correlation='reader_harris')` → `OrificeFlowResult` ✅
  - **Returns:** Dataclass with `mdot`, `v`, `Re_D`, `Re_d`, `Cd`, `epsilon`, `rho_corrected`
  - **Computes:** All orifice flow properties in one call with real gas correction
  - **Real gas:** Uses `rho_corrected = rho_ideal / Z` for accurate mass flow
  - **Benefits:** Single function call, IDE autocomplete, consistent with Phase 8 pattern

### Utility Functions (simple calculations):
- [x] `orifice_velocity_from_mdot(mdot, rho, d, Z=1.0)` - Velocity through orifice [m/s] ✅
  - Uses `rho_corrected = rho / Z` for real gas
- [x] `orifice_area_from_beta(D, beta)` - Orifice area from beta ratio [m²] ✅
- [x] `beta_from_diameters(d, D)` - Beta ratio from diameters [-] ✅
- [x] `orifice_Re_d_from_mdot(mdot, d, mu)` - Orifice Reynolds number [-] ✅

**OrificeFlowResult Properties (7):**
- [x] `mdot` - Mass flow rate [kg/s] ✅
- [x] `v` - Velocity through orifice throat [m/s] ✅
- [x] `Re_D` - Pipe Reynolds number (based on D) [-] ✅
- [x] `Re_d` - Orifice Reynolds number (based on d) [-] ✅
- [x] `Cd` - Discharge coefficient [-] ✅
- [x] `epsilon` - Expansibility factor [-] (1.0 for incompressible) ✅
- [x] `rho_corrected` - Density corrected for compressibility [kg/m³] ✅

**Implementation:**
- [x] Create `struct OrificeFlowResult` in `include/orifice.h` (67 lines with documentation)
- [x] Implement `orifice_flow()` composite function in `src/orifice.cpp` (106 lines)
  - Computes composition from (T, P) using ideal gas law
  - Applies Z correction: `rho_corrected = rho_ideal / Z`
  - Calls `solve_orifice_mdot()` with corrected density
  - Computes all derived properties (v, Re_D, Re_d)
  - Returns populated `OrificeFlowResult` struct
- [x] Implement 4 utility functions in `src/orifice.cpp` (71 lines total)
  - Added Z parameter to density-dependent functions
  - Simple, machine-precision calculations
- [x] Add pybind11 class binding for `OrificeFlowResult` in `python/combaero/_core.cpp` (140 lines)
  - Used `py::class_<OrificeFlowResult>` with `.def_readonly()` for each property ✅
  - Added docstrings with units for each property ✅
  - Added `__repr__` for nice debugging output ✅
- [x] Add pybind11 bindings for `orifice_flow()` and utility functions ✅
- [x] Export `OrificeFlowResult` class and functions in `python/combaero/__init__.py` ✅
- [x] Write comprehensive unit tests in `python/tests/test_orifice_flow.py` (29 tests, 11 classes)
- [x] **CRITICAL TESTS:**
  - Verified `mdot` matches `solve_orifice_mdot()` with Z correction ✅
  - Verified `v = mdot / (rho_corrected * A)` within machine precision ✅
  - Verified `Re_d = 4*mdot / (π*d*mu)` within machine precision ✅
  - Verified `rho_corrected = rho_ideal / Z` for Z ≠ 1 ✅
  - Tested with Z=1.0 (ideal gas) matches existing functions ✅
  - Tested with Z=0.8 (high-pressure natural gas) shows correct density correction ✅
  - Tested all utility functions match manual calculations ✅
  - Tested round-trip: area → beta → area ✅
- [x] Update `docs/API_REFERENCE.md` with struct definition and usage examples ✅
- [x] **Add unit entries to `include/units_data.h` (5 entries: 1 composite + 4 utilities)** ✅
- [x] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`** ✅
- [x] Run all tests: `pytest python/tests/` - **280 tests pass** ✅
- [x] Commit: "Complete Phase 10: Add orifice flow utilities with real gas support" ✅
- [x] **Update TODO.md to mark phase complete** ✅

**Testing results:**
- All 7 properties match individual function calls within machine precision ✅
- Tested with ideal gas (Z=1.0) and real gas (Z=0.7-0.9) ✅
- Tested with realistic orifice geometries (beta = 0.2 to 0.7) ✅
- Tested at various pressures: low (1 bar), medium (5 bar), high (100 bar) ✅
- Verified derived properties computed correctly ✅
- Edge cases and error handling verified ✅
- Typical engineering applications tested (natural gas metering, air flow) ✅

**Achievements:**
- IDE autocomplete works perfectly for all properties
- Type-safe attribute access (read-only)
- Nice `__repr__` for debugging
- Comprehensive error messages for invalid inputs
- Real gas support via Z factor (critical for high-pressure applications)
- Single function call replaces multiple calculations
- Consistent with AirProperties pattern from Phase 8
- Backward compatible (Z=1.0 default)

**Status:** ✅ COMPLETE (Commit: 663501b)

---

## Phase 11: Acoustics Utilities ✅ COMPLETE

**Priority:** MEDIUM - Useful for combustor design and instability analysis
**Difficulty:** EASY - Simple frequency and wavelength calculations + dataclass pattern

**Rationale:** The acoustics module is comprehensive but lacks simple helper functions and a convenient bundle for acoustic properties. Following the dataclass pattern from Phases 8 and 10.

**Functions implemented (1 composite + 5 utilities):**

### Composite Function (returns dataclass):
- [x] `acoustic_properties(f, rho, c, p_rms, p_ref)` → `AcousticProperties` ✅
  - **Returns:** Dataclass with `wavelength`, `frequency`, `impedance`, `particle_velocity`, `spl`
  - **Benefits:** Single function call, IDE autocomplete, consistent with AirProperties pattern

### Utility Functions:
- [x] `wavelength(f, c)` - Acoustic wavelength [m] ✅
- [x] `frequency_from_wavelength(lambda, c)` - Frequency from wavelength [Hz] ✅
- [x] `acoustic_impedance(rho, c)` - Characteristic acoustic impedance [Pa·s/m] ✅
- [x] `sound_pressure_level(p_rms, p_ref)` - SPL in dB [dB] ✅
- [x] `particle_velocity(p, rho, c)` - Particle velocity from pressure [m/s] ✅

**Implementation:**
- [x] Add `AcousticProperties` struct to `include/acoustics.h` (56 lines with documentation) ✅
- [x] Add functions to `src/acoustics.cpp` (98 lines total) ✅
- [x] Add Python bindings to `python/combaero/_core.cpp` (136 lines) ✅
  - Used `py::class_<AcousticProperties>` with `.def_readonly()` for each property ✅
  - Added docstrings with units for each property ✅
  - Added `__repr__` for nice debugging output ✅
- [x] Export in `python/combaero/__init__.py` ✅
- [x] Write unit tests in `python/tests/test_acoustic_properties.py` (42 tests, 12 classes) ✅
- [x] **CRITICAL TESTS:**
  - Verified wavelength = c / f within machine precision ✅
  - Verified frequency_from_wavelength = c / lambda within machine precision ✅
  - Verified acoustic_impedance = rho * c within machine precision ✅
  - Verified SPL = 20*log10(p_rms/p_ref) within machine precision ✅
  - Verified particle_velocity = p / (rho*c) within machine precision ✅
  - Tested consistency: wavelength * frequency = c ✅
  - Tested consistency: impedance * particle_velocity = p_rms ✅
- [x] Update `docs/API_REFERENCE.md` with struct definition and usage examples ✅
- [x] **Add unit entries to `include/units_data.h` (6 entries: 1 composite + 5 utilities)** ✅
- [x] **Run `python scripts/generate_units_md.py` to regenerate `docs/UNITS.md`** ✅
- [x] Run all tests: `pytest python/tests/` - **322 tests pass** ✅
- [x] Commit: "Complete Phase 11: Add acoustic properties utilities with dataclass" ✅
- [x] **Update TODO.md to mark phase complete** ✅

**Testing results:**
- All 5 properties match individual function calls within machine precision ✅
- Tested with standard air conditions (20°C, 1 atm) ✅
- Tested with combustor conditions (1500K, 15 bar) ✅
- Tested with water acoustics ✅
- Tested frequency range: 20 Hz to 20 kHz ✅
- Tested SPL range: 0 dB (threshold of hearing) to 130 dB (threshold of pain) ✅
- Edge cases and error handling verified ✅
- Physical relationships verified (λ·f=c, Z·u=p) ✅

**Achievements:**
- IDE autocomplete works perfectly for all properties
- Type-safe attribute access (read-only)
- Nice `__repr__` for debugging
- Comprehensive error messages for invalid inputs
- Single function call replaces 5+ individual calculations
- Consistent with AirProperties and OrificeFlowResult patterns
- Default parameters for standard reference pressure (20 μPa)

**Status:** ✅ COMPLETE (Commit: [pending])

---

## Completion Status

### Completed:
- ✅ Phase 1: Geometry helpers (pipe_area, annular_area, pipe_volume)
- ✅ Phase 3: Velocity from mass flow (already existed as pipe_velocity)
- ✅ Phase 4: Pressure drop composite (pressure_drop_pipe)
- ✅ Phase 6: Pipe roughness database
- ✅ Phase 7: LMTD convenience wrappers
- ✅ Phase 9: Examples refactored to use convenience functions

### Pending:
- [ ] Phase 2: Mass-specific thermodynamic properties (cp_mass, cv_mass, h_mass, s_mass)
- [ ] Phase 5: Heat transfer coefficient composite (htc_pipe)
- [ ] Phase 8: Air properties bundle (optional)
- [ ] Phase 10: Orifice flow utilities + example
- [ ] Phase 11: Acoustics utilities + example

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

## Estimated Effort (Remaining)

### Completed: ~8-10 hours
- Phase 1, 3, 4, 6, 7, 9 ✅

### Remaining:
- **Phase 2:** 1-1.5 hours (4 simple wrappers)
- **Phase 5:** 2-3 hours (htc_pipe composite)
- **Phase 8:** 2 hours (air_properties - optional)
- **Phase 10:** 1.5 hours (orifice utilities + example)
- **Phase 11:** 1.5 hours (acoustics utilities + example)

**Remaining Core (2, 5):** ~3-4.5 hours
**Remaining Extended (2, 5, 10, 11):** ~6-8 hours
**Remaining Complete (2, 5, 8, 10, 11):** ~8-10 hours

---

## Priority Order (Remaining Work)

### **High Priority (Core):**
1. **Phase 2** - Mass-specific thermo (symmetric interface, frequently needed)
2. **Phase 5** - htc_pipe (high value composite for heat transfer)

### **Medium Priority (Extended):**
3. **Phase 10** - Orifice utilities + example (flow metering applications)
4. **Phase 11** - Acoustics utilities + example (combustor design)

### **Low Priority (Optional):**
5. **Phase 8** - Air properties bundle (nice to have, but not critical)

**Recommended next steps:** Complete Phase 2 and Phase 5 for core functionality.
Add Phases 10-11 for domain-specific applications if needed.
Phase 8 can be deferred or skipped.

---

## Notes

- Do NOT commit this TODO.md file itself
- Each phase must be fully tested before moving to next
- Composite functions are the most valuable but require careful testing
- Start with simple geometry helpers to build confidence
- Machine precision testing is critical for composite functions
- Update examples in a separate task after all functions are implemented
