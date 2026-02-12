# Python Bindings Implementation Plan

## Overview
Add Python bindings for friction, heat transfer, and geometry utility functions.

**Requirements:**
1. Every binding needs a test verifying implementation is correct
2. Units always commented and checked
3. Documentation updated with new bindings
4. Commit after each phase completion (tests must pass)
5. Continue to next phase only if all tests pass

---

## Phase 1: Friction Factor Correlations ✅

**Functions to bind (4):**
- [x] `friction_haaland(Re, e_D)` - Haaland correlation [-]
- [x] `friction_serghides(Re, e_D)` - Serghides correlation [-]
- [x] `friction_colebrook(Re, e_D, tol, max_iter)` - Colebrook-White [-]
- [x] `friction_petukhov(Re)` - Petukhov smooth pipe [-]

**Tasks:**
- [x] Add bindings to `python/combaero/_core.cpp`
- [x] Write unit tests in `python/tests/test_friction.py` (12 tests)
- [x] Verify units in docstrings (dimensionless friction factor)
- [x] Update `docs/API_REFERENCE.md` with friction functions
- [x] Run all tests: `pytest python/tests/` - **12/12 PASSED**
- [x] Commit: "Add Python bindings for friction factor correlations"

**Status:** ✅ COMPLETED

---

## Phase 2: Basic Heat Transfer Correlations ✅

**Functions to bind (6):**
- [x] `nusselt_dittus_boelter(Re, Pr, heating)` - Dittus-Boelter Nu [-]
- [x] `nusselt_gnielinski(Re, Pr, f)` - Gnielinski Nu with friction [-]
- [x] `nusselt_gnielinski(Re, Pr)` - Gnielinski Nu auto-friction [-]
- [x] `nusselt_sieder_tate(Re, Pr, mu_ratio)` - Sieder-Tate Nu [-]
- [x] `htc_from_nusselt(Nu, k, L)` - Heat transfer coefficient [W/(m²·K)]
- [x] `lmtd(dT1, dT2)` - Log mean temperature difference [K]

**Tasks:**
- [x] Add bindings to `python/combaero/_core.cpp`
- [x] Write unit tests in `python/tests/test_heat_transfer.py` (18 tests)
- [x] Verify units in docstrings (Nu dimensionless, h in W/(m²·K), LMTD in K)
- [x] Update `docs/API_REFERENCE.md` with heat transfer functions
- [x] Run all tests: `pytest python/tests/` - **18/18 PASSED**
- [x] Commit: "Add Python bindings for heat transfer correlations"

**Status:** ✅ COMPLETED

---

## Phase 3: Geometry Utilities ⏳

**Functions to bind (6):**
- [ ] `hydraulic_diameter(A, P_wetted)` - Generic hydraulic diameter [m]
- [ ] `hydraulic_diameter_rect(a, b)` - Rectangular duct [m]
- [ ] `hydraulic_diameter_annulus(D_outer, D_inner)` - Annular duct [m]
- [ ] `residence_time(V, Q)` - Residence time from volume flow [s]
- [ ] `residence_time_mdot(V, mdot, rho)` - Residence time from mass flow [s]
- [ ] `space_velocity(Q, V)` - Space velocity [1/s]

**Tasks:**
- [ ] Add bindings to `python/combaero/_core.cpp`
- [ ] Write unit tests in `python/tests/test_geometry.py`
- [ ] Verify units in docstrings (m, s, 1/s)
- [ ] Update `docs/API_REFERENCE.md` with geometry functions
- [ ] Run all tests: `pytest python/tests/`
- [ ] Commit: "Add Python bindings for geometry utilities"

**Status:** Not started

---

## Completion Checklist

- [ ] All Phase 1 tests passing
- [ ] All Phase 2 tests passing
- [ ] All Phase 3 tests passing
- [ ] Documentation updated for all new functions
- [ ] All commits pushed to main
- [ ] Total functions added: 16

---

## Notes

- Do NOT commit this TODO.md file itself
- Each phase must be fully tested before moving to next
- Units must be explicitly documented in all docstrings
- Follow existing binding patterns in `_core.cpp`
