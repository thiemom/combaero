fix: correct heat transfer Jacobian bugs (106% → 0% FD error)

Fixed multiple critical bugs in heat transfer Jacobian calculations that caused
the wall coupling Jacobian validation test to fail with 106.2% relative error.

## C++ Changes (src/heat_transfer.cpp + include/heat_transfer.h)

1. **f_multiplier ordering**: Apply `f *= f_multiplier` BEFORE Nu computation
   so Gnielinski/Petukhov correlations (which use f in their formulas) are
   consistent with the dNu/dRe finite difference stencils.

2. **Laminar df_dRe**: Remove double-application of f_multiplier. Since f
   already includes the multiplier, `df_dRe = -f/Re` is correct.

3. **dNu/dPr FD**: Revert commit 7240bb2 — now that f is multiplied before
   Nu computation, the FD stencil should use f directly (not f_raw).

4. **Petukhov df_dRe**: Analytical formula applies to raw f, not multiplied f.
   Restored `pow(f/f_multiplier, 1.5) * f_multiplier`.

5. **Add dT_aw derivatives**: Added `dT_aw_dT` and `dT_aw_dmdot` fields to
   ChannelResult struct. Computed unconditionally (moved outside T_wall guard)
   so Python relay can use them. Applied to all 5 channel functions.

## Python Changes (solver.py)

6. **ROOT CAUSE — Unit mismatch**: `dT_mix_d_delta_h = 1/cp` is dT/d(specific_h)
   in K/(J/kg), but wall Q is total heat rate in Watts. Fixed by dividing by
   total mass flow: `dT_mix_dQ = dT_mix_d_delta_h / m_dot_total`.

7. **Use actual dT_aw derivatives**: Relay now uses `ch.dT_aw_dT` and
   `ch.dT_aw_dmdot` from ChannelResult instead of hardcoded 1.0 and 0.0.

## Other Changes

- Added pybind11 bindings for new fields (_core.cpp)
- Added units metadata for dT_aw_dT and dT_aw_dmdot (units_data.h)
- Added integration test (test_wall_coupling_integration.py)

## Test Results

- C++ tests: 15/15 pass
- Python tests: 1013/1014 pass (1 pre-existing units_sync failure)
- **Jacobian validation: 0.0% max relative error** (was 106.2%)
