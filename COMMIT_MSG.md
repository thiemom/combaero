test: add FD validation for adiabatic and stagnation Jacobians

Added C++ finite difference tests for previously untested Jacobian functions
used in combustor networks and compressible flow calculations.

New tests (all pass with 0.1% tolerance):
- AdiabaticTCompleteDerivatives: validates adiabatic_T_complete_and_jacobian_T
- AdiabaticTEquilibriumDerivatives: validates adiabatic_T_equilibrium_and_jacobians
- T0FromStaticDerivatives: validates T0_from_static_and_jacobian_M
- P0FromStaticDerivatives: validates P0_from_static_and_jacobian_M

These functions are critical for combustor modeling and compressible flow
but had no direct C++ validation. The new tests use central finite differences
and tight tolerances to ensure Jacobian accuracy.

Test results:
- C++ tests: 9/9 pass (5 original + 4 new)
- All tests enabled with tight 0.1% relative tolerance
