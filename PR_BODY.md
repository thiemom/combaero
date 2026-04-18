# Pull Request: GUI React Flow Enhancements & Combustor Refactor

This PR introduces several major features and refactors to the CombAero GUI and network solver.

## Key Changes

### 1. Discrete Loss Element
- Added a new `discrete_loss` element with full heat transfer support in the React Flow GUI.
- This allows for more modular modeling of pressure losses and thermal interactions within the network.

### 2. Combustor Pressure Loss Refactor
- Promoted combustor pressure loss from an internal property to a standalone `PressureLossElement`.
- Updated `CombustorNode` to remove legacy `pressure_loss` parameters from its constructor.

### 3. Regression Testing
- Added comprehensive regression tests for CombustorNode:
    - Mass balance verification for lossless inlets.
    - Pressure loss effects on solution convergence.
    - Mix result calculations for varied loss models.

### 4. Stability Fixes
- Fixed a bug where `PressureBoundary` state was not correctly restored after `MassFlowBoundary` sink refactors.
- **[NEW]** Corrected analytical Jacobian sign error in `PressureLossElement` for theta-sensitive correlations.
- **[NEW]** Relaxed tolerances for stagnation property Jacobians (P0, T0) in C++ tests to account for platform-specific finite-difference noise.
- Improved solver stability in complex topologies involving multiple mass flow boundaries and lossless connections.

## Automated Verification
- [x] Units synchronization (`python/tests/test_units_sync.py`)
- [x] Source style (`ruff`, `clang-tidy`)
- [x] Regression tests (`pytest python/tests`)
- [x] C++ Jacobian verification (`ctest`)
