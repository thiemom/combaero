# Pull Request: Network GUI (React Flow) + Combustor Pressure Loss Architecture

This PR delivers the complete web-based node-editor GUI for CombAero network modelling, built on
React Flow, together with a FastAPI backend and a restructured combustor pressure loss API.

## Primary Change: Network GUI

A brand-new React Flow frontend and FastAPI backend allow users to build, edit, and
solve flow networks interactively in the browser.

**Frontend (`gui/`)**
- Node-based editor for all component types (Plenum, Combustor, Channel, Orifice, …)
- Live telemetry panel: displays converged state (T, P, m, ϕ, …) overlaid on nodes
- Probe node for point-in-network thermodynamic monitoring
- Multi-handle plenums and thermal connection routing
- User-defined node labels; file open/save/export to DataFrame

**Backend**
- FastAPI server serialises `FlowNetwork` to/from JSON (Pydantic schemas)
- `/solve` endpoint returns flat solution dict; `/state` streams diagnostics
- Results-to-DataFrame utility (`results_to_dataframe`)

**Biome** is used for frontend linting/formatting; pre-commit hook enforces style.

---

## Combustor Pressure Loss Refactor

Combustor pressure loss is now a first-class network element instead of a node-internal
parameter, enabling modular loss modelling and correct Jacobian coupling.

- New `PressureLossElement` + `PressureLossCallable` protocol
- Correlation library: `ConstantFractionLoss`, `ConstantHeadLoss`,
  `LinearThetaFractionLoss`, `LinearThetaHeadLoss`
- `CombustorNode` delegates loss via `_zero_loss` when no element is attached
- 15 regression tests cover mass balance, exact pressure ratios, cold-flow
  fallback, ambiguous-source detection, and override behaviour

---

## Stability Fixes & Analytical Jacobians

- **Root-cause fix for CI failures:** `T0_from_static_and_jacobian_M` and
  `P0_from_static_and_jacobian_M` previously used an asymmetric finite-difference
  (smooth-floored M_minus vs. forward M_plus), producing O(h) error ≈ 0.05 Pa accuracy.
  Replaced with **exact analytical chain-rule derivatives**:

  - `dT0/dM = M · a(T)² / cp_mass(T₀)`
  - `dP0/dM = P₀ · cp_mass(T₀) / (T₀ · R_specific) · dT0/dM`

  Tests now pass with relative error < 1e-9 instead of ~2e-4.

- **PressureLossElement Jacobian sign correction:** `jac[dres/dT_burned]` had an
  inverted sign for linear-theta correlations; corrected via chain-rule analysis.

---

## Automated Verification

- [x] C++ unit tests — all 13 Jacobian tests pass with strict original tolerances
- [x] Python unit tests — all network solver / combustor regression tests pass
- [x] Cantera validation — 602 tests pass
- [x] Units synchronisation (`test_units_sync.py`)
- [x] Source style (`ruff`, clang-format via pre-commit)
