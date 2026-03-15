# Network Solver Review — Follow-Up

**Date:** 2026-03-15
**Scope:** Re-review of implementation against `SOLVER_REVIEW.md` plan
**Baseline:** 258/258 C++ tests, 916/916 Python tests passing

---

## Executive Summary

The derived-state node architecture (Phases 2–3 of SOLVER_REVIEW.md) has been
successfully implemented. All nodes return only `[P, P_total]` as solver unknowns,
with `T` and `Y` derived via forward propagation with chain-rule relay. The
combustor network test that was previously deselected now passes.

Key gaps remain: C++ smoothness fixes were not applied, the `Q` heat exchange
parameter was not added, combustor pressure loss is dead code, and most of the
prescribed test plan was not implemented.

---

## Implementation Status vs SOLVER_REVIEW.md

### Fully Implemented ✅

| Item | Reference |
|------|-----------|
| **P0-1:** `_default_Y` uses mass fractions | `solver.py:84` |
| **P0-2:** `_build_x0` matches `.Y[` pattern | `solver.py:112` |
| **P1-4:** Y clamping before thermo conversion | `solver.py:311-317` |
| **P1-5:** Debug exception printing removed | `solver.py:558-560` |
| **P1-6:** `strict=False` in hot path | `solver.py:270` |
| **P2-9:** `HeatExchangerElement` removed | — |
| **Phase 2:** All nodes `[P, P_total]` only | `components.py:121,182,332` |
| **Phase 2:** `compute_derived_state()` on all node types | `components.py:123,184,334,244,289` |
| **Phase 3:** Topological ordering | `solver.py:143-170` |
| **Phase 3:** Forward state propagation | `solver.py:172-235` |
| **Phase 3:** Chain-rule relay in Jacobian assembly | `solver.py:355-371, 434-449` |
| **Phase 3:** `_get_node_state` uses cached derived state | `solver.py:292-297` |
| **Phase 3:** Derived states exported in solution dict | `solver.py:591-597` |
| **Orifice Python:** `if/else` branching removed | `components.py:402-438` |
| **Pipe Python:** `if/else` branching removed | `components.py:689-735` |

### Not Implemented ❌

| Item | SOLVER_REVIEW Section | Impact |
|------|----------------------|--------|
| C++ regularized `√(dP)` for orifice | Smoothness Audit | Orifice Jacobian wrong for dP < 0 |
| C++ bidirectional orifice interface | Smoothness Audit | Single-sided thermo only |
| C++ `√(Re²+1)` for pipe Re floor | Smoothness Audit | Kink at Re = 1 |
| Phase 1 C++ pressure-only residual functions | Implementation Plan §1 | Stagnation handled in Python instead |
| Node-level `Q` parameter | Design §Universal Node Model | No heat exchange capability |
| Phase 4 cleanup (dead code removal) | Implementation Plan §4 | Stale code remains |
| 17 of 19 prescribed tests | Test Plan | See test gap table below |

---

## New Issues Found

### 1. `abs(s.m_dot) + 1e-10` workaround (P0)

All three node types use this in `compute_derived_state`:

```python
# components.py:133, 194, 345
streams = [
    cb.MassStream(abs(s.m_dot) + 1e-10, s.T_total, s.P_total, s.Y) for s in upstream_states
]
```

**Problems:**
- `abs()` destroys flow direction — reverse flow is invisible to the mixer
- `1e-10` bias shifts results for small flows (e.g., pilot fuel at 0.001 kg/s)
- The guard should be inside the C++ mixer, not applied by every caller

**Fix:** Add zero-flow handling in C++ `mixer_from_streams_and_jacobians` (return
first stream's properties when `mdot_tot ≈ 0`), then remove `abs()` and bias
from Python callers.

### 2. `CombustorNode.pressure_loss_frac` is dead code (P0)

Defined at `components.py:322` but **never used**. The combustor residual is just
`P_total - P = 0` (same as plenum) — there is no pressure loss equation.

Users setting `pressure_loss_frac=0.05` get zero effect. The old C++
`combustor_residuals_and_jacobian` had pressure loss, but it's no longer called.

**Fix:** Either:
- (A) Add a pressure loss residual: `P_total_out = P_total_in_avg * (1 - loss_frac)`,
  making the combustor have 3 residuals (stagnation, pressure loss, mass conservation)
  and 3 unknowns (`P, P_total, P_total_in` or similar), OR
- (B) Remove `pressure_loss_frac` and document that combustor pressure loss should
  be modeled as a separate element (pipe or orifice upstream/downstream)

Option (B) is cleaner for the pure pressure-flow architecture.

### 3. Debug equation name block in hot path (P1)

```python
# solver.py:453-473
# DEBUG ATTACH NAMES
if not hasattr(self, "_debug_eqn_names") or len(self._debug_eqn_names) != len(res):
    names = []
    for node in self.network.nodes.values():
        names.append(f"{node.id} mass_con")
    for node in self.network.nodes.values():
        if isinstance(node, cb.network.components.CombustorNode):
            names.append(f"{node.id} t")
            for i in range(14):
                names.append(f"{node.id} X[{i}]")
            ...
```

**Problems:**
- Runs on every residual evaluation (first call per size change)
- References old T/Y/P_total equation ordering (pre-refactor)
- Imports `CombustorNode` via string path in hot path

**Fix:** Remove entirely or move to a `debug=True` solver option.

### 4. `PipeElement.residuals()` dead code (P1)

```python
# components.py:694-697
dP_fwd = state_in.P_total - state_out.P

if dP_fwd > 0 or (dP_fwd == 0 and m_dot > 0):
    import combaero as cb
```

`dP_fwd` is computed but the `if` branch only contains a redundant `import`.
Leftover from removed branching logic.

### 5. Stale `enable_mixing` flag (P1)

`PlenumNode.__init__` and `MomentumChamberNode.__init__` still accept
`enable_mixing` (`components.py:114, 163`). This flag is vestigial — all nodes
now derive state uniformly via `compute_derived_state()`. It only affects
`resolve_topology` (whether to store upstream elements), but the solver's
`_propagate_states` discovers upstream elements independently.

**Fix:** Remove the parameter. Always resolve upstream elements in
`resolve_topology`.

### 6. Stale `NetworkSolver` docstring (P1)

```python
# solver.py:30-60
#    CRITICAL TODOs for Phase 2 Completion (Combustion & Mixing):
#    ================================================================
#    1. IMPLEMENT PROPER COMBUSTION RESIDUALS in CombustorNode:
#    ...
```

Phase 2 is complete. This 30-line TODO block is misleading.

### 7. `MomentumChamberNode.residuals()` stale TODO (P1)

```python
# components.py:199-208
"""
TODO: This is a placeholder implementation.

For Phase 2+, this should implement proper energy and species conservation
using cb.solver.enthalpy_and_jacobian() and other solver_interface.h functions.
"""
```

Phase 2+ is complete. The momentum chamber currently uses the same `P_total = P`
stagnation as a plenum. If this is intentional (momentum handled elsewhere),
remove the TODO. If not, implement `P_total = P + 0.5 * rho * v^2`.

### 8. Relay missing `d_P_total` sensitivity (P2)

`_propagate_states` (solver.py:220-233) propagates `T` and `Y` sensitivities
through `d_T` and `d_Y`, but **never propagates `d_P_total`** from the
`StreamJacobian`. For most networks this is irrelevant (mixer T/Y don't depend
on P_total), but equilibrium combustion where `T_ad` depends on pressure has
a missing chain-rule term.

### 9. Unused C++ functions still compiled (P3)

`combustor_residuals_and_jacobian` and `plenum_residuals_and_jacobian` in
`solver_interface.cpp` are no longer called from Python. They consume compile
time and maintenance attention.

---

## Test Gap Analysis

### SOLVER_REVIEW.md Prescribed Tests

| Test | Status | Notes |
|------|--------|-------|
| `test_mixer_derived_state_values` | ❌ | Verify T_mix, Y_mix for 2-stream mix |
| `test_combustor_derived_state_values` | ❌ | Verify T_ad, Y_burned match standalone |
| `test_single_stream_passthrough` | ❌ | 1 stream → T_out=T_in, Y_out=Y_in |
| `test_derived_state_jacobian_fd` | ❌ | Python FD validation of MixerResult |
| `test_plenum_with_Q` | ❌ | Q not implemented |
| `test_combustor_with_Q` | ❌ | Q not implemented |
| `test_Q_zero_recovers_adiabatic` | ❌ | Q not implemented |
| `test_chain_rule_orifice_through_combustor` | ❌ | FD validate chain-rule relay |
| `test_chain_rule_pipe_through_plenum` | ❌ | FD validate chain-rule relay |
| `test_chain_rule_single_stream_passthrough` | ❌ | dT/dT_in = 1 through non-mixing plenum |
| `test_chain_rule_cascaded_derived_nodes` | ❌ | Recursive chain through 2 derived nodes |
| `test_all_nodes_only_P_Ptotal_unknowns` | ❌ | Assert system size = 12 not 28 |
| `test_no_negative_Y_ever` | ❌ | Y >= 0 at every iteration |
| `test_plenum_temperature_passthrough` | ❌ | Plenum T ≈ 500K not default 300K |
| `test_mixing_plenum_convergence` | ❌ | 2× MassFlow → Plenum → Orifice |
| `test_equilibrium_combustor_convergence` | ❌ | method="equilibrium" network |
| `test_heat_exchange_node` | ❌ | Q ≠ 0 test (blocked on Q) |
| `test_combustor_network_convergence` | ✅ | `test_combustor_network_conservation` |
| Smoothness tests (5 prescribed) | ❌ | Blocked on C++ regularization |

### Existing Tests That Cover the Refactor

| Test | What It Proves |
|------|----------------|
| `test_combustor_network_conservation` | Combustor derived state + chain-rule works end-to-end |
| `test_jacobian_accuracy` | Analytical Jacobian matches FD for pressure-only network |
| `test_flow_reversal_simple_orifice` | Negative m_dot converges |
| `test_flow_reversal_momentum_chamber` | MomentumChamber with reversal converges |
| `test_step_1..4_*` scenarios | Various topologies with derived state |

### Recommended Test Priority

**High — should exist before further refactoring:**
1. `test_chain_rule_orifice_through_combustor` — proves the relay is correct
2. `test_plenum_temperature_passthrough` — proves Problem 2 is fixed
3. `test_all_nodes_only_P_Ptotal_unknowns` — structural assertion
4. `test_single_stream_passthrough` — proves identity chain rule

**Medium — after cleanup:**
5. `test_mixing_plenum_convergence`
6. `test_equilibrium_combustor_convergence`
7. `test_derived_state_jacobian_fd`

**Blocked on features:**
8–12. Q tests (blocked on Q implementation)
13–17. Smoothness tests (blocked on C++ regularization)

---

## Prioritized Action Items

### P0 — Correctness

| # | Item | Effort | Files |
|---|------|--------|-------|
| 1 | Fix `abs(s.m_dot) + 1e-10` — handle zero flow in C++ mixer | S | `solver_interface.cpp`, `components.py` |
| 2 | Remove dead `pressure_loss_frac` or implement combustor ΔP | S | `components.py` |

### P1 — Dead Code & Stale Docs

| # | Item | Effort | Files |
|---|------|--------|-------|
| 3 | Remove debug equation name block from hot path | S | `solver.py:453-473` |
| 4 | Remove dead `dP_fwd` + conditional import in pipe | S | `components.py:694-697` |
| 5 | Remove vestigial `enable_mixing` flag | S | `components.py` |
| 6 | Update `NetworkSolver` docstring (remove Phase 2 TODOs) | S | `solver.py:30-60` |
| 7 | Update `MomentumChamberNode.residuals()` TODO | S | `components.py:199-208` |
| 8 | Add high-priority tests (chain-rule, passthrough, system size) | M | `python/tests/` |

### P2 — Architecture Improvements

| # | Item | Effort | Files |
|---|------|--------|-------|
| 9 | C++ regularized `√(dP)` for orifice | M | `solver_interface.cpp` |
| 10 | C++ bidirectional orifice interface | M | `solver_interface.cpp`, `solver_interface.h`, `_core.cpp` |
| 11 | Node-level `Q` parameter for heat exchange | M | `components.py`, `solver.py` |
| 12 | Add `d_P_total` to relay propagation | S | `solver.py` |

### P3 — Cleanup

| # | Item | Effort | Files |
|---|------|--------|-------|
| 13 | Remove unused C++ `combustor_residuals_and_jacobian`, `plenum_residuals_and_jacobian` | S | `solver_interface.cpp`, `solver_interface.h`, `_core.cpp` |
| 14 | Implement proper `MomentumChamberNode` with `P_total ≠ P` | M | `components.py` |

---

**Review completed:** 2026-03-15
**Next steps:** Execute P0 items, then P1 cleanup + tests, then P2 architecture
