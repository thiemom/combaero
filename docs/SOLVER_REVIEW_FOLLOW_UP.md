# Network Solver Review — Follow-Up

**Date:** 2026-03-15 (initial), **updated** 2026-03-15 (rev 3)
**Scope:** Re-review of implementation against `SOLVER_REVIEW.md` plan
**Baseline:** 258/258 C++ tests, 916/916 Python tests passing

---

## Executive Summary

The derived-state node architecture (Phases 2–3 of SOLVER_REVIEW.md), the C++
smoothness audit, and all Python-side cleanup items have been implemented. All
nodes return only `[P, P_total]` as solver unknowns, with `T` and `Y` derived
via forward propagation with chain-rule relay (including `P_total` sensitivity).
The C++ orifice uses regularized sqrt `dP / (dP² + ε²)^0.25`, all friction
models use regularized `Re_eff = √(Re² + ε²)`, and unused C++ functions have
been removed. Python dead code (`abs+bias`, `pressure_loss_frac`,
`enable_mixing`, debug block, stale TODOs/docstrings, pipe dead code) has been
cleaned up.

**One accuracy gap remains:** the velocity-of-approach factor `E = 1/√(1 − β⁴)`
is missing from all orifice mass flow paths (solver and standalone). Also, most
of the prescribed test plan has not been implemented, and the bidirectional
orifice interface and `Q` parameter are deferred.

---

## Implementation Status vs SOLVER_REVIEW.md

### Fully Implemented ✅

| Item | Reference |
|------|-----------|
| **P0-1:** `_default_Y` uses mass fractions | `solver.py:55` |
| **P0-2:** `_build_x0` matches `.Y[` pattern | `solver.py:83` |
| **P1-4:** Y clamping before thermo conversion | `solver.py:302-308` |
| **P1-5:** Debug exception printing removed | `solver.py:527-529` |
| **P1-6:** `strict=False` in hot path | `solver.py:261,290` |
| **P2-9:** `HeatExchangerElement` removed | — |
| **Phase 2:** All nodes `[P, P_total]` only | `components.py:120,176,312` |
| **Phase 2:** `compute_derived_state()` on all node types | `components.py:122,178,314,226,271` |
| **Phase 3:** Topological ordering | `solver.py:114-141` |
| **Phase 3:** Forward state propagation | `solver.py:143-226` |
| **Phase 3:** Chain-rule relay incl. `d_P_total` | `solver.py:172-214` |
| **Phase 3:** `_get_node_state` uses cached derived state | `solver.py:284-288` |
| **Phase 3:** Derived states exported in solution dict | `solver.py:563-569` |
| **Orifice Python:** `if/else` branching removed | `components.py:382-420` |
| **Pipe Python:** `if/else` branching removed, dead code cleaned | `components.py:676-718` |
| **C++ Smoothness:** Regularized `√(dP)` for orifice | `solver_interface.cpp:31-51` |
| **C++ Smoothness:** Regularized Re in all friction models | `solver_interface.cpp:157-158` (haaland), `:200-201` (serghides), `:265-266` (colebrook) |
| **C++ Smoothness:** Orifice accepts signed dP (no throw) | `solver_interface.cpp:31-51` |
| **C++ Cleanup:** `combustor_residuals_and_jacobian` removed | `solver_interface.cpp`, `solver_interface.h` |
| **C++ Cleanup:** `plenum_residuals_and_jacobian` removed | `solver_interface.cpp`, `solver_interface.h` |
| **Python Cleanup:** `abs(m_dot) + 1e-10` removed | `components.py:131,187,324` — now `cb.MassStream(s.m_dot, ...)` |
| **Python Cleanup:** `pressure_loss_frac` removed from `CombustorNode` | `components.py:296-300` |
| **Python Cleanup:** `enable_mixing` flag removed | `components.py:114,154` |
| **Python Cleanup:** Debug equation name block removed | `solver.py` — no longer present |
| **Python Cleanup:** Stale `NetworkSolver` docstring updated | `solver.py:25-32` |
| **Python Cleanup:** `MomentumChamberNode` stale TODO removed | `components.py:191-195` |

### Not Implemented ❌

| Item | SOLVER_REVIEW Section | Impact |
|------|----------------------|--------|
| Velocity-of-approach factor `E = 1/√(1-β⁴)` | Orifice physics | Up to 14% error at β=0.7 |
| C++ bidirectional orifice interface (both T/Y) | Smoothness Audit | Single-sided thermo only |
| Node-level `Q` parameter | Design §Universal Node Model | No heat exchange capability |
| 17 of 19 prescribed tests | Test Plan | See test gap table below |

---

## C++ Smoothness Audit — Verified

### Orifice: Regularized `dP / (dP² + ε²)^0.25`

The old code threw on `dP < 0` and had a `{0.0, 1e12}` discontinuity at `dP = 0`.
The new code (`solver_interface.cpp:31-51`) implements:

```
f_reg(dP) = C · dP / (dP² + ε²)^0.25
```

with `ε² = 1.0` (1 Pa²). This is:
- **C∞ smooth** everywhere, including at `dP = 0`
- **Sign-preserving**: `f_reg(dP)` has the same sign as `dP`
- **Asymptotically correct**: for `|dP| >> 1`, converges to `C · sign(dP) · √|dP|`

The analytical derivative is also correctly implemented:
```
df_reg/d(dP) = C · (0.5·dP² + ε²) / (dP² + ε²)^1.25
```

The `orifice_residuals_and_jacobian` function (`solver_interface.cpp:973-1009`)
now passes signed `dP` directly to `orifice_mdot_and_jacobian` — no `abs()`,
no `sign()` multiply, no derivative sign flip. This eliminates the Jacobian
inconsistency documented in SOLVER_REVIEW.md §Smoothness Layer 2.

**Remaining issue (minor):** The function still only accepts upstream thermo state
(`T_up, Y_up`). A fully bidirectional interface would accept both sides and
blend density. This is a minor accuracy issue (not a smoothness issue) since
the regularized sqrt already handles the sign correctly.

### Missing Velocity-of-Approach Factor `E = 1/√(1 - β⁴)` (P0-accuracy)

The ISO 5167 orifice mass flow equation is:

```
q_m = C · E · ε · A · √(2·ρ·ΔP)
```

where **E = 1/√(1 - β⁴)** is the velocity-of-approach factor and **β = d/D**
(orifice bore diameter / pipe diameter). Both orifice code paths omit E:

| Code path | Formula | Missing |
|-----------|---------|---------|
| `solver_interface.cpp:orifice_mdot_and_jacobian` | `Cd · A · f_reg(dP, ρ)` | E |
| `orifice.cpp:orifice_mdot` | `Cd · A · √(2ρ·dP)` | E (and ε) |
| `orifice.cpp:solve_orifice_mdot` | `Cd · ε · A · √(2ρ·dP)` | E |

The Cd correlations (Reader-Harris/Gallagher etc.) return the ISO discharge
coefficient **C**, which does **not** absorb E.

**Impact by β:**

| β | E | Error |
|---|---|-------|
| ≤ 0.3 | ≤ 1.004 | negligible |
| 0.5 | 1.033 | 3.3% |
| 0.6 | 1.072 | **7.2%** |
| 0.7 | 1.143 | **14.3%** |

**`OrificeElement.resolve_topology()` discovers pipe diameters but never uses them:**

```python
# components.py:355-356 — discovered but dead
self.upstream_diameter: float | None = None
self.downstream_diameter: float | None = None
```

`residuals()` (line 382) never references these. β could be computed as
`√(4·A_bore / (π·D²))` but isn't.

**Interaction with regularization:** The regularization `dP/(dP²+ε²)^0.25` is
orthogonal to β. E is a constant geometry multiplier (β doesn't change during
a solve), so all derivatives simply scale by E — no new Jacobian columns needed.

**Recommended fix:** Add optional `beta` parameter to `orifice_mdot_and_jacobian`:

```cpp
std::tuple<double, double> orifice_mdot_and_jacobian(
    double dP, double rho, double Cd, double area, double beta = 0.0) {
  double E = (beta > 0.0) ? 1.0 / std::sqrt(1.0 - std::pow(beta, 4.0)) : 1.0;
  double constant = Cd * E * area * std::sqrt(2.0 * rho);
  // ... rest unchanged
```

Then have `OrificeElement.residuals()` compute β from `self.upstream_diameter`
and pass it through. This also fixes `orifice_residuals_and_jacobian` and the
`dmdot_drho` derivative chain. The same fix should be applied to
`orifice.cpp:orifice_mdot` and `solve_orifice_mdot`.

**Note:** `EffectiveAreaConnectionElement` and `AreaDischargeCoefficientConnectionElement`
skip topology resolution and document that no β correction is needed — they
assume the user provides an already-corrected effective area. This is correct.
The `AreaDischargeCoefficientConnectionElement` relationship `ζ = 1/Cd² − 1` is
the β→0 limit; for non-zero β the correct form is `K = (1/Cd² − 1)(1 − β⁴)`
(see `orifice.cpp:K_from_Cd`). This should be documented but is not a bug since
those elements don't claim to handle β.

### Friction: Regularized Reynolds Number

All explicit friction models (haaland, serghides, colebrook) now use:
```
Re_eff = √(Re² + ε_Re²)    with ε_Re = 10
```

This replaces the old `max(Re, 1.0)` hard floor. The derivative
`dRe_eff/dRe = Re / Re_eff` is correctly propagated through the chain rule
in all three models. The `1/Re` singularity at zero flow is eliminated smoothly.

**Verified:** C++ Jacobian FD tests (`test_solver_jacobians.cpp`) pass for both
orifice and pipe, confirming analytical derivatives match finite differences.

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
| Smoothness tests (5 prescribed) | ⚠️ | C++ regularization done; Python-level smoothness tests not yet written |

### Existing Tests That Cover the Refactor

| Test | What It Proves |
|------|----------------|
| `test_combustor_network_conservation` | Combustor derived state + chain-rule works end-to-end |
| `test_jacobian_accuracy` | Analytical Jacobian matches FD for pressure-only network |
| `test_flow_reversal_simple_orifice` | Negative m_dot converges (now via smooth C++ orifice) |
| `test_flow_reversal_momentum_chamber` | MomentumChamber with reversal converges |
| `test_step_1..4_*` scenarios | Various topologies with derived state |
| C++ `SolverJacobianTest.OrificeDerivatives` | Regularized orifice Jacobian matches FD |
| C++ `SolverJacobianTest.PipeDerivatives` | Pipe Jacobian with regularized friction matches FD |

### Recommended Test Priority

**High — should exist before further refactoring:**
1. `test_chain_rule_orifice_through_combustor` — proves the relay is correct
2. `test_plenum_temperature_passthrough` — proves derived-state passthrough works
3. `test_all_nodes_only_P_Ptotal_unknowns` — structural assertion
4. `test_single_stream_passthrough` — proves identity chain rule
5. `test_orifice_smoothness_at_zero_dP` — proves regularization works end-to-end

**Medium — after cleanup:**
6. `test_mixing_plenum_convergence`
7. `test_equilibrium_combustor_convergence`
8. `test_derived_state_jacobian_fd`

**Blocked on features:**
9–11. Q tests (blocked on Q implementation)

---

## Prioritized Action Items

### P0 — Accuracy

| # | Item | Effort | Files |
|---|------|--------|-------|
| 1 | Add velocity-of-approach factor `E = 1/√(1-β⁴)` to orifice — solver and standalone paths | M | `solver_interface.cpp`, `solver_interface.h`, `orifice.cpp`, `components.py`, `_core.cpp` |

### P1 — Tests

| # | Item | Effort | Files |
|---|------|--------|-------|
| 2 | Add high-priority tests (chain-rule, passthrough, system size, smoothness, β-correction) | M | `python/tests/` |

### P2 — Architecture Improvements

| # | Item | Effort | Files |
|---|------|--------|-------|
| 3 | C++ bidirectional orifice interface (both T/Y) | M | `solver_interface.cpp`, `solver_interface.h`, `_core.cpp` |
| 4 | Node-level `Q` parameter for heat exchange | M | `components.py`, `solver.py` |

### P3 — Future

| # | Item | Effort | Files |
|---|------|--------|-------|
| 5 | Implement proper `MomentumChamberNode` with `P_total ≠ P` | M | `components.py` |

---

**Initial review:** 2026-03-15
**Rev 2:** 2026-03-15 (C++ smoothness verified, β-correction gap identified)
**Rev 3:** 2026-03-15 (all Python cleanup verified done; P0 #2–3, P1 #3–8 resolved)
**Next steps:** P0 #1 (β-correction), then P1 #2 (tests), then P2 architecture
