feat(network): Step 4 - Solver wall coupling via EnergyBoundary + Jacobian relay

Architecture (per design document):
- Evaluate walls BEFORE compute_derived_state() so Q flows through
  the existing EnergyBoundary abstraction (single mixing pass)
- Create dedicated wall EnergyBoundary per affected node during setup
- Update wall EB Q each Newton iteration before mixing

Wall evaluation (_evaluate_walls_for_node):
- Calls htc_and_T() on both coupled elements
- Calls wall_coupling_and_jacobian() (C++) for Q and Jacobians
- Accumulates Q from multiple walls before setting on EnergyBoundary
- Skips same-node walls (net Q = 0 by energy conservation)
- Deduplicates walls via processed_walls set

Back-edge handling (_get_node_state_with_prev):
- Preserves _prev_derived_states from previous Newton iteration
- Uses lagged T/Y for nodes not yet visited in current topological pass
- Standard Newton behaviour for cross-stream coupling

Jacobian relay (first-order temperature coupling):
- dT_node/dx += dT_mix/dQ * sign * dQ/dT_aw_a * dT_A/dx
- dT_node/dx += dT_mix/dQ * sign * dQ/dT_aw_b * dT_B/dx
- Uses dQ/dT_aw ~ dQ/dT approximation (recovery factor ~ 1)
- Mass-flow coupling (dQ/dmdot) deferred until channel Jacobians available

All 1019 tests pass.
