refactor(gui): consume solver diagnostics directly and eliminate redundant evaluations

- Integrated extract_complete_states into the post-solve path of NetworkSolver.
- Optimized extract_complete_states to reuse internal solver state and avoid redundant mass-to-mole conversions.
- Refactored GUI backend (_solve_sync) to consume pre-calculated diagnostics directly from sol_dict.
- Eliminated redundant cb.MixtureState construction and duplicate node/element diagnostic calls in the backend.
- verified with full test suite and scratch verification scripts.

Phase 2 completion (Items 3 & 4) of REVIEW_NETWORK_GUI_INTEGRATION.md.
