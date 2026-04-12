fix(robustness): implement safe density barrier, standardized pi, and species guards

- Added _safe_rho utility with a softplus barrier (and symmetric underflow guard) to ensure numerical robustness in heat transfer models.
- Integrated _safe_rho into PipeElement and MomentumChamberNode htc_and_T methods.
- Standardized PI constants across all components to use math.pi for precision and consistency.
- Implemented defensive clamping guard in MixtureState.X property to handle unphysical mass fractions gracefully during solver iteration.
- updated REVIEW_NETWORK_GUI_INTEGRATION.md (Phase 3 completion).
