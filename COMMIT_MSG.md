feat: expose solver method selection in the GUI sidebar

Exposed the `method` parameter of `NetworkSolver.solve()` (e.g., `hybr`, `lm`, `krylov`) to the user via the "Solver Configuration" sidebar.
- Added `method` to `SolverSettings` Pydantic model with Literal validation.
- Updated FastAPI solve handler to pass the selected method to the solver.
- Integrated `method` into the frontend Zustand store.
- Added a dropdown menu in the Sidebar component for algorithm selection.
