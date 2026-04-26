feat: implement continuation solver and solver method selection

- Exposed solver method options (lm, krylov, etc.) in backend and GUI.
- Implemented 'continuation' initialization strategy to reuse converged solutions.
- Added module-level persistence in FastAPI backend for solver warm-starts.
- Implemented safety resets for topology mismatches and initial guess changes.
- Added automated tests for continuation logic and fingerprinting.
- Updated Sidebar with dynamic continuation availability and reset recovery.
