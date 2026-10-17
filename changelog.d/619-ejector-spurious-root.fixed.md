- A compressible `EjectorElement` no longer converges to a spurious zero-flow state with its mixing-plane static at ~1e13 Pa. The stress sweep found it in 9 of 11 "converged" cold draws, all reported as success. The C-D nozzle closure behind R0/R1 is now strictly monotone in `P_py`:
  - above the supply it passes reverse flow, where it used to return zero flow flat in `P_py`;
  - below the critical ratio the exit Mach is held at 1, where the supersonic-exit flux used to dip under the choked cap.

  Converged suite states moved by at most 4e-7 in primary flow (#619).
