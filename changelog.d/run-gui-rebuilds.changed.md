- **`scripts/run-gui.sh` rebuilds a stale C++ extension before launching.**
  The install is editable, so Python changes are live, but the compiled
  `combaero._core` is not. After a pull that touched C++, the backend died at
  import against the old build. The new `scripts/ensure-gui-build.sh`
  rebuilds when a C++ input is newer than the build, then checks that the
  backend imports cleanly, counting combaero's swallowed "Import failed"
  fallback as a failure. It stops with the error instead of launching a broken
  backend. `--check` reports without rebuilding.
