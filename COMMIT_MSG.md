fix: state property units mismatch and example runner automation

- Fixed critical units mismatch in `_core.cpp` where `State` properties (SP, HP, etc.) were double-converting mass-based units.
- Created `scripts/run_examples.py` for automated discovery and verification of Python examples.
- Simplified `network_from_grid.py` and moved intensive Mach sweep to `benchmarks/mach_sweep_network.py`.
- Fixed multiple API regressions in examples and corrected orifice dimensions.
