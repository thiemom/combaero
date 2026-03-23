perf: optimize pytest suite by simplifying slow network tests

- Moved intensive high-Mach convergence benchmark to `benchmarks/convergence_benchmark.py`.
- Reduced `test_fully_coupled_...` in `test_network_compressible.py` to a single point at Mach 0.2.
- Reduced `time.sleep` in solver timeout tests from 0.2s to 0.05s.
- Reduced total `pytest` execution time from ~57s to **6.66s** on macOS.
- Verified all 1030+ tests still pass (or xfail as expected).
