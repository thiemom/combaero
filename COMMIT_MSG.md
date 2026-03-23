per: optimize pytest suite by simplifying slow network tests

- Reduced Mach sweep resolution in `test_network_compressible.py` (benchmark: 14 -> 3 points, comparison: 3 -> 2 points).
- Reduced `time.sleep` in solver timeout tests from 0.2s to 0.05s.
- Reduced total `pytest` execution time from ~57s to ~24s on macOS.
