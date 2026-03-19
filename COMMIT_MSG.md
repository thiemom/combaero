feat: improve network solver robustness and validation

- Add residual norm guard to catch scipy false-positive convergence
- Add automatic fallback to 'hybr' method when requested method fails
- Add network solver validation example (flow_regime_comparison.py Example 5)
- Fix compressible_vs_incompressible_network.py with realistic parameters
- Add high-Mach convergence benchmark test (xfail for future improvement)

Validation:
- Network solver matches direct brentq reference to 0.05-0.13%
- All compressible network tests pass (7 passed, 1 xfailed benchmark)
- Coupled sweep shows 0-4.82% PR difference with pipe-dominated physics
