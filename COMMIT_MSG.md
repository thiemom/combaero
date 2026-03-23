fix: compressible momentum chambers and robust friction blending

Implemented compressible momentum chambers and a robust friction model transition between laminar and turbulent regimes.

- Added MomentumBoundary node and compressible support for MomentumChamberNode.
- Implemented smooth C1 blending between laminar (64/Re) and turbulent friction models.
- Added analytical Jacobians for blended friction to solver interface.
- Optimized benchmark suite with dynamic repeats and parameterized grid tests.
- Consolidated grid examples in python/examples.
