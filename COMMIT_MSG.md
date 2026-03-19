feat: reorganize examples and create paired combustor cooling implementations

Phase 2 preparation: move examples to python/examples and create network vs standalone comparison.

## Changes

**1. Example file organization**
- Moved `examples/combustor_network_example.py` → `python/examples/`
- Moved `examples/wall_coupling_simple.py` → `python/examples/`
- Renamed `combustor_liner_cooling.py` → `combustor_liner_cooling_standalone.py`

**2. Paired combustor cooling implementations**
- **Standalone version** (`combustor_liner_cooling_standalone.py`):
  - Uses direct combaero API calls (channel_smooth, channel_ribbed, etc.)
  - Manual iteration for hot/cold side coupling
  - Bypass sweep to find minimum cooling flow for T_metal < 800°C

- **Network version** (`combustor_liner_cooling_network.py`):
  - Uses FlowNetwork with ConvectiveSurface and PipeElement
  - Demonstrates network solver with thermal coupling
  - Should produce identical results to standalone for validation

**3. Comparison script**
- Created `compare_combustor_implementations.py`
- Runs both implementations with identical inputs
- Verifies results match within 1% tolerance
- Tests multiple f_cool values and channel types (Smooth, Ribbed, Dimpled)

## Purpose

This paired implementation approach:
- Validates network heat transfer against known standalone calculations
- Provides users with both simple (standalone) and advanced (network) examples
- Demonstrates that FlowNetwork produces correct physics
- Serves as regression test for network solver thermal coupling

## Next Steps

Run comparison script to verify implementations match, investigate any discrepancies.
