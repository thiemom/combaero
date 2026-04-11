fix: resolve combustor equivalence ratio bug and enhance probe diagnostics

- Refactored `src/composition.cpp` to correct plenum mixing logic and eliminate negative mass fractions.
- Updated `src/combustion.cpp` for robust thermodynamic coupling and diagnostic extraction.
- Enhanced network component logic in `python/combaero/network/components.py` for backend serialization.
- Updated frontend `ProbeNode` and `Inspector` to visualize equivalence ratios and temperature rise.
- Synchronized `include/units_data.h` and `docs/UNITS.md` with updated physical quantities.
- Added `test_phi_cpp.py` for backend verification.
