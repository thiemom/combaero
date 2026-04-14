refactor: finalize project-wide nomenclature transition from 'Pipe' to 'Channel'

Unified 'Channel' nomenclature across the entire project to improve codebase
consistency and alignment with high-level physical modeling.

Key changes:
- Renamed all 'Pipe' related classes and methods in `_core.cpp` and public Python API.
- Refactored `combaero.incompressible` and `combaero.compressible` submodules.
- Standardized `ChannelElement` as the primary network element for fluid transport.
- Updated `OrificeElement` constructor to support both 'diameter' and 'area' for improved API flexibility.
- Fixed critical attribute resolution bug in `python/combaero/__init__.py`.
- Modernized GUI backend schemas and graph builder to utilize `Channel` terminology.
- Updated comprehensive examples and docstrings to reflect the new standardized nomenclature.
- Verified fix with full test suite (1051/1051 items passed).
