feat: stabilize compressible network solver and unify nomenclature

- Renamed PipeElement/pipe_flow to ChannelElement/channel_flow project-wide
- Provided aliases for backward compatibility (PipeElement, pipe_flow, etc.)
- Relaxed strict Mach validation in fanno_channel_rough (src/compressible.cpp)
- Fixed thermodynamic consistency by using T_total in element residuals
- Updated unit metadata and documentation (API_PYTHON.md, API_CPP.md, UNITS.md)
- Verified all regressions pass in test_network_compressible.py领域,Description:
