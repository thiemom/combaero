feat: redesign orifice element and harden solver logic

1. Resolved Python-C++ binding compatibility by strictly using C++ fallback implementation logic within python component space for Stolz and Miller correlations.
2. Restored internal exception raising logic for T <= 0 and P <= 0 within transport_state C++ bindings, preventing failure silently defaulting to previous states.
3. Adjusted network element resolution to gracefully skip correlation fallback evaluation for non-aligned data parameters and removed "Auto Cd" usage fallback tests.
4. Cleared out formatting and strict typing anomalies caught by validation suites. All Cantera validation checks pass.
