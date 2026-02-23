# CombAero AI/Agent Development Guide

This document outlines the hard rules, coding styles, and project conventions for AI agents (like Gemini, Claude, or Cursor) working on the CombAero project. **These rules must be strictly followed.**

## 1. Hard Rules & Constraints

*   **No system Python installs:** Agents **must never** install packages into the system Python or use global `pip`.
*   **Virtual Environment Only:** Always use the project-local virtual environment located at `.venv`. If activating a shell, use `source .venv/bin/activate`.
*   **Do not bypass checks:** Never bypass Git pre-commit hooks or GitHub Actions runner checks. All code must pass these checks cleanly.
*   **Auto-generated files:** Do not manually edit auto-generated files. They will be overwritten. Use the designated scripts to regenerate them.

## 2. Coding Style

### C++ (Modern C++17)

*   **Standard:** Use C++17 features (the project baseline).
*   **Memory Management:** Prefer smart pointers (`std::unique_ptr`, `std::shared_ptr`) over raw pointers. Never use `new` or `delete` manually. Be strictly RAII compliant.
*   **Header Hygiene:** Use `#pragma once` in all header files. Keep includes sorted and minimal (include what you use).
*   **Type Safety:** Use `auto` for complex iterator types, but keep variable types explicit where it aids readability.
*   **Comments:** Use line comments (`//`) only. Block comments (`/* ... */`) are flagged by the style checker.
*   **Math Constants:** Use `math_constants.h` instead of macros like `M_PI`.
*   **Encoding:** Only ASCII characters are allowed outside of strings and comments.

### Python (Type-Safe & Clean)

*   **Type Hinting:** Mandatory type annotations for all function signatures (arguments and return types). Avoid `Any` unless absolutely necessary.
*   **Modern Typing:** Use the `|` operator for unions (e.g., `str | None`) instead of `Union` or `Optional` (project targets Python 3.12).
*   **Explicit Returns:** Always include `-> None` for functions that do not return a value.
*   **LBYL (Look Before You Leap):** Proactively check conditions rather than relying on `try/except` for control flow.
*   **Formatting/Linting:** The project uses `ruff` for both linting and formatting. Ensure code complies by running the formatters locally.
*   **Encoding:** Strictly no non-ASCII characters allowed anywhere in Python files (including docstrings and comments).

## 3. Auto-Generated Files & Scripts

Several critical files in this project are auto-generated. **Do not modify the targets directly.** Modify the source data or generation script, then run the proper updater.

### Units Documentation (`docs/UNITS.md`)
*   **Source:** `include/units_data.h`
*   **Script:** `python scripts/generate_units_md.py`
*   **Action:** If fixing a unit description or adding a section, modify the C++ header, then run the python script to regenerate the markdown.

### Thermodynamic & Transport Data (`include/thermo_transport_data.h`)
*   **Source:** NASA CEA output text files and Cantera YAML mechanisms (external), combined into a local JSON cache.
*   **Generator Directory:** `thermo_data_generator/`
*   **Action:** To update the C++ thermo data header, use the scripts inside `thermo_data_generator/` (e.g., `generate_thermo_data.py`). See `thermo_data_generator/README.md` for the exact multi-step process.

## 4. Environment & Tooling Usage

### The Virtual Environment (`.venv`)
The root `.venv` is the single source of truth for all Python operations in this repository:
*   Building the C++ Python wheel binding (`scikit-build-core`)
*   Running C++/Python tests via `pytest`
*   Running the thermo data generators
*   Code style checks (`ruff`, `mypy`) via `./scripts/check-python-style.sh`

**To ensure the environment is correct:**
```bash
./scripts/bootstrap.sh
source .venv/bin/activate
```

### Script Execution and Verification
Before committing code, verify it locally using the provided scripts:

*   **Check C++ Style:** `./scripts/check-source-style.sh`
*   **Check Python Style (and format):** `./scripts/check-python-style.sh --fix`
*   **Run Clang-Tidy:** `./scripts/run-clang-tidy.sh` (requires a CMake build directory with `compile_commands.json`)
*   **Run Test Suites:**
    *   C++: `cd build && ctest`
    *   Python: `./.venv/bin/python -m pytest python/tests`
    *   Cantera Validation: `CANTERA_DATA="$PWD/cantera_validation_tests" ./.venv/bin/python -m pytest cantera_validation_tests -v`

### Version Synchronization
When updating tool versions or Python targets, ensure the following files remain synchronized:
1.  `.python-version`
2.  `ruff.toml` (`target-version`)
3.  `.github/workflows/ci.yml`
4.  `.github/workflows/validation-tests.yml`
5.  `.pre-commit-config.yaml`
6.  `Makefile` and the `.sh` scripts in `scripts/`
