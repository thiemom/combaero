# CombAero Python Tests

This directory contains the `pytest` suite for the `combaero` Python bindings and API.

The tests are run via the Github Actions CI pipeline across all major operating systems. You can run them locally via the root virtual environment:
```bash
./scripts/bootstrap.sh
source .venv/bin/activate
python -m pytest python/tests -v
```

## Unit Synchronization (`test_units_sync.py`)

A critical requirement of CombAero is that **every single mathematical calculation, function, or class property exposed to the end-user MUST have an accompanying unit definition** in `include/units_data.h`. This ensures that downstream tools can always query `combaero.has_units("...")` and receive the output unit (e.g., `J/kg` versus `J/mol`).

The `test_units_sync.py` script automatically utilizes Python introspection (`inspect`) on `combaero.__all__` to guarantee this sync. **If the test fails, it means you have exposed a new API endpoint in `combaero` without defining its units in C++.**

### Rules for AI Agents: How to Fix a Failing Sync Test

If you are an AI assistant and `test_units_sync.py` is failing after you implemented a request:

1. **Did you add a mathematical calculation or property?**
   If yes, you **MUST** update `include/units_data.h` with the appropriate unit strings and recompile the C++ module (`pip install -e .`). Do not bypass the test.

2. **Did you add a structural module, helper, or non-physical enum?**
   If you added an entire submodule (e.g., `combaero.new_math_module`), a utility tool (`combaero.suppress_warnings`), or an enum (`combaero.Flag`), these do not generate math outputs and therefore do not have physical units.
   In this case, you **MUST** open `python/tests/test_units_sync.py` and append your new entity to the `IGNORE_LIST` dictionary to explicitly allow it to exist without units.
