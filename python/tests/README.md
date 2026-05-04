# CombAero Python Tests

This directory contains the `pytest` suite for the `combaero` Python bindings and API.

Run locally from the repo root:

```bash
uv run pytest python/tests/ -v
```

CI runs these tests across Linux, macOS, and Windows on every push.

## Unit Synchronization (`test_units_sync.py`)

Every calculation, function, or class property exposed to the end-user **must** have a corresponding unit definition in `include/units_data.h`. `test_units_sync.py` enforces this via Python introspection on `combaero.__all__`.

**If this test fails after adding an API symbol:**

1. **Mathematical calculation or property** — update `include/units_data.h` with the correct unit string and rebuild (`uv pip install -e .`). Do not bypass the test.

2. **Structural module, utility, or enum without physical units** — add the new symbol to the `IGNORE_LIST` in `test_units_sync.py` with a short explanation.
