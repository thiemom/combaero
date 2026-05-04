# CombAero Python Bindings

Python bindings for the CombAero C++ library, built with [pybind11](https://pybind11.readthedocs.io/) and [scikit-build-core](https://scikit-build-core.readthedocs.io/).

## Installation

```bash
uv pip install combaero
```

For an editable dev install from the repo root:

```bash
uv pip install -e .
```

## API Reference

Full usage, examples, and function signatures:

- **[API Reference (Python)](../docs/API_PYTHON.md)** — thermodynamics, combustion, flow, cooling, acoustics
- **[Units Guide](../docs/UNITS.md)** — SI units for every output quantity

## Running Tests

```bash
uv run pytest python/tests/
```
