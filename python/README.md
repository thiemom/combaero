# CombAero Python bindings

This package provides Python bindings for a subset of the CombAero C++ library:

- Thermodynamic properties for gas mixtures
- Selected combustion utilities
- Humid air properties

The bindings are implemented with [pybind11](https://pybind11.readthedocs.io/) and built via [scikit-build-core](https://scikit-build-core.readthedocs.io/).

## Installation (from source)

From the `python/` subdirectory:

```bash
python -m pip install -U pip build
python -m build
```

This will produce wheels and an sdist under `python/dist/`. To install a built wheel:

```bash
python -m pip install dist/combaero-*.whl
```

## Usage

Example using NumPy arrays and the high-level `combaero` package:

```python
import numpy as np
import combaero as ca

T = 300.0
X = np.array([1.0, 0.0, ...])  # your species vector

print("h =", ca.mixture_h(T, X))
print("T_ad =", ca.adiabatic_T_wgs(T, X))
```
