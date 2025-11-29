# CombAero Python bindings

This package provides Python bindings for the core CombAero C++ library:

- Thermodynamic and transport properties for gas mixtures
- Combustion helpers (oxygen requirement, complete combustion)
- Equivalence ratio and Bilger mixture-fraction utilities
- Humid air properties and standard dry-air composition
- Species common-name mapping

The bindings are implemented with [pybind11](https://pybind11.readthedocs.io/) and built via [scikit-build-core](https://scikit-build-core.readthedocs.io/).

## Installation (from source)

The recommended workflow is from the repository root; see the top-level README for details. In short:

```bash
python -m build --wheel
python -m pip install dist/combaero-0.0.1-*.whl
```

## Usage

Example using NumPy arrays and the high-level `combaero` package:

```python
import numpy as np
import combaero as ca

# Get species index for N2 and O2
i_N2 = ca.species_index_from_name("N2")
i_O2 = ca.species_index_from_name("O2")

# Create air composition (79% N2, 21% O2)
X = np.zeros(ca.num_species(), dtype=float)
X[i_N2] = 0.79
X[i_O2] = 0.21

T = 300.0   # K
P = 101325.0  # Pa

print(f"cp = {ca.cp(T, X):.2f} J/(mol·K)")
print(f"density = {ca.density(T, P, X):.4f} kg/m³")
print(f"viscosity = {ca.viscosity(T, P, X):.2e} Pa·s")
```
