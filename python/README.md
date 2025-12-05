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

### State-based API

The `State` class provides property getters and setters with chaining:

```python
import numpy as np
import combaero as ca

# Create a State with setters
s = ca.State()
s.set_T(300.0).set_P(101325.0).set_X(X_air)

# Access properties directly
print(f"cp = {s.cp():.2f} J/(mol·K)")
print(f"h = {s.h():.0f} J/mol")
print(f"rho = {s.rho():.4f} kg/m³")
print(f"gamma = {s.gamma():.3f}")
print(f"speed of sound = {s.a():.1f} m/s")
```

### Combustion and Equilibrium

```python
import numpy as np
import combaero as ca

# Create a CH4 + air mixture
X = np.zeros(ca.num_species(), dtype=float)
X[ca.species_index_from_name("CH4")] = 0.095
X[ca.species_index_from_name("O2")] = 0.19
X[ca.species_index_from_name("N2")] = 0.715

# Adiabatic complete combustion
burned = ca.complete_combustion(T=300.0, X=X)
print(f"Adiabatic flame T: {burned.T:.0f} K")
print(f"Flame density: {burned.rho():.4f} kg/m³")

# WGS equilibrium (isothermal or adiabatic)
eq = ca.wgs_equilibrium_adiabatic(T=1500.0, X=burned.X)
print(f"Equilibrium T: {eq.T:.0f} K")
```

### Stream Mixing

Mix multiple streams with mass and enthalpy balance:

```python
import numpy as np
import combaero as ca

# Create streams
air = ca.Stream()
air.set_T(400.0).set_P(101325.0).set_X(X_air).set_mdot(10.0)

fuel = ca.Stream()
fuel.set_T(300.0).set_P(101325.0).set_X(X_fuel).set_mdot(0.5)

# Mix streams (uses minimum inlet pressure by default)
mixed = ca.mix([air, fuel])
print(f"Mixed T: {mixed.T():.1f} K")
print(f"Mixed mdot: {mixed.mdot:.2f} kg/s")

# Or specify output pressure explicitly
mixed2 = ca.mix([air, fuel], P_out=150000.0)
```
