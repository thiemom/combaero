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

The `State` class uses Pythonic properties for all attribute access:

```python
import numpy as np
import combaero as ca

# Create a State - direct property assignment
s = ca.State()
s.T = 300.0
s.P = 101325.0
s.X = X_air

# Or use fluent setters for chaining
s = ca.State()
s.set_T(300.0).set_P(101325.0).set_X(X_air)

# Access thermodynamic properties (all are properties, not methods)
print(f"cp = {s.cp:.2f} J/(mol·K)")
print(f"h = {s.h:.0f} J/mol")
print(f"rho = {s.rho:.4f} kg/m³")
print(f"gamma = {s.gamma:.3f}")
print(f"speed of sound = {s.a:.1f} m/s")

# Access transport properties
print(f"viscosity = {s.mu:.2e} Pa·s")
print(f"thermal conductivity = {s.k:.4f} W/(m·K)")
print(f"Prandtl = {s.Pr:.3f}")
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
print(f"Flame density: {burned.rho:.4f} kg/m³")

# WGS equilibrium (isothermal or adiabatic)
eq = ca.wgs_equilibrium_adiabatic(T=1500.0, X=burned.X)
print(f"Equilibrium T: {eq.T:.0f} K")
```

### Reforming Equilibrium (Rich Combustion)

For rich combustion products with unburned hydrocarbons, use reforming equilibrium
to compute the steam reforming + water-gas shift equilibrium:

```python
import numpy as np
import combaero as ca

# Rich combustion products (with unburned CH4, C2H6, etc.)
X = np.zeros(ca.num_species(), dtype=float)
X[ca.species_index_from_name("CH4")] = 0.02    # Unburned methane
X[ca.species_index_from_name("C2H6")] = 0.003  # Unburned ethane
X[ca.species_index_from_name("H2O")] = 0.20
X[ca.species_index_from_name("CO2")] = 0.08
X[ca.species_index_from_name("N2")] = 0.70

# General reforming equilibrium (handles all hydrocarbons)
# CnHm + n*H2O <-> n*CO + (n + m/2)*H2
# CO + H2O <-> CO2 + H2 (WGS)
eq = ca.reforming_equilibrium_adiabatic(T=2200.0, X=X)
print(f"Equilibrium T: {eq.T:.0f} K")
print(f"CO: {eq.X[ca.species_index_from_name('CO')]*100:.2f}%")
print(f"H2: {eq.X[ca.species_index_from_name('H2')]*100:.2f}%")

# Isothermal version (fixed temperature)
eq_iso = ca.reforming_equilibrium(T=2000.0, X=X)

# SMR+WGS (CH4 only, for backward compatibility)
eq_smr = ca.smr_wgs_equilibrium_adiabatic(T=2200.0, X=X)
```

**Available equilibrium functions:**
- `wgs_equilibrium(T, X, P)` - Isothermal WGS (CO + H2O ⇌ CO2 + H2)
- `wgs_equilibrium_adiabatic(T, X, P)` - Adiabatic WGS
- `smr_wgs_equilibrium(T, X, P)` - Isothermal SMR+WGS (CH4 only)
- `smr_wgs_equilibrium_adiabatic(T, X, P)` - Adiabatic SMR+WGS (CH4 only)
- `reforming_equilibrium(T, X, P)` - Isothermal reforming (all hydrocarbons)
- `reforming_equilibrium_adiabatic(T, X, P)` - Adiabatic reforming (all hydrocarbons)
- `combustion_equilibrium(T, X, P)` - **One-step**: combustion + reforming + WGS

**Tip:** Use `combustion_equilibrium` when starting from an unburned fuel+air mixture.
It combines complete combustion with reforming equilibrium in one call:

```python
# Unburned fuel + air -> equilibrium products in one step
X_unburned = ...  # CH4 + O2 + N2
result = ca.combustion_equilibrium(T=300.0, X=X_unburned)
print(f"Equilibrium T: {result.T:.0f} K")
```

### Stream Mixing

Mix multiple streams with mass and enthalpy balance:

```python
import numpy as np
import combaero as ca

# Create streams - use property assignment (Pythonic)
air = ca.Stream()
air.T = 400.0
air.P = 101325.0
air.X = X_air
air.mdot = 10.0

fuel = ca.Stream()
fuel.T = 300.0
fuel.P = 101325.0
fuel.X = X_fuel
fuel.mdot = 0.5

# Or use fluent setters for one-liners
fuel.set_T(300.0).set_P(101325.0).set_X(X_fuel).set_mdot(0.5)

# Mix streams (uses minimum inlet pressure by default)
mixed = ca.mix([air, fuel])
print(f"Mixed T: {mixed.T:.1f} K")
print(f"Mixed mdot: {mixed.mdot:.2f} kg/s")

# Or specify output pressure explicitly
mixed2 = ca.mix([air, fuel], P_out=150000.0)
```
