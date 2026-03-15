# CombAero Python API Reference

This document provides the high-level reference for the `combaero` Python package. For the C++ technical reference, see [API_CPP.md](API_CPP.md).

## Table of Contents
- [Core Thermodynamics](#core-thermodynamics)
- [Combustion & Equilibrium](#combustion--equilibrium)
- [Flow Regimes (Symmetric API)](#flow-regimes-symmetric-api)
- [Heat Transfer](#heat-transfer)
- [Acoustics](#acoustics)
- [Network Solver](#network-solver)
- [Psychrometrics (Humid Air)](#psychrometrics-humid-air)

---

## Core Thermodynamics

### State Object
The `State` object is the central data structure, managing T, P, and composition.

```python
import combaero as cb

state = cb.State()
state.TPX = (300.0, 101325.0, "O2:0.21, N2:0.79")
# OR setting mass fractions directly
state.Y = y_vector

print(f"Mole fractions: {state.X}")
print(f"Mass fractions: {state.Y}")

print(f"Density: {state.rho} kg/m³")
print(f"Enthalpy: {state.h} J/mol")
print(f"Viscosity: {state.mu} Pa·s")
```

### MixtureState (Network Data)
The `MixtureState` struct is used for network node properties and carries both static and stagnation conditions.

```python
# Constructor: MixtureState(P, P_total, T, T_total, m_dot, Y)
node_state = cb.MixtureState(1e5, 1.1e5, 300, 300, 1.5, y_vec)

print(node_state.P)       # Static [Pa]
print(node_state.m_dot)   # Mass flow [kg/s]
print(node_state.density()) # Static density [kg/m³]
```

### Species Lookup
```python
# Metadata
num = cb.num_species()
name = cb.species_name(0)
idx = cb.species_index_from_name("H2O")
mw = cb.species_molar_mass_from_name("CH4")
```

---

## Combustion & Equilibrium

### Complete Combustion
```python
# ADIABATIC flame temperature and products
# smooth=True ensures a continuous Jacobian for root-finding near Phi=0
burned = cb.complete_combustion(state, smooth=True)
print(f"Adiabatic T: {burned.T} K")

# Isothermal products
products = cb.complete_combustion_isothermal(state, smooth=True)
```

### Chemical Equilibrium
```python
# T, P constant equilibrium
eq_state = cb.combustion_equilibrium(state, smooth=True)

# Water-gas shift
wgs = cb.wgs_equilibrium(state)
```

---

## Flow Regimes (Symmetric API)

CombAero provides symmetric submodules for incompressible and compressible flow.

```python
from combaero import incompressible as flow   # OR: from combaero import compressible as flow

sol = flow.pipe_flow(T=400, P=2e5, X=air, u=10.0, L=2.0, D=0.05, f=0.02)
print(f"Pressure drop: {sol.dP} Pa")
print(f"Outlet Mach: {sol.M}") # nan for incompressible
```

### FlowSolution Fields
- `mdot`: Mass flow rate [kg/s]
- `v`: bulk velocity [m/s]
- `dP`: Pressure drop [Pa]
- `M`: Mach number
- `choked`: Boolean flag

---

## Heat Transfer

### Correlations
```python
Nu = cb.nusselt_gnielinski(Re=1e5, Pr=0.7)
h = cb.htc_from_nusselt(Nu, k=0.026, L=0.05)

# Multi-layer wall temperature profile
temps, q = cb.wall_temperature_profile(T_hot=1200, T_cold=300, h_hot=200, h_cold=20, t_over_k=[0.01/50, 0.05/0.5])
```

---

## Acoustics

```python
# Acoustic properties bundle
props = cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=1.0)
print(f"SPL: {props.spl} dB")

# Cavity resonators
f_helm = cb.helmholtz_frequency(V=0.001, A_neck=1e-4, L_neck=0.01, c=340)

# Duct modes
tube = cb.Tube(L=1.0, D=0.1)
modes = cb.tube_axial_modes(tube, c=340, bc1=cb.BoundaryCondition.Closed, bc2=cb.BoundaryCondition.Open)
```

---

## Network Solver

```python
from combaero.network import FlowNetwork, NetworkSolver, OrificeElement

graph = FlowNetwork()
graph.add_element(OrificeElement("ori1", "nodeA", "nodeB", Cd=0.6, area=1e-4))

solver = NetworkSolver(graph)
results = solver.solve(timeout=5.0)
```

### Advanced (f, J) Interface
High-accuracy analytical Jacobians are available for performance-critical solver loops.
```python
# Exact residuals and derivatives wrt (P_tot, P_static, T, Y)
res_ori = cb.orifice_residuals_and_jacobian(m_dot, P_tot, P_stat, T, Y, P_down, Cd, area)
res_comb = cb.combustor_residuals_and_jacobian(m_dot, P_tot, P_stat, T, Y, Q_comb, method, smooth)
res_plen = cb.plenum_residuals_and_jacobian(m_dot_vec, P_target, T_target, Y_target)
```

Combustor results return mapping specific to `(m_dot, P_total, T, Y)` state vectors.

---

## Psychrometrics (Humid Air)

```python
air = cb.HumidAir()
air.set_TP_RH(300.0, 101325.0, 0.5)

print(f"Humidity ratio: {air.humidity_ratio} kg/kg")
print(f"Dewpoint: {air.dewpoint} K")
print(f"Underlying state P: {air.state.P} Pa")
```
