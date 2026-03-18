# CombAero Python API Reference

This document provides the high-level reference for the `combaero` Python package. For the C++ technical reference, see [API_CPP.md](API_CPP.md).

## Table of Contents
- [Core Thermodynamics](#core-thermodynamics)
- [Combustion & Equilibrium](#combustion--equilibrium)
- [Flow Regimes (Symmetric API)](#flow-regimes-symmetric-api)
- [Advanced Flow Functions](#advanced-flow-functions)
- [Comprehensive Orifice Functions](#comprehensive-orifice-functions)
- [Heat Transfer](#heat-transfer)
- [Acoustics](#acoustics)
- [Network Solver](#network-solver)
- [Geometry & Materials](#geometry--materials)
- [Advanced Thermodynamics](#advanced-thermodynamics)
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

## Advanced Flow Functions

### Inverse Solvers
```python
# Solve for unknown given mass flow rate
A_eff = cb.solve_A_eff_from_mdot(T0=300, P0=2e5, P_back=1e5, mdot_target=0.1, X=air)
P_back = cb.solve_P_back_from_mdot(T0=300, P0=2e5, A_eff=1e-4, mdot_target=0.1, X=air)
P0 = cb.solve_P0_from_mdot(T0=300, P_back=1e5, A_eff=1e-4, mdot_target=0.1, X=air)
```

### Utility Functions
```python
# Critical pressure ratio and Mach number
P_crit_ratio = cb.critical_pressure_ratio(T0=300, P0=2e5, X=air)
Mach = cb.mach_from_pressure_ratio(T0=300, P0=2e5, P=1.5e5, X=air)

# Mass flux and flow analysis
G = cb.mass_flux_isentropic(T0=300, P0=2e5, P=1.5e5, X=air)
L_max = cb.fanno_max_length(T_in=300, P_in=2e5, u_in=50, D=0.05, f=0.02, X=air)
```

### Rocket Nozzle Functions
```python
# Converging-diverging nozzle analysis
sol = cb.nozzle_cd(T0=300, P0=2e5, P_design=1e5, P_amb=101325,
                   A_inlet=1e-3, A_throat=5e-4, A_exit=2e-3,
                   x_throat=0.1, x_exit=0.2, X=air)

# Thrust calculations
thrust = cb.nozzle_thrust_cd(T0=300, P0=2e5, P_design=1e5, P_amb=101325,
                             A_inlet=1e-3, A_throat=5e-4, A_exit=2e-3,
                             x_throat=0.1, x_exit=0.2, X=air)
print(f"Thrust: {thrust.thrust} N")
print(f"Specific impulse: {thrust.specific_impulse} s")
```

---

## Comprehensive Orifice Functions

### Discharge Coefficient Correlations
```python
# Individual correlations
Cd = cb.Cd_sharp_thin_plate(geom, state)
Cd = cb.Cd_thick_plate(geom, state)
Cd = cb.Cd_rounded_entry(geom, state)

# Auto-selection based on geometry
Cd = cb.Cd_orifice(geom, state)
```

### Geometry and Flow Analysis
```python
# Area calculations
area = cb.orifice_area(d=0.01)
area = cb.orifice_area_from_beta(beta=0.5, D=0.02)

# Cd-K conversions
K = cb.orifice_K_from_Cd(Cd=0.65)
Cd = cb.orifice_Cd_from_K(K=0.5)

# Flow state analysis
state = cb.orifice_flow_state(P1=2e5, P2=1e5, rho=1.2, mu=1.8e-5, d=0.01, D=0.02)
Re_d = cb.orifice_Re_d_from_mdot(mdot=0.1, d=0.01, mu=1.8e-5)
```

### Advanced Flow Functions
```python
# Thermodynamic orifice flow
sol = cb.orifice_flow_thermo(T=300, P=2e5, X=air, m_dot=0.1, area=1e-4, Cd=0.65)

# Impedance with flow effects
Z = cb.orifice_impedance_with_flow(mdot=0.1, area=1e-4, Cd=0.65, rho=1.2, c=340)

# Thickness corrections
Cd_corrected = cb.orifice_thickness_correction(Cd=0.65, t_over_d=0.1)
```

### Pressure and Flow Calculations
```python
# Pressure drop
dP = cb.orifice_dP(mdot=0.1, area=1e-4, Cd=0.65, rho=1.2)
dP = cb.orifice_dP_Cd(mdot=0.1, area=1e-4, Cd=0.65, rho=1.2)

# Mass flow
mdot = cb.orifice_mdot(P1=2e5, P2=1e5, area=1e-4, Cd=0.65, rho=1.2)
mdot = cb.orifice_mdot_Cd(P1=2e5, P2=1e5, area=1e-4, Cd=0.65, rho=1.2)

# Velocity and heat transfer
v = cb.orifice_velocity(mdot=0.1, area=1e-4, rho=1.2)
v = cb.orifice_velocity_from_mdot(mdot=0.1, area=1e-4, rho=1.2)
Q = cb.orifice_Q(mdot=0.1, h1=3e5, h2=2.8e5)
```

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

### Basic Usage
```python
from combaero.network import FlowNetwork, NetworkSolver, OrificeElement

graph = FlowNetwork()
graph.add_element(OrificeElement("ori1", "nodeA", "nodeB", Cd=0.6, area=1e-4))

solver = NetworkSolver(graph)
results = solver.solve(timeout=5.0)
```

### Network Nodes
```python
from combaero.network import (
    PressureBoundary, MassFlowBoundary, PlenumNode,
    CombustorNode, MomentumChamberNode, EnergyBoundary
)

# Boundary conditions
inlet = PressureBoundary("inlet", P_total=2e5, T_total=300)
outlet = PressureBoundary("outlet", P_total=1e5, T_total=300)
mass_in = MassFlowBoundary("mass_in", m_dot=0.1, T_total=400, Y=air)

# Internal nodes
junction = PlenumNode("junction", P=1.5e5, T=350, Y=air)
combustor = CombustorNode("combustor", P=2e5, T=1800, Y=products)
momentum = MomentumChamberNode("momentum", P=2e5, T=1800, Y=products)
energy = EnergyBoundary("energy", Q=50000)  # Heat addition [W]
```

### Network Elements
```python
from combaero.network import (
    OrificeElement, PipeElement, EffectiveAreaConnectionElement,
    LosslessConnectionElement, AreaDischargeCoefficientConnectionElement
)

# Flow elements
orifice = OrificeElement("orifice", "node1", "node2", Cd=0.65, area=1e-4, regime="compressible")
pipe = PipeElement("pipe", "node2", "node3", length=2.0, diameter=0.05, roughness=1e-4, regime="compressible_fanno")

# Connection elements
effective_area = EffectiveAreaConnectionElement("ea", "node3", "node4", area=2e-4)
lossless = LosslessConnectionElement("lossless", "node4", "node5")
area_cd = AreaDischargeCoefficientConnectionElement("area_cd", "node5", "node6", area=1e-4, Cd=0.7)
```

### Combustion Integration
```python
from combaero.network.combustion import (
    combustion_from_streams, combustion_from_phi, mix_streams, stoichiometric_products
)

# Create combustion from streams
fuel_stream = cb.Stream(m_dot=0.01, T=300, P_total=2e5, Y=fuel_y)
oxidizer_stream = cb.Stream(m_dot=0.2, T=300, P_total=2e5, Y=air_y)

result = combustion_from_streams([fuel_stream, oxidizer_stream], phi=1.0)
print(f"Products: {result.products.Y}")
print(f"Adiabatic T: {result.T_ad} K")

# Or from equivalence ratio
result = combustion_from_phi(phi=0.8, T=300, P=2e5, fuel="CH4", oxidizer="air")
```

### MixtureState for Network Data
```python
from combaero.network import MixtureState

# Constructor: MixtureState(P, P_total, T, T_total, m_dot, Y)
node_state = MixtureState(1e5, 1.1e5, 300, 300, 1.5, air_vec)

print(f"Static pressure: {node_state.P} Pa")
print(f"Total pressure: {node_state.P_total} Pa")
print(f"Mass flow: {node_state.m_dot} kg/s")
print(f"Static density: {node_state.density()} kg/m³")
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

## Geometry & Materials

### Geometry Classes
```python
# Basic geometry
tube = cb.Tube(L=1.0, D=0.1)
annulus = cb.Annulus(R_outer=0.1, R_inner=0.05)

# Can-annular geometry
can_geo = cb.CanAnnularFlowGeometry(
    N_cans=12, R_can=0.05, L_can=0.3,
    R_annulus_inner=0.06, R_annulus_outer=0.1
)
```

### Material Properties
```python
# List available materials
materials = cb.list_materials()

# Thermal conductivity [W/(m·K)]
k_al = cb.k_aluminum_6061
k_ss = cb.k_stainless_steel_316
k_inconel = cb.k_inconel718
k_haynes = cb.k_haynes230
k_tbc = cb.k_tbc_ysz
```

### Advanced Geometry Functions
```python
# Annular areas
area_ann = cb.annular_area(R_outer=0.1, R_inner=0.05)

# Residence times
t_res = cb.residence_time(V=0.01, Q=0.1)
t_res_ann = cb.residence_time_annulus(V=0.01, Q=0.1, R_outer=0.1, R_inner=0.05)
t_res_can = cb.residence_time_can_annular(V=0.01, Q=0.1, N_cans=12, R_can=0.05)
```

### Pin Fin Heat Transfer
```python
# Pin fin correlations
Nu_pin = cb.pin_fin_nusselt(Re=1e4, Pr=0.7, height=0.01, diameter=0.001)
f_pin = cb.pin_fin_friction(Re=1e4, height=0.01, diameter=0.001, spacing=0.02)
```

---

## Advanced Thermodynamics

### Complete State Functions
```python
# Air properties bundle
air = cb.air_properties(T=300, P=101325, humidity=0.0)
print(f"Density: {air.density} kg/m³")
print(f"Viscosity: {air.viscosity} Pa·s")
print(f"Thermal conductivity: {air.thermal_conductivity} W/(m·K)")

# Complete thermodynamic state
thermo = cb.thermo_state(T=300, P=101325, X=air)
complete = cb.complete_state(T=300, P=101325, X=air)
```

### Derivatives
```python
# Temperature derivatives
dh_dT = cb.dh_dT(T=300, X=air)
ds_dT = cb.ds_dT(T=300, X=air)
dcp_dT = cb.dcp_dT(T=300, X=air)
dg_dT = cb.dg_over_RT_dT(T=300, X=air)
```

### Advanced Inverse Solvers
```python
# Temperature from various properties
T_from_h = cb.calc_T_from_h(h_target=3e5, X=air)
T_from_s = cb.calc_T_from_s(s_target=7000, P=101325, X=air)
T_from_cp = cb.calc_T_from_cp(cp_target=1000, X=air)
T_from_u = cb.calc_T_from_u(u_target=2.2e5, X=air)
```

### Flow Analysis
```python
# Dimensionless numbers
Re = cb.reynolds(rho=1.2, v=10, D=0.05, mu=1.8e-5)
Pe = cb.peclet(rho=1.2, v=10, D=0.05, cp=1000, k=0.026)
KE = cb.kinetic_energy(rho=1.2, v=10)
```

---

## Psychrometrics (Humid Air)

```python
air = cb.HumidAir()
air.set_TP_RH(300.0, 101325.0, 0.5)

print(f"Humidity ratio: {air.humidity_ratio} kg/kg")
print(f"Dewpoint: {air.dewpoint} K")
print(f"Underlying state P: {air.state.P} Pa")
```
