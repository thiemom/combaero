# CombAero Python API Reference

This document provides the high-level reference for the `combaero` Python package. For the C++ technical reference, see [API_CPP.md](API_CPP.md).

## Table of Contents
- [Core Thermodynamics](#core-thermodynamics)
- [Combustion & Equilibrium](#combustion--equilibrium)
- [Flow Regimes (Symmetric API)](#flow-regimes-symmetric-api)
- [Advanced Flow Functions](#advanced-flow-functions)
- [Comprehensive Orifice Functions](#comprehensive-orifice-functions)
- [Heat Transfer](#heat-transfer)
  - [Correlations](#correlations)
  - [Network Heat Transfer](#network-heat-transfer)
- [Acoustics](#acoustics)
- [Network Solver](#network-solver)
- [NetworkRunner (GUI-JSON Programmatic Driver)](#networkrunner-gui-json-programmatic-driver)
- [Geometry & Materials](#geometry--materials)
- [Advanced Thermodynamics](#advanced-thermodynamics)
- [Psychrometrics (Humid Air)](#psychrometrics-humid-air)

---

## Core Thermodynamics

### State Object
The `State` object is the central data structure, managing T [K], P [Pa], and composition.

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
**Units**: P [Pa], T [K], m_dot [kg/s].

```python
# Constructor: MixtureState(P, P_total, T, T_total, m_dot, Y)
node_state = cb.MixtureState(1e5, 1.1e5, 300, 300, 1.5, y_vec)

print(node_state.P)       # Static [Pa]
print(node_state.m_dot)   # Mass flow [kg/s]
print(node_state.density()) # Static density [kg/m³]
```

### Species Module
The `combaero.species` module provides high-level helpers for composition vectors.

```python
import combaero as cb

# Get standard compositions (returns numpy array)
X_air = cb.species.dry_air()        # Mole fractions
Y_air = cb.mole_to_mass(X_air)   # Mass fractions (standard for network)

# Humid air at T [K], P [Pa], RH [0-1]
Y_humid = cb.species.humid_air_mass(300, 101325, 0.5)

# Pure species or from mapping
Y_ch4 = cb.species.pure_species("CH4")
Y_mix = cb.species.from_mapping({"O2": 0.20, "N2": 78, "AR": 0.02})

# Explicit conversion
Y = cb.species.to_mass(X_air)
X = cb.species.to_mole(Y_air)
```

### Species Lookup (Legacy)
```python
# Low-level metadata from core
num = cb.num_species()
name = cb.species_name(0)
idx = cb.species_index_from_name("H2O")
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

sol = flow.channel_flow(T=400, P=2e5, X=air, u=10.0, L=2.0, D=0.05, f=0.02)
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

#### Plenum-to-plenum holes: McGreehan and Schotsch (1988)

The correlations above are ISO 5167 pipe metering -- a single hole IN a pipe
run. For a cooling transfer hole or a jet-plate hole discharging
plenum-to-plenum, use the composite chain of McGreehan and Schotsch (1988),
ASME J. Turbomachinery 110(2), 213-217, which covers Reynolds number, inlet
corner radius, orifice length and inlet crossflow:

```python
# A bare, sharp-edged jet plate, plate thickness = hole diameter
cb.mcgreehan_schotsch_1988_cd(Re=1e4, r_over_d=0.0, L_over_d=1.0)   # 0.785

# A radiused, long cooling hole
cb.mcgreehan_schotsch_1988_cd(Re=3e4, r_over_d=0.1, L_over_d=3.0)   # 0.889

# Fed from a duct rather than a plenum: supply-side crossflow
cb.mcgreehan_schotsch_1988_cd(Re=3e4, r_over_d=0.0, L_over_d=1.0,
                              U1_over_Vi=0.5)                       # 0.770
```

`U1_over_Vi` is the **inlet (approach, supply-side)** tangential velocity over
the ideal through-flow velocity, and defaults to 0 for a plenum-fed plate. It
is *not* a discharge-side crossflow ratio -- do not pass Florschuetz's
`Gc/Gj`, which is spent air in the impingement channel on the far face of the
plate. `V_i` must be built from static inlet conditions, not from a total
pressure that already contains the tangential velocity head.

Two behaviours that look like bugs and are not:

- **Not monotonic in `U1_over_Vi`.** `Cd` rises up to ~5.7% above its
  zero-crossflow value near `U1_over_Vi ~ 0.09` before falling away. This is
  drawn in the source's Fig. 4 and supported by Rohde (NASA TN D-5467), whose
  result 2 is that slanting an orifice into the flow increases `Cd`.
- **Constant below `Re = 1e4`.** That is the correlation's stated floor, and
  below it the underlying equation diverges (it reaches `Cd = 1.0` at
  `Re = 904`), so the value is held at the floor rather than extrapolated.

When the plate's zero-crossflow `Cd` is already known -- measured, or a
literature value such as Florschuetz's per-configuration Table 1 -- apply only
the crossflow correction:

```python
cb.mcgreehan_schotsch_1988_crossflow_cd(cd_base=0.79, U1_over_Vi=0.5)  # 0.794
cb.mcgreehan_schotsch_1988_crossflow_cd(cd_base=0.79, U1_over_Vi=2.0)  # 0.541
```

Note the first: at a 0.79 baseline, `U1/Vi = 0.5` still sits inside the rise,
because a higher baseline pushes the peak to larger `U1/Vi` (`R_v` carries a
`(Cd/0.6)^-3` factor). The decay is well established by `U1/Vi = 2`.

This is how the source uses Eq. (17) in its own validation. Its Figs. 5 and 6
anchor to "a set baseline point at `U1/Vi = 0`" taken from Rohde's
measurements (0.64, 0.73, 0.88), not to the chain's prediction for the same
geometry -- the chain runs 3-12% higher, which is the paper's own remark that
Rohde's "basic values are lower".

**Measured agreement.** Against Rohde's data as replotted in Fig. 6
(`t/d = 0.51`, baseline 0.64, 11 points), Eq. (17) scores bias `+5.95%`,
RMS `7.54%` -- reading high, and increasingly so with crossflow (`+3%` below
`U1/Vi = 0.4`, `+16%` at 1.41). Scored by
`validation/cooling/orifice_runner.py`.

Provenance, the ten checks behind every constant, and the two errata found in
the paper are in
`validation/cooling/extractions/orifice_discharge_coefficient.md`.

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

### Network Heat Transfer

#### ConvectiveSurface
Defines convective heat transfer properties on network elements. Supports different channel models.

```python
from combaero.heat_transfer import ConvectiveSurface, SmoothModel, RibbedModel, DimpledModel, PinFinModel, ImpingementModel

# Smooth channel with Gnielinski correlation (default)
surface = ConvectiveSurface(
    area=np.pi * 0.04 * 1.0,  # Surface area [m²]
    model=SmoothModel(correlation="gnielinski")
)

# Ribbed surface with geometry parameters
ribbed = ConvectiveSurface(
    area=2.5,
    model=RibbedModel(
        e_D=0.05,      # Rib height / hydraulic diameter
        p_e=10.0,     # Pitch / height
        w_e=0.5,      # Rib width / height
        correlation="gnielinski"
    )
)

# Pin fin array
pin_fin = ConvectiveSurface(
    area=3.0,
    model=PinFinModel(
        L_H=1.0,       # Fin height / hydraulic diameter
        S_H=2.0,       # Spanwise spacing / height
        S_L=2.0,       # Streamwise spacing / height
        t_D=0.1,       # Fin thickness / diameter
        correlation="gnielinski"
    )
)
```

#### WallConnection
Thermal coupling between two network elements through a shared wall.

```python
from combaero.heat_transfer import WallConnection

# Couple hot and cold channels through a wall
wall = WallConnection(
    id="coupling_wall",
    element_a="hot_channel",      # First element ID
    element_b="cold_channel",     # Second element ID
    wall_thickness=0.002,      # Wall thickness [m]
    wall_conductivity=25.0,    # Wall thermal conductivity [W/(m·K)]
    contact_area=None          # Optional: override surface area
)

# Add to network
network.add_wall(wall)
```

#### Element Integration
Elements with convective surfaces support heat transfer calculations.

```python
from combaero.network import ChannelElement

# Channel with convective heat transfer
channel = ChannelElement(
    id="heated_channel",
    from_node="inlet",
    to_node="outlet",
    diameter=0.04,
    length=1.0,
    roughness=0.0,
    surface=surface  # ConvectiveSurface instance
)

# Get heat transfer coefficient and temperature
h, T = channel.htc_and_T(state)  # Returns (htc [W/(m²·K)], T [K])
```

#### Thermal Coupling Toggle
Enable/disable thermal coupling globally for debugging or fast solves.

```python
network = FlowNetwork()
network.thermal_coupling_enabled = True   # Enable wall heat transfer
network.thermal_coupling_enabled = False  # Disable (faster, no heat transfer)
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

## Programmatic Network Design

For automated design and agents, `combaero.network` provides a pure-Python interface to build and solve networks without the GUI.

```python
from combaero.network import FlowNetwork, NetworkSolver, OrificeElement, PressureBoundary, PlenumNode
import combaero.species as species

# 1. Initialize Network
net = FlowNetwork()

# 2. Define Boundaries & Nodes
net.add_node(PressureBoundary("inlet", P_total=5e5, T_total=400, Y=species.dry_air_mass()))
net.add_node(PlenumNode("plenum", P=4e5, T=400))
net.add_node(PressureBoundary("outlet", P_total=1e5, T_total=300))

# 3. Add Elements
net.add_element(OrificeElement("feed", "inlet", "plenum", diameter=0.01, Cd=0.65))
net.add_element(OrificeElement("exhaust", "plenum", "outlet", diameter=0.015, Cd=0.62))

# 4. Solve
solver = NetworkSolver(net)
results = solver.solve(method="hybr")

# 5. Access Results
print(f"Plenum Pressure: {results.nodes['plenum'].P} Pa")
print(f"System Mass Flow: {results.elements['feed'].m_dot} kg/s")
```

> [!TIP]
> **Agent Hint**: Use `net.to_dict()` to get a JSON-serializable representation of the network, which can be saved or sent to the GUI.

---

## Network Solver

### What a solve reports

`NetworkSolver.solve` returns the solved unknowns plus a set of `__dunder__`
bookkeeping keys. Five of them say what happened:

| key | type | meaning |
|---|---|---|
| `__success__` | `bool` | converged **and** consistent. Unchanged; keep using it for "can I trust this result". |
| `__converged__` | `bool` | the root finder reached the residual tolerance, whatever the consistency checks then said |
| `__consistent__` | `bool \| None` | the elements' own physical-consistency verdict. **`None` means not checked** -- either the solve never converged, or no element in this network has a verifier. It is not a pass. |
| `__inconsistent_elements__` | `list[str]` | ids of the elements that rejected the solution |
| `__outcome__` | `SolveOutcome` | why it ended, in one machine-readable value |
| `__worst_residuals__` | `list[dict]` | the rows carrying the residual, largest first |

`__success__` is False in two quite different situations, and reading it alone
cannot tell them apart: Newton never got there, or Newton got there and a
junction rejected the root as unphysical. The second reports a *small*
`__final_norm__` next to `success=False`, which looks like a contradiction
until `__converged__` and `__consistent__` are read separately.

```python
from combaero.network import SolveOutcome

result = solver.solve()
if result["__outcome__"] == SolveOutcome.INCONSISTENT:
    print("converged, then rejected by", result["__inconsistent_elements__"])
elif not result["__success__"]:
    worst = result["__worst_residuals__"][0]
    print(f"{result['__outcome__']}: {worst['name']} carries {worst['residual']:.3e}")
```

`SolveOutcome` is a `StrEnum`, so it compares and serialises as a plain
string: `CONVERGED`, `INCONSISTENT`, `NO_PROGRESS`, `RESIDUAL_TOO_LARGE`,
`TIMEOUT`, `ERROR`, `NOT_CONVERGED`, `NO_UNKNOWNS`.

> [!TIP]
> Prefer `__outcome__` over matching on `__message__`. Part of that text comes
> from SciPy and can be reworded without notice; `__outcome__` is classified
> once, by the code that knows the answer.

Note there is no `__residual_norm__`; the norm is `__final_norm__`.


### Basic Usage
```python
from combaero.network import FlowNetwork, NetworkSolver, OrificeElement

graph = FlowNetwork()
graph.add_element(OrificeElement("ori1", "nodeA", "nodeB", Cd=0.6, diameter=0.011284))

solver = NetworkSolver(graph)
# init_strategy options: "default", "analytical_pt_prop", "homotopy",
# "continuation", "outletref_warmstart", "incompressible_warmstart"
# (deprecated).
# "default" auto-upgrades to "analytical_pt_prop" for networks that
# contain a MultiPortChamberBase; pass another strategy or an
# explicit x0 to opt out.
# "outletref_warmstart" solves the outlet-referenced incompressible
# proxy (densities at the downstream static) and warm-starts the
# compressible solve from it directly -- the auto-retry's seed as a
# primary strategy, for networks known to stall the cold path.
# auto_retry (default True): failed cold solves on compressible
# MPCE networks are retried once from an outlet-referenced
# incompressible warm start (densities at the downstream static).
# The primary attempt gets 40% of `timeout`, the retry the rest.
# Stall detection ends doomed phases early: when hybr's best |F|
# plateaus far from tolerance it hands over to the LM fallback, and
# when the LM fallback plateaus too the attempt fails fast so the
# auto-retry gets the budget. Cold attempts only -- seeded solves
# (explicit x0 or a warm-start seed) are the last rung and keep
# their full budget.
results = solver.solve(method="hybr", init_strategy="homotopy")
```

### Network Nodes
```python
from combaero.network import (
    PressureBoundary, MassFlowBoundary, WallNode, PlenumNode,
    CombustorNode, MomentumChamberNode, MomentumBoundary, EnergyBoundary
)

# Boundary conditions
inlet = PressureBoundary("inlet", P_total=2e5, T_total=300)
outlet = PressureBoundary("outlet", P_total=1e5, T_total=300)
mass_in = MassFlowBoundary("mass_in", m_dot=0.1, T_total=400, Y=air)
mom_in = MomentumBoundary("mom_in", P_total=2e5, T_total=300, area=1e-4)

# Closed-end boundary (zero mass flow)
# Typical use: symmetry plane of a ring manifold, or any dead-end branch.
wall = WallNode("wall")

# Internal nodes
junction = PlenumNode("junction", P=1.5e5, T=350, Y=air)
combustor = CombustorNode("combustor", P=2e5, T=1800, Y=products)
momentum = MomentumChamberNode("momentum", P=2e5, T=1800, Y=products, regime="compressible")
energy = EnergyBoundary("energy", Q=50000)  # Heat addition [W]
```

### Network Elements
```python
from combaero.network import (
    OrificeElement, ChannelElement, EffectiveAreaConnectionElement,
    LosslessConnectionElement, DiameterDischargeCoefficientConnectionElement,
    TeeJunctionElement, VortexElement,
    BorderCarnotLossElement,
)
from combaero.network.mpce_element import ConstantKTeeElement, MultiPortChamberElement

# Flow elements
orifice = OrificeElement("orifice", "node1", "node2", Cd=0.65, diameter=0.011284, regime="compressible")
channel = ChannelElement("channel", "node2", "node3", length=2.0, diameter=0.05, roughness=1e-4, regime="compressible")

# Connection elements
effective_area = EffectiveAreaConnectionElement("ea", "node3", "node4", diameter=0.015958)
lossless = LosslessConnectionElement("lossless", "node4", "node5")
area_cd = DiameterDischargeCoefficientConnectionElement("area_cd", "node5", "node6", diameter=0.011284, Cd=0.7)

# Three-port tee junction (Unified0D compressible model, Mynard & Valen-Sendstad 2015)
import math
# Merging tee: two inlets (straight_node, branch_node) -> one outlet (common_node)
tee_merge = TeeJunctionElement(
    "tee", common_node="outlet", straight_node="inlet_s", branch_node="inlet_b",
    theta=math.pi / 2.0, F_C=0.01, psi=1.0, tee_type="merging",
)
# Branching tee: one inlet (common_node) -> two outlets (straight_node, branch_node)
tee_branch = TeeJunctionElement(
    "tee", common_node="inlet", straight_node="outlet_s", branch_node="outlet_b",
    theta=math.pi / 2.0, F_C=0.01, psi=1.0, tee_type="branching",
)
# Solved unknowns: tee.m_dot_com (total), tee.m_dot_branch (branch arm)
# Straight flow is implicit: m_dot_straight = m_dot_com - m_dot_branch

# Momentum-CV junction (PDF spec, supersedes K-closure for n>3 manifolds and
# high-Mach / ejector behaviour; see docs/junction/momentum cv implementation guide.pdf).
# N port-MCNs feed a single junction; each lateral port carries a turning loss.
mc_com = MomentumChamberNode("mc_com", area=0.01)
mc_str = MomentumChamberNode("mc_str", area=0.01)
mc_bra = MomentumChamberNode("mc_bra", area=0.008)

# MultiPortChamberBase is the ABSTRACT BASE since 0.6.0 -- it owns the port
# machinery and no longer carries a residual. Use MultiPortChamberElement (the Mynard
# closure) or ConstantKTeeElement (fixed handbook K).
jct = MultiPortChamberElement(
    "jct",
    inlet_nodes=["mc_com"],
    outlet_nodes=["mc_str", "mc_bra"],
    inlet_angles_deg=[0.0],
    outlet_angles_deg=[0.0, 90.0],   # geometric branch angles per outlet port
    flow_direction="branch",         # 1 supplier + N-1 collectors
    # port_areas (length = N_in + N_out) inherited from connecting channels if not given
)
# Lateral port turning loss (sharp-edged 90 deg branch). Straight ports
# (delta_geom = 0) do not need a loss element.
loss_bra = BorderCarnotLossElement(
    "loss_bra", from_node="mc_bra", to_node="ch_bra_in",
    delta_geom_deg=90.0,
    # area inherited from neighbouring channel if not given
)
# Solved unknowns: jct.P_jct (junction static pressure). Per-port mass flows
# live on the connecting channels / loss elements; the junction reads them via
# the graph. The N+1 residuals = N port total-pressure relations
# (Pt_i - Pt_jct +/- K_i * q_dyn_com) + global sum-of-port-mdots = 0, computed
# in C++ (see mpce_residuals_and_jacobian in API_CPP.md).

# Constant-K junction (the "simplest model" tier): fixed handbook loss
# coefficients instead of the Mynard closure. K is referenced to the
# common-leg dynamic head (Idelchik convention) and does not vary with the
# flow split. K_ports keys are port indices (inlets first, then outlets);
# the common port (single inlet for 'branch', single outlet for 'merge')
# is ignored. Exposed in the GUI as the tee's "Constant K (handbook)"
# junction model.
from combaero.network.mpce_element import ConstantKTeeElement
jct_k = ConstantKTeeElement(
    "jct",
    inlet_nodes=["mc_com"],
    outlet_nodes=["mc_str", "mc_bra"],
    flow_direction="branch",
    K_ports={1: 0.4, 2: 1.0},   # K_straight, K_branch
)

# Rotating cavity / disc-pump pressure rise (Vatistas n-vortex model)
vortex = VortexElement(
    "vortex", "node_in", "node_out",
    r_c=0.02,       # vortex core radius [m]
    r_out=0.10,     # outer evaluation radius [m]
    r_in=0.0,       # inner evaluation radius [m] (0 = on-axis)
    omega_rpm=3000, # shaft speed [rpm]
    n=2.0,          # Vatistas shape parameter (>= 1, default n=2)
)
# Residual: Pt_out - Pt_in - dP_vortex = 0
# dP_vortex = vatistas_delta_p(r_out, ...) - vatistas_delta_p(r_in, ...)
# Gamma is derived from omega so that V_theta(r_c) = omega * r_c exactly (solid-body core).
# Diagnostics include: m_dot, omega_rpm, dP_vortex [Pa], Pt_in, Pt_out.
# If omega_rpm=0 the element is lossless (transparent); r_in=0 omits the on-axis term.

# Supersonic ejector across all operating regimes -- critical (double-choked),
# subcritical droop, and unchoked-primary subsonic jet pump -- in one C1 residual
# system (Huang 1999 entrainment ratio + Kracik & Dvorak 2016 mixing closure;
# see validation/ejector/README.md and OPERATING_REGIMES_DESIGN.md).
from combaero.network.ejector_element import EjectorElement
ejector = EjectorElement(
    "ej1", primary_node="mc_primary", secondary_node="mc_secondary", outlet_node="mc_outlet",
    throat_area=3.14e-5,      # primary nozzle throat area A_t [m^2]
    nozzle_exit_area=1.0e-4,  # primary nozzle exit area A_p1 [m^2], must exceed throat_area
    mixing_area=8.0e-4,       # constant-area mixing section A_3 [m^2], must exceed nozzle_exit_area
    recovery_efficiency=1.0,  # multiplies the lossless mixed stagnation pressure; not fitted to data
)
# Solved unknowns: the three port mass flows + ej1.P_jct, the single owned scalar,
# repurposed as the mixing-plane static pressure P_py. Four residual rows: primary
# C-D nozzle (choked or unchoked), blended entrainment, mass conservation, and a
# regime-dependent outlet/discharge closure. The regime is picked smoothly by two
# smootherstep weights (s_choke on mp/mdot_choked, s_sub on outlet.Pt/P_c*): in
# critical mode the outlet floats free below P_c* and P_jct is the diagnostic
# recovered stagnation; in the subcritical/jet-pump regimes the outlet is pinned
# to the ejector's discharge and P_py floats. Working fluid: combaero's real-fluid
# EOS (combustion/air species only, no halocarbon refrigerants); gamma/R evaluated
# live at the entrained-flow choking plane. The whole 4-row (f, J) is assembled in
# C++ (full analytic Jacobian, no FD fallback). diagnostics() reports the actual
# omega = ms/mp, a critical_mode flag, and (in critical mode) p_c_star_pa.
# verify_solution_consistent(sol) returns True unconditionally now that every
# regime is modeled.
```

### Flow-Area Inference

Components that need a reference flow area infer one from the topology when you
do not supply it. Inference searches the neighbourhood of the attached nodes:
channels first (so networks that resolved before resolve identically), then any
other area-bearing element (`AreaChangeElement`, `TeeJunctionElement`,
`BorderCarnotLossElement`, `PressureLossElement`), then adjacent nodes whose
area you set explicitly (`MomentumChamberNode`, `CombustorNode`). An area that
was itself auto-sized is never used as a source -- inheriting a placeholder
would launder one guess into another.

What happens when nothing can be inferred depends on what the area means:

| Component | Parameter | Unresolvable |
|---|---|---|
| `AreaChangeElement` | `F0`, `F1` | **raises** |
| `ChannelElement` | `diameter` | **raises** |
| `TeeJunctionElement` | `F_C` | **raises** |
| `BorderCarnotLossElement` | `area` | **raises** |
| `MomentumChamberNode` | `area` | nominal default |
| `CombustorNode` | `area` | nominal default |
| `PressureLossElement` | `area` | nominal default, warns if a head-loss correlation reads it |

The first group raises because the area *is* the geometry being modelled -- a
default there invents a contraction or expansion the network does not contain,
and the solver either stalls on it or converges to a confidently wrong answer.
The second group defaults because the area is only a dynamic-head reference and
a network can legitimately never read it.

```python
# Unresolvable -> a named error saying which parameter clears it
AreaChangeElement("expansion", "plenum_a", "plenum_b")
# ValueError: AreaChangeElement 'expansion': cannot determine upstream area F0
# -- no neighbouring channel, area-bearing element or node with a known area
# was found at 'plenum_a'. Set upstream area F0 explicitly ...

# Either fix clears it: set the value, or give a neighbour a known area
AreaChangeElement("expansion", "plenum_a", "plenum_b", F0=0.05, F1=0.08)
MomentumChamberNode("chamber", area=0.15)  # neighbour the area change can read
```

### Tuned Constants in the Junction Closure

`MultiPortChamberElement`'s Mynard closure carries three tuned constants, each named
and documented at its definition in `combaero.network._mynard2010`:

| constant | value | what it is | switch |
|---|---|---|---|
| `MYNARD_ETA_A0`, `MYNARD_ETA_A1` | 0.8, -0.2 | Mynard 2015 Eq 36 energy-transfer factor, CFD-fitted | `eta_scale` |
| `FLOW_RATIO_DAMPING` | 0.02 | regulariser from the Matlab reference, not physics | -- |
| `MultiPortChamberElement.DEFAULT_JOINING_ETRANSFER_ALPHA` | 0.2 | combaero's joining-side correction, fitted by `validation/junction/calibrate_etransfer.py` | pass `joining_etransfer_alpha=0.0` |

Under this repo's policy a tuned constant earns its place only by improving
agreement on the digitised validation data; the on/off tables live on issue
#271. `eta_scale` exists for that measurement:

```python
MultiPortChamberElement(..., eta_scale=1.0)   # default: the faithful Mynard port
MultiPortChamberElement(..., eta_scale=0.0)   # energy-transfer term off, for scoring
junction_loss_coefficient(U, A, theta, eta_scale=0.0)
```

Changing any of these values is a retune and needs a before/after table
before it lands, not a silent edit -- `python/tests/test_junction_tuned_constants.py`
pins the documented values.

### The Junction Soft-Barrier Weight

Distinct from the tuned constants above: this is a **numerical** parameter, not
physics, and it is derived rather than declared.

When `strict=False` and a port flows against its declared direction,
`MultiPortChamberElement` replaces the physics with a continuity residual plus a
one-sided quadratic penalty, `alpha * max(0, -e_i * mdot_i)^2`. The penalty
shares its row with the continuity relation, so it balances against a pressure
error rather than driving the offending flow to zero, and has a fixed point at
`slack* = sqrt(dP / alpha)`. A solve that reaches it parks there.

`alpha` therefore carries `Pa/(kg/s)^2` and cannot be a constant: the weight
needed scales as `1/m_ref^2`, two decades per decade of network size.
`NetworkSolver` derives it before each solve from the reference state it
already computes for seeding:

```python
alpha = P_ref / (BARRIER_SLACK_FRACTION * m_ref) ** 2   # f = 0.005
```

placing the fixed point at 0.5% of the reference mass flow whatever the
network's size. It is frozen for the solve, so the residual and Jacobian are
unchanged in form.

| name | where | meaning |
|---|---|---|
| `BARRIER_SLACK_FRACTION` | `combaero.network.mpce_element` | where the fixed point is placed, as a fraction of `m_ref` |
| `DEFAULT_SOFT_PENALTY_ALPHA` | same | fallback when no solver has supplied a weight, or the reference state is degenerate |
| `scaled_penalty_alpha(P_ref, m_ref)` | same | the derivation, exposed for testing |
| `MultiPortChamberElement.effective_penalty_alpha()` | element | the weight actually used, after precedence |

Precedence: an explicitly set `soft_penalty_alpha` always wins, then the
solver-supplied scale-aware weight, then the fallback. So the tuning knob keeps
working:

```python
element.soft_penalty_alpha = 5.0e7   # explicit: overrides the derived weight
```

Raising `alpha` shrinks the fixed point and never destabilises the solve --
the response is monotone and saturates -- so when in doubt, larger.

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
res_chan = cb.channel_residuals_and_jacobian(m_dot, P_tot, P_stat, T, Y, P_down, L, D, roughness, correlation)
res_comb = cb.combustor_residuals_and_jacobians(m_dot, P_tot, P_stat, T, Y, Q_comb, method, smooth, pressure_loss_func)
res_plen = cb.plenum_residuals_and_jacobian(m_dot_vec, P_target, T_target, Y_target)
```

Combustor results return mapping specific to `(m_dot, P_total, T, Y)` state vectors.

### Tee Junction Interface

Two solver-facing APIs are available:

**Unified0D compressible model** (current, used by `TeeJunctionElement`): stagnation-pressure
residuals with mass-flow continuity, consistent with compressible duct elements. Supports
arbitrary per-branch gamma and R_gas.

**Legacy Bassett 2001 model**: empirical incompressible K tables retained for direct
low-level access; not used by the network solver.

#### Unified0D compressible interface

`BranchInput` holds per-branch thermodynamic state. `CompressibleTeeResult` holds the
two residuals and the full Jacobian.

```python
import combaero._core as _core
import math

# Build per-branch state (P_static [Pa], Pt [Pa], T [K], m_dot [kg/s],
#                         A [m^2], theta [rad], gamma_eff [-], R_gas [J/kg/K])
com = _core.BranchInput(P_static=1.95e5, Pt=2.0e5, T=400.0, m_dot=0.5,
                         A=0.01, theta=0.0, gamma_eff=1.4, R_gas=287.0)
str_b = _core.BranchInput(P_static=1.88e5, Pt=1.95e5, T=400.0, m_dot=0.3,
                           A=0.01, theta=0.0, gamma_eff=1.4, R_gas=287.0)
bra = _core.BranchInput(P_static=1.82e5, Pt=1.90e5, T=400.0, m_dot=0.2,
                         A=0.01, theta=math.pi/2, gamma_eff=1.4, R_gas=287.0)

# Branching: common=supplier, straight+branch=collectors
res = _core.compressible_branching_tee_rj(com=com, str=str_b, bra=bra)
# res.R_0, res.R_1           -- residuals [Pa]
# res.dR0_dPt_com, res.dR0_dPt_str, ...  -- full Jacobian fields

# Merging: straight+branch=suppliers, common=collector
res = _core.compressible_merging_tee_rj(com=com, str=str_b, bra=bra)
```

#### Legacy Bassett 2001 interface

```python
import combaero._core as _core
import math

Y = [0.0] * 15
Y[0] = 0.767  # N2
Y[1] = 0.233  # O2

# Check input validity (non-throwing)
status = _core.tee_check_inputs(q=0.5, psi=1.0, theta=math.pi/2)
# status.valid, status.q_in_range, status.psi_valid, status.theta_valid

# Raw K-coefficient functions (Bassett 2001 Table 2 + Eq 33/34 angle corrections).
# K2 and K5 take only q (psi/theta-independent per Eq 15); the rest take (q, psi, theta).
# Provided for diagnostics and validation against measured data.
k1  = _core.tee_K1(q, psi, theta)
k2  = _core.tee_K2(q)
k3  = _core.tee_K3(q, psi, theta)
k4  = _core.tee_K4(q, psi, theta)
k5  = _core.tee_K5(q)               # straight arm, separating tee
k6  = _core.tee_K6(q, psi, theta)   # branch arm, separating tee
k7  = _core.tee_K7(q, psi, theta)
k8  = _core.tee_K8(q, psi, theta)
k9  = _core.tee_K9(q, psi, theta)
k10 = _core.tee_K10(q, psi, theta)
k11 = _core.tee_K11(q, psi, theta)  # straight arm, joining tee
k12 = _core.tee_K12(q, psi, theta)  # branch arm, joining tee

# Blended K functions (always finite, smooth at q=0 and across topology reversal)
ks = _core.merging_tee_K_straight(q, psi, theta)
kb = _core.merging_tee_K_branch(q, psi, theta)

# Legacy solver residuals and full Jacobians
res = _core.merging_tee_residuals_and_jacobian(
    m_dot_com=0.5, m_dot_branch=0.2,
    dP0_straight=150.0, dP0_branch=200.0,
    P_static_com=2e5, T_com=600.0, Y_com=Y,
    theta=math.pi/2, psi=1.0, F_C=0.01,
)
# res.R_straight, res.R_branch
# res.dR_straight_d_mdot_com, res.dR_straight_d_mdot_branch
# res.dR_straight_dP_static_com, res.dR_straight_dT_com
# res.topology_valid, res.status  (CorrelationValidity.VALID or EXTRAPOLATED)

res = _core.branching_tee_residuals_and_jacobian(
    m_dot_com=0.5, m_dot_branch=0.15,
    dP0_straight=100.0, dP0_branch=180.0,
    P_static_com=1.5e5, T_com=500.0, Y_com=Y,
    theta=math.pi/3, psi=1.2, F_C=0.008,
)
```

---

## Vatistas n-Vortex Model

Reference: Vatistas, Kozel, Mih (1991), *Exp. Fluids* 11, 73–76.

Shape parameter `n >= 1` (default 2.0). `n=2` gives the best fit to most
experimental swirl data.  Closed-form antiderivatives for the pressure integral
are used for `n=1` and `n=2`; composite Simpson quadrature for general `n`.
All Jacobians are analytical.

### Free functions

```python
import combaero as cb

# Normalised shape functions (no units, r_bar = r / r_c)
v0   = cb.vatistas_v0_bar(r_bar, n)           # tangential velocity shape
dv0  = cb.vatistas_dv0_bar_drbar(r_bar, n)    # d(V0_bar)/d(r_bar)
vr   = cb.vatistas_vr_bar(r_bar, n)           # radial velocity shape (< 0)
dvr  = cb.vatistas_dvr_bar_drbar(r_bar, n)    # d(Vr_bar)/d(r_bar)
I    = cb.vatistas_pressure_integral(r_bar, n) # integral_0^r V0^2/r' dr'
dI   = cb.vatistas_d_pressure_integral_drbar(r_bar, n)  # V0_bar^2/r_bar (FTC)

# Dimensional functions (Gamma [m^2/s], r_c [m], rho [kg/m^3])
Vt   = cb.vatistas_v_theta(r, Gamma, r_c, n)  # [m/s]
Vt, dVt_dr, dVt_dG, dVt_drc = cb.vatistas_v_theta_and_jacobians(r, Gamma, r_c, n)

dP   = cb.vatistas_delta_p(r, rho, Gamma, r_c, n)  # P(r)-P(0) [Pa]
dP, ddP_dr, ddP_dG, ddP_drc = cb.vatistas_delta_p_and_jacobians(r, rho, Gamma, r_c, n)
```

### `VatistasVortex` class

```python
import numpy as np
import combaero as cb

vortex = cb.VatistasVortex(Gamma=0.5, r_c=0.05, n=2.0)

r = np.linspace(0, 0.2, 200)
V  = vortex.V_theta(r)          # [m/s], vectorised
dV = vortex.dV_theta_dr(r)      # [1/s], vectorised

Vmax = vortex.V_theta_max()     # peak tangential velocity at r = r_c

dP   = vortex.delta_P(r=0.1, rho=1.2)    # [Pa], scalar
dPbar = vortex.delta_P_bar(r=0.1)        # normalised [0, 1)
```

---

## NetworkRunner (GUI-JSON Programmatic Driver)

`NetworkRunner` loads a network saved from the GUI and runs it from Python
scripts or Jupyter notebooks — no web server required.

### Loading a network

```python
from gui.backend.runner import NetworkRunner

# From a file downloaded from the GUI
runner = NetworkRunner.from_file("my_network.json")

# From an already-loaded dict (normalises GUI camelCase keys automatically)
import json
with open("my_network.json") as f:
    runner = NetworkRunner.from_dict(json.load(f))
```

Node labels are read from the GUI's label field.  If a node has no label the
node ID is used as the fallback key for overrides and result lookup.

### Single solve with boundary-condition overrides

```python
# Override format: "<label>.<attribute>": value
result = runner.solve({
    "air_inlet.m_dot": 1.2,   # kg/s
    "air_inlet.Tt": 650.0,    # K
    "outlet.Pt": 200_000.0,   # Pa
})

print(result.success, result.final_norm)

# Scalar result by label.quantity (resolves label → node ID automatically)
T_combustor = result.get("combustor.T")         # K
phi          = result.get("combustor.phi")
m_dot_out    = result.get("hot_channel.m_dot")  # kg/s

# Full thermodynamic state dict for a node
state = result.node_state("combustor")
# {"T": 1738.0, "P": 195000.0, "Pt": 200000.0, "Tt": 1760.0, ...}

# Full-detail DataFrame (matches GUI CSV export, unit-annotated columns)
df = result.to_dataframe()
```

### Parametric sweep

```python
import pandas as pd

params = pd.DataFrame({
    "fuel_inlet.m_dot": [0.020, 0.025, 0.030, 0.035],
})

# Compact: one row per solve, only requested metrics
sweep_df = runner.sweep(params, metrics=["combustor.T", "combustor.phi"])
# columns: fuel_inlet.m_dot | combustor.T | combustor.phi | success | final_norm

# Full detail: complete to_dataframe() output for each solve, stacked
full_df = runner.sweep(params)
# columns: _sweep_index | fuel_inlet.m_dot | <all entity columns> | success | final_norm
```

### Injecting labels for unlabelled networks

GUI exports may omit node labels.  Inject them before constructing the runner:

```python
import json

with open("network.json") as f:
    schema = json.load(f)

label_map = {
    "node_1778698958490": "air_inlet",
    "node_1778698961265": "fuel_inlet",
    "node_1778699014636": "combustor",
}
for node in schema["nodes"]:
    if node["id"] in label_map:
        node["data"]["label"] = label_map[node["id"]]

runner = NetworkRunner.from_dict(schema)
```

### NetworkResult attributes

| Attribute | Type | Description |
|---|---|---|
| `success` | `bool` | True when solver converged |
| `message` | `str` | Human-readable solver status |
| `final_norm` | `float \| None` | Residual norm at convergence |

| Method | Returns | Description |
|---|---|---|
| `get(key)` | `float` | Scalar by `<id>.<qty>` or `<label>.<qty>` |
| `node_state(label)` | `dict` | Full state dict for a node |
| `to_dataframe()` | `DataFrame` | Unit-annotated full-detail export |
| `swap_boundary(label, new_type)` | `NetworkRunner` | Retype a boundary node (mass ↔ pressure) and return a new runner pre-seeded with the full solved state |

### Boundary condition swap

```python
result = runner.solve({"air_inlet.m_dot": 1.0})

# Re-run the same point with a pressure BC (Pt auto-populated from solved state)
p_runner = result.swap_boundary("air_inlet", "pressure_boundary")
result2 = p_runner.solve()  # starts from solved state, converges immediately

# Sweep inlet pressure around the design point
import pandas as pd, numpy as np
design_Pt = result.node_state("air_inlet")["Pt"]
params = pd.DataFrame({"air_inlet.Pt": design_Pt * np.linspace(0.90, 1.10, 11)})
sweep_df = p_runner.sweep(params, metrics=["combustor.T"])
```

Only `mass_boundary` ↔ `pressure_boundary` swaps are supported.  The returned
runner carries the full solved state (pressures, flows, junction states) as a
warm start so `solve()` requires no additional `init_strategy`.

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

### Rib Correlations (parametrised)

Rib correlations are **data, not code**. A parameter set carries the
coefficients, their normalisers, the validity band, the stated accuracy, and
where all of it came from -- so a built-in set and a set you tuned on your own
rig are structurally distinguishable.

```python
import combaero as cb

s = cb.han_1988_orthogonal()          # Han (1988), 90 deg orthogonal ribs
g = cb.RibGeometry(e_D=0.047, p_e=10.0, W_H=1.0, alpha_deg=90.0)

r = cb.evaluate_rib(s, g, Re=10_000)
r.R, r.f, r.e_plus, r.G, r.St_r       # 3.2000, 0.04576, 71.1, 12.21, 0.00968
r.extrapolated                        # outside the set's advisory validity
```

**Two named sets, chosen explicitly.** They cover disjoint Reynolds-number
regimes reported by different papers, not one superseding the other -- so
there is no auto-switching between them:

```python
s90 = cb.han_1988_orthogonal()        # 90 deg, Re 10,000-60,000
s45 = cb.rallabandi_2009_high_re()    # 45 deg sharp ribs, Re 30,000-400,000

g = cb.RibGeometry(e_D=0.14, p_e=7.5, W_H=1.0, alpha_deg=45.0)
r = cb.evaluate_rib(s45, g, Re=100_000)
r.G                                    # 35.40 -- see han_ribbed_high_re.md
```

`rallabandi_2009_high_re` is 45 deg and sharp-edged ribs only: the source
reports round-edged ribs at the same conditions instead following Han's
correlation, but calls the agreement "coincidental" rather than physically
grounded, so it is not wired as an automatic edge-profile switch.

**A third set covers angled ribs generally**, `han_park_1988_angled` (Eq.
4.17/4.18, `alpha` 30-90 deg, `W/H` 1-4):

```python
angled = cb.han_park_1988_angled()
g = cb.RibGeometry(e_D=0.06, p_e=15.0, W_H=2.0, alpha_deg=60.0)
r = cb.evaluate_rib(angled, g, Re=30_000)
r.R, r.G                              # 3.2332, 18.4408
```

Its `R` and `G` carry genuine, unsmoothed discontinuities the source
states rather than a numerical artefact: `R`'s `(W/H)^m` term switches at
`alpha == 90 deg` (up to 62.5% at `W/H = 4`), and `G`'s `alpha`/`p_e`
exponents switch on whether the channel is square (`W/H == 1`, up to 27%
at low `alpha`). Neither the printed equations nor the extraction record a
smooth transition between the branches, so none is invented -- a caller
whose solver traverses either boundary exactly needs its own guard.

`RAlphaShape` and `GShapeModel` (both importable from `combaero`) are what
make this representable as data: setting `R_alpha_shape` to
`QuadraticAlpha` and populating `R_quad_c0/c1/c2` plus the two
`R_quad_WH_exponent_*` fields lets a user-supplied set describe its own
angle-dependent `R`, the same way `RibTerm`'s normaliser lets one describe
its own geometry dependence.

**Supply your own.** Real hardware needs it -- no published correlation is
precise enough for a specific rig:

```python
mine = cb.han_1988_orthogonal()
mine.name = "rig_3"
mine.source = "measured 2026-09"
mine.provenance = cb.RibProvenance.User   # not Extracted -- the claim differs
mine.C_G = 4.1
cb.validate_rib_set(mine)                 # rejects malformed sets, loudly
```

Three things worth knowing:

- **The normaliser is data.** `RibTerm(exponent, reference)` divides by
  `reference` before applying the exponent. `3.2 (p/e/10)^0.35` and
  `1.4294 (p/e)^0.35` are the same function; using one constant under the
  other convention is wrong by `10^0.35 = 2.24x`, uniformly, which never looks
  like a trend.
- **Validity is advisory.** `r.extrapolated` reports; nothing refuses. A band
  belongs to the source's rig, not to yours.
- **`St_r` is the ribbed side.** Combining it with the smooth walls is the
  caller's job, and depends on how many walls are ribbed.

Bad *parameters* raise from `validate_rib_set`. Bad *operating points* never
raise: reverse flow, zero flow and extreme values are guarded smoothly, because
a solver probes states that are not physical and a throw inside a residual
kills the solve.

### Ribbed Channels

```python
from combaero.network import ChannelElement, ConvectiveSurface, RibbedModel

surface = ConvectiveSurface(
    area=0.1,
    model=RibbedModel(
        e_D=0.06, p_e=10.0, alpha_deg=90.0,
        W_H=1.0,                  # aspect ratio: a Dh does not determine it
        n_ribbed_walls=2,         # 2 = two OPPOSITE walls, as Han measured
    ),
)
```

`n_ribbed_walls` accepts 1, 2 or 4. Three is rejected -- it has no unambiguous
geometry.

**Two asymmetries, both from the source rather than convenience:**

- **Friction needs no wall weighting.** The correlation's `f` is already the
  four-sided channel value, so the element uses it directly. Nothing multiplies
  pipe friction.
- **Heat transfer does.** The correlation gives the ribbed side; the smooth
  walls come from the base correlation and the channel average is the
  area-weighted combination. The result exposes `h_ribbed` and `h_smooth`
  separately, because the average hides a modelling choice worth seeing.

**A documented gap.** Plain smooth walls give `h_s/h_r` around 0.42 where Han's
own channel average implies 0.70 -- ribs enhance the adjacent smooth wall by
10-50% too, which no correlation here covers. The channel average is therefore
about **20% below** Han's measurement for the two-ribbed-wall square case.

```python
model = RibbedModel(..., smooth_wall_Nu_multiplier=1.67)   # reproduces Han
```

That knob exists because the ribbed side already has a better one: change `C_G`
on the parameter set, which records what you changed and why. The smooth walls
come from Gnielinski and have no set of their own.

### Impingement Channels

```python
from combaero.network import ConvectiveSurface, ImpingementModel, SingleJetImpingementModel

surface = ConvectiveSurface(
    area=0.01,   # this ROW's own target-plate footprint
    model=ImpingementModel(
        d_jet=0.002, xn_d=8.0, yn_d=6.0, z_d=2.0,
        row=3,       # 1-indexed, counting from upstream
    ),
)
```

**One element models one spanwise row**, not a whole array -- there is no
single "channel Nu" for a jet array the way there is for a smooth or ribbed
duct, since downstream rows see progressively more crossflow than upstream
ones. Model a real array by chaining `n_rows` separate elements, each with
its own `row`; row 1 sees zero crossflow by definition (`Gc/Gj = 0`,
Florschuetz's own `Nu1`).

**Mass flow is recovered from `area`, the same convention every
`ConvectiveSurface` model uses** -- `ImpingementModel` does not introduce a
new one. The element treats `area` as this row's own footprint, divides by
`xn_d * d_jet * yn_d * d_jet` to get the row's hole count, and gets this
row's own jet velocity from the upstream state's total mass flow through
that area. Whether every row gets the same total flow (a uniform-supply
approximation) or a row-dependent one (Florschuetz's own Eq. 7, deliberately
not implemented -- see `impingement_correlation.h`'s module comment) is a
choice for whoever assembles the chain, not something this element decides.

**Pressure drop is not modelled here.** Model the jet plate's own orifice
loss with a proper `OrificeElement` upstream, using a real
discharge-coefficient correlation -- duplicating it inside this heat-transfer
model would risk double-counting it. `f`/`dP` reported by this element are
the crossflow's own plain smooth-duct values (Han's book gives no
impingement-specific friction correlation, unlike ribs' own `R`/`f`
relationship), borrowed to compute a `T_aw` Eq. 4.2/4.3's own recovery-factor
model would otherwise have to supply.

**A single free jet** (Goldstein, Behbahani and Heppelmann 1986) is a
separate model, since there is no array and no crossflow:

```python
surface = ConvectiveSurface(
    area=0.001,   # this jet's own target patch
    model=SingleJetImpingementModel(d_jet=0.003, L_D=7.75, R_D=5.0),
)
```

`R_D` is a real modelling choice, not a detail: the correlation's `Nu` is
LOCAL to the radial position, so this reports one representative value
rather than an area-averaged profile.

Both models expose the same wall-coupling derivatives (`dh_dmdot`, `dh_dT`,
`dT_aw_dmdot`, `dT_aw_dT`) and an `extrapolated` flag (array only -- a single
jet has no stated validity range to be outside of), matching `RibbedModel`'s
contract.

### Jet Impingement Correlations (parametrised)

Two independent regimes, matching
`validation/cooling/extractions/han_impingement.md`'s own split, wired into
the network elements above (`ImpingementModel`, `SingleJetImpingementModel`).

**Single jet** (Goldstein, Behbahani and Heppelmann, 1986): one free round
jet on a flat plate, no crossflow, no array.

```python
import combaero as cb

s = cb.goldstein_1986_single_jet()
nu_q = cb.single_jet_impingement_nu(
    s, cb.ImpingementThermalBC.ConstantHeatFlux, Re=25_000, L_D=7.75, R_D=5.0)
nu_t = cb.single_jet_impingement_nu(
    s, cb.ImpingementThermalBC.ConstantWallTemperature, Re=25_000, L_D=7.75, R_D=5.0)
nu_q, nu_t                            # 59.93, 55.71 -- Han's own check point (60, 56)
```

`L_D = 7.75` is the optimum spacing -- structural, not a separately fitted
fact: `A` itself is the numerator's value there, for any `Re` or `R_D`.

**Jet array with crossflow** (Florschuetz, Truman and Metzger, 1981): a
staggered or inline array fed from a common plenum, where every row's Nu is
degraded by the crossflow accumulated from every row upstream of it.

```python
inline = cb.florschuetz_1981_inline()

result = cb.jet_array_impingement_nu(
    inline, Re_j=10_000, Gc_Gj=0.3, Pr=0.7, xn_d=10.0, yn_d=6.0, z_d=2.0)
result.Nu, result.extrapolated        # 27.93, False
```

`Re_j` and `Gc_Gj` are that ROW's own local values, not an array mean --
Florschuetz correlates "the individual spanwise row jet Reynolds number".
`Gc_Gj` at a row has its own closed form (the paper's own Eq. 8), needing
only geometry and a jet-plate discharge coefficient:

```python
cb.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=cb.FLORSCHUETZ_1981_DEFAULT_CD, row=1)   # 0.0 -- Nu1's own definition
cb.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=cb.FLORSCHUETZ_1981_DEFAULT_CD, row=10)  # 0.404
```

It depends on `(yn/d)(z/d)` only, not `xn/d` -- the source states the flow
distribution is independent of streamwise hole spacing and hole pattern.
`FLORSCHUETZ_1981_DEFAULT_CD` (0.79) is the paper's own recommended default
absent a measured value, and remains the default. To predict `C_D` from the
plate's own geometry instead, pass `mcgreehan_schotsch_1988_cd` (see above);
for a plenum-fed plate leave its `U1_over_Vi` at 0, since Florschuetz's
`Gc/Gj` is discharge-side crossflow and must not be fed to it:

```python
C_D = cb.mcgreehan_schotsch_1988_cd(Re=1e4, r_over_d=0.0, L_over_d=t_over_d)
cb.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=C_D, row=10)
```

At `t/d = 1` the two agree to 0.6%, so this changes little; at `t/d = 0.5` the
correlation gives 0.688 and the predicted `Gc/Gj` rises by up to 12% at
downstream rows (about 1% in `Nu`). The ISO 5167 `Cd_sharp_thin_plate` family
in `orifice.h` remains inapplicable here -- it models a hole in a pipe run.

**Two named sets**, chosen explicitly by hole pattern -- their coefficients
differ, not just their validity:

```python
inl = cb.florschuetz_1981_inline()      # xn/d 5-15
stg = cb.florschuetz_1981_staggered()   # xn/d 5-10, genuinely tighter
```

Bad *parameters* raise from `validate_single_jet_set`/`validate_jet_array_set`.
Bad *operating points* never raise: reverse flow and out-of-range crossflow
ratios are guarded smoothly, matching the ribs' policy above.

Leading-edge/curved-surface impingement (Section 4.1.4) and the simpler,
looser forms (Eq. 4.6 Kercher-Tabakoff -- graphical, not closed-form anyway;
Eq. 4.7/4.8 -- Florschuetz's own less-tight alternate) are out of scope,
deferred per the extraction's I3.

### Enhanced Cooling Surfaces

Removed in 0.7.0. The pin-fin, dimple, rib and impingement correlations could
not be traced to their cited sources -- the rib friction multiplier was 4-5x
below the only rib datum in the repository. Provenanced rib and jet
impingement correlations are re-added above (issues #334, #337); see issue
#339 for the rebuild.

`channel_smooth` and the base convective correlations (Gnielinski,
Dittus-Boelter, Sieder-Tate, Petukhov) are unaffected, as are the user-set
`Nu_multiplier` and `f_multiplier` knobs on `ConvectiveSurface`.

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
# Molar-basis inverses (targets in J/mol or J/(mol*K))
T_from_h   = cb.calc_T_from_h(h_target=3e5, X=air)            # J/mol -> K
T_from_s   = cb.calc_T_from_s(s_target=7000, P=101325, X=air) # J/(mol*K), Pa -> K
T_from_cp  = cb.calc_T_from_cp(cp_target=29, X=air)           # J/(mol*K) -> K
T_from_u   = cb.calc_T_from_u(u_target=2.2e5, X=air)          # J/mol -> K

# Mass-specific inverses (targets in J/kg or J/(kg*K))
T_from_h_m = cb.calc_T_from_h_mass(h_mass_target=3e5, X=air)
T_from_s_m = cb.calc_T_from_s_mass(s_mass_target=7000, P=101325, X=air)
T_from_u_m = cb.calc_T_from_u_mass(u_mass_target=2.2e5, X=air)

# Flash calculations — solve T from two independent properties
T_from_sv  = cb.calc_T_from_sv_mass(s_mass_target=7000, v_mass_target=0.8, X=air)
T_from_sh  = cb.calc_T_from_sh_mass(s_mass_target=7000, h_mass_target=3e5, X=air)
```

### Flow Analysis
```python
# Bulk velocity and kinetic energy
v = cb.bulk_velocity(m_dot=0.5, rho=1.2, area=1e-3)   # m/s
KE = cb.kinetic_energy(v=v)                             # J/kg

# Dimensionless numbers
# From full mixture state (recomputes rho, mu internally):
Re = cb.reynolds(T=300, P=101325, X=air, V=v, L=0.05)
# From pre-computed state properties (use when CompleteState is available):
Re = cb.reynolds_from_state(rho=1.2, v=v, L=0.05, mu=1.8e-5)

# Stagnation conditions from static state and velocity
Tt, Pt = cb.stagnation_from_static(T=300, P=101325, v=v, X=air)  # (K, Pa)
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
