# Network Convective Heat Transfer — Design Document

This document describes the architecture and implementation plan for adding
convective heat transfer coupling to the CombAero network solver.

---

## 1. Motivation

The current network solver treats elements as **pressure-loss devices** and
nodes as **mixing / energy-balance points**.  Temperature is forward-propagated
through the graph via `compute_derived_state()` at each node; elements never
modify the fluid temperature.

For convective heat transfer (turbine blade cooling, combustor liner cooling,
heat exchangers) an element must be able to **transfer heat between two flow
paths** that share a wall.  The transferred heat rate Q depends on the flow
state of *both* sides (Re, Pr → h) and on wall geometry and material.

### Target use-case

Combustor with counterflow liner cooling:

```
  Hot gas path:     [Combustor] ── hot liner ──── [Turbine Inlet]
                                      ║ Q (wall)
  Cooling path:     [Supply]    ── cool passage ── [Cool Exit] ── dilution ── [Combustor]
```

The cooling air picks up heat along the passage, then feeds back into the
combustor as dilution air — creating a **thermally-coupled loop** that the
solver must converge simultaneously with pressure and mass flow.

---

## 2. Architecture Overview

Three new concepts, layered on the existing solver without breaking current
behaviour:

| Concept | Role |
|---------|------|
| **`ConvectiveSurface`** | Per-element dataclass describing the wetted surface (area, model, correlation, empirical multipliers) |
| **`htc_and_T()`** | Method on `NetworkElement` returning (h, T_aw, A_conv) from the current flow state |
| **`WallConnection`** | Coupler object connecting two elements through a shared wall; computes Q and its Jacobians |

### Operating modes

| `A_conv` | `WallConnection` | Behaviour |
|----------|-------------------|-----------|
| 0 (default) | — | No HTC computed.  Fast pressure-flow only. |
| > 0 | absent | HTC computed but **post-processing only** — no solver feedback. |
| > 0 | present | **Fully coupled** thermal solve via Newton. |

A global toggle `FlowNetwork.thermal_coupling_enabled` (default `True`)
allows disabling all wall couplings for debugging without removing them from
the graph.

---

## 3. `ConvectiveSurface` and Model Subclasses

Each convective surface carries:

- **`area`** `[m²]` — wetted convective area.  Default `0.0` (disabled).
- **`model`** — a model-specific subclass holding the geometry parameters that
  map 1-to-1 to a C++ `channel_*` function.
- **`correlation`** — Nusselt correlation selector (smooth model only).
- **`heating`** — `bool | None`.  `None` = auto-detect from sign of
  `T_wall − T_fluid`.
- **`Nu_multiplier`** — empirical correction factor on Nu (default `1.0`).
- **`f_multiplier`** — empirical correction factor on friction (default `1.0`).

### Model subclasses

```python
from dataclasses import dataclass

@dataclass
class SmoothModel:
    """Parameters for channel_smooth."""
    correlation: str = "gnielinski"      # "gnielinski" | "dittus_boelter" | "sieder_tate" | "petukhov"
    mu_ratio: float = 1.0               # μ_bulk / μ_wall (Sieder-Tate)
    roughness: float = 0.0              # absolute roughness [m]

@dataclass
class RibbedModel:
    """Parameters for channel_ribbed."""
    e_D: float = 0.0                    # rib height / hydraulic diameter
    pitch_to_height: float = 0.0        # rib pitch / rib height
    alpha_deg: float = 90.0             # rib angle [deg]

@dataclass
class DimpledModel:
    """Parameters for channel_dimpled."""
    d_Dh: float = 0.0                   # dimple diameter / hydraulic diameter
    h_d: float = 0.0                    # dimple depth / dimple diameter
    S_d: float = 0.0                    # dimple pitch / dimple diameter

@dataclass
class PinFinModel:
    """Parameters for channel_pin_fin."""
    pin_diameter: float = 0.0           # pin diameter [m]
    channel_height: float = 0.0         # channel height [m]
    S_D: float = 2.0                    # transverse pitch / pin diameter
    X_D: float = 2.0                    # streamwise pitch / pin diameter
    N_rows: int = 1                     # number of pin rows
    is_staggered: bool = True

@dataclass
class ImpingementModel:
    """Parameters for channel_impingement."""
    d_jet: float = 0.0                  # jet hole diameter [m]
    z_D: float = 0.0                    # jet-to-target distance / d_jet
    x_D: float = 0.0                    # streamwise pitch / d_jet
    y_D: float = 0.0                    # spanwise pitch / d_jet
    A_target: float = 0.0              # target area [m²]
    Cd_jet: float = 0.8                 # jet discharge coefficient


ChannelModel = SmoothModel | RibbedModel | DimpledModel | PinFinModel | ImpingementModel


@dataclass
class ConvectiveSurface:
    """Convective surface description attached to a NetworkElement."""
    area: float = 0.0                            # A_conv [m²] — 0 disables
    model: ChannelModel = field(default_factory=SmoothModel)
    heating: bool | None = None                  # None = auto-detect
    Nu_multiplier: float = 1.0                   # empirical correction on Nu
    f_multiplier: float = 1.0                    # empirical correction on f
```

Each subclass maps directly to the corresponding C++ `channel_*` function,
making misconfiguration (e.g. setting rib parameters on a smooth model) a
type-level error.

---

## 4. Element Integration

### Base class extension

```python
class NetworkElement(ABC):
    ...
    def htc_and_T(self, state: MixtureState) -> tuple[float, float, float] | None:
        """Return (h [W/(m²·K)], T_aw [K], A_conv [m²]) or None."""
        return None   # default: no convective surface
```

### `PipeElement` override

```python
class PipeElement(NetworkElement):
    def __init__(self, ..., surface: ConvectiveSurface | None = None):
        ...
        self.surface = surface or ConvectiveSurface()  # default: area=0

    def htc_and_T(self, state: MixtureState) -> tuple[float, float, float] | None:
        if self.surface.area == 0.0:
            return None

        rho = density(state.T, state.P, state.Y)
        u = state.m_dot / (rho * self.area)  # element cross-section area

        # Dispatch to appropriate C++ channel_* function based on model type
        result = _dispatch_channel(
            self.surface.model, state.T, state.P, state.Y,
            u=u, diameter=self.diameter, length=self.length,
            heating=self.surface.heating,
            Nu_multiplier=self.surface.Nu_multiplier,
            f_multiplier=self.surface.f_multiplier,
        )
        return result.h, result.T_aw, self.surface.area
```

The velocity `u = mdot / (ρ · A_cross)` comes from the element's own
hydraulic cross-section — independent of node type.  This works correctly
whether the element connects to a plenum or a momentum chamber.

---

## 5. `WallConnection`

```python
@dataclass
class WallConnection:
    """Thermal coupling between two elements through a shared wall."""
    id: str
    element_a: str                       # element ID — side A
    element_b: str                       # element ID — side B
    wall_thickness: float                # [m]
    wall_conductivity: float             # [W/(m·K)]
    contact_area: float | None = None    # override A_conv if different
```

### Heat transfer calculation

Given the current flow state of both elements:

1. Call `htc_and_T()` on each element → `(h_a, T_aw_a, A_a)` and
   `(h_b, T_aw_b, A_b)`.
2. Compute overall HTC through the wall:
   `U = overall_htc_wall(h_a, h_b, t/k)` — already in C++ core.
3. Use the smaller of `A_a`, `A_b`, or `contact_area` as the effective area.
4. Compute `Q = U · A · (T_aw_a − T_aw_b)` (positive = heat flows A→B).
5. Apply `−Q` to element A's downstream node, `+Q` to element B's downstream
   node via the existing `EnergyBoundary` mechanism.

### `heating` auto-detection

When `heating` is `None` on a `ConvectiveSurface`, the wall temperature
`T_wall` is needed to determine the Dittus-Boelter exponent.  `T_wall` is
computed by the `WallConnection` from the balance of both sides:

```
T_wall = (h_a · T_aw_a + h_b · T_aw_b) / (h_a + h_b)  (simplified, no wall resistance)
```

Each element then checks `T_wall > T_fluid` to set `heating = True` or
`False`.  The user can override by setting `heating` explicitly on the
`ConvectiveSurface`.

### Registration

```python
class FlowNetwork:
    def __init__(self):
        ...
        self.walls: dict[str, WallConnection] = {}
        self.thermal_coupling_enabled: bool = True

    def add_wall(self, wall: WallConnection) -> None:
        ...
```

---

## 6. Solver Integration — Proper Newton Coupling

### Approach: flow-dependent `EnergyBoundary`

The solver evaluates `WallConnection`s **during** `_propagate_states()`, so
that `Q` participates in the node mixing and its derivatives appear in the
Jacobian.  This gives **quadratic convergence** for the thermal coupling.

#### Mechanism

In `NetworkSolver._propagate_states()`, after constructing `state_up` for
each upstream element of a node:

1. Check if the element has a `WallConnection`.
2. If so, read the partner element's upstream state (which must already be
   propagated — see ordering below).
3. Call `htc_and_T()` on both elements.
4. Compute `Q` and its Jacobians via `wall_coupling_and_jacobian()` (C++).
5. Set `Q` on a flow-dependent `EnergyBoundary` attached to this node.
6. The mixer's existing `Q` parameter in `mixer_from_streams_and_jacobians()`
   receives this Q.

#### Topological ordering and cross-stream coupling

For a wall connecting element_hot (A→B) and element_cold (C→D):

- When processing node B: need `T_A` (upstream of hot, already propagated)
  and `T_C` (upstream of cold — may or may not be propagated yet).
- If `T_C` is not yet available (topological ordering issue), use the
  **previous Newton iteration's** value.  This is standard Newton behaviour
  for back-edges in cyclic graphs.
- The Jacobian includes `dQ/dT_A`, `dQ/dT_C`, `dQ/dmdot_hot`,
  `dQ/dmdot_cold` — chained into the sensitivity relay.

For acyclic topologies (the common case), both upstream states are available
and Q is exact at every Newton step.

#### Sensitivity relay extension

The existing relay tracks `dT/dx`, `dY/dx`, `dP_total/dx` for each node.
The wall coupling adds:

```
dT_node_B/dx = dT_mix/dQ · dQ/dx + (existing terms)

where  dQ/dx = dQ/dT_A · dT_A/dx + dQ/dmdot_hot · dmdot_hot/dx
             + dQ/dT_C · dT_C/dx + dQ/dmdot_cold · dmdot_cold/dx
```

The `dT_mix/dQ` term already exists — it is the `dT_mix_d_delta_h` field in
`MixerResult`.

#### Global toggle

```python
if not self.network.thermal_coupling_enabled:
    # Skip all wall evaluations — pure pressure-flow solve
    pass
```

### Rejected alternative: separate thermal pass (lagged Q)

An alternative is to evaluate walls **after** full state propagation and use
the resulting Q in the **next** Newton iteration.  This is simpler but yields:

- **Linear convergence** for the thermal coupling (Q lags by one iteration).
- Risk of oscillation or divergence for thermally-coupled loops (combustor ↔
  cooling passage) where the loop gain exceeds 1.
- Effectively removes `dQ/dx` from the Jacobian.

The proper Newton approach is preferred because:

1. Heat transfer in cooling passages is moderately coupled — the loop gain is
   usually < 1 but not negligible.
2. The combustor feedback loop (coolant → dilution → combustion → hot gas →
   wall → coolant) requires tight coupling to converge reliably.
3. The existing `EnergyBoundary` and `MixerResult.dT_mix_d_delta_h` mechanism
   already supports Q injection with analytical Jacobians — the infrastructure
   exists.

---

## 7. C++ Core Changes

### 7.1 Extended `ChannelResult`

Add Jacobian fields to the existing struct in `heat_transfer.h`:

```cpp
struct ChannelResult {
  // Existing fields
  double h;       // convective HTC [W/(m²·K)]
  double Nu;      // Nusselt number [-]
  double Re;      // Reynolds number [-]
  double Pr;      // Prandtl number [-]
  double f;       // friction / loss coefficient [-]
  double dP;      // pressure drop [Pa]
  double M;       // Mach number [-]
  double T_aw;    // adiabatic wall temperature [K]
  double q;       // heat flux [W/m²]

  // Jacobians (new) — derivatives w.r.t. mass flow and temperature
  double dh_dmdot;     // ∂h/∂ṁ  [W/(m²·K) / (kg/s)]
  double dh_dT;        // ∂h/∂T  [W/(m²·K) / K]
  double ddP_dmdot;    // ∂dP/∂ṁ [Pa / (kg/s)]
  double ddP_dT;       // ∂dP/∂T [Pa / K]
  double dq_dmdot;     // ∂q/∂ṁ  [W/m² / (kg/s)]
  double dq_dT;        // ∂q/∂T  [W/m² / K]
  double dq_dT_wall;   // ∂q/∂T_wall [W/m² / K]  = −h (once h is known)
};
```

The Jacobians compose existing building blocks from `solver_interface.h`:

- `nusselt_and_jacobian_gnielinski(Re, Pr, f)` → `(Nu, dNu/dRe)`
- `friction_and_jacobian(tag, Re, e_D)` → `(f, df/dRe)`
- `density_and_jacobians(T, P, X)` → `(ρ, dρ/dT, dρ/dP)`
- `enthalpy_and_jacobian(T, X)` → `(h_mass, cp)`

Chain rule through `Re = ρ · u · D / μ` and `u = mdot / (ρ · A)`:

```
dRe/dmdot = Re / mdot
dRe/dT    = Re · (1/ρ · dρ/dT − 1/μ · dμ/dT)
dNu/dmdot = dNu/dRe · dRe/dmdot
dh/dmdot  = dNu/dmdot · k / D
```

### 7.2 `Nu_multiplier` and `f_multiplier`

Add to all `channel_*` function signatures:

```cpp
ChannelResult
channel_smooth(double T, double P, const std::vector<double> &X,
               double velocity, double diameter, double length,
               double T_wall = NaN,
               const std::string &correlation = "gnielinski",
               bool heating = true, double mu_ratio = 1.0,
               double roughness = 0.0,
               double Nu_multiplier = 1.0,    // new
               double f_multiplier = 1.0);    // new
```

Applied inside the function **before** computing derived quantities:

```cpp
Nu *= Nu_multiplier;
f  *= f_multiplier;
h   = htc_from_nusselt(Nu, k, diameter);
dP  = f * (L / D) * 0.5 * rho * v * v;
```

Since the multipliers are constants, all Jacobians simply scale:

```
dh_corrected/dmdot = Nu_multiplier · dh_base/dmdot
ddP_corrected/dmdot = f_multiplier · ddP_base/dmdot
```

The same pattern applies to `channel_ribbed`, `channel_dimpled`,
`channel_pin_fin`, and `channel_impingement`.

### 7.3 `WallCouplingResult`

New struct and function:

```cpp
struct WallCouplingResult {
  double Q;             // heat transfer rate [W] (positive A→B)
  double T_wall;        // wall temperature [K] (hot-side surface)
  double dQ_dh_a;       // ∂Q/∂h_a
  double dQ_dh_b;       // ∂Q/∂h_b
  double dQ_dT_aw_a;    // ∂Q/∂T_aw_a
  double dQ_dT_aw_b;    // ∂Q/∂T_aw_b
};

WallCouplingResult wall_coupling_and_jacobian(
    double h_a, double T_aw_a,
    double h_b, double T_aw_b,
    double t_over_k,        // wall_thickness / wall_conductivity [m²·K/W]
    double A);              // effective heat transfer area [m²]
```

Uses `overall_htc_wall(h_a, h_b, t_over_k)` internally.  Derivatives are
analytical since `U = 1 / (1/h_a + t/k + 1/h_b)` is a simple rational
function of h_a, h_b.

### 7.4 Pybind11 Bindings

- Bind new `ChannelResult` Jacobian fields as `.def_readonly()`.
- Bind `Nu_multiplier`, `f_multiplier` as `py::arg()` with default `1.0`.
- Bind `WallCouplingResult` and `wall_coupling_and_jacobian`.

---

## 8. User-Facing API Examples

### Simple pipe-to-pipe wall coupling

```python
import combaero as cb
from combaero.network import (
    FlowNetwork, PressureBoundary, PlenumNode,
    PipeElement, WallConnection, ConvectiveSurface,
    SmoothModel, RibbedModel,
)
from math import pi

net = FlowNetwork()

# Hot side
net.add_node(PressureBoundary("hot_in",  P_total=500e3, T_total=1200.0))
net.add_node(PlenumNode("hot_out_plenum"))
net.add_node(PressureBoundary("hot_out", P_total=480e3, T_total=1200.0))

net.add_element(PipeElement(
    "hot_pipe", "hot_in", "hot_out_plenum",
    length=0.3, diameter=0.05, roughness=1e-4,
    surface=ConvectiveSurface(
        area=pi * 0.05 * 0.3,
        model=SmoothModel(correlation="gnielinski"),
    ),
))

# Cold side
net.add_node(PressureBoundary("cold_in",  P_total=600e3, T_total=400.0))
net.add_node(PlenumNode("cold_out_plenum"))
net.add_node(PressureBoundary("cold_out", P_total=580e3, T_total=400.0))

net.add_element(PipeElement(
    "cold_pipe", "cold_in", "cold_out_plenum",
    length=0.3, diameter=0.005, roughness=1e-5,
    surface=ConvectiveSurface(
        area=pi * 0.005 * 0.3,
        model=RibbedModel(e_D=0.05, pitch_to_height=10, alpha_deg=60),
        Nu_multiplier=1.15,    # 15% enhancement from test data
    ),
))

# Wall coupling
net.add_wall(WallConnection(
    "liner_wall",
    element_a="hot_pipe",
    element_b="cold_pipe",
    wall_thickness=0.002,
    wall_conductivity=25.0,   # Hastelloy X
))

result = net.solve()
```

### Combustor with counterflow liner cooling

```python
from combaero.network import (
    FlowNetwork, PressureBoundary, PlenumNode, MomentumChamberNode,
    CombustionNode, PipeElement, OrificeElement,
    WallConnection, ConvectiveSurface, SmoothModel, RibbedModel,
)
from math import pi

net = FlowNetwork()

# --- Boundaries ---
net.add_node(PressureBoundary("air_supply", P_total=3e6, T_total=700.0))
net.add_node(PressureBoundary("fuel_supply", P_total=3.5e6, T_total=300.0))
net.add_node(PressureBoundary("turbine_inlet", P_total=2.8e6))

# --- Hot gas path ---
net.add_node(CombustionNode("combustor", method="equilibrium"))
net.add_element(OrificeElement("primary_air", "air_supply", "combustor",
                               Cd=0.7, diameter=0.035682))
net.add_element(OrificeElement("fuel_injector", "fuel_supply", "combustor",
                               Cd=0.6, diameter=0.015958))
net.add_element(PipeElement(
    "hot_liner", "combustor", "turbine_inlet",
    length=0.4, diameter=0.15, roughness=1e-4,
    surface=ConvectiveSurface(
        area=pi * 0.15 * 0.4,
        model=SmoothModel(correlation="gnielinski"),
    ),
))

# --- Cooling path ---
net.add_node(PlenumNode("cool_plenum"))
net.add_element(OrificeElement("bleed_orifice", "air_supply", "cool_plenum",
                               Cd=0.65, diameter=0.025231))
net.add_node(MomentumChamberNode("cool_exit", area=pi/4 * 0.005**2))
net.add_element(PipeElement(
    "cooling_passage", "cool_plenum", "cool_exit",
    length=0.4, diameter=0.005, roughness=1e-5,
    surface=ConvectiveSurface(
        area=pi * 0.005 * 0.4,
        model=RibbedModel(e_D=0.05, pitch_to_height=10, alpha_deg=60),
        Nu_multiplier=1.1,
    ),
))

# Heated coolant feeds back as dilution air
net.add_element(OrificeElement("dilution_holes", "cool_exit", "combustor",
                               Cd=0.7, diameter=0.031915))

# --- Wall coupling ---
net.add_wall(WallConnection(
    "combustor_liner",
    element_a="hot_liner",
    element_b="cooling_passage",
    wall_thickness=0.002,
    wall_conductivity=25.0,
))

result = net.solve()

# Post-processing: inspect thermal margins
for wall_id, wall_result in result.wall_results.items():
    print(f"{wall_id}: Q={wall_result.Q:.1f} W, T_wall={wall_result.T_wall:.1f} K")
```

---

## 9. Phased Implementation Plan

### Phase 1 — C++ `channel_smooth` Jacobians + multipliers

**Scope**: Extend `ChannelResult` struct, implement analytical Jacobians for
`channel_smooth`, add `Nu_multiplier` / `f_multiplier` parameters, pybind11
bindings, C++ and Python unit tests.

**Deliverables**:
- `heat_transfer.h`: extended `ChannelResult` with Jacobian fields
- `heat_transfer.cpp`: `channel_smooth` computes Jacobians via chain rule
  through existing `nusselt_and_jacobian_*` and `friction_and_jacobian`
- `_core.cpp`: bind new fields and parameters
- `tests/test_heat_transfer_jacobians.cpp`: validate vs finite differences
- `python/tests/test_heat_transfer_jacobians.py`: validate from Python

### Phase 2 — C++ remaining channel variants + wall coupling

**Scope**: Extend Jacobians to `channel_ribbed`, `channel_dimpled`,
`channel_pin_fin`, `channel_impingement`.  Implement `WallCouplingResult` and
`wall_coupling_and_jacobian`.

**Deliverables**:
- All `channel_*` functions compute Jacobians and accept multipliers
- `WallCouplingResult` struct and `wall_coupling_and_jacobian` in
  `solver_interface.h` / `solver_interface.cpp`
- Pybind11 bindings for all new structs/functions
- C++ and Python tests

### Phase 3 — Python `ConvectiveSurface` + `htc_and_T()`

**Scope**: Implement `ConvectiveSurface`, model subclasses, `htc_and_T()`
method on `NetworkElement` and `PipeElement`.  Post-processing mode only (no
solver coupling yet).

**Deliverables**:
- `components.py`: `ConvectiveSurface`, model subclasses, `htc_and_T()`
- Python tests: verify `htc_and_T()` returns correct values for each model
- Python tests: verify `A_conv = 0` returns `None`

### Phase 4 — Python `WallConnection` + Newton solver integration

**Scope**: Implement `WallConnection`, `FlowNetwork.add_wall()`, solver
integration with flow-dependent `EnergyBoundary`, Jacobian relay extension.

**Deliverables**:
- `components.py`: `WallConnection` class
- `graph.py`: `add_wall()`, serialization support
- `solver.py`: wall evaluation during `_propagate_states()`, relay extension
- `thermal_coupling_enabled` global toggle
- Integration tests: two-pipe wall coupling converges to correct Q
- Integration tests: verify Jacobian accuracy (analytical vs FD)

### Phase 5 — Combustor counterflow demo + documentation

**Scope**: End-to-end combustor liner cooling example, performance
benchmarking, API documentation update.

**Deliverables**:
- Example script / notebook: combustor with counterflow cooling
- Verify thermal loop convergence (coolant → dilution → combustion → wall)
- Update `API_PYTHON.md` with new classes and usage
- Update `API_CPP.md` with new structs and functions

---

## 10. Testing Strategy

| Level | What | How |
|-------|------|-----|
| **C++ unit** | Each `channel_*` Jacobian | Compare analytical vs central FD (ε = 1e-6), require < 0.1% error |
| **C++ unit** | `wall_coupling_and_jacobian` | Validate against known analytical solution |
| **Python unit** | `ConvectiveSurface` + `htc_and_T()` | Verify correct dispatch to C++ channel functions |
| **Python unit** | `A_conv = 0` returns `None` | No C++ call made |
| **Python integration** | Two-pipe wall coupling | Converged Q matches analytical solution within tolerance |
| **Python integration** | Combustor thermal loop | Newton converges, T_wall in expected range |
| **Python integration** | `thermal_coupling_enabled = False` | Identical to no-wall solve |
| **Jacobian validation** | Full network Jacobian | Compare assembled J vs column-wise FD perturbation |

---

## 11. Existing C++ Building Blocks

Functions already in the codebase that this design composes:

| Function | Location | Role |
|----------|----------|------|
| `nusselt_and_jacobian_gnielinski` | `solver_interface.h` | `(Nu, dNu/dRe)` |
| `nusselt_and_jacobian_dittus_boelter` | `solver_interface.h` | `(Nu, dNu/dRe)` |
| `nusselt_and_jacobian_sieder_tate` | `solver_interface.h` | `(Nu, dNu/dRe)` |
| `nusselt_and_jacobian_petukhov` | `solver_interface.h` | `(Nu, dNu/dRe)` |
| `friction_and_jacobian` | `solver_interface.h` | `(f, df/dRe)` dispatcher |
| `density_and_jacobians` | `solver_interface.h` | `(ρ, dρ/dT, dρ/dP)` |
| `enthalpy_and_jacobian` | `solver_interface.h` | `(h_mass, cp)` |
| `overall_htc_wall` | `heat_transfer.h` | `U` from inner/outer HTC + wall |
| `wall_temperature_profile` | `heat_transfer.h` | Multi-layer T_wall |
| `pin_fin_nusselt_and_jacobian` | `solver_interface.h` | `(Nu, dNu/dRe)` |
| `dimple_nusselt_enhancement_and_jacobian` | `solver_interface.h` | `(enh, denh/dRe)` |
| `rib_enhancement_factor_high_re_and_jacobian` | `solver_interface.h` | `(enh, denh/dRe)` |
| `impingement_nusselt_and_jacobian` | `solver_interface.h` | `(Nu, dNu/dRe)` |
| `mixer_from_streams_and_jacobians` | `solver_interface.h` | Mixing with Q input |
