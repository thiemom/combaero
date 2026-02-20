# CombAero Network Solver — Project Roadmap

## Vision

A Python-based thermodynamic network solver with a React GUI, built on top of
CombAero's C++ physics kernels. Users draw networks of combustion system
elements, specify boundary conditions, and solve for pressures, temperatures,
mass flows, and compositions throughout the system.

---

## Tech Stack

| Layer | Technology | License | Rationale |
|---|---|---|---|
| Physics kernels | C++ / CombAero | yours | Thermodynamics, transport, acoustics — already exists |
| Python bindings | pybind11 | BSD | Already exists |
| Network solver | Python + SciPy `root` | BSD | Phase 1–2; upgrade to CasADi for exact Jacobians in Phase 3+ |
| Auto-differentiation | CasADi (Phase 3+) | LGPL | Exact sparse Jacobians through element residuals |
| Graph topology | NetworkX | BSD | Incidence matrix, traversal, validation |
| Regression elements | scikit-learn / PyTorch | BSD / BSD | Custom element pressure loss maps |
| API backend | FastAPI | MIT | Exposes solver as `POST /solve` |
| GUI frontend | React + React Flow | MIT | Node-graph editor |
| UI components | Tailwind + shadcn/ui | MIT | Controls, panels, results display |
| Packaging | Docker or PyInstaller | — | Standalone app or server deployment |

---

## Scope Separation: C++ vs Python

```
┌──────────────────────────────────────────────────────┐
│  PYTHON — network layer                              │
│                                                      │
│  NetworkSolver      Newton iteration, convergence    │
│  NetworkGraph       topology, validation, assembly   │
│  NetworkNode ABC    PlenumNode, MomentumChamber,     │
│                     JunctionNode, BoundaryNode       │
│  NetworkElement ABC OrificeElement, PipeElement,     │
│                     CombustionElement, HeatExchanger,│
│                     RegressionElement, MixingElement │
│  BoundaryCondition  PressureBC, MassFlowBC           │
│  MixtureState       P, T, m_dot, X[14]              │
└────────────────────────┬─────────────────────────────┘
                         │ pybind11
┌────────────────────────▼─────────────────────────────┐
│  C++ — CombAero physics kernels                      │
│                                                      │
│  h(), cp(), cv()         mixture enthalpy/heat cap   │
│  density()               ideal gas law               │
│  speed_of_sound()        isentropic                  │
│  viscosity(), prandtl()  transport properties        │
│  calc_T_from_h()         Newton inversion            │
│  mole_to_mass()          composition conversion      │
│  orifice_mdot()          incompressible Cd·A          │
│  pipe_flow()             Darcy-Weisbach + thermo       │
│  fanno_pipe()            compressible Fanno flow       │
│  cooling_correlations    Nusselt, friction factor    │
│  can_annular_eigenmodes  acoustics (Phase 4)         │
└──────────────────────────────────────────────────────┘
```

**Rule:** C++ owns property evaluation and correlations. Python owns equation
assembly, graph logic, and solver orchestration. The solver never calls C++
directly — always through element `residuals()`.

---

## Node Types

### Interior Nodes (unknowns solved by Newton)

| Node | Unknowns | Physics |
|---|---|---|
| `PlenumNode` | P, T, X | P_total = P_static (v ≈ 0), energy balance |
| `MomentumChamberNode` | P_static, P_total, v | Isentropic total/static relation, momentum flux, requires flow area |
| `JunctionNode` | P, T | Massless split/merge, Σṁ = 0 |

### Boundary Nodes (fixed values, contribute no residuals)

| BC | Fixed inputs | Free (solved) | Typical use |
|---|---|---|---|
| `PressureBoundary` | P_total, T_total, X | ṁ | Inlet reservoir, atmosphere exit |
| `MassFlowBoundary` | ṁ, T_total, X | P | Fuel injector, compressor delivery |

**Constraint:** every network requires at least one `PressureBoundary` for
well-posedness. The solver validates this at setup.

---

## Element Catalogue

### Flow Elements

| Element | Residual | CombAero call |
|---|---|---|
| `OrificeElement` | `ṁ = Cd·A·√(2ρΔP)` | `orifice_mdot()` or `orifice_flow_thermo()` |
| `PipeElement` (incompressible) | `ΔP = f·(L/D)·½ρv²` | `pipe_flow()` / `pipe_flow_rough()` |
| `PipeElement` (compressible) | Fanno adiabatic friction | `fanno_pipe()` / `fanno_pipe_rough()` |
| `MomentumChamberElement` | Momentum + area change | `density()`, `speed_of_sound()` |

### Thermal Elements

| Element | Residual | CombAero call |
|---|---|---|
| `HeatExchangerElement` | `Q = h·A·(T_aw − T_wall)`, solves T_wall | `htc_pipe()`, `T_adiabatic_wall()` from `stagnation.h` |
| `AdiabaticWallElement` | `Q = 0` | — |

### Reaction / Mixing Elements

| Element | Residual | CombAero call |
|---|---|---|
| `MixingElement` | Mass + species + enthalpy balance | `h()`, `mole_to_mass()` |
| `CombustionElement` | Atom balance + energy | `h()`, `calc_T_from_h()` — see `COMBUSTION_ELEMENTS.md` |

### Custom / Data-Driven Elements

| Element | Residual | Backend |
|---|---|---|
| `RegressionElement` | `ΔP = model(ṁ, T, ...)` | sklearn / PyTorch |

---

## MixtureState Dataclass

The common currency between all nodes and elements. This is the **canonical
definition** — `COMBUSTION_ELEMENTS.md` uses a simplified subset for Phase 1.

```python
from dataclasses import dataclass, field

@dataclass
class MixtureState:
    P: float                    # static pressure [Pa]
    P_total: float              # total pressure [Pa]  (= P for plenum nodes)
    T: float                    # static temperature [K]
    T_total: float              # total temperature [K]  (= T for plenum nodes)
    m_dot: float                # mass flow rate [kg/s]
    X: list[float]              # mole fractions [14 species]

    # Derived — computed on demand via CombAero
    def density(self) -> float: ...        # combaero.density(T, P, X)
    def enthalpy(self) -> float: ...       # combaero.h(T, X)
    def speed_of_sound(self) -> float: ... # combaero.speed_of_sound(T, X)
```

For a plenum node `P = P_total` and `T = T_total`. For a momentum chamber
node they differ and the isentropic relation closes the system.

**Phase 1** uses only `P` and `m_dot` (T and X are fixed, no composition
tracking). `P_total` and `T_total` are introduced in Phase 2.

---

## Static vs Stagnation Convention

`MixtureState` carries both static (`T`, `P`) and total (`T_total`, `P_total`).
The distinction matters for each part of the network:

| Context | Temperature to use | Notes |
|---|---|---|
| Element residuals (pipe ΔP, orifice ṁ) | static `T`, `P` | Darcy-Weisbach, Nusselt/Re/Pr all at bulk static conditions |
| Boundary conditions | total `T_total`, `P_total` | Reservoir / compressor delivery are naturally total conditions |
| Heat flux driving temperature | `T_aw` (adiabatic wall T) | **Not** T_static or T_total; use `combaero.T_adiabatic_wall()` |
| Combustion / mixing enthalpy balance | static `h(T)` + v²/2 = h₀ | Stagnation enthalpy h₀ is conserved, not static h |
| Isentropic elements (nozzle, momentum chamber) | convert via `combaero.T0_from_static()` / `P0_from_static()` | Requires known M at the element |

**Why T_aw for heat transfer?**

```
q = h_conv · (T_aw − T_wall)          ← correct for all Mach numbers

T_aw = T_static + r · v² / (2·cp)     r = Pr^(1/3) turbulent
                                           r = Pr^(1/2) laminar

T_static < T_aw < T_total  (since r < 1 for air)
```

- Using `T_static`: underestimates heat load at M > 0.3 (2–15% error for combustor liner / turbine cooling)
- Using `T_total`: overcorrects — `T_aw < T_total` always since r < 1
- Re and Pr for the Nusselt correlation are always evaluated at static `T` (correct)

**CombAero utilities** (`stagnation.h`, available as `combaero.*`):

```python
T_aw = combaero.T_adiabatic_wall(T_static, v, T, P, X)      # from velocity
T_aw = combaero.T_adiabatic_wall_mach(T_static, M, T, P, X) # from Mach number
T0   = combaero.T0_from_static(T, M, X)
P0   = combaero.P0_from_static(P, T, M, X)
T    = combaero.T_from_stagnation(T0, M, X)
P    = combaero.P_from_stagnation(P0, T0, M, X)
```

---

## Abstract Interfaces

```python
from abc import ABC, abstractmethod

class NetworkNode(ABC):
    @abstractmethod
    def unknowns(self) -> list[str]:
        """Names of unknowns this node contributes to the solver."""
        ...

    @abstractmethod
    def residuals(self, state: MixtureState) -> list[float]:
        """Returns zero when node equations are satisfied."""
        ...

class NetworkElement(ABC):
    @abstractmethod
    def residuals(self,
                  state_in: MixtureState,
                  state_out: MixtureState) -> list[float]:
        """Returns zero when element equations are satisfied."""
        ...

    @abstractmethod
    def n_equations(self) -> int:
        """Number of residual equations contributed."""
        ...
```

---

## Network JSON Format

```json
{
  "nodes": [
    {
      "id": "inlet",
      "type": "pressure_bc",
      "P_total": 500000,
      "T_total": 700,
      "X": {"N2": 0.76, "O2": 0.21, "AR": 0.01, "CO2": 0.02}
    },
    {
      "id": "plenum_1",
      "type": "plenum"
    },
    {
      "id": "chamber_1",
      "type": "momentum_chamber",
      "area": 0.05
    },
    {
      "id": "exit",
      "type": "pressure_bc",
      "P_total": 101325,
      "T_total": 300,
      "X": {"N2": 0.79, "O2": 0.21}
    }
  ],
  "edges": [
    {
      "id": "orif_1",
      "type": "orifice",
      "from": "inlet",
      "to": "plenum_1",
      "Cd": 0.7,
      "area": 0.001
    },
    {
      "id": "burner_1",
      "type": "combustion",
      "from": "plenum_1",
      "to": "chamber_1",
      "eta": 0.99,
      "delta_P_frac": 0.04,
      "fuel": {
        "type": "mass_flow_bc",
        "m_dot": 0.05,
        "T_total": 300,
        "X": {"CH4": 1.0}
      }
    },
    {
      "id": "orif_2",
      "type": "orifice",
      "from": "chamber_1",
      "to": "exit",
      "Cd": 0.65,
      "area": 0.002
    }
  ]
}
```

---

## Phased Roadmap

### Phase 1 — Isothermal Flow Network

**Goal:** solve P and ṁ throughout a network of orifices, pipes, and plenums.

- `MixtureState`: P and ṁ only (T and X fixed, no composition tracking)
- Nodes: `PlenumNode`, `JunctionNode`, `PressureBoundary`, `MassFlowBoundary`
- Elements: `OrificeElement`, `PipeElement`
- Solver: `scipy.optimize.root` with finite-difference Jacobian
- Graph: NetworkX incidence matrix for equation assembly
- Network described and validated via JSON
- **No GUI** — JSON in, JSON out
- Validation: orifices in series/parallel, known analytical solutions

### Phase 2 — Combustion and Composition

**Goal:** track mixture composition and temperature; add combustion element.

- `MixtureState` gains T and X[14]
- Nodes: `MomentumChamberNode` with total/static distinction
- Elements: `CombustionElement`, `MixingElement`, `RegressionElement`
- Segregated solve: flow first, then energy/composition propagation
- See `COMBUSTION_ELEMENTS.md` for combustion function design
- Validation: adiabatic flame temperature against Cantera

### Phase 3 — Heat Exchange

**Goal:** conjugate heat transfer between streams through walls.

- `HeatExchangerElement` with `T_wall` as additional unknown per segment
- Nusselt correlations via `cooling_correlations` (already in CombAero)
- Upgrade solver to CasADi for exact sparse Jacobians
- Fully coupled Newton solve over [P, ṁ, T, T_wall, X]
- BC types: fixed T_wall, fixed Q, adiabatic, convective both sides

### Phase 4 — GUI

**Goal:** visual network editor. Can develop in parallel with Phase 2–3.

- React + React Flow node-graph editor
- Nodes and edges map 1:1 to Python element classes
- FastAPI backend: `POST /solve` accepts network JSON, returns solution JSON
- Results overlay: colour-coded P, T, ṁ per element
- Export: network JSON, CSV results

### Phase 5 — Acoustics

**Goal:** linearise mean flow solution → acoustic transfer matrix propagation.

- Mean flow from Phase 1–3 provides the operating point
- Acoustic transfer matrices propagate on top of mean flow
- `can_annular_eigenmodes`, `annular_duct_modes` already in CombAero
- Assemble global transfer matrix from element acoustic matrices
- GUI: display mode shapes and frequency response

---

## Well-Posedness Rules

1. Every network requires ≥1 `PressureBoundary`
2. Number of pressure BCs + mass flow BCs = number of independent flow paths
3. Every node must be reachable from a boundary node
4. No isolated subgraphs

The solver validates all four at setup before attempting Newton iteration.
