# GUI Technical Reference

This document provides the technical specification for the JSON schemas and data flow used by the CombAero Network Designer. This is intended for developers extending the GUI and agents performing automated network generation or inspection.

## Data Flow Architecture

The GUI follows a decoupled architecture:
1.  **Frontend (React Flow)**: Manages the visual graph, nodes, and edges.
2.  **State Persistence**: The graph is serialized to a JSON format (see `NetworkGraphSchema` below).
3.  **Backend (FastAPI)**: Receives the JSON, maps React Flow nodes/edges to C++ physics objects, and invokes the `NetworkSolver`.
4.  **Results Propagation**: The solver returns a `SolveResponse` containing the converged state of every node and element.

## Network Schema

The primary data structure for saving/loading networks is the `NetworkGraphSchema`.

### Root Structure
```json
{
  "nodes": [],           // List of ReactFlowNode objects
  "edges": [],           // List of ReactFlowEdge objects
  "solver_settings": {}  // SolverSettings object
}
```

### Node Data Schemas
Each node in the `nodes` list carries a `data` object corresponding to its type.

#### Pressure Boundary (`pressure_boundary`)
Used for inlets/outlets with fixed total conditions.
- `Pt`: Total Pressure [Pa] (Default: 101325.0)
- `Tt`: Total Temperature [K] (Default: 300.0)
- `composition`: Composition object (Mode: mole/mass, Source: dry_air/fuel/custom)

#### Plenum Node (`plenum`)
Internal mixing node where pressure and temperature are solved.
- `initial_guess`: `{"P": float, "T": float}` (Optional)

#### Combustor Node (`combustor`)
Internal node with heat release logic.
- `method`: "complete" or "equilibrium"
- `area`: Cross-sectional area [m²]

### Element Data Schemas
Elements are defined by edges in the graph but carry their own physical parameters.

#### Orifice Element (`orifice`)
- `diameter`: Orifice diameter [m]
- `Cd`: Discharge coefficient (if using "fixed" correlation)
- `correlation`: "ReaderHarrisGallagher", "Stolz", "Miller", "ThickPlate", "RoundedEntry", or "fixed".
- `regime`: "default", "incompressible", or "compressible".

#### Channel Element (`channel`)
- `L`: Length [m]
- `D`: Hydraulic diameter [m]
- `roughness`: Surface roughness [m]
- `surface`: Surface model (smooth, ribbed, dimpled, pin_fin, impingement).

## Solver Settings
```json
{
  "global_regime": "incompressible" | "compressible",
  "init_strategy": "default" | "incompressible_warmstart" | "homotopy",
  "method": "hybr" | "lm" | "broyden1" | ...,
  "timeout": 180.0
}
```

## Solve Results Schema
After a solve, the backend returns a `SolveResponse`.

- `success`: boolean
- `message`: Error message if failed.
- `node_results`: Map of node IDs to results.
  - `state`: `{"P": float, "T": float, "Pt": float, "Tt": float, "Y": list[float], ...}`
- `element_results`: Map of element IDs to results.
  - `m_dot`: Mass flow rate [kg/s]
