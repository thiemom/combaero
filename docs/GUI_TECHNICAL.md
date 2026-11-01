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
Internal mixing node; its pressure is solved (P and Pt are the unknowns, and
a plenum's gas is at rest, so Pt = P). Temperature is not an unknown: it is
propagated from the inflows at every evaluation.
- `initial_guess`: `{"P": float}` (Optional; seeds both P and Pt). The
  Inspector offers P only; a `T` key in a hand-written file has nothing to seed.

#### Combustor Node (`combustor`)
Internal node with heat release logic.
- `method`: "complete" or "equilibrium"
- `area`: Cross-sectional area [m²]

#### Heat Source (`heat_source`)
A heat source or sink (#792): an `EnergyBoundary` on the plenum, momentum
chamber or combustor that its ONE edge (handle `heat-out`, edge
`data.type = "heat"`) points at. Not a flow path. The backend refuses an
unwired source, an edge to an element or boundary, and edges to two nodes.
- `mode`: `"Q"` (default), `"fraction"` or `"Tt_target"`
- `Q`: Heat flow [W], > 0 heats the fluid (Default: 1000.0)
- `fraction`: Fraction of the inflow's sensible enthalpy [-] (Default: 0.01)
- `Tt_target`: The node's outlet TOTAL temperature [K]; the heat is then
  solved for (Default: 600.0)

Its result (`node_results[<id>]`) carries `Q` [W], the heat it applied in
every mode (NaN when no through-flow determines it), `Tt_target` in that mode,
`heated_node`, and the heated node's `state`.

#### Fuel Boundary (`fuel_boundary`)
A fuel controller (#892): a `FuelBoundary` whose fuel flow is SOLVED so that
the combustor named by `combustor_id` is at `Tt_target`. It only supplies
(one outlet handle, `flow-source`). The canvas draws a dashed yellow link to
its combustor, derived from `combustor_id` like a probe's from `target_id`
(never stored as an edge); drag from its link handle (`probe-out`) onto a
combustor to choose one, or connect the fuel straight to a combustor, which
links it when it names none. Changing the combustor is a structural edit.
The backend refuses a fuel boundary with no `combustor_id`.
- `combustor_id`: the combustor it controls (required)
- `Tt_target`: Combustor exit total temperature [K] (Default: 1500.0)
- `branch`: `"lean"` (default) or `"rich"`
- `Tt`: Fuel total temperature [K] (Default: 300.0)
- `composition`: as a boundary's (Default: source `fuel`, CH4)
- `m_dot_guess`: optional start [kg/s]; otherwise sized on the target

Its result `state` carries `m_dot`, `phi`, `Tt_target`, `Tt_without_fuel`
and `rich_branch`.

### Element Data Schemas
Elements are defined by edges in the graph but carry their own physical parameters.

#### Orifice Element (`orifice`)
- `diameter`: Orifice diameter [m]
- `Cd`: Discharge coefficient (if using "fixed" correlation)
- `correlation`: "ReaderHarrisGallagher", "Stolz", "IdelchikThick", "IdelchikBeveled", "IdelchikRounded", "McGreehanSchotsch", "Lichtarowicz", or "fixed".
  - "Stolz" is deprecated (#876), removed in 0.9.0. A saved file with it still loads and solves; the Inspector labels it deprecated and `/solve` returns a note for each such orifice in `SolveResponse.deprecations`, shown on the solver status badge.
- `regime`: "default", "incompressible", or "compressible".

#### Compressor and Turbine (`compressor`, `turbine`)
A `CompressorElement` / `TurbineElement` at a fixed total-to-total pressure
ratio (#893): Pt_out = PR x Pt_in (compressor) or Pt_in / PR (turbine); the
network sets the flow. One inlet (`flow-target`) and one outlet
(`flow-source`); no surface, so no thermal handles.
- `pressure_ratio`: >= 1 (Default: 4 compressor, 3 turbine)
- `efficiency_type`: `"polytropic"` (default) or `"isentropic"`
- `efficiency`: in (0, 1] (Default: 0.88 compressor, 0.9 turbine)

Its result carries `power` [W] > 0 (absorbed / delivered), `dh` [J/kg] signed
into the flow, `Tt_in`, `Tt_out`, `Tt_out_isentropic`, `eta_poly`, `eta_is`,
`pressure_ratio`. A reversed-flow root is refused; a machine between two
pressure boundaries is refused (nothing sets its flow).

#### Channel Element (`channel`)
- `L`: Length [m]
- `D`: Hydraulic diameter [m]
- `roughness`: Surface roughness [m]
- `surface`: Surface model (smooth, ribbed, dimpled, pin_fin, impingement).

#### Ejector Element (`ejector`)
Critical-mode supersonic ejector (3 ports: `primary` + `secondary` inlets,
one `outlet`; handle ids `port-{primary,secondary,outlet}-{target,source}`).
- `throat_area`: Primary nozzle throat area A_t [m^2].
- `nozzle_exit_area`: Primary nozzle exit area A_p1 [m^2] (must exceed throat).
- `mixing_area`: Constant-area mixing-section area A_3 [m^2] (must exceed nozzle exit).
- `recovery_efficiency`: Multiplies the lossless mixed stagnation pressure to
  give the critical back pressure P_c* (default 1.0 = no artificial loss).

Diagnostics returned on solve include `entrainment_ratio`,
`P_c_star` (critical back pressure), and `gamma`.

## Solver Settings
```json
{
  "global_regime": "incompressible" | "compressible",
  "init_strategy": "default" | "analytical_pt_prop" | "outletref_warmstart" | "homotopy" | "continuation",
  "method": "hybr" | "lm" | "broyden1" | ...,
  "timeout": 30.0
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

## CSV Export
`POST /export` solves the network and returns `NetworkResult.to_dataframe()`
as CSV (ASCII, UTF-8 without BOM). The format is one row per entity:

| Columns | Rows | Content |
|---|---|---|
| `type`, `kind`, `is_boundary`, `combustion_method`, `id`, `label` | all | what the row is |
| `success`, `converged`, `consistent`, `outcome`, `final_norm` | all | the solve's status |
| `<key> [<unit>]` (`T [K]`, `m_dot [kg/s]`, `Y[N2] [-]`, ...) | result rows | results (#774) |
| `input:<key> [<unit>]` (`input:L [m]`, `input:composition.source`) | entity rows | the case's inputs, nested fields dotted (#784) |
| `initial_guess:<key>` | entity rows | the entity's initial guess |
| `solver:<key>` (`solver:method`, `solver:timeout`, ...) | network row | how the case was solved |

Result rows are `node`, `element` and `thermal_wall`, as before. Since 0.8.0
the input columns follow every result column, and two row types are added:
one `network` row (id `_network`) with the network-wide settings
(`input:global_regime`, ambient conditions, multipliers) and the solver
options, and an `edge` row for an edge with inputs and no result row. A flow
edge's inputs go on its auto-linked element row. Every numeric input carries
its unit (`INPUT_UNIT_MAP` in `gui/backend/units.py`); text inputs are bare.
A text cell starting with `= + - @`, a tab or a carriage return gets a
leading `'`.

`gui.backend.results.load_results("combaero_results.csv")` reads the file back
into the one-row-per-case results table, inputs included; see
[API_PYTHON.md](API_PYTHON.md#results-table-one-row-per-case).
