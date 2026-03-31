feat: enable thermal coupling for nodes and enhance heat flow telemetry

This commit enhances the CombAero network solver and GUI to support thermal
coupling between nodes (e.g., MomentumChambers) and elements.

Key Changes:
- Modified `FlowNetwork.add_wall` to allow `Node` IDs as endpoints.
- Implemented `htc_and_T` in `MomentumChamberNode` for convection calculations.
- Integrated `edge_results` ($Q, T_{wall}$) into the `NetworkSolver` output.
- Updated FastAPI backend to relay heat flow telemetry to the GUI.
- Enhanced `ThermalEdge` in React Flow to display heat flux labels (W/kW).
- Restored thermal handles on `MomentumChamberNode` and fixed SVG accessibility.
- Automatic calculation of pipe convective areas in GUI graph builder.

Closes thermal coupling visualization issues.
