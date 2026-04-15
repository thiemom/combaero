feat(gui): ensure thermal edge data consistency and refine typography

- Initializes new thermal edges with a default layers array in useStore.ts to resolve missing label issues.
- Implements fallback layer logic in ThermalEdge.tsx and Inspector.tsx to handle missing data.layers gracefully.
- Refines ThermalEdge styling:
  - Restores 2px orange border (matching boundary node stroke).
  - Isolates font-weight: bold for primary Q-label, normal (light) for diagnostic wall temperature.
  - Standardizes label size to 12px (text-xs).
- Updates wall temperature probe text to remove redundant "WALL" prefix.
