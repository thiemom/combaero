from __future__ import annotations

import inspect
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .components import NetworkElement, NetworkNode, ThermalWall


class FlowNetwork:
    """
    Manages the topology of the flow network graph, storing nodes and elements.
    Provides methods for elements to discover their upstream/downstream neighbors.
    """

    def __init__(self) -> None:
        self.nodes: dict[str, NetworkNode] = {}
        self.elements: dict[str, NetworkElement] = {}

        # Thermal coupling walls dictionary
        self.walls: dict[str, ThermalWall] = {}
        self.thermal_coupling_enabled: bool = True

        # Directed graph mapping: node_id -> list of element_ids
        self._upstream_of_node: dict[str, list[str]] = {}
        self._downstream_of_node: dict[str, list[str]] = {}

    def add_node(self, node: NetworkNode) -> None:
        """Register a new node in the network."""
        if node.id in self.nodes:
            raise ValueError(f"Node '{node.id}' already exists in topology.")
        self.nodes[node.id] = node
        self._upstream_of_node[node.id] = []
        self._downstream_of_node[node.id] = []

    def add_element(self, element: NetworkElement) -> None:
        """
        Register a new flow element. Connects the element to its from_node and
        to_node, ensuring they exist in the network.
        """
        if element.id in self.elements:
            raise ValueError(f"Element '{element.id}' already exists in topology.")

        if element.from_node not in self.nodes:
            raise ValueError(f"Unknown from_node '{element.from_node}' for element {element.id}.")

        if element.to_node not in self.nodes:
            raise ValueError(f"Unknown to_node '{element.to_node}' for element {element.id}.")

        if element.from_node == element.to_node:
            raise ValueError(
                f"Element '{element.id}' has from_node == to_node ('{element.from_node}'). "
                "Self-loops are not permitted."
            )

        self.elements[element.id] = element

        # Element 'from_node' -> implies element is downstream of the from_node
        self._downstream_of_node[element.from_node].append(element.id)

        # Element 'to_node' -> implies element is upstream of the to_node
        self._upstream_of_node[element.to_node].append(element.id)

    def add_wall(self, wall: ThermalWall) -> None:
        """Register a thermal coupling wall between two elements.

        Parameters
        ----------
        wall : ThermalWall
            Wall connection defining thermal coupling between two elements.

        Raises
        ------
        ValueError
            If wall ID already exists or if element_a/element_b don't exist.
        """
        if wall.id in self.walls:
            raise ValueError(f"Wall '{wall.id}' already exists in network.")

        if wall.element_a not in self.elements and wall.element_a not in self.nodes:
            raise ValueError(f"Unknown element_a/node_a '{wall.element_a}' for wall '{wall.id}'.")

        if wall.element_b not in self.elements and wall.element_b not in self.nodes:
            raise ValueError(f"Unknown element_b/node_b '{wall.element_b}' for wall '{wall.id}'.")

        self.walls[wall.id] = wall

    def get_upstream_elements(self, node_id: str) -> list[NetworkElement]:
        """Return all elements feeding into the specified node."""
        if node_id not in self.nodes:
            raise KeyError(f"Node '{node_id}' not found.")
        return [self.elements[eid] for eid in self._upstream_of_node[node_id]]

    def get_downstream_elements(self, node_id: str) -> list[NetworkElement]:
        """Return all elements flowing out of the specified node."""
        if node_id not in self.nodes:
            raise KeyError(f"Node '{node_id}' not found.")
        return [self.elements[eid] for eid in self._downstream_of_node[node_id]]

    def resolve_all_topology(self) -> None:
        """
        Trigger the auto-discovery phase where all components inspect the graph
        to extract their needed operational parameters (e.g. Orifice discovering
        upstream pipe diameter).
        """
        for node in self.nodes.values():
            node.resolve_topology(self)
        for element in self.elements.values():
            element.resolve_topology(self)

    def _reachable_from(self, start_id: str) -> set[str]:
        """BFS over undirected adjacency to find all node IDs reachable from start_id."""
        visited: set[str] = set()
        queue = [start_id]
        while queue:
            nid = queue.pop()
            if nid in visited:
                continue
            visited.add(nid)
            for eid in self._upstream_of_node[nid] + self._downstream_of_node[nid]:
                elem = self.elements[eid]
                neighbour = elem.from_node if elem.to_node == nid else elem.to_node
                if neighbour not in visited:
                    queue.append(neighbour)
        return visited

    def validate(self) -> None:
        """
        Validate the network topology.
        - Must have at least one node and one element.
        - Must have at least one PressureBoundary.
        - Must contain at least one element with losses (to prevent infinite flow).
        - No isolated subgraphs (all nodes reachable from first node).
        - Every interior node must be connected to at least one element.
        """
        from .components import PressureBoundary

        if not self.nodes:
            raise ValueError("FlowNetwork validation failed: the network contains no nodes.")

        if len(self.elements) == 0:
            raise ValueError("FlowNetwork validation failed: the network contains no elements.")

        has_pressure_bc = any(isinstance(node, PressureBoundary) for node in self.nodes.values())
        if not has_pressure_bc:
            raise ValueError(
                "FlowNetwork validation failed: the network must contain at least "
                "one PressureBoundary to serve as a reference."
            )

        from .components import LosslessConnectionElement

        all_lossless = all(isinstance(e, LosslessConnectionElement) for e in self.elements.values())
        if all_lossless:
            raise ValueError(
                "Network contains only lossless connection elements. "
                "There must be at least one pressure drop element to solve flow."
            )

        # Isolated subgraph check: all nodes must be reachable from *some* PressureBoundary.
        # This prevents islands that have elements but no reference pressure.
        reachable_from_bc = set()
        for nid, node in self.nodes.items():
            if isinstance(node, PressureBoundary):
                reachable_from_bc.update(self._reachable_from(nid))

        isolated = set(self.nodes) - reachable_from_bc
        if isolated:
            raise ValueError(
                f"FlowNetwork validation failed: isolated node(s) detected: {sorted(isolated)}. "
                "Every node must be reachable from a PressureBoundary."
            )

        # Every interior node must have at least 2 connected elements to satisfy continuity (in and out)
        # Boundary nodes only need at least 1.
        for node_id, node in self.nodes.items():
            connected = self._upstream_of_node[node_id] + self._downstream_of_node[node_id]
            is_boundary = type(node).__name__ in ("PressureBoundary", "MassFlowBoundary")

            if not connected:
                raise ValueError(
                    f"FlowNetwork validation failed: node '{node_id}' has no connected elements."
                )

            if not is_boundary and len(connected) < 2:
                raise ValueError(
                    f"FlowNetwork validation failed: interior node '{node_id}' has only {len(connected)} "
                    "connected element(s). Interior nodes must have at least 2 connections to satisfy "
                    "mass conservation."
                )

    def to_dict(self) -> dict:
        """
        Serialize the entire network topology and physics settings into a declarative dictionary suitable for JSON.
        Uses inspect.signature to dynamically capture constructor arguments from component instances.
        """
        nodes_data = {}
        for nid, node in self.nodes.items():
            sig = inspect.signature(node.__init__)
            # Extract arguments identically matching the __init__ signature keys
            args = {
                k: getattr(node, k, getattr(node, f"_{k}", None))
                for k in sig.parameters
                if k != "self" and k != "args" and k != "kwargs"
            }
            nodes_data[nid] = {"type": type(node).__name__, "kwargs": args}

        elements_data = {}
        for eid, element in self.elements.items():
            sig = inspect.signature(element.__init__)
            args = {
                k: getattr(element, k, getattr(element, f"_{k}", None))
                for k in sig.parameters
                if k != "self" and k != "args" and k != "kwargs"
            }
            elements_data[eid] = {"type": type(element).__name__, "kwargs": args}

        walls_data = {}
        for wid, wall in self.walls.items():
            sig = inspect.signature(wall.__init__)
            args = {
                k: getattr(wall, k, getattr(wall, f"_{k}", None))
                for k in sig.parameters
                if k != "self" and k != "args" and k != "kwargs"
            }
            walls_data[wid] = {"type": type(wall).__name__, "kwargs": args}

        return {
            "nodes": nodes_data,
            "elements": elements_data,
            "walls": walls_data,
            "thermal_coupling_enabled": self.thermal_coupling_enabled,
        }

    @classmethod
    def from_dict(cls, data: dict) -> FlowNetwork:
        """
        Deserialize a flow network from a declarative dictionary generated by to_dict().
        """
        from . import components

        network = cls()

        # Restore thermal coupling toggle
        network.thermal_coupling_enabled = data.get("thermal_coupling_enabled", True)

        for _nid, n_data in data.get("nodes", {}).items():
            NodeClass = getattr(components, n_data["type"])
            network.add_node(NodeClass(**n_data["kwargs"]))

        for _eid, e_data in data.get("elements", {}).items():
            ElemClass = getattr(components, e_data["type"])
            network.add_element(ElemClass(**e_data["kwargs"]))

        for _wid, w_data in data.get("walls", {}).items():
            WallClass = getattr(components, w_data["type"])
            network.add_wall(WallClass(**w_data["kwargs"]))

        return network
