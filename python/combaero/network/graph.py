from __future__ import annotations

import inspect
from typing import TYPE_CHECKING, Self

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
        Register a new flow element. Connects the element to its source and sink
        nodes, ensuring they exist in the network.
        """
        if element.id in self.elements:
            raise ValueError(f"Element '{element.id}' already exists in topology.")

        all_nodes = element.all_source_nodes() + element.all_sink_nodes()
        for nid in all_nodes:
            if nid not in self.nodes:
                raise ValueError(f"Unknown node '{nid}' for element '{element.id}'.")

        if len(all_nodes) != len(set(all_nodes)):
            raise ValueError(
                f"Element '{element.id}' has duplicate nodes among its source/sink nodes."
            )

        self.elements[element.id] = element

        for src in element.all_source_nodes():
            self._downstream_of_node[src].append(element.id)

        for snk in element.all_sink_nodes():
            self._upstream_of_node[snk].append(element.id)

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
        # Second pass: elements that inherit from other elements (e.g. TeeJunction.F_C from a
        # channel whose D was itself inherited) get a second chance after the first pass resolved
        # all direct values.
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
                if nid in elem.all_sink_nodes():
                    neighbours = elem.all_source_nodes()
                else:
                    neighbours = elem.all_sink_nodes()
                for nb in neighbours:
                    if nb not in visited:
                        queue.append(nb)
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

        from .components import LosslessConnectionElement, PressureLossElement

        all_lossless = all(isinstance(e, LosslessConnectionElement) for e in self.elements.values())
        if all_lossless:
            raise ValueError(
                "Network contains only lossless connection elements. "
                "There must be at least one pressure drop element to solve flow."
            )

        # Per-element specific validation
        for _, elem in self.elements.items():
            elem.validate()

        # PressureLossElement ambiguity: both endpoints have has_theta=True and no
        # explicit theta_source override.  Requires user disambiguation.
        for eid, elem in self.elements.items():
            if isinstance(elem, PressureLossElement) and elem._both_endpoints_theta:
                raise ValueError(
                    f"FlowNetwork validation failed: PressureLossElement '{eid}' has "
                    f"has_theta=True nodes on both endpoints "
                    f"('{elem.from_node}' and '{elem.to_node}'). "
                    f"Pass theta_source='<node_id>' to disambiguate."
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
            upstream = self._upstream_of_node[node_id]
            downstream = self._downstream_of_node[node_id]
            connected = upstream + downstream
            is_boundary = type(node).__name__ in (
                "PressureBoundary",
                "MassFlowBoundary",
                "WallNode",
            )

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

            # WallNode is a dead end - no flow may exit it.
            if type(node).__name__ == "WallNode" and downstream:
                raise ValueError(
                    f"FlowNetwork validation failed: WallNode '{node_id}' has downstream "
                    "connections. A wall is a closed end - connect only one element leading "
                    "into the wall, not out of it."
                )

            # Every interior node needs at least one element delivering flow in and one taking
            # flow out. All-upstream or all-downstream means mass cannot be conserved --
            # the most common cause is a tee junction wired with the wrong port as inlet.
            if not is_boundary and (not upstream or not downstream):
                direction = "upstream" if not upstream else "downstream"
                raise ValueError(
                    f"FlowNetwork validation failed: interior node '{node_id}' has no {direction} "
                    "connections -- flow has nowhere to go. This often means a tee junction port "
                    "is wired backwards (e.g. C on a branching tee should be the inlet)."
                )

            # MomentumChamberNode topology guard (#174). MCN's scalar closure
            # Pt = P + 0.5 rho v^2 only models a single bulk velocity, so it
            # cannot represent merging or splitting streams in the chamber
            # itself. >1 incoming or >1 outgoing edges puts MCN outside its
            # valid envelope; use MultiPortChamberElement (momentum-CV
            # junction with one MCN per port) for those topologies.
            if type(node).__name__ == "MomentumChamberNode" and (
                len(upstream) > 1 or len(downstream) > 1
            ):
                direction = "incoming" if len(upstream) > 1 else "outgoing"
                edges = upstream if len(upstream) > 1 else downstream
                raise ValueError(
                    f"FlowNetwork validation failed: MomentumChamberNode "
                    f"'{node_id}' has multiple {direction} edges ({sorted(edges)}). "
                    "MCN's scalar Pt = P + 0.5 rho v^2 closure cannot represent "
                    "merging or splitting streams. Use a MultiPortChamberElement "
                    "junction (one MCN per port) instead -- see "
                    "python/tests/test_multi_port_chamber.py for the pattern. "
                    "Tracked in #174."
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
    def from_dict(cls, data: dict) -> Self:
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
