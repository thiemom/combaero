from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .components import NetworkElement, NetworkNode


class FlowNetwork:
    """
    Manages the topology of the flow network graph, storing nodes and elements.
    Provides methods for elements to discover their upstream/downstream neighbors.
    """

    def __init__(self) -> None:
        self.nodes: dict[str, NetworkNode] = {}
        self.elements: dict[str, NetworkElement] = {}

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

        self.elements[element.id] = element

        # Element 'from_node' -> implies element is downstream of the from_node
        self._downstream_of_node[element.from_node].append(element.id)

        # Element 'to_node' -> implies element is upstream of the to_node
        self._upstream_of_node[element.to_node].append(element.id)

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

    def validate(self) -> None:
        """
        Validate the network topology.
        - Must have at least one PressureBoundary.
        - Must contain at least one element with losses (to prevent infinite flow).
        - No isolated subgraphs (implicitly checked if Newton solver runs, but good practice).
        """
        from .components import LosslessConnectionElement, PressureBoundary

        has_pressure_bc = any(isinstance(node, PressureBoundary) for node in self.nodes.values())
        if not has_pressure_bc:
            raise ValueError("Network must contain at least one PressureBoundary node.")

        all_lossless = all(isinstance(e, LosslessConnectionElement) for e in self.elements.values())
        if self.elements and all_lossless:
            raise ValueError(
                "Network contains only lossless connection elements. "
                "There must be at least one pressure drop element to solve flow."
            )
