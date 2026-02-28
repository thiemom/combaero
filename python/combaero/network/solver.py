import numpy as np
from scipy.optimize import root

from .components import (
    MassFlowBoundary,
    MixtureState,
    NetworkNode,
    PressureBoundary,
)
from .graph import FlowNetwork


class NetworkSolver:
    """
    Orchestrates the numerical solution of a fluid flow network using scipy.optimize.root.
    Translates the graph topology into a flat array of unknowns and residuals.
    """

    def __init__(self, network: FlowNetwork):
        self.network = network
        self.unknown_names: list[str] = []
        # Mapping from node/element ID to index in the flattened x array
        self._unknown_indices: dict[str, list[int]] = {}

    def _build_x0(self) -> np.ndarray:
        """
        Constructs the initial guess vector by gathering all unknowns.
        Returns a flat numpy array.
        """
        self.unknown_names = []
        self._unknown_indices = {}
        x0_list = []

        # Iterate through interior nodes and gather unknowns
        for node_id, node in self.network.nodes.items():
            if not isinstance(node, (PressureBoundary, MassFlowBoundary)):
                unknowns = node.unknowns()
                if unknowns:
                    start_idx = len(x0_list)
                    for unk in unknowns:
                        self.unknown_names.append(unk)
                        # Very simple initial guess logic based on variable name
                        if unk.endswith(".P") or unk.endswith(".P_total"):
                            x0_list.append(101325.0)  # Default 1 atm
                        elif unk.endswith(".T") or unk.endswith(".T_total"):
                            x0_list.append(300.0)
                        else:
                            x0_list.append(1.0)
                    self._unknown_indices[node_id] = list(range(start_idx, len(x0_list)))

        # Iterate through elements and gather unknowns (m_dot)
        for elem_id, element in self.network.elements.items():
            unknowns = element.unknowns()
            if unknowns:
                start_idx = len(x0_list)
                for unk in unknowns:
                    self.unknown_names.append(unk)
                    if unk.endswith(".m_dot"):
                        x0_list.append(0.1)  # Default 0.1 kg/s
                    else:
                        x0_list.append(1.0)
                self._unknown_indices[elem_id] = list(range(start_idx, len(x0_list)))

        return np.array(x0_list, dtype=float)

    def _get_node_state(self, node: NetworkNode, x: np.ndarray) -> MixtureState:
        """Constructs a MixtureState for a given node based on the current solver vector x."""
        # species_names order: N2=0, O2=1, AR=2, CO2=3, H2O=4, ...
        # species names are exposed to python in python/combaero/species.py
        X_air = [0.0] * 14
        X_air[0] = 0.79  # N2
        X_air[1] = 0.21  # O2

        # Determine boundaries
        if isinstance(node, PressureBoundary):
            # For PressureBoundary, defaults.
            # We first check if the test suite provided an initial_guess dict for external states
            guess = getattr(node, "initial_guess", {})

            P_tot = guess.get(f"{node.id}.P_total", getattr(node, "P_total", 101325.0))
            P_stat = guess.get(f"{node.id}.P", P_tot)
            T_tot = guess.get(f"{node.id}.T_total", getattr(node, "T_total", 300.0))
            T_stat = guess.get(f"{node.id}.T", T_tot)
            X = getattr(node, "X", X_air)
            return MixtureState(P=P_stat, P_total=P_tot, T=T_stat, T_total=T_tot, m_dot=0.0, X=X)

        if isinstance(node, MassFlowBoundary):
            guess = getattr(node, "initial_guess", {})
            m = getattr(node, "m_dot", 0.1)
            T_tot = guess.get(f"{node.id}.T_total", getattr(node, "T_total", 300.0))
            T_stat = guess.get(f"{node.id}.T", T_tot)
            X = getattr(node, "X", X_air)
            return MixtureState(P=101325.0, P_total=101325.0, T=T_stat, T_total=T_tot, m_dot=m, X=X)

        # For interior nodes, unpack from x
        indices = self._unknown_indices.get(node.id, [])
        unknowns = node.unknowns()

        # Default state
        state_dict = {
            "P": 101325.0,
            "P_total": 101325.0,
            "T": 300.0,
            "T_total": 300.0,
            "m_dot": 0.0,
            "X": X_air,
        }

        for i, unk in zip(indices, unknowns, strict=False):
            var_name = unk.split(".")[-1]
            if var_name in state_dict:
                state_dict[var_name] = x[i]

        return MixtureState(**state_dict)

    def _residuals(self, x: np.ndarray) -> np.ndarray:
        """
        Evaluates the global residual vector given the current unknown vector x.
        """
        res = []

        # 1. Node Residuals (P_total constraint + Mass Conservation)
        for node_id, node in self.network.nodes.items():
            if isinstance(node, (PressureBoundary, MassFlowBoundary)):
                continue

            # Node thermodynamic constraints
            state = self._get_node_state(node, x)
            node_res = node.residuals(state)
            res.extend(node_res)

            # Mass Conservation: Sum(m_dot_in) - Sum(m_dot_out) = 0
            # Includes MassFlowBoundary contributions from adjacent boundary nodes.
            upstream_elems = self.network.get_upstream_elements(node_id)
            downstream_elems = self.network.get_downstream_elements(node_id)

            m_dot_in = 0.0
            m_dot_out = 0.0

            for elem in upstream_elems:
                indices = self._unknown_indices.get(elem.id)
                if indices:
                    m_dot_in += x[indices[0]]
                else:
                    # Zero-DOF element: read m_dot from its from_node if it is a boundary
                    from_node = self.network.nodes[elem.from_node]
                    # Upstream boundary element: its flow goes into this node, so it adds to m_dot_in
                    m_dot_in += getattr(from_node, "m_dot", 0.0)

            for elem in downstream_elems:
                indices = self._unknown_indices.get(elem.id)
                if indices:
                    m_dot_out += x[indices[0]]
                else:
                    # Downstream boundary element: its flow leaves this node, so it adds to m_dot_out
                    to_node = self.network.nodes[elem.to_node]
                    m_dot_out += getattr(to_node, "m_dot", 0.0)

            res.append(m_dot_in - m_dot_out)

        # 2. Element Residuals (Momentum Conservation)
        for elem_id, element in self.network.elements.items():
            state_in = self._get_node_state(self.network.nodes[element.from_node], x)
            state_out = self._get_node_state(self.network.nodes[element.to_node], x)

            # Insert the element's current m_dot from x (if it has a DOF)
            indices = self._unknown_indices.get(elem_id)
            if indices:
                m_dot = x[indices[0]]
                state_in.m_dot = m_dot
                state_out.m_dot = m_dot

            elem_res = element.residuals(state_in, state_out)
            res.extend(elem_res)

        return np.array(res, dtype=float)

    def solve(self, method="hybr") -> dict[str, float]:
        """
        Executes the root finding algorithm to solve the network.
        Returns a dictionary mapping unknown names to their solved values.
        """
        # Ensure topology is resolved before solving
        self.network.resolve_all_topology()
        self.network.validate()

        x0 = self._build_x0()

        # Dimension check
        res0 = self._residuals(x0)
        if len(x0) != len(res0):
            raise ValueError(f"System not square: {len(x0)} unknowns vs {len(res0)} equations.")

        # Solve
        sol = root(self._residuals, x0, method=method)

        if not sol.success:
            import warnings

            warnings.warn(
                f"NetworkSolver did not converge: {sol.message} "
                f"(|F|={float(np.linalg.norm(sol.fun)):.3e}). "
                "Returning best iterate.",
                stacklevel=2,
            )

        return dict(zip(self.unknown_names, sol.x, strict=False))
