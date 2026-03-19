import time
import warnings
from typing import Any

import numpy as np
from scipy.optimize import root

import combaero as cb

from .components import (
    EnergyBoundary,
    MassFlowBoundary,
    MixtureState,
    NetworkNode,
    PressureBoundary,
)
from .graph import FlowNetwork


class SolverTimeoutError(Exception):
    """Raised when the solver exceeds the specified wall-clock timeout."""

    pass


class NetworkSolver:
    """
    Orchestrates the numerical solution of a fluid flow network using scipy.optimize.root.
    Translates the graph topology into a flat array of unknowns and residuals.
    Nodes define their Pressure (P) and Total Pressure (P_total) as unknowns,
    while Temperature (T) and Composition (Y) are forward-propagated using
    topological ordering and automatic differentiation (chain-rule relay).
    """

    SUPPORTED_METHODS = [
        "hybr",
        "lm",
        "broyden1",
        "broyden2",
        "anderson",
        "linearmixing",
        "diagbroyden",
        "excitingmixing",
        "krylov",
        "df-sane",
    ]

    def __init__(self, network: FlowNetwork):
        self.network = network
        self.unknown_names: list[str] = []
        # Mapping from node/element ID to index in the flattened x array
        self._unknown_indices: dict[str, list[int]] = {}
        # Mapping from unknown name to global index
        self._name_to_index: dict[str, int] = {}
        # Canonical default composition (Mass Fractions)
        self._default_Y: list[float] = list(cb.mole_to_mass(cb.standard_dry_air_composition()))
        self._topological_order: list[str] = []
        self._derived_states: dict[str, tuple[float, list[float], Any]] = {}

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
            if not isinstance(node, PressureBoundary):
                unknowns = node.unknowns()
                if unknowns:
                    start_idx = len(x0_list)
                    for unk in unknowns:
                        self.unknown_names.append(unk)
                        guess_dict = getattr(node, "initial_guess", {})
                        if unk in guess_dict:
                            x0_list.append(guess_dict[unk])
                        elif unk.endswith(".P") or unk.endswith(".P_total"):
                            x0_list.append(101325.0)  # Default 1 atm
                        elif unk.endswith(".T") or unk.endswith(".T_total"):
                            x0_list.append(300.0)
                        elif ".Y[" in unk:
                            idx = int(unk.split(".Y[")[1].replace("]", ""))
                            x0_list.append(self._default_Y[idx])
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
                    guess_dict = getattr(element, "initial_guess", {})
                    if unk in guess_dict:
                        x0_list.append(guess_dict[unk])
                    elif unk.endswith(".m_dot"):
                        x0_list.append(0.1)  # Default 0.1 kg/s
                    else:
                        x0_list.append(1.0)
                self._unknown_indices[elem_id] = list(range(start_idx, len(x0_list)))

        # Build name to index mapping
        self._name_to_index = {name: i for i, name in enumerate(self.unknown_names)}

        # Pre-compute topological order for forward state propagation
        self._topological_order = self._compute_topological_order()

        # Set up flow-dependent wall EnergyBoundaries for thermal coupling
        self._wall_ebs: dict[str, EnergyBoundary] = {}
        if self.network.thermal_coupling_enabled and self.network.walls:
            self._setup_wall_energy_boundaries()

        return np.array(x0_list, dtype=float)

    def _compute_topological_order(self) -> list[str]:
        """Simple topological sort starting from boundary nodes."""
        order = []
        visited = set()
        queue = []
        # Boundaries are starting points
        for nid, node in self.network.nodes.items():
            if isinstance(node, (PressureBoundary, MassFlowBoundary)):
                queue.append(nid)
                visited.add(nid)

        while queue:
            nid = queue.pop(0)
            order.append(nid)
            for elem in self.network.get_downstream_elements(nid):
                to_id = elem.to_node
                if to_id not in visited:
                    up_elems = self.network.get_upstream_elements(to_id)
                    # A node is ready if all its upstream 'from' nodes are visited
                    if all(e.from_node in visited for e in up_elems):
                        visited.add(to_id)
                        queue.append(to_id)

        # Cleanup: append any unvisited (cycles or islands)
        for nid in self.network.nodes:
            if nid not in visited:
                order.append(nid)
        return order

    def _setup_wall_energy_boundaries(self) -> None:
        """Create flow-dependent EnergyBoundary objects for wall-coupled nodes.

        For each node that has upstream elements participating in a WallConnection,
        a dedicated EnergyBoundary is created and attached to the node.  Its Q is
        updated every Newton iteration inside ``_propagate_states``.
        """
        # Build a mapping: node_id -> set of wall IDs that affect it
        affected_nodes: dict[str, set[str]] = {}
        for wall in self.network.walls.values():
            elem_a = self.network.elements[wall.element_a]
            elem_b = self.network.elements[wall.element_b]

            # Skip walls where both sides feed the same node (net Q = 0)
            if elem_a.to_node == elem_b.to_node:
                continue

            affected_nodes.setdefault(elem_a.to_node, set()).add(wall.id)
            affected_nodes.setdefault(elem_b.to_node, set()).add(wall.id)

        for nid in affected_nodes:
            node = self.network.nodes[nid]
            if not hasattr(node, "energy_boundaries"):
                continue
            eb = EnergyBoundary(id=f"_wall_{nid}", Q=0.0)
            self._wall_ebs[nid] = eb
            node.add_energy_boundary(eb)

    def _evaluate_walls_for_node(
        self,
        nid: str,
        up_elems: list,
        x: np.ndarray,
    ) -> list[tuple]:
        """Evaluate all wall couplings affecting *nid* and return contributions.

        Returns a list of ``(wall, Q_inject, wall_result, elem_a_id, elem_b_id,
        is_side_a, ch_result_a, ch_result_b)`` tuples.  Also updates the wall
        ``EnergyBoundary`` Q on the node so that ``compute_derived_state`` picks
        it up.
        """
        wall_contributions: list[tuple] = []
        processed_walls: set[str] = set()
        up_elem_ids = {e.id for e in up_elems}

        for wall in self.network.walls.values():
            if wall.id in processed_walls:
                continue

            # Does this wall involve an upstream element of nid?
            elem_a_id = wall.element_a
            elem_b_id = wall.element_b
            if elem_a_id not in up_elem_ids and elem_b_id not in up_elem_ids:
                continue

            elem_a = self.network.elements[elem_a_id]
            elem_b = self.network.elements[elem_b_id]

            # Skip if both elements feed the same node (net Q = 0)
            if elem_a.to_node == elem_b.to_node:
                continue

            # Determine which side feeds the current node
            is_side_a = elem_a_id in up_elem_ids

            # Get states for both sides (previous-iteration fallback for back-edges)
            state_a = self._get_node_state_with_prev(self.network.nodes[elem_a.from_node], x)
            state_b = self._get_node_state_with_prev(self.network.nodes[elem_b.from_node], x)

            # Set mass flow rates from solver vector
            m_indices_a = self._unknown_indices.get(elem_a_id, [])
            if m_indices_a:
                state_a.m_dot = x[m_indices_a[0]]
            m_indices_b = self._unknown_indices.get(elem_b_id, [])
            if m_indices_b:
                state_b.m_dot = x[m_indices_b[0]]

            # Call htc_and_T() on both elements
            ch_result_a = elem_a.htc_and_T(state_a)
            ch_result_b = elem_b.htc_and_T(state_b)

            # Only proceed if both elements have heat transfer surfaces
            if ch_result_a is not None and ch_result_b is not None:
                # Extract scalars from ChannelResult
                h_a = ch_result_a.h
                T_aw_a = ch_result_a.T_aw
                h_b = ch_result_b.h
                T_aw_b = ch_result_b.T_aw

                # Get convective areas from element surfaces (safe: only elements
                # with a ConvectiveSurface can return non-None ChannelResult)
                A_conv_a = elem_a.surface.area
                A_conv_b = elem_b.surface.area

                A_eff = (
                    wall.contact_area if wall.contact_area is not None else min(A_conv_a, A_conv_b)
                )
                t_over_k = wall.wall_thickness / wall.wall_conductivity

                wall_result = cb.wall_coupling_and_jacobian(
                    h_a, T_aw_a, h_b, T_aw_b, t_over_k, A_eff
                )

                # Convention: Q > 0 means heat flows A->B
                Q_inject = -wall_result.Q if is_side_a else wall_result.Q

                # Store ChannelResults for Jacobian relay
                wall_contributions.append(
                    (
                        wall,
                        Q_inject,
                        wall_result,
                        elem_a_id,
                        elem_b_id,
                        is_side_a,
                        ch_result_a,
                        ch_result_b,
                    )
                )
                processed_walls.add(wall.id)

        # Set accumulated wall Q on the node's wall EnergyBoundary
        if nid in self._wall_ebs:
            self._wall_ebs[nid].Q = sum(c[1] for c in wall_contributions)

        return wall_contributions

    def _propagate_states(self, x: np.ndarray) -> dict:
        """
        Forward-propagate T and Y through the network graph in topological order.
        Also computes global sensitivities (relay) for the chain-rule.
        """
        # Preserve previous iteration's derived states for back-edge wall coupling
        self._prev_derived_states = getattr(self, "_derived_states", {})
        self._derived_states = {}
        # relay[node_id][global_unknown_index] = {'T': dT/dx, 'Y': [dY/dx]}
        relay = {nid: {} for nid in self.network.nodes}

        has_walls = self.network.thermal_coupling_enabled and bool(self.network.walls)

        for nid in self._topological_order:
            node = self.network.nodes[nid]

            if isinstance(node, (PressureBoundary, MassFlowBoundary)):
                T, Y, _ = node.compute_derived_state([])
                self._derived_states[nid] = (T, Y, None)
                continue

            up_elems = self.network.get_upstream_elements(nid)
            up_states = []
            for elem in up_elems:
                state_up = self._get_node_state(self.network.nodes[elem.from_node], x)
                indices = self._unknown_indices.get(elem.id, [])
                if indices:
                    state_up.m_dot = x[indices[0]]
                    # Tag state with element ID for momentum chamber Jacobian
                    state_up._element_id = elem.id
                up_states.append(state_up)

            # Wall Coupling: evaluate walls BEFORE compute_derived_state
            # so that Q flows through the EnergyBoundary abstraction.
            wall_contributions: list[tuple] = []
            if has_walls:
                wall_contributions = self._evaluate_walls_for_node(nid, up_elems, x)

            T, Y, mix_res = node.compute_derived_state(up_states)
            self._derived_states[nid] = (T, Y, mix_res)

            # Store wall coupling debug info on the node
            if wall_contributions:
                m_dot_total_dbg = sum(s.m_dot for s in up_states) or 1.0
                dT_mix_dQ = (mix_res.dT_mix_d_delta_h / m_dot_total_dbg) if mix_res else 0.0
                node._wall_couplings = []
                for (
                    wall,
                    Q_inject,
                    wall_result,
                    ea_id,
                    eb_id,
                    is_side_a,
                    ch_a,
                    ch_b,
                ) in wall_contributions:
                    # dQ_inject/dQ_wall = -1 if is_side_a, +1 otherwise
                    dT_dQ_signed = -dT_mix_dQ if is_side_a else dT_mix_dQ
                    node._wall_couplings.append(
                        {
                            "wall": wall,
                            "wall_result": wall_result,
                            "Q_inject": Q_inject,
                            "is_side_a": is_side_a,
                            "elem_a_id": ea_id,
                            "elem_b_id": eb_id,
                            "dT_mix_dQ": dT_dQ_signed,
                            "ch_result_a": ch_a,
                            "ch_result_b": ch_b,
                        }
                    )

            # Sensitivity Relay (Chain-Rule)
            if mix_res:
                n_species = cb.num_species()
                for i, elem in enumerate(up_elems):
                    from_nid = elem.from_node
                    t_jac = mix_res.dT_mix_d_stream[i]
                    pt_jac = mix_res.dP_total_mix_d_stream[i]
                    # dY_mix_d_stream is indexed by [species][stream]
                    y_jacs = [mix_res.dY_mix_d_stream[k][i] for k in range(n_species)]

                    # 1. Direct dependency on upstream mass flow (if it's a solver unknown)
                    m_indices = self._unknown_indices.get(elem.id)
                    if m_indices:
                        idx = m_indices[0]
                        node_relay = relay[nid].setdefault(
                            idx, {"T": 0.0, "Y": np.zeros(n_species), "P_total": 0.0}
                        )
                        node_relay["T"] += t_jac.d_mdot
                        node_relay["P_total"] += pt_jac.d_mdot
                        for k in range(n_species):
                            node_relay["Y"][k] += y_jacs[k].d_mdot

                    # 2. Recursive dependency on upstream P/T/Y via from_node relay
                    if from_nid in relay:
                        for idx, sens_up in relay[from_nid].items():
                            node_relay = relay[nid].setdefault(
                                idx, {"T": 0.0, "Y": np.zeros(n_species), "P_total": 0.0}
                            )
                            # dT/dx = sum_i( dT/dT_in_i * dT_in_i/dx + dT/dP_in_i * dP_in_i/dx + dT/dY_in_i * dY_in_i/dx )
                            node_relay["T"] += t_jac.d_T * sens_up["T"]
                            node_relay["T"] += t_jac.d_P_total * sens_up["P_total"]
                            node_relay["T"] += np.dot(t_jac.d_Y, sens_up["Y"])

                            # dP_total/dx
                            node_relay["P_total"] += pt_jac.d_T * sens_up["T"]
                            node_relay["P_total"] += pt_jac.d_P_total * sens_up["P_total"]
                            node_relay["P_total"] += np.dot(pt_jac.d_Y, sens_up["Y"])

                            # dY_k/dx
                            for k in range(n_species):
                                node_relay["Y"][k] += y_jacs[k].d_T * sens_up["T"]
                                node_relay["Y"][k] += y_jacs[k].d_P_total * sens_up["P_total"]
                                node_relay["Y"][k] += np.dot(y_jacs[k].d_Y, sens_up["Y"])

                # Wall coupling relay extension: dT_node/dx += dT_mix/dQ * dQ/dx
                # Full chain rule with channel Jacobians:
                #   dQ/dmdot = dQ/dh_a * dh_a/dmdot + dQ/dT_aw_a * dT_aw_a/dmdot + (side B terms)
                #   dQ/dT    = dQ/dh_a * dh_a/dT    + dQ/dT_aw_a * dT_aw_a/dT    + (side B terms)
                if wall_contributions:
                    # dT_mix_d_delta_h = 1/cp_mix is dT/d(specific_h) [K/(J/kg)].
                    # Wall Q is total heat rate [W], so dT/dQ = (1/cp) / m_dot_total.
                    m_dot_total = sum(s.m_dot for s in up_states) or 1.0
                    dT_mix_dQ = mix_res.dT_mix_d_delta_h / m_dot_total
                    for (
                        _wall,
                        _Qi,
                        wall_result,
                        ea_id,
                        eb_id,
                        is_side_a,
                        ch_a,
                        ch_b,
                    ) in wall_contributions:
                        sign = -1.0 if is_side_a else 1.0

                        # Side A contributions
                        from_a = self.network.elements[ea_id].from_node
                        m_indices_a = self._unknown_indices.get(ea_id, [])

                        # Temperature coupling: dT_node/dT_upstream_a
                        if from_a in relay:
                            # dQ/dT_aw_a * dT_aw_a/dT + dQ/dh_a * dh_a/dT
                            factor_T_a = (
                                dT_mix_dQ
                                * sign
                                * (
                                    wall_result.dQ_dT_aw_a * ch_a.dT_aw_dT
                                    + wall_result.dQ_dh_a * ch_a.dh_dT
                                )
                            )
                            for idx, sens_up in relay[from_a].items():
                                nr = relay[nid].setdefault(
                                    idx, {"T": 0.0, "Y": np.zeros(n_species), "P_total": 0.0}
                                )
                                nr["T"] += factor_T_a * sens_up["T"]

                        # Mass-flow coupling: dT_node/dmdot_a (direct)
                        if m_indices_a:
                            idx_mdot_a = m_indices_a[0]
                            # dQ/dh_a * dh_a/dmdot + dQ/dT_aw_a * dT_aw_a/dmdot
                            factor_mdot_a = (
                                dT_mix_dQ
                                * sign
                                * (
                                    wall_result.dQ_dh_a * ch_a.dh_dmdot
                                    + wall_result.dQ_dT_aw_a * ch_a.dT_aw_dmdot
                                )
                            )
                            nr = relay[nid].setdefault(
                                idx_mdot_a, {"T": 0.0, "Y": np.zeros(n_species), "P_total": 0.0}
                            )
                            nr["T"] += factor_mdot_a

                        # Side B contributions
                        from_b = self.network.elements[eb_id].from_node
                        m_indices_b = self._unknown_indices.get(eb_id, [])

                        # Temperature coupling: dT_node/dT_upstream_b
                        if from_b in relay:
                            factor_T_b = (
                                dT_mix_dQ
                                * sign
                                * (
                                    wall_result.dQ_dT_aw_b * ch_b.dT_aw_dT
                                    + wall_result.dQ_dh_b * ch_b.dh_dT
                                )
                            )
                            for idx, sens_up in relay[from_b].items():
                                nr = relay[nid].setdefault(
                                    idx, {"T": 0.0, "Y": np.zeros(n_species), "P_total": 0.0}
                                )
                                nr["T"] += factor_T_b * sens_up["T"]

                        # Mass-flow coupling: dT_node/dmdot_b (direct)
                        if m_indices_b:
                            idx_mdot_b = m_indices_b[0]
                            factor_mdot_b = (
                                dT_mix_dQ
                                * sign
                                * (
                                    wall_result.dQ_dh_b * ch_b.dh_dmdot
                                    + wall_result.dQ_dT_aw_b * ch_b.dT_aw_dmdot
                                )
                            )
                            nr = relay[nid].setdefault(
                                idx_mdot_b, {"T": 0.0, "Y": np.zeros(n_species), "P_total": 0.0}
                            )
                            nr["T"] += factor_mdot_b

            # 3. Add own unknowns (P_total) to the relay
            node_unks = self._unknown_indices.get(nid, [])
            for unk_idx in node_unks:
                unk_name = self.unknown_names[unk_idx]
                if unk_name.endswith(".P_total"):
                    node_relay = relay[nid].setdefault(
                        unk_idx, {"T": 0.0, "Y": np.zeros(cb.num_species()), "P_total": 0.0}
                    )
                    node_relay["P_total"] = 1.0

        return relay

    def _get_node_state_with_prev(self, node: NetworkNode, x: np.ndarray) -> MixtureState:
        """
        Constructs a MixtureState for a node, using previous iteration's derived state
        for back-edges (nodes not yet visited in current topological pass).
        This ensures wall coupling uses lagged values for cross-stream coupling.
        """
        default_Y = self._default_Y

        # Boundaries always use current values
        if isinstance(node, (PressureBoundary, MassFlowBoundary)):
            return self._get_node_state(node, x)

        # For interior nodes, check if already propagated in current iteration
        if node.id in self._derived_states:
            # Use current iteration's values
            return self._get_node_state(node, x)

        # Back-edge: use previous iteration's derived state
        if node.id in self._prev_derived_states:
            T_prev, Y_prev, _ = self._prev_derived_states[node.id]
        else:
            # First iteration: use defaults
            T_prev = 300.0
            Y_prev = default_Y

        # Get pressure unknowns from current x vector
        indices = self._unknown_indices.get(node.id, [])
        unknowns = node.unknowns()

        state_dict = {
            "P": 101325.0,
            "P_total": 101325.0,
            "T": T_prev,
            "T_total": T_prev,
            "m_dot": 0.0,
            "Y": Y_prev,
        }

        for i, unk in zip(indices, unknowns, strict=False):
            var_name = unk.split(".")[-1]
            if var_name in state_dict:
                state_dict[var_name] = x[i]

        return MixtureState(**state_dict)

    def _get_node_state(self, node: NetworkNode, x: np.ndarray) -> MixtureState:
        """Constructs a MixtureState for a given node based on the current solver vector x."""
        default_Y = self._default_Y

        # Determine boundaries
        if isinstance(node, PressureBoundary):
            # For PressureBoundary, defaults.
            # We first check if the test suite provided an initial_guess dict for external states
            guess = getattr(node, "initial_guess", {})

            P_tot = guess.get(f"{node.id}.P_total", getattr(node, "P_total", 101325.0))
            P_stat = guess.get(f"{node.id}.P", P_tot)
            T_tot = guess.get(f"{node.id}.T_total", getattr(node, "T_total", 300.0))
            T_stat = guess.get(f"{node.id}.T", T_tot)
            Y = getattr(node, "Y", None)
            if Y is None:
                Y = default_Y
            return MixtureState(P=P_stat, P_total=P_tot, T=T_stat, T_total=T_tot, m_dot=0.0, Y=Y)

        if isinstance(node, MassFlowBoundary):
            guess = getattr(node, "initial_guess", {})
            m = getattr(node, "m_dot", 0.1)
            T_tot = guess.get(f"{node.id}.T_total", getattr(node, "T_total", 300.0))
            T_stat = guess.get(f"{node.id}.T", T_tot)
            Y = getattr(node, "Y", None)
            if Y is None:
                Y = default_Y

            p_val = guess.get(f"{node.id}.P", 101325.0)
            ptot_val = guess.get(f"{node.id}.P_total", 101325.0)

            indices = self._unknown_indices.get(node.id, [])
            unknowns = node.unknowns()
            for i, unk in zip(indices, unknowns, strict=False):
                if unk.endswith(".P"):
                    p_val = x[i]
                elif unk.endswith(".P_total"):
                    ptot_val = x[i]

            return MixtureState(P=p_val, P_total=ptot_val, T=T_stat, T_total=T_tot, m_dot=m, Y=Y)

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
            "Y": default_Y,
        }

        # Override with derived state if available
        if node.id in self._derived_states:
            T, Y, _ = self._derived_states[node.id]
            state_dict["T"] = T
            state_dict["T_total"] = T
            state_dict["Y"] = Y

        for i, unk in zip(indices, unknowns, strict=False):
            var_name = unk.split(".")[-1]
            if var_name in state_dict:
                state_dict[var_name] = x[i]
            elif var_name.startswith("Y[") and var_name.endswith("]"):
                idx_str = var_name[2:-1]
                idx = int(idx_str)
                # default_Y is shared, so we must copy it if we intend to mutate it.
                if state_dict["Y"] is default_Y:
                    state_dict["Y"] = list(default_Y)
                state_dict["Y"][idx] = x[i]

        # P1-4: Clamp fractions to [0, 1] before creating state to avoid NaNs in thermo conversion
        Y_raw = state_dict["Y"]
        Y_safe = [max(1e-12, min(1.0, float(yi))) for yi in Y_raw]
        y_sum = sum(Y_safe)
        if y_sum > 0:
            Y_safe = [yi / y_sum for yi in Y_safe]
        state_dict["Y"] = Y_safe

        return MixtureState(**state_dict)

    def _residuals_and_jacobian(self, x: np.ndarray) -> tuple[np.ndarray, Any]:
        """
        Evaluates the global residual vector and its sparse Jacobian.
        """
        import scipy.sparse as sp

        # Pure Pressure-Flow: Propagate states forward first
        relay = self._propagate_states(x)

        res = []
        rows = []
        cols = []
        data = []

        # 1. Node Residuals
        for node_id, node in self.network.nodes.items():
            if isinstance(node, PressureBoundary):
                continue

            start_res_idx = len(res)

            # Node thermodynamic/stagnation constraints
            state = self._get_node_state(node, x)
            node_res, node_jac = node.residuals(state)
            res.extend(node_res)

            # Assemble node Jacobian entries
            for eq_idx, var_derivs in node_jac.items():
                row = start_res_idx + eq_idx
                for unk_name, deriv in var_derivs.items():
                    if unk_name in self._name_to_index:
                        rows.append(row)
                        cols.append(self._name_to_index[unk_name])
                        data.append(deriv)
                    elif "." in unk_name:
                        # Relay for derived properties (T, Y) in node residuals
                        parts = unk_name.split(".")
                        nid, prop = parts[0], parts[1]
                        if nid in relay:
                            for unk_idx, sens_pkg in relay[nid].items():
                                if prop == "T" and "T" in sens_pkg:
                                    rows.append(row)
                                    cols.append(unk_idx)
                                    data.append(deriv * sens_pkg["T"])
                                elif prop.startswith("Y") and "Y" in sens_pkg:
                                    # Handle species composition relay
                                    if "[" in prop:
                                        s_idx = int(prop.split("[")[1].replace("]", ""))
                                        rows.append(row)
                                        cols.append(unk_idx)
                                        data.append(deriv * sens_pkg["Y"][s_idx])

            # Mass Conservation: Sum(m_dot_in) - Sum(m_dot_out) = 0
            mass_res_idx = len(res)
            upstream_elems = self.network.get_upstream_elements(node_id)
            downstream_elems = self.network.get_downstream_elements(node_id)

            m_dot_in = 0.0
            m_dot_out = 0.0

            if isinstance(node, MassFlowBoundary):
                m_dot_in += getattr(node, "m_dot", 0.0)

            for elem in upstream_elems:
                indices = self._unknown_indices.get(elem.id)
                if indices:
                    m_dot_in += x[indices[0]]
                    rows.append(mass_res_idx)
                    cols.append(indices[0])
                    data.append(1.0)
                else:
                    from_node = self.network.nodes[elem.from_node]
                    m_dot_in += getattr(from_node, "m_dot", 0.0)

            for elem in downstream_elems:
                indices = self._unknown_indices.get(elem.id)
                if indices:
                    m_dot_out += x[indices[0]]
                    rows.append(mass_res_idx)
                    cols.append(indices[0])
                    data.append(-1.0)
                else:
                    to_node = self.network.nodes[elem.to_node]
                    m_dot_out += getattr(to_node, "m_dot", 0.0)

            res.append(m_dot_in - m_dot_out)

        # 2. Element Residuals
        for elem_id, element in self.network.elements.items():
            start_res_idx = len(res)

            node_in = self.network.nodes[element.from_node]
            node_out = self.network.nodes[element.to_node]

            state_in = self._get_node_state(node_in, x)
            state_out = self._get_node_state(node_out, x)

            # Unpack element-specific unknowns (like m_dot) into the state
            m_indices = self._unknown_indices.get(elem_id)
            if m_indices:
                state_in.m_dot = x[m_indices[0]]
                state_out.m_dot = x[m_indices[0]]

            elem_res, elem_jac = element.residuals(state_in, state_out)
            res.extend(elem_res)

            for eq_idx, var_derivs in elem_jac.items():
                row = start_res_idx + eq_idx
                for unk_name, deriv in var_derivs.items():
                    if unk_name in self._name_to_index:
                        rows.append(row)
                        cols.append(self._name_to_index[unk_name])
                        data.append(deriv)
                    elif "." in unk_name:
                        # Relay for derived properties (T, Y) in element residuals
                        parts = unk_name.split(".")
                        nid, prop = parts[0], parts[1]
                        if nid in relay:
                            for unk_idx, sens_pkg in relay[nid].items():
                                if prop == "T" and "T" in sens_pkg:
                                    rows.append(row)
                                    cols.append(unk_idx)
                                    data.append(deriv * sens_pkg["T"])
                                elif prop.startswith("Y") and "Y" in sens_pkg:
                                    if "[" in prop:
                                        s_idx = int(prop.split("[")[1].replace("]", ""))
                                        rows.append(row)
                                        cols.append(unk_idx)
                                        data.append(deriv * sens_pkg["Y"][s_idx])

        jac_sparse = sp.coo_matrix((data, (rows, cols)), shape=(len(res), len(x))).tocsr()

        return np.array(res, dtype=float), jac_sparse

    def _residuals(self, x: np.ndarray) -> np.ndarray:
        """
        Evaluates the global residual vector given the current unknown vector x.
        """
        res, _ = self._residuals_and_jacobian(x)
        return res

    def solve(
        self,
        method: str = "hybr",
        timeout: float | None = None,
        options: dict[str, Any] | None = None,
        use_jac: bool = True,
    ) -> dict[str, float]:
        """
        Executes the root finding algorithm to solve the network.
        Returns a dictionary mapping unknown names to their solved values.

        Args:
            method: The scipy.optimize.root method to use.
            timeout: Maximum wall-clock time in seconds before stopping and returning the best iterate.
            options: Dictionary of solver options (e.g., maxiter, maxfev, xtol).
        """
        if method not in self.SUPPORTED_METHODS:
            raise ValueError(
                f"Method '{method}' is not supported. Supported methods are: {', '.join(self.SUPPORTED_METHODS)}"
            )

        if options is None:
            options = {}

        # Default iteration limits to prevent infinite loops in C extensions
        # Note: different methods in root use different names for iteration limits
        if method in ("hybr", "lm"):
            if "maxfev" not in options:
                options["maxfev"] = 500
        else:
            if "maxiter" not in options:
                options["maxiter"] = 500

        # Ensure topology is resolved before solving
        self.network.resolve_all_topology()
        self.network.validate()

        x0 = self._build_x0()

        # Dimension check
        res0 = self._residuals(x0)
        if len(x0) != len(res0):
            raise ValueError(f"System not square: {len(x0)} unknowns vs {len(res0)} equations.")

        # State for timeout and best iterate tracking
        start_time = time.perf_counter()
        best_x = x0.copy()
        best_res_norm = np.linalg.norm(res0)

        def residuals_wrapper(x: np.ndarray) -> Any:
            nonlocal best_x, best_res_norm

            # Check timeout
            if timeout is not None and (time.perf_counter() - start_time) > timeout:
                raise SolverTimeoutError(f"Solver timed out after {timeout} seconds.")

            try:
                # Evaluate residuals and jacobian
                res, jac = self._residuals_and_jacobian(x)

                # Track best iterate
                res_norm = np.linalg.norm(res)
                if res_norm < best_res_norm:
                    best_res_norm = res_norm
                    best_x = x.copy()

                if use_jac and method in ("hybr", "lm"):
                    # Scipy's root (hybr, lm) does not support sparse Jacobians directly.
                    # Use toarray() to provide the expected dense format while still
                    # leveraging analytical calculation.
                    return res, jac.toarray()
                else:
                    # Other methods might not support (f, J) or sparse J properly via the jac=True flag
                    return res
            except Exception:
                # P1-5: Avoid printing in hot path unless absolutely necessary
                raise

        # Solve
        try:
            sol = root(residuals_wrapper, x0, method=method, options=options, jac=use_jac)
            final_x = sol.x
            success = sol.success
            message = sol.message
            final_norm = float(np.linalg.norm(sol.fun))
        except SolverTimeoutError as e:
            final_x = best_x
            success = False
            message = str(e)
            final_norm = float(best_res_norm)
        except Exception as e:
            # Fallback for unexpected errors during residuals evaluation
            final_x = best_x
            success = False
            message = f"Unexpected error during residual evaluation: {e}"
            final_norm = float(best_res_norm)

        if not success:
            pass

        if not success:
            warnings.warn(
                f"NetworkSolver did not converge: {message} "
                f"(|F|={final_norm:.3e}). "
                "Returning best iterate.",
                stacklevel=2,
            )

        sol_dict = dict(zip(self.unknown_names, final_x, strict=False))

        # Include derived states (T, Y) for consistency and state inspection
        self._propagate_states(final_x)
        for nid, (T, Y, _) in self._derived_states.items():
            sol_dict[f"{nid}.T"] = T
            sol_dict[f"{nid}.T_total"] = T  # Canonical total T
            for i, yi in enumerate(Y):
                sol_dict[f"{nid}.Y[{i}]"] = float(yi)

        sol_dict["__success__"] = success
        sol_dict["__message__"] = message
        sol_dict["__final_norm__"] = final_norm
        return sol_dict
