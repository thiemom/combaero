import math
import time
import warnings
from typing import Any, Literal

import numpy as np
from scipy.optimize import root

import combaero as cb

from .components import (
    EnergyBoundary,
    LosslessConnectionElement,
    MassFlowBoundary,
    MixtureState,
    NetworkNode,
    PipeElement,
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

    @staticmethod
    def _is_pressure_unknown(name: str) -> bool:
        return name.endswith(".P") or name.endswith(".P_total")

    @staticmethod
    def _is_mdot_unknown(name: str) -> bool:
        return name.endswith(".m_dot")

    @staticmethod
    def _to_incompressible_regime(regime: str) -> str | None:
        if regime in ("compressible", "compressible_fanno"):
            return "incompressible"
        return None

    def _compressible_element_overrides(self) -> dict[str, str]:
        """Map compressible element IDs to their incompressible fallback regimes."""
        overrides: dict[str, str] = {}
        for elem_id, element in self.network.elements.items():
            regime = getattr(element, "regime", None)
            if not isinstance(regime, str):
                continue
            fallback = self._to_incompressible_regime(regime)
            if fallback is not None:
                overrides[elem_id] = fallback
        return overrides

    def _reference_scales(self, x0_use: np.ndarray) -> tuple[float, float]:
        pressure_vals = [
            abs(float(v))
            for n, v in zip(self.unknown_names, x0_use, strict=False)
            if self._is_pressure_unknown(n)
        ]
        mdot_vals = [
            abs(float(v))
            for n, v in zip(self.unknown_names, x0_use, strict=False)
            if self._is_mdot_unknown(n)
        ]

        ref_p_abs = float(np.median(pressure_vals)) if pressure_vals else 1.0e5
        ref_mdot = float(np.median(mdot_vals)) if mdot_vals else 1.0

        # Pressure-residual scale: use boundary DP (not absolute P)
        # to match actual residual magnitudes in pressure equations
        if len(pressure_vals) >= 2:
            p_spread = max(pressure_vals) - min(pressure_vals)
            ref_p_residual = max(p_spread, ref_p_abs * 1e-4)
        else:
            ref_p_residual = ref_p_abs * 0.01

        ref_p_residual = max(ref_p_residual, 1.0)
        ref_mdot = max(ref_mdot, 1e-6)
        return ref_p_residual, ref_mdot

    def _build_residual_scales(self, x0_use: np.ndarray) -> np.ndarray:
        """Build characteristic residual scales by equation type.

        Pressure-like residuals use ``ref_p_residual`` [Pa] (based on boundary DP).
        Mass-flow-like residuals use ``ref_mdot`` [kg/s].
        """
        ref_p, ref_mdot = self._reference_scales(x0_use)

        scales: list[float] = []

        for _node_id, node in self.network.nodes.items():
            if isinstance(node, PressureBoundary):
                continue

            state = self._get_node_state(node, x0_use)
            node_res, _ = node.residuals(state)
            scales.extend([ref_p] * len(node_res))

            # Node mass conservation equation appended in _residuals_and_jacobian
            scales.append(ref_mdot)

        for element in self.network.elements.values():
            elem_scale = (
                ref_p if isinstance(element, (PipeElement, LosslessConnectionElement)) else ref_mdot
            )
            scales.extend([elem_scale] * element.n_equations())

        return np.asarray(scales, dtype=float)

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

    def _infer_reference_state(self) -> dict[str, Any]:
        """Derive sensible default values for unknowns from boundary nodes.

        Uses the **minimum** ``PressureBoundary`` pressure so that
        interior nodes start near the downstream end of the network,
        consistent with the previous hard-coded 1 atm default for
        low-pressure networks while giving a physically meaningful
        starting point for high-pressure ones.  Temperature and
        composition are averaged
        across all boundary nodes that carry them.  Mass-flow defaults
        to the total prescribed mass flow divided by the number of flow
        elements.

        Returns a dict with keys ``'P'``, ``'T'``, ``'Y'``, ``'m_dot'``.
        """
        pressures: list[float] = []
        temperatures: list[float] = []
        compositions: list[list[float]] = []
        total_mdot = 0.0

        for node in self.network.nodes.values():
            if isinstance(node, PressureBoundary):
                pressures.append(getattr(node, "P_total", 101325.0))
                T = getattr(node, "T_total", None)
                if T is not None:
                    temperatures.append(T)
                Y = getattr(node, "Y", None)
                if Y is not None:
                    compositions.append(list(Y))

            elif isinstance(node, MassFlowBoundary):
                T = getattr(node, "T_total", None)
                if T is not None:
                    temperatures.append(T)
                Y = getattr(node, "Y", None)
                if Y is not None:
                    compositions.append(list(Y))
                total_mdot += getattr(node, "m_dot", 0.0)

        ref_P = float(min(pressures)) if pressures else 101325.0
        ref_T = float(np.mean(temperatures)) if temperatures else 300.0
        if compositions:
            ref_Y = [float(v) for v in np.mean(compositions, axis=0)]
        else:
            ref_Y = list(self._default_Y)

        n_elems = max(len(self.network.elements), 1)
        ref_mdot = total_mdot / n_elems if total_mdot > 0 else 0.1

        return {"P": ref_P, "T": ref_T, "Y": ref_Y, "m_dot": ref_mdot}

    def _propagate_pressure_guess(self, ref: dict[str, Any]) -> dict[str, float]:
        """Estimate per-node pressures by walking the graph from boundaries.

        Starting from every ``PressureBoundary``, pressures are propagated
        both downstream (subtracting dP) and upstream (adding dP) per
        element via BFS.  ``MassFlowBoundary`` nodes are *not* seeded -
        their pressure is discovered by reverse propagation from known
        outlets, ensuring the pressure gradient points in the correct
        flow direction.
        """
        p_guess: dict[str, float] = {}

        # Seed only PressureBoundary nodes (known pressures).
        # MassFlowBoundary nodes have unknown pressure - let BFS
        # discover them by propagating upstream from known outlets.
        for nid, node in self.network.nodes.items():
            if isinstance(node, PressureBoundary):
                p_guess[nid] = getattr(node, "P_total", ref["P"])

        # Estimated per-element dP: spread across elements, or 0.1 % of ref_P
        seed_pressures = list(p_guess.values())
        if len(seed_pressures) >= 2:
            dp_est = (max(seed_pressures) - min(seed_pressures)) / max(
                len(self.network.elements), 1
            )
        else:
            dp_est = ref["P"] * 0.001

        # BFS propagation from known nodes
        visited = set(p_guess.keys())
        queue = list(visited)
        while queue:
            nid = queue.pop(0)
            for elem in self.network.get_downstream_elements(nid):
                to_id = elem.to_node
                if to_id not in visited:
                    p_guess[to_id] = p_guess[nid] - dp_est
                    visited.add(to_id)
                    queue.append(to_id)
            for elem in self.network.get_upstream_elements(nid):
                from_id = elem.from_node
                if from_id not in visited:
                    p_guess[from_id] = p_guess[nid] + dp_est
                    visited.add(from_id)
                    queue.append(from_id)

        return p_guess

    def _propagate_mdot_guess(self, ref: dict[str, Any]) -> dict[str, float]:
        """Estimate per-element mass flow by propagating boundary flows.

        Mass flow is seeded from ``MassFlowBoundary`` nodes and propagated
        in topological order. Merge nodes sum incoming flows; split nodes
        divide the available flow equally among downstream elements.

        For compressible elements, the propagated value is capped using
        an isentropic estimate at a conservative target Mach number
        (``M=0.3``), based on the element area and reference state.
        This keeps startup states safely subsonic while preserving
        topology-derived split behavior.

        Returns a dict mapping element ID to estimated ``m_dot`` [kg/s].
        """
        mach_target = 0.3
        mdot_guess: dict[str, float] = {}
        node_flow: dict[str, float] = {
            nid: float(getattr(node, "m_dot", 0.0))
            for nid, node in self.network.nodes.items()
            if isinstance(node, MassFlowBoundary)
        }

        ref_X = list(cb.mass_to_mole(ref["Y"]))
        gamma_ref = float(cb.isentropic_expansion_coefficient(ref["T"], ref_X))
        rho_ref = float(cb.density(ref["T"], ref["P"], ref_X))
        a_ref = float(cb.speed_of_sound(ref["T"], ref_X))

        term = 1.0 + 0.5 * (gamma_ref - 1.0) * mach_target**2
        exponent = -(gamma_ref + 1.0) / (2.0 * (gamma_ref - 1.0))
        mdot_ratio = mach_target * term**exponent

        for nid in self._compute_topological_order():
            node = self.network.nodes[nid]

            if isinstance(node, MassFlowBoundary):
                available = node_flow[nid]
            else:
                available = node_flow.get(nid, ref["m_dot"])

            down_elems = self.network.get_downstream_elements(nid)
            if not down_elems:
                continue

            share = available / len(down_elems)
            if abs(share) < 1e-12:
                share = ref["m_dot"]

            for elem in down_elems:
                elem_share = share
                regime = getattr(elem, "regime", None)
                if regime in ("compressible", "compressible_fanno"):
                    area = getattr(elem, "area", None)
                    if area is not None and area > 0.0:
                        mdot_cap = rho_ref * a_ref * area * mdot_ratio
                        elem_share = min(elem_share, mdot_cap)

                mdot_guess[elem.id] = elem_share
                node_flow[elem.to_node] = node_flow.get(elem.to_node, 0.0) + elem_share

        return mdot_guess

    def _build_x0(self) -> np.ndarray:
        """
        Constructs the initial guess vector by gathering all unknowns.

        Pressure unknowns are initialised via topological propagation
        from boundary nodes with an estimated per-element drop.
        Temperature and composition are averaged across boundaries.
        Mass-flow is propagated from ``MassFlowBoundary`` nodes through
        the topology, with conservative capping for compressible elements.

        Returns a flat numpy array.
        """
        self.unknown_names = []
        self._unknown_indices = {}
        x0_list = []

        # --- Infer sensible defaults from boundary nodes ----------------
        ref = self._infer_reference_state()
        p_guess = self._propagate_pressure_guess(ref)
        mdot_guess = self._propagate_mdot_guess(ref)

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
                            x0_list.append(p_guess.get(node_id, ref["P"]))
                        elif unk.endswith(".T") or unk.endswith(".T_total"):
                            x0_list.append(ref["T"])
                        elif ".Y[" in unk:
                            idx = int(unk.split(".Y[")[1].replace("]", ""))
                            x0_list.append(ref["Y"][idx])
                        else:
                            x0_list.append(ref["m_dot"])
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
                        x0_list.append(mdot_guess.get(elem_id, ref["m_dot"]))
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
        x0: np.ndarray | None = None,
        init_strategy: Literal["default", "incompressible_warmstart"] = "default",
        warmstart_maxfev: int = 150,
    ) -> dict[str, float]:
        """
        Executes the root finding algorithm to solve the network.
        Returns a dictionary mapping unknown names to their solved values.

        Args:
            method: The scipy.optimize.root method to use.
            timeout: Maximum wall-clock time in seconds before stopping
                and returning the best iterate.
            options: Dictionary of solver options (e.g., maxiter, maxfev,
                xtol).
            use_jac: Whether to supply the analytical Jacobian to the
                solver (default True).
            x0: Optional initial guess vector.  When ``None`` (default),
                an initial guess is built automatically.  Pass a
                previous solution's ``_last_x`` for warm-starting
                parameter sweeps.
            init_strategy: Initialization strategy. ``default`` uses
                direct x0 construction; ``incompressible_warmstart``
                first solves an incompressible-regime proxy network and
                uses that solution as x0.
            warmstart_maxfev: Maximum function evaluations for the
                incompressible proxy solve when
                ``init_strategy='incompressible_warmstart'``.
        """
        if method not in self.SUPPORTED_METHODS:
            raise ValueError(
                f"Method '{method}' is not supported. Supported methods are: {', '.join(self.SUPPORTED_METHODS)}"
            )

        if init_strategy not in ("default", "incompressible_warmstart"):
            raise ValueError("init_strategy must be one of: 'default', 'incompressible_warmstart'.")
        if warmstart_maxfev <= 0:
            raise ValueError("warmstart_maxfev must be > 0.")

        options = {} if options is None else dict(options)

        # Normalise caller option names across methods.
        # scipy.root uses different names for the iteration limit:
        # - hybr: maxfev
        # - lm and quasi-Newton/fixed-point methods: maxiter
        # Translate before applying defaults to avoid stale duplicates.
        if method != "hybr" and "maxfev" in options:
            options.setdefault("maxiter", options.pop("maxfev"))

        if method == "hybr":
            options.setdefault("maxfev", 500)
        else:
            options.setdefault("maxiter", 500)

        # Ensure topology is resolved before solving
        self.network.resolve_all_topology()
        self.network.validate()

        warmstart_x0: np.ndarray | None = None
        if x0 is None and init_strategy == "incompressible_warmstart":
            overrides = self._compressible_element_overrides()
            if overrides:
                original_regimes = {
                    elem_id: getattr(self.network.elements[elem_id], "regime", None)
                    for elem_id in overrides
                }
                try:
                    for elem_id, fallback_regime in overrides.items():
                        self.network.elements[elem_id].regime = fallback_regime

                    self.network.resolve_all_topology()

                    warm_options = {"maxfev": int(warmstart_maxfev)}
                    warm_sol = self.solve(
                        method="hybr",
                        timeout=timeout,
                        options=warm_options,
                        use_jac=use_jac,
                        x0=None,
                        init_strategy="default",
                    )

                    if bool(warm_sol.get("__success__", False)):
                        candidate = getattr(self, "_last_x", None)
                        if candidate is not None:
                            warmstart_x0 = np.asarray(candidate, dtype=float)
                except Exception as e:
                    warnings.warn(
                        "Incompressible warmstart failed. Proceeding with default initialization: "
                        f"{e}",
                        RuntimeWarning,
                        stacklevel=2,
                    )
                finally:
                    for elem_id, regime in original_regimes.items():
                        self.network.elements[elem_id].regime = regime
                    self.network.resolve_all_topology()

        # --- Initial guess: user-supplied (warm-start) or auto-built ----
        x0_auto = self._build_x0()
        if x0 is not None:
            if len(x0) != len(x0_auto):
                raise ValueError(
                    f"Warm-start x0 length ({len(x0)}) does not match "
                    f"number of unknowns ({len(x0_auto)})."
                )
            x0_use = np.asarray(x0, dtype=float)
        elif warmstart_x0 is not None:
            x0_use = warmstart_x0 if len(warmstart_x0) == len(x0_auto) else x0_auto
        else:
            x0_use = x0_auto

        # Seed forward-propagated thermo/composition states once before the
        # first residual/scaling pass. This lets the initial evaluation use
        # derived T/Y (e.g., combustion-updated temperatures) rather than
        # purely default placeholders.
        try:
            self._propagate_states(x0_use)
        except Exception as e:
            # pragma: no cover - defensive initialization fallback
            warnings.warn(
                "Initial state propagation failed. Proceeding with direct residual "
                f"evaluation, but convergence may be affected: {e}",
                category=RuntimeWarning,
                stacklevel=2,
            )

        # Dimension check
        res0 = self._residuals(x0_use)
        if len(x0_use) != len(res0):
            raise ValueError(f"System not square: {len(x0_use)} unknowns vs {len(res0)} equations.")

        # --- Variable & residual scaling --------------------------------
        # Scale unknowns so they are O(1): x_real = x_scaled * D_x
        # Scale residuals so they are O(1): F_scaled = F_real * D_f
        # This dramatically improves Jacobian conditioning when unknowns
        # span different orders of magnitude (Pa vs kg/s).
        D_x = np.abs(x0_use).copy()
        D_x[D_x < 1e-12] = 1.0  # avoid division by zero for near-zero guesses
        D_f = self._build_residual_scales(x0_use)
        inv_D_x = 1.0 / D_x
        inv_D_f = 1.0 / D_f

        x0_scaled = x0_use * inv_D_x  # should be ~1.0

        # State for timeout and best iterate tracking (in real space)
        start_time = time.perf_counter()
        best_x = x0_use.copy()
        best_res_norm = np.linalg.norm(res0)

        def residuals_wrapper(x_scaled: np.ndarray) -> Any:
            nonlocal best_x, best_res_norm

            # Check timeout
            if timeout is not None and (time.perf_counter() - start_time) > timeout:
                raise SolverTimeoutError(f"Solver timed out after {timeout} seconds.")

            # Un-scale to real space
            x_real = x_scaled * D_x

            try:
                # Evaluate residuals and jacobian
                res, jac = self._residuals_and_jacobian(x_real)

                # Track best iterate (in real space, unscaled residual)
                res_norm = np.linalg.norm(res)
                if res_norm < best_res_norm:
                    best_res_norm = res_norm
                    best_x = x_real.copy()

                # Scale residuals
                res_scaled = res * inv_D_f

                if use_jac and method in ("hybr", "lm"):
                    # J_scaled = diag(inv_D_f) @ J @ diag(D_x)
                    jac_dense = jac.toarray()
                    jac_scaled = (inv_D_f[:, None] * jac_dense) * D_x[None, :]
                    return res_scaled, jac_scaled
                else:
                    return res_scaled
            except SolverTimeoutError:
                raise
            except Exception:
                # Return a large penalty residual so the trust-region
                # solver can shrink its step instead of terminating.
                # Warn once so developers can trace unexpected failures.
                if not getattr(self, "_penalty_warning_issued", False):
                    warnings.warn(
                        "NetworkSolver hit an unphysical intermediate state; "
                        "returning penalty residual to allow step reduction. "
                        "Enable debug logging to see the original exception.",
                        RuntimeWarning,
                        stacklevel=2,
                    )
                    self._penalty_warning_issued = True
                penalty = np.full(len(x_scaled), 10.0, dtype=float)
                if use_jac and method in ("hybr", "lm"):
                    return penalty, np.eye(len(x_scaled))
                return penalty

        # Solve with timing
        _RESIDUAL_TOL = 1e-3
        solve_start_time = time.time()
        try:
            sol = root(residuals_wrapper, x0_scaled, method=method, options=options, jac=use_jac)
            solve_end_time = time.time()
            self._wall_time = solve_end_time - solve_start_time
            final_x = sol.x * D_x  # un-scale
            success = sol.success
            message = sol.message
            # Compute unscaled residual norm for reporting
            final_res = self._residuals(final_x)
            final_norm = float(np.linalg.norm(final_res))
            # Guard against methods (e.g. lm) that report success at a
            # non-zero local minimum of ||F||^2.
            if success and final_norm > _RESIDUAL_TOL:
                success = False
                message = (
                    f"Solver reported success but |F|={final_norm:.3e} "
                    f"exceeds residual tolerance ({_RESIDUAL_TOL})."
                )
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

        # Automatic fallback: if the chosen method failed and isn't hybr,
        # retry with hybr (Powell's method) which is more robust for square
        # nonlinear systems.
        if not success and method != "hybr":
            try:
                sol2 = root(residuals_wrapper, x0_scaled, method="hybr", jac=use_jac)
                fallback_x = sol2.x * D_x
                fallback_res = self._residuals(fallback_x)
                fallback_norm = float(np.linalg.norm(fallback_res))
                if sol2.success and fallback_norm < _RESIDUAL_TOL:
                    final_x = fallback_x
                    success = True
                    message = f"Converged after fallback to hybr (original {method} failed)."
                    final_norm = fallback_norm
            except Exception:
                pass  # keep original failure

        if not success:
            warnings.warn(
                f"NetworkSolver did not converge: {message} "
                f"(|F|={final_norm:.3e}). "
                "Returning best iterate.",
                stacklevel=2,
            )

        # Store solution for warm-start reuse
        self._last_x = final_x.copy()

        # Store diagnostic data
        self._diagnostic_data = {
            "x0": x0_use.copy(),
            "solution": final_x.copy(),
            "converged": success,
            "final_norm": final_norm,
            "message": message,
            "unknown_names": self.unknown_names.copy(),
            "solver_method": method,
            "solver_options": options.copy() if options else {},
            "wall_time": getattr(self, "_wall_time", None),
        }

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

    def extract_complete_states(self, result: dict[str, float]) -> dict[str, Any]:
        """
        Extract CompleteState for all network nodes regardless of convergence status.

        This method provides comprehensive thermodynamic and transport properties
        for every node in the network, leveraging the existing C++ CompleteState
        structure for maximum accuracy and performance.

        Args:
            result: Dictionary returned by NetworkSolver.solve() containing
                node temperatures, pressures, and compositions.

        Returns:
            Dictionary mapping node IDs to CompleteState objects containing
            28 properties each (16 thermodynamic + 12 transport).

        Notes:
            - Works regardless of solver convergence status
            - Uses C++ CompleteState for fast, accurate property calculations
            - All properties are computed from the same (T, P, X) state
            - Returns nan for nodes with missing/invalid data
        """
        # CompleteState is accessed via complete_state function

        complete_states = {}

        # Get composition from result or use default air
        # Try to get composition from first node that has Y data
        default_X = None
        for node_id in self.network.nodes:
            if f"{node_id}.Y[0]" in result:
                # Extract composition array
                Y = []
                i = 0
                while f"{node_id}.Y[{i}]" in result:
                    Y.append(result[f"{node_id}.Y[{i}]"])
                    i += 1
                if Y:
                    # Convert mass fractions to mole fractions
                    from combaero import mass_to_mole

                    default_X = mass_to_mole(np.array(Y))
                    break

        # Fallback to air if no composition found
        if default_X is None:
            default_X = cb.humid_air_composition(300.0, 101325.0, 0.0)  # Dry air at STP

        # Extract CompleteState for each node
        for node_id in self.network.nodes:
            try:
                # Get temperature and pressure from result
                T = result.get(f"{node_id}.T", float("nan"))
                P = result.get(f"{node_id}.P", float("nan"))

                # Get composition for this node (fallback to default)
                X = default_X.copy()
                if f"{node_id}.Y[0]" in result:
                    Y = []
                    i = 0
                    while f"{node_id}.Y[{i}]" in result:
                        Y.append(result[f"{node_id}.Y[{i}]"])
                        i += 1
                    if Y:
                        X = mass_to_mole(np.array(Y))

                # Create CompleteState if we have valid T and P
                if not (math.isnan(T) or math.isnan(P)):
                    complete_state_obj = cb.complete_state(T, P, X)
                    complete_states[node_id] = complete_state_obj
                else:
                    # Create a placeholder with nan values
                    complete_states[node_id] = None

            except Exception:
                # Store None for nodes that fail state calculation
                complete_states[node_id] = None
                # Could add logging here for debugging

        return complete_states

    def get_diagnostic_report(self) -> dict[str, Any]:
        """
        Generate comprehensive diagnostic report comparing initial guess to final solution.

        Returns:
            Dictionary containing detailed analysis of solver behavior including
            x0 vs solution comparison, convergence metrics, and per-variable analysis.
        """
        if not hasattr(self, "_diagnostic_data"):
            raise RuntimeError("No diagnostic data available. Run solve() first.")

        data = self._diagnostic_data
        x0 = data["x0"]
        solution = data["solution"]
        unknown_names = data["unknown_names"]

        # Calculate differences and statistics
        differences = solution - x0
        abs_differences = np.abs(differences)
        rel_differences = np.abs(differences / (x0 + 1e-10))  # Avoid division by zero

        # Per-variable analysis
        variable_analysis = {}
        for i, name in enumerate(unknown_names):
            variable_analysis[name] = {
                "x0": float(x0[i]),
                "solution": float(solution[i]),
                "difference": float(differences[i]),
                "abs_difference": float(abs_differences[i]),
                "rel_difference": float(rel_differences[i]),
            }

        # Overall statistics
        overall_stats = {
            "max_abs_difference": float(np.max(abs_differences)),
            "mean_abs_difference": float(np.mean(abs_differences)),
            "median_abs_difference": float(np.median(abs_differences)),
            "max_rel_difference": float(np.max(rel_differences)),
            "mean_rel_difference": float(np.mean(rel_differences)),
            "median_rel_difference": float(np.median(rel_differences)),
            "std_difference": float(np.std(differences)),
        }

        # Convergence analysis
        convergence_analysis = {
            "converged": data["converged"],
            "final_norm": data["final_norm"],
            "message": data["message"],
            "wall_time": data["wall_time"],
            "solver_method": data["solver_method"],
            "solver_options": data["solver_options"],
            "num_unknowns": len(unknown_names),
        }

        # Identify most changed variables
        most_changed_indices = np.argsort(abs_differences)[-5:]  # Top 5
        most_changed = [
            {
                "name": unknown_names[i],
                "abs_difference": float(abs_differences[i]),
                "rel_difference": float(rel_differences[i]),
                "x0": float(x0[i]),
                "solution": float(solution[i]),
            }
            for i in reversed(most_changed_indices)
        ]

        return {
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "variable_analysis": variable_analysis,
            "overall_statistics": overall_stats,
            "convergence_analysis": convergence_analysis,
            "most_changed_variables": most_changed,
            "unknown_names": unknown_names,
        }
