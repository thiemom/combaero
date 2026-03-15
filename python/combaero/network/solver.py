import time
import warnings
from typing import Any

import numpy as np
from scipy.optimize import root

import combaero as cb

from .components import (
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

    CRITICAL TODOs for Phase 2 Completion (Combustion & Mixing):
    ================================================================

    1. IMPLEMENT PROPER COMBUSTION RESIDUALS in CombustorNode:
       - Replace placeholder with cb.solver.enthalpy_and_jacobian() calls
       - Add species composition residuals using atom balancing
       - Implement pressure drop residuals with analytical Jacobians
       - Handle combustion efficiency properly

    2. FIX STOICHIOMETRY IMPLEMENTATION:
       - stoichiometric_products() needs access to molecular_structures
       - Must properly handle rich/lean mixtures (φ ≠ 1)
       - Implement water-gas shift for rich combustion

    4. COMPLETE CHAMBER NODES:
       - PlenumNode: Move mixing residuals from solver to node
       - MomentumChamberNode: Add momentum conservation
       - Both: Use cb.solver.* functions for analytical Jacobians

    5. ADD COMPREHENSIVE TESTS:
       - Jacobian accuracy tests (finite difference vs analytical)
       - Cantera validation for combustion cases
       - Energy conservation verification
       - Species tracking tests

    6. RESOLVE TECHNICAL DEBT:
       - mix_streams() temperature calculation issues
       - Circular import issues in combustion module
       - Proper error handling for edge cases

    See COMBUSTION_ELEMENTS.md and SOLVER.md for detailed requirements.
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

    def _propagate_states(self, x: np.ndarray) -> dict:
        """
        Forward-propagate T and Y through the network graph in topological order.
        Also computes global sensitivities (relay) for the chain-rule.
        """
        self._derived_states = {}
        # relay[node_id][global_unknown_index] = {'T': dT/dx, 'Y': [dY/dx]}
        relay = {nid: {} for nid in self.network.nodes}

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
                up_states.append(state_up)

            T, Y, mix_res = node.compute_derived_state(up_states)
            self._derived_states[nid] = (T, Y, mix_res)

            if mix_res is None:
                continue

            # Sensitivity Relay (Chain-Rule)
            for i, elem in enumerate(up_elems):
                from_nid = elem.from_node
                t_jac = mix_res.dT_mix_d_stream[i]
                y_jac_list = mix_res.dY_mix_d_stream[i]
                # Direct dependency on upstream mass flow (if it's an unknown)
                m_indices = self._unknown_indices.get(elem.id)
                if m_indices:
                    idx = m_indices[0]
                    node_relay = relay[nid].setdefault(
                        idx, {"T": 0.0, "Y": np.zeros(cb.num_species())}
                    )
                    node_relay["T"] += t_jac.d_mdot
                    for k in range(len(y_jac_list)):
                        node_relay["Y"][k] += y_jac_list[k].d_mdot

                # Recursive dependency on upstream P/T/Y via from_node relay
                if from_nid in relay:
                    for idx, sens_up in relay[from_nid].items():
                        node_relay = relay[nid].setdefault(
                            idx, {"T": 0.0, "Y": np.zeros(cb.num_species())}
                        )
                        # dT/dx = sum_i( (dT/dT_in_i * dT_in_i/dx) + (dT/dY_in_i * dY_in_i/dx) )
                        node_relay["T"] += t_jac.d_T * sens_up["T"]
                        node_relay["T"] += np.dot(t_jac.d_Y, sens_up["Y"])

                        # dY_k/dx = sum_i( (dY_k/dT_in_i * dT_in_i/dx) + sum_j(dY_k/dY_in_j * dY_in_j/dx) )
                        for k in range(len(y_jac_list)):
                            node_relay["Y"][k] += y_jac_list[k].d_T * sens_up["T"]
                            node_relay["Y"][k] += np.dot(y_jac_list[k].d_Y, sens_up["Y"])

        return relay

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

        # DEBUG ATTACH NAMES
        if not hasattr(self, "_debug_eqn_names") or len(self._debug_eqn_names) != len(res):
            names = []
            for node in self.network.nodes.values():
                names.append(f"{node.id} mass_con")
            for node in self.network.nodes.values():
                if isinstance(node, cb.network.components.CombustorNode):
                    names.append(f"{node.id} t")
                    for i in range(14):
                        names.append(f"{node.id} X[{i}]")
                    names.append(f"{node.id} p_tot")
                    names.append(f"{node.id} p_stat")
                else:
                    s = self._get_node_state(node, x)
                    node_res, _ = node.residuals(s)
                    for k in range(len(node_res)):
                        names.append(f"{node.id} intrinsic[{k}]")
            for elem_id, element in self.network.elements.items():
                for k in range(element.n_equations()):
                    names.append(f"{elem_id} eq[{k}]")
            self._debug_eqn_names = names

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
