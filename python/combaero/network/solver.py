import math
import time
import warnings
from typing import Any, Literal

import numpy as np
import scipy.sparse as sp
from scipy.optimize import root

import combaero as cb

from .components import (
    ChannelElement,
    EnergyBoundary,
    LosslessConnectionElement,
    MassFlowBoundary,
    NetworkMixtureState,
    NetworkNode,
    PressureBoundary,
    WallNode,
)
from .graph import FlowNetwork


class SolverTimeoutError(Exception):
    """Raised when the solver exceeds the specified wall-clock timeout."""

    pass


class NetworkSolver:
    """
    Orchestrates the numerical solution of a fluid flow network using scipy.optimize.root.
    Translates the graph topology into a flat array of unknowns and residuals.
    Nodes define their Pressure (P) and Total Pressure (Pt) as unknowns,
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
        return name.endswith(".P") or name.endswith(".Pt")

    @staticmethod
    def _is_mdot_unknown(name: str) -> bool:
        return (
            name.endswith(".m_dot") or name.endswith(".m_dot_com") or name.endswith(".m_dot_branch")
        )

    @staticmethod
    def _to_incompressible_regime(regime: str) -> str | None:
        if regime == "compressible":
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

            # Node mass conservation equation appended in _residuals_and_jacobian,
            # unless the node is a junction port (then the MultiPortChamberElement's
            # sum-mass residual covers conservation; mass row is skipped).
            if not getattr(node, "_is_junction_port", False):
                scales.append(ref_mdot)

        from .components import TeeJunctionElement

        for element in self.network.elements.values():
            elem_scale = (
                ref_p
                if isinstance(
                    element,
                    (ChannelElement, LosslessConnectionElement, TeeJunctionElement),
                )
                else ref_mdot
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
        self._default_Y: list[float] = list(cb.species.dry_air_mass())
        self._n_species: int = len(self._default_Y)
        self._topological_order: list[str] = []
        self._derived_states: dict[str, tuple[float, list[float], Any]] = {}
        # Optional initial-guess overrides keyed by unknown name; consulted
        # by _build_x0 after user-provided initial_guess dicts and before
        # the topological p_guess / mdot_guess propagators. Populated by
        # init_strategy branches (e.g. analytical_pt_prop) and cleared
        # once the solve returns.
        self._init_overrides: dict[str, float] = {}

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
                pressures.append(getattr(node, "Pt", 101325.0))
                T = getattr(node, "Tt", None)
                if T is not None:
                    temperatures.append(T)
                Y = getattr(node, "Y", None)
                if Y is not None:
                    compositions.append(list(Y))

            elif isinstance(node, MassFlowBoundary):
                T = getattr(node, "Tt", None)
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
                p_guess[nid] = getattr(node, "Pt", ref["P"])

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

        # WallNode dead-ends carry no flow, so the correct pressure equals the
        # upstream node pressure (no friction drop).  Override the BFS estimate.
        for nid, node in self.network.nodes.items():
            if isinstance(node, WallNode):
                for elem in self.network.get_upstream_elements(nid):
                    p_guess[nid] = p_guess.get(elem.from_node, p_guess.get(nid, ref["P"]))
                    break

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
            # Only substitute a nonzero reference guess for elements whose
            # source is unknown.  If the share is zero because an upstream
            # MassFlowBoundary is genuinely OFF (m_dot=0), keep zero.
            if abs(share) < 1e-12 and not isinstance(node, MassFlowBoundary):
                share = ref["m_dot"]

            from_mass_boundary = isinstance(node, MassFlowBoundary)
            for elem in down_elems:
                elem_share = share
                regime = getattr(elem, "regime", None)

                # Dead-end channels (to_node is a WallNode) carry zero flow.
                # Start Newton at m_dot=0 so the mass-balance constraint is
                # trivially satisfied from the first iteration.
                if isinstance(self.network.nodes.get(elem.to_node), WallNode):
                    elem_share = 0.0

                # Only apply the Mach cap when the flow is estimated (pressure-
                # boundary seeded).  When a MassFlowBoundary seeds the element
                # directly, the mass flow is a known constraint and capping it
                # produces an initial guess that is too far from the solution,
                # causing hybr to converge to a spurious local minimum.
                elif regime == "compressible" and not from_mass_boundary:
                    area = getattr(elem, "area", None)
                    if area is not None and area > 0.0:
                        mdot_cap = rho_ref * a_ref * area * mdot_ratio
                        elem_share = min(abs(elem_share), mdot_cap)
                        if share < 0:
                            elem_share = -elem_share

                # Accumulate rather than overwrite: multi-source elements (e.g.
                # merging tee) are visited once per source node and each visit
                # contributes its share, so the total equals the sum of all inflows.
                mdot_guess[elem.id] = mdot_guess.get(elem.id, 0.0) + elem_share
                node_flow[elem.to_node] = node_flow.get(elem.to_node, 0.0) + elem_share

        return mdot_guess

    def _propagate_analytical_pt_prop(self, ref: dict[str, Any]) -> dict[str, float]:
        """Compute topology-aware initial guesses for MPCE compressible cold-start.

        For each ``ChannelElement`` with a signed pressure gradient across
        it, seed the m_dot unknown via a Bernoulli estimate
        ``m_dot = sign(dP) * A * sqrt(2 * rho_ref * |dP|)``.

        For each ``MultiPortChamberElement`` (and subclasses), seed the
        ``P_jct`` unknown from the propagated Pt at the junction's
        "common" port (the arm on the single-inlet side for a branch or
        the single-outlet side for a merge). The default
        ``ref["P"] = min(boundary Pt)`` is almost never near the true
        junction static, which is what motivated this strategy (see
        LHS-32 audit, ``tmp/mpce_cold_start_audit.py``).

        Returns a dict mapping unknown name to seed value. Consumed by
        ``_build_x0`` via ``self._init_overrides``.
        """
        from .components import MultiPortChamberElement

        p_guess = self._propagate_pressure_guess(ref)
        if not p_guess:
            return {}
        # Reference density from Pt/RT at the propagated mean pressure and
        # ambient reference temperature. R_air is a stand-in for the local
        # gas; MPCE audit shows this level of approximation is enough to
        # get Newton into the correct basin.
        R_air = 287.0
        Tt = float(ref.get("T", 300.0))
        Pt_ref = float(np.mean(list(p_guess.values())))
        rho_ref = max(Pt_ref / (R_air * Tt), 1e-3)

        overrides: dict[str, float] = {}
        for elem_id, elem in self.network.elements.items():
            if isinstance(elem, ChannelElement):
                pt_from = p_guess.get(elem.from_node)
                pt_to = p_guess.get(elem.to_node)
                if pt_from is None or pt_to is None:
                    continue
                dp = pt_from - pt_to
                if dp == 0.0:
                    continue
                diameter = getattr(elem, "diameter", None)
                if diameter is None or diameter <= 0.0:
                    continue
                area = math.pi * (float(diameter) / 2.0) ** 2
                mdot_mag = area * math.sqrt(2.0 * rho_ref * abs(dp))
                overrides[f"{elem_id}.m_dot"] = math.copysign(mdot_mag, dp)
            elif isinstance(elem, MultiPortChamberElement):
                inlets = list(getattr(elem, "inlet_nodes", []))
                outlets = list(getattr(elem, "outlet_nodes", []))
                if len(inlets) == 1:
                    common = inlets[0]
                elif len(outlets) == 1:
                    common = outlets[0]
                else:
                    continue
                pt = p_guess.get(common)
                if pt is None:
                    continue
                overrides[f"{elem_id}.P_jct"] = pt
        return overrides

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

        # Bernoulli-based initial guess for TeeJunctionElement unknowns.
        # The generic mdot_guess propagator seeds from ref["m_dot"]=0.1 (arbitrary
        # when no MassFlowBoundary exists), but tee flows scale as A*sqrt(2*rho*dP)
        # which can be O(1) kg/s for typical pressure differences.  A factor-of-10
        # error in the starting mass flow prevents Newton from converging.
        from .components import TeeJunctionElement as _TeeJE

        _ref_X_tee = list(cb.mass_to_mole(ref["Y"]))
        _rho_tee = float(cb.density(ref["T"], ref["P"], _ref_X_tee))
        _a_tee = float(cb.speed_of_sound(ref["T"], _ref_X_tee))
        _gamma_tee = float(cb.isentropic_expansion_coefficient(ref["T"], _ref_X_tee))
        _mach_cap_target = 0.3
        _term_tee = 1.0 + 0.5 * (_gamma_tee - 1.0) * _mach_cap_target**2
        _exp_tee = -(_gamma_tee + 1.0) / (2.0 * (_gamma_tee - 1.0))
        _mdot_ratio_tee = _mach_cap_target * _term_tee**_exp_tee

        def _node_pt(n: str) -> float:
            node = self.network.nodes[n]
            if isinstance(node, PressureBoundary):
                return float(node.Pt)
            return p_guess.get(n, ref["P"])

        _tee_branch_guess: dict[str, float] = {}
        for _eid, _elem in self.network.elements.items():
            if not isinstance(_elem, _TeeJE):
                continue
            _src = _elem.all_source_nodes()
            _snk = _elem.all_sink_nodes()
            _pt_src = max(_node_pt(n) for n in _src)
            _pt_snk = min(_node_pt(n) for n in _snk)
            _dP_total = max(_pt_src - _pt_snk, 1.0)
            # F_C may be None before topology is resolved (e.g. when _build_x0
            # is called directly before NetworkSolver.solve resolves topology).
            # Skip the Bernoulli estimate and fall back to the propagated guess.
            if not _elem.F_C:
                continue
            # Mach cap: the tee duct cannot carry more than M=0.3 at reference
            # conditions.  Without this cap the Bernoulli estimate can overestimate
            # the tee flow by 3x when downstream channel losses are included in the
            # boundary-to-boundary dP, causing Newton to start at a wildly
            # inconsistent mass balance.
            _mdot_cap_tee = _rho_tee * _a_tee * _elem.F_C * _mdot_ratio_tee
            _bernoulli_com = min(
                _elem.F_C * math.sqrt(2.0 * _rho_tee * _dP_total),
                _mdot_cap_tee,
            )
            # Prefer a propagated mdot (e.g. from MassFlowBoundary) over the
            # Bernoulli estimate when the former is larger -- the propagated value
            # is exact when a mass-flow BC seeds the network, while the Bernoulli
            # estimate can be a factor of ~6 too small when only one pressure
            # boundary exists (dP from p_guess is only 0.1% of ref_P).
            mdot_guess[_eid] = max(_bernoulli_com, mdot_guess.get(_eid, 0.0))
            # Branch guess: Bernoulli estimate for the branch arm pressure drop.
            _pt_br = _node_pt(_elem.branch_node)
            if _elem.tee_type == "merging":
                _dP_br = max(_pt_br - _pt_snk, 1.0)
            else:
                _dP_br = max(_pt_src - _pt_br, 1.0)
            _bernoulli_br = min(
                _elem.F_C * math.sqrt(2.0 * _rho_tee * _dP_br),
                mdot_guess[_eid] * 0.75,
            )
            _tee_branch_guess[_eid] = max(_bernoulli_br, mdot_guess[_eid] * 0.25)

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
                        elif unk in self._init_overrides:
                            x0_list.append(self._init_overrides[unk])
                        elif unk.endswith(".P") or unk.endswith(".Pt"):
                            x0_list.append(p_guess.get(node_id, ref["P"]))
                        elif unk.endswith(".T") or unk.endswith(".Tt"):
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
                    elif unk in self._init_overrides:
                        x0_list.append(self._init_overrides[unk])
                    elif unk.endswith(".m_dot") or unk.endswith(".m_dot_com"):
                        x0_list.append(mdot_guess.get(elem_id, ref["m_dot"]))
                    elif unk.endswith(".m_dot_branch"):
                        x0_list.append(
                            _tee_branch_guess.get(
                                elem_id, mdot_guess.get(elem_id, ref["m_dot"]) * 0.5
                            )
                        )
                    elif unk.endswith(".P_jct"):
                        # MultiPortChamberElement: P_jct should be close to the
                        # ambient port static pressure. Use the max boundary Pt
                        # as a stand-in (P_jct >= Pt at every port since the
                        # impulse equation gives P_jct = Pt + q with q >= 0).
                        x0_list.append(ref["P"])
                    else:
                        x0_list.append(1.0)
                self._unknown_indices[elem_id] = list(range(start_idx, len(x0_list)))

        # Build name to index mapping
        self._name_to_index = {name: i for i, name in enumerate(self.unknown_names)}

        # MultiPortChamberElement ports carry their throughflow on the outer
        # connecting elements' m_dot unknowns. Stash the resolved global
        # indices on each junction so flow_at_node / flow_jac_at_node can
        # report real port flows during state propagation: collector-port
        # MCNs need them for the Pt = P + 0.5*rho*v^2 closure and for
        # mass-weighted mixing of merging streams.
        from .components import MultiPortChamberElement as _MPCElem

        for element in self.network.elements.values():
            if isinstance(element, _MPCElem):
                element._port_outer_mdot_idx = {
                    outer_id: self._name_to_index[f"{outer_id}.m_dot"]
                    for outer_id in element._port_element_ids
                    if f"{outer_id}.m_dot" in self._name_to_index
                }

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
                for to_id in elem.all_sink_nodes():
                    if to_id not in visited:
                        up_elems = self.network.get_upstream_elements(to_id)
                        # Ready when all source nodes of all upstream elements are visited
                        if all(src in visited for e in up_elems for src in e.all_source_nodes()):
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

        def _endpoint_to_node_id(endpoint_id: str) -> str | None:
            if endpoint_id in self.network.elements:
                return self.network.elements[endpoint_id].to_node
            if endpoint_id in self.network.nodes:
                return endpoint_id
            return None

        # Build a mapping: node_id -> set of wall IDs that affect it
        affected_nodes: dict[str, set[str]] = {}
        for wall in self.network.walls.values():
            node_a = _endpoint_to_node_id(wall.element_a)
            node_b = _endpoint_to_node_id(wall.element_b)

            if node_a is None or node_b is None:
                continue

            # Skip walls where both sides feed the same node (net Q = 0)
            if node_a == node_b:
                continue

            affected_nodes.setdefault(node_a, set()).add(wall.id)
            affected_nodes.setdefault(node_b, set()).add(wall.id)

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

        for wall in self.network.walls.values():
            if wall.id in processed_walls:
                continue

            # Endpoints may be elements or nodes
            elem_a_id = wall.element_a
            elem_b_id = wall.element_b

            endpoint_a_is_node = elem_a_id in self.network.nodes
            endpoint_b_is_node = elem_b_id in self.network.nodes

            obj_a = self.network.elements.get(elem_a_id) or self.network.nodes.get(elem_a_id)
            obj_b = self.network.elements.get(elem_b_id) or self.network.nodes.get(elem_b_id)

            if obj_a is None or obj_b is None:
                continue

            # Find the actual target node for side A and side B
            target_node_a = obj_a.id if endpoint_a_is_node else obj_a.to_node
            target_node_b = obj_b.id if endpoint_b_is_node else obj_b.to_node

            # This wall does not affect the current node
            if nid not in (target_node_a, target_node_b):
                continue

            # Skip if both sides feed/are the same node (net Q = 0)
            if target_node_a == target_node_b:
                continue

            # Determine which side feeds the current node
            is_side_a = target_node_a == nid

            # Get states for both sides (previous-iteration fallback for back-edges)
            node_feeding_a = (
                self.network.nodes[obj_a.id]
                if endpoint_a_is_node
                else self.network.nodes[obj_a.from_node]
            )
            node_feeding_b = (
                self.network.nodes[obj_b.id]
                if endpoint_b_is_node
                else self.network.nodes[obj_b.from_node]
            )

            state_a = self._get_node_state_with_prev(node_feeding_a, x)
            state_b = self._get_node_state_with_prev(node_feeding_b, x)

            # Set mass flow rates from solver vector for elements
            if not endpoint_a_is_node:
                m_indices_a = self._unknown_indices.get(elem_a_id, [])
                if m_indices_a:
                    state_a.m_dot = x[m_indices_a[0]]
            if not endpoint_b_is_node:
                m_indices_b = self._unknown_indices.get(elem_b_id, [])
                if m_indices_b:
                    state_b.m_dot = x[m_indices_b[0]]

            # Call htc_and_T() on both objects
            ch_result_a = obj_a.htc_and_T(state_a) if obj_a.has_convective_surface else None
            ch_result_b = obj_b.htc_and_T(state_b) if obj_b.has_convective_surface else None

            # Only proceed if both elements have heat transfer surfaces
            if ch_result_a is not None and ch_result_b is not None:
                # Extract scalars from ChannelResult
                h_a = ch_result_a.h
                T_aw_a = ch_result_a.T_aw
                h_b = ch_result_b.h
                T_aw_b = ch_result_b.T_aw

                # Call multi-layer coupling logic
                wall_result = wall.compute_coupling(
                    h_a, T_aw_a, obj_a.surface.area, h_b, T_aw_b, obj_b.surface.area
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
        # relay[node_id][global_unknown_index] = {'T': dT/dx, 'Y': [dY/dx], 'Pt_mix': dP_total_mix/dx}
        relay = {nid: {} for nid in self.network.nodes}

        has_walls = self.network.thermal_coupling_enabled and bool(self.network.walls)

        for nid in self._topological_order:
            node = self.network.nodes[nid]

            up_elems = self.network.get_upstream_elements(nid)

            # Boundaries with NO upstream connections are pure sources
            if isinstance(node, (PressureBoundary, MassFlowBoundary)) and not up_elems:
                T, Y, _ = node.compute_derived_state([])
                self._derived_states[nid] = (T, Y, None)
                # Seed relay for boundary's own unknowns
                node_indices = self._unknown_indices.get(nid, [])
                unknown_names = self.unknown_names
                for global_idx in node_indices:
                    name = unknown_names[global_idx]
                    if name.endswith(".Pt"):
                        relay[nid].setdefault(
                            global_idx, {"T": 0.0, "Y": np.zeros(self._n_species), "Pt": 0.0}
                        )["Pt"] = 1.0
                continue

            # stream_info: (state, elem, src_nid) - one entry per upstream stream.
            # Multi-port elements (e.g. TeeJunctionElement) contribute multiple entries.
            stream_info: list[tuple] = []
            for elem in up_elems:
                indices = self._unknown_indices.get(elem.id, [])
                # Compute named m_dot Jacobian for this element once (flow into nid).
                # Stored on each upstream state so MomentumChamberNode can build the
                # correct d_res/dm_dot Jacobian even when the element has named unknowns
                # like tee.m_dot_com / tee.m_dot_branch instead of a plain elem.m_dot.
                elem_jac_names: dict[str, float] = {}
                if indices:
                    for col, coeff in elem.flow_jac_at_node(nid, indices).items():
                        if col < len(self.unknown_names):
                            elem_jac_names[self.unknown_names[col]] = coeff
                _single_source = len(elem.all_source_nodes()) == 1
                for src_nid in elem.all_source_nodes():
                    state_up = self._get_node_state(self.network.nodes[src_nid], x)
                    if indices:
                        # For single-source elements (channels, branching tees) use the
                        # flow INTO the current sink node (nid) so that each sink's
                        # _total_m_dot reflects its actual partial flow, not the
                        # supplier's total.  For multi-source elements (merging tees) each
                        # source contributes its own partial flow, so keep flow_at_node
                        # called with src_nid.
                        if _single_source:
                            state_up.m_dot = elem.flow_at_node(nid, x, indices)
                        else:
                            state_up.m_dot = elem.flow_at_node(src_nid, x, indices)
                        state_up._element_id = elem.id
                        state_up._m_dot_jac_names = elem_jac_names
                    stream_info.append((state_up, elem, src_nid))
            up_states = [s for s, _, _ in stream_info]

            # Wall Coupling: evaluate walls BEFORE compute_derived_state
            # so that Q flows through the EnergyBoundary abstraction.
            wall_contributions: list[tuple] = []
            if has_walls:
                wall_contributions = self._evaluate_walls_for_node(nid, up_elems, x)

            T, Y, mix_res = node.compute_derived_state(up_states)
            T = min(max(float(T), 200.0), 5000.0)
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
                n_species = self._n_species
                # Iterate over stream_info (one entry per upstream stream, including
                # multiple entries for multi-port elements like TeeJunctionElement).
                for i, (_, elem, src_nid) in enumerate(stream_info):
                    t_jac = mix_res.dT_mix_d_stream[i]
                    pt_jac = mix_res.dP_total_mix_d_stream[i]
                    # dY_mix_d_stream is indexed by [species][stream]
                    y_jacs = [mix_res.dY_mix_d_stream[k][i] for k in range(n_species)]

                    # 1. Direct dependency on upstream mass flow (if it's a solver unknown)
                    m_indices = self._unknown_indices.get(elem.id)
                    if m_indices:
                        # flow_jac_at_node maps each unknown index to its coefficient
                        # in d(stream_i.m_dot)/d(x). For 2-port: {idx: 1.0}. For tee:
                        # may involve 2 unknowns (m_dot_com and m_dot_branch).
                        for idx, coeff in elem.flow_jac_at_node(src_nid, m_indices).items():
                            node_relay = relay[nid].setdefault(
                                idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt_mix": 0.0}
                            )
                            node_relay["T"] += t_jac.d_mdot * coeff
                            node_relay["Pt_mix"] += pt_jac.d_mdot * coeff
                            for k in range(n_species):
                                node_relay["Y"][k] += y_jacs[k].d_mdot * coeff

                    # 2. Recursive dependency on upstream P/T/Y via src_nid relay
                    if src_nid in relay:
                        for idx, sens_up in relay[src_nid].items():
                            node_relay = relay[nid].setdefault(
                                idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt_mix": 0.0}
                            )
                            # dT/dx = sum_i( dT/dT_in_i * dT_in_i/dx + dT/dP_in_i * dP_in_i/dx + dT/dY_in_i * dY_in_i/dx )
                            node_relay["T"] += t_jac.d_T * sens_up["T"]
                            # Note: stagnation Pt doesn't usually have T/Y dependency, but we relay it if present
                            node_relay["T"] += t_jac.d_P_total * sens_up.get("Pt", 0.0)
                            node_relay["T"] += np.dot(t_jac.d_Y, sens_up["Y"])

                            # dP_total_mix/dx
                            node_relay["Pt_mix"] += pt_jac.d_T * sens_up["T"]
                            node_relay["Pt_mix"] += pt_jac.d_P_total * sens_up.get("Pt", 0.0)
                            node_relay["Pt_mix"] += np.dot(pt_jac.d_Y, sens_up["Y"])

                            # dY_k/dx
                            for k in range(n_species):
                                node_relay["Y"][k] += y_jacs[k].d_T * sens_up["T"]
                                node_relay["Y"][k] += y_jacs[k].d_P_total * sens_up.get("Pt", 0.0)
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
                        obj_a = self.network.elements.get(ea_id) or self.network.nodes.get(ea_id)
                        from_a = ea_id if ea_id in self.network.nodes else obj_a.from_node
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
                                    idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0}
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
                                idx_mdot_a, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0}
                            )
                            nr["T"] += factor_mdot_a

                        # Side B contributions
                        obj_b = self.network.elements.get(eb_id) or self.network.nodes.get(eb_id)
                        from_b = eb_id if eb_id in self.network.nodes else obj_b.from_node
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
                                    idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0}
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
                                idx_mdot_b, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0}
                            )
                            nr["T"] += factor_mdot_b

            # 3. Add own unknowns (Pt) to the relay
            node_unks = self._unknown_indices.get(nid, [])
            for unk_idx in node_unks:
                unk_name = self.unknown_names[unk_idx]
                if unk_name.endswith(".Pt"):
                    node_relay = relay[nid].setdefault(
                        unk_idx, {"T": 0.0, "Y": np.zeros(self._n_species), "Pt": 0.0}
                    )
                    node_relay["Pt"] = 1.0

        return relay

    def _get_node_state_with_prev(self, node: NetworkNode, x: np.ndarray) -> NetworkMixtureState:
        """
        Constructs a NetworkMixtureState for a node, using previous iteration's derived state
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
            "Pt": 101325.0,
            "T": T_prev,
            "Tt": T_prev,
            "m_dot": 0.0,
            "Y": Y_prev,
        }

        for i, unk in zip(indices, unknowns, strict=False):
            var_name = unk.split(".")[-1]
            if var_name in state_dict:
                state_dict[var_name] = x[i]

        return NetworkMixtureState(
            float(state_dict["P"]),
            float(state_dict["Pt"]),
            float(state_dict["T"]),
            float(state_dict["Tt"]),
            float(state_dict["m_dot"]),
            state_dict["Y"],
        )

    def _get_node_state(self, node: NetworkNode, x: np.ndarray) -> NetworkMixtureState:
        """Constructs a NetworkMixtureState for a given node based on the current solver vector x."""
        default_Y = self._default_Y

        # Determine boundaries
        if isinstance(node, PressureBoundary):
            guess = getattr(node, "initial_guess", {})
            P_tot = guess.get(f"{node.id}.Pt", getattr(node, "Pt", 101325.0))
            P_stat = guess.get(f"{node.id}.P", P_tot)
            T_tot = guess.get(f"{node.id}.Tt", getattr(node, "Tt", 300.0))
            T_stat = guess.get(f"{node.id}.T", T_tot)
            Y = getattr(node, "Y", None)
            if Y is None:
                Y = default_Y
            return NetworkMixtureState(
                float(P_stat), float(P_tot), float(T_stat), float(T_tot), 0.0, Y
            )

        if isinstance(node, MassFlowBoundary):
            guess = getattr(node, "initial_guess", {})
            m = getattr(node, "m_dot", 0.0)
            T_tot = guess.get(f"{node.id}.Tt", getattr(node, "Tt", 300.0))
            Y = getattr(node, "Y", None)

            # If this boundary has upstream connections (sink), use derived mixed state
            if node.id in self._derived_states:
                T_tot_derived, Y_derived, _ = self._derived_states[node.id]
                T_tot = T_tot_derived
                Y = Y_derived

            if Y is None:
                Y = default_Y

            T_stat = guess.get(f"{node.id}.T", T_tot)
            p_val = guess.get(f"{node.id}.P", 101325.0)
            ptot_val = guess.get(f"{node.id}.Pt", 101325.0)

            indices = self._unknown_indices.get(node.id, [])
            unknowns = node.unknowns()
            for i, unk in zip(indices, unknowns, strict=False):
                if unk.endswith(".P"):
                    p_val = x[i]
                elif unk.endswith(".Pt"):
                    ptot_val = x[i]

            return NetworkMixtureState(
                float(p_val), float(ptot_val), float(T_stat), float(T_tot), float(m), Y
            )

        # For interior nodes, unpack from x
        indices = self._unknown_indices.get(node.id, [])
        unknowns = node.unknowns()

        # Default state
        state_dict = {
            "P": 101325.0,
            "Pt": 101325.0,
            "T": 300.0,
            "Tt": 300.0,
            "m_dot": 0.0,
            "Y": default_Y,
        }

        # Override with derived state if available
        if node.id in self._derived_states:
            T, Y, _ = self._derived_states[node.id]
            state_dict["T"] = T
            state_dict["Tt"] = T
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

        # P5: Floor on static pressure to prevent NaN in thermo calls at near-zero P.
        P_safe = max(float(state_dict["P"]), 1000.0)

        return NetworkMixtureState(
            P_safe,
            float(state_dict["Pt"]),
            min(max(float(state_dict["T"]), 200.0), 5000.0),
            min(max(float(state_dict["Tt"]), 200.0), 5000.0),
            float(state_dict["m_dot"]),
            Y_safe,
        )

    def _residuals_and_jacobian(
        self, x: np.ndarray, *, compute_jacobian: bool = True
    ) -> tuple[np.ndarray, Any]:
        """
        Evaluates the global residual vector and its sparse Jacobian.
        """
        # Pure Pressure-Flow: Propagate states forward first
        relay = self._propagate_states(x)

        res = []
        rows = [] if compute_jacobian else None
        cols = [] if compute_jacobian else None
        data = [] if compute_jacobian else None

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
            if compute_jacobian:
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
                                    elif prop == "Pt_mix" and "Pt_mix" in sens_pkg:
                                        rows.append(row)
                                        cols.append(unk_idx)
                                        data.append(deriv * sens_pkg["Pt_mix"])

            # Mass Conservation: Sum(m_dot_in) - Sum(m_dot_out) = 0
            #
            # Skip for port-MCNs of a MultiPortChamberElement: the junction's
            # sum-mdot residual replaces the per-port mass row. Without this
            # skip the row degenerates to 0=0 (both the channel and the junction
            # report the same flow at the node with opposite signs, producing a
            # zero row in the Jacobian).
            if getattr(node, "_is_junction_port", False):
                continue
            mass_res_idx = len(res)
            upstream_elems = self.network.get_upstream_elements(node_id)
            downstream_elems = self.network.get_downstream_elements(node_id)

            m_dot_in = 0.0
            m_dot_out = 0.0

            if isinstance(node, MassFlowBoundary):
                # If connected as a sink (has upstream, no downstream), subtract from balance.
                # Otherwise (source or injection) add to balance.
                if upstream_elems and not downstream_elems:
                    m_dot_out += getattr(node, "m_dot", 0.0)
                else:
                    m_dot_in += getattr(node, "m_dot", 0.0)

            for elem in upstream_elems:
                indices = self._unknown_indices.get(elem.id)
                if indices:
                    m_dot_in += elem.flow_at_node(node_id, x, indices)
                    if compute_jacobian:
                        for col, coeff in elem.flow_jac_at_node(node_id, indices).items():
                            rows.append(mass_res_idx)
                            cols.append(col)
                            data.append(coeff)
                else:
                    from_node = self.network.nodes[elem.from_node]
                    m_dot_in += getattr(from_node, "m_dot", 0.0)

            for elem in downstream_elems:
                indices = self._unknown_indices.get(elem.id)
                if indices:
                    m_dot_out += elem.flow_at_node(node_id, x, indices)
                    if compute_jacobian:
                        for col, coeff in elem.flow_jac_at_node(node_id, indices).items():
                            rows.append(mass_res_idx)
                            cols.append(col)
                            data.append(-coeff)
                else:
                    to_node = self.network.nodes[elem.to_node]
                    m_dot_out += getattr(to_node, "m_dot", 0.0)

            res.append(m_dot_in - m_dot_out)

        # 2. Element Residuals
        from .components import MultiPortChamberElement, TeeJunctionElement

        for elem_id, element in self.network.elements.items():
            start_res_idx = len(res)

            m_indices = self._unknown_indices.get(elem_id)

            if isinstance(element, TeeJunctionElement):
                nodes = self.network.nodes
                state_com = self._get_node_state(nodes[element.common_node], x)
                state_straight = self._get_node_state(nodes[element.straight_node], x)
                state_branch = self._get_node_state(nodes[element.branch_node], x)
                if m_indices:
                    state_com.m_dot = float(x[m_indices[0]])
                    state_branch.m_dot = float(x[m_indices[1]])
                    state_straight.m_dot = float(x[m_indices[0]]) - float(x[m_indices[1]])
                elem_res, elem_jac = element.residuals(state_com, state_straight, state_branch)
            elif isinstance(element, MultiPortChamberElement):
                # N-port momentum-CV junction: one P_jct unknown plus N+1 residuals.
                # Port mdots are sourced from each connecting element's m_dot
                # unknown, sign-mapped to junction convention (positive = out).
                nodes = self.network.nodes
                port_states = [self._get_node_state(nodes[pid], x) for pid in element.port_nodes]
                P_jct_val = float(x[m_indices[0]]) if m_indices else 0.0
                port_mdots: list[float] = []
                for i, outer_id in enumerate(element._port_element_ids):
                    outer_indices = self._unknown_indices.get(outer_id, [])
                    outer_mdot = float(x[outer_indices[0]]) if outer_indices else 0.0
                    port_mdots.append(element._port_signs[i] * outer_mdot)
                elem_res, elem_jac = element.residuals(port_states, P_jct_val, port_mdots)
            else:
                node_in = self.network.nodes[element.from_node]
                node_out = self.network.nodes[element.to_node]

                state_in = self._get_node_state(node_in, x)
                state_out = self._get_node_state(node_out, x)

                # Unpack element-specific unknowns (like m_dot) into the state
                if m_indices:
                    state_in.m_dot = x[m_indices[0]]
                    state_out.m_dot = x[m_indices[0]]

                elem_res, elem_jac = element.residuals(state_in, state_out)
            res.extend(elem_res)

            if compute_jacobian:
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
                                    elif prop == "Pt_mix" and "Pt_mix" in sens_pkg:
                                        rows.append(row)
                                        cols.append(unk_idx)
                                        data.append(deriv * sens_pkg["Pt_mix"])
                                    elif prop == "Pt" and "Pt" in sens_pkg:
                                        rows.append(row)
                                        cols.append(unk_idx)
                                        data.append(deriv * sens_pkg["Pt"])

        jac_sparse = None
        if compute_jacobian:
            jac_sparse = sp.coo_matrix((data, (rows, cols)), shape=(len(res), len(x))).tocsr()

        return np.array(res, dtype=float), jac_sparse

    def _residuals(self, x: np.ndarray) -> np.ndarray:
        """
        Evaluates the global residual vector given the current unknown vector x.
        """
        res, _ = self._residuals_and_jacobian(x, compute_jacobian=False)
        return res

    def solve(
        self,
        method: str = "hybr",
        timeout: float | None = None,
        options: dict[str, Any] | None = None,
        use_jac: bool = True,
        x0: np.ndarray | None = None,
        init_strategy: Literal[
            "default", "incompressible_warmstart", "homotopy", "analytical_pt_prop"
        ] = "default",
        warmstart_maxfev: int = 150,
        lambda_steps: list[float] | None = None,
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
                direct x0 construction; for networks containing a
                ``MultiPortChamberElement`` it auto-upgrades to
                ``analytical_pt_prop`` (32/32 vs 28/32 certified-root
                convergence on the 2026-07 inverse-design audit,
                ``tmp/mpce_audit_v2_runner.py``) -- pass another
                strategy or an explicit ``x0`` to opt out.
                ``analytical_pt_prop`` seeds topology-aware initial
                guesses (channel Bernoulli m_dot + MPCE junction P_jct
                at the common port). ``incompressible_warmstart``
                (DEPRECATED: 10/32 on the same audit at ~10x the wall
                time) first solves an incompressible-regime proxy
                network and uses that solution as x0.
            warmstart_maxfev: Maximum function evaluations for the
                incompressible proxy solve when
                ``init_strategy='incompressible_warmstart'``.
            lambda_steps: Optional list of mass flow scaling factors for
                ``init_strategy='homotopy'``. Defaults to 10 steps
                from 0.1 to 1.0.
        """
        if method not in self.SUPPORTED_METHODS:
            raise ValueError(
                f"Method '{method}' is not supported. Supported methods are: {', '.join(self.SUPPORTED_METHODS)}"
            )

        if init_strategy not in (
            "default",
            "incompressible_warmstart",
            "homotopy",
            "continuation",
            "analytical_pt_prop",
        ):
            raise ValueError(
                "init_strategy must be one of: 'default', 'incompressible_warmstart', "
                "'homotopy', 'continuation', 'analytical_pt_prop'."
            )
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

        if init_strategy == "incompressible_warmstart":
            warnings.warn(
                "init_strategy='incompressible_warmstart' is deprecated. "
                "On the certified MPCE cold-start audit (2026-07) it "
                "converged 10/32 cases at ~10x the wall time of "
                "'analytical_pt_prop' (32/32). Use 'analytical_pt_prop' "
                "or leave init_strategy='default'.",
                DeprecationWarning,
                stacklevel=2,
            )

        # Ensure topology is resolved before solving
        self.network.resolve_all_topology()
        self.network.validate()

        # MPCE networks auto-upgrade the default initialization to the
        # analytical Pt-propagation seeds: on the certified inverse-design
        # audit (tmp/mpce_audit_v2_*, 2026-07-04) analytical_pt_prop
        # converged 32/32 to the certified root vs 28/32 for the plain
        # cold start, and was the fastest strategy overall. Explicitly
        # passing any other init_strategy (or an x0) opts out.
        # The incompressible_warmstart proxy solve passes "default" and
        # must keep the legacy plain cold start (flag below), otherwise
        # the deprecated strategy's behavior silently changes.
        if (
            x0 is None
            and init_strategy == "default"
            and not getattr(self, "_in_warmstart_proxy", False)
        ):
            from .components import MultiPortChamberElement as _MPCElem

            if any(isinstance(e, _MPCElem) for e in self.network.elements.values()):
                init_strategy = "analytical_pt_prop"

        # analytical_pt_prop seeds channel m_dot + MPCE P_jct through the
        # regular _build_x0 path via self._init_overrides. _build_x0 reads
        # the dict once and clears it, so no cross-solve leakage.
        if x0 is None and init_strategy == "analytical_pt_prop":
            _ref = self._infer_reference_state()
            self._init_overrides = self._propagate_analytical_pt_prop(_ref)
        else:
            self._init_overrides = {}

        # Build initial guess to populate unknown_names early
        x0_auto = self._build_x0()
        self._init_overrides = {}
        if not self.unknown_names:
            return {
                "__success__": True,
                "__message__": "Network is empty or fully constrained (no unknowns).",
                "__iterations__": 0,
            }

        # Set default iteration limit based on problem size (matches hybr's own default
        # of 200*(n+1)) now that we know n.  Applied after the early-return check so
        # the network can't be empty when we read len(unknown_names).
        _n_unk = len(self.unknown_names)
        _default_iters = max(500, 200 * (_n_unk + 1))
        if method == "hybr":
            options.setdefault("maxfev", _default_iters)
            # Reduce initial trust-region step bound (hybr default is 100).
            # Smaller steps prevent large oscillating Newton steps near
            # ill-conditioned Jacobians (e.g., compressible flow near choking).
            options.setdefault("factor", 10)
        else:
            options.setdefault("maxiter", _default_iters)

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
                    self._in_warmstart_proxy = True
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", DeprecationWarning)
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
                    self._in_warmstart_proxy = False
                    for elem_id, regime in original_regimes.items():
                        self.network.elements[elem_id].regime = regime
                    self.network.resolve_all_topology()

        if x0 is None and init_strategy == "homotopy":
            from .components import MassFlowBoundary

            # Identify all MassFlowBoundary nodes and their target m_dot
            targets = {
                nid: node.m_dot
                for nid, node in self.network.nodes.items()
                if isinstance(node, MassFlowBoundary)
            }

            if targets:
                # Sequential solve loop with increasing mass flow lambda.
                # Use more steps for better robustness near choking.
                if lambda_steps is None:
                    lambda_steps = [0.2, 0.4, 0.6, 0.8, 0.9, 1.0]

                current_x0 = None
                last_sol = {"__success__": False}

                for lam in lambda_steps:
                    # Update boundaries
                    for nid, target_mdot in targets.items():
                        self.network.nodes[nid].m_dot = lam * target_mdot

                    # Solve intermediate step.
                    # Limit maxfev/maxiter for intermediate steps to stay fast.
                    # This is decremental to convergence: trade speed versus robustness near choking.
                    step_options = dict(options)
                    # options already has maxfev/_default_iters set above; honour
                    # it.  The setdefault here is only reached if the caller
                    # explicitly zeroed out maxfev/maxiter, which should not happen.

                    last_sol = self.solve(
                        method=method,
                        timeout=timeout,
                        options=step_options,
                        use_jac=use_jac,
                        x0=current_x0,
                        init_strategy="default",
                    )

                    # Update x0 for next step anyway (best guess propagation)
                    current_x0 = self._last_x

                    if not last_sol["__success__"]:
                        # If an intermediate step fails badly, we might still
                        # want to continue if progress was made, but for now
                        # we stop to avoid compounding errors.
                        # Note: we still use the failed result as x0 if it's the best we have.
                        break

                # Restore original m_dot (lam=1.0 loop already did this if it reached it)
                for nid, target_mdot in targets.items():
                    self.network.nodes[nid].m_dot = target_mdot

                # If homotopy reached lam=1.0 and succeeded, return that solution
                if last_sol["__success__"]:
                    return last_sol

                # If it failed at lam=1.0 but current_x0 is much better than cold start,
                # we will naturally use it as warmstart_x0 below.
                warmstart_x0 = current_x0
        elif x0 is None and init_strategy == "continuation":
            raise ValueError(
                "continuation init: no prior converged solution provided via 'x0'. "
                "The caller is responsible for persisting and passing the previous 'x' vector."
            )

        # --- Initial guess: user-supplied (warm-start) or auto-built ----
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
        # first residual/scaling pass by ensuring the first residual evaluation
        # triggers the propagation internally.

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

        # For compressible networks hybr can oscillate near choking.  Reserve
        # 35% of the wall-clock budget for an LM fallback that has better
        # global convergence via (J^TJ + lI) regularisation.
        # Incompressible inner solves (regime overrides applied) have no
        # compressible elements and keep the full budget.
        # Tee junctions couple two supplier stagnation pressures whose
        # common mode leaves hybr's Jacobian near-singular; LM's regularised
        # step resolves it, so route tee networks through the same fallback.
        # Same near-singularity affects MultiPortChamberElement networks (the
        # P_jct unknown couples N port-pressure rows that share a common-mode
        # direction at low-Mach).
        from .components import MultiPortChamberElement as _MPCElem
        from .components import TeeJunctionElement as _TeeJE

        _has_tee = any(isinstance(e, _TeeJE) for e in self.network.elements.values())
        _has_mpce = any(isinstance(e, _MPCElem) for e in self.network.elements.values())
        _has_compressible = method == "hybr" and (
            bool(self._compressible_element_overrides()) or _has_tee or _has_mpce
        )
        _hybr_budget = timeout * 0.65 if (_has_compressible and timeout is not None) else timeout
        # Mutable so the LM fallback block can extend the effective limit.
        _timeout_eff: list[float | None] = [_hybr_budget]

        # Convergence history: record on every ||F|| improvement and at most
        # once per 0.1 s to show oscillation without excessive overhead.
        # Cap at 1 000 entries (~50 KB JSON) regardless of problem size.
        _history: list[dict] = []
        _last_record_t: list[float] = [0.0]
        _RECORD_INTERVAL = 0.1
        _eval_count: list[int] = [0]
        _lm_started_at_eval: list[int | None] = [None]

        # State for timeout and best iterate tracking (in real space)
        start_time = time.perf_counter()
        best_x = x0_use.copy()
        best_res_norm = np.linalg.norm(res0)
        last_exception: Exception | None = None

        def residuals_wrapper(x_scaled: np.ndarray) -> Any:
            nonlocal best_x, best_res_norm, last_exception

            # Check timeout
            if _timeout_eff[0] is not None and (time.perf_counter() - start_time) > _timeout_eff[0]:
                raise SolverTimeoutError(f"Solver timed out after {_timeout_eff[0]:.1f} seconds.")

            # Un-scale to real space
            x_real = x_scaled * D_x

            try:
                # Evaluate residuals and jacobian
                res, jac = self._residuals_and_jacobian(x_real)

                # Track best iterate (in real space, unscaled residual)
                res_norm = np.linalg.norm(res)
                _improved = res_norm < best_res_norm
                if _improved:
                    best_res_norm = res_norm
                    best_x = x_real.copy()

                # Record convergence history point
                _t = time.perf_counter() - start_time
                if (_improved or _t - _last_record_t[0] >= _RECORD_INTERVAL) and len(
                    _history
                ) < 1000:
                    _history.append({"eval": _eval_count[0], "t": round(_t, 3), "norm": res_norm})
                    _last_record_t[0] = _t
                _eval_count[0] += 1

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
            except Exception as e:
                last_exception = e
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
                # Penalty that pulls Newton back toward the best physical point
                # seen so far.  F = x_scaled - x_best means Newton step = x_best
                # exactly (J=I), so hybr retreats to a known-safe iterate rather
                # than stepping deeper into unphysical territory.
                x_best_scaled = best_x * inv_D_x
                penalty = x_scaled - x_best_scaled
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
            final_x = best_x
            success = sol.success
            message = sol.message
            # Efficiently use tracked norm instead of re-evaluating
            final_norm = float(best_res_norm)
            # Guard against methods (e.g. lm) that report success at a
            # non-zero local minimum of ||F||^2.
            if success and final_norm > _RESIDUAL_TOL:
                success = False
                message = (
                    f"Solver reported success but |F|={final_norm:.3e} "
                    f"exceeds residual tolerance ({_RESIDUAL_TOL})."
                )
            if not success and last_exception:
                raise last_exception
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
                fallback_x = best_x
                fallback_norm = float(best_res_norm)
                if sol2.success and fallback_norm < _RESIDUAL_TOL:
                    final_x = fallback_x
                    success = True
                    message = f"Converged after fallback to hybr (original {method} failed)."
                    final_norm = fallback_norm
            except Exception:
                pass  # keep original failure

        # LM fallback for hybr oscillation on compressible networks.
        # Levenberg-Marquardt's (J^TJ + lI) regularisation prevents the
        # ill-conditioned Newton step that causes hybr to bounce near a
        # choked solution without making progress.  Start from the best
        # iterate hybr found; restore the full timeout so LM gets the
        # remaining 35% of the budget.
        if not success and _has_compressible:
            _elapsed = time.perf_counter() - start_time
            _remaining = (timeout - _elapsed) if timeout is not None else None
            if _remaining is None or _remaining > 2.0:
                _lm_started_at_eval[0] = _eval_count[0]
                _timeout_eff[0] = timeout
                # For MPCE networks the impulse + sum-mass residual structure
                # leaves hybr near-singular at low Mach; its "best iterate" is
                # often a worse LM starting point than the original x0. Restart
                # LM from x0 in that case to avoid getting trapped in hybr's
                # poor basin. For tee / compressible networks the existing
                # warm-start-from-hybr path is fine.
                _lm_x0_scaled = x0_scaled if _has_mpce else best_x * inv_D_x
                try:
                    root(
                        residuals_wrapper,
                        _lm_x0_scaled,
                        method="lm",
                        jac=use_jac,
                        options={"maxiter": _default_iters},
                    )
                    _lm_norm = float(best_res_norm)
                    if _lm_norm < _RESIDUAL_TOL:
                        final_x = best_x
                        success = True
                        message = (
                            f"Converged with LM fallback "
                            f"(hybr |F|={final_norm:.3e}, LM |F|={_lm_norm:.3e})."
                        )
                        final_norm = _lm_norm
                except SolverTimeoutError:
                    pass
                except Exception:
                    pass

        # Post-solve direction verification for constrained junctions.
        # strict=False soft-barrier mode can manufacture EXACT artifact
        # roots: the one-sided penalty alpha * mdot^2 on a slightly
        # reversed port can cancel that row's Pt-continuity mismatch, so
        # the solver converges (|F| ~ 1e-10) to a wrong-direction state
        # (exposed by constant-K warm starts on the certified audit,
        # 2026-07-04: a merge feed reversed to -0.0035 kg/s while the
        # strict residual raises at the same state). Demote such
        # solutions to honest failures instead of reporting them.
        if success:
            _sol_names = dict(zip(self.unknown_names, final_x, strict=False))
            _bad_junctions = [
                eid
                for eid, element in self.network.elements.items()
                if hasattr(element, "verify_solution_consistent")
                and not element.verify_solution_consistent(_sol_names)
            ]
            if _bad_junctions:
                success = False
                message = (
                    "Converged to a wrong-direction artifact root at "
                    f"junction(s) {_bad_junctions}: the soft-barrier "
                    "penalty balanced the residual with reversed or zero "
                    "port flow. Treating as non-converged; try "
                    "init_strategy='analytical_pt_prop' or a different "
                    "initial guess."
                )

        if not success:
            warnings.warn(
                f"NetworkSolver did not converge: {message} "
                f"(|F|={final_norm:.3e}). "
                "Returning best iterate.",
                stacklevel=2,
            )

        # Store solution for warm-start reuse
        self._last_x = final_x.copy()

        # Decompose final residual vector into per-unknown contributions.
        # One extra _residuals() call; same cost as a single Newton step.
        # Guard with try-except: mocked or patched _residuals_and_jacobian in
        # tests may raise, and diagnostic collection must never break the caller.
        try:
            _res_final = self._residuals(final_x)
            _breakdown = {
                name: float(abs(r)) for name, r in zip(self.unknown_names, _res_final, strict=False)
            }
            _worst = sorted(_breakdown.items(), key=lambda kv: kv[1], reverse=True)[:10]
        except Exception:
            _breakdown = {}
            _worst = []

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
            "convergence_history": _history,
            "residual_breakdown": _breakdown,
            "worst_residuals": [{"name": n, "residual": v} for n, v in _worst],
            "lm_started_at_eval": _lm_started_at_eval[0],
            "solver_settings_used": {
                "method": method,
                "init_strategy": init_strategy,
                "timeout": timeout,
                **options,
            },
        }

        sol_dict = dict(zip(self.unknown_names, final_x, strict=False))

        # Ensure _derived_states matches final_x (best iterate), not the
        # last trial point the solver happened to evaluate.  This is cheap:
        # forward T/Y propagation only, no residuals or Jacobian.
        self._propagate_states(final_x)
        from combaero import mass_to_mole

        for nid, (T, Y, _) in self._derived_states.items():
            sol_dict[f"{nid}.T"] = T
            for i, yi in enumerate(Y):
                sol_dict[f"{nid}.Y[{i}]"] = float(yi)

            # Store mole fractions as well (cheap, prevents redundant conversions)
            X = mass_to_mole(np.array(Y))
            for i, xi in enumerate(X):
                sol_dict[f"{nid}.X[{i}]"] = float(xi)

            # Node-specific diagnostics (e.g. mach, transport/thermo properties)
            node = self.network.nodes[nid]
            state = self._get_node_state(node, final_x)
            sol_dict[f"{nid}.mach"] = node.mach(state)
            if hasattr(node, "diagnostics"):
                diag = node.diagnostics(state)
                for key, val in diag.items():
                    sol_dict[f"{nid}.{key}"] = val

        # Element-specific diagnostics (e.g. throat Mach, P-ratio)
        from .components import (
            MultiPortChamberElement as _MultiPortChamberElement,
        )
        from .components import (
            TeeJunctionElement as _TeeJunctionElement,
        )

        sol_dict["__element_diag__"] = {}
        for eid, element in self.network.elements.items():
            m_indices = self._unknown_indices.get(eid, [])
            if isinstance(element, _TeeJunctionElement):
                nodes = self.network.nodes
                state_com = self._get_node_state(nodes[element.common_node], final_x)
                state_straight = self._get_node_state(nodes[element.straight_node], final_x)
                state_branch = self._get_node_state(nodes[element.branch_node], final_x)
                if m_indices:
                    state_com.m_dot = float(final_x[m_indices[0]])
                    state_branch.m_dot = float(final_x[m_indices[1]])
                    state_straight.m_dot = float(final_x[m_indices[0]]) - float(
                        final_x[m_indices[1]]
                    )
                diag = element.diagnostics(state_com, state_straight, state_branch)
            elif isinstance(element, _MultiPortChamberElement):
                nodes = self.network.nodes
                port_states = [
                    self._get_node_state(nodes[pid], final_x) for pid in element.port_nodes
                ]
                P_jct_val = float(final_x[m_indices[0]]) if m_indices else 0.0
                # Per-port mdots in junction convention (positive = out of
                # junction), same sign-mapping the residual call uses.
                port_mdots: list[float] = []
                for i, outer_id in enumerate(element._port_element_ids):
                    outer_indices = self._unknown_indices.get(outer_id, [])
                    outer_mdot = float(final_x[outer_indices[0]]) if outer_indices else 0.0
                    port_mdots.append(element._port_signs[i] * outer_mdot)
                diag = element.diagnostics(port_states, P_jct_val, port_mdots)
            else:
                state_in = self._get_node_state(self.network.nodes[element.from_node], final_x)
                state_out = self._get_node_state(self.network.nodes[element.to_node], final_x)
                # Inject solved m_dot so velocity-dependent diagnostics (e.g. Mach) are correct
                if m_indices:
                    m_solved = float(final_x[m_indices[0]])
                    state_in.m_dot = m_solved
                    state_out.m_dot = m_solved
                diag = element.diagnostics(state_in, state_out)
            sol_dict["__element_diag__"][eid] = diag
            for key, val in diag.items():
                sol_dict[f"{eid}.{key}"] = val

        # Wall-specific thermodynamics
        import combaero as cb

        for wid, wall in self.network.walls.items():
            obj_a = self.network.elements.get(wall.element_a) or self.network.nodes.get(
                wall.element_a
            )
            obj_b = self.network.elements.get(wall.element_b) or self.network.nodes.get(
                wall.element_b
            )

            node_feeding_a = (
                obj_a
                if wall.element_a in self.network.nodes
                else self.network.nodes[obj_a.from_node]
            )
            node_feeding_b = (
                obj_b
                if wall.element_b in self.network.nodes
                else self.network.nodes[obj_b.from_node]
            )

            state_a = self._get_node_state(node_feeding_a, final_x)
            state_b = self._get_node_state(node_feeding_b, final_x)

            if wall.element_a in self.network.elements:
                m_indices_a = self._unknown_indices.get(wall.element_a, [])
                if m_indices_a:
                    state_a.m_dot = final_x[m_indices_a[0]]
            if wall.element_b in self.network.elements:
                m_indices_b = self._unknown_indices.get(wall.element_b, [])
                if m_indices_b:
                    state_b.m_dot = final_x[m_indices_b[0]]

            ch_result_a = obj_a.htc_and_T(state_a) if obj_a.has_convective_surface else None
            ch_result_b = obj_b.htc_and_T(state_b) if obj_b.has_convective_surface else None

            if ch_result_a and ch_result_b:
                # Call multi-layer coupling logic (same as in residual evaluation)
                wall_res = wall.compute_coupling(
                    ch_result_a.h,
                    ch_result_a.T_aw,
                    obj_a.surface.area,
                    ch_result_b.h,
                    ch_result_b.T_aw,
                    obj_b.surface.area,
                )

                sol_dict[f"{wid}.Q"] = float(wall_res.Q)
                sol_dict[f"{wid}.T_hot"] = float(wall_res.T_hot)
                sol_dict[f"{wid}.h_a"] = float(ch_result_a.h)
                sol_dict[f"{wid}.h_b"] = float(ch_result_b.h)

                # Detailed temperature profile [K] (diagnostics)
                # Ensure profile uses correct side-A/side-B ordering
                t_over_k_layers = [L.r_val for L in wall.layers]
                profile, _q = cb.wall_temperature_profile(
                    ch_result_a.T_aw,
                    ch_result_b.T_aw,
                    ch_result_a.h,
                    ch_result_b.h,
                    t_over_k_layers,
                    wall.R_fouling,
                )
                sol_dict[f"{wid}.T_interface"] = [float(tp) for tp in profile]

        sol_dict["__complete_states__"] = self.extract_complete_states(sol_dict)
        sol_dict["__success__"] = success
        sol_dict["__message__"] = message
        sol_dict["__final_norm__"] = final_norm
        sol_dict["__unknown_names__"] = list(self.unknown_names)
        sol_dict["__x_solution__"] = list(final_x)
        # Convergence history: per-eval (eval_idx, t_elapsed_s, residual_norm)
        # triples accumulated by the residuals wrapper. Throttled at 0.1s
        # intervals (or whenever the iterate improves), capped at 1000 points.
        # Used by GUI diagnostics to plot the |F| vs iteration trace -- helps
        # diagnose slow-converging networks ("Newton overshoots 5 steps then
        # damps" vs "monotone but many iterations" tell different stories).
        sol_dict["__convergence_history__"] = list(_history)
        return sol_dict

    def extract_complete_states(self, result: dict[str, float] | None = None) -> dict[str, Any]:
        """
        Extract CompleteState for all network nodes.

        Leverages the existing C++ CompleteState structure for maximum accuracy
        and performance. By default, it uses the solver's internal cached derived
        states populated during solve().

        Args:
            result: Optional dictionary containing node T, P, X and Y. If omitted,
                uses internal solver state from the last solve() call.

        Returns:
            Dictionary mapping node IDs to CompleteState objects.
        """
        import combaero as cb
        from combaero import mass_to_mole

        complete_states = {}

        if result is None:
            # OPTIMIZED: Use cached derived states from _propagate_states
            # Note: This is primarily for manual diagnostic exploration after solve()
            if not self._derived_states:
                return {}

            for nid, (T, Y, _) in self._derived_states.items():
                # For P, we rely on node.P or total P if defined in the network
                node = self.network.nodes[nid]
                P = getattr(node, "P", getattr(node, "Pt", 101325.0))
                X = mass_to_mole(np.array(Y))
                complete_states[nid] = cb.complete_state(T, P, X)
            return complete_states

        # DICT PATH: Consumption from solve() result (as used in the GUI)
        # Try to get composition from first node that has data to use as default
        default_X = None
        for node_id in self.network.nodes:
            # Prefer pre-computed X if available
            if f"{node_id}.X[0]" in result:
                X_list = []
                idx = 0
                while f"{node_id}.X[{idx}]" in result:
                    X_list.append(result[f"{node_id}.X[{idx}]"])
                    idx += 1
                if X_list:
                    default_X = np.array(X_list)
                    break
            elif f"{node_id}.Y[0]" in result:
                Y_list = []
                idx = 0
                while f"{node_id}.Y[{idx}]" in result:
                    Y_list.append(result[f"{node_id}.Y[{idx}]"])
                    idx += 1
                if Y_list:
                    default_X = mass_to_mole(np.array(Y_list))
                    break

        if default_X is None:
            default_X = cb.humid_air_composition(350.0, 101325.0, 0.0)

        for node_id in self.network.nodes:
            try:
                T = result.get(f"{node_id}.T", result.get(f"{node_id}.Tt", 300.0))
                P = result.get(f"{node_id}.P", result.get(f"{node_id}.Pt", 101325.0))

                # Use pre-computed X if available, otherwise fallback to Y-parsing
                X = None
                if f"{node_id}.X[0]" in result:
                    X_list = []
                    idx = 0
                    while f"{node_id}.X[{idx}]" in result:
                        X_list.append(result[f"{node_id}.X[{idx}]"])
                        idx += 1
                    if X_list:
                        X = np.array(X_list)

                if X is None and f"{node_id}.Y[0]" in result:
                    Y_list = []
                    idx = 0
                    while f"{node_id}.Y[{idx}]" in result:
                        Y_list.append(result[f"{node_id}.Y[{idx}]"])
                        idx += 1
                    if Y_list:
                        X = mass_to_mole(np.array(Y_list))

                X = X if X is not None else default_X

                if not (math.isnan(T) or math.isnan(P)):
                    complete_states[node_id] = cb.complete_state(T, P, X)
                else:
                    complete_states[node_id] = None
            except Exception:
                complete_states[node_id] = None

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
