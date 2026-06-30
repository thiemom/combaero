"""
MPCE-v2: Mynard Unified0D residual structure on top of MPCE-v1's topology.

Subclasses :class:`MultiPortChamberElement` so all the topology resolution
(port-MCN flagging, port area inheritance, sign convention validation) is
inherited unchanged. Only the residuals method differs: instead of the
impulse-CV formula with the empirical cross-coupling correction
(``2*sin(theta)*cos((3/4)theta)*rho*u_0*u_i``), we use Mynard's
pseudosupplier-based smooth K relation.

Residual structure:

    For each port i:
        R_i = Pt_port_i - Pt_jct + K_i * q_dyn_pseudo = 0
    Plus the sum-mass conservation residual:
        R_mass = sum(port_mdots) = 0

where:
  Pt_jct          = junction chamber STAGNATION pressure (the lone unknown,
                    re-purposed from MPCE-v1's static P_jct slot)
  Pt_port_i       = port-MCN stagnation pressure (from the MCN's converged
                    state propagation)
  K_i             = Mynard per-port loss coefficient. For supplier ports
                    (flow into junction) K_i = 0 (total-pressure continuity
                    per Mynard Eq 14). For collector ports K_i comes from
                    the pseudosupplier-reorientation C value.
  q_dyn_pseudo    = 0.5 * rho_ref * U_pseudo^2 -- reference dynamic head at
                    the pseudosupplier velocity computed by Mynard.

This is a *prototype*. The Jacobian is left as an empty dict so the solver
falls back to numerical differencing. Once the residual structure is
validated empirically (does it converge under all three BC topologies?
Does it match Mynard's reference K at simple cases? Does theta=120 no
longer fail?), we move to either analytical Jacobian or a C++ port.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

from combaero.network._mpce_v2_jacobian import dKQ_dmdot_separating_T
from combaero.network.components import (
    MultiPortChamberElement,
    NetworkMixtureState,
)
from validation.junction.models.mynard2010 import junction_loss_coefficient

FlowDirection = Literal["merge", "branch"]


class MPCEv2Element(MultiPortChamberElement):
    """Mynard Unified0D residual on MPCE-v1's topology framework.

    The ``jacobian_method`` attribute controls how the K*q_dyn term's
    Jacobian is computed:
      - ``"sympy"`` (default): analytical sympy-derived Jacobian for the
        canonical 3-port separating T (1 supplier on port 0, straight at
        0 deg, branch at theta_branch). Uses 1 Mynard evaluation per
        residual call. Non-canonical topologies (joining, theta != 0 on
        the straight arm, N != 3) auto-fall-back to FD because the
        sympy derivation is specific to the canonical case. Validated on
        the separating K6 / K5 audit: identical convergence + accuracy
        vs FD with ~1.14x wall-time speedup.
      - ``"fd"``: numerical finite-difference on Mynard. Costs N+1
        Mynard evaluations per residual call. Always available across
        all topologies. The previous default; flipped to sympy in PR
        adding GUI solver tweaks once the empirical "fd is more robust
        on mfb_two_pb low-q" concern was found to no longer apply at
        the post-soft-barrier audit point.

    ``flow_direction`` constrains which physical flow regime this element
    represents:
      - ``"merge"``: 2 suppliers + 1 collector (joining flow). The
        residual asserts at solve time and raises if the observed mdot
        signs imply separating flow instead.
      - ``"branch"``: 1 supplier + 2 collectors (separating flow). Errors
        on observed joining flow.
    """

    jacobian_method: str = "sympy"

    # Soft-barrier penalty scale used when ``strict=False`` and the observed
    # flow direction disagrees with the declared one. Wraps the
    # one-sided quadratic ``alpha * max(0, -expected_sign * mdot)^2``.
    # Calibrated so that a wrong-sign mdot of 0.1 kg/s contributes roughly
    # 1e5 Pa to the residual -- the natural Pt scale. Tunable via attribute.
    soft_penalty_alpha: float = 1.0e7

    # Joining-side etransfer correction (combaero extension to faithful
    # Mynard 2010). Default 0.2 calibrated against Bassett K11_corr/K12_corr
    # + Idelchik 1966 tabulated values at psi in [1.25, 3.33], theta in
    # {30, 45, 90} (analytical-only anchors; measured points held out as
    # independent validation). Independent validation on Bassett-measured
    # K11/K12 (145 points across 3 topologies) shows -10% mean error and
    # -14% max error on the imposed_q topology (the cleanest test);
    # essentially tied on three_pb / mfb_two_pb. Reduces to identity at
    # psi=1 by construction so the equal-area baseline is unchanged. See
    # the calibration write-up in the iteration-8 commit and the
    # validation script in scripts/calibrate_mpce_v2_etransfer.py.
    DEFAULT_JOINING_ETRANSFER_ALPHA: float = 0.2

    def __init__(
        self,
        id: str,
        inlet_nodes: list[str],
        outlet_nodes: list[str],
        inlet_angles_deg: list[float] | None = None,
        outlet_angles_deg: list[float] | None = None,
        port_areas: list[float] | None = None,
        flow_direction: FlowDirection = "branch",
        strict: bool = True,
        joining_etransfer_alpha: float | None = None,
    ):
        super().__init__(
            id=id,
            inlet_nodes=inlet_nodes,
            outlet_nodes=outlet_nodes,
            inlet_angles_deg=inlet_angles_deg,
            outlet_angles_deg=outlet_angles_deg,
            port_areas=port_areas,
        )
        if flow_direction not in ("merge", "branch"):
            raise ValueError(
                f"MPCEv2Element '{id}': flow_direction must be 'merge' or "
                f"'branch', got {flow_direction!r}."
            )
        self.flow_direction: FlowDirection = flow_direction
        # strict=True: raise on direction mismatch (default, safe).
        # strict=False: return a soft-barrier residual that pulls Newton
        # back toward the correct sign. Use for hard cases where the
        # solver starts in the wrong basin but a physical root exists.
        self.strict: bool = strict
        # joining_etransfer_alpha: pass None to use the calibrated default
        # (DEFAULT_JOINING_ETRANSFER_ALPHA = 0.2). Pass 0.0 to disable the
        # correction entirely (faithful-port Mynard). Pass a custom value
        # if you've re-calibrated against a different anchor set.
        self.joining_etransfer_alpha: float = (
            self.DEFAULT_JOINING_ETRANSFER_ALPHA
            if joining_etransfer_alpha is None
            else float(joining_etransfer_alpha)
        )

    def verify_solution_consistent(
        self,
        sol: dict[str, float],
        eps: float = 1e-6,
    ) -> bool:
        """Post-solve check that converged mdots match the declared direction.

        For each port-connecting element, the canonical-direction unknown
        ``{elem_id}.m_dot`` should be > 0 at convergence (the sign mapping
        ``port_mdots[i] = port_signs[i] * outer_mdot[i]`` guarantees that
        outer_mdot > 0 implies the correct in/out direction at the port).
        Returns False when soft mode landed at a wrong-basin fixed point
        with at least one outer mdot at or below ``eps``.
        """
        for elem_id in self._port_element_ids:
            if not elem_id:
                continue
            key = f"{elem_id}.m_dot"
            if key not in sol:
                continue  # connecting element has no m_dot unknown
            if sol[key] < eps:
                return False
        return True

    def diagnostics(  # type: ignore[override]
        self,
        states: list[NetworkMixtureState],
        Pt_jct: float,
        port_mdots: list[float] | None = None,
    ) -> dict[str, float]:
        """Emit MPCE-v2 diagnostics: parent fields + Mynard K + named aliases.

        Beyond the parent class output (``P_jct``, ``n_ports``, per-port
        ``P/T/area/sign``), this re-evaluates Mynard at the converged
        operating point to expose the K coefficient used in the residual,
        plus the topology-aware aliases that the GUI displays:

          - separating (flow_direction='branch'): ``K_straight``,
            ``K_branch``, ``mass_flow_ratio = m_dot_branch / m_dot_com``
          - joining (flow_direction='merge'): ``K11``, ``K12``,
            ``mass_flow_ratio = m_dot_branch / m_dot_com`` (same convention
            as Bassett's joining-T q)

        Re-evaluation rather than caching keeps diagnostics consistent with
        a converged state and avoids the residual call having to stash
        per-element state.
        """
        diag = super().diagnostics(states, Pt_jct, port_mdots)
        if port_mdots is None or len(port_mdots) != self.N:
            return diag

        # Per-port mdot in junction convention is already in `port_mdots`.
        # Convert to Mynard's U (positive = supplier into junction).
        rho_port = np.array([float(s.density()) for s in states])
        A = np.array([float(a or 0.0) for a in self.port_areas])
        if (A <= 0.0).any():
            return diag  # underspecified geometry; skip K emission
        U_mynard = -np.array(port_mdots) / (rho_port * A)

        # Reproduce the residual's angle-handling: axial-back at pi for the
        # lone opposite-sign port.
        theta_rad = np.array([math.radians(float(t)) for t in self.port_angles_deg])
        sup_count = int(np.sum(U_mynard > 1e-9))
        col_count = int(np.sum(U_mynard < -1e-9))
        if sup_count == 1 and col_count == self.N - 1:
            theta_rad[int(np.argmax(U_mynard))] = math.pi
        elif col_count == 1 and sup_count == self.N - 1:
            theta_rad[int(np.argmin(U_mynard))] = math.pi
        else:
            return diag  # degenerate / multi-supplier-multi-collector

        try:
            mynard = junction_loss_coefficient(
                U_mynard,
                A,
                theta_rad,
                joining_etransfer_alpha=self.joining_etransfer_alpha,
            )
        except Exception:
            return diag

        if mynard.K is None:
            return diag

        # Mynard K is per non-common port (collectors for separating, suppliers
        # for joining). non_common_idxs preserves Si/Ci order, which for the
        # 3-port canonical setup matches (straight, branch) for separating and
        # (straight, branch) for joining (the lone-opposite port is common).
        sup_mask = U_mynard > 1e-9
        col_mask = U_mynard < -1e-9
        non_common_idxs = np.where(col_mask)[0] if sup_count == 1 else np.where(sup_mask)[0]
        if len(non_common_idxs) != len(mynard.K):
            return diag

        K_per_port = np.zeros(self.N)
        for j, port_idx in enumerate(non_common_idxs):
            K_per_port[port_idx] = float(mynard.K[j])
            diag[f"port_{port_idx}_K"] = float(mynard.K[j])

        # Topology-aware aliases. port ordering:
        #   "branch" (separating): port 0 = com (inlet), 1 = straight, 2 = branch
        #   "merge"  (joining):    port 0 = straight, 1 = branch, 2 = com (outlet)
        if self.N == 3:
            if self.flow_direction == "branch":
                diag["K_straight"] = float(K_per_port[1])
                diag["K_branch"] = float(K_per_port[2])
                m_com = abs(float(port_mdots[0]))
                m_branch = abs(float(port_mdots[2]))
            else:  # "merge"
                # K11 (straight->common), K12 (branch->common)
                diag["K11"] = float(K_per_port[0])
                diag["K12"] = float(K_per_port[1])
                m_com = abs(float(port_mdots[2]))
                m_branch = abs(float(port_mdots[1]))
            if m_com > 1e-12:
                diag["mass_flow_ratio"] = m_branch / m_com

        # Convenience: Pt at the chamber (equals P_jct in MPCE-v2 naming, since
        # `P_jct` re-purposes the MPCE-v1 static slot for Pt).
        diag["Pt_jct"] = float(Pt_jct)
        return diag

    def _soft_barrier_residual(
        self,
        states: list[NetworkMixtureState],
        Pt_jct: float,
        port_mdots: list[float],
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """Continuity + sign-coercing penalty when in the wrong basin.

        Returns ``R_i = (Pt_i - Pt_jct) + alpha * max(0, -e_i * mdot_i)^2``
        where ``e_i = self._port_signs[i]`` is the canonical sign for that
        port given the declared ``flow_direction``. The penalty is C^1
        smooth, zero in the correct-direction half-space, and has gradient
        ``2 * alpha * (-e_i) * max(0, -e_i * mdot_i)`` in the wrong half.
        Newton's step on this residual is pulled toward ``mdot_i = 0``,
        from which a sign flip restores the strict-physics residual on the
        next iteration.
        """
        N = self.N
        alpha = float(self.soft_penalty_alpha)
        residuals: list[float] = []
        jac: dict[int, dict[str, float]] = {}
        for i in range(N):
            e_i = float(self._port_signs[i])
            mdot_i = float(port_mdots[i])
            slack = max(0.0, -e_i * mdot_i)  # positive when wrong-direction
            pen = alpha * slack * slack
            Pt_i = float(states[i].Pt)
            residuals.append(Pt_i - Pt_jct + pen)

            row: dict[str, float] = {
                f"{self.port_nodes[i]}.Pt": 1.0,
                f"{self.id}.P_jct": -1.0,
            }
            if slack > 0.0:
                # d(pen)/d(mdot_i) = 2 * alpha * slack * (-e_i)
                # Chain through outer-mdot sign: d(mdot_i)/d(outer_i) = sign_i = e_i.
                # So d(pen)/d(outer_i) = 2 * alpha * slack * (-e_i) * e_i = -2*alpha*slack.
                mdot_var = f"{self._port_element_ids[i]}.m_dot"
                row[mdot_var] = row.get(mdot_var, 0.0) - 2.0 * alpha * slack
            jac[i] = row

        residuals.append(sum(port_mdots))
        mass_row: dict[str, float] = {}
        for i in range(N):
            mass_var = f"{self._port_element_ids[i]}.m_dot"
            mass_row[mass_var] = mass_row.get(mass_var, 0.0) + self._port_signs[i]
        jac[N] = mass_row
        return residuals, jac

    def residuals(  # type: ignore[override]
        self,
        states: list[NetworkMixtureState],
        Pt_jct: float,
        port_mdots: list[float],
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        N = self.N
        A = np.array([float(a) for a in self.port_areas])

        # Convert: port_mdots are in junction convention (positive = out of
        # junction). For Mynard: positive U = supplier (into junction).
        # So Mynard's U = -port_mdots / (rho * A).
        rho_port = np.array([float(s.density()) for s in states])
        U_mynard = -np.array(port_mdots) / (rho_port * A)

        # Convert MPCE-v1 branch angles (angle from main axis) to Mynard
        # vessel-direction convention. The "axial back" port -- the single
        # port whose flow direction opposes the other N-1 ports (the lone
        # supplier in separating flow, or lone collector in joining flow) --
        # gets angle pi (pointing back along the main duct axis). All other
        # ports keep their MPCE-v1 angles. This is flow-direction-aware so
        # the same element handles separating and joining cleanly.
        theta_rad = np.array([math.radians(float(t)) for t in self.port_angles_deg])
        sup_count = int(np.sum(U_mynard > 1e-9))
        col_count = int(np.sum(U_mynard < -1e-9))
        if sup_count == 1 and col_count == N - 1:
            axial_idx = int(np.argmax(U_mynard))
            theta_rad[axial_idx] = math.pi
        elif col_count == 1 and sup_count == N - 1:
            axial_idx = int(np.argmin(U_mynard))
            theta_rad[axial_idx] = math.pi

        # Mynard requires at least one supplier (U > 0) AND at least one
        # collector (U < 0). When the imposed BCs imply a degenerate state
        # (all U near zero, or all same sign), fall back to a continuity-only
        # residual: Pt_i = Pt_jct, mass conservation.
        sup_mask = U_mynard > 1e-9
        col_mask = U_mynard < -1e-9
        if not sup_mask.any() or not col_mask.any():
            residuals = [float(s.Pt) - Pt_jct for s in states]
            residuals.append(sum(port_mdots))
            jac: dict[int, dict[str, float]] = {}
            return residuals, jac

        # Constrained topology: refuse if observed direction does not match
        # the declared one. Per-port check: port_mdots[i] should have the
        # canonical sign self._port_signs[i] (e.g., port_signs[i]=-1 for an
        # inlet -> port_mdots[i] < 0). The earlier count-based check
        # (1 supplier vs N-1) was insufficient because it accepted
        # configurations with the right count but wrong port distribution
        # (e.g., for "merge" port_mdots=[-, +, -] has 2 suppliers but at
        # the wrong ports). Per-port catches this exactly.
        # Correct direction at port i: sign(port_mdots[i]) == sign(port_signs[i]),
        # i.e., port_signs[i] * port_mdots[i] > 0. Wrong when the product is
        # strictly negative (a tolerance keeps the boundary mdot=0 quiet).
        wrong_ports = [i for i in range(N) if self._port_signs[i] * port_mdots[i] < -1e-9]
        if wrong_ports:
            if self.strict:
                expected = [
                    "into junction" if self._port_signs[i] < 0 else "out of junction"
                    for i in wrong_ports
                ]
                raise ValueError(
                    f"MPCEv2Element '{self.id}': declared "
                    f"flow_direction={self.flow_direction!r} but observed "
                    f"wrong flow direction at port(s) {wrong_ports} "
                    f"(expected {expected})."
                )
            # Soft-barrier fallback. Replace the physics residual with
            # continuity (Pt_i = Pt_jct) plus a one-sided quadratic penalty
            # on each port mdot whose sign is wrong, scaled by
            # ``soft_penalty_alpha``. Pulls Newton back toward the boundary
            # mdot=0 from which a sign flip can restore the declared regime.
            return self._soft_barrier_residual(states, Pt_jct, port_mdots)

        try:
            mynard = junction_loss_coefficient(
                U_mynard,
                A,
                theta_rad,
                joining_etransfer_alpha=self.joining_etransfer_alpha,
            )
        except Exception:
            # Numerical failure (e.g. singular pseudosupplier). Same fallback.
            residuals = [float(s.Pt) - Pt_jct for s in states]
            residuals.append(sum(port_mdots))
            return residuals, {}

        # ITERATION-2: use mynard.K (Matlab line 73) with common-side q_dyn.
        # This is the physically correct normalization: K is defined such
        # that Pt_common - Pt_other = K * 0.5 * rho_com * u_com^2.
        # The common branch is the single supplier (diverging) or single
        # collector (converging); mynard.K is per non-common port.
        #
        # Sign of the K*q_dyn term in the residual:
        #   Separating: common is supplier (higher Pt), collectors have
        #               LOWER Pt -> Pt_col = Pt_jct - K*q_dyn -> +1 in residual
        #   Joining:    common is collector (lower Pt), suppliers have
        #               HIGHER Pt -> Pt_sup = Pt_jct + K*q_dyn -> -1 in residual
        if int(np.sum(sup_mask)) == 1:
            common_mask = sup_mask
            non_common_idxs = np.where(col_mask)[0]
            K_term_sign = +1.0
        else:
            common_mask = col_mask
            non_common_idxs = np.where(sup_mask)[0]
            K_term_sign = -1.0

        K_per_port = np.zeros(N)
        if mynard.K is not None and len(mynard.K) == len(non_common_idxs):
            for j, port_idx in enumerate(non_common_idxs):
                K_per_port[port_idx] = float(mynard.K[j])

        # Common-side dynamic head: q_dyn at the single common port.
        common_idx = int(np.where(common_mask)[0][0])
        u_com = abs(float(U_mynard[common_idx]))
        rho_com = float(rho_port[common_idx])
        q_dyn_com = 0.5 * rho_com * u_com * u_com

        # Residual: Pt_i - Pt_jct + K_i * q_dyn_com = 0
        # Common port: K_i = 0, so Pt_common = Pt_jct (continuity).
        # Other ports: K * q_dyn_com is the loss term per Mynard's
        # stagnation-pressure relation (Eq 15).
        residuals: list[float] = []
        for i in range(N):
            Pt_i = float(states[i].Pt)
            R_i = Pt_i - Pt_jct + K_term_sign * K_per_port[i] * q_dyn_com
            residuals.append(R_i)
        residuals.append(sum(port_mdots))

        # Jacobian: linear pieces explicit + analytical (sympy-derived) for
        # the K*q_dyn term in the canonical 3-port separating case; FD
        # fallback for everything else.
        KQ_base = K_per_port * q_dyn_com  # per-port loss term (collectors nonzero)
        dKQ_dmdot = np.zeros((N, N))

        is_canonical_separating_T = (
            N == 3
            and int(np.sum(sup_mask)) == 1
            and bool(sup_mask[0])  # port 0 is the supplier
            and abs(self.port_angles_deg[1]) < 1e-9  # straight at 0
        )
        if is_canonical_separating_T and self.jacobian_method == "sympy":
            dKQ_dmdot = dKQ_dmdot_separating_T(
                np.asarray(port_mdots, dtype=float),
                rho_port,
                A,
                math.radians(float(self.port_angles_deg[2])),
            )
        else:
            # FD fallback: N+1 Mynard calls.
            eps_scale = 1e-4
            for j in range(N):
                mdot_eps = max(abs(port_mdots[j]) * eps_scale, 1e-7)
                mdot_pert = list(port_mdots)
                mdot_pert[j] = port_mdots[j] + mdot_eps
                U_pert = -np.array(mdot_pert) / (rho_port * A)
                if (np.sign(U_pert) != np.sign(U_mynard)).any():
                    continue
                try:
                    m_pert = junction_loss_coefficient(
                        U_pert,
                        A,
                        theta_rad,
                        joining_etransfer_alpha=self.joining_etransfer_alpha,
                    )
                except Exception:
                    continue
                if m_pert.K is None or len(m_pert.K) != len(non_common_idxs):
                    continue
                K_pert = np.zeros(N)
                for k, port_idx in enumerate(non_common_idxs):
                    K_pert[port_idx] = float(m_pert.K[k])
                u_com_pert = abs(float(U_pert[common_idx]))
                q_dyn_pert = 0.5 * rho_com * u_com_pert * u_com_pert
                KQ_pert = K_pert * q_dyn_pert
                dKQ_dmdot[:, j] = (KQ_pert - KQ_base) / mdot_eps

        jac: dict[int, dict[str, float]] = {}
        for i in range(N):
            row: dict[str, float] = {
                f"{self.port_nodes[i]}.Pt": 1.0,
                f"{self.id}.P_jct": -1.0,
            }
            for j in range(N):
                if abs(dKQ_dmdot[i, j]) > 0.0:
                    mdot_var = f"{self._port_element_ids[j]}.m_dot"
                    # Outer mdot has been sign-mapped: port_mdots[j] = sign_j * outer.
                    # Chain rule: d(KQ_i)/d(outer_j) = dKQ_dmdot[i,j] * sign_j.
                    # K_term_sign carries the joining/separating residual flip.
                    row[mdot_var] = row.get(mdot_var, 0.0) + (
                        K_term_sign * dKQ_dmdot[i, j] * self._port_signs[j]
                    )
            jac[i] = row

        mass_row: dict[str, float] = {}
        for i in range(N):
            mass_var = f"{self._port_element_ids[i]}.m_dot"
            mass_row[mass_var] = mass_row.get(mass_var, 0.0) + self._port_signs[i]
        jac[N] = mass_row

        return residuals, jac
