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
from combaero.network._mynard2010 import junction_loss_coefficient
from combaero.network.components import (
    MultiPortChamberElement,
    NetworkMixtureState,
)

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

    #: TUNED CONSTANT -- combaero's joining-side etransfer correction, an
    #: extension to Mynard 2015 (not in the paper). Vanishes at psi = 1 by
    #: construction.
    #:
    #: Source of the value: ``validation/junction/calibrate_etransfer.py``,
    #: anchored on Bassett K11_corr/K12_corr + Idelchik tables over
    #: psi in [1.25, 3.33], theta in {30, 45, 90}, q in [0.1, 0.9].
    #:
    #: Provenance caveats (issue #271):
    #:   - 0.2 was fitted with the Bassett K11 anchor MIRRORED (lateral
    #:     fraction passed where Table 1 wants the straight one) and on
    #:     pre-#212 plumbing. The corrected-axis optimum is 0.28-0.31 and the
    #:     fit is shallow (#283). The value here has NOT been moved yet.
    #:   - The "145-point independent validation" once cited here is not in
    #:     the script that produced the value (scripts/
    #:     calibrate_mpce_v2_etransfer.py from #197: mirrored axis, no
    #:     validation code, stale import path, absolute DATA_DIR -- retired
    #:     and superseded by validation/junction/calibrate_etransfer.py). If
    #:     it was measured, it was outside the script and is not
    #:     reproducible; the in-network scorecard (#282) is the first
    #:     validation this constant has had that can be re-run.
    #:   - It is measured together with Mynard's eta (both live in
    #:     ``etransfer``; see ``eta_scale``), never one half of the pair.
    #:
    #: PROVEN 2026-09-05 (#271 step 1.3, imposed_q, eta on, alpha swept):
    #:   alpha        0.0     0.1     0.15    0.20    0.25    0.30
    #:   Idelchik all 0.512   0.404   0.380   0.377   0.385   0.411   (MAE)
    #:   Idelchik K11 bias    -0.309  ...    -0.025  +0.046  +0.116
    #: In-network optimum is 0.2, shallow over [0.15, 0.25]; sources pull
    #: slightly apart (Bassett K11 prefers ~0.25, Bassett K12 and Idelchik
    #: K11 ~0.15). The corrected-axis ANCHOR refit (#283) said 0.28-0.31 and
    #: is worse here (Idelchik +0.034): the anchors are Bassett's analytical
    #: curves, this table is measured/tabulated data, and measured wins.
    #: No alpha fixes Idelchik K12 (bias -0.27 at the MAE minimum) -- that
    #: is the closure's K12@90 shape limit, not a tuning question.
    #: No effect on any dividing cell or any psi = 1 cell, verified.
    #: Switch: pass 0.0.
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
        eta_scale: float = 1.0,
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
        # Measurement switch on Mynard's eta (MYNARD_ETA_A0/A1); 1.0 = faithful
        # port. Exposed so the term can be scored on and off (#271 step 1.3).
        self.eta_scale: float = float(eta_scale)
        # Mynard's K conversion (Eq 18) is defined for three branches only
        # (the reference code computes K for n <= 3). A larger junction used
        # to return finite residuals with an entirely zero loss-term Jacobian
        # -- silently (issue #271). Refuse at construction; the native C-form
        # residual is the route to N > 3 (design doc sec 8 step 3b).
        if self.N > 3:
            raise ValueError(
                f"MPCEv2Element '{id}': {self.N} ports, but the Mynard K closure is "
                f"defined for 3-branch junctions only. Split the junction, or use "
                f"MultiPortChamberElement."
            )

    def verify_solution_consistent(
        self,
        sol: dict[str, float],
        eps: float = 1e-6,
        energy_rel_tol: float = 1e-9,
    ) -> bool:
        """Post-solve check: declared flow directions, and an energy balance.

        **Direction.** For each port-connecting element, the canonical-direction
        unknown ``{elem_id}.m_dot`` should be > 0 at convergence (the sign
        mapping ``port_mdots[i] = port_signs[i] * outer_mdot[i]`` guarantees
        that outer_mdot > 0 implies the correct in/out direction at the port).
        Returns False when soft mode landed at a wrong-basin fixed point with
        at least one outer mdot at or below ``eps``.

        **Energy.** A passive junction cannot create flow work, so the
        mass-weighted total pressure must not increase across it:

            dissipation = sum_in |m_i| Pt_i - sum_out |m_i| Pt_i >= 0

        v1 has an energy check (`MultiPortChamberElement`); v2 had none, which
        is the gap this closes (issue #271, defect 10). The v1 form bounds each
        collector's Pt by the single supplier's, which is deliberately weak and
        single-supplier only. The mass-weighted balance is both stricter and
        more permissive in the right places: it **allows** one branch to gain
        total pressure at another's expense -- real physics, and the reason
        Bassett's K5 and Hager's xi_t go negative in dividing flow -- while
        forbidding a net gain, and it applies to joining flow too.

        Equivalently, in coefficient terms the condition is that the
        mass-weighted mean K is non-negative. Audited before being imposed:
        Bassett's analytical pair satisfies it everywhere and tends to exactly
        zero as the lateral flow vanishes, which is what a junction that
        diverts nothing must dissipate. The MODEL violates it at low lateral
        fraction (down to -0.06 near q = 0.1), which is the K_straight defect
        of #272 showing up as a hard physical violation rather than an accuracy
        gap. The mechanism is Mynard's energy-transfer factor (Eq 35-36): it is
        written as an exchange between collectors, but nothing constrains the
        credit to equal the debit, and measured on the closure it lowers the
        flow-weighted sum by up to 0.27.

        ``energy_rel_tol`` is relative to ``sum |m_i| * max |Pt_i|``. A
        converged solve leaves a residual of order 1e-10 Pa against dynamic
        heads of order 10 Pa, so the numerical floor on this quantity is
        ~1e-12 relative; 1e-9 sits three decades above it and still catches
        every measured violation, the smallest of which is ~1e-3 in
        weighted-mean-K units.

        SCOPE. This is the model's own dissipation identity, not a general
        entropy statement: it says the flow-weighted sum of the closure's loss
        terms may not be net negative. Where the port total temperatures differ
        the true entropy production also carries a mixing term, which is
        positive and is not counted here -- so on a non-isothermal junction the
        check is conservative in the wrong direction and could in principle
        reject a solution whose full entropy budget is fine. The junction is
        documented for low Mach and the validation fixtures are isothermal, so
        that case does not arise today; revisit it if the element is ever used
        across a large temperature difference.

        Missing keys mean True: incomplete solution dicts are not policed, the
        same convention v1 uses.
        """
        for elem_id in self._port_element_ids:
            if not elem_id:
                continue
            key = f"{elem_id}.m_dot"
            if key not in sol:
                continue  # connecting element has no m_dot unknown
            if sol[key] < eps:
                return False

        flows: list[float] = []
        pts: list[float] = []
        for i in range(self.N):
            elem_id = self._port_element_ids[i]
            m_key = f"{elem_id}.m_dot"
            pt_key = f"{self.port_nodes[i]}.Pt"
            if not elem_id or m_key not in sol or pt_key not in sol:
                return True
            # Junction convention: positive = flow OUT of the junction.
            flows.append(float(self._port_signs[i]) * float(sol[m_key]))
            pts.append(float(sol[pt_key]))

        net_outflow_of_energy = sum(f * pt for f, pt in zip(flows, pts, strict=True))
        scale = sum(abs(f) for f in flows) * max(abs(pt) for pt in pts)
        return net_outflow_of_energy <= energy_rel_tol * scale

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
                eta_scale=self.eta_scale,
            )
        except (IndexError, ValueError) as exc:
            # The two exceptions the closure can raise on a degenerate flow
            # split. Post-solve, so report rather than raise -- but never
            # silently: an absent K field must be attributable.
            diag["closure_error"] = type(exc).__name__
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
            diag["m_dot_com"] = m_com
            diag["m_dot_branch"] = m_branch
            diag["m_dot_straight"] = max(m_com - m_branch, 0.0)

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
        if not sup_mask.any() and not col_mask.any():
            # Every port ~zero: a cold-start iterate. The soft barrier at zero
            # slack IS the continuity residual -- with its Jacobian. This used
            # to return jac = {}: N+1 zero rows, and every solve that reached
            # it stalled (26 of 26 on the scorecard, issue #271).
            return self._soft_barrier_residual(states, Pt_jct, port_mdots)
        # Every port flowing the same way is not handled here: it violates
        # mass conservation and is a wrong-direction state for at least one
        # port, which the check below routes to the strict raise or the soft
        # barrier -- both of which carry a Jacobian that pulls the offending
        # port back toward zero. Returning a lossless residual without one
        # was where those solves went to die.

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

        # A port excluded from both masks (|U| <= 1e-9, typically an exact
        # zero in an initial guess) must not be left for the closure to
        # classify on its own: Mynard puts Q = -0.0 with the suppliers, the two
        # classifications disagree, K comes back mis-shaped and was silently
        # zeroed -- a lossless junction 384 times per scorecard (issue #271).
        # Snap it to a tiny flow in its DECLARED direction and rebuild the
        # masks from the snapped values, so both sides agree by construction
        # and the FlowRatio -> 0 limit (K -> 1 for a dead outlet, step 1.4)
        # takes over.
        excluded = ~(sup_mask | col_mask)
        if excluded.any():
            U_mynard = U_mynard.copy()
            for i in np.where(excluded)[0]:
                # _port_signs: -1 = inlet (Mynard supplier, U > 0), +1 = outlet.
                U_mynard[i] = -float(self._port_signs[i]) * 1e-9
            sup_mask = U_mynard > 0.0
            col_mask = U_mynard < 0.0

        try:
            mynard = junction_loss_coefficient(
                U_mynard,
                A,
                theta_rad,
                joining_etransfer_alpha=self.joining_etransfer_alpha,
                eta_scale=self.eta_scale,
            )
        except (IndexError, ValueError) as exc:
            # The closure can raise only on a degenerate flow split (empty
            # supplier or collector mask), and the guard above already routes
            # those to the continuity fallback before this call. Reaching
            # here therefore means the two classifications disagree -- an
            # anomaly worth seeing, not one to paper over. This used to
            # return a lossless residual with an EMPTY Jacobian for *any*
            # exception, which turned a plumbing error into a plausible-
            # looking lossless junction mid-solve (issue #271). Programming
            # errors now propagate untouched.
            raise RuntimeError(
                f"MPCEv2Element '{self.id}': Mynard closure failed on a "
                f"non-degenerate state (port_mdots={list(port_mdots)}, "
                f"suppliers={sup_mask.tolist()}, collectors={col_mask.tolist()}). "
                f"The degenerate-split guard should have caught this first."
            ) from exc

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
                        eta_scale=self.eta_scale,
                    )
                except (IndexError, ValueError) as exc:
                    # The sign-flip guard above already skips perturbations
                    # that cross zero; a raise here is a state the closure
                    # cannot classify at all. Leaving the column zero was a
                    # silent wrong derivative (issue #271, #280).
                    raise RuntimeError(
                        f"MPCEv2Element '{self.id}': Mynard closure failed on the "
                        f"FD perturbation of port {j} (perturbed mdots={mdot_pert})."
                    ) from exc
                if m_pert.K is None or len(m_pert.K) != len(non_common_idxs):
                    continue
                K_pert = np.zeros(N)
                for k, port_idx in enumerate(non_common_idxs):
                    K_pert[port_idx] = float(m_pert.K[k])
                u_com_pert = abs(float(U_pert[common_idx]))
                q_dyn_pert = 0.5 * rho_com * u_com_pert * u_com_pert
                KQ_pert = K_pert * q_dyn_pert
                dKQ_dmdot[:, j] = (KQ_pert - KQ_base) / mdot_eps

        # Sensitivity of the loss term to the common port's STATIC pressure.
        # This column was absent, leaving a silently zero Jacobian entry -- the
        # same defect PR #230 fixed in ConstantKTeeElement (issue #271).
        # Temperature is derived forward for a MomentumChamberNode rather than
        # solved, so there is no .T unknown to differentiate against.
        #
        # Two paths carry the dependence, and taking only the first is wrong by
        # a few per cent:
        #   (a) q_dyn_com = mdot^2 / (2*rho*A^2) ~ 1/rho, and rho ~ P/T for a
        #       near-ideal gas, so d(q_dyn)/dP = -q_dyn / P.
        #   (b) Mynard's K is a function of the port velocities U = -mdot/(rho*A),
        #       so it moves with rho too.
        # For (b), U depends on P and on mdot through the same 1/rho factor:
        #       dU/dP = -U/P     and     dU/dmdot = U/mdot
        #   =>  dK/dP = -(mdot_com / P_com) * dK/dmdot_com,
        # which lets the existing dKQ/dmdot column supply it. Substituting and
        # using d(q_dyn)/dmdot_com = 2*q_dyn/mdot_com, the two paths combine to
        #       dR_i/dP_com = sign * [ K_i*q_dyn/P - (mdot_com/P) * dKQ_i/dmdot_com ]
        # Note the leading term ends up POSITIVE: path (b) more than cancels the
        # naive -K*q/P. Pinned by test_mpce_v2_jacobian_rows.py.
        P_com = float(states[common_idx].P)
        m_com = float(port_mdots[common_idx])
        p_var_com = f"{self.port_nodes[common_idx]}.P"

        jac: dict[int, dict[str, float]] = {}
        for i in range(N):
            row: dict[str, float] = {
                f"{self.port_nodes[i]}.Pt": 1.0,
                f"{self.id}.P_jct": -1.0,
            }
            if P_com > 0.0:
                dR_dP_com = float(K_per_port[i]) * q_dyn_com / P_com - (m_com / P_com) * float(
                    dKQ_dmdot[i, common_idx]
                )
                if dR_dP_com != 0.0:
                    row[p_var_com] = row.get(p_var_com, 0.0) + K_term_sign * dR_dP_com
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


class ConstantKTeeElement(MPCEv2Element):
    """Junction with FIXED per-port loss coefficients (the "simplest
    model" tier): handbook/datasheet K values instead of the Mynard
    closure.

    Residual on each non-common port i:

        Pt_i - Pt_jct + sign * K_i * q_dyn_com = 0

    with ``sign = +1`` separating (flow_direction='branch') and ``-1``
    joining ('merge'); Pt continuity (K = 0) on the common port; the
    usual junction sum-mass row. ``q_dyn_com`` is the common port's own
    dynamic head, so K follows the Idelchik convention (coefficients
    referenced to the combined-leg velocity head). K does NOT depend on
    the flow split: the junction rows depend on m_com only and the
    network nearly decouples.

    The rows are smooth and even in m_com (no strict/soft barrier
    machinery; the ``strict`` flag is ignored). Like every even-form
    junction, the sign-flipped image of a root is also a root; the
    inherited MPCEv2 ``verify_solution_consistent`` direction guard
    demotes wrong-direction results post-solve.

    Promoted from the certified-audit prototype
    (tmp/mpce_audit_v2_runner.py ConstantKTee, 2026-07-04: 28/32
    converged in 4-14 evaluations).
    """

    def __init__(
        self,
        id: str,
        inlet_nodes: list[str],
        outlet_nodes: list[str],
        inlet_angles_deg: list[float] | None = None,
        outlet_angles_deg: list[float] | None = None,
        port_areas: list[float] | None = None,
        flow_direction: FlowDirection = "branch",
        strict: bool = False,
        K_ports: dict[int, float] | None = None,
    ):
        super().__init__(
            id=id,
            inlet_nodes=inlet_nodes,
            outlet_nodes=outlet_nodes,
            inlet_angles_deg=inlet_angles_deg,
            outlet_angles_deg=outlet_angles_deg,
            port_areas=port_areas,
            flow_direction=flow_direction,
            strict=strict,
        )
        if flow_direction == "branch" and len(inlet_nodes) != 1:
            raise ValueError(
                f"ConstantKTeeElement '{id}': flow_direction='branch' "
                f"requires exactly 1 inlet port, got {len(inlet_nodes)}."
            )
        if flow_direction == "merge" and len(outlet_nodes) != 1:
            raise ValueError(
                f"ConstantKTeeElement '{id}': flow_direction='merge' "
                f"requires exactly 1 outlet port, got {len(outlet_nodes)}."
            )
        # Per-port loss coefficient, keyed by port index (inlets first,
        # then outlets); the common port's entry is ignored (K = 0).
        self.K_ports: dict[int, float] = dict(K_ports) if K_ports else {}

    def _common_port_index(self) -> int:
        # Port order is inlet_nodes + outlet_nodes: the single inlet of a
        # branch tee is port 0; the single outlet of a merge tee is the
        # last port.
        return 0 if self.flow_direction == "branch" else self.N - 1

    def residuals(  # type: ignore[override]
        self,
        states: list[NetworkMixtureState],
        Pt_jct: float,
        port_mdots: list[float],
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        N = self.N
        common = self._common_port_index()
        sign = +1.0 if self.flow_direction == "branch" else -1.0
        A_c = float(self.port_areas[common])
        rho_c = float(states[common].density())
        P_c = float(states[common].P)
        T_c = float(states[common].T)
        m_c = float(port_mdots[common])
        q_dyn = m_c * m_c / (2.0 * rho_c * A_c * A_c)
        dq_dm = m_c / (rho_c * A_c * A_c)
        # q ~ 1/rho with rho ~ P/T (near-ideal gas): d q/dP = -q/P and
        # d q/dT = +q/T at the common port's static state.
        dq_dP = -q_dyn / P_c if P_c > 0.0 else 0.0
        dq_dT = q_dyn / T_c if T_c > 0.0 else 0.0

        residuals: list[float] = []
        jac: dict[int, dict[str, float]] = {}
        m_var_c = f"{self._port_element_ids[common]}.m_dot"
        p_var_c = f"{self.port_nodes[common]}.P"
        t_var_c = f"{self.port_nodes[common]}.T"
        for i in range(N):
            K_i = 0.0 if i == common else float(self.K_ports.get(i, 0.0))
            residuals.append(float(states[i].Pt) - Pt_jct + sign * K_i * q_dyn)
            row: dict[str, float] = {
                f"{self.port_nodes[i]}.Pt": 1.0,
                f"{self.id}.P_jct": -1.0,
            }
            if K_i != 0.0:
                # port_mdots[common] = sign_c * outer_c; q is quadratic in
                # m_c, so chain d q / d outer_c through the port sign.
                row[m_var_c] = row.get(m_var_c, 0.0) + (
                    sign * K_i * dq_dm * self._port_signs[common]
                )
                row[p_var_c] = row.get(p_var_c, 0.0) + sign * K_i * dq_dP
                row[t_var_c] = row.get(t_var_c, 0.0) + sign * K_i * dq_dT
            jac[i] = row

        residuals.append(sum(port_mdots))
        mass_row: dict[str, float] = {}
        for i in range(N):
            mv = f"{self._port_element_ids[i]}.m_dot"
            mass_row[mv] = mass_row.get(mv, 0.0) + self._port_signs[i]
        jac[N] = mass_row
        return residuals, jac

    def diagnostics(  # type: ignore[override]
        self,
        states: list[NetworkMixtureState],
        Pt_jct: float,
        port_mdots: list[float] | None = None,
    ) -> dict[str, float]:
        """Parent per-port fields plus the FIXED K values actually used.

        Deliberately skips MPCEv2's diagnostics (which re-evaluate the
        Mynard closure -- not the model in use here). For the 3-port tee
        the two non-common ports are aliased ``K_straight`` / ``K_branch``
        in port order, matching the GUI wiring.
        """
        diag = MultiPortChamberElement.diagnostics(self, states, Pt_jct, port_mdots)
        common = self._common_port_index()
        noncommon = [i for i in range(self.N) if i != common]
        for i in noncommon:
            diag[f"port_{i}_K"] = float(self.K_ports.get(i, 0.0))
        if self.N == 3:
            diag["K_straight"] = float(self.K_ports.get(noncommon[0], 0.0))
            diag["K_branch"] = float(self.K_ports.get(noncommon[1], 0.0))
        if port_mdots is not None and len(port_mdots) == self.N:
            m_c = float(port_mdots[common])
            if abs(m_c) > 1e-12 and self.N == 3:
                diag["mass_flow_ratio"] = abs(float(port_mdots[noncommon[1]])) / abs(m_c)
        return diag
