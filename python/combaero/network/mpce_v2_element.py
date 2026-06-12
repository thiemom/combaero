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

import numpy as np

from combaero.network._mpce_v2_jacobian import dKQ_dmdot_separating_T
from combaero.network.components import (
    MultiPortChamberElement,
    NetworkMixtureState,
)
from validation.junction.models.mynard2010 import junction_loss_coefficient


class MPCEv2Element(MultiPortChamberElement):
    """Mynard Unified0D residual on MPCE-v1's topology framework."""

    def residuals(  # type: ignore[override]
        self,
        states: list[NetworkMixtureState],
        Pt_jct: float,
        port_mdots: list[float],
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        N = self.N
        A = np.array([float(a) for a in self.port_areas])
        # Convert MPCE-v1 branch angles (angle from main axis) to Mynard
        # convention (vessel direction). Inlet ports have +pi added because
        # they point away from the junction in the opposite direction from
        # the canonical inlet flow.
        theta_rad = np.array(
            [
                math.radians(float(t)) + (math.pi if self._port_signs[i] < 0 else 0.0)
                for i, t in enumerate(self.port_angles_deg)
            ]
        )

        # Convert: port_mdots are in junction convention (positive = out of
        # junction). For Mynard: positive U = supplier (into junction).
        # So Mynard's U = -port_mdots / (rho * A).
        rho_port = np.array([float(s.density()) for s in states])
        U_mynard = -np.array(port_mdots) / (rho_port * A)

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

        try:
            mynard = junction_loss_coefficient(U_mynard, A, theta_rad)
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
        if int(np.sum(sup_mask)) == 1:
            common_mask = sup_mask
            non_common_idxs = np.where(col_mask)[0]
        else:
            common_mask = col_mask
            non_common_idxs = np.where(sup_mask)[0]

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
            R_i = Pt_i - Pt_jct + K_per_port[i] * q_dyn_com
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
        if is_canonical_separating_T:
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
                    m_pert = junction_loss_coefficient(U_pert, A, theta_rad)
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
                    row[mdot_var] = row.get(mdot_var, 0.0) + (dKQ_dmdot[i, j] * self._port_signs[j])
            jac[i] = row

        mass_row: dict[str, float] = {}
        for i in range(N):
            mass_var = f"{self._port_element_ids[i]}.m_dot"
            mass_row[mass_var] = mass_row.get(mass_var, 0.0) + self._port_signs[i]
        jac[N] = mass_row

        return residuals, jac
