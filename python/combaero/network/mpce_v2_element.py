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

        # ITERATION-1 PROTOTYPE: use Mynard's C per branch with the global
        # supplier-side q_dyn. This converges across the slices MPCE-v1
        # converges at (and at theta=120 where v1 fails to converge), but
        # the extracted K is biased because C and q_dyn here are
        # mismatched in normalization (C is collector-relative, q_dyn is
        # supplier-relative). Iteration-2 will either:
        #   - switch to mynard.K with common-side q_dyn (correct
        #     normalization but currently fails convergence -- needs
        #     residual scaling or analytical Jacobian),
        #   - or keep C and switch to per-port collector-side q_dyn (also
        #     currently fails convergence -- to be investigated).
        # The runner scorecard quantifies the bias so we can iterate
        # against measured data.
        C_per_port = mynard.C.copy()

        # Supplier-side dynamic head (mass-flow-weighted; equals u_sup^2
        # for the single-supplier case).
        sup_idxs = np.where(sup_mask)[0]
        u_sup_sq_avg = float(
            np.sum(U_mynard[sup_idxs] ** 2 * rho_port[sup_idxs] * A[sup_idxs])
            / np.sum(rho_port[sup_idxs] * A[sup_idxs])
        )
        rho_sup_avg = float(
            np.sum(rho_port[sup_idxs] * A[sup_idxs] * U_mynard[sup_idxs])
            / np.sum(A[sup_idxs] * U_mynard[sup_idxs])
        )
        q_dyn_ref = 0.5 * rho_sup_avg * u_sup_sq_avg

        residuals: list[float] = []
        for i in range(N):
            Pt_i = float(states[i].Pt)
            R_i = Pt_i - Pt_jct + C_per_port[i] * q_dyn_ref
            residuals.append(R_i)
        residuals.append(sum(port_mdots))

        # Stub Jacobian -- solver falls back to numerical differencing.
        return residuals, {}
