"""
CorrelationModel adapter scoring Mynard 2010 Unified0D against the dataset.

Maps Bassett K_id + (q, psi, theta) into Mynard's (U, A, theta) inputs
and reads back the appropriate per-branch K. Covers separating K5/K6 and
joining K11/K12; other K_id's involve non-canonical topologies and are
left as None (the scorecard renders them as '-').

This adapter is positioned alongside `PaperCeiling` so the scorecard can
compare a candidate junction model against BOTH the Bassett analytical
reference AND the Mynard Unified0D reference. The latter is the published
ceiling for the model class our `tee_junction.h::K_dat_j_closed` aims to
implement.
"""

from __future__ import annotations

import math

import numpy as np

from validation.junction.models.mynard2010 import junction_loss_coefficient


class MynardAnalytical:
    """Mynard & Valen-Sendstad 2010 Unified0D loss coefficients."""

    name = "mynard_unified0d_analytical"

    def evaluate(
        self,
        paper: str,
        K_id: str,
        q: float,
        psi: float | None,
        theta_rad: float | None,
        M_3: float | None = None,
        **kwargs: float,
    ) -> float | None:
        if paper != "bassett2001":
            return None
        psi = psi if psi is not None else 1.0
        theta_rad = theta_rad if theta_rad is not None else math.pi / 2.0

        # Mynard's K is per-branch from a 3-branch (U, A, theta) input.
        # Map Bassett K_id -> (flow direction, which-K-to-return).
        # Sign convention: positive U = supplier (into junction).
        #
        # Areas: Bassett psi = F_C / F_B, so A_branch = A_com / psi.
        A_com, A_str = 1.0, 1.0
        A_bra = 1.0 / psi
        # Theta convention: common at pi (flow direction into junction from "left"),
        # straight at 0 (flow leaves "right"), branch at theta (lateral).
        theta_arr = np.array([math.pi, 0.0, theta_rad])

        if K_id in {"K5", "K2", "K6"}:
            # Separating: common is supplier, straight + branch are collectors.
            U_arr = np.array([1.0, -(1.0 - q), -q])
            A_arr = np.array([A_com, A_str, A_bra])
            try:
                r = junction_loss_coefficient(U_arr, A_arr, theta_arr)
            except Exception:
                return None
            if r.K is None or len(r.K) < 2:
                return None
            # K[0] = straight collector (K5/K2 analog); K[1] = branch (K6 analog)
            return float(r.K[0]) if K_id in {"K5", "K2"} else float(r.K[1])

        if K_id in {"K11", "K12"}:
            # Joining: common is collector, straight + branch are suppliers.
            # Bassett q convention for K11/K12 = m_branch / m_common
            # (Table 1). So branch supplier flow rate is q*Q_total, straight
            # is (1-q)*Q_total. Common collector takes all flow (Q = -1).
            U_arr = np.array([-1.0, (1.0 - q), q])
            A_arr = np.array([A_com, A_str, A_bra])
            try:
                r = junction_loss_coefficient(U_arr, A_arr, theta_arr)
            except Exception:
                return None
            if r.K is None or len(r.K) < 2:
                return None
            # K[0] = straight supplier (K11 analog); K[1] = branch (K12 analog)
            return float(r.K[0]) if K_id == "K11" else float(r.K[1])

        # K1/K3/K4/K7/K8/K9/K10 use specific Bassett control volumes that
        # don't map 1:1 onto Mynard's symmetric 3-branch geometry. Defer.
        return None
