"""
Network-mode adapter for `MultiPortChamberElement` (MPCE-v1).

Mirrors `TeeJunctionElementNetwork` but uses the momentum-CV junction
landed in PR #181. Same imposed-q topology so a head-to-head comparison
against the K-closure is apples-to-apples; only the junction element
differs.

Coverage limitations match the MPCE-v1 milestone scope (see
`python/tests/test_momentum_cv_tier1_bassett.py`):
  - K6 (separating lateral): MPCE+cross-coupling reproduces Bassett K6
    within ~15% under bounded LM solver. Default NetworkSolver lands on
    a different basin for imposed-q topology; expect non-convergence or
    biased K extraction.
  - K5 (separating straight): not reproduced -- sin^2(theta=0) loses
    axial dynamic-head coupling.
  - K12 (joining): wrong cross-coupling sign for joining-flow direction
    convention. Tracked as MPCE-v2 follow-up.

The adapter reports converged=False with a diagnostic message for K_id's
known not to converge cleanly, so the scorecard shows convergence rate
rather than fabricated bad K values.
"""

from __future__ import annotations

import math

from combaero.network import (
    BorderCarnotLossElement,
    LosslessConnectionElement,
    MomentumChamberNode,
    MultiPortChamberElement,
)

from validation.junction.models._network_builder import (
    _F_C,
    _M_DOT_REF,
    NetworkResult,
    build_separating_network_skeleton,
    solve_and_extract,
)


class MPCEv1Network:
    """CorrelationModel wrapping MultiPortChamberElement in a real network."""

    name = "mpce_v1_network"

    def evaluate_network(
        self,
        paper: str,
        K_id: str,
        q: float,
        psi: float | None,
        theta_rad: float | None,
        **kwargs: float,
    ) -> NetworkResult:
        if paper != "bassett2001":
            return NetworkResult(converged=False, message="non-Bassett papers unsupported")
        if K_id == "K6":
            return self._separating_lateral(q, psi or 1.0, theta_rad or math.pi / 2.0)
        if K_id in {"K5", "K2"}:
            return NetworkResult(
                converged=False,
                message="K_straight: MPCE-v1 known limitation (sin^2(0)=0 loses coupling)",
            )
        if K_id in {"K11", "K12"}:
            return NetworkResult(
                converged=False,
                message="joining flow: MPCE-v1 cross-coupling sign issue (MPCE-v2 task)",
            )
        return NetworkResult(converged=False, message=f"K_id {K_id} not in MPCE-v1 scope")

    def _separating_lateral(self, q: float, psi: float, theta_rad: float) -> NetworkResult:
        """Imposed-q separating network with MPCE + per-port BorderCarnotLoss.

        Mirrors test_momentum_cv_tier1_bassett._build_imposed_q_network: MPCE
        as the junction, BorderCarnotLossElement on the lateral branch with
        the geometric angle, lossless straight branch.
        """
        m_in = _M_DOT_REF
        m_lat = q * m_in
        net = build_separating_network_skeleton(m_in=m_in, m_lateral=m_lat)
        # Add a post-loss MCN on the lateral branch (MPCE pattern).
        net.add_node(MomentumChamberNode("port_bra_post", area=_F_C))
        net.add_element(
            BorderCarnotLossElement(
                "loss_bra",
                from_node="port_bra",
                to_node="port_bra_post",
                delta_geom_deg=math.degrees(theta_rad),
                area=_F_C,
            )
        )
        net.add_element(LosslessConnectionElement("lc_bra", "port_bra_post", "mb_lat"))
        net.add_element(
            MultiPortChamberElement(
                id="jct",
                inlet_nodes=["port_com"],
                outlet_nodes=["port_str", "port_bra"],
                inlet_angles_deg=[0.0],
                outlet_angles_deg=[0.0, math.degrees(theta_rad)],
                port_areas=[_F_C, _F_C, _F_C],
            )
        )
        # Extract K against the POST-LOSS node so the loss element is captured.
        return solve_and_extract(
            net,
            common_node="port_com",
            straight_node="port_str",
            lateral_node="port_bra_post",
            m_dot_ref=m_in,
            area=_F_C,
        )
