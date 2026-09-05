"""
Network-mode adapter for `MultiPortChamberElement` (MPCE-v1).

Mirrors `TeeJunctionElementNetwork` but uses the momentum-CV junction
landed in PR #181. Same imposed-q topology so a head-to-head comparison
against the K-closure is apples-to-apples; only the junction element
differs.

Coverage limitations match the MPCE-v1 milestone scope (see
`python/tests/test_momentum_cv_tier1_bassett.py`):
  - K5/K2 (separating straight) and K6 (separating lateral) are both
    extracted from the same separating network and scored against the
    Bassett measured curves. K_straight was previously declared
    unsupported because the sin^2(theta_i) gate collapsed the collinear
    port to static pressure equality; the gate was removed in issue #272.
  - K12 (joining): wrong cross-coupling sign for joining-flow direction
    convention. Tracked as MPCE-v2 follow-up.

The adapter reports converged=False with a diagnostic message for K_id's
known not to converge cleanly, so the scorecard shows convergence rate
rather than fabricated bad K values.
"""

from __future__ import annotations

import math
from dataclasses import replace

from combaero.network import (
    BorderCarnotLossElement,
    LosslessConnectionElement,
    MomentumChamberNode,
    MultiPortChamberElement,
)

from validation.junction.models import bassett2001
from validation.junction.models._network_builder import (
    _F_C,
    _M_DOT_REF,
    ALL_TOPOLOGIES,
    NetworkResult,
    Topology,
    build_separating_mfb_two_pb_skeleton,
    build_separating_network_skeleton,
    build_separating_three_pb_skeleton,
    solve_and_extract,
)


class MPCEv1Network:
    """CorrelationModel wrapping MultiPortChamberElement in a real network."""

    name = "mpce_v1_network"

    SUPPORTED_TOPOLOGIES: tuple[Topology, ...] = ALL_TOPOLOGIES

    def evaluate_network(
        self,
        paper: str,
        K_id: str,
        q: float,
        psi: float | None,
        theta_rad: float | None,
        topology: Topology = "imposed_q",
        **kwargs: float,
    ) -> NetworkResult:
        if paper != "bassett2001":
            return NetworkResult(converged=False, message="non-Bassett papers unsupported")
        if topology not in self.SUPPORTED_TOPOLOGIES:
            return NetworkResult(
                converged=False, message=f"topology {topology!r} not wired for this adapter"
            )
        # K5/K2 (straight) and K6 (lateral) come from the same separating
        # network; the runner selects which coefficient it reads. K_straight
        # used to be declared unsupported because the sin^2(theta_i) gate
        # collapsed the collinear port to static equality; with the gate
        # removed (issue #272) it is a real prediction and is scored.
        if K_id in {"K5", "K2", "K6"}:
            # Bassett indexes each coefficient by the fraction in ITS OWN leg
            # (equivalences.py), so a K5/K2 file's q is the STRAIGHT fraction
            # while K6's is the lateral one. The network is built from the
            # lateral fraction, so straight-leg coefficients need 1-q. Feeding
            # the raw q builds a mirrored operating point (issue #272).
            q_lateral = 1.0 - q if K_id in {"K5", "K2"} else q
            result = self._separating(
                q_lateral, psi or 1.0, theta_rad or math.pi / 2.0, topology
            )
            if K_id in {"K5", "K2"} and result.q_converged is not None:
                # The achieved operating point takes the same transform as the
                # target did on the way in, or the drift is measured on a
                # mirrored axis. Only reachable since q_converged was populated
                # (#290), which is why the input half predates it.
                result = replace(result, q_converged=1.0 - result.q_converged)
            return result
        if K_id in {"K11", "K12"}:
            return NetworkResult(
                converged=False,
                message="joining flow: MPCE-v1 cross-coupling sign issue (MPCE-v2 task)",
            )
        return NetworkResult(converged=False, message=f"K_id {K_id} not in MPCE-v1 scope")

    def _separating(
        self, q: float, psi: float, theta_rad: float, topology: Topology
    ) -> NetworkResult:
        """MPCE + per-port BorderCarnotLoss in any of the 3 separating topologies."""
        if topology == "imposed_q":
            m_in = _M_DOT_REF
            net = build_separating_network_skeleton(m_in=m_in, m_lateral=q * m_in)
            lateral_terminal = "mb_lat"
        else:
            K_lat = bassett2001.K6(q, psi, theta_rad)
            K_str = bassett2001.K5(q)
            if topology == "three_pb":
                net = build_separating_three_pb_skeleton(
                    K_lateral_target=K_lat, K_straight_target=K_str
                )
            else:  # mfb_two_pb
                net = build_separating_mfb_two_pb_skeleton(
                    K_lateral_target=K_lat, K_straight_target=K_str
                )
            lateral_terminal = "pb_bra"
        # Add post-loss MCN on the lateral branch (MPCE pattern).
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
        net.add_element(LosslessConnectionElement("lc_bra", "port_bra_post", lateral_terminal))
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
        return solve_and_extract(
            net,
            common_node="port_com",
            straight_node="port_str",
            lateral_node="port_bra_post",
            m_dot_ref=_M_DOT_REF,
            area=_F_C,
        )
