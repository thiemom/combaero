"""
Network-mode adapter for `TeeJunctionElement` (the v3 K-closure).

Wraps the Python `TeeJunctionElement` inside an imposed-q FlowNetwork and
extracts K from the converged state. This tests the WRAPPED element
behaviour (solver coupling, residual scaling, jacobian propagation), not
just the raw K math which `TeeJunctionRaw` already covers.

Coverage: only separating-flow K6 (lateral) and K2/K5 (straight) at
psi=1.0. The branching-tee API is symmetric so the same network builder
works for both; we report K_lateral = K6 and K_straight = K5.

Joining-flow K11/K12 requires a different topology (two MFB inlets, one
PB outlet). This is implemented similarly but skipped for K cases that
don't match the expected (paper, K_id).
"""

from __future__ import annotations

import math

from combaero.network import (
    BorderCarnotLossElement,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    PressureBoundary,
    TeeJunctionElement,
)

from validation.junction.models import bassett2001
from validation.junction.models._network_builder import (
    _F_C,
    _M_DOT_REF,
    _PT_REF,
    _TT_REF,
    _DRY_AIR_Y,
    ALL_TOPOLOGIES,
    NetworkResult,
    Topology,
    build_separating_mfb_two_pb_skeleton,
    build_separating_network_skeleton,
    build_separating_three_pb_skeleton,
    solve_and_extract,
)


class TeeJunctionElementNetwork:
    """CorrelationModel wrapping TeeJunctionElement in a real network."""

    name = "tee_junction_element_network"

    # Topologies the adapter has wiring for. Add more here as joining-flow
    # topologies are implemented.
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
        if K_id not in {"K6", "K5", "K2"}:
            return NetworkResult(
                converged=False,
                message=f"K_id {K_id} requires non-separating topology, not yet wired",
            )
        return self._separating(q, psi or 1.0, theta_rad or math.pi / 2.0, topology)

    def _separating(
        self, q: float, psi: float, theta_rad: float, topology: Topology
    ) -> NetworkResult:
        if topology == "imposed_q":
            return self._separating_imposed_q(q, psi, theta_rad)
        # For PB-based topologies, use Bassett analytical K to set up BCs.
        K_lat = bassett2001.K6(q, psi, theta_rad)
        K_str = bassett2001.K5(q)
        if topology == "three_pb":
            net = build_separating_three_pb_skeleton(
                K_lateral_target=K_lat, K_straight_target=K_str
            )
            common_in = ("pb_com", "port_com", "lc_com")
        else:  # mfb_two_pb
            net = build_separating_mfb_two_pb_skeleton(
                K_lateral_target=K_lat, K_straight_target=K_str
            )
            common_in = None  # no extra wiring; mb_in already connected by skeleton
        net.add_element(LosslessConnectionElement("lc_bra", "port_bra", "pb_bra"))
        net.add_element(
            TeeJunctionElement(
                id="tee",
                common_node="port_com",
                straight_node="port_str",
                branch_node="port_bra",
                theta=theta_rad,
                F_C=_F_C,
                psi=psi,
                tee_type="branching",
            )
        )
        return solve_and_extract(
            net,
            common_node="port_com",
            straight_node="port_str",
            lateral_node="port_bra",
            m_dot_ref=_M_DOT_REF,
            area=_F_C,
        )

    def _separating_imposed_q(self, q: float, psi: float, theta_rad: float) -> NetworkResult:
        m_in = _M_DOT_REF
        m_lat = q * m_in
        net = build_separating_network_skeleton(m_in=m_in, m_lateral=m_lat)
        net.add_element(LosslessConnectionElement("lc_bra", "port_bra", "mb_lat"))
        net.add_element(
            TeeJunctionElement(
                id="tee",
                common_node="port_com",
                straight_node="port_str",
                branch_node="port_bra",
                theta=theta_rad,
                F_C=_F_C,
                psi=psi,
                tee_type="branching",
            )
        )
        return solve_and_extract(
            net,
            common_node="port_com",
            straight_node="port_str",
            lateral_node="port_bra",
            m_dot_ref=m_in,
            area=_F_C,
        )
