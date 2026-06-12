"""
Network-mode adapter for MPCE-v2 (Mynard Unified0D residual).

Mirrors `MPCEv1Network` but builds the junction with `MPCEv2Element`
(subclass of MPCE-v1 with Mynard's K-based residual replacing the
cross-coupling impulse formula). No BorderCarnotLossElement on the
lateral branch because Mynard's K already captures the full loss.
"""

from __future__ import annotations

import math

from combaero.network import LosslessConnectionElement
from combaero.network.mpce_v2_element import MPCEv2Element
from validation.junction.models import bassett2001
from validation.junction.models._network_builder import (
    _F_C,
    _M_DOT_REF,
    ALL_TOPOLOGIES,
    NetworkResult,
    Topology,
    build_joining_imposed_q_skeleton,
    build_joining_mfb_two_pb_skeleton,
    build_joining_three_pb_skeleton,
    build_separating_mfb_two_pb_skeleton,
    build_separating_network_skeleton,
    build_separating_three_pb_skeleton,
    solve_and_extract,
)


class MPCEv2Network:
    """MPCE-v2 (Mynard residual) in the three separating topologies."""

    name = "mpce_v2_network"

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
            return NetworkResult(converged=False, message=f"topology {topology!r} not wired")
        if K_id in {"K6", "K5", "K2"}:
            return self._separating(q, psi or 1.0, theta_rad or math.pi / 2.0, topology)
        if K_id in {"K11", "K12"}:
            return self._joining(q, psi or 1.0, theta_rad or math.pi / 2.0, topology)
        return NetworkResult(
            converged=False,
            message=f"K_id {K_id} not yet wired for MPCE-v2",
        )

    def _separating(
        self, q: float, psi: float, theta_rad: float, topology: Topology
    ) -> NetworkResult:
        m_in = _M_DOT_REF
        A_bra = _F_C / psi  # Bassett psi = F_C/F_B, so A_branch = A_com/psi

        if topology == "imposed_q":
            net = build_separating_network_skeleton(m_in=m_in, m_lateral=q * m_in)
            lateral_terminal = "mb_lat"
        else:
            K_lat = bassett2001.K6(q, psi, theta_rad)
            K_str = bassett2001.K5(q)
            if topology == "three_pb":
                net = build_separating_three_pb_skeleton(
                    K_lateral_target=K_lat, K_straight_target=K_str
                )
            else:
                net = build_separating_mfb_two_pb_skeleton(
                    K_lateral_target=K_lat, K_straight_target=K_str
                )
            lateral_terminal = "pb_bra"

        # Replace the default branch port area to match psi.
        net.nodes["port_bra"].area = A_bra

        net.add_element(LosslessConnectionElement("lc_bra", "port_bra", lateral_terminal))
        net.add_element(
            MPCEv2Element(
                id="jct",
                inlet_nodes=["port_com"],
                outlet_nodes=["port_str", "port_bra"],
                inlet_angles_deg=[0.0],
                outlet_angles_deg=[0.0, math.degrees(theta_rad)],
                port_areas=[_F_C, _F_C, A_bra],
                flow_direction="branch",
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

    def _joining(
        self, q: float, psi: float, theta_rad: float, topology: Topology
    ) -> NetworkResult:
        """Build a joining-flow MPCE-v2 network.

        Geometry mirrors the separating case (str at 0deg, bra at theta) but
        flow reverses: str and bra are inlets, com is the outlet. The
        ``MPCEv2Element`` residual reads the flow direction from the
        signed mass flows at runtime, so the same element handles both.
        """
        m_in = _M_DOT_REF
        A_bra = _F_C / psi
        if topology == "imposed_q":
            net = build_joining_imposed_q_skeleton(m_in=m_in, m_lateral=q * m_in)
        else:
            K11 = bassett2001.K11_corr(q, psi, theta_rad)
            K12 = bassett2001.K12_corr(q, psi, theta_rad)
            if topology == "three_pb":
                net = build_joining_three_pb_skeleton(K11_target=K11, K12_target=K12)
            else:
                net = build_joining_mfb_two_pb_skeleton(
                    m_in=m_in, K11_target=K11, K12_target=K12
                )

        net.nodes["port_bra"].area = A_bra
        # Branch supplier connects to its MFB/PB terminal.
        if topology == "imposed_q":
            net.add_element(LosslessConnectionElement("lc_bra", "mb_bra", "port_bra"))
        else:
            net.add_element(LosslessConnectionElement("lc_bra", "pb_bra", "port_bra"))

        # MPCE-v1 inlet/outlet declaration: for joining, str and bra are
        # inlets (flow into junction from supply), com is the outlet.
        net.add_element(
            MPCEv2Element(
                id="jct",
                inlet_nodes=["port_str", "port_bra"],
                outlet_nodes=["port_com"],
                inlet_angles_deg=[0.0, math.degrees(theta_rad)],
                outlet_angles_deg=[0.0],
                port_areas=[_F_C, A_bra, _F_C],
                flow_direction="merge",
            )
        )
        return solve_and_extract(
            net,
            common_node="port_com",
            straight_node="port_str",
            lateral_node="port_bra",
            m_dot_ref=m_in,
            area=_F_C,
            flow_direction="joining",
        )
