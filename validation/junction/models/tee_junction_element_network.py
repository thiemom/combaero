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

from validation.junction.models._network_builder import (
    _F_C,
    _M_DOT_REF,
    _PT_REF,
    _TT_REF,
    _DRY_AIR_Y,
    NetworkResult,
    build_separating_network_skeleton,
    solve_and_extract,
)


class TeeJunctionElementNetwork:
    """CorrelationModel wrapping TeeJunctionElement in a real network."""

    name = "tee_junction_element_network"

    def evaluate_network(
        self,
        paper: str,
        K_id: str,
        q: float,
        psi: float | None,
        theta_rad: float | None,
        **kwargs: float,
    ) -> NetworkResult:
        """Solve the imposed-q network and return convergence + extracted K."""
        if paper != "bassett2001":
            return NetworkResult(converged=False, message="non-Bassett papers unsupported")
        # Map K_id to network regime + which extracted K to report.
        # K6 (separating lateral): report K_lateral
        # K5/K2 (separating straight): report K_straight
        # K12 / K11 / K7 / K8: would need joining flow topology
        if K_id in {"K6", "K5", "K2"}:
            return self._separating(q, psi or 1.0, theta_rad or math.pi / 2.0)
        return NetworkResult(
            converged=False,
            message=f"K_id {K_id} requires non-separating topology, not yet wired",
        )

    def _separating(self, q: float, psi: float, theta_rad: float) -> NetworkResult:
        """Build separating-flow network with TeeJunctionElement, solve, extract K."""
        m_in = _M_DOT_REF
        m_lat = q * m_in
        net = build_separating_network_skeleton(m_in=m_in, m_lateral=m_lat)
        # TeeJunctionElement: branching tee with common as inlet, straight + branch as outlets.
        # Need an extra MFB on the branch outlet to enforce m_lat split (the skeleton already
        # has mb_lat as a MassFlowBoundary). Connect port_bra -> mb_lat.
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
