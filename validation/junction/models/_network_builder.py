"""
Shared helpers for the network-mode validation adapters.

The network mode tests the *wrapped* junction element (TeeJunctionElement,
MultiPortChamberElement) as it behaves inside a NetworkSolver, not the raw
K math. Both adapters use the same imposed-q topology so the comparison is
apples-to-apples; only the junction element itself differs.

Imposed-q topology (separating flow, single-inlet common):

    mb_in (MFB, m_in)  ->  port_com  ->  junction  ->  port_str -> pb_str (PB)
                                                  \\->  port_bra -> mb_lat (MFB, m_lat)

Mass balance forces m_str = m_in - m_lat; PB at the straight end sets the
pressure datum. K6 = (Pt_com - Pt_bra) / q_dyn_com, K2 = (Pt_com - Pt_str) /
q_dyn_com (with q_dyn_com computed from converged port_com state).

For joining flow the topology is mirrored: two MFB inlets feed the junction,
one PB outlet on the common branch.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Literal

import combaero as cb
from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_F_C = 0.01  # 100 cm^2 -- standard area for the validation networks
_M_DOT_REF = 0.1  # kg/s -- low-Mach inlet flow
_TT_REF = 300.0  # K
_PT_REF = 1.0e5  # Pa


@dataclass
class NetworkResult:
    """Output of a single network-mode evaluation."""

    converged: bool
    K_lateral: float | None = None
    K_straight: float | None = None
    residual_norm: float = math.inf
    wall_time_s: float = 0.0
    message: str = ""


def _x_air() -> list[float]:
    """Dry-air mole fractions (cached call to cb.mass_to_mole)."""
    return list(cb.mass_to_mole(_DRY_AIR_Y))


def _extract_K(
    sol: dict[str, float],
    *,
    common_node: str,
    straight_node: str,
    lateral_node: str,
    m_dot_ref: float,
    area: float,
) -> tuple[float, float]:
    """Extract K_lateral and K_straight from a converged solution.

    K = (Pt_common - Pt_branch) / (0.5 * rho_common * u_common^2),
    where Pt is the stagnation pressure at the junction face of each branch.
    """
    P_com = sol[f"{common_node}.P"]
    Pt_com = sol[f"{common_node}.Pt"]
    Pt_str = sol[f"{straight_node}.Pt"]
    Pt_lat = sol[f"{lateral_node}.Pt"]
    X = _x_air()
    rho_com = float(cb.density(_TT_REF, P_com, X))
    u_com = m_dot_ref / (rho_com * area)
    q_dyn_com = 0.5 * rho_com * u_com * u_com
    return (
        (Pt_com - Pt_lat) / q_dyn_com,
        (Pt_com - Pt_str) / q_dyn_com,
    )


# ---------------------------------------------------------------------------
# Network builders
# ---------------------------------------------------------------------------


def build_separating_network_skeleton(
    *,
    m_in: float = _M_DOT_REF,
    m_lateral: float | None = None,
    Pt_ref: float = _PT_REF,
) -> FlowNetwork:
    """Build the boundary + port-node skeleton for separating flow.

    Returns a FlowNetwork with the standard MFB+MFB+PB boundary set and
    the three port MCNs (common, straight, branch). The caller adds the
    junction element + any per-port loss elements + LosslessConnectionElement
    pieces to close the topology.

    All port areas are equal to `_F_C`.
    """
    if m_lateral is None:
        m_lateral = 0.5 * m_in
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(MassFlowBoundary("mb_in", m_dot=m_in, Tt=_TT_REF, Y=Y))
    net.add_node(MassFlowBoundary("mb_lat", m_dot=m_lateral, Tt=_TT_REF, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_ref, Tt=_TT_REF, Y=Y))
    net.add_node(MomentumChamberNode("port_com", area=_F_C))
    net.add_node(MomentumChamberNode("port_str", area=_F_C))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C))
    net.add_element(LosslessConnectionElement("lc_in", "mb_in", "port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    return net


def solve_and_extract(
    net: FlowNetwork,
    *,
    common_node: str,
    straight_node: str,
    lateral_node: str,
    m_dot_ref: float,
    area: float = _F_C,
) -> NetworkResult:
    """Run NetworkSolver, return convergence + K diagnostics."""
    solver = NetworkSolver(net)
    t0 = time.perf_counter()
    try:
        sol = solver.solve()
    except Exception as e:
        return NetworkResult(
            converged=False, message=str(e)[:120], wall_time_s=time.perf_counter() - t0
        )
    wall_time = time.perf_counter() - t0
    converged = bool(sol.get("__success__", False))
    res_norm = float(sol.get("__residual_norm__", math.inf))
    if not converged:
        return NetworkResult(
            converged=False,
            residual_norm=res_norm,
            wall_time_s=wall_time,
            message=sol.get("__message__", "")[:120],
        )
    try:
        K_lat, K_str = _extract_K(
            sol,
            common_node=common_node,
            straight_node=straight_node,
            lateral_node=lateral_node,
            m_dot_ref=m_dot_ref,
            area=area,
        )
    except KeyError as e:
        return NetworkResult(
            converged=False,
            residual_norm=res_norm,
            wall_time_s=wall_time,
            message=f"missing solution key {e}",
        )
    return NetworkResult(
        converged=True,
        K_lateral=K_lat,
        K_straight=K_str,
        residual_norm=res_norm,
        wall_time_s=wall_time,
    )
