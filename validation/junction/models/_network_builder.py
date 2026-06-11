"""
Shared helpers for the network-mode validation adapters.

The network mode tests the *wrapped* junction element (TeeJunctionElement,
MultiPortChamberElement) as it behaves inside a NetworkSolver, not the raw
K math.

Three separating-flow topologies are supported. The same junction element
is exercised under all three; only the boundary set changes.

  IMPOSED_Q  (2x MFB + 1x PB)
        mb_in (MFB, m_in)  ->  port_com -> junction -> port_str -> pb_str (PB)
                                                    \\-> port_bra -> mb_lat (MFB, m_lat)
        Mass flows imposed; pressure datum at the straight outlet.

  THREE_PB  (3x PB)
        pb_com (PB) -> port_com -> junction -> port_str -> pb_str (PB)
                                            \\-> port_bra -> pb_bra (PB)
        Pressures imposed at every boundary; junction determines split.
        BC pressures are derived from the analytical K so the target q is
        a converged root (tests solver robustness, not K accuracy).

  MFB_TWO_PB  (1x MFB + 2x PB)
        mb_in (MFB) -> port_com -> junction -> port_str -> pb_str (PB)
                                            \\-> port_bra -> pb_bra (PB)
        Inlet flow imposed; outlet pressures set the split.
        Production-realistic BC set.
"""

from __future__ import annotations

import math
import time
import warnings
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


Topology = Literal["imposed_q", "three_pb", "mfb_two_pb"]
ALL_TOPOLOGIES: tuple[Topology, ...] = ("imposed_q", "three_pb", "mfb_two_pb")


@dataclass
class NetworkResult:
    """Output of a single network-mode evaluation."""

    converged: bool
    K_lateral: float | None = None
    K_straight: float | None = None
    residual_norm: float = math.inf
    wall_time_s: float = 0.0
    message: str = ""
    q_converged: float | None = None  # actual q the solver found, may drift from target


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


def _q_dyn_ref(m_dot: float = _M_DOT_REF, area: float = _F_C, Pt: float = _PT_REF) -> float:
    """Reference dynamic head at the low-Mach inlet (ρ·u²/2)."""
    X = _x_air()
    rho = float(cb.density(_TT_REF, Pt, X))
    u = m_dot / (rho * area)
    return 0.5 * rho * u * u


def build_separating_network_skeleton(
    *,
    m_in: float = _M_DOT_REF,
    m_lateral: float | None = None,
    Pt_ref: float = _PT_REF,
) -> FlowNetwork:
    """IMPOSED_Q skeleton: MFB on common inlet + MFB on lateral outlet + PB on straight.

    Mass flows imposed; PB at the straight sets the pressure datum.
    Caller adds the junction element + per-port loss + lossless wiring.
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


def build_separating_three_pb_skeleton(
    *,
    K_lateral_target: float,
    K_straight_target: float,
    Pt_com: float = _PT_REF,
) -> FlowNetwork:
    """THREE_PB skeleton: PB on every branch.

    Pressure drops imply the target K relations: Pt_com - Pt_lat = K_lat * q_dyn
    and Pt_com - Pt_str = K_str * q_dyn. The network is solved blind; if the
    junction element's K differs from the target, the solver finds the q where
    the model's K matches the imposed Pt drops.

    BC pressures are based on the reference dynamic head, which is fixed
    rather than re-derived per (m, A) -- the absolute value drops out of the
    K extraction.
    """
    q_dyn = _q_dyn_ref(Pt=Pt_com)
    Pt_lat = Pt_com - K_lateral_target * q_dyn
    Pt_str = Pt_com - K_straight_target * q_dyn
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(PressureBoundary("pb_com", Pt=Pt_com, Tt=_TT_REF, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_str, Tt=_TT_REF, Y=Y))
    net.add_node(PressureBoundary("pb_bra", Pt=Pt_lat, Tt=_TT_REF, Y=Y))
    net.add_node(MomentumChamberNode("port_com", area=_F_C))
    net.add_node(MomentumChamberNode("port_str", area=_F_C))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C))
    net.add_element(LosslessConnectionElement("lc_com", "pb_com", "port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    return net


def build_separating_mfb_two_pb_skeleton(
    *,
    m_in: float = _M_DOT_REF,
    K_lateral_target: float,
    K_straight_target: float,
    Pt_com: float = _PT_REF,
) -> FlowNetwork:
    """MFB_TWO_PB skeleton: MFB on common inlet + PB on each outlet.

    Inlet flow imposed at m_in; outlet pressures sized from analytical K
    so the target q is a converged root. The most production-realistic BC
    set for steady-state design.
    """
    q_dyn = _q_dyn_ref(m_dot=m_in, Pt=Pt_com)
    Pt_lat = Pt_com - K_lateral_target * q_dyn
    Pt_str = Pt_com - K_straight_target * q_dyn
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(MassFlowBoundary("mb_in", m_dot=m_in, Tt=_TT_REF, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_str, Tt=_TT_REF, Y=Y))
    net.add_node(PressureBoundary("pb_bra", Pt=Pt_lat, Tt=_TT_REF, Y=Y))
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
        # Suppress the per-call non-convergence UserWarnings -- the scorecard
        # already captures success/residual norm explicitly.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
