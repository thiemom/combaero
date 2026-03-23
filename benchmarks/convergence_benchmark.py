#!/usr/bin/env python3
"""Benchmark: Solver convergence for high-Mach coupled compressible networks.

This benchmark targets the regime where the Fanno pipe approaches choking
and the initial-guess propagation becomes critical.
"""

import numpy as np
import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
    WallConnection,
)

def _area(diameter: float) -> float:
    return float(np.pi * (diameter / 2.0) ** 2)

def _mdot_for_mach(mach: float, T: float, P_ref: float, X: list[float], area: float) -> float:
    rho = cb.density(T, P_ref, X)
    a = cb.speed_of_sound(T, X)
    return float(mach * rho * a * area)

_T_HOT = 740.0
_T_COLD = 320.0
_P_OUT = 2.0e5

def _build_fully_coupled_network(mach_target: float) -> FlowNetwork:
    x_air = cb.standard_dry_air_composition()
    y_air = list(cb.mole_to_mass(x_air))

    t_hot = _T_HOT
    t_cold = _T_COLD
    p_out_hot = _P_OUT
    p_out_cold = _P_OUT
    d_hot = 0.018
    d_cold = 0.014
    length = 2.0

    mdot_hot = _mdot_for_mach(mach_target, t_hot, p_out_hot, x_air, _area(d_hot))
    mdot_cold = 0.65 * _mdot_for_mach(mach_target, t_cold, p_out_cold, x_air, _area(d_cold))

    net = FlowNetwork()
    net.add_node(MassFlowBoundary("hot_inlet", m_dot=mdot_hot, T_total=t_hot, Y=y_air))
    net.add_node(PlenumNode("hot_plenum"))
    net.add_node(PressureBoundary("hot_outlet", P_total=p_out_hot, T_total=t_hot, Y=y_air))

    net.add_node(MassFlowBoundary("cold_inlet", m_dot=mdot_cold, T_total=t_cold, Y=y_air))
    net.add_node(PlenumNode("cold_plenum"))
    net.add_node(PressureBoundary("cold_outlet", P_total=p_out_cold, T_total=t_cold, Y=y_air))

    hot_surface = ConvectiveSurface(area=np.pi * d_hot * length, model=SmoothModel())
    cold_surface = ConvectiveSurface(area=np.pi * d_cold * length, model=SmoothModel())

    net.add_element(PipeElement(id="hot_pipe", from_node="hot_inlet", to_node="hot_plenum", length=length, diameter=d_hot, roughness=5.0e-5, regime="compressible_fanno", surface=hot_surface))
    net.add_element(OrificeElement(id="hot_orifice", from_node="hot_plenum", to_node="hot_outlet", Cd=0.72, area=_area(0.030), regime="compressible"))
    net.add_element(PipeElement(id="cold_pipe", from_node="cold_inlet", to_node="cold_plenum", length=length, diameter=d_cold, roughness=5.0e-5, regime="compressible_fanno", surface=cold_surface))
    net.add_element(OrificeElement(id="cold_orifice", from_node="cold_plenum", to_node="cold_outlet", Cd=0.72, area=_area(0.025), regime="compressible"))

    net.add_wall(WallConnection(id="coupling_wall", element_a="hot_pipe", element_b="cold_pipe", wall_thickness=0.002, wall_conductivity=20.0))
    net.thermal_coupling_enabled = True
    return net

if __name__ == "__main__":
    mach_values = np.linspace(0.20, 1.50, 14)
    print(f"Running high-Mach convergence benchmark ({len(mach_values)} points)...")

    x0_prev = None
    for mach_target in mach_values:
        net = _build_fully_coupled_network(float(mach_target))
        solver = NetworkSolver(net)

        # Try fresh, then warm-start from previous point
        sol = solver.solve(method="hybr")
        if not sol["__success__"] and x0_prev is not None:
            sol = solver.solve(method="hybr", x0=x0_prev)

        status = "SUCCESS" if sol["__success__"] else "FAILED"
        print(f"M_target={mach_target:.3f}: {status} (|F|={sol['__final_norm__']:.3e})")
        if sol["__success__"]:
            x0_prev = solver._last_x

    print("Benchmark complete.")
