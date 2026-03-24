#!/usr/bin/env python3
"""Benchmark: Continuation (Warm-Start) vs Cold Start initialization.

This script compares the solver performance when starting from a previous
solution (Continuation) versus starting from scratch (Cold Start) for
a parameter sweep across Mach numbers.
"""

import time
import numpy as np
import combaero as cb
from dataclasses import dataclass
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

@dataclass(frozen=True)
class SweepConfig:
    t_hot_in: float = 750.0
    t_cold_in: float = 320.0
    p_out_hot: float = 2.0e5
    p_out_cold: float = 2.0e5
    d_pipe_hot: float = 0.018
    d_pipe_cold: float = 0.014
    l_pipe: float = 1.2
    roughness: float = 1.0e-5
    cd: float = 0.72
    d_orifice_hot: float = 0.012
    d_orifice_cold: float = 0.010
    wall_thickness: float = 0.002
    wall_k: float = 20.0

def _area_from_diameter(diameter: float) -> float:
    return float(np.pi * (diameter / 2.0) ** 2)

def _mdot_for_target_mach(mach: float, T: float, P_ref: float, X_mole: list[float], area: float) -> float:
    rho = cb.density(T, P_ref, X_mole)
    a = cb.speed_of_sound(T, X_mole)
    return float(mach * rho * a * area)

def build_network(regime: str, mach_target: float, cfg: SweepConfig) -> FlowNetwork:
    x_mole = cb.species.dry_air()
    y_mass = cb.mole_to_mass(x_mole)
    area_hot = _area_from_diameter(cfg.d_pipe_hot)
    area_cold = _area_from_diameter(cfg.d_pipe_cold)

    mdot_hot = _mdot_for_target_mach(mach_target, cfg.t_hot_in, cfg.p_out_hot, x_mole, area_hot)
    mdot_cold = 0.65 * _mdot_for_target_mach(mach_target, cfg.t_cold_in, cfg.p_out_cold, x_mole, area_cold)

    net = FlowNetwork()
    net.add_node(MassFlowBoundary("hot_inlet", m_dot=mdot_hot, T_total=cfg.t_hot_in, Y=y_mass))
    net.add_node(PlenumNode("hot_plenum"))
    net.add_node(PressureBoundary("hot_outlet", P_total=cfg.p_out_hot, T_total=cfg.t_hot_in, Y=y_mass))

    net.add_node(MassFlowBoundary("cold_inlet", m_dot=mdot_cold, T_total=cfg.t_cold_in, Y=y_mass))
    net.add_node(PlenumNode("cold_plenum"))
    net.add_node(PressureBoundary("cold_outlet", P_total=cfg.p_out_cold, T_total=cfg.t_cold_in, Y=y_mass))

    pipe_regime = "compressible" if regime == "compressible" else "incompressible"
    orifice_regime = "compressible" if regime == "compressible" else "incompressible"

    net.add_element(PipeElement(id="hot_pipe", from_node="hot_inlet", to_node="hot_plenum", length=cfg.l_pipe, diameter=cfg.d_pipe_hot, roughness=cfg.roughness, regime=pipe_regime, surface=ConvectiveSurface(area=np.pi * cfg.d_pipe_hot * cfg.l_pipe, model=SmoothModel())))
    net.add_element(OrificeElement(id="hot_orifice", from_node="hot_plenum", to_node="hot_outlet", Cd=cfg.cd, area=_area_from_diameter(cfg.d_orifice_hot), regime=orifice_regime))
    net.add_element(PipeElement(id="cold_pipe", from_node="cold_inlet", to_node="cold_plenum", length=cfg.l_pipe, diameter=cfg.d_pipe_cold, roughness=cfg.roughness, regime=pipe_regime, surface=ConvectiveSurface(area=np.pi * cfg.d_pipe_cold * cfg.l_pipe, model=SmoothModel())))
    net.add_element(OrificeElement(id="cold_orifice", from_node="cold_plenum", to_node="cold_outlet", Cd=cfg.cd, area=_area_from_diameter(cfg.d_orifice_cold), regime=orifice_regime))
    net.add_wall(WallConnection(id="coupling_wall", element_a="hot_pipe", element_b="cold_pipe", wall_thickness=cfg.wall_thickness, wall_conductivity=cfg.wall_k))
    net.thermal_coupling_enabled = True
    return net

def run_benchmark(strategy: str, mach_values: np.ndarray, cfg: SweepConfig):
    print(f"\nRunning benchmark with strategy: {strategy}")
    total_time = 0.0
    success_count = 0
    x_prev = None

    times = []

    for m in mach_values:
        net = build_network("compressible", m, cfg)
        solver = NetworkSolver(net)

        start = time.perf_counter()

        if strategy == "continuation" and x_prev is not None:
            sol = solver.solve(method="hybr", x0=x_prev)
        else:
            sol = solver.solve(method="hybr", x0=None)

        end = time.perf_counter()
        dt = end - start

        if sol["__success__"]:
            success_count += 1
            x_prev = solver._last_x
            times.append(dt)
            total_time += dt
            # print(f"M={m:.2f}: SUCCESS in {dt:.4f}s")
        else:
            # print(f"M={m:.2f}: FAILED")
            x_prev = None # Reset on failure

    avg_time = np.mean(times) if times else 0.0
    return {
        "strategy": strategy,
        "total_time": total_time,
        "avg_time": avg_time,
        "success_rate": success_count / len(mach_values)
    }

if __name__ == "__main__":
    cfg = SweepConfig()
    # High resolution sweep for timing
    mach_values = np.linspace(0.1, 1.0, 50)

    results_cold = run_benchmark("cold_start", mach_values, cfg)
    results_cont = run_benchmark("continuation", mach_values, cfg)

    print("\n--- BENCHMARK RESULTS ---")
    print(f"{'Strategy':<15} | {'Success Rate':<15} | {'Total Time (s)':<15} | {'Avg Time/Point (s)':<20}")
    print("-" * 75)
    for r in [results_cold, results_cont]:
        print(f"{r['strategy']:<15} | {r['success_rate']:<15.2%} | {r['total_time']:<15.4f} | {r['avg_time']:<20.4f}")

    speedup = results_cold['total_time'] / results_cont['total_time'] if results_cont['total_time'] > 0 else 0
    print(f"\nContinuation Speedup: {speedup:.2f}x")
