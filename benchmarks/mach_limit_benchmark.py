#!/usr/bin/env python3
"""Benchmark: Solver convergence limit for high-Mach simple networks.

Network: MassFlowBoundary -> PipeElement -> OrificeElement -> PressureBoundary
Includes comparison of Cold Start vs Continuation initialization.
"""

import time

import numpy as np

import combaero as cb
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
)


def _area(diameter: float) -> float:
    return float(np.pi * (diameter / 2.0) ** 2)


def _mdot_for_mach(mach: float, T: float, P_ref: float, X: list[float], area: float) -> float:
    rho = cb.density(T, P_ref, X)
    a = cb.speed_of_sound(T, X)
    return float(mach * rho * a * area)


def build_network(regime: str, mach_target: float) -> FlowNetwork:
    x_air = cb.species.dry_air()
    y_air = cb.mole_to_mass(x_air)

    t_in = 300.0
    p_out = 101325.0
    d_pipe = 0.2  # Larger diameter for better supersonic stability
    length = 1.0

    mdot = _mdot_for_mach(mach_target, t_in, p_out, x_air, _area(d_pipe))

    net = FlowNetwork()
    net.add_node(MassFlowBoundary("inlet", m_dot=mdot, T_total=t_in, Y=y_air))
    net.add_node(PlenumNode("plenum"))
    net.add_node(PressureBoundary("outlet", P_total=p_out, T_total=t_in, Y=y_air))

    pipe_regime = "compressible" if regime == "compressible" else "incompressible"
    orifice_regime = "compressible" if regime == "compressible" else "incompressible"

    net.add_element(
        PipeElement(
            "pipe",
            "inlet",
            "plenum",
            length=length,
            diameter=d_pipe,
            roughness=1e-5,
            regime=pipe_regime,
        )
    )
    net.add_element(
        OrificeElement(
            "orifice", "plenum", "outlet", Cd=0.6, area=_area(d_pipe * 0.8), regime=orifice_regime
        )
    )

    return net


def run_limit_sweep(regime: str, strategy: str = "auto"):
    # Resolution for the sweep
    mach_values = np.linspace(0.1, 5.0, 50)
    results = []
    x_prev = None

    print(f"\nStarting {regime} sweep (Strategy: {strategy})...")
    total_start = time.perf_counter()

    for m in mach_values:
        net = build_network(regime, m)
        solver = NetworkSolver(net)

        # Determine x0 strategy
        x0_use = None
        if strategy == "continuation" and x_prev is not None:
            x0_use = x_prev

        start = time.perf_counter()
        sol = solver.solve(method="hybr", x0=x0_use)

        # If auto-mode and failed, try warm-start fallback
        if strategy == "auto" and not sol["__success__"] and x_prev is not None:
            sol = solver.solve(method="hybr", x0=x_prev)

        dt = time.perf_counter() - start

        # Strict success check (norm must be low)
        is_converged = sol["__success__"] and (sol["__final_norm__"] < 1e-3)

        if not is_converged:
            print(f"M={m:.2f}: FAILED. Stopping sweep.")
            break

        report = solver.get_diagnostic_report()
        stats = report["overall_statistics"]

        results.append(
            {
                "mach": m,
                "success": True,
                "dt": dt,
                "norm": sol["__final_norm__"],
                "max_rel_diff_x0": stats["max_rel_difference"],
                "p_plenum": sol["plenum.P"],
            }
        )
        x_prev = solver._last_x
        # print(f"M={m:.2f}: SUCCESS ({dt:.4f}s)")

    total_time = time.perf_counter() - total_start
    return results, total_time


if __name__ == "__main__":
    # Compare strategies on compressible regime
    res_cold, time_cold = run_limit_sweep("compressible", strategy="cold")
    res_cont, time_cont = run_limit_sweep("compressible", strategy="continuation")

    print("\n--- STRATEGY COMPARISON (Compressible) ---")
    print(
        f"{'Strategy':<15} | {'Max Mach':<10} | {'Total Time (s)':<15} | {'Avg Time/Point (s)':<20}"
    )
    print("-" * 65)
    print(
        f"{'Cold Start':<15} | {res_cold[-1]['mach'] if res_cold else 0:<10.2f} | {time_cold:<15.4f} | {time_cold / len(res_cold) if res_cold else 0:<20.4f}"
    )
    print(
        f"{'Continuation':<15} | {res_cont[-1]['mach'] if res_cont else 0:<10.2f} | {time_cont:<15.4f} | {time_cont / len(res_cont) if res_cont else 0:<20.4f}"
    )

    if time_cont > 0:
        print(f"\nContinuation Speedup: {time_cold / time_cont:.2f}x")

    # Generate Markdown Report
    with open("benchmarks/mach_limit_report.md", "w") as f:
        f.write("# Mach Limit Benchmark Report\n\n")
        f.write("## Strategy Comparison\n\n")
        f.write("| Strategy | Max Mach | Total Time [s] | Avg Time/Point [s] |\n")
        f.write("| :--- | :--- | :--- | :--- |\n")
        f.write(
            f"| Cold Start | {res_cold[-1]['mach'] if res_cold else 0:.2f} | {time_cold:.4f} | {time_cold / len(res_cold) if res_cold else 0:.4f} |\n"
        )
        f.write(
            f"| Continuation | {res_cont[-1]['mach'] if res_cont else 0:.2f} | {time_cont:.4f} | {time_cont / len(res_cont) if res_cont else 0:.4f} |\n\n"
        )

        f.write("## Detail: Continuation Data\n\n")
        f.write("| Mach | dt [s] | norm | max_x0_diff | P_plenum |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- |\n")
        for r in res_cont:
            f.write(
                f"| {r['mach']:.2f} | {r['dt']:.4f} | {r['norm']:.2e} | {r['max_rel_diff_x0']:.2e} | {r['p_plenum']:.1f} |\n"
            )

    print("\nReport updated: benchmarks/mach_limit_report.md")
