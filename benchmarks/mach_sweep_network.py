"""Benchmark: Fully coupled Mach sweep (0.05 to 1.0).

This script performs an intensive Mach sweep calculation for both
compressible and incompressible coupled network models.
"""

import os
import sys

from dataclasses import dataclass
from typing import Literal

import matplotlib.pyplot as plt
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

# Reuse logic from example but with full resolution
Regime = Literal["compressible", "incompressible"]


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


def _mdot_for_target_mach(
    mach: float, T: float, P_ref: float, X_mole: list[float], area: float
) -> float:
    rho = cb.density(T, P_ref, X_mole)
    a = cb.speed_of_sound(T, X_mole)
    return float(mach * rho * a * area)


def _build_coupled_network(
    regime: Regime, mach_target: float, cfg: SweepConfig
) -> tuple[FlowNetwork, float]:
    x_mole = cb.species.dry_air()
    y_mass = cb.mole_to_mass(x_mole)
    area_hot = _area_from_diameter(cfg.d_pipe_hot)
    area_cold = _area_from_diameter(cfg.d_pipe_cold)
    mdot_hot = _mdot_for_target_mach(mach_target, cfg.t_hot_in, cfg.p_out_hot, x_mole, area_hot)
    mdot_cold = 0.65 * _mdot_for_target_mach(
        mach_target, cfg.t_cold_in, cfg.p_out_cold, x_mole, area_cold
    )
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("hot_inlet", m_dot=mdot_hot, T_total=cfg.t_hot_in, Y=y_mass))
    net.add_node(PlenumNode("hot_plenum"))
    net.add_node(
        PressureBoundary("hot_outlet", P_total=cfg.p_out_hot, T_total=cfg.t_hot_in, Y=y_mass)
    )
    net.add_node(MassFlowBoundary("cold_inlet", m_dot=mdot_cold, T_total=cfg.t_cold_in, Y=y_mass))
    net.add_node(PlenumNode("cold_plenum"))
    net.add_node(
        PressureBoundary("cold_outlet", P_total=cfg.p_out_cold, T_total=cfg.t_cold_in, Y=y_mass)
    )
    pipe_regime = "compressible" if regime == "compressible" else "incompressible"
    orifice_regime = "compressible" if regime == "compressible" else "incompressible"
    hot_surface = ConvectiveSurface(area=np.pi * cfg.d_pipe_hot * cfg.l_pipe, model=SmoothModel())
    cold_surface = ConvectiveSurface(area=np.pi * cfg.d_pipe_cold * cfg.l_pipe, model=SmoothModel())
    net.add_element(
        PipeElement(
            id="hot_pipe",
            from_node="hot_inlet",
            to_node="hot_plenum",
            length=cfg.l_pipe,
            diameter=cfg.d_pipe_hot,
            roughness=cfg.roughness,
            regime=pipe_regime,
            surface=hot_surface,
        )
    )
    net.add_element(
        OrificeElement(
            id="hot_orifice",
            from_node="hot_plenum",
            to_node="hot_outlet",
            Cd=cfg.cd,
            area=_area_from_diameter(cfg.d_orifice_hot),
            regime=orifice_regime,
        )
    )
    net.add_element(
        PipeElement(
            id="cold_pipe",
            from_node="cold_inlet",
            to_node="cold_plenum",
            length=cfg.l_pipe,
            diameter=cfg.d_pipe_cold,
            roughness=cfg.roughness,
            regime=pipe_regime,
            surface=cold_surface,
        )
    )
    net.add_element(
        OrificeElement(
            id="cold_orifice",
            from_node="cold_plenum",
            to_node="cold_outlet",
            Cd=cfg.cd,
            area=_area_from_diameter(cfg.d_orifice_cold),
            regime=orifice_regime,
        )
    )
    net.add_wall(
        WallConnection(
            id="coupling_wall",
            element_a="hot_pipe",
            element_b="cold_pipe",
            wall_thickness=cfg.wall_thickness,
            wall_conductivity=cfg.wall_k,
        )
    )
    net.thermal_coupling_enabled = True
    return net, mdot_hot


def run_sweep(mach_values: np.ndarray, regime: Regime, cfg: SweepConfig) -> dict[str, np.ndarray]:
    pr_hot, pr_cold, mach_solved_hot, delta_t_hot, delta_t_cold = [], [], [], [], []
    x_mole = cb.species.dry_air()
    area_hot = _area_from_diameter(cfg.d_pipe_hot)
    x0_prev = None
    for m_target in mach_values:
        net, mdot_hot = _build_coupled_network(regime=regime, mach_target=float(m_target), cfg=cfg)
        solver = NetworkSolver(net)
        sol = solver.solve(method="hybr")
        if not sol.get("__success__", False) and x0_prev is not None:
            sol = solver.solve(method="hybr", x0=x0_prev)
        if not sol.get("__success__", False):
            print(f"WARNING: {regime} solve failed at target Mach={m_target:.3f}")
            continue
        x0_prev = solver._last_x
        p_in_hot, p_in_cold = float(sol["hot_inlet.P_total"]), float(sol["cold_inlet.P_total"])
        pr_hot.append(cfg.p_out_hot / p_in_hot)
        pr_cold.append(cfg.p_out_cold / p_in_cold)
        rho_hot = cb.density(cfg.t_hot_in, p_in_hot, x_mole)
        a_hot = cb.speed_of_sound(cfg.t_hot_in, x_mole)
        mach_solved_hot.append(mdot_hot / (rho_hot * a_hot * area_hot))
        t_hot_out = float(solver._derived_states["hot_plenum"][0])
        t_cold_out = float(solver._derived_states["cold_plenum"][0])
        delta_t_hot.append(cfg.t_hot_in - t_hot_out)
        delta_t_cold.append(t_cold_out - cfg.t_cold_in)

    return {
        k: np.array(v)
        for k, v in zip(
            ["mach_solved_hot", "pr_hot", "pr_cold", "delta_t_hot", "delta_t_cold"],
            [mach_solved_hot, pr_hot, pr_cold, delta_t_hot, delta_t_cold],
            strict=True,
        )
    }


# Results processing
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "python", "examples"))
)

try:
    from plot_utils import show_or_save
except ImportError:
    # Fallback if invoked from a different location
    def show_or_save(fig, filename):
        fig.savefig(filename)
        print(f"Saved: {filename}")


if __name__ == "__main__":
    cfg = SweepConfig()
    mach_values = np.linspace(0.05, 1.0, 30)
    print(f"Running intensive Mach sweep ({len(mach_values)} points)...")

    data_comp = run_sweep(mach_values, "compressible", cfg)
    data_incomp = run_sweep(mach_values, "incompressible", cfg)

    # Plotting
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True)

    ax[0].plot(data_comp["mach_solved_hot"], data_comp["pr_hot"], "-", label="Compressible")
    ax[0].plot(data_incomp["mach_solved_hot"], data_incomp["pr_hot"], "--", label="Incompressible")
    ax[0].set_ylabel("Pressure Ratio (P_out / P_in)")
    ax[0].set_title("Benchmark: Coupled Network Mach Sweep")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()

    pr_diff = 100.0 * (data_incomp["pr_hot"] - data_comp["pr_hot"]) / data_comp["pr_hot"]
    ax[1].plot(data_comp["mach_solved_hot"], pr_diff, "r-")
    ax[1].set_ylabel("PR Difference [%]\n(Incomp - Comp)")
    ax[1].set_xlabel("Mach Number")
    ax[1].grid(True, alpha=0.3)

    fig.tight_layout()
    show_or_save(fig, "benchmark_mach_sweep.png")
    print("Done. Benchmark complete.")
