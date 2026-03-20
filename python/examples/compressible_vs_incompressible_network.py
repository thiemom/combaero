#!/usr/bin/env python3
"""Fully coupled compressible vs incompressible pressure-ratio sweep.

This example solves a wall-coupled two-stream network for a Mach sweep
and compares the pressure ratio required by:

- compressible elements (Fanno pipe + compressible orifice)
- incompressible elements (Darcy/Bernoulli variants)

Both models use the same fully coupled topology (hot/cold streams linked
by ``WallConnection`` with ``thermal_coupling_enabled=True``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from plot_utils import show_or_save

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
    d_orifice_hot: float = 0.030
    d_orifice_cold: float = 0.025
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
    x_mole = cb.standard_dry_air_composition()
    y_mass = list(cb.mole_to_mass(x_mole))

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

    pipe_regime = "compressible_fanno" if regime == "compressible" else "incompressible"
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
    pr_hot: list[float] = []
    pr_cold: list[float] = []
    mach_solved_hot: list[float] = []
    delta_t_hot: list[float] = []
    delta_t_cold: list[float] = []

    x_mole = cb.standard_dry_air_composition()
    area_hot = _area_from_diameter(cfg.d_pipe_hot)
    x0_prev: np.ndarray | None = None

    for m_target in mach_values:
        net, mdot_hot = _build_coupled_network(regime=regime, mach_target=float(m_target), cfg=cfg)
        solver = NetworkSolver(net)

        # Try fresh solve; fall back to warm-start from the previous point.
        sol = solver.solve(method="hybr")
        if not sol.get("__success__", False) and x0_prev is not None:
            sol = solver.solve(method="hybr", x0=x0_prev)
        if not sol.get("__success__", False):
            raise RuntimeError(
                f"{regime} solve failed at target Mach={m_target:.3f}: "
                f"|F|={sol.get('__final_norm__', '?'):.3e}  "
                f"{sol.get('__message__', 'unknown')}"
            )
        x0_prev = solver._last_x

        p_in_hot = float(sol["hot_inlet.P_total"])
        p_in_cold = float(sol["cold_inlet.P_total"])
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
        "mach_target": np.asarray(mach_values, dtype=float),
        "mach_solved_hot": np.asarray(mach_solved_hot, dtype=float),
        "pr_hot": np.asarray(pr_hot, dtype=float),
        "pr_cold": np.asarray(pr_cold, dtype=float),
        "delta_t_hot": np.asarray(delta_t_hot, dtype=float),
        "delta_t_cold": np.asarray(delta_t_cold, dtype=float),
    }


def main() -> None:
    cfg = SweepConfig()
    mach_values = np.linspace(0.05, 1.0, 25)

    data_incomp = run_sweep(mach_values, regime="incompressible", cfg=cfg)
    data_comp = run_sweep(mach_values, regime="compressible", cfg=cfg)

    pr_rel_diff_pct = 100.0 * (data_incomp["pr_hot"] - data_comp["pr_hot"]) / data_comp["pr_hot"]

    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10, 10), sharex=True)

    ax[0].plot(
        data_comp["mach_solved_hot"], data_comp["pr_hot"], "-", label="PR hot (compressible)"
    )
    ax[0].plot(
        data_incomp["mach_solved_hot"],
        data_incomp["pr_hot"],
        "--",
        label="PR hot (incompressible)",
    )
    ax[0].set_ylabel("P_out / P_in")
    ax[0].set_title("Fully coupled network: pressure-ratio sweep")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()

    ax[1].plot(data_comp["mach_solved_hot"], pr_rel_diff_pct, "-")
    ax[1].axhline(0.0, color="k", linestyle="--", linewidth=1.0)
    ax[1].set_ylabel("PR difference [%]\n(incomp - comp)")
    ax[1].grid(True, alpha=0.3)

    ax[2].plot(data_comp["mach_solved_hot"], data_comp["delta_t_hot"], "-", label="Hot cooling ΔT")
    ax[2].plot(
        data_comp["mach_solved_hot"], data_comp["delta_t_cold"], "-", label="Cold heating ΔT"
    )
    ax[2].set_xlabel("Mach number (hot stream)")
    ax[2].set_ylabel("Thermal coupling ΔT [K]")
    ax[2].grid(True, alpha=0.3)
    ax[2].legend()

    fig.tight_layout()
    show_or_save(fig, "compressible_vs_incompressible_coupled_sweep.png")

    print("Rendered: compressible_vs_incompressible_coupled_sweep.png")
    print(
        "Hot-stream PR difference range [%]: "
        f"{float(np.min(pr_rel_diff_pct)):.2f} .. {float(np.max(pr_rel_diff_pct)):.2f}"
    )


if __name__ == "__main__":
    main()
