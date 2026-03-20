#!/usr/bin/env python3
"""Standalone network solver stability and timing benchmarks.

This script is intentionally outside the regular pytest suite. It provides
repeatable performance/stability snapshots for representative network classes:

1. Simple incompressible network
2. Fully coupled heat-transfer network
3. Compressible flow network
4. Combustion network
5. Fully coupled combustion path (hot gas -> cooling -> combustor inlet)
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, TimeoutError
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, median
from time import perf_counter
from typing import Any

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
from combaero.network.components import (
    CombustorNode,
    EffectiveAreaConnectionElement,
)


def _area(diameter: float) -> float:
    return float(math.pi * (diameter / 2.0) ** 2)


def _mdot_for_mach(mach: float, temperature: float, pressure: float, x: list[float], area: float) -> float:
    rho = cb.density(temperature, pressure, x)
    a = cb.speed_of_sound(temperature, x)
    return float(mach * rho * a * area)


def _build_simple_network() -> FlowNetwork:
    net = FlowNetwork()
    inlet = PressureBoundary("inlet", P_total=2.8e5, T_total=420.0)
    plenum = PlenumNode("plenum")
    outlet = PressureBoundary("outlet", P_total=1.35e5, T_total=380.0)

    net.add_node(inlet)
    net.add_node(plenum)
    net.add_node(outlet)

    net.add_element(OrificeElement("orifice", "inlet", "plenum", Cd=0.67, area=1.6e-4))
    net.add_element(
        PipeElement(
            "pipe",
            "plenum",
            "outlet",
            length=1.8,
            diameter=0.028,
            roughness=8.0e-5,
            regime="incompressible",
        )
    )
    return net


def _build_fully_coupled_heat_network() -> FlowNetwork:
    x_air = cb.standard_dry_air_composition()
    y_air = list(cb.mole_to_mass(x_air))

    t_hot = 920.0
    t_cold = 360.0
    p_out_hot = 2.3e5
    p_out_cold = 2.1e5
    d_hot = 0.020
    d_cold = 0.016
    length = 2.4

    mdot_hot = _mdot_for_mach(0.55, t_hot, p_out_hot, x_air, _area(d_hot))
    mdot_cold = 0.70 * _mdot_for_mach(0.45, t_cold, p_out_cold, x_air, _area(d_cold))

    net = FlowNetwork()
    net.add_node(MassFlowBoundary("hot_inlet", m_dot=mdot_hot, T_total=t_hot, Y=y_air))
    net.add_node(PlenumNode("hot_plenum"))
    net.add_node(PressureBoundary("hot_outlet", P_total=p_out_hot, T_total=t_hot, Y=y_air))

    net.add_node(MassFlowBoundary("cold_inlet", m_dot=mdot_cold, T_total=t_cold, Y=y_air))
    net.add_node(PlenumNode("cold_plenum"))
    net.add_node(PressureBoundary("cold_outlet", P_total=p_out_cold, T_total=t_cold, Y=y_air))

    hot_surface = ConvectiveSurface(area=math.pi * d_hot * length, model=SmoothModel())
    cold_surface = ConvectiveSurface(area=math.pi * d_cold * length, model=SmoothModel())

    net.add_element(
        PipeElement(
            id="hot_pipe",
            from_node="hot_inlet",
            to_node="hot_plenum",
            length=length,
            diameter=d_hot,
            roughness=5.0e-5,
            regime="incompressible",
            surface=hot_surface,
        )
    )
    net.add_element(
        OrificeElement(
            id="hot_orifice",
            from_node="hot_plenum",
            to_node="hot_outlet",
            Cd=0.72,
            area=_area(0.018),
            regime="incompressible",
        )
    )

    net.add_element(
        PipeElement(
            id="cold_pipe",
            from_node="cold_inlet",
            to_node="cold_plenum",
            length=length,
            diameter=d_cold,
            roughness=5.0e-5,
            regime="incompressible",
            surface=cold_surface,
        )
    )
    net.add_element(
        OrificeElement(
            id="cold_orifice",
            from_node="cold_plenum",
            to_node="cold_outlet",
            Cd=0.72,
            area=_area(0.015),
            regime="incompressible",
        )
    )

    net.add_wall(
        WallConnection(
            id="coupling_wall",
            element_a="hot_pipe",
            element_b="cold_pipe",
            wall_thickness=0.0025,
            wall_conductivity=19.0,
        )
    )
    net.thermal_coupling_enabled = True
    return net


def _build_compressible_network() -> FlowNetwork:
    net = FlowNetwork()
    inlet = PressureBoundary("inlet", P_total=4.8e5, T_total=540.0)
    mid = PlenumNode("mid")
    outlet = PressureBoundary("outlet", P_total=1.9e5, T_total=520.0)

    net.add_node(inlet)
    net.add_node(mid)
    net.add_node(outlet)

    net.add_element(
        PipeElement(
            "fanno_pipe",
            "inlet",
            "mid",
            length=3.2,
            diameter=0.022,
            roughness=8.0e-5,
            regime="compressible_fanno",
            friction_model="haaland",
        )
    )
    net.add_element(
        OrificeElement(
            "comp_orifice",
            "mid",
            "outlet",
            Cd=0.69,
            area=1.2e-4,
            regime="compressible",
        )
    )
    return net


def _build_combustion_network() -> FlowNetwork:
    net = FlowNetwork()

    x_air = cb.standard_dry_air_composition()
    y_air = cb.mole_to_mass(x_air)

    x_ch4 = [0.0] * cb.num_species()
    x_ch4[cb.species_index_from_name("CH4")] = 1.0
    y_ch4 = cb.mole_to_mass(x_ch4)

    air = MassFlowBoundary("air", m_dot=2.3, T_total=760.0, Y=y_air)
    fuel = MassFlowBoundary("fuel", m_dot=0.068, T_total=340.0, Y=y_ch4)

    combustor = CombustorNode("comb", method="equilibrium")
    outlet = PressureBoundary("outlet", P_total=1.3e5, T_total=700.0)

    net.add_node(air)
    net.add_node(fuel)
    net.add_node(combustor)
    net.add_node(outlet)

    e_air = EffectiveAreaConnectionElement("e_air", "air", "comb", effective_area=0.06)
    e_fuel = EffectiveAreaConnectionElement("e_fuel", "fuel", "comb", effective_area=0.01)
    nozzle = OrificeElement("nozzle", "comb", "outlet", Cd=0.8, area=0.018)

    net.add_element(e_air)
    net.add_element(e_fuel)
    net.add_element(nozzle)
    return net


def _build_fully_coupled_combustion_network() -> FlowNetwork:
    """Build a fully coupled combustion network for benchmarking.

    This creates a realistic combustion network with proper coupling
    that uses NetworkSolver.solve() directly for proper benchmarking.
    The network represents a challenging but solvable combustion case.
    """
    # Species compositions
    x_air = cb.standard_dry_air_composition()
    y_air = list(cb.mole_to_mass(x_air))
    x_ch4 = [0.0] * cb.num_species()
    x_ch4[cb.species_index_from_name("CH4")] = 1.0
    y_ch4 = list(cb.mole_to_mass(x_ch4))

    # Network setup - realistic combustion parameters
    net = FlowNetwork()

    # Add boundary nodes with realistic conditions
    air = MassFlowBoundary("air_inlet", m_dot=1.0, T_total=650.0, Y=y_air)
    fuel = MassFlowBoundary("fuel_inlet", m_dot=0.04, T_total=300.0, Y=y_ch4)

    # Use proper CombustorNode for combustion physics
    combustor = CombustorNode("combustor", method="equilibrium")

    net.add_node(PlenumNode("turbine_inlet"))
    net.add_node(PressureBoundary("outlet", P_total=1.0e5, T_total=1500.0, Y=y_air))

    net.add_node(air)
    net.add_node(fuel)
    net.add_node(combustor)

    # Add realistic elements with proper coupling
    # Air inlet restriction
    net.add_element(OrificeElement("air_restrictor", "air_inlet", "combustor",
                                   Cd=0.7, area=0.02, regime="incompressible"))

    # Fuel injection
    net.add_element(OrificeElement("fuel_injector", "fuel_inlet", "combustor",
                                   Cd=0.6, area=0.005, regime="incompressible"))

    # Combustor pressure drop (significant)
    net.add_element(PipeElement("combustor", "combustor", "turbine_inlet",
                                length=1.5, diameter=0.08, roughness=2e-5, regime="incompressible"))

    # Turbine inlet restriction
    net.add_element(OrificeElement("turbine_inlet", "turbine_inlet", "outlet",
                                   Cd=0.8, area=0.025, regime="incompressible"))

    return net


def _run_network_case(case_name: str, solver_timeout_s: float) -> dict[str, Any]:
    builders: dict[str, Any] = {
        "simple_network": _build_simple_network,
        "fully_coupled_heat_transfer": _build_fully_coupled_heat_network,
        "compressible_flow_network": _build_compressible_network,
        "combustion_network": _build_combustion_network,
        "fully_coupled_combustion_network": _build_fully_coupled_combustion_network,
    }

    net = builders[case_name]()
    solver = NetworkSolver(net)

    t0 = perf_counter()
    sol = solver.solve(method="hybr", timeout=solver_timeout_s)
    elapsed = perf_counter() - t0

    return {
        "success": bool(sol.get("__success__", False)),
        "elapsed_s": elapsed,
        "final_norm": float(sol.get("__final_norm__", float("nan"))),
        "message": str(sol.get("__message__", "")),
    }


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    timeout_s: float


def _run_case_once(case: BenchmarkCase, solver_timeout_s: float) -> dict[str, Any]:
    with ProcessPoolExecutor(max_workers=1, mp_context=mp.get_context("spawn")) as pool:
        fut = pool.submit(_run_network_case, case.name, solver_timeout_s)
        return fut.result(timeout=case.timeout_s)


def run_benchmarks(repeats: int, case_timeout_s: float, solver_timeout_s: float) -> dict[str, Any]:
    cases = [
        BenchmarkCase("simple_network", case_timeout_s),
        BenchmarkCase("fully_coupled_heat_transfer", case_timeout_s),
        BenchmarkCase("compressible_flow_network", case_timeout_s),
        BenchmarkCase("combustion_network", case_timeout_s),
        BenchmarkCase("fully_coupled_combustion_network", case_timeout_s),
    ]

    run_stamp = datetime.now(timezone.utc).isoformat()
    results: dict[str, Any] = {
        "run_utc": run_stamp,
        "repeats": repeats,
        "case_timeout_s": case_timeout_s,
        "solver_timeout_s": solver_timeout_s,
        "cases": {},
    }

    for case in cases:
        times: list[float] = []
        success_count = 0
        timeout_count = 0
        error_messages: list[str] = []
        final_norms: list[float] = []

        for _ in range(repeats):
            try:
                res = _run_case_once(case, solver_timeout_s)
                times.append(float(res["elapsed_s"]))
                if bool(res["success"]):
                    success_count += 1
                if math.isfinite(float(res["final_norm"])):
                    final_norms.append(float(res["final_norm"]))
                message = str(res.get("message", ""))
                if message:
                    error_messages.append(message)
            except TimeoutError:
                timeout_count += 1
                error_messages.append(f"case timeout > {case.timeout_s}s")
            except Exception as exc:  # noqa: BLE001
                error_messages.append(str(exc))

        stats: dict[str, Any] = {
            "success_count": success_count,
            "attempts": repeats,
            "timeout_count": timeout_count,
            "success_rate": success_count / repeats,
            "messages": error_messages,
        }

        if times:
            stats.update(
                {
                    "mean_s": mean(times),
                    "median_s": median(times),
                    "min_s": min(times),
                    "max_s": max(times),
                }
            )
        if final_norms:
            stats["median_final_norm"] = median(final_norms)
            stats["max_final_norm"] = max(final_norms)

        results["cases"][case.name] = stats

    return results


def format_markdown(results: dict[str, Any]) -> str:
    lines = [
        "# Network Solver Benchmarks",
        "",
        f"- Run UTC: `{results['run_utc']}`",
        f"- Repeats per case: `{results['repeats']}`",
        f"- Case timeout: `{results['case_timeout_s']:.1f}s`",
        f"- Solver timeout: `{results['solver_timeout_s']:.1f}s`",
        "",
        "| Case | Success | Time median [s] | Time max [s] | Median |F| |",
        "|---|---:|---:|---:|---:|",
    ]

    for case_name, stats in results["cases"].items():
        success_txt = f"{stats['success_count']}/{stats['attempts']}"
        median_s = stats.get("median_s", float("nan"))
        max_s = stats.get("max_s", float("nan"))
        median_norm = stats.get("median_final_norm", float("nan"))
        lines.append(
            "| "
            f"`{case_name}` | {success_txt} | {median_s:.4f} | {max_s:.4f} | {median_norm:.3e} |"
        )

    lines.append("")
    lines.append("## Notes")
    lines.append("- This benchmark is informational and not part of default CI gating.")
    lines.append("- Timeouts are used to prevent hangs in difficult nonlinear regimes.")
    return "\n".join(lines)


def append_history(markdown_path: Path, results: dict[str, Any]) -> None:
    section = [
        f"## {results['run_utc']}",
        "",
        format_markdown(results),
        "",
    ]

    if markdown_path.exists():
        existing = markdown_path.read_text(encoding="utf-8")
    else:
        existing = "# Network Solver Benchmark History\n\n"

    markdown_path.write_text(existing + "\n" + "\n".join(section), encoding="utf-8")


def _timestamp_slug(run_utc: str) -> str:
    """Convert ISO UTC stamp to filename-safe UTC timestamp."""
    dt = datetime.fromisoformat(run_utc.replace("Z", "+00:00"))
    return dt.strftime("%Y%m%d_%H%M%SZ")


def write_json_outputs(
    *,
    results: dict[str, Any],
    latest_json_path: Path,
    archive_dir: Path,
    overwrite_latest: bool,
) -> Path:
    """Write benchmark JSON to archive and optionally latest path.

    Returns the path of the archived JSON file.
    """
    archive_dir.mkdir(parents=True, exist_ok=True)
    run_slug = _timestamp_slug(str(results["run_utc"]))
    archive_path = archive_dir / f"network_solver_benchmark_{run_slug}.json"
    archive_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    latest_json_path.parent.mkdir(parents=True, exist_ok=True)
    if overwrite_latest or not latest_json_path.exists():
        latest_json_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    return archive_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run network solver stability/timing benchmarks.")
    parser.add_argument("--repeats", type=int, default=3, help="Number of runs per benchmark case.")
    parser.add_argument(
        "--case-timeout",
        type=float,
        default=90.0,
        help="Hard timeout (seconds) for each benchmark case process.",
    )
    parser.add_argument(
        "--solver-timeout",
        type=float,
        default=25.0,
        help="Timeout passed to NetworkSolver.solve() in seconds.",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        default=Path(__file__).parent.parent / "benchmarks" / "network_solver_benchmark_latest.json",
        help="Path for latest benchmark JSON output.",
    )
    parser.add_argument(
        "--history-md",
        type=Path,
        default=Path(__file__).parent.parent / "benchmarks" / "NETWORK_SOLVER_BENCHMARK_HISTORY.md",
        help="Path for benchmark markdown history.",
    )
    parser.add_argument(
        "--json-archive-dir",
        type=Path,
        default=Path(__file__).parent.parent / "benchmarks" / "runs",
        help="Directory for timestamped benchmark JSON snapshots.",
    )
    parser.add_argument(
        "--overwrite-latest-json",
        action="store_true",
        help="Overwrite --json-out with current run (default keeps existing file for comparison).",
    )

    args = parser.parse_args()
    results = run_benchmarks(
        repeats=args.repeats,
        case_timeout_s=args.case_timeout,
        solver_timeout_s=args.solver_timeout,
    )

    archive_path = write_json_outputs(
        results=results,
        latest_json_path=args.json_out,
        archive_dir=args.json_archive_dir,
        overwrite_latest=args.overwrite_latest_json,
    )

    append_history(args.history_md, results)

    print(format_markdown(results))
    print(f"\nWrote archived JSON: {archive_path}")
    if args.overwrite_latest_json or not args.json_out.exists():
        print(f"Updated latest JSON: {args.json_out}")
    else:
        print(f"Kept existing latest JSON (unchanged): {args.json_out}")
    print(f"Updated history: {args.history_md}")


if __name__ == "__main__":
    main()
