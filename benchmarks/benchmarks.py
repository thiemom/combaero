#!/usr/bin/env python3
"""Standalone network solver stability and timing benchmarks.

This script replaces the old bespoke benchmark networks with the highly
parameterized "v3" grid logic. It evaluates solver robustness and speed
across complex grids (e.g. 1x1 up to 10x10) tracking both fully
incompressible boundary assumptions against highly unstable compressible bounds.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, TimeoutError
from dataclasses import dataclass
from datetime import UTC, datetime
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
    MomentumChamberNode,
)


def _build_fully_coupled_grid(
    m_dot_air: float,
    T_ad: float,
    n_serial: int,
    n_parallel: int,
    total_area: float,
    total_ht_area: float,
    use_compressible: bool,
) -> FlowNetwork:
    """Builds a highly parameterized recuperative combustor flow network."""
    T_air = 560.0
    P_ref = 101325.0
    X_air = cb.humid_air_composition(288.16, P_ref, 0.6)
    Y_air = cb.mole_to_mass(X_air)

    Y_fuel = [0.0] * cb.num_species()
    Y_fuel[cb.species_index_from_name("H2")] = 1.0
    X_fuel = cb.mass_to_mole(Y_fuel)

    air_stream = cb.Stream()
    air_stream.set_T(T_air).set_P(P_ref).set_X(X_air).set_mdot(m_dot_air)

    fuel_stream = cb.Stream()
    fuel_stream.set_T(300.0).set_P(P_ref).set_X(X_fuel)

    fuel_stream = cb.set_fuel_stream_for_Tad(T_ad, fuel_stream, air_stream)
    m_dot_fuel = fuel_stream.mdot

    net = FlowNetwork()

    inlet_air = MassFlowBoundary("inlet_air", m_dot=m_dot_air, T_total=T_air, Y=Y_air)
    inlet_fuel = MassFlowBoundary("inlet_fuel", m_dot=m_dot_fuel, T_total=300.0, Y=Y_fuel)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0, Y=Y_air)

    net.add_node(inlet_air)
    net.add_node(inlet_fuel)
    net.add_node(outlet)

    distributor = PlenumNode("distributor")
    collector = PlenumNode("collector")
    net.add_node(distributor)
    net.add_node(collector)

    pipe_regime = "compressible" if use_compressible else "incompressible"
    ori_regime = "compressible" if use_compressible else "incompressible"

    net.add_element(
        OrificeElement(
            "ori_air_in", "inlet_air", "distributor", Cd=0.8, area=total_area, regime=ori_regime
        )
    )

    A_branch = total_area / n_parallel
    D_branch = math.sqrt(4 * A_branch / math.pi)

    A_ht_pipe = (total_ht_area / n_parallel) / max(1, n_serial)
    L_pipe = 1.0

    distributor_fuel = PlenumNode("distributor_fuel")
    net.add_node(distributor_fuel)
    net.add_element(
        OrificeElement(
            "ori_fuel_in",
            "inlet_fuel",
            "distributor_fuel",
            Cd=0.6,
            area=A_branch * 0.1 * n_parallel,
            regime=ori_regime,
        )
    )

    for p in range(n_parallel):
        mixing = PlenumNode(f"mixing_plenum_{p}")
        combustor = CombustorNode(f"combustor_{p}", method="complete")
        net.add_node(mixing)
        net.add_node(combustor)

        net.add_element(
            OrificeElement(
                f"ori_fuel_{p}",
                "distributor_fuel",
                f"mixing_plenum_{p}",
                Cd=0.6,
                area=A_branch * 0.1,
                regime=ori_regime,
            )
        )

        for i in range(n_serial + 1):
            net.add_node(MomentumChamberNode(f"cold_chamber_{p}_{i}", area=A_branch))
            net.add_node(MomentumChamberNode(f"hot_chamber_{p}_{i}", area=A_branch))

        net.add_element(
            OrificeElement(
                f"ori_cold_in_{p}",
                "distributor",
                f"cold_chamber_{p}_0",
                Cd=0.8,
                area=A_branch,
                regime=ori_regime,
            )
        )
        for i in range(n_serial):
            net.add_element(
                PipeElement(
                    f"cold_pipe_{p}_{i}",
                    f"cold_chamber_{p}_{i}",
                    f"cold_chamber_{p}_{i + 1}",
                    length=L_pipe,
                    diameter=D_branch,
                    roughness=2e-5,
                    regime=pipe_regime,
                    surface=ConvectiveSurface(area=A_ht_pipe, model=SmoothModel())
                    if A_ht_pipe > 0
                    else None,
                )
            )

        net.add_element(
            OrificeElement(
                f"ori_comb_in_{p}",
                f"cold_chamber_{p}_{n_serial}",
                f"mixing_plenum_{p}",
                Cd=0.8,
                area=A_branch,
                regime=ori_regime,
            )
        )
        net.add_element(
            OrificeElement(
                f"ori_mix_to_comb_{p}",
                f"mixing_plenum_{p}",
                f"combustor_{p}",
                Cd=0.8,
                area=A_branch,
                regime=ori_regime,
            )
        )
        net.add_element(
            OrificeElement(
                f"ori_comb_out_{p}",
                f"combustor_{p}",
                f"hot_chamber_{p}_0",
                Cd=0.8,
                area=A_branch,
                regime=ori_regime,
            )
        )

        for i in range(n_serial):
            net.add_element(
                PipeElement(
                    f"hot_pipe_{p}_{i}",
                    f"hot_chamber_{p}_{i}",
                    f"hot_chamber_{p}_{i + 1}",
                    length=L_pipe,
                    diameter=D_branch,
                    roughness=2e-5,
                    regime=pipe_regime,
                    surface=ConvectiveSurface(area=A_ht_pipe, model=SmoothModel())
                    if A_ht_pipe > 0
                    else None,
                )
            )

            if A_ht_pipe > 0:
                wall_cold_pipe_idx = n_serial - 1 - i
                net.add_wall(
                    WallConnection(
                        id=f"wall_{p}_{i}",
                        element_a=f"hot_pipe_{p}_{i}",
                        element_b=f"cold_pipe_{p}_{wall_cold_pipe_idx}",
                        wall_thickness=0.005,
                        wall_conductivity=20.0,
                        contact_area=A_ht_pipe,
                    )
                )

        net.add_element(
            OrificeElement(
                f"ori_hot_out_{p}",
                f"hot_chamber_{p}_{n_serial}",
                "collector",
                Cd=0.8,
                area=A_branch,
                regime=ori_regime,
            )
        )

    net.add_element(
        OrificeElement(
            "ori_exhaust", "collector", "outlet", Cd=0.8, area=total_area, regime=ori_regime
        )
    )

    net.thermal_coupling_enabled = True
    return net


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    timeout_s: float
    m_dot_air: float
    T_ad: float
    n_serial: int
    n_parallel: int
    total_area: float
    total_ht_area: float
    use_compressible: bool


def _run_network_case(case: BenchmarkCase, solver_timeout_s: float) -> dict[str, Any]:
    net = _build_fully_coupled_grid(
        m_dot_air=case.m_dot_air,
        T_ad=case.T_ad,
        n_serial=case.n_serial,
        n_parallel=case.n_parallel,
        total_area=case.total_area,
        total_ht_area=case.total_ht_area,
        use_compressible=case.use_compressible,
    )

    solver = NetworkSolver(net)

    solver_options = {
        "maxfev": 10000,
        "xtol": 1e-4,
        "ftol": 1e-6,
    }

    t0 = perf_counter()
    sol = solver.solve(
        method="lm",
        timeout=solver_timeout_s,
        options=solver_options,
        init_strategy="incompressible_warmstart",
    )
    elapsed = perf_counter() - t0

    return {
        "success": bool(sol.get("__success__", False)),
        "elapsed_s": elapsed,
        "final_norm": float(sol.get("__final_norm__", float("nan"))),
        "message": str(sol.get("__message__", "")),
        "iterations": int(sol.get("__nfev__", 0)),
    }


def _run_case_once(case: BenchmarkCase, solver_timeout_s: float) -> dict[str, Any]:
    with ProcessPoolExecutor(max_workers=1, mp_context=mp.get_context("spawn")) as pool:
        fut = pool.submit(_run_network_case, case, solver_timeout_s)
        return fut.result(timeout=case.timeout_s)


def run_benchmarks(repeats: int, case_timeout_s: float, solver_timeout_s: float) -> dict[str, Any]:
    total_area = 0.25 * math.pi * (0.8) ** 2

    grids = [(1, 1), (2, 2), (5, 5), (10, 10)]

    cases = []
    for s, p in grids:
        ht_area = max(50.0, float(s * p * 2.0))
        for comp in [False, True]:
            reg_str = "comp" if comp else "incomp"
            cases.append(
                BenchmarkCase(
                    name=f"grid_{s}x{p}_{reg_str}",
                    timeout_s=case_timeout_s,
                    m_dot_air=10.0,
                    T_ad=1700.0,
                    n_serial=s,
                    n_parallel=p,
                    total_area=total_area,
                    total_ht_area=ht_area,
                    use_compressible=comp,
                )
            )

    run_stamp = datetime.now(UTC).isoformat()
    results: dict[str, Any] = {
        "run_utc": run_stamp,
        "repeats": repeats,
        "case_timeout_s": case_timeout_s,
        "solver_timeout_s": solver_timeout_s,
        "cases": {},
    }

    print(f"--- Running Network Solver Benchmarks ({len(cases)} cases, max {repeats} repeats) ---")

    for case in cases:
        times: list[float] = []
        success_count = 0
        timeout_count = 0
        error_messages: list[str] = []
        final_norms: list[float] = []

        # Speed up heavy cases
        case_repeats = 1 if case.n_serial * case.n_parallel >= 25 else repeats

        print(f"Testing {case.name} ({case_repeats} repeats)...", flush=True)

        for _ in range(case_repeats):
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
            "attempts": case_repeats,
            "timeout_count": timeout_count,
            "success_rate": success_count / case_repeats,
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
        f"- Repeats per case: `{results['repeats']}` (hard grids: `1`)",
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
            f"| `{case_name}` | {success_txt} | {median_s:.4f} | {max_s:.4f} | {median_norm:.3e} |"
        )

    lines.append("")
    lines.append("## Notes")
    lines.append(
        "- Tested heavily non-linear parameterized combustor grid arrays (from 1x1 up to 10x10)."
    )
    lines.append(
        '- "comp" forces strict Fanno/Orifice compressibility, which frequently crashes generic non-damped Newton solvers.'
    )
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
        default=300.0,
        help="Hard timeout (seconds) for each benchmark case process.",
    )
    parser.add_argument(
        "--solver-timeout",
        type=float,
        default=120.0,
        help="Timeout passed to NetworkSolver.solve() in seconds.",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        default=Path(__file__).parent.parent
        / "benchmarks"
        / "network_solver_benchmark_latest.json",
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
