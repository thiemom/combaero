"""
Diagnostic utilities for NetworkSolver analysis and visualization.
"""

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


def export_diagnostic_report(report: dict[str, Any], filepath: str | Path) -> None:
    """
    Export diagnostic report to JSON file.

    Args:
        report: Diagnostic report from NetworkSolver.get_diagnostic_report()
        filepath: Output file path for JSON report
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, "w") as f:
        json.dump(report, f, indent=2)

    print(f"Diagnostic report exported to: {filepath}")


def create_grid_heatmaps(
    solver, report: dict[str, Any], output_dir: str | Path = "diagnostic_plots"
) -> None:
    """
    Create heatmap visualizations for grid network analysis.

    Args:
        solver: NetworkSolver instance with diagnostic data
        report: Diagnostic report from get_diagnostic_report()
        output_dir: Directory to save plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract grid data from unknown names
    unknown_names = report["unknown_names"]
    x0 = solver._diagnostic_data["x0"]
    solution = solver._diagnostic_data["solution"]

    # Parse grid structure from variable names
    grid_data = parse_grid_variables(unknown_names, x0, solution)

    if not grid_data:
        print("No grid structure detected in variable names")
        return

    # Create heatmaps
    create_heatmap(
        grid_data["x0"], "Initial Guess (x0)", output_dir / "initial_guess.png", grid_data["shape"]
    )
    create_heatmap(
        grid_data["solution"],
        "Final Solution",
        output_dir / "final_solution.png",
        grid_data["shape"],
    )
    create_heatmap(
        grid_data["difference"],
        "Solution - Initial Guess",
        output_dir / "difference.png",
        grid_data["shape"],
    )
    create_heatmap(
        grid_data["relative_difference"],
        "Relative Difference",
        output_dir / "relative_difference.png",
        grid_data["shape"],
    )

    print(f"Grid heatmaps saved to: {output_dir}")


def parse_grid_variables(
    unknown_names: list[str], x0: np.ndarray, solution: np.ndarray
) -> dict[str, Any]:
    """
    Parse grid variables from unknown names to extract 2D structure.

    Args:
        unknown_names: List of variable names
        x0: Initial guess vector
        solution: Final solution vector

    Returns:
        Dictionary with parsed grid data or empty dict if no grid found
    """
    # Look for momentum chamber pressure variables (grid structure)
    chamber_vars = [
        name for name in unknown_names if name.startswith("momentum_chamber_") and ".P" in name
    ]

    if not chamber_vars:
        return {}

    # Extract grid indices
    grid_indices = []
    grid_values_x0 = []
    grid_values_solution = []

    for var in chamber_vars:
        # Parse indices like "momentum_chamber_i_j.P"
        parts = var.split("_")
        if len(parts) >= 4:
            try:
                i = int(parts[2])
                j = int(parts[3].split(".")[0])

                idx = unknown_names.index(var)
                grid_indices.append((i, j))
                grid_values_x0.append(x0[idx])
                grid_values_solution.append(solution[idx])
            except (ValueError, IndexError):
                continue

    if not grid_indices:
        return {}

    # Determine grid shape
    i_indices, j_indices = zip(*grid_indices, strict=True)
    n_rows = max(i_indices) + 1
    n_cols = max(j_indices) + 1

    # Create 2D arrays
    x0_grid = np.full((n_rows, n_cols), np.nan)
    solution_grid = np.full((n_rows, n_cols), np.nan)

    for (i, j), val_x0, val_sol in zip(
        grid_indices, grid_values_x0, grid_values_solution, strict=True
    ):
        x0_grid[i, j] = val_x0
        solution_grid[i, j] = val_sol

    difference = solution_grid - x0_grid
    relative_difference = np.abs(difference / (x0_grid + 1e-10))

    return {
        "x0": x0_grid,
        "solution": solution_grid,
        "difference": difference,
        "relative_difference": relative_difference,
        "shape": (n_rows, n_cols),
        "variable_type": "pressure",
    }


def create_heatmap(data: np.ndarray, title: str, filepath: Path, shape: tuple[int, int]) -> None:
    """
    Create and save heatmap plot.

    Args:
        data: 2D array to plot
        title: Plot title
        filepath: Output file path
        shape: Grid shape (rows, cols)
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create heatmap
    im = ax.imshow(data, cmap="viridis", aspect="auto")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Value")

    # Set labels and title
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Column Index", fontsize=12)
    ax.set_ylabel("Row Index", fontsize=12)

    # Add grid
    ax.set_xticks(range(shape[1]))
    ax.set_yticks(range(shape[0]))
    ax.grid(True, alpha=0.3)

    # Add text annotations for small grids
    if shape[0] * shape[1] <= 100:
        for i in range(shape[0]):
            for j in range(shape[1]):
                if not np.isnan(data[i, j]):
                    ax.text(
                        j,
                        i,
                        f"{data[i, j]:.2e}",
                        ha="center",
                        va="center",
                        color="white",
                        fontsize=8,
                    )

    plt.tight_layout()
    plt.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close()


def print_diagnostic_summary(report: dict[str, Any]) -> None:
    """
    Print a concise summary of the diagnostic report.

    Args:
        report: Diagnostic report from get_diagnostic_report()
    """
    conv = report["convergence_analysis"]
    stats = report["overall_statistics"]

    print("\n" + "=" * 60)
    print("NETWORK SOLVER DIAGNOSTIC SUMMARY")
    print("=" * 60)

    print("\n[Convergence] Status:")
    print(f"   Converged: {conv['converged']}")
    print(f"   Method: {conv['solver_method']}")
    print(f"   Final Norm: {conv['final_norm']:.3e}")
    print(f"   Wall Time: {conv['wall_time']:.2f} s")
    print(f"   Unknowns: {conv['num_unknowns']}")

    print("\n[Solution] Changes:")
    print(f"   Max Absolute Change: {stats['max_abs_difference']:.3e}")
    print(f"   Mean Absolute Change: {stats['mean_abs_difference']:.3e}")
    print(f"   Max Relative Change: {stats['max_rel_difference']:.3%}")
    print(f"   Mean Relative Change: {stats['mean_rel_difference']:.3%}")

    print("\n[Most] Changed Variables:")
    for i, var in enumerate(report["most_changed_variables"][:3], 1):
        print(
            f"   {i}. {var['name']}: Delta={var['abs_difference']:.3e} ({var['rel_difference']:.1%})"
        )

    print("=" * 60)
