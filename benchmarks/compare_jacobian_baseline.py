#!/usr/bin/env python3
"""
Jacobian Baseline Comparison Tool

This script compares current Jacobian test results against the baseline
to detect accuracy regressions in the solver interface functions.

It integrates functionality for detailed printing and record saving.

Usage:
    python benchmarks/compare_jacobian_baseline.py
    python benchmarks/compare_jacobian_baseline.py --verbose  # Show detailed diffs
    python benchmarks/compare_jacobian_baseline.py --save     # Persist results to benchmarks/runs/
"""

import argparse
import datetime
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


def load_csv(filepath):
    """Load CSV file with error handling."""
    try:
        df = pd.read_csv(filepath)
        return df
    except FileNotFoundError:
        print(f"Error: File {filepath} not found")
        return None
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def compare_differences(baseline_df, current_df):
    """Compare Jacobian differences between baseline and current results."""
    # Merge on test_name and derivative
    merged = pd.merge(
        baseline_df, current_df, on=["test_name", "derivative"], suffixes=("_baseline", "_current")
    )

    if len(merged) == 0:
        print("No matching derivatives found between baseline and current")
        return None

    # Calculate changes
    merged["abs_diff_change"] = merged["abs_diff_current"] - merged["abs_diff_baseline"]
    merged["rel_diff_change"] = merged["rel_diff_current"] - merged["rel_diff_baseline"]

    # Find regressions (increased differences)
    regressions = merged[merged["abs_diff_change"] > 1e-10].copy()
    improvements = merged[merged["abs_diff_change"] < -1e-10].copy()

    return {
        "merged": merged,
        "regressions": regressions,
        "improvements": improvements,
        "total_compared": len(merged),
    }


def print_summary(results):
    """Print comparison summary."""
    if not results:
        return

    merged = results["merged"]
    regressions = results["regressions"]
    improvements = results["improvements"]

    print("\n=== JACOBIAN BASELINE COMPARISON ===")
    print(f"Total derivatives compared: {results['total_compared']}")
    print(f"Regressions (worse accuracy): {len(regressions)}")
    print(f"Improvements (better accuracy): {len(improvements)}")
    print(f"Unchanged: {results['total_compared'] - len(regressions) - len(improvements)}")

    # Statistics
    max_abs_change = merged["abs_diff_change"].abs().max()
    max_rel_change = merged["rel_diff_change"].abs().max()

    print(f"\nMaximum absolute difference change: {max_abs_change:.2e}")
    print(f"Maximum relative difference change: {max_rel_change:.2e}")

    # Show top regressions
    if len(regressions) > 0:
        print("\n=== TOP 5 REGRESSIONS ===")
        worst_regressions = regressions.nlargest(5, "abs_diff_change")
        for _, row in worst_regressions.iterrows():
            print(
                f"{row['test_name']}.{row['derivative']}: "
                f"{row['abs_diff_baseline']:.2e} → {row['abs_diff_current']:.2e} "
                f"(+{row['abs_diff_change']:.2e})"
            )

    # Show top improvements
    if len(improvements) > 0:
        print("\n=== TOP 5 IMPROVEMENTS ===")
        best_improvements = improvements.nsmallest(5, "abs_diff_change")
        for _, row in best_improvements.iterrows():
            print(
                f"{row['test_name']}.{row['derivative']}: "
                f"{row['abs_diff_baseline']:.2e} → {row['abs_diff_current']:.2e} "
                f"({row['abs_diff_change']:.2e})"
            )


def main():
    parser = argparse.ArgumentParser(description="Compare Jacobian test results against baseline")
    root_dir = Path(__file__).resolve().parent.parent
    benchmarks_dir = root_dir / "benchmarks"
    runs_dir = benchmarks_dir / "runs"

    parser.add_argument(
        "--baseline",
        default=str(benchmarks_dir / "jacobian_baseline_current.csv"),
        help="Baseline CSV file",
    )
    parser.add_argument(
        "--current", help="Current results CSV file (if not provided, will run tests)"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=1e-10,
        help="Minimum change to consider as regression (default: 1e-10)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable detailed difference reporting (PRINT_JACOBIAN_DIFFERENCES=1)",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save current results to benchmarks/runs/ instead of a temporary file",
    )

    args = parser.parse_args()

    # Find baseline file
    baseline_path = Path(args.baseline)
    if not baseline_path.exists():
        print(f"Baseline file not found: {args.baseline}")
        return 1

    # Load baseline
    print(f"Loading baseline: {baseline_path}")
    baseline_df = load_csv(baseline_path)
    if baseline_df is None:
        return 1

    # Get current results
    if args.current:
        current_df = load_csv(args.current)
        if current_df is None:
            return 1
    else:
        print("Running current tests...")

        timestamp = datetime.datetime.now(datetime.UTC).strftime("%Y%m%d_%H%M%SZ")
        if args.save:
            runs_dir.mkdir(parents=True, exist_ok=True)
            output_file_name = f"jacobian_test_{timestamp}.csv"
            output_path = runs_dir / output_file_name
        else:
            output_file_name = f"jacobian_temp_{timestamp}.csv"
            output_path = root_dir / "build" / output_file_name

        env = {**os.environ}
        if args.verbose:
            env["PRINT_JACOBIAN_DIFFERENCES"] = "1"
        env["JACOBIAN_OUTPUT_FILE"] = str(output_path)

        try:
            # Ensure we are in build directory or find the executable
            build_dir = root_dir / "build"
            exe_path = build_dir / "test_solver_jacobians"

            if not exe_path.exists():
                print(f"Error: Executable not found at {exe_path}")
                return 1

            result = subprocess.run(
                [str(exe_path), "--gtest_print_time=0"],
                env=env,
                capture_output=not args.verbose,
                text=True,
                cwd=str(build_dir),
            )

            if result.returncode != 0 and not args.verbose:
                print(f"Error running tests: {result.stderr}")
                return 1

            # If verbose, the output already went to stdout.
            # If not verbose, we might want to see if there were any errors.

            # The output file is written relative to the CWD of the process if it's just a filename.
            # Our C++ code uses: char* val = getenv("JACOBIAN_OUTPUT_FILE"); if(val) ...
            # I'll check both root and build just in case.

            final_output_path = output_path
            if not final_output_path.exists():
                # Check relative to build
                alt_path = build_dir / output_file_name
                if alt_path.exists():
                    final_output_path = alt_path
                else:
                    print(f"Error: Output file not produced at {output_path}")
                    return 1

            current_df = load_csv(final_output_path)
            if current_df is None:
                return 1

            if args.save:
                print(f"Saved current results to: {final_output_path}")

        finally:
            # Clean up temp file if not saving
            if not args.save and output_path.exists():
                output_path.unlink()
            elif not args.save:
                alt_path = build_dir / output_file_name
                if alt_path.exists():
                    alt_path.unlink()

    # Compare results
    results = compare_differences(baseline_df, current_df)
    if results is None:
        return 1

    print_summary(results)

    # Return error code if there are regressions
    if len(results["regressions"]) > 0:
        print(f"⚠️  DETECTED {len(results['regressions'])} REGRESSIONS")
        return 1
    else:
        print("\n✅ NO REGRESSIONS DETECTED")
        return 0


if __name__ == "__main__":
    sys.exit(main())
