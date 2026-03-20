#!/usr/bin/env python3
"""
Jacobian Baseline Comparison Tool

This script compares current Jacobian test results against the baseline
to detect accuracy regressions in the solver interface functions.

Usage:
    python scripts/compare_jacobian_baseline.py
    python scripts/compare_jacobian_baseline.py --current new_results.csv
    python scripts/compare_jacobian_baseline.py --baseline baseline_20260320.csv
"""

import argparse
import pandas as pd
import sys
import os
from pathlib import Path

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
    merged = pd.merge(baseline_df, current_df,
                     on=['test_name', 'derivative'],
                     suffixes=('_baseline', '_current'))

    if len(merged) == 0:
        print("No matching derivatives found between baseline and current")
        return None

    # Calculate changes
    merged['abs_diff_change'] = merged['abs_diff_current'] - merged['abs_diff_baseline']
    merged['rel_diff_change'] = merged['rel_diff_current'] - merged['rel_diff_baseline']

    # Find regressions (increased differences)
    regressions = merged[merged['abs_diff_change'] > 1e-10].copy()
    improvements = merged[merged['abs_diff_change'] < -1e-10].copy()

    return {
        'merged': merged,
        'regressions': regressions,
        'improvements': improvements,
        'total_compared': len(merged)
    }

def print_summary(results):
    """Print comparison summary."""
    if not results:
        return

    merged = results['merged']
    regressions = results['regressions']
    improvements = results['improvements']

    print("\n=== JACOBIAN BASELINE COMPARISON ===")
    print(f"Total derivatives compared: {results['total_compared']}")
    print(f"Regressions (worse accuracy): {len(regressions)}")
    print(f"Improvements (better accuracy): {len(improvements)}")
    print(f"Unchanged: {results['total_compared'] - len(regressions) - len(improvements)}")

    # Statistics
    max_abs_change = merged['abs_diff_change'].abs().max()
    max_rel_change = merged['rel_diff_change'].abs().max()

    print(f"\nMaximum absolute difference change: {max_abs_change:.2e}")
    print(f"Maximum relative difference change: {max_rel_change:.2e}")

    # Show top regressions
    if len(regressions) > 0:
        print(f"\n=== TOP 5 REGRESSIONS ===")
        worst_regressions = regressions.nlargest(5, 'abs_diff_change')
        for _, row in worst_regressions.iterrows():
            print(f"{row['test_name']}.{row['derivative']}: "
                  f"{row['abs_diff_baseline']:.2e} → {row['abs_diff_current']:.2e} "
                  f"(+{row['abs_diff_change']:.2e})")

    # Show top improvements
    if len(improvements) > 0:
        print(f"\n=== TOP 5 IMPROVEMENTS ===")
        best_improvements = improvements.nsmallest(5, 'abs_diff_change')
        for _, row in best_improvements.iterrows():
            print(f"{row['test_name']}.{row['derivative']}: "
                  f"{row['abs_diff_baseline']:.2e} → {row['abs_diff_current']:.2e} "
                  f"({row['abs_diff_change']:.2e})")

def main():
    parser = argparse.ArgumentParser(description='Compare Jacobian test results against baseline')
    parser.add_argument('--baseline', default='jacobian_baseline_current.csv',
                       help='Baseline CSV file (default: jacobian_baseline_current.csv)')
    parser.add_argument('--current',
                       help='Current results CSV file (if not provided, will run tests)')
    parser.add_argument('--threshold', type=float, default=1e-10,
                       help='Minimum change to consider as regression (default: 1e-10)')

    args = parser.parse_args()

    # Find baseline file
    baseline_path = Path(args.baseline)
    if not baseline_path.exists():
        # Try in project root
        baseline_path = Path(__file__).parent.parent / args.baseline
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
        print("Running current tests to get baseline comparison...")
        import subprocess
        import tempfile
        import datetime

        # Run tests and capture output
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        temp_file = f"jacobian_temp_{timestamp}.csv"

        try:
            result = subprocess.run([
                './test_solver_jacobians',
                '--gtest_print_time=0'
            ], env={**os.environ, 'JACOBIAN_OUTPUT_FILE': temp_file},
               capture_output=True, text=True, cwd='build')

            if result.returncode != 0:
                print(f"Error running tests: {result.stderr}")
                return 1

            current_df = load_csv(f"build/{temp_file}")
            if current_df is None:
                return 1

        finally:
            # Clean up temp file
            temp_path = Path(f"build/{temp_file}")
            if temp_path.exists():
                temp_path.unlink()

    # Compare results
    results = compare_differences(baseline_df, current_df)
    if results is None:
        return 1

    print_summary(results)

    # Return error code if there are regressions
    if len(results['regressions']) > 0:
        print(f"⚠️  DETECTED {len(results['regressions'])} REGRESSIONS")
        return 1
    else:
        print(f"\n✅ NO REGRESSIONS DETECTED")
        return 0

if __name__ == "__main__":
    sys.exit(main())
