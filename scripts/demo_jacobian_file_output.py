#!/usr/bin/env python3
"""
Demonstration of Jacobian difference file output capability.

This script shows how to create persistent records of Jacobian test accuracy
using the new file output feature.
"""

import os
import subprocess
import sys
import pandas as pd
from datetime import datetime

def run_file_output_demo():
    print("=== Jacobian File Output Demo ===\n")

    # Change to build directory
    build_dir = "build"
    if not os.path.exists(build_dir):
        print("Error: build directory not found. Run 'cmake -B build -S .' first.")
        return 1

    os.chdir(build_dir)

    # Demo 1: Create timestamped record
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"jacobian_record_{timestamp}.csv"

    print(f"1. Creating persistent record: {output_file}")
    print("   Command: JACOBIAN_OUTPUT_FILE=jacobian_record.csv ./test_solver_jacobians --gtest_filter='*OrificeDerivatives*'")

    env = os.environ.copy()
    env["JACOBIAN_OUTPUT_FILE"] = "jacobian_record.csv"
    result = subprocess.run(
        ["./test_solver_jacobians", "--gtest_filter=*OrificeDerivatives*", "--gtest_print_time=0"],
        capture_output=True, text=True, env=env
    )

    if os.path.exists("jacobian_record.csv"):
        print("   ✓ File created successfully")

        # Show file contents
        with open("jacobian_record.csv", "r") as f:
            lines = f.readlines()
            print(f"   File contains {len(lines)} lines (including header)")
            print("   First few lines:")
            for i, line in enumerate(lines[:3]):
                print(f"     {line.strip()}")

        # Rename with timestamp
        os.rename("jacobian_record.csv", output_file)
        print(f"   ✓ Renamed to: {output_file}")
    else:
        print("   ✗ File was not created")

    print()

    # Demo 2: Append to daily log
    daily_file = f"jacobian_daily_{datetime.now().strftime('%Y%m%d')}.csv"
    print(f"2. Appending to daily log: {daily_file}")
    print("   Command: JACOBIAN_OUTPUT_FILE=daily_file.csv ./test_solver_jacobians --gtest_filter='*AdiabaticWallDerivatives*'")

    env["JACOBIAN_OUTPUT_FILE"] = daily_file
    result = subprocess.run(
        ["./test_solver_jacobians", "--gtest_filter=*AdiabaticWallDerivatives*", "--gtest_print_time=0"],
        capture_output=True, text=True, env=env
    )

    if os.path.exists(daily_file):
        with open(daily_file, "r") as f:
            lines = f.readlines()
            print(f"   ✓ Daily log now contains {len(lines)} total entries")
    print()

    # Demo 3: Show analysis potential
    print("3. Analysis capabilities with CSV data:")
    print("   The CSV files can be analyzed with:")
    print("   - pandas: pd.read_csv('jacobian_record.csv')")
    print("   - Excel: Direct import")
    print("   - R: read.csv()")
    print("   - Shell tools: grep, awk, sort")
    print()

    # Demo 4: Show sample analysis if pandas is available
    try:
        if os.path.exists(output_file):
            df = pd.read_csv(output_file)
            print("4. Sample pandas analysis:")
            print(f"   Records: {len(df)}")
            print(f"   Tests: {df['test_name'].nunique()}")
            print(f"   Max absolute difference: {df['abs_diff'].max():.2e}")
            print(f"   Max relative difference: {df['rel_diff'].max():.2e}")
            print()

            # Group by test
            print("   Differences by test:")
            grouped = df.groupby('test_name').agg({
                'abs_diff': ['max', 'mean'],
                'rel_diff': ['max', 'mean']
            }).round(10)
            print(grouped)

    except ImportError:
        print("4. Install pandas for analysis: pip install pandas")
    except Exception as e:
        print(f"4. Analysis failed: {e}")

    print()
    print("=== File Output Benefits ===")
    print("✓ Persistent records for historical tracking")
    print("✓ Easy data analysis with standard tools")
    print("✓ Time series analysis over git history")
    print("✓ Integration with CI/CD pipelines")
    print("✓ Statistical analysis and visualization")
    print()

    return 0

if __name__ == "__main__":
    sys.exit(run_file_output_demo())
