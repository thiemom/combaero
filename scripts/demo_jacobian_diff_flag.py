#!/usr/bin/env python3
"""
Demonstration of the PRINT_JACOBIAN_DIFFERENCES flag for Jacobian tests.

This script shows how to use the environment variable to collect and analyze
Jacobian test differences across all C++ tests.

Usage:
    # Run tests normally (no detailed output, just summary)
    cd build && ./test_solver_jacobians

    # Run with detailed difference reporting
    cd build && PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians

    # Run specific tests with detailed reporting
    cd build && PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians --gtest_filter="*Orifice*"
"""

import os
import subprocess
import sys

def run_demo():
    print("=== Jacobian Difference Reporting Demo ===\n")

    # Change to build directory
    build_dir = "build"
    if not os.path.exists(build_dir):
        print("Error: build directory not found. Run 'cmake -B build -S .' first.")
        return 1

    os.chdir(build_dir)

    # Demo 1: Normal run (no detailed output)
    print("1. Normal test run (no detailed output):")
    print("   Command: ./test_solver_jacobians --gtest_filter='*OrificeDerivatives*'")
    result = subprocess.run(
        ["./test_solver_jacobians", "--gtest_filter=*OrificeDerivatives*", "--gtest_print_time=0"],
        capture_output=True, text=True
    )
    print("   Output (last few lines):")
    lines = result.stdout.strip().split('\n')
    for line in lines[-10:]:
        print(f"   {line}")
    print()

    # Demo 2: With detailed reporting
    print("2. With detailed difference reporting:")
    print("   Command: PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians --gtest_filter='*OrificeDerivatives*'")
    env = os.environ.copy()
    env["PRINT_JACOBIAN_DIFFERENCES"] = "1"
    result = subprocess.run(
        ["./test_solver_jacobians", "--gtest_filter=*OrificeDerivatives*", "--gtest_print_time=0"],
        capture_output=True, text=True, env=env
    )
    print("   Output (showing JACOBIAN_DIFF lines):")
    for line in result.stdout.split('\n'):
        if "JACOBIAN_DIFF" in line or "===" in line or "Test" in line:
            print(f"   {line}")
    print()

    # Demo 3: All tests with summary
    print("3. All tests with summary statistics:")
    print("   Command: PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians")
    print("   (This would show patterns across ALL Jacobian tests)")
    print("   Use this to identify which functions have the largest differences")
    print()

    print("=== How to Add Reporting to New Tests ===")
    print("To add difference reporting to a new test:")
    print("1. Call report_jacobian_difference() for each derivative:")
    print("   report_jacobian_difference('TestName', 'derivative_name', analytical, fd)")
    print("2. The function automatically collects stats and prints summary")
    print("3. Set PRINT_JACOBIAN_DIFFERENCES=1 to see detailed output")
    print()

    print("=== Pattern Analysis ===")
    print("The summary helps identify patterns like:")
    print("- Which functions have the largest absolute differences")
    print("- Which have the largest relative differences")
    print("- Whether certain types of derivatives are consistently less accurate")
    print()

    return 0

if __name__ == "__main__":
    sys.exit(run_demo())
