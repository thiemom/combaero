#!/usr/bin/env python3
"""Compare standalone vs network implementations of combustor liner cooling.

Runs both implementations with identical inputs and verifies results match.
"""

from __future__ import annotations

import sys
from pathlib import Path

import combaero as cb

# Import solve functions from both implementations
sys.path.insert(0, str(Path(__file__).resolve().parent))
from combustor_liner_cooling_network import solve_operating_point_network as solve_network
from combustor_liner_cooling_standalone import solve_operating_point as solve_standalone


def compare_results(standalone: dict, network: dict, tol_rel: float = 0.01) -> bool:
    """Compare results from both implementations.

    Args:
        standalone: Results from standalone implementation
        network: Results from network implementation
        tol_rel: Relative tolerance for comparison (default 1%)

    Returns:
        True if all key results match within tolerance
    """
    keys_to_compare = [
        "u_cool",
        "h_cool",
        "h_hot",
        "Re",
        "Nu",
        "dP_cool",
        "T_cool_exit",
        "T_hot",
        "FAR",
        "q_wall",
        "T_hot_surf",
        "T_metal_hot",
        "eta_thp",
    ]

    all_match = True
    max_error = 0.0

    print("\n" + "=" * 80)
    print("COMPARISON: Standalone vs Network Implementation")
    print("=" * 80)
    print(f"{'Parameter':20} {'Standalone':>15} {'Network':>15} {'Rel Error %':>15} {'Status':>10}")
    print("-" * 80)

    for key in keys_to_compare:
        val_s = standalone[key]
        val_n = network[key]

        # Compute relative error
        rel_error = abs(val_n - val_s) / abs(val_s) if abs(val_s) > 1e-12 else abs(val_n - val_s)

        max_error = max(max_error, rel_error)
        status = "✓ PASS" if rel_error < tol_rel else "✗ FAIL"

        if rel_error >= tol_rel:
            all_match = False

        print(f"{key:20} {val_s:15.6g} {val_n:15.6g} {rel_error * 100:15.4f} {status:>10}")

    print("-" * 80)
    print(f"Maximum relative error: {max_error * 100:.4f}%")
    print(f"Overall: {'✓ ALL PASS' if all_match else '✗ SOME FAILURES'}")
    print("=" * 80)

    return all_match


def main() -> int:

    print("=" * 80)
    print("Combustor Liner Cooling — Implementation Comparison")
    print("=" * 80)

    # Test conditions (same as both implementations)
    PR = 20.0
    P2 = PR * 101325.0
    T2 = 667.0  # K
    RH_in = 0.60
    COT_target = 1758.0  # K
    dP_burner_frac = 0.02
    D_ann = 0.30  # m
    u_hot_ann = 25.0  # m/s

    X_cool = cb.species.humid_air(288.15, 101325.0, RH_in)
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.empty()
    X_ch4[cb.species.indices["CH4"]] = 1.0

    mdot_total = 1.0  # kg/s

    # Channel geometry
    D_ch = 0.008  # m
    L_ch = 0.35  # m
    N_ch = 8

    # Materials
    k_inconel = cb.k_inconel718(900.0)
    t_wall = 0.003
    wall_layers = [(t_wall, k_inconel)]

    # Test cases: different f_cool and channel types
    test_cases = [
        (0.10, "Smooth"),
        (0.20, "Smooth"),
        (0.15, "Ribbed"),
        (0.15, "Dimpled"),
    ]

    all_tests_pass = True

    for f_cool, ch_name in test_cases:
        print(f"\n{'=' * 80}")
        print(f"Test Case: f_cool={f_cool:.2f}, channel={ch_name}")
        print(f"{'=' * 80}")

        # Run standalone
        print("\nRunning standalone implementation...")
        result_standalone = solve_standalone(
            f_cool,
            mdot_total,
            T2,
            P2,
            X_cool,
            X_air,
            X_ch4,
            COT_target,
            dP_burner_frac,
            D_ch,
            L_ch,
            N_ch,
            D_ann,
            u_hot_ann,
            ch_name,
            wall_layers,
        )

        # Run network
        print("Running network implementation...")
        result_network = solve_network(
            f_cool,
            mdot_total,
            T2,
            P2,
            X_cool,
            X_air,
            X_ch4,
            COT_target,
            dP_burner_frac,
            D_ch,
            L_ch,
            N_ch,
            D_ann,
            u_hot_ann,
            ch_name,
            wall_layers,
        )

        # Compare
        test_pass = compare_results(result_standalone, result_network, tol_rel=0.01)

        if not test_pass:
            all_tests_pass = False
            print(f"\n⚠️  Test case FAILED: f_cool={f_cool:.2f}, channel={ch_name}")
        else:
            print(f"\n✓ Test case PASSED: f_cool={f_cool:.2f}, channel={ch_name}")

    # Final summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    if all_tests_pass:
        print("✓ ALL TEST CASES PASSED")
        print("Both implementations produce matching results within 1% tolerance.")
        return 0
    else:
        print("✗ SOME TEST CASES FAILED")
        print("Discrepancies found between implementations. Investigation needed.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
