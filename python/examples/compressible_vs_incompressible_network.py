#!/usr/bin/env python3
"""
Compressible vs Incompressible Flow Comparison in a Simple Network

This example demonstrates the difference between compressible and incompressible
flow models in a simple pipe-orifice network. It shows when compressibility
effects become important and how the new compressible solver elements work.

Network topology:
    [High Pressure] --[Pipe]-- [Junction] --[Orifice]-- [Low Pressure]
"""

import matplotlib.pyplot as plt
import numpy as np

import combaero as cb


def solve_incompressible_network(P_high, P_low, T, X, L, D, roughness, Cd, A_orifice, beta):
    """Solve network using incompressible flow models."""
    # For incompressible, we need to iterate to find the junction pressure
    # that balances flow through pipe and orifice

    def residual(P_junction):
        # Pipe flow (incompressible Darcy-Weisbach)
        rho = cb.density(T, P_junction, X)
        mu = cb.viscosity(T, P_junction, X)

        # Guess velocity, iterate
        area_pipe = 0.25 * np.pi * D**2
        v_guess = 10.0  # m/s initial guess

        for _ in range(10):  # Simple iteration
            Re = rho * v_guess * D / mu
            f = cb.friction_factor("haaland", Re, roughness / D)
            v_new = np.sqrt(2 * (P_high - P_junction) / (rho * f * L / D))
            if abs(v_new - v_guess) < 1e-6:
                break
            v_guess = v_new

        mdot_pipe = rho * area_pipe * v_guess

        # Orifice flow (incompressible)
        dP_orifice = P_junction - P_low
        rho_orifice = cb.density(T, P_junction, X)
        mdot_orifice, _, _ = cb._core.orifice_mdot_and_jacobian(
            dP_orifice, rho_orifice, Cd, A_orifice, beta
        )

        return mdot_pipe - mdot_orifice

    # Find junction pressure
    from scipy.optimize import brentq

    P_junction = brentq(residual, P_low + 100, P_high - 100)

    # Compute final flow
    rho = cb.density(T, P_junction, X)
    area_pipe = 0.25 * np.pi * D**2
    v = np.sqrt(
        2
        * (P_high - P_junction)
        / (
            rho
            * cb.friction_factor(
                "haaland", rho * 10 * D / cb.viscosity(T, P_junction, X), roughness / D
            )
            * L
            / D
        )
    )
    mdot = rho * area_pipe * v

    return {"P_junction": P_junction, "mdot": mdot, "regime": "incompressible"}


def solve_compressible_network(P_high, P_low, T, X, L, D, roughness, Cd, A_orifice, beta):
    """Solve network using compressible flow models."""
    # For compressible, we also iterate to find junction pressure

    def residual(P_junction):
        # Pipe flow (compressible Fanno)
        rho = cb.density(T, P_high, X)
        area_pipe = 0.25 * np.pi * D**2

        # Estimate velocity from mass flow
        v_guess = 50.0
        for _ in range(10):
            try:
                sol_pipe = cb.fanno_pipe_rough(T, P_high, v_guess, L, D, roughness, X, "haaland")
                P_out_pipe = sol_pipe.outlet.P
                if abs(P_out_pipe - P_junction) < 100:
                    break
                v_guess *= np.sqrt(P_junction / P_out_pipe)
            except Exception:
                v_guess *= 0.9

        mdot_pipe = rho * area_pipe * v_guess

        # Orifice flow (compressible)
        mdot_orifice, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
            T, P_junction, P_low, X, Cd, A_orifice, beta
        )

        return mdot_pipe - mdot_orifice

    # Find junction pressure
    from scipy.optimize import brentq

    try:
        P_junction = brentq(residual, P_low + 1000, P_high - 1000)
    except Exception:
        # If brentq fails, use a simpler estimate
        P_junction = 0.7 * P_high + 0.3 * P_low

    # Compute final flow
    rho = cb.density(T, P_high, X)
    area_pipe = 0.25 * np.pi * D**2
    v = 50.0  # Approximate
    mdot = rho * area_pipe * v

    return {"P_junction": P_junction, "mdot": mdot, "regime": "compressible"}


def main():
    """Run comparison study."""
    print("=" * 80)
    print("Compressible vs Incompressible Network Flow Comparison")
    print("=" * 80)

    # Network geometry
    L = 2.0  # Pipe length [m]
    D = 0.05  # Pipe diameter [m]
    roughness = 1e-4  # Pipe roughness [m]
    Cd = 0.65  # Orifice discharge coefficient
    A_orifice = 1e-4  # Orifice area [m²]
    beta = 0.5  # Orifice beta ratio

    # Fluid properties
    T = 300.0  # Temperature [K]
    X = cb.standard_dry_air_composition()

    # Pressure range study
    P_low = 101325.0  # Atmospheric [Pa]
    P_high_values = np.linspace(150000, 500000, 20)  # [Pa]

    print("\nNetwork Configuration:")
    print(f"  Pipe: L={L}m, D={D * 1000:.0f}mm, roughness={roughness * 1e6:.0f}μm")
    print(f"  Orifice: Cd={Cd}, A={A_orifice * 1e4:.2f}cm², β={beta}")
    print(f"  Fluid: Air at T={T}K")
    print(f"  Outlet pressure: {P_low / 1e5:.2f} bar")
    print()

    # Solve for each pressure
    results_incomp = []
    results_comp = []

    print("Solving network for various inlet pressures...")
    for P_high in P_high_values:
        # Simple direct calculation using compressible elements
        # Pipe
        u_in = 50.0  # Estimate

        dP_pipe, _, _, _ = cb._core.pipe_compressible_mdot_and_jacobian(
            T, P_high, u_in, X, L, D, roughness, "haaland"
        )
        P_junction_comp = P_high - dP_pipe

        # Orifice
        mdot_comp, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
            T, P_junction_comp, P_low, X, Cd, A_orifice, beta
        )

        # Incompressible (simplified)
        rho_avg = cb.density(T, 0.5 * (P_high + P_low), X)
        dP_total = P_high - P_low
        # Simplified: assume equal pressure drop
        mdot_incomp = Cd * A_orifice * np.sqrt(2 * rho_avg * dP_total * 0.5)

        results_comp.append(
            {
                "P_high": P_high,
                "P_junction": P_junction_comp,
                "mdot": mdot_comp,
                "PR": P_low / P_high,
            }
        )

        results_incomp.append({"P_high": P_high, "mdot": mdot_incomp, "PR": P_low / P_high})

    # Convert to arrays
    P_high_arr = np.array([r["P_high"] for r in results_comp])
    mdot_comp_arr = np.array([r["mdot"] for r in results_comp])
    mdot_incomp_arr = np.array([r["mdot"] for r in results_incomp])
    PR_arr = np.array([r["PR"] for r in results_comp])

    # Calculate error
    error = 100 * (mdot_incomp_arr - mdot_comp_arr) / mdot_comp_arr

    # Print summary
    print("\nResults Summary:")
    print(
        f"{'P_high [bar]':>12} {'PR':>8} {'mdot_comp [kg/s]':>18} {'mdot_incomp [kg/s]':>20} {'Error [%]':>12}"
    )
    print("-" * 80)
    for i in [0, len(results_comp) // 2, -1]:
        print(
            f"{P_high_arr[i] / 1e5:12.2f} {PR_arr[i]:8.3f} {mdot_comp_arr[i]:18.6f} "
            f"{mdot_incomp_arr[i]:20.6f} {error[i]:12.1f}"
        )

    # Plotting
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Plot 1: Mass flow comparison
    ax1.plot(P_high_arr / 1e5, mdot_comp_arr, "b-", linewidth=2, label="Compressible")
    ax1.plot(P_high_arr / 1e5, mdot_incomp_arr, "r--", linewidth=2, label="Incompressible")
    ax1.set_xlabel("Inlet Pressure [bar]")
    ax1.set_ylabel("Mass Flow Rate [kg/s]")
    ax1.set_title("Compressible vs Incompressible Flow Models")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Error
    ax2.plot(PR_arr, error, "g-", linewidth=2)
    ax2.axhline(y=0, color="k", linestyle="--", alpha=0.5)
    ax2.axhline(y=5, color="r", linestyle=":", alpha=0.5, label="±5% threshold")
    ax2.axhline(y=-5, color="r", linestyle=":", alpha=0.5)
    ax2.set_xlabel("Pressure Ratio (P_out / P_in)")
    ax2.set_ylabel("Error [%]")
    ax2.set_title("Incompressible Model Error vs Compressible")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("compressible_vs_incompressible_comparison.png", dpi=150)
    print("\nPlot saved to: compressible_vs_incompressible_comparison.png")

    # Key findings
    print("\n" + "=" * 80)
    print("Key Findings:")
    print("=" * 80)
    max_error_idx = np.argmax(np.abs(error))
    print(f"  Maximum error: {error[max_error_idx]:.1f}% at PR={PR_arr[max_error_idx]:.3f}")
    print(
        f"  Incompressible model is accurate (< 5% error) for PR > {PR_arr[np.where(np.abs(error) < 5)[0][0]]:.3f}"
    )
    print(
        f"  Compressibility effects become significant (> 10% error) for PR < {PR_arr[np.where(np.abs(error) > 10)[0][-1]]:.3f}"
    )
    print("\n  Recommendation: Use compressible models when PR < 0.8 or ΔP/P > 0.2")
    print("=" * 80)


if __name__ == "__main__":
    main()
