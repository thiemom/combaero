#!/usr/bin/env python3
"""Friction factor and pressure drop calculations.

This example demonstrates:
- Friction factor correlations (Haaland, Serghides, Colebrook, Petukhov)
- Pressure drop calculations in pipes
- Effect of roughness on friction
- Comparison of different correlations
- Practical pipe flow design
"""

from __future__ import annotations

import numpy as np

import combaero as ca


def main() -> None:
    print("=" * 80)
    print("FRICTION FACTOR AND PRESSURE DROP CALCULATIONS")
    print("=" * 80)

    # =========================================================================
    # Example 1: Smooth pipe friction factors
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 1: Smooth Pipe Friction Factors")
    print("=" * 80)

    Re = 1e5  # Reynolds number
    e_D = 0.0  # Smooth pipe (zero roughness)

    print(f"\nReynolds number: {Re:.0e}")
    print(f"Relative roughness ε/D: {e_D}")

    # Compare different correlations
    f_haaland = ca.friction_haaland(Re, e_D)
    f_serghides = ca.friction_serghides(Re, e_D)
    f_colebrook = ca.friction_colebrook(Re, e_D)
    f_petukhov = ca.friction_petukhov(Re)

    print("\nFriction factors (smooth pipe):")
    print(f"  Haaland:    f = {f_haaland:.6f}")
    print(f"  Serghides:  f = {f_serghides:.6f}")
    print(f"  Colebrook:  f = {f_colebrook:.6f}")
    print(f"  Petukhov:   f = {f_petukhov:.6f}")

    # =========================================================================
    # Example 2: Effect of roughness
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 2: Effect of Pipe Roughness")
    print("=" * 80)

    Re = 1e5
    roughness_values = {
        "Smooth (drawn tubing)": 0.0,
        "Commercial steel": 0.00045 / 0.1,  # ε=0.045mm, D=100mm
        "Galvanized iron": 0.0015 / 0.1,  # ε=0.15mm, D=100mm
        "Cast iron": 0.0026 / 0.1,  # ε=0.26mm, D=100mm
    }

    print(f"\nReynolds number: {Re:.0e}")
    print("\nFriction factors for different pipe materials:")
    for material, e_D in roughness_values.items():
        f = ca.friction_haaland(Re, e_D)
        print(f"  {material:25s}: f = {f:.6f}  (ε/D = {e_D:.6f})")

    # =========================================================================
    # Example 3: Pressure drop calculation
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 3: Pressure Drop in Pipe")
    print("=" * 80)

    # Pipe geometry
    D = 0.1  # Diameter [m]
    L = 100.0  # Length [m]
    e = 0.045e-3  # Roughness (commercial steel) [m]
    e_D = e / D

    # Air properties at 300K
    T = 300.0  # Temperature [K]
    P = 101325.0  # Pressure [Pa]
    X_air = ca.standard_dry_air_composition()

    # Flow conditions
    v = 10.0  # Velocity [m/s]

    # Calculate Reynolds number
    rho = ca.density(T, P, X_air)
    mu = ca.viscosity(T, P, X_air)
    Re = rho * v * D / mu

    # Calculate friction factor
    f = ca.friction_haaland(Re, e_D)

    # Pressure drop: ΔP = f * (L/D) * (ρ*v²/2)
    dP = f * (L / D) * (rho * v**2 / 2.0)

    print("\nPipe specifications:")
    print(f"  Diameter:        D = {D*1000:.1f} mm")
    print(f"  Length:          L = {L:.1f} m")
    print(f"  Roughness:       ε = {e*1e6:.1f} μm (commercial steel)")
    print(f"  Relative rough:  ε/D = {e_D:.6f}")

    print("\nFlow conditions:")
    print(f"  Fluid:           Air at {T:.0f} K")
    print(f"  Velocity:        v = {v:.1f} m/s")
    print(f"  Density:         ρ = {rho:.3f} kg/m³")
    print(f"  Viscosity:       μ = {mu*1e6:.2f} μPa·s")
    print(f"  Reynolds number: Re = {Re:.2e}")

    print("\nResults:")
    print(f"  Friction factor: f = {f:.6f}")
    print(f"  Pressure drop:   ΔP = {dP:.1f} Pa = {dP/1000:.3f} kPa")
    print(f"  Pressure drop:   ΔP/L = {dP/L:.2f} Pa/m")

    # =========================================================================
    # Example 4: Reynolds number sweep
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 4: Friction Factor vs Reynolds Number")
    print("=" * 80)

    e_D = 0.001  # Moderate roughness
    Re_values = np.logspace(3.5, 6, 10)  # 3162 to 1e6

    print(f"\nRelative roughness: ε/D = {e_D:.4f}")
    print(f"\n{'Re':>12s}  {'f (Haaland)':>12s}  {'f (Serghides)':>12s}  {'Diff %':>10s}")
    print("-" * 50)

    for Re in Re_values:
        f_h = ca.friction_haaland(Re, e_D)
        f_s = ca.friction_serghides(Re, e_D)
        diff = abs(f_h - f_s) / f_s * 100
        print(f"{Re:12.2e}  {f_h:12.6f}  {f_s:12.6f}  {diff:10.4f}")

    # =========================================================================
    # Example 5: Practical design - sizing a pipe
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 5: Pipe Sizing for Maximum Pressure Drop")
    print("=" * 80)

    # Design requirements
    mdot = 0.5  # Mass flow rate [kg/s]
    L = 50.0  # Pipe length [m]
    dP_max = 5000.0  # Maximum allowable pressure drop [Pa]
    e = 0.045e-3  # Commercial steel roughness [m]

    # Air properties
    T = 300.0
    P = 101325.0
    X_air = ca.standard_dry_air_composition()
    rho = ca.density(T, P, X_air)
    mu = ca.viscosity(T, P, X_air)

    print("\nDesign requirements:")
    print(f"  Mass flow rate:  ṁ = {mdot:.2f} kg/s")
    print(f"  Pipe length:     L = {L:.1f} m")
    print(f"  Max pressure drop: ΔP_max = {dP_max/1000:.1f} kPa")
    print(f"  Roughness:       ε = {e*1e6:.1f} μm")

    # Try different diameters
    print(
        f"\n{'D [mm]':>10s}  {'v [m/s]':>10s}  {'Re':>12s}  {'f':>10s}  {'ΔP [kPa]':>12s}  {'OK?':>6s}"
    )
    print("-" * 70)

    for D in [0.05, 0.075, 0.1, 0.125, 0.15]:
        A = np.pi * (D / 2) ** 2
        v = mdot / (rho * A)
        Re = rho * v * D / mu
        e_D = e / D
        f = ca.friction_haaland(Re, e_D)
        dP = f * (L / D) * (rho * v**2 / 2.0)
        ok = "✓" if dP <= dP_max else "✗"
        print(f"{D*1000:10.1f}  {v:10.2f}  {Re:12.2e}  {f:10.6f}  {dP/1000:12.3f}  {ok:>6s}")

    print("\n✓ = Meets pressure drop requirement")
    print("✗ = Exceeds pressure drop limit")

    # =========================================================================
    # Example 6: Laminar vs Turbulent
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 6: Laminar vs Turbulent Flow")
    print("=" * 80)

    print(f"\n{'Re':>10s}  {'Regime':>12s}  {'f (formula)':>15s}  {'f (correlation)':>18s}")
    print("-" * 60)

    for Re in [1000, 2000, 2300, 3000, 5000, 10000]:
        if Re < 2300:
            regime = "Laminar"
            f_formula = 64.0 / Re
            f_corr = "N/A"
        else:
            regime = "Turbulent"
            f_formula = "N/A"
            f_corr = f"{ca.friction_haaland(Re, 0.0):.6f}"

        f_formula_str = f"{f_formula:.6f}" if isinstance(f_formula, float) else f_formula
        print(f"{Re:10.0f}  {regime:>12s}  {f_formula_str:>15s}  {f_corr:>18s}")

    print("\nNote: Correlations are valid for turbulent flow (Re > ~2300)")
    print("For laminar flow, use: f = 64/Re")


if __name__ == "__main__":
    main()
