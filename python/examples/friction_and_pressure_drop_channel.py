#!/usr/bin/env python3
"""Friction factor and pressure drop calculations.

This example demonstrates:
- Friction factor correlations (Haaland, Serghides, Colebrook, Petukhov)
- Pressure drop calculations in channels
- Effect of roughness on friction
- Comparison of different correlations
- Practical channel flow design
"""

from __future__ import annotations

import numpy as np

import combaero as cb


def main() -> None:
    print("=" * 80)
    print("FRICTION FACTOR AND PRESSURE DROP CALCULATIONS")
    print("=" * 80)

    # =========================================================================
    # Example 1: Smooth channel friction factors
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 1: Smooth Channel Friction Factors")
    print("=" * 80)

    Re = 1e5  # Reynolds number
    e_D = 0.0  # Smooth channel (zero roughness)

    print(f"\nReynolds number: {Re:.0e}")
    print(f"Relative roughness eps/D: {e_D}")

    # Compare different correlations
    f_haaland = cb.friction_haaland(Re, e_D)
    f_serghides = cb.friction_serghides(Re, e_D)
    f_colebrook = cb.friction_colebrook(Re, e_D)
    f_petukhov = cb.friction_petukhov(Re)

    print("\nFriction factors (smooth channel):")
    print(f"  Haaland:    f = {f_haaland:.6f}")
    print(f"  Serghides:  f = {f_serghides:.6f}")
    print(f"  Colebrook:  f = {f_colebrook:.6f}")
    print(f"  Petukhov:   f = {f_petukhov:.6f}")

    # =========================================================================
    # Example 2: Effect of roughness
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 2: Effect of Channel Roughness")
    print("=" * 80)

    Re = 1e5
    D = 0.1  # Channel diameter [m]

    # Use channel_roughness database for standard materials
    materials = ["drawn_tubing", "commercial_steel", "galvanized_iron", "cast_iron"]

    print(f"\nReynolds number: {Re:.0e}")
    print(f"Channel diameter: {D * 1000:.0f} mm")
    print("\nFriction factors for different channel materials:")
    for material in materials:
        eps = cb.channel_roughness(material)  # Get roughness from database
        e_D = eps / D
        f = cb.friction_haaland(Re, e_D)
        print(f"  {material:25s}: f = {f:.6f}  (eps = {eps * 1e6:.1f} mum, eps/D = {e_D:.6f})")

    # =========================================================================
    # Example 3: Pressure drop calculation
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 3: Pressure Drop in Channel")
    print("=" * 80)

    # Channel geometry
    D = 0.1  # Diameter [m]
    L = 100.0  # Length [m]
    eps = cb.channel_roughness("commercial_steel")  # Get roughness from database
    e_D = eps / D

    # Air properties at 300K
    T = 300.0  # Temperature [K]
    P = 101325.0  # Pressure [Pa]
    X_air = cb.species.dry_air()

    # Flow conditions
    v = 10.0  # Velocity [m/s]

    # Use pressure_drop_channel composite function (calculates Re, f, and DeltaP)
    dP, Re, f = cb.pressure_drop_channel(T, P, X_air, v, D, L, eps, "haaland")

    # Also get density and viscosity for display
    rho = cb.density(T, P, X_air)
    mu = cb.viscosity(T, P, X_air)

    print("\nChannel specifications:")
    print(f"  Diameter:        D = {D * 1000:.1f} mm")
    print(f"  Length:          L = {L:.1f} m")
    print(f"  Roughness:       eps = {eps * 1e6:.1f} mum (commercial steel)")
    print(f"  Relative rough:  eps/D = {e_D:.6f}")

    print("\nFlow conditions:")
    print(f"  Fluid:           Air at {T:.0f} K")
    print(f"  Velocity:        v = {v:.1f} m/s")
    print(f"  Density:         rho = {rho:.3f} kg/m^3")
    print(f"  Viscosity:       mu = {mu * 1e6:.2f} muPa*s")
    print(f"  Reynolds number: Re = {Re:.2e}")

    print("\nResults:")
    print(f"  Friction factor: f = {f:.6f}")
    print(f"  Pressure drop:   DeltaP = {dP:.1f} Pa = {dP / 1000:.3f} kPa")
    print(f"  Pressure drop:   DeltaP/L = {dP / L:.2f} Pa/m")

    # =========================================================================
    # Example 4: Reynolds number sweep
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 4: Friction Factor vs Reynolds Number")
    print("=" * 80)

    e_D = 0.001  # Moderate roughness
    Re_values = np.logspace(3.5, 6, 10)  # 3162 to 1e6

    print(f"\nRelative roughness: eps/D = {e_D:.4f}")
    print(f"\n{'Re':>12s}  {'f (Haaland)':>12s}  {'f (Serghides)':>12s}  {'Diff %':>10s}")
    print("-" * 50)

    for Re in Re_values:
        f_h = cb.friction_haaland(Re, e_D)
        f_s = cb.friction_serghides(Re, e_D)
        diff = abs(f_h - f_s) / f_s * 100
        print(f"{Re:12.2e}  {f_h:12.6f}  {f_s:12.6f}  {diff:10.4f}")

    # =========================================================================
    # Example 5: Practical design - sizing a channel
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 5: Channel Sizing for Maximum Pressure Drop")
    print("=" * 80)

    # Design requirements
    mdot = 0.5  # Mass flow rate [kg/s]
    L = 50.0  # Channel length [m]
    dP_max = 5000.0  # Maximum allowable pressure drop [Pa]
    eps = cb.channel_roughness("commercial_steel")  # Get roughness from database

    # Air properties
    T = 300.0
    P = 101325.0
    X_air = cb.species.dry_air()
    rho = cb.density(T, P, X_air)
    mu = cb.viscosity(T, P, X_air)

    print("\nDesign requirements:")
    print(f"  Mass flow rate:  mdot = {mdot:.2f} kg/s")
    print(f"  Channel length:     L = {L:.1f} m")
    print(f"  Max pressure drop: DeltaP_max = {dP_max / 1000:.1f} kPa")
    print(f"  Roughness:       eps = {eps * 1e6:.1f} mum (commercial steel)")

    # Try different diameters
    print(
        f"\n{'D [mm]':>10s}  {'v [m/s]':>10s}  {'Re':>12s}  {'f':>10s}  {'DeltaP [kPa]':>12s}  {'OK?':>6s}"
    )
    print("-" * 70)

    for D in [0.05, 0.075, 0.1, 0.125, 0.15]:
        # Use channel_area helper
        A = cb.channel_area(D)
        v = mdot / (rho * A)
        # Use pressure_drop_channel composite function
        dP, Re, f = cb.pressure_drop_channel(T, P, X_air, v, D, L, eps, "haaland")
        ok = "OK" if dP <= dP_max else "FAIL"
        print(f"{D * 1000:10.1f}  {v:10.2f}  {Re:12.2e}  {f:10.6f}  {dP / 1000:12.3f}  {ok:>6s}")

    print("\nOK = Meets pressure drop requirement")
    print("FAIL = Exceeds pressure drop limit")

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
            f_corr = f"{cb.friction_haaland(Re, 0.0):.6f}"

        f_formula_str = f"{f_formula:.6f}" if isinstance(f_formula, float) else f_formula
        print(f"{Re:10.0f}  {regime:>12s}  {f_formula_str:>15s}  {f_corr:>18s}")

    print("\nNote: Correlations are valid for turbulent flow (Re > ~2300)")
    print("For laminar flow, use: f = 64/Re")


if __name__ == "__main__":
    main()
