#!/usr/bin/env python3
"""Combustion example: Natural gas + humid air with combustion_equilibrium.

This example demonstrates:
- Creating fuel and air streams
- Mixing streams with enthalpy balance
- combustion_equilibrium: one-step combustion + reforming + WGS equilibrium
- Equivalence ratio sweep from lean to rich
"""
from __future__ import annotations

import numpy as np

import combaero as ca
from combaero.species import SpeciesLocator


def main() -> None:
    sp = SpeciesLocator.from_core()
    n = len(sp.names)

    # =========================================================================
    # Define fuel: Natural gas (90% CH4, 5% C2H6, 2% C3H8, 2% N2, 1% CO2)
    # =========================================================================
    X_fuel = sp.empty()
    X_fuel[sp.indices["CH4"]] = 0.90
    X_fuel[sp.indices["C2H6"]] = 0.05
    X_fuel[sp.indices["C3H8"]] = 0.02
    X_fuel[sp.indices["N2"]] = 0.02
    X_fuel[sp.indices["CO2"]] = 0.01

    fuel = ca.Stream()
    fuel.T = 300.0
    fuel.P = 101325.0
    fuel.X = X_fuel
    fuel.mdot = 1.0

    # =========================================================================
    # Define air: Humid air at 60% RH
    # =========================================================================
    X_air = np.array(ca.humid_air_composition(300.0, 101325.0, 0.6))

    air = ca.Stream()
    air.T = 300.0
    air.P = 101325.0
    air.X = X_air
    air.mdot = 10.0  # kg/s

    # =========================================================================
    # Equivalence ratio sweep
    # =========================================================================
    print("=" * 75)
    print("Combustion Example: Natural Gas + Humid Air")
    print("=" * 75)
    print(f"\nFuel: Natural gas at {fuel.T:.1f} K")
    print(
        f"      CH4: {X_fuel[sp.indices['CH4']]*100:.0f}%, "
        f"C2H6: {X_fuel[sp.indices['C2H6']]*100:.0f}%, "
        f"C3H8: {X_fuel[sp.indices['C3H8']]*100:.0f}%"
    )
    print(f"Air:  Humid air at {air.T:.0f} K, 60% RH, {air.mdot:.1f} kg/s")

    print("\nEquivalence Ratio Sweep")
    print("-" * 75)
    print(
        f"{'phi':>6s} {'mdot_f':>10s} {'T_mix':>10s} {'T_eq':>10s} " f"{'X_CO':>10s} {'X_H2':>10s}"
    )
    print(f"{'[-]':>6s} {'[kg/s]':>10s} {'[K]':>10s} {'[K]':>10s} " f"{'[%]':>10s} {'[%]':>10s}")
    print("-" * 65)

    # Pre-calculate stoichiometric ratio
    Y_fuel = np.array(ca.mole_to_mass(list(X_fuel)))
    Y_air = np.array(ca.mole_to_mass(list(X_air)))
    stoich_ratio = ca.bilger_stoich_mixture_fraction_mass(list(Y_fuel), list(Y_air))

    for phi in np.arange(0.5, 1.25, 0.1):
        # Calculate fuel mass flow for target equivalence ratio
        # Z = phi * Z_st / (1 + (phi - 1) * Z_st)
        # Z = mdot_fuel / (mdot_fuel + mdot_air)
        Z_target = phi * stoich_ratio / (1 + (phi - 1) * stoich_ratio)
        mdot_fuel = Z_target * air.mdot / (1 - Z_target)

        fuel.mdot = mdot_fuel

        # Mix streams
        mixed = ca.mix([fuel, air])

        # One-step: combustion + reforming + WGS equilibrium
        eq = ca.combustion_equilibrium(mixed.T, mixed.X, mixed.P)

        # Get CO and H2 mole fractions
        X_CO = eq.X[sp.indices["CO"]] * 100
        X_H2 = eq.X[sp.indices["H2"]] * 100

        print(
            f"{phi:6.2f} {mdot_fuel:10.4f} {mixed.T:10.2f} {eq.T:10.2f} "
            f"{X_CO:10.4f} {X_H2:10.4f}"
        )

    # =========================================================================
    # Detailed results at stoichiometric
    # =========================================================================
    print("\n" + "=" * 75)
    print("Detailed Results at Stoichiometric (phi = 1.0)")
    print("=" * 75)

    phi = 1.0
    Z_target = phi * stoich_ratio / (1 + (phi - 1) * stoich_ratio)
    fuel.mdot = Z_target * air.mdot / (1 - Z_target)

    mixed = ca.mix([fuel, air])

    # One-step: combustion + reforming + WGS equilibrium
    eq = ca.combustion_equilibrium(mixed.T, mixed.X, mixed.P)

    # For comparison, also show complete combustion (before equilibrium)
    burned = ca.complete_combustion(mixed.T, mixed.X, mixed.P)

    print(f"\nMixed Stream (before combustion):")
    print(f"  T = {mixed.T:.2f} K")
    print(f"  P = {mixed.P:.0f} Pa")
    print(f"  mdot = {mixed.mdot:.4f} kg/s")

    print(f"\nComplete Combustion (CO2 + H2O only):")
    print(f"  T_ad = {burned.T:.2f} K")

    print(f"\nCombustion + Equilibrium (one-step):")
    print(f"  T_eq = {eq.T:.2f} K")

    print(f"\nProduct Composition (mole fractions > 0.1%):")
    for i, name in enumerate(sp.names):
        if eq.X[i] > 0.001:
            print(f"  {name:>8s}: {eq.X[i]*100:6.2f}%")

    # =========================================================================
    # Rich case (phi = 1.2) - Reforming converts hydrocarbons to CO + H2
    # =========================================================================
    print("\n" + "=" * 75)
    print("Rich Case (phi = 1.2) - Reforming + WGS Equilibrium")
    print("=" * 75)

    phi = 1.2
    Z_target = phi * stoich_ratio / (1 + (phi - 1) * stoich_ratio)
    fuel.mdot = Z_target * air.mdot / (1 - Z_target)

    mixed = ca.mix([fuel, air])

    # One-step: combustion + reforming + WGS equilibrium
    eq = ca.combustion_equilibrium(mixed.T, mixed.X, mixed.P)

    # For comparison, also show complete combustion (before equilibrium)
    burned = ca.complete_combustion(mixed.T, mixed.X, mixed.P)

    print(f"\nComplete Combustion (before equilibrium):")
    print(f"  T_ad = {burned.T:.2f} K")
    print(f"  Unburned hydrocarbons:")
    for hc in ["CH4", "C2H6", "C3H8"]:
        if burned.X[sp.indices[hc]] > 1e-6:
            print(f"    {hc}: {burned.X[sp.indices[hc]]*100:.4f}%")

    print(f"\nCombustion + Equilibrium (one-step, hydrocarbons reformed to CO + H2):")
    print(f"  T_eq = {eq.T:.2f} K (dropped {burned.T - eq.T:.1f} K)")
    print(f"  Product composition:")
    for i, name in enumerate(sp.names):
        if eq.X[i] > 0.001:
            print(f"    {name:>8s}: {eq.X[i]*100:6.2f}%")


if __name__ == "__main__":
    main()
