#!/usr/bin/env python3
"""Thermodynamic and transport properties of air, fuel, and combustion products.

Demonstrates:
- Standard dry air and humid air composition helpers
- Temperature sweep of cp, gamma, speed of sound, viscosity, Pr
- Fuel (CH4) properties and stoichiometric oxygen requirement
- Combustion products via mix + complete_combustion
- Inverse solvers: find T from h, s, cp
"""

from __future__ import annotations

import combaero as ca
from combaero.species import SpeciesLocator


def main() -> None:
    sp = SpeciesLocator.from_core()
    P0 = 101325.0

    # -------------------------------------------------------------------------
    # 1. Standard dry air
    # -------------------------------------------------------------------------
    print("=" * 65)
    print("1. Standard Dry Air -- Temperature Sweep")
    print("=" * 65)

    X_air = ca.standard_dry_air_composition()
    print(
        f"\n{'T [K]':>8} {'cp [J/mol·K]':>14} {'gamma':>8} {'a [m/s]':>10} {'mu [µPa·s]':>12} {'Pr':>8}"
    )
    print("-" * 65)
    for T in (300, 600, 900, 1200, 1500, 2000):
        print(
            f"{T:8.0f}"
            f" {ca.cp(T, X_air):14.3f}"
            f" {ca.isentropic_expansion_coefficient(T, X_air):8.4f}"
            f" {ca.speed_of_sound(T, X_air):10.2f}"
            f" {ca.viscosity(T, P0, X_air) * 1e6:12.3f}"
            f" {ca.prandtl(T, P0, X_air):8.4f}"
        )

    # -------------------------------------------------------------------------
    # 2. Humid air
    # -------------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("2. Humid Air (60% RH, 300 K)")
    print("=" * 65)

    X_humid = ca.humid_air_composition(300.0, P0, 0.60)
    idx_h2o = sp.indices["H2O"]
    print(f"\n  H2O mole fraction : {X_humid[idx_h2o] * 100:.3f} %")
    print(f"  cp                : {ca.cp(300.0, X_humid):.3f} J/(mol·K)")
    print(f"  rho               : {ca.density(300.0, P0, X_humid):.4f} kg/m³")

    # -------------------------------------------------------------------------
    # 3. Fuel: pure CH4
    # -------------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("3. Fuel: Pure CH4")
    print("=" * 65)

    X_ch4 = sp.empty()
    X_ch4[sp.indices["CH4"]] = 1.0

    print(f"\n  LHV (mass basis)  : {ca.fuel_lhv_mass(X_ch4) / 1e6:.3f} MJ/kg")
    print(
        f"  O2 required       : {ca.oxygen_required_per_kg_fuel(sp.indices['CH4']):.4f} kg O2/kg fuel"
    )
    print(f"  cp at 300 K       : {ca.cp(300.0, X_ch4):.3f} J/(mol·K)")
    print(f"  cp at 800 K       : {ca.cp(800.0, X_ch4):.3f} J/(mol·K)")

    # -------------------------------------------------------------------------
    # 4. Combustion products: CH4/air at phi=0.8, mix + combust
    # -------------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("4. Combustion Products (CH4/air, phi=0.8)")
    print("=" * 65)

    fuel = ca.Stream()
    fuel.T = 300.0
    fuel.P = P0
    fuel.X = X_ch4

    air = ca.Stream()
    air.T = 300.0
    air.P = P0
    air.X = X_air
    air.mdot = 10.0

    fuel_phi = ca.set_fuel_stream_for_phi(0.8, fuel, air)
    mixed = ca.mix([fuel_phi, air])
    burned = ca.complete_combustion(mixed.T, mixed.X, mixed.P)

    print(f"\n  T_ad              : {burned.T:.1f} K")
    print(f"  cp                : {ca.cp(burned.T, burned.X):.3f} J/(mol·K)")
    print(f"  gamma             : {ca.isentropic_expansion_coefficient(burned.T, burned.X):.4f}")
    print(f"  a                 : {ca.speed_of_sound(burned.T, burned.X):.2f} m/s")
    print(f"  mu                : {ca.viscosity(burned.T, P0, burned.X) * 1e6:.3f} µPa·s")
    print(f"  Pr                : {ca.prandtl(burned.T, P0, burned.X):.4f}")
    print("\n  Major species (> 0.1 %):")
    for i, name in enumerate(sp.names):
        if burned.X[i] > 0.001:
            print(f"    {name:>6s}: {burned.X[i] * 100:.2f} %")

    # -------------------------------------------------------------------------
    # 5. Inverse solvers: recover T from h, s, cp
    # -------------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("5. Inverse Solvers (air, T_ref = 1200 K)")
    print("=" * 65)

    T_ref = 1200.0
    h_ref = ca.h(T_ref, X_air)
    s_ref = ca.s(T_ref, X_air, P0)
    cp_ref = ca.cp(T_ref, X_air)

    T_from_h = ca.calc_T_from_h(h_ref, X_air)
    T_from_s = ca.calc_T_from_s(s_ref, P0, X_air)
    T_from_cp = ca.calc_T_from_cp(cp_ref, X_air)

    print(f"\n  T_ref             : {T_ref:.2f} K")
    print(f"  T from h          : {T_from_h:.6f} K  (err {abs(T_from_h - T_ref):.2e} K)")
    print(f"  T from s          : {T_from_s:.6f} K  (err {abs(T_from_s - T_ref):.2e} K)")
    print(f"  T from cp         : {T_from_cp:.6f} K  (err {abs(T_from_cp - T_ref):.2e} K)")


if __name__ == "__main__":
    main()
