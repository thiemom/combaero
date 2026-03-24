#!/usr/bin/env python3
"""Thermodynamic and transport properties of air, fuel, and combustion products.

Demonstrates:
- Standard dry air and humid air composition helpers
- Temperature sweep of cp, gamma, speed of sound, viscosity, Pr
- Fuel (CH4) properties and stoichiometric oxygen requirement
- Combustion products via mix + complete_combustion
- Inverse solvers: find T from h, s, cp
- Unit introspection: query expected API units
"""

from __future__ import annotations

import combaero as cb


def main() -> None:

    P0 = 101325.0

    # -------------------------------------------------------------------------
    # 1. Standard dry air
    # -------------------------------------------------------------------------
    print("=" * 65)
    print("1. Standard Dry Air -- Temperature Sweep")
    print("=" * 65)

    X_air = cb.species.dry_air()
    print(
        f"\n{'T [K]':>8} {'cp [J/mol K]':>14} {'gamma':>8} {'a [m/s]':>10} {'mu [uPa s]':>12} {'Pr':>8}"
    )
    print("-" * 65)
    for T in (300, 600, 900, 1200, 1500, 2000):
        print(
            f"{T:8.0f}"
            f" {cb.cp(T, X_air):14.3f}"
            f" {cb.isentropic_expansion_coefficient(T, X_air):8.4f}"
            f" {cb.speed_of_sound(T, X_air):10.2f}"
            f" {cb.viscosity(T, P0, X_air) * 1e6:12.3f}"
            f" {cb.prandtl(T, P0, X_air):8.4f}"
        )

    # -------------------------------------------------------------------------
    # 2. Humid air
    # -------------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("2. Humid Air (60% RH, 300 K)")
    print("=" * 65)

    X_humid = cb.species.humid_air(300.0, P0, 0.60)
    idx_h2o = cb.species.indices["H2O"]
    print(f"\n  H2O mole fraction : {X_humid[idx_h2o] * 100:.3f} %")
    print(f"  cp                : {cb.cp(300.0, X_humid):.3f} J/(mol K)")
    print(f"  rho               : {cb.density(300.0, P0, X_humid):.4f} kg/m3")

    # -------------------------------------------------------------------------
    # 3. Fuel: pure CH4
    # -------------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("3. Fuel: Pure CH4")
    print("=" * 65)

    X_ch4 = cb.species.empty()
    X_ch4[cb.species.indices["CH4"]] = 1.0

    print(f"\n  LHV (mass basis)  : {cb.fuel_lhv_mass(X_ch4) / 1e6:.3f} MJ/kg")
    print(
        f"  O2 required       : {cb.oxygen_required_per_kg_fuel(cb.species.indices['CH4']):.4f} kg O2/kg fuel"
    )
    print(f"  cp at 300 K       : {cb.cp(300.0, X_ch4):.3f} J/(mol K)")
    print(f"  cp at 800 K       : {cb.cp(800.0, X_ch4):.3f} J/(mol K)")

    # -------------------------------------------------------------------------
    # 4. Combustion products: CH4/air at phi=0.8, mix + combust
    # -------------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("4. Combustion Products (CH4/air, phi=0.8)")
    print("=" * 65)

    fuel = cb.Stream()
    fuel.state.TPX = 300.0, P0, X_ch4

    air = cb.Stream()
    air.state.TPX = 300.0, P0, X_air
    air.mdot = 10.0

    fuel_phi = cb.set_fuel_stream_for_phi(0.8, fuel, air)
    mixed = cb.mix([fuel_phi, air])
    burned = cb.complete_combustion(mixed.T, mixed.X, mixed.P)

    print(f"\n  T_ad              : {burned.T:.1f} K")
    print(f"  cp                : {cb.cp(burned.T, burned.X):.3f} J/(mol K)")
    print(f"  gamma             : {cb.isentropic_expansion_coefficient(burned.T, burned.X):.4f}")
    print(f"  a                 : {cb.speed_of_sound(burned.T, burned.X):.2f} m/s")
    print(f"  mu                : {cb.viscosity(burned.T, P0, burned.X) * 1e6:.3f} uPa s")
    print(f"  Pr                : {cb.prandtl(burned.T, P0, burned.X):.4f}")
    print("\n  Major species (> 0.1 %):")
    for i, name in enumerate(cb.species.names):
        if burned.X[i] > 0.001:
            print(f"    {name:>6s}: {burned.X[i] * 100:.2f} %")

    # -------------------------------------------------------------------------
    # 5. Inverse solvers: recover T from h, s, cp
    # -------------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("5. Inverse Solvers (air, T_ref = 1200 K)")
    print("=" * 65)

    T_ref = 1200.0
    h_ref = cb.h(T_ref, X_air)
    s_ref = cb.s(T_ref, X_air, P0)
    cp_ref = cb.cp(T_ref, X_air)

    T_from_h = cb.calc_T_from_h(h_ref, X_air)
    T_from_s = cb.calc_T_from_s(s_ref, P0, X_air)
    T_from_cp = cb.calc_T_from_cp(cp_ref, X_air)

    print(f"\n  T_ref             : {T_ref:.2f} K")
    print(f"  T from h          : {T_from_h:.6f} K  (err {abs(T_from_h - T_ref):.2e} K)")
    print(f"  T from s          : {T_from_s:.6f} K  (err {abs(T_from_s - T_ref):.2e} K)")
    print(f"  T from cp         : {T_from_cp:.6f} K  (err {abs(T_from_cp - T_ref):.2e} K)")

    # -------------------------------------------------------------------------
    # 6. Unit Introspection API
    # -------------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("6. Unit Introspection API")
    print("=" * 65)

    print("\n  Every API parameter and output is explicitly mapped:")
    print(f"  cb.h outputs                      : {cb.output_units('h')} (molar basis)")
    print(f"  cb.density outputs                : {cb.output_units('density')}")
    print(f"  cb.dynamic_viscosity outputs      : {cb.output_units('viscosity')}")
    print(f"  cb.fuel_lhv_mass outputs          : {cb.output_units('fuel_lhv_mass')}")
    print()
    print("  We can also query function input signatures:")
    print(f"  cb.cp required inputs             : {cb.input_units('cp')}")
    print(f"  cb.species.humid_air inputs   : {cb.input_units('humid_air_composition')}")


if __name__ == "__main__":
    main()
