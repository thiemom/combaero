#!/usr/bin/env python3
"""Reactor design with geometry utilities.

This example demonstrates:
- Hydraulic diameter calculations (generic, rectangular, annular)
- Residence time calculations
- Space velocity
- Reactor sizing and design
- Combustor design applications
"""

from __future__ import annotations

import numpy as np

import combaero as ca
from combaero.species import SpeciesLocator


def main() -> None:
    print("=" * 80)
    print("REACTOR DESIGN WITH GEOMETRY UTILITIES")
    print("=" * 80)

    # =========================================================================
    # Example 1: Hydraulic diameter calculations
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 1: Hydraulic Diameter for Different Geometries")
    print("=" * 80)

    # Circular pipe
    D_pipe = 0.1  # Diameter [m]
    A_pipe = ca.pipe_area(D_pipe)  # Use pipe_area helper
    P_pipe = np.pi * D_pipe
    Dh_pipe = ca.hydraulic_diameter(A_pipe, P_pipe)

    print("\nCircular pipe:")
    print(f"  Diameter:           D = {D_pipe * 1000:.1f} mm")
    print(f"  Area:               A = {A_pipe * 1e6:.2f} mm2")
    print(f"  Perimeter:          P = {P_pipe * 1000:.2f} mm")
    print(f"  Hydraulic diameter: Dh = {Dh_pipe * 1000:.2f} mm")
    print("  (Should equal D for circular pipe)")

    # Square duct
    a = 0.1  # Side length [m]
    Dh_square = ca.hydraulic_diameter_rect(a, a)

    print("\nSquare duct:")
    print(f"  Side length:        a = {a * 1000:.1f} mm")
    print(f"  Hydraulic diameter: Dh = {Dh_square * 1000:.2f} mm")
    print("  (Should equal a for square)")

    # Rectangular duct
    a_rect = 0.2  # Width [m]
    b_rect = 0.1  # Height [m]
    Dh_rect = ca.hydraulic_diameter_rect(a_rect, b_rect)

    print("\nRectangular duct:")
    print(f"  Width:              a = {a_rect * 1000:.1f} mm")
    print(f"  Height:             b = {b_rect * 1000:.1f} mm")
    print(f"  Hydraulic diameter: Dh = {Dh_rect * 1000:.2f} mm")
    print(f"  Dh = 2ab/(a+b) = {2 * a_rect * b_rect / (a_rect + b_rect) * 1000:.2f} mm")

    # Annular duct
    D_outer = 0.15  # Outer diameter [m]
    D_inner = 0.10  # Inner diameter [m]
    Dh_annulus = ca.hydraulic_diameter_annulus(D_outer, D_inner)

    print("\nAnnular duct:")
    print(f"  Outer diameter:     D_o = {D_outer * 1000:.1f} mm")
    print(f"  Inner diameter:     D_i = {D_inner * 1000:.1f} mm")
    print(f"  Hydraulic diameter: Dh = {Dh_annulus * 1000:.2f} mm")
    print(f"  Dh = D_o - D_i = {(D_outer - D_inner) * 1000:.2f} mm")

    # =========================================================================
    # Example 2: Residence time in a reactor
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 2: Residence Time in a Tubular Reactor")
    print("=" * 80)

    # Reactor geometry
    D = 0.2  # Diameter [m]
    L = 2.0  # Length [m]
    V = ca.pipe_volume(D, L)  # Use pipe_volume helper

    # Flow conditions
    Q = 0.01  # Volumetric flow rate [m3/s]

    # Calculate residence time
    tau = ca.residence_time(V, Q)
    SV = ca.space_velocity(Q, V)

    print("\nReactor geometry:")
    print(f"  Diameter:        D = {D * 1000:.0f} mm")
    print(f"  Length:          L = {L:.1f} m")
    print(f"  Volume:          V = {V * 1000:.2f} L")

    print("\nFlow conditions:")
    print(f"  Volumetric flow: Q = {Q * 1000:.1f} L/s")

    print("\nResults:")
    print(f"  Residence time:  tau = {tau:.2f} s = {tau / 60:.3f} min")
    print(f"  Space velocity:  SV = {SV:.4f} s^-^1 = {SV * 3600:.1f} h^-^1")

    # =========================================================================
    # Example 3: Residence time from mass flow rate
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 3: Residence Time from Mass Flow Rate")
    print("=" * 80)

    # Reactor volume
    V = 0.05  # Volume [m3]

    # Air flow
    mdot = 0.5  # Mass flow rate [kg/s]
    T = 300.0  # Temperature [K]
    P = 101325.0  # Pressure [Pa]
    X_air = ca.standard_dry_air_composition()
    rho = ca.density(T, P, X_air)

    # Calculate residence time
    tau_mdot = ca.residence_time_mdot(V, mdot, rho)

    # Also calculate from volumetric flow for comparison
    Q = mdot / rho
    tau_Q = ca.residence_time(V, Q)

    print(f"\nReactor volume:  V = {V * 1000:.1f} L")
    print("\nAir flow:")
    print(f"  Mass flow rate:  mdot = {mdot:.2f} kg/s")
    print(f"  Temperature:     T = {T:.0f} K")
    print(f"  Density:         rho = {rho:.3f} kg/m3")
    print(f"  Volumetric flow: Q = {Q * 1000:.2f} L/s")

    print("\nResidence time:")
    print(f"  From mass flow:  tau = {tau_mdot:.3f} s")
    print(f"  From vol. flow:  tau = {tau_Q:.3f} s")
    print("  (Should be identical)")

    # =========================================================================
    # Example 4: Combustor design
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 4: Gas Turbine Combustor Design")
    print("=" * 80)

    # Design requirements
    tau_target = 0.015  # Target residence time [s] = 15 ms
    mdot_air = 10.0  # Air mass flow [kg/s]
    P_comb = 15e5  # Pressure [Pa] = 15 bar
    phi_comb = 0.6  # Lean primary zone

    # Build burned gas composition for realistic density
    sp = SpeciesLocator.from_core()
    X_air_c = ca.standard_dry_air_composition()
    X_ch4 = sp.empty()
    X_ch4[sp.indices["CH4"]] = 1.0

    air_c = ca.Stream()
    air_c.T, air_c.P, air_c.X, air_c.mdot = 700.0, P_comb, X_air_c, mdot_air

    fuel_c = ca.Stream()
    fuel_c.T, fuel_c.P, fuel_c.X = 300.0, P_comb, X_ch4

    fuel_c = ca.set_fuel_stream_for_phi(phi_comb, fuel_c, air_c)
    mixed_c = ca.mix([fuel_c, air_c])
    burned_c = ca.complete_combustion(mixed_c.T, mixed_c.X, mixed_c.P)

    mdot_total = mixed_c.mdot
    rho = ca.density(burned_c.T, P_comb, burned_c.X)

    # Calculate required volume
    V_required = tau_target * mdot_total / rho

    # Annular combustor geometry
    D_outer = 0.5  # Outer diameter [m]
    D_inner = 0.35  # Inner diameter [m]
    A_annular = np.pi * ((D_outer / 2) ** 2 - (D_inner / 2) ** 2)

    # Required length
    L_required = V_required / A_annular

    # Hydraulic diameter
    Dh = ca.hydraulic_diameter_annulus(D_outer, D_inner)

    print("\nDesign requirements:")
    print(f"  Target residence time: tau = {tau_target * 1000:.1f} ms")
    print(f"  Air mass flow:         mdot_air = {mdot_air:.1f} kg/s")
    print(f"  Fuel mass flow:        mdot_fuel = {fuel_c.mdot:.3f} kg/s  (phi={phi_comb:.1f})")
    print(f"  Total mass flow:       mdot = {mdot_total:.3f} kg/s")
    print(f"  Operating pressure:    P = {P_comb / 1e5:.0f} bar")
    print(f"  T_ad (burned):         T = {burned_c.T:.0f} K")

    print("\nBurned gas properties:")
    print(f"  Density:               rho = {rho:.4f} kg/m3  (burned gas at T_ad)")

    print("\nAnnular combustor:")
    print(f"  Outer diameter:        D_o = {D_outer * 1000:.0f} mm")
    print(f"  Inner diameter:        D_i = {D_inner * 1000:.0f} mm")
    print(f"  Hydraulic diameter:    Dh = {Dh * 1000:.0f} mm")
    print(f"  Annular area:          A = {A_annular * 1e4:.2f} cm2")

    print("\nResults:")
    print(f"  Required volume:       V = {V_required * 1000:.2f} L")
    print(f"  Required length:       L = {L_required * 1000:.0f} mm")

    # Verify
    tau_actual = ca.residence_time_mdot(V_required, mdot_total, rho)
    print(f"  Actual residence time: tau = {tau_actual * 1000:.2f} ms OK")

    # =========================================================================
    # Example 5: Reactor sizing for different flow rates
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 5: Reactor Sizing for Different Flow Rates")
    print("=" * 80)

    # Fixed reactor volume
    V = 0.1  # Volume [m3]

    # Air properties
    T = 300.0
    P = 101325.0
    X_air = ca.standard_dry_air_composition()
    rho = ca.density(T, P, X_air)

    print(f"\nFixed reactor volume: V = {V * 1000:.1f} L")
    print(f"Air at {T:.0f} K, {P / 1000:.0f} kPa (rho = {rho:.3f} kg/m3)")

    print(f"\n{'mdot [kg/s]':>12s}  {'Q [L/s]':>12s}  {'tau [s]':>12s}  {'SV [h^-^1]':>12s}")
    print("-" * 52)

    for mdot in [0.1, 0.5, 1.0, 2.0, 5.0]:
        Q = mdot / rho
        tau = ca.residence_time_mdot(V, mdot, rho)
        SV = ca.space_velocity(Q, V)
        print(f"{mdot:12.2f}  {Q * 1000:12.2f}  {tau:12.3f}  {SV * 3600:12.1f}")

    # =========================================================================
    # Example 6: Damkohler number (reaction vs residence time)
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 6: Damkohler Number Analysis")
    print("=" * 80)

    # Reactor parameters
    V = 0.05  # Volume [m3]
    Q = 0.01  # Flow rate [m3/s]
    tau = ca.residence_time(V, Q)

    # Reaction time scale (example: first-order reaction)
    k = 10.0  # Reaction rate constant [1/s]
    tau_reaction = 1.0 / k  # Reaction time scale [s]

    # Damkohler number
    Da = tau / tau_reaction

    print("\nReactor:")
    print(f"  Volume:          V = {V * 1000:.1f} L")
    print(f"  Flow rate:       Q = {Q * 1000:.1f} L/s")
    print(f"  Residence time:  tau = {tau:.2f} s")

    print("\nReaction:")
    print(f"  Rate constant:   k = {k:.1f} s^-^1")
    print(f"  Reaction time:   tau_rxn = {tau_reaction:.3f} s")

    print("\nDamkohler number:")
    print(f"  Da = tau/tau_rxn = {Da:.2f}")

    if Da > 10:
        print("  -> Reaction-limited (fast flow, slow reaction)")
    elif Da < 0.1:
        print("  -> Mixing-limited (slow flow, fast reaction)")
    else:
        print("  -> Intermediate regime")

    # =========================================================================
    # Example 7: Pipe flow residence time
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 7: Residence Time in Pipe Flow")
    print("=" * 80)

    # Pipe geometry
    D = 0.05  # Diameter [m]
    L = 10.0  # Length [m]
    A = np.pi * (D / 2) ** 2
    V = A * L

    # Flow velocity
    v = 5.0  # Velocity [m/s]
    Q = A * v

    # Residence time
    tau = ca.residence_time(V, Q)

    # Simple calculation: tau = L/v
    tau_simple = L / v

    print("\nPipe:")
    print(f"  Diameter:        D = {D * 1000:.0f} mm")
    print(f"  Length:          L = {L:.1f} m")
    print(f"  Volume:          V = {V * 1000:.3f} L")

    print("\nFlow:")
    print(f"  Velocity:        v = {v:.1f} m/s")
    print(f"  Volumetric flow: Q = {Q * 1000:.3f} L/s")

    print("\nResidence time:")
    print(f"  tau = V/Q = {tau:.3f} s")
    print(f"  tau = L/v = {tau_simple:.3f} s")
    print("  (Should be identical for constant area)")


if __name__ == "__main__":
    main()
