#!/usr/bin/env python3
"""De Laval (converging-diverging) nozzle example.

Classical supersonic nozzle problem demonstrating:
- Subsonic flow in converging section
- Sonic conditions at throat (M=1)
- Supersonic expansion in diverging section

Chamber gas: real H2/air combustion products at phi=1.0.
"""

import matplotlib.pyplot as plt
import numpy as np
from plot_utils import show_or_save

import combaero as cb
from combaero.species import SpeciesLocator


def main() -> None:
    # =========================================================================
    # Nozzle geometry (small rocket motor)
    # =========================================================================
    D_inlet = 0.05  # m (inlet diameter)
    D_throat = 0.025  # m (throat diameter)
    D_exit = 0.04  # m (exit diameter)

    A_inlet = np.pi * (D_inlet / 2) ** 2  # 19.6 cm^2
    A_throat = np.pi * (D_throat / 2) ** 2  # 4.9 cm^2
    A_exit = np.pi * (D_exit / 2) ** 2  # 12.6 cm^2

    x_throat = 0.05  # m (throat at 5 cm)
    x_exit = 0.12  # m (exit at 12 cm)

    # Area ratio
    AR = A_exit / A_throat
    print("=" * 60)
    print("De Laval Nozzle Analysis")
    print("=" * 60)
    print("\nGeometry:")
    print(f"  Inlet diameter:  {D_inlet * 100:.1f} cm  (A = {A_inlet * 1e4:.2f} cm^2)")
    print(f"  Throat diameter: {D_throat * 100:.1f} cm  (A = {A_throat * 1e4:.2f} cm^2)")
    print(f"  Exit diameter:   {D_exit * 100:.1f} cm  (A = {A_exit * 1e4:.2f} cm^2)")
    print(f"  Area ratio A_exit/A_throat: {AR:.2f}")

    # =========================================================================
    # Chamber gas: H2/air combustion products at phi=1.0, P0=20 bar
    # =========================================================================
    sp = SpeciesLocator.from_core()
    P0 = 2_000_000.0  # Pa  (20 bar chamber pressure)

    X_air = cb.standard_dry_air_composition()
    X_h2 = sp.empty()
    X_h2[sp.indices["H2"]] = 1.0

    fuel = cb.Stream()
    fuel.T, fuel.P, fuel.X = 300.0, P0, X_h2

    air = cb.Stream()
    air.T, air.P, air.X, air.mdot = 300.0, P0, X_air, 10.0

    fuel_phi = cb.set_fuel_stream_for_phi(1.0, fuel, air)
    mixed = cb.mix([fuel_phi, air])
    burned = cb.complete_combustion(mixed.T, mixed.X, mixed.P)

    T0 = burned.T  # adiabatic flame temperature = stagnation T
    X = burned.X  # combustion products composition

    print("\nChamber gas (H2/air, phi=1.0):")
    print(f"  T0 = {T0:.0f} K  (adiabatic flame temperature)")
    print(f"  P0 = {P0 / 1e5:.1f} bar")

    # Get gas properties at stagnation
    gamma = cb.isentropic_expansion_coefficient(T0, X)
    a0 = cb.speed_of_sound(T0, X)
    print(f"  gamma = {gamma:.3f}")
    print(f"  a0 = {a0:.1f} m/s")

    # Critical pressure ratio
    P_ratio_crit = cb.critical_pressure_ratio(T0, P0, X)
    P_crit = P0 * P_ratio_crit
    print("\nCritical (sonic) conditions:")
    print(f"  P*/P0 = {P_ratio_crit:.4f}")
    print(f"  P* = {P_crit / 1e5:.2f} bar")

    # =========================================================================
    # Case 1: Fully expanded supersonic flow (design condition)
    # =========================================================================
    # For supersonic exit, P_exit must be low enough
    P_exit_design = 100000.0  # Pa (1 bar, atmospheric)

    print("\n" + "-" * 60)
    print("Case 1: Design condition (supersonic exit)")
    print("-" * 60)

    sol_design = cb.nozzle_cd(
        T0, P0, P_exit_design, A_inlet, A_throat, A_exit, x_throat, x_exit, X, n_stations=100
    )

    print(f"Exit pressure: {P_exit_design / 1e5:.2f} bar")
    print(f"Choked: {sol_design.choked}")
    print(f"Mass flow: {sol_design.mdot:.4f} kg/s")

    # Find exit conditions from profile
    exit_station = sol_design.profile[-1]
    throat_idx = np.argmin([abs(st.x - x_throat) for st in sol_design.profile])
    throat_station = sol_design.profile[throat_idx]

    print(f"\nThroat conditions (x = {throat_station.x:.3f} m):")
    print(f"  M = {throat_station.M:.3f}")
    print(f"  T = {throat_station.T:.1f} K")
    print(f"  P = {throat_station.P / 1e5:.2f} bar")
    print(f"  u = {throat_station.u:.1f} m/s")

    print(f"\nExit conditions (x = {exit_station.x:.3f} m):")
    print(f"  M = {exit_station.M:.3f}")
    print(f"  T = {exit_station.T:.1f} K")
    print(f"  P = {exit_station.P / 1e5:.2f} bar")
    print(f"  u = {exit_station.u:.1f} m/s")

    # Thrust at sea level: P_amb = 1 bar
    thrust_sl = cb.nozzle_thrust(sol_design, P_exit_design)
    print("\nThrust performance (sea level, P_amb=1 bar):")
    print(f"  Thrust: {thrust_sl.thrust:.1f} N")
    print(f"  Specific impulse: {thrust_sl.specific_impulse:.1f} s")
    print(f"  Thrust coefficient: {thrust_sl.thrust_coefficient:.3f}")

    # Vacuum thrust: reuse same nozzle solution, just change P_amb
    thrust_vac = cb.nozzle_thrust(sol_design, 0.0)
    print("\nVacuum performance (P_amb = 0):")
    print(f"  Thrust: {thrust_vac.thrust:.1f} N")
    print(f"  Specific impulse: {thrust_vac.specific_impulse:.1f} s")

    # =========================================================================
    # Case 2: Subsonic flow (high back pressure)
    # =========================================================================
    P_exit_subsonic = 1500000.0  # Pa (15 bar)

    print("\n" + "-" * 60)
    print("Case 2: Subsonic flow (high back pressure)")
    print("-" * 60)

    sol_subsonic = cb.nozzle_cd(
        T0, P0, P_exit_subsonic, A_inlet, A_throat, A_exit, x_throat, x_exit, X, n_stations=100
    )

    print(f"Exit pressure: {P_exit_subsonic / 1e5:.2f} bar")
    print(f"Choked: {sol_subsonic.choked}")
    print(f"Mass flow: {sol_subsonic.mdot:.4f} kg/s")

    exit_sub = sol_subsonic.profile[-1]
    print(f"Exit Mach: {exit_sub.M:.3f}")

    # =========================================================================
    # Plot results
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Extract profile data for design case
    x = np.array([st.x for st in sol_design.profile]) * 100  # cm
    A = np.array([st.A for st in sol_design.profile]) * 1e4  # cm^2
    M = np.array([st.M for st in sol_design.profile])
    T = np.array([st.T for st in sol_design.profile])
    P = np.array([st.P for st in sol_design.profile]) / 1e5  # bar

    # Area profile
    ax = axes[0, 0]
    ax.plot(x, A, "b-", linewidth=2)
    ax.axvline(x_throat * 100, color="r", linestyle="--", label="Throat")
    ax.set_xlabel("Axial position [cm]")
    ax.set_ylabel("Area [cm^2]")
    ax.set_title("Nozzle Area Profile")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Mach number
    ax = axes[0, 1]
    ax.plot(x, M, "g-", linewidth=2)
    ax.axhline(1.0, color="r", linestyle="--", label="M = 1 (sonic)")
    ax.axvline(x_throat * 100, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Axial position [cm]")
    ax.set_ylabel("Mach number [-]")
    ax.set_title("Mach Number Distribution")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Temperature
    ax = axes[1, 0]
    ax.plot(x, T, "r-", linewidth=2)
    ax.axhline(T0, color="gray", linestyle="--", label=f"T0 = {T0:.0f} K")
    ax.axvline(x_throat * 100, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Axial position [cm]")
    ax.set_ylabel("Temperature [K]")
    ax.set_title("Temperature Distribution")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Pressure
    ax = axes[1, 1]
    ax.plot(x, P, "m-", linewidth=2)
    ax.axhline(P0 / 1e5, color="gray", linestyle="--", label=f"P0 = {P0 / 1e5:.0f} bar")
    ax.axvline(x_throat * 100, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Axial position [cm]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("Pressure Distribution")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle("De Laval Nozzle - Design Condition (Supersonic Exit)", fontsize=14)
    plt.tight_layout()
    show_or_save(fig, "de_laval_nozzle.png")


if __name__ == "__main__":
    main()
