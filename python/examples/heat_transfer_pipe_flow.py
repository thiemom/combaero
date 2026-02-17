#!/usr/bin/env python3
"""Heat transfer in pipe flow.

This example demonstrates:
- Nusselt number correlations (Dittus-Boelter, Gnielinski, Sieder-Tate)
- Heat transfer coefficient calculations
- Pipe heating/cooling design
- LMTD for heat exchangers
- Integration with friction factors
"""

from __future__ import annotations

import numpy as np

import combaero as ca


def main() -> None:
    print("=" * 80)
    print("HEAT TRANSFER IN PIPE FLOW")
    print("=" * 80)

    # =========================================================================
    # Example 1: Nusselt number correlations
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 1: Nusselt Number Correlations")
    print("=" * 80)

    Re = 5e4  # Reynolds number
    Pr = 0.7  # Prandtl number (air)

    print(f"\nReynolds number: Re = {Re:.0e}")
    print(f"Prandtl number:  Pr = {Pr:.2f}")

    # Compare correlations
    Nu_db_heat = ca.nusselt_dittus_boelter(Re, Pr, heating=True)
    Nu_db_cool = ca.nusselt_dittus_boelter(Re, Pr, heating=False)
    Nu_gn = ca.nusselt_gnielinski(Re, Pr)
    Nu_st = ca.nusselt_sieder_tate(Re, Pr, mu_ratio=1.0)

    print("\nNusselt numbers:")
    print(f"  Dittus-Boelter (heating): Nu = {Nu_db_heat:.2f}")
    print(f"  Dittus-Boelter (cooling): Nu = {Nu_db_cool:.2f}")
    print(f"  Gnielinski:               Nu = {Nu_gn:.2f}")
    print(f"  Sieder-Tate:              Nu = {Nu_st:.2f}")

    # =========================================================================
    # Example 2: Heat transfer coefficient for air in a pipe
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 2: Heat Transfer Coefficient - Air in Pipe")
    print("=" * 80)

    # Pipe and flow conditions
    D = 0.05  # Diameter [m]
    v = 10.0  # Velocity [m/s]
    T = 300.0  # Air temperature [K]
    P = 101325.0  # Pressure [Pa]

    # Air properties
    X_air = ca.standard_dry_air_composition()
    rho = ca.density(T, P, X_air)
    mu = ca.viscosity(T, P, X_air)
    k = ca.thermal_conductivity(T, P, X_air)
    cp = ca.cp(T, X_air)
    Pr = ca.prandtl(T, P, X_air)

    # Calculate Reynolds number
    Re = rho * v * D / mu

    # Calculate Nusselt number
    Nu = ca.nusselt_dittus_boelter(Re, Pr, heating=True)

    # Calculate heat transfer coefficient
    h = ca.htc_from_nusselt(Nu, k, D)

    print(f"\nPipe diameter:   D = {D * 1000:.1f} mm")
    print(f"Air velocity:    v = {v:.1f} m/s")
    print(f"Air temperature: T = {T:.0f} K")

    print("\nAir properties:")
    print(f"  Density:              ρ = {rho:.3f} kg/m3")
    print(f"  Viscosity:            μ = {mu * 1e6:.2f} μPa*s")
    print(f"  Thermal conductivity: k = {k * 1000:.2f} mW/(m*K)")
    print(f"  Specific heat:        cp = {cp:.1f} J/(mol*K)")
    print(f"  Prandtl number:       Pr = {Pr:.3f}")

    print("\nResults:")
    print(f"  Reynolds number:      Re = {Re:.2e}")
    print(f"  Nusselt number:       Nu = {Nu:.2f}")
    print(f"  Heat transfer coeff:  h = {h:.2f} W/(m2*K)")

    # =========================================================================
    # Example 3: Pipe heating design
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 3: Pipe Heating Design")
    print("=" * 80)

    # Design problem: Heat air from 300K to 400K
    T_in = 300.0  # Inlet temperature [K]
    T_out = 400.0  # Outlet temperature [K]
    T_wall = 450.0  # Wall temperature [K]
    mdot = 0.1  # Mass flow rate [kg/s]
    D = 0.05  # Pipe diameter [m]

    # Air properties (average temperature)
    T_avg = (T_in + T_out) / 2
    X_air = ca.standard_dry_air_composition()
    rho = ca.density(T_avg, P, X_air)
    mu = ca.viscosity(T_avg, P, X_air)
    k = ca.thermal_conductivity(T_avg, P, X_air)
    cp_mass = ca.cp_mass(T_avg, X_air)  # Use cp_mass convenience function
    Pr = ca.prandtl(T_avg, P, X_air)

    # Flow velocity
    A = ca.pipe_area(D)  # Use pipe_area helper
    v = mdot / (rho * A)
    Re = rho * v * D / mu

    # Heat transfer coefficient
    Nu = ca.nusselt_dittus_boelter(Re, Pr, heating=True)
    h = ca.htc_from_nusselt(Nu, k, D)

    # Heat transfer rate required
    Q = mdot * cp_mass * (T_out - T_in)

    # LMTD
    dT1 = T_wall - T_in  # Inlet end
    dT2 = T_wall - T_out  # Outlet end
    LMTD = ca.lmtd(dT1, dT2)

    # Required length: Q = h * A_surface * LMTD
    # A_surface = π * D * L
    L = Q / (h * np.pi * D * LMTD)

    print("\nDesign requirements:")
    print(f"  Heat air from {T_in:.0f} K to {T_out:.0f} K")
    print(f"  Mass flow rate:  ṁ = {mdot:.3f} kg/s")
    print(f"  Wall temperature: T_wall = {T_wall:.0f} K")
    print(f"  Pipe diameter:    D = {D * 1000:.1f} mm")

    print(f"\nFlow conditions (at T_avg = {T_avg:.0f} K):")
    print(f"  Velocity:         v = {v:.2f} m/s")
    print(f"  Reynolds number:  Re = {Re:.2e}")
    print(f"  Prandtl number:   Pr = {Pr:.3f}")

    print("\nHeat transfer:")
    print(f"  Nusselt number:   Nu = {Nu:.2f}")
    print(f"  HTC:              h = {h:.2f} W/(m2*K)")
    print(f"  Heat duty:        Q = {Q / 1000:.2f} kW")
    print(f"  LMTD:             LMTD = {LMTD:.2f} K")

    print("\nResult:")
    print(f"  Required length:  L = {L:.2f} m")

    # =========================================================================
    # Example 4: Heat exchanger LMTD
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 4: Heat Exchanger LMTD")
    print("=" * 80)

    # Counter-flow heat exchanger
    print("\nCounter-flow heat exchanger:")
    T_hot_in = 373.15  # 100degC
    T_hot_out = 333.15  # 60degC
    T_cold_in = 293.15  # 20degC
    T_cold_out = 323.15  # 50degC

    # Use lmtd_counterflow convenience function
    LMTD_counter = ca.lmtd_counterflow(T_hot_in, T_hot_out, T_cold_in, T_cold_out)
    dT1 = T_hot_in - T_cold_out  # For display
    dT2 = T_hot_out - T_cold_in

    print(f"  Hot side:  {T_hot_in - 273.15:.0f}degC → {T_hot_out - 273.15:.0f}degC")
    print(f"  Cold side: {T_cold_in - 273.15:.0f}degC → {T_cold_out - 273.15:.0f}degC")
    print(f"  DeltaT1 = {dT1:.2f} K (hot in - cold out)")
    print(f"  DeltaT2 = {dT2:.2f} K (hot out - cold in)")
    print(f"  LMTD = {LMTD_counter:.2f} K")

    # Parallel-flow heat exchanger
    print("\nParallel-flow heat exchanger (same temperatures):")
    # Use lmtd_parallelflow convenience function
    LMTD_parallel = ca.lmtd_parallelflow(T_hot_in, T_hot_out, T_cold_in, T_cold_out)
    dT1_par = T_hot_in - T_cold_in  # For display
    dT2_par = T_hot_out - T_cold_out

    print(f"  DeltaT1 = {dT1_par:.2f} K (both inlets)")
    print(f"  DeltaT2 = {dT2_par:.2f} K (both outlets)")
    print(f"  LMTD = {LMTD_parallel:.2f} K")

    print("\nCounter-flow is more efficient:")
    print(f"  LMTD_counter / LMTD_parallel = {LMTD_counter / LMTD_parallel:.3f}")

    # =========================================================================
    # Example 5: Viscosity correction (Sieder-Tate)
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 5: Viscosity Correction for Heating")
    print("=" * 80)

    Re = 5e4
    Pr = 0.7

    # Heating: wall cooler than bulk, so μ_wall > μ_bulk
    mu_ratios = [0.6, 0.8, 1.0, 1.2, 1.4]

    print(f"\nReynolds number: Re = {Re:.0e}")
    print(f"Prandtl number:  Pr = {Pr:.2f}")
    print(f"\n{'μ_bulk/μ_wall':>15s}  {'Nu (Sieder-Tate)':>20s}  {'Effect':>10s}")
    print("-" * 50)

    Nu_base = ca.nusselt_sieder_tate(Re, Pr, mu_ratio=1.0)

    for mu_ratio in mu_ratios:
        Nu = ca.nusselt_sieder_tate(Re, Pr, mu_ratio=mu_ratio)
        effect = (Nu - Nu_base) / Nu_base * 100
        print(f"{mu_ratio:15.2f}  {Nu:20.2f}  {effect:+9.1f}%")

    print("\nNote: μ_bulk/μ_wall < 1 → heating (wall hotter)")
    print("      μ_bulk/μ_wall > 1 → cooling (wall cooler)")

    # =========================================================================
    # Example 6: Combined friction and heat transfer
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 6: Gnielinski with Friction Factor")
    print("=" * 80)

    Re = 5e4
    Pr = 0.7
    D = 0.1  # Pipe diameter [m]
    eps = ca.pipe_roughness("commercial_steel")  # Use pipe_roughness database
    e_D = eps / D

    # Get friction factor
    f = ca.friction_haaland(Re, e_D)

    # Use Gnielinski with explicit friction factor
    Nu_with_f = ca.nusselt_gnielinski(Re, Pr, f)

    # Compare with auto friction (smooth pipe)
    Nu_auto = ca.nusselt_gnielinski(Re, Pr)

    print(f"\nReynolds number: Re = {Re:.0e}")
    print(f"Prandtl number:  Pr = {Pr:.2f}")
    print(f"Pipe diameter:   D = {D * 1000:.0f} mm")
    print(f"Roughness:       ε = {eps * 1e6:.1f} μm (commercial steel)")
    print(f"Relative rough:  ε/D = {e_D:.6f}")

    print("\nResults:")
    print(f"  Friction factor:          f = {f:.6f}")
    print(f"  Nu (with friction):       Nu = {Nu_with_f:.2f}")
    print(f"  Nu (auto, smooth pipe):   Nu = {Nu_auto:.2f}")
    print(f"  Difference:               {(Nu_with_f - Nu_auto) / Nu_auto * 100:+.2f}%")

    print("\nNote: Roughness affects both friction and heat transfer!")


if __name__ == "__main__":
    main()
