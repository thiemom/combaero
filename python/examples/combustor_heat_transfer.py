#!/usr/bin/env python3
"""Integrated combustor design with heat transfer.

This example demonstrates:
- Complete combustor design workflow
- Combustion + heat transfer integration
- Pressure drop and heat transfer combined
- Realistic gas turbine combustor analysis
- Using multiple CombAero functions together
"""

from __future__ import annotations

import numpy as np

import combaero as ca
from combaero.species import SpeciesLocator


def main() -> None:
    print("=" * 80)
    print("INTEGRATED COMBUSTOR DESIGN WITH HEAT TRANSFER")
    print("=" * 80)

    sp = SpeciesLocator.from_core()

    # =========================================================================
    # Example 1: Combustor design workflow
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 1: Gas Turbine Combustor Design")
    print("=" * 80)

    # Operating conditions
    P_combustor = 15e5  # Combustor pressure [Pa] = 15 bar
    T_air_in = 700.0  # Compressor exit temperature [K]
    mdot_air = 10.0  # Air mass flow [kg/s]
    mdot_fuel = 0.2  # Fuel mass flow [kg/s]
    phi = 0.8  # Equivalence ratio (lean burn)

    # Fuel composition (natural gas)
    X_fuel = sp.empty()
    X_fuel[sp.indices["CH4"]] = 0.95
    X_fuel[sp.indices["C2H6"]] = 0.03
    X_fuel[sp.indices["N2"]] = 0.02

    # Air composition
    X_air = ca.standard_dry_air_composition()

    print("\nOperating conditions:")
    print(f"  Pressure:          P = {P_combustor/1e5:.0f} bar")
    print(f"  Air inlet temp:    T = {T_air_in:.0f} K")
    print(f"  Air mass flow:     ṁ_air = {mdot_air:.1f} kg/s")
    print(f"  Fuel mass flow:    ṁ_fuel = {mdot_fuel:.2f} kg/s")
    print(f"  Equivalence ratio: φ = {phi:.2f}")

    # Create streams
    air = ca.Stream()
    air.T = T_air_in
    air.P = P_combustor
    air.X = X_air
    air.mdot = mdot_air

    fuel = ca.Stream()
    fuel.T = 300.0  # Fuel at ambient
    fuel.P = P_combustor
    fuel.X = X_fuel
    fuel.mdot = mdot_fuel

    # Mix fuel and air
    mixed = ca.mix([fuel, air], P_out=P_combustor)

    print("\nMixed stream before combustion:")
    print(f"  Temperature:       T = {mixed.T:.1f} K")
    print(f"  Total mass flow:   ṁ = {mixed.mdot:.2f} kg/s")

    # Combustion with equilibrium
    burned = ca.combustion_equilibrium(mixed.T, mixed.X, mixed.P)

    print("\nAfter combustion:")
    print(f"  Temperature:       T = {burned.T:.1f} K (adiabatic flame temp)")
    print(f"  Enthalpy:          h = {burned.h:.1f} J/mol")

    # =========================================================================
    # Example 2: Combustor geometry and residence time
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 2: Combustor Geometry and Residence Time")
    print("=" * 80)

    # Annular combustor geometry
    D_outer = 0.6  # Outer diameter [m]
    D_inner = 0.4  # Inner diameter [m]
    L_combustor = 0.4  # Length [m]

    # Calculate volume using annular_area helper
    A_annular = ca.annular_area(D_outer, D_inner)
    V_combustor = A_annular * L_combustor

    # Hydraulic diameter
    Dh = ca.hydraulic_diameter_annulus(D_outer, D_inner)

    # Average properties (use average of inlet and outlet)
    T_avg = (mixed.T + burned.T) / 2
    rho_avg = ca.density(T_avg, P_combustor, burned.X)

    # Residence time
    tau = ca.residence_time_mdot(V_combustor, mixed.mdot, rho_avg)
    SV = ca.space_velocity(mixed.mdot / rho_avg, V_combustor)

    print("\nCombustor geometry:")
    print(f"  Outer diameter:    D_o = {D_outer*1000:.0f} mm")
    print(f"  Inner diameter:    D_i = {D_inner*1000:.0f} mm")
    print(f"  Hydraulic diameter: Dh = {Dh*1000:.0f} mm")
    print(f"  Length:            L = {L_combustor*1000:.0f} mm")
    print(f"  Annular area:      A = {A_annular*1e4:.2f} cm2")
    print(f"  Volume:            V = {V_combustor*1000:.2f} L")

    print("\nFlow characteristics:")
    print(f"  Average temp:      T_avg = {T_avg:.0f} K")
    print(f"  Average density:   ρ_avg = {rho_avg:.3f} kg/m3")
    print(f"  Residence time:    τ = {tau*1000:.1f} ms")
    print(f"  Space velocity:    SV = {SV*3600:.0f} h⁻¹")

    # =========================================================================
    # Example 3: Pressure drop through combustor
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 3: Pressure Drop Through Combustor")
    print("=" * 80)

    # Flow velocity
    v_avg = mixed.mdot / (rho_avg * A_annular)

    # Use pressure_drop_pipe composite function
    # Note: Using cast_iron roughness as proxy for combustor liner
    eps = ca.pipe_roughness("cast_iron")
    dP_friction, Re, f = ca.pressure_drop_pipe(
        T_avg, P_combustor, burned.X, v_avg, Dh, L_combustor, eps, "haaland"
    )

    # Also get viscosity for display
    mu_avg = ca.viscosity(T_avg, P_combustor, burned.X)

    # Total pressure drop (friction + other losses, typically 3-5%)
    dP_total = dP_friction * 1.5  # Account for additional losses

    print("\nFlow conditions:")
    print(f"  Average velocity:  v = {v_avg:.1f} m/s")
    print(f"  Reynolds number:   Re = {Re:.2e}")
    print(f"  Viscosity:         μ = {mu_avg*1e6:.2f} μPa*s")

    print("\nPressure drop:")
    print(f"  Friction factor:   f = {f:.6f}")
    print(f"  Friction loss:     DeltaP_f = {dP_friction/1000:.2f} kPa")
    print(f"  Total loss:        DeltaP = {dP_total/1000:.2f} kPa")
    print(f"  Pressure ratio:    DeltaP/P = {dP_total/P_combustor*100:.2f}%")

    # =========================================================================
    # Example 4: Heat transfer to combustor liner
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 4: Heat Transfer to Combustor Liner")
    print("=" * 80)

    # Hot gas side (inside combustor)
    T_gas = burned.T  # Hot combustion products
    k_gas = ca.thermal_conductivity(T_gas, P_combustor, burned.X)
    Pr_gas = ca.prandtl(T_gas, P_combustor, burned.X)

    # Calculate Nusselt number for hot gas
    Nu_gas = ca.nusselt_dittus_boelter(Re, Pr_gas, heating=False)  # Cooling the gas
    h_gas = ca.htc_from_nusselt(Nu_gas, k_gas, Dh)

    # Cooling air side (outside liner, simplified)
    T_cool = 600.0  # Cooling air temperature [K]
    v_cool = 50.0  # High velocity cooling air [m/s]
    D_cool = 0.01  # Cooling passage hydraulic diameter [m]

    # Cooling air properties
    rho_cool = ca.density(T_cool, P_combustor, X_air)
    mu_cool = ca.viscosity(T_cool, P_combustor, X_air)
    k_cool = ca.thermal_conductivity(T_cool, P_combustor, X_air)
    Pr_cool = ca.prandtl(T_cool, P_combustor, X_air)

    Re_cool = rho_cool * v_cool * D_cool / mu_cool
    Nu_cool = ca.nusselt_dittus_boelter(Re_cool, Pr_cool, heating=True)  # Heating the cooling air
    h_cool = ca.htc_from_nusselt(Nu_cool, k_cool, D_cool)

    # Liner wall (assume thermal barrier coating + metal)
    t_tbc = 0.001  # TBC thickness [m]
    k_tbc = 1.0  # TBC thermal conductivity [W/(m*K)]
    t_metal = 0.003  # Metal thickness [m]
    k_metal = 20.0  # Metal thermal conductivity [W/(m*K)]

    # Overall heat transfer coefficient
    R_tbc = t_tbc / k_tbc
    R_metal = t_metal / k_metal
    U = 1.0 / (1.0 / h_gas + R_tbc + R_metal + 1.0 / h_cool)

    # Heat flux
    q = U * (T_gas - T_cool)

    # Total heat transfer (approximate, using mean diameter)
    D_mean = (D_outer + D_inner) / 2
    A_surface = np.pi * D_mean * L_combustor
    Q_total = q * A_surface

    print("\nHot gas side (combustion products):")
    print(f"  Temperature:       T = {T_gas:.0f} K")
    print(f"  Reynolds number:   Re = {Re:.2e}")
    print(f"  Prandtl number:    Pr = {Pr_gas:.3f}")
    print(f"  Nusselt number:    Nu = {Nu_gas:.1f}")
    print(f"  HTC:               h = {h_gas:.1f} W/(m2*K)")

    print("\nCooling air side:")
    print(f"  Temperature:       T = {T_cool:.0f} K")
    print(f"  Velocity:          v = {v_cool:.0f} m/s")
    print(f"  Reynolds number:   Re = {Re_cool:.2e}")
    print(f"  Nusselt number:    Nu = {Nu_cool:.1f}")
    print(f"  HTC:               h = {h_cool:.1f} W/(m2*K)")

    print("\nLiner wall:")
    print(f"  TBC thickness:     t = {t_tbc*1000:.1f} mm (k = {k_tbc:.1f} W/(m*K))")
    print(f"  Metal thickness:   t = {t_metal*1000:.1f} mm (k = {k_metal:.0f} W/(m*K))")

    print("\nHeat transfer:")
    print(f"  Overall HTC:       U = {U:.1f} W/(m2*K)")
    print(f"  Heat flux:         q = {q/1000:.1f} kW/m2")
    print(f"  Total heat loss:   Q = {Q_total/1e6:.2f} MW")
    print("  (Fraction of combustion power)")

    # =========================================================================
    # Example 5: Combustor performance summary
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 5: Combustor Performance Summary")
    print("=" * 80)

    # Combustion efficiency (simplified)
    LHV_CH4 = 50e6  # Lower heating value of CH4 [J/kg]
    Q_combustion = mdot_fuel * LHV_CH4  # Approximate

    # Temperature rise
    dT_combustion = burned.T - mixed.T

    # Pressure loss percentage
    pressure_loss_pct = dP_total / P_combustor * 100

    # Pattern factor (temperature non-uniformity, typical 0.2-0.3)
    pattern_factor = 0.25  # Assumed

    print("\nCombustor Performance:")
    print(f"  Inlet temperature:     T_in = {mixed.T:.0f} K")
    print(f"  Outlet temperature:    T_out = {burned.T:.0f} K")
    print(f"  Temperature rise:      DeltaT = {dT_combustion:.0f} K")
    print(f"  Pressure loss:         DeltaP/P = {pressure_loss_pct:.2f}%")
    print(f"  Residence time:        τ = {tau*1000:.1f} ms")
    print(f"  Heat release:          Q = {Q_combustion/1e6:.1f} MW")
    print(
        f"  Heat loss to walls:    Q_loss = {Q_total/1e6:.2f} MW ({Q_total/Q_combustion*100:.1f}%)"
    )
    print(f"  Pattern factor:        PF = {pattern_factor:.2f}")

    print("\nDesign metrics:")
    print(f"  Loading parameter:     ṁ√T/P = {mdot_air*np.sqrt(T_air_in)/P_combustor:.4f}")
    print(f"  Volumetric heat release: {Q_combustion/V_combustor/1e9:.1f} GW/m3")

    # =========================================================================
    # Example 6: Parametric study - effect of equivalence ratio
    # =========================================================================
    print("\n" + "=" * 80)
    print("Example 6: Effect of Equivalence Ratio on Combustor Performance")
    print("=" * 80)

    print(f"\n{'φ':>6s}  {'T_out [K]':>12s}  {'τ [ms]':>10s}  {'DeltaP/P [%]':>12s}")
    print("-" * 45)

    for phi_test in [0.6, 0.7, 0.8, 0.9, 1.0]:
        # Adjust fuel flow for equivalence ratio
        # (simplified - would need proper stoichiometry calculation)
        mdot_fuel_test = mdot_fuel * phi_test / phi

        # Create streams
        fuel_test = ca.Stream()
        fuel_test.T = 300.0
        fuel_test.P = P_combustor
        fuel_test.X = X_fuel
        fuel_test.mdot = mdot_fuel_test

        # Mix and burn
        mixed_test = ca.mix([fuel_test, air], P_out=P_combustor)
        burned_test = ca.combustion_equilibrium(mixed_test.T, mixed_test.X, mixed_test.P)

        # Recalculate residence time and pressure drop
        T_avg_test = (mixed_test.T + burned_test.T) / 2
        rho_avg_test = ca.density(T_avg_test, P_combustor, burned_test.X)
        tau_test = ca.residence_time_mdot(V_combustor, mixed_test.mdot, rho_avg_test)

        v_avg_test = mixed_test.mdot / (rho_avg_test * A_annular)
        # Use pressure_drop_pipe composite function
        dP_friction_test, Re_test, f_test = ca.pressure_drop_pipe(
            T_avg_test, P_combustor, burned_test.X, v_avg_test, Dh, L_combustor, eps, "haaland"
        )
        dP_test = dP_friction_test * 1.5  # Account for additional losses

        print(
            f"{phi_test:6.2f}  {burned_test.T:12.1f}  {tau_test*1000:10.2f}  {dP_test/P_combustor*100:12.3f}"
        )

    print("\nNote: Leaner mixtures → lower temperature, longer residence time")


if __name__ == "__main__":
    main()
