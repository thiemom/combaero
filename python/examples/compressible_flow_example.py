#!/usr/bin/env python3
"""Compressible flow example demonstrating nozzle and pipe flow calculations.

This example shows:
1. Isentropic nozzle flow (subsonic and choked)
2. Critical pressure ratio calculation
3. Converging-diverging nozzle with axial profile
4. Fanno flow (adiabatic pipe with friction)
5. Friction factor correlations
"""

import combaero as cb


def main() -> None:
    print("=" * 60)
    print("Compressible Flow Examples")
    print("=" * 60)

    # Use standard dry air composition
    X = cb.standard_dry_air_composition()

    # -------------------------------------------------------------------------
    # 1. Isentropic Nozzle Flow
    # -------------------------------------------------------------------------
    print("\n1. Isentropic Nozzle Flow")
    print("-" * 40)

    T0 = 500.0  # K (stagnation temperature)
    P0 = 500000.0  # Pa (stagnation pressure, 5 bar)
    A_eff = 0.001  # m^2 (effective flow area)

    # Calculate critical pressure ratio
    P_ratio_crit = cb.critical_pressure_ratio(T0, P0, X)
    P_crit = P0 * P_ratio_crit
    print(f"Stagnation: T0 = {T0:.1f} K, P0 = {P0 / 1000:.1f} kPa")
    print(f"Critical pressure ratio: P*/P0 = {P_ratio_crit:.4f}")
    print(f"Critical pressure: P* = {P_crit / 1000:.1f} kPa")

    # Subsonic flow (P_back > P_crit)
    P_back_sub = 400000.0  # Pa
    sol_sub = cb.nozzle_flow(T0, P0, P_back_sub, A_eff, X)
    print(f"\nSubsonic case (P_back = {P_back_sub / 1000:.0f} kPa):")
    print(f"  Mach number: {sol_sub.M:.3f}")
    print(f"  Mass flow: {sol_sub.mdot:.4f} kg/s")
    print(f"  Outlet T: {sol_sub.outlet.T:.1f} K")
    print(f"  Choked: {sol_sub.choked}")

    # Choked flow (P_back < P_crit)
    P_back_choked = 100000.0  # Pa (atmospheric)
    sol_choked = cb.nozzle_flow(T0, P0, P_back_choked, A_eff, X)
    print(f"\nChoked case (P_back = {P_back_choked / 1000:.0f} kPa):")
    print(f"  Mach number at throat: {sol_choked.M:.3f}")
    print(f"  Mass flow: {sol_choked.mdot:.4f} kg/s")
    print(f"  Outlet T: {sol_choked.outlet.T:.1f} K")
    print(f"  Choked: {sol_choked.choked}")

    # -------------------------------------------------------------------------
    # 2. Inverse Problem: Find Area for Target Mass Flow
    # -------------------------------------------------------------------------
    print("\n2. Inverse Problem: Find Area for Target Mass Flow")
    print("-" * 40)

    mdot_target = 0.5  # kg/s
    P_back = 300000.0  # Pa

    A_required = cb.solve_A_eff_from_mdot(T0, P0, P_back, mdot_target, X)
    print(f"Target mass flow: {mdot_target:.2f} kg/s")
    print(f"Required effective area: {A_required * 1e4:.2f} cm²")

    # Verify
    sol_verify = cb.nozzle_flow(T0, P0, P_back, A_required, X)
    print(f"Verification: mdot = {sol_verify.mdot:.4f} kg/s")

    # -------------------------------------------------------------------------
    # 3. Converging-Diverging Nozzle
    # -------------------------------------------------------------------------
    print("\n3. Converging-Diverging Nozzle")
    print("-" * 40)

    # Nozzle geometry
    A_inlet = 0.01  # m^2 (inlet area)
    A_throat = 0.005  # m^2 (throat area)
    A_exit = 0.008  # m^2 (exit area)
    x_throat = 0.1  # m (throat position)
    x_exit = 0.2  # m (exit position)

    T0_cd = 600.0  # K
    P0_cd = 1000000.0  # Pa (10 bar)
    P_exit = 200000.0  # Pa

    sol_cd = cb.nozzle_cd(
        T0_cd, P0_cd, P_exit, A_inlet, A_throat, A_exit, x_throat, x_exit, X, n_stations=50
    )

    print(f"Inlet: A = {A_inlet * 1e4:.1f} cm², Throat: A = {A_throat * 1e4:.1f} cm²")
    print(f"Exit: A = {A_exit * 1e4:.1f} cm²")
    print(f"Stagnation: T0 = {sol_cd.T0:.1f} K, P0 = {sol_cd.P0 / 1000:.0f} kPa")
    print(f"Mass flow: {sol_cd.mdot:.4f} kg/s")
    print(f"Choked: {sol_cd.choked}")
    print(f"Throat location: x = {sol_cd.x_throat:.3f} m")

    # Print a few stations from the profile
    print("\nAxial profile (selected stations):")
    print(f"{'x [m]':>8} {'A [cm²]':>10} {'M [-]':>8} {'T [K]':>8} {'P [kPa]':>10}")
    for i in range(0, len(sol_cd.profile), len(sol_cd.profile) // 5):
        st = sol_cd.profile[i]
        print(f"{st.x:8.3f} {st.A * 1e4:10.2f} {st.M:8.3f} {st.T:8.1f} {st.P / 1000:10.1f}")

    # -------------------------------------------------------------------------
    # 4. Fanno Flow (Adiabatic Pipe with Friction)
    # -------------------------------------------------------------------------
    print("\n4. Fanno Flow (Adiabatic Pipe with Friction)")
    print("-" * 40)

    T_in = 400.0  # K
    P_in = 500000.0  # Pa
    u_in = 100.0  # m/s
    D = 0.05  # m (pipe diameter)
    L = 5.0  # m (pipe length)

    # Calculate friction factor using Colebrook
    rho = cb.density(T_in, P_in, X)
    mu = cb.viscosity(T_in, P_in, X)
    Re = rho * u_in * D / mu
    e_D = 0.0001  # relative roughness (smooth pipe)
    f = cb.friction_colebrook(Re, e_D)

    print(f"Inlet: T = {T_in:.0f} K, P = {P_in / 1000:.0f} kPa, u = {u_in:.0f} m/s")
    print(f"Pipe: D = {D * 100:.1f} cm, L = {L:.1f} m")
    print(f"Reynolds number: Re = {Re:.0f}")
    print(f"Friction factor (Colebrook): f = {f:.5f}")

    # Calculate maximum length before choking
    L_star = cb.fanno_max_length(T_in, P_in, u_in, D, f, X)
    print(f"Maximum length before choking: L* = {L_star:.2f} m")

    # Solve Fanno flow
    sol_fanno = cb.fanno_pipe(T_in, P_in, u_in, L, D, f, X, store_profile=True)

    print(f"\nOutlet conditions (L = {L:.1f} m):")
    print(f"  T = {sol_fanno.outlet.T:.1f} K")
    print(f"  P = {sol_fanno.outlet.P / 1000:.1f} kPa")
    print(f"  Choked: {sol_fanno.choked}")
    print(f"  Mass flow: {sol_fanno.mdot:.4f} kg/s")

    # -------------------------------------------------------------------------
    # 5. Friction Factor Comparison
    # -------------------------------------------------------------------------
    print("\n5. Friction Factor Correlations")
    print("-" * 40)

    Re_test = 100000.0
    e_D_test = 0.001

    f_haaland = cb.friction_haaland(Re_test, e_D_test)
    f_serghides = cb.friction_serghides(Re_test, e_D_test)
    f_colebrook = cb.friction_colebrook(Re_test, e_D_test)

    print(f"Re = {Re_test:.0f}, ε/D = {e_D_test}")
    print(f"  Haaland:   f = {f_haaland:.6f}")
    print(f"  Serghides: f = {f_serghides:.6f}")
    print(f"  Colebrook: f = {f_colebrook:.6f} (reference)")
    print(f"  Haaland error:   {100 * (f_haaland - f_colebrook) / f_colebrook:+.2f}%")
    print(f"  Serghides error: {100 * (f_serghides - f_colebrook) / f_colebrook:+.3f}%")

    print("\n" + "=" * 60)
    print("Done!")


if __name__ == "__main__":
    main()
