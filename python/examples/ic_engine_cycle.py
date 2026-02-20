#!/usr/bin/env python3
"""Otto cycle (spark-ignition IC engine) simulation.

Models an ideal Otto cycle using real gas properties:

  State 1 -> 2 : Isentropic compression  (s = const, V decreases)
  State 2 -> 3 : Constant-volume heat addition (combustion at TDC)
  State 3 -> 4 : Isentropic expansion    (s = const, V increases)
  State 4 -> 1 : Constant-volume heat rejection (exhaust / BDC)

Key combaero tools used:
  - set_fuel_stream_for_phi + mix + complete_combustion  (heat addition)
  - s / s_mass / calc_T_from_s  (isentropic path with variable gamma)
  - u_mass  (constant-volume heat rejection: q_out = u4 - u1 at fixed composition)
  - cv_mass  (constant-volume heat addition: q_in via integral, composition changes)
  - density / specific_gas_constant  (thermodynamic state)
  - isentropic_expansion_coefficient  (local gamma along path)
"""

from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np

import combaero as ca
from combaero.species import SpeciesLocator

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def adiabatic_path(
    T_start: float, P_start: float, P_end: float, X: list, n: int = 200
) -> tuple[np.ndarray, np.ndarray]:
    """Pressure-volume path for isentropic process (variable gamma)."""
    # Use entropy conservation: at each P along the path, find T
    s_ref = ca.s(T_start, X, P_start)  # molar entropy [J/(mol K)]
    P_path = np.linspace(P_start, P_end, n)
    V_path = np.empty(n)
    for i, P in enumerate(P_path):
        T = ca.calc_T_from_s(s_ref, P, X)
        rho = ca.density(T, P, X)
        V_path[i] = 1.0 / rho  # specific volume [m3/kg]
    return V_path, P_path


def main() -> None:
    sp = SpeciesLocator.from_core()

    # =========================================================================
    # Engine parameters
    # =========================================================================
    CR = 10.0  # compression ratio V1/V2
    phi = 0.9  # equivalence ratio (slightly lean)
    T1 = 350.0  # K  intake temperature (slightly above ambient)
    P1 = 101325.0  # Pa intake pressure (naturally aspirated)

    X_air = ca.standard_dry_air_composition()
    X_ch4 = sp.empty()
    X_ch4[sp.indices["CH4"]] = 1.0

    print("=" * 65)
    print("Otto Cycle -- Spark-Ignition IC Engine")
    print("=" * 65)
    print(f"\n  Compression ratio : CR  = {CR:.1f}")
    print(f"  Equivalence ratio : phi = {phi:.2f}")
    print(f"  Intake conditions : T1  = {T1:.0f} K,  P1 = {P1 / 1e3:.1f} kPa")

    # =========================================================================
    # State 1: BDC after intake (fresh charge = premixed fuel + air)
    # =========================================================================
    air = ca.Stream()
    air.T, air.P, air.X, air.mdot = T1, P1, X_air, 1.0

    fuel = ca.Stream()
    fuel.T, fuel.P, fuel.X = T1, P1, X_ch4

    fuel = ca.set_fuel_stream_for_phi(phi, fuel, air)
    charge = ca.mix([fuel, air])

    X1 = charge.X  # unburned mixture composition
    rho1 = ca.density(T1, P1, X1)
    v1 = 1.0 / rho1  # specific volume [m3/kg]

    print("\n--- State 1 (BDC, unburned charge) ---")
    print(f"  T1  = {T1:.1f} K")
    print(f"  P1  = {P1 / 1e3:.2f} kPa")
    print(f"  v1  = {v1:.4f} m3/kg")
    print(f"  FAR = {fuel.mdot / air.mdot:.5f}  (fuel/air mass ratio)")

    # =========================================================================
    # State 2: TDC after isentropic compression
    # =========================================================================
    # P2/P1 = CR^gamma  (first estimate with local gamma, then exact via s=const)
    gamma1 = ca.isentropic_expansion_coefficient(T1, X1)
    P2_est = P1 * CR**gamma1
    s1 = ca.s(T1, X1, P1)  # molar entropy [J/(mol K)]
    T2, P2 = ca.calc_T_from_s(s1, P2_est, X1), P2_est
    # Refine: actual v2 = v1/CR, so iterate P2 until v2 = v1/CR at s=const
    for _ in range(10):
        rho2 = ca.density(T2, P2, X1)
        v2_actual = 1.0 / rho2
        v2_target = v1 / CR
        P2 = P2 * (v2_actual / v2_target) ** ca.isentropic_expansion_coefficient(T2, X1)
        T2 = ca.calc_T_from_s(s1, P2, X1)

    rho2 = ca.density(T2, P2, X1)
    v2 = 1.0 / rho2

    print("\n--- State 2 (TDC, end of compression) ---")
    print(f"  T2  = {T2:.1f} K")
    print(f"  P2  = {P2 / 1e3:.2f} kPa")
    print(f"  v2  = {v2:.4f} m3/kg  (v1/v2 = {v1 / v2:.2f})")

    # =========================================================================
    # State 3: TDC after constant-volume combustion
    # =========================================================================
    # Burn the compressed charge; volume (and thus density) stays constant.
    # complete_combustion gives T_ad at the post-compression pressure.
    burned = ca.complete_combustion(T2, X1, P2)
    T3 = burned.T
    X3 = burned.X

    # At constant volume: rho3 = rho2, so P3 = rho2 * R_specific(X3) * T3
    # Use density identity: rho = P / (R_sp * T)  =>  P3 = P2 * (T3/T2) * (R3/R2)
    R2 = ca.specific_gas_constant(X1)
    R3 = ca.specific_gas_constant(X3)
    P3 = P2 * (T3 / T2) * (R3 / R2)
    v3 = v2  # constant volume

    # Heat added per kg of charge: constant-volume process -> q = integral cv dT
    T_cv_in = np.linspace(T2, T3, 200)
    q_in = np.trapz([ca.cv_mass(T, X1) for T in T_cv_in], T_cv_in)  # J/kg

    print("\n--- State 3 (TDC, after combustion) ---")
    print(f"  T3  = {T3:.1f} K")
    print(f"  P3  = {P3 / 1e3:.2f} kPa")
    print(f"  v3  = {v3:.4f} m3/kg")
    print(f"  q_in  ~ {q_in / 1e3:.1f} kJ/kg")

    # =========================================================================
    # State 4: BDC after isentropic expansion
    # =========================================================================
    # Expand back to v4 = v1 (same specific volume as State 1)
    # Mirror of compression: find P4 such that v4 = v1 at constant s
    gamma3 = ca.isentropic_expansion_coefficient(T3, X3)
    P4_est = P3 / CR**gamma3
    s3 = ca.s(T3, X3, P3)  # molar entropy [J/(mol K)]
    P4 = P4_est
    T4 = ca.calc_T_from_s(s3, P4, X3)
    for _ in range(10):
        rho4 = ca.density(T4, P4, X3)
        v4_actual = 1.0 / rho4
        P4 = P4 * (v4_actual / v1) ** ca.isentropic_expansion_coefficient(T4, X3)
        T4 = ca.calc_T_from_s(s3, P4, X3)

    rho4 = ca.density(T4, P4, X3)
    v4 = 1.0 / rho4

    print("\n--- State 4 (BDC, end of expansion) ---")
    print(f"  T4  = {T4:.1f} K")
    print(f"  P4  = {P4 / 1e3:.2f} kPa")
    print(f"  v4  = {v4:.4f} m3/kg  (v4/v1 = {v4 / v1:.3f})")

    # =========================================================================
    # Cycle performance
    # =========================================================================
    # Heat rejected: constant-volume process 4->1, composition fixed at X3
    # u_mass is valid here since both states share the same composition
    q_out = ca.u_mass(T4, X3) - ca.u_mass(T1, X3)  # J/kg (positive)

    w_net = q_in - q_out  # net work [J/kg]
    eta_th = w_net / q_in  # thermal efficiency

    # Ideal Otto efficiency for comparison (constant gamma)
    gamma_avg = 0.5 * (gamma1 + gamma3)
    eta_otto_ideal = 1.0 - CR ** (1.0 - gamma_avg)

    print(f"\n{'=' * 65}")
    print("Cycle Performance")
    print(f"{'=' * 65}")
    print(f"  q_in            : {q_in / 1e3:.2f} kJ/kg")
    print(f"  q_out           : {q_out / 1e3:.2f} kJ/kg")
    print(f"  w_net           : {w_net / 1e3:.2f} kJ/kg")
    print(f"  eta_thermal     : {eta_th * 100:.2f} %")
    print(f"  eta_Otto(ideal) : {eta_otto_ideal * 100:.2f} %  (avg gamma={gamma_avg:.3f})")

    # =========================================================================
    # P-V and T-s diagrams
    # =========================================================================
    # --- isentropic paths (P-V) ---
    v_comp, P_comp = adiabatic_path(T1, P1, P2, X1)
    v_exp, P_exp = adiabatic_path(T3, P3, P4, X3)

    # --- entropy at each state [J/(kg K)] ---
    s1 = ca.s_mass(T1, X1, P1)
    s2 = ca.s_mass(T2, X1, P2)  # == s1 (isentropic)
    s3 = ca.s_mass(T3, X3, P3)
    s4 = ca.s_mass(T4, X3, P4)  # == s3 (isentropic)

    # isochoric heat addition 2->3: ds = cv/T dT  (variable composition approx: use X1)
    n_ts = 200
    T_23 = np.linspace(T2, T3, n_ts)
    s_23 = s2 + np.array(
        [
            np.trapz([ca.cv_mass(t, X1) / t for t in T_23[: i + 1]], T_23[: i + 1])
            for i in range(n_ts)
        ]
    )

    # isochoric heat rejection 4->1: ds = cv/T dT (burned gas X3)
    T_41 = np.linspace(T4, T1, n_ts)
    s_41 = s4 + np.array(
        [
            np.trapz([ca.cv_mass(t, X3) / t for t in T_41[: i + 1]], T_41[: i + 1])
            for i in range(n_ts)
        ]
    )

    fig, (ax_pv, ax_ts) = plt.subplots(1, 2, figsize=(13, 5))
    title = f"Otto Cycle  (CR={CR:.0f}, phi={phi:.2f}, eta={eta_th * 100:.1f} %)"
    fig.suptitle(title, fontsize=13)

    # --- P-V ---
    ax_pv.plot(v_comp * 1e3, P_comp / 1e3, label="1-2 Compression")
    ax_pv.plot([v2 * 1e3, v3 * 1e3], [P2 / 1e3, P3 / 1e3], label="2-3 Combustion")
    ax_pv.plot(v_exp * 1e3, P_exp / 1e3, label="3-4 Expansion")
    ax_pv.plot([v4 * 1e3, v1 * 1e3], [P4 / 1e3, P1 / 1e3], label="4-1 Exhaust")
    for lbl, v, Ppt in [("1", v1, P1), ("2", v2, P2), ("3", v3, P3), ("4", v4, P4)]:
        ax_pv.plot(v * 1e3, Ppt / 1e3, "ko", markersize=5)
        ax_pv.annotate(
            lbl,
            xy=(v * 1e3, Ppt / 1e3),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=11,
            fontweight="bold",
        )
    ax_pv.set_xlabel("Specific volume [L/kg]")
    ax_pv.set_ylabel("Pressure [kPa]")
    ax_pv.set_title("P-V Diagram")
    ax_pv.legend()
    ax_pv.grid(True, alpha=0.3)

    # --- T-s ---
    ax_ts.plot([s1, s2], [T1, T2], label="1-2 Compression")  # vertical line (s=const)
    ax_ts.plot(s_23, T_23, label="2-3 Combustion")
    ax_ts.plot([s3, s4], [T3, T4], label="3-4 Expansion")  # vertical line (s=const)
    ax_ts.plot(s_41, T_41, label="4-1 Exhaust")
    for lbl, s_pt, T_pt in [("1", s1, T1), ("2", s2, T2), ("3", s3, T3), ("4", s4, T4)]:
        ax_ts.plot(s_pt, T_pt, "ko", markersize=5)
        ax_ts.annotate(
            lbl,
            xy=(s_pt, T_pt),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=11,
            fontweight="bold",
        )
    ax_ts.set_xlabel("Specific entropy [J/(kg K)]")
    ax_ts.set_ylabel("Temperature [K]")
    ax_ts.set_title("T-s Diagram")
    ax_ts.legend()
    ax_ts.grid(True, alpha=0.3)

    plt.tight_layout()
    if plt.get_backend().lower() == "agg":
        out_path = pathlib.Path(__file__).parent / "ic_engine_cycle.png"
        plt.savefig(out_path, dpi=150)
        print(f"\nPlot saved to '{out_path}'")
    else:
        plt.show()


if __name__ == "__main__":
    main()
