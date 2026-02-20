#!/usr/bin/env python3
"""Brayton cycle (open-cycle gas turbine) simulation.

Models an ideal open Brayton cycle with real gas properties:

  State 1 -> 2 : Isentropic compression   (compressor, s = const)
  State 2 -> 3 : Isobaric heat addition   (combustor, P = const)
  State 3 -> 4 : Isentropic expansion     (turbine, s = const)
  State 4 -> 1 : Isobaric heat rejection  (exhaust to atmosphere)

Key combaero tools used:
  - set_fuel_stream_for_phi + mix + complete_combustion  (combustor)
  - s / calc_T_from_s  (isentropic compression/expansion)
  - h_mass  (open-system energy balance: w_net = h_in - h_out per kg air)
  - cp_mass  (compressor work: same composition, sensible only)
  - s_mass   (entropy at each state for T-s diagram)
  - density / isentropic_expansion_coefficient  (state properties)

Energy balance (per kg air):
  w_net = h_air(T1) + FAR*h_fuel(T_fuel) - (1+FAR)*h_exhaust(T4)
  q_in  = (1+FAR)*h_exhaust(T3) - h_air(T2) - FAR*h_fuel(T_fuel)
  w_c   = h_air(T2) - h_air(T1)   [compressor, same composition]
  w_t   = w_net + w_c              [turbine, by difference]
"""

from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np

import combaero as ca
from combaero.species import SpeciesLocator


def dh_same_X(T_a: float, T_b: float, X: list, n: int = 300) -> float:
    """Isobaric enthalpy change [J/kg] via cp_mass integral (same composition)."""
    Ts = np.linspace(T_a, T_b, n)
    return float(np.trapz([ca.cp_mass(T, X) for T in Ts], Ts))


def isentropic_T(T_start: float, P_start: float, P_end: float, X: list) -> float:
    """Temperature after isentropic process from P_start to P_end."""
    return ca.calc_T_from_s(ca.s(T_start, X, P_start), P_end, X)


def main() -> None:
    sp = SpeciesLocator.from_core()

    # =========================================================================
    # Engine parameters
    # =========================================================================
    PR = 20.0  # pressure ratio P2/P1
    phi = 0.30  # equivalence ratio (lean, typical gas turbine)
    T1 = 288.15  # K  compressor inlet (ISA sea level)
    P1 = 101325.0  # Pa

    X_air = ca.standard_dry_air_composition()
    X_ch4 = sp.empty()
    X_ch4[sp.indices["CH4"]] = 1.0

    print("=" * 65)
    print("Brayton Cycle -- Open-Cycle Gas Turbine")
    print("=" * 65)
    print(f"\n  Pressure ratio    : PR  = {PR:.0f}")
    print(f"  Equivalence ratio : phi = {phi:.2f}")
    print(f"  Inlet conditions  : T1  = {T1:.2f} K,  P1 = {P1 / 1e3:.3f} kPa")

    # =========================================================================
    # State 1: Compressor inlet
    # =========================================================================
    rho1 = ca.density(T1, P1, X_air)
    gamma1 = ca.isentropic_expansion_coefficient(T1, X_air)

    print("\n--- State 1 (Compressor inlet) ---")
    print(f"  T1    = {T1:.2f} K")
    print(f"  P1    = {P1 / 1e3:.3f} kPa")
    print(f"  rho1  = {rho1:.4f} kg/m3")
    print(f"  gamma = {gamma1:.4f}")

    # =========================================================================
    # State 2: After isentropic compression
    # =========================================================================
    P2 = P1 * PR
    T2 = isentropic_T(T1, P1, P2, X_air)
    rho2 = ca.density(T2, P2, X_air)

    print("\n--- State 2 (Compressor exit) ---")
    print(f"  T2    = {T2:.1f} K")
    print(f"  P2    = {P2 / 1e3:.1f} kPa")
    print(f"  rho2  = {rho2:.4f} kg/m3")

    # =========================================================================
    # State 3: After isobaric combustion
    # =========================================================================
    air = ca.Stream()
    air.T, air.P, air.X, air.mdot = T2, P2, X_air, 1.0

    fuel = ca.Stream()
    fuel.T, fuel.P, fuel.X = 300.0, P2, X_ch4

    fuel = ca.set_fuel_stream_for_phi(phi, fuel, air)
    mixed = ca.mix([fuel, air])
    burned = ca.complete_combustion(mixed.T, mixed.X, mixed.P)

    T3 = burned.T
    X3 = burned.X
    P3 = P2  # isobaric combustor
    rho3 = ca.density(T3, P3, X3)
    gamma3 = ca.isentropic_expansion_coefficient(T3, X3)

    print("\n--- State 3 (Turbine inlet) ---")
    print(f"  T3    = {T3:.1f} K  (TIT)")
    print(f"  P3    = {P3 / 1e3:.1f} kPa")
    print(f"  rho3  = {rho3:.4f} kg/m3")
    print(f"  gamma = {gamma3:.4f}")
    print(f"  FAR   = {fuel.mdot / air.mdot:.5f}")

    # =========================================================================
    # State 4: After isentropic expansion
    # =========================================================================
    P4 = P1
    T4 = isentropic_T(T3, P3, P4, X3)
    rho4 = ca.density(T4, P4, X3)

    print("\n--- State 4 (Turbine exit / exhaust) ---")
    print(f"  T4    = {T4:.1f} K")
    print(f"  P4    = {P4 / 1e3:.3f} kPa")
    print(f"  rho4  = {rho4:.4f} kg/m3")

    # =========================================================================
    # Cycle performance
    # =========================================================================
    far = fuel.mdot / air.mdot  # fuel/air mass ratio

    # Open-system energy balance (per kg air): w_net = h_in - h_out
    # h_mass uses NASA9 absolute enthalpies (incl. formation), valid across compositions
    w_net = ca.h_mass(T1, X_air) + far * ca.h_mass(fuel.T, X_ch4) - (1.0 + far) * ca.h_mass(T4, X3)
    # q_in = fuel chemical energy (LHV basis); complete_combustion is adiabatic so
    # the combustor enthalpy difference is zero by construction
    q_in = far * ca.fuel_lhv_mass(X_ch4)

    # Compressor and turbine work by difference (same composition each)
    w_c = dh_same_X(T1, T2, X_air)  # compressor: air only, sensible
    w_t = w_net + w_c  # turbine: by energy balance

    eta_th = w_net / q_in
    bwr = w_c / w_t  # back-work ratio

    # Ideal Brayton efficiency for comparison (avg gamma)
    gamma_avg = 0.5 * (gamma1 + gamma3)
    eta_brayton_ideal = 1.0 - PR ** ((1.0 - gamma_avg) / gamma_avg)

    print(f"\n{'=' * 65}")
    print("Cycle Performance")
    print(f"{'=' * 65}")
    print(f"  w_compressor    : {w_c / 1e3:.2f} kJ/kg")
    print(f"  w_turbine       : {w_t / 1e3:.2f} kJ/kg")
    print(f"  w_net           : {w_net / 1e3:.2f} kJ/kg")
    print(f"  q_in            : {q_in / 1e3:.2f} kJ/kg")
    print(f"  eta_thermal     : {eta_th * 100:.2f} %")
    print(f"  eta_Brayton(id) : {eta_brayton_ideal * 100:.2f} %  (avg gamma={gamma_avg:.3f})")
    print(f"  back-work ratio : {bwr * 100:.1f} %")

    # =========================================================================
    # Pressure ratio sweep
    # =========================================================================
    print(f"\n{'=' * 65}")
    print("Pressure Ratio Sweep (phi=0.30, T1=288 K)")
    print(f"{'=' * 65}")
    print(
        f"  {'PR':>5}  {'T2 [K]':>8}  {'T3 [K]':>8}  {'T4 [K]':>8}  {'eta [%]':>9}  {'BWR [%]':>8}"
    )
    print(f"  {'-' * 5}  {'-' * 8}  {'-' * 8}  {'-' * 8}  {'-' * 9}  {'-' * 8}")

    pr_vals, eta_vals, bwr_vals = [], [], []
    for pr in [5, 8, 10, 15, 20, 25, 30, 40]:
        p2_ = P1 * pr
        t2_ = isentropic_T(T1, P1, p2_, X_air)
        air_ = ca.Stream()
        air_.T, air_.P, air_.X, air_.mdot = t2_, p2_, X_air, 1.0
        fuel_ = ca.Stream()
        fuel_.T, fuel_.P, fuel_.X = 300.0, p2_, X_ch4
        fuel_ = ca.set_fuel_stream_for_phi(phi, fuel_, air_)
        mixed_ = ca.mix([fuel_, air_])
        burned_ = ca.complete_combustion(mixed_.T, mixed_.X, mixed_.P)
        t3_ = burned_.T
        x3_ = burned_.X
        t4_ = isentropic_T(t3_, p2_, P1, x3_)
        far_ = fuel_.mdot / air_.mdot
        wnet_ = (
            ca.h_mass(T1, X_air)
            + far_ * ca.h_mass(fuel_.T, X_ch4)
            - (1.0 + far_) * ca.h_mass(t4_, x3_)
        )
        qin_ = far_ * ca.fuel_lhv_mass(X_ch4)
        wc_ = dh_same_X(T1, t2_, X_air)
        wt_ = wnet_ + wc_
        eta_ = wnet_ / qin_
        bwr_ = wc_ / wt_
        pr_vals.append(pr)
        eta_vals.append(eta_ * 100)
        bwr_vals.append(bwr_ * 100)
        print(
            f"  {pr:5d}  {t2_:8.1f}  {t3_:8.1f}  {t4_:8.1f}  {eta_ * 100:9.2f}  {bwr_ * 100:8.1f}"
        )

    # =========================================================================
    # Plots: T-s diagram + eta vs PR
    # =========================================================================
    # Entropy at each state
    s1 = ca.s_mass(T1, X_air, P1)
    s2 = ca.s_mass(T2, X_air, P2)  # == s1 (isentropic)
    s3 = ca.s_mass(T3, X3, P3)
    s4 = ca.s_mass(T4, X3, P4)  # == s3 (isentropic)

    # Isobaric paths for T-s: interpolate s linearly between endpoint states.
    # Composition changes along both paths (combustion, exhaust mixing with atmosphere),
    # so anchoring to the exact s values at each state guarantees the cycle closes.
    T_23 = np.linspace(T2, T3, 200)
    s_23 = s2 + (s3 - s2) * (T_23 - T2) / (T3 - T2)
    T_41 = np.linspace(T4, T1, 200)
    s_41 = s4 + (s1 - s4) * (T_41 - T4) / (T1 - T4)

    fig, (ax_ts, ax_eta) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"Brayton Cycle  (PR={PR:.0f}, phi={phi:.2f}, eta={eta_th * 100:.1f} %)", fontsize=13
    )

    # --- T-s diagram ---
    ax_ts.plot([s1, s2], [T1, T2], label="1-2 Compression")
    ax_ts.plot(s_23, T_23, label="2-3 Combustion")
    ax_ts.plot([s3, s4], [T3, T4], label="3-4 Expansion")
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

    # --- eta and BWR vs PR ---
    ax_eta.plot(pr_vals, eta_vals, marker="o", label="Thermal efficiency")
    ax_eta.plot(pr_vals, bwr_vals, marker="s", linestyle="--", label="Back-work ratio")
    ax_eta.axvline(PR, color="gray", linestyle=":", alpha=0.6, label=f"Design PR={PR:.0f}")
    ax_eta.set_xlabel("Pressure ratio [-]")
    ax_eta.set_ylabel("[%]")
    ax_eta.set_title("Performance vs Pressure Ratio")
    ax_eta.legend()
    ax_eta.grid(True, alpha=0.3)

    plt.tight_layout()
    if plt.get_backend().lower() == "agg":
        out_path = pathlib.Path(__file__).parent / "gas_turbine_cycle.png"
        plt.savefig(out_path, dpi=150)
        print(f"\nPlot saved to '{out_path}'")
    else:
        plt.show()


if __name__ == "__main__":
    main()
