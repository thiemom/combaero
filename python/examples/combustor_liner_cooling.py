#!/usr/bin/env python3
"""Combustor liner convective cooling design.

Models a gas turbine combustor liner wall cooled by compressor delivery air
flowing through internal channels. Compares three surface geometries:

  - Smooth channel
  - Rib-roughened channel  (Han et al. 1988)
  - Dimpled channel        (Chyu et al. 1997)

And two wall configurations:

  - Bare Inconel 718 liner
  - Inconel 718 + YSZ thermal barrier coating (APS, fresh and aged 1000 h)

Operating point: modern 20-bar gas turbine combustor.

  Coolant : humid compressor delivery air (RH=60% at inlet, T2~667 K, P=20 bar)
  Hot gas : combustion products at COT ~1758 K (phi~0.50, CH4/air)
  Hot-side HTC : estimated from Dittus-Boelter on annular combustor passage
  Burner dP/P  : 4 % (typical modern combustor)

Key combaero tools used:
  - humid_air_composition          (realistic coolant composition)
  - complete_combustion            (hot-gas composition and temperature)
  - channel_smooth / channel_ribbed / channel_dimpled  (coolant-side HTC + dP)
  - wall_temperature_profile       (temperature through each wall layer)
  - overall_htc_wall_multilayer    (U-value for wall stack)
  - k_inconel718 / k_tbc_ysz      (temperature-dependent material properties)
  - list_materials                 (material database)
  - nusselt_dittus_boelter / htc_from_nusselt  (hot-side HTC)
"""

from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np

import combaero as ca
from combaero.species import SpeciesLocator

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def hot_side_htc(T_hot: float, P: float, X: list, D_h: float, u: float) -> float:
    """Estimate hot-gas-side HTC via Dittus-Boelter on the combustor annulus."""
    rho = ca.density(T_hot, P, X)
    mu = ca.viscosity(T_hot, P, X)
    k = ca.thermal_conductivity(T_hot, P, X)
    Re = rho * u * D_h / mu
    Nu = ca.nusselt_dittus_boelter(Re, ca.prandtl(T_hot, P, X), heating=False)
    return ca.htc_from_nusselt(Nu, k, D_h)


def wall_stack_temps(
    T_hot: float,
    T_cool: float,
    h_hot: float,
    h_cool: float,
    layers: list[tuple[float, float]],  # [(t, k), ...]
) -> tuple[list[float], float]:
    """Temperatures at each interface through the wall stack.

    Returns (interface_temps, q) where interface_temps includes
    T_hot_side_surface and T_cool_side_surface.
    """
    t_over_k = np.array([t / k for t, k in layers])
    temps, q = ca.wall_temperature_profile(T_hot, T_cool, h_hot, h_cool, t_over_k)
    return list(temps), float(q)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    sp = SpeciesLocator.from_core()

    print("=" * 70)
    print("Combustor Liner Convective Cooling Design")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Materials
    # -----------------------------------------------------------------------
    print(f"\nAvailable materials: {ca.list_materials()}")

    T_metal_limit = 1073.15  # K = 800 °C  (Inconel 718 oxidation limit)
    T_ref_metal = 900.0  # K  representative Inconel temperature for k lookup
    T_ref_tbc = 1200.0  # K  TBC outer surface temperature (hotter than metal)

    k_inconel = ca.k_inconel718(T_ref_metal)
    k_ysz_fresh = ca.k_tbc_ysz(T_ref_tbc, hours=0.0, is_ebpvd=False)  # APS as-sprayed
    k_ysz_aged = ca.k_tbc_ysz(T_ref_tbc, hours=1000.0, is_ebpvd=False)  # APS after 1000 h
    k_ysz_ebpvd = ca.k_tbc_ysz(T_ref_tbc, hours=0.0, is_ebpvd=True)

    t_wall = 0.003  # m  Inconel liner thickness
    t_tbc = 0.0003  # m  TBC thickness (300 µm typical APS)

    print("\n--- Wall Materials ---")
    print(f"  Inconel 718        : k = {k_inconel:.2f} W/(m·K)  @ {T_ref_metal:.0f} K")
    print(f"  YSZ APS  (fresh)   : k = {k_ysz_fresh:.3f} W/(m·K)  @ {T_ref_tbc:.0f} K")
    print(f"  YSZ APS  (1000 h)  : k = {k_ysz_aged:.3f} W/(m·K)  @ {T_ref_tbc:.0f} K")
    print(f"  YSZ EB-PVD (fresh) : k = {k_ysz_ebpvd:.3f} W/(m·K)  @ {T_ref_tbc:.0f} K")
    print(f"  Liner thickness    : {t_wall * 1e3:.1f} mm")
    print(f"  TBC thickness      : {t_tbc * 1e6:.0f} µm")

    # -----------------------------------------------------------------------
    # Operating conditions
    # -----------------------------------------------------------------------
    PR = 20.0
    P_cool = PR * 101325.0  # Pa  compressor delivery pressure
    T_cool = 667.0  # K   compressor exit temperature (isentropic, PR=20)
    RH_in = 0.60  # relative humidity at compressor inlet
    T_inlet = 288.15  # K   ISA sea level
    P_inlet = 101325.0  # Pa

    X_cool = ca.humid_air_composition(T_inlet, P_inlet, RH_in)

    # Hot gas: CH4/air combustion at phi=0.50 -> COT ~1758 K
    # (COT = Combustor Outlet Temperature; TIT is lower after dilution cooling)
    phi_comb = 0.50
    X_air = ca.standard_dry_air_composition()
    X_ch4 = sp.empty()
    X_ch4[sp.indices["CH4"]] = 1.0

    air = ca.Stream()
    air.T, air.P, air.X, air.mdot = T_cool, P_cool, X_air, 1.0
    fuel = ca.Stream()
    fuel.T, fuel.P, fuel.X = 300.0, P_cool, X_ch4
    fuel = ca.set_fuel_stream_for_phi(phi_comb, fuel, air)
    mixed = ca.mix([fuel, air])
    burned = ca.complete_combustion(mixed.T, mixed.X, mixed.P)

    T_hot = burned.T
    X_hot = burned.X
    P_hot = P_cool * (1.0 - 0.04)  # 4% burner pressure loss

    print("\n--- Operating Conditions ---")
    print(f"  Compressor exit  : T = {T_cool:.1f} K,  P = {P_cool / 1e5:.1f} bar")
    print(f"  Coolant          : humid air (RH={RH_in:.0%} at inlet)")
    print(f"  Combustor exit   : T = {T_hot:.1f} K  (COT),  phi = {phi_comb:.2f}")
    print(f"  Burner dP/P      : 4 %   ->  P_hot = {P_hot / 1e5:.2f} bar")
    print(f"  Metal limit      : {T_metal_limit - 273.15:.0f} °C")

    # -----------------------------------------------------------------------
    # Geometry
    # -----------------------------------------------------------------------
    # Combustor annulus (hot side)
    D_annulus = 0.30  # m  hydraulic diameter of combustor annulus
    u_hot = 25.0  # m/s  mean hot-gas velocity in combustor

    # Cooling channel geometry (per liner panel)
    D_ch = 0.008  # m  hydraulic diameter of cooling channel
    L_ch = 0.35  # m  channel length (liner axial extent)
    u_cool = 60.0  # m/s  coolant velocity

    # -----------------------------------------------------------------------
    # Hot-side HTC
    # -----------------------------------------------------------------------
    h_hot = hot_side_htc(T_hot, P_hot, X_hot, D_annulus, u_hot)
    print("\n--- Hot-Side Convection ---")
    print(f"  D_annulus = {D_annulus * 1e3:.0f} mm,  u_hot = {u_hot:.0f} m/s")
    print(f"  h_hot     = {h_hot:.1f} W/(m²·K)")

    # -----------------------------------------------------------------------
    # Coolant-side channels
    # -----------------------------------------------------------------------
    channels = {
        "Smooth": ca.channel_smooth(T_cool, P_cool, X_cool, u_cool, D_ch, L_ch),
        "Ribbed (e/D=0.07, P/e=8, 60°)": ca.channel_ribbed(
            T_cool, P_cool, X_cool, u_cool, D_ch, L_ch, e_D=0.07, P_e=8.0, alpha=60.0
        ),
        "Dimpled (d/Dh=0.2, h/d=0.2, S/d=2)": ca.channel_dimpled(
            T_cool, P_cool, X_cool, u_cool, D_ch, L_ch, d_Dh=0.20, h_d=0.20, S_d=2.0
        ),
    }

    print("\n--- Coolant Channel Performance ---")
    print(f"  D_h = {D_ch * 1e3:.1f} mm,  L = {L_ch * 1e3:.0f} mm,  u = {u_cool:.0f} m/s")
    print(f"\n  {'Channel':40s}  {'Re':>8}  {'Nu':>6}  {'h [W/m²K]':>11}  {'dP [kPa]':>9}")
    print(f"  {'-' * 40}  {'-' * 8}  {'-' * 6}  {'-' * 11}  {'-' * 9}")
    for name, r in channels.items():
        print(f"  {name:40s}  {r.Re:8.0f}  {r.Nu:6.1f}  {r.h:11.1f}  {r.dP / 1e3:9.2f}")

    # -----------------------------------------------------------------------
    # Wall temperature analysis: 4 configurations per channel type
    # -----------------------------------------------------------------------
    wall_configs = {
        "Inconel only": [(t_wall, k_inconel)],
        "Inconel + YSZ APS (fresh)": [(t_tbc, k_ysz_fresh), (t_wall, k_inconel)],
        "Inconel + YSZ APS (1000 h)": [(t_tbc, k_ysz_aged), (t_wall, k_inconel)],
        "Inconel + YSZ EB-PVD": [(t_tbc, k_ysz_ebpvd), (t_wall, k_inconel)],
    }

    print(f"\n--- Wall Temperature Analysis (T_hot={T_hot:.0f} K, T_cool={T_cool:.0f} K) ---")
    print("  Layers listed hot-side first.")
    print()

    # Store results for plotting
    results: dict[str, dict] = {}

    for ch_name, ch_res in channels.items():
        h_cool = ch_res.h
        results[ch_name] = {}
        print(f"  Channel: {ch_name}")
        print(
            f"  {'Wall config':35s}  {'T_hot-surf [°C]':>16}  {'T_metal [°C]':>13}  {'q [kW/m²]':>10}  {'OK?':>5}"
        )
        print(f"  {'-' * 35}  {'-' * 16}  {'-' * 13}  {'-' * 10}  {'-' * 5}")
        for cfg_name, layers in wall_configs.items():
            temps, q = wall_stack_temps(T_hot, T_cool, h_hot, h_cool, layers)
            # temps[0] = hot-gas-side wall surface, temps[-1] = coolant-side surface
            T_hot_surf = temps[0]  # hot face (TBC outer surface if TBC present)
            # temps layout: [T_hot_wall, ...interfaces..., T_cool_wall]
            # With TBC: [T_tbc_hot, T_tbc_cool=T_metal_hot, T_metal_cool]
            # Without:  [T_metal_hot, T_metal_cool]
            T_metal_hot = temps[1] if len(temps) > 2 else temps[0]
            ok = "OK" if T_metal_hot < T_metal_limit else "FAIL"
            results[ch_name][cfg_name] = {
                "temps": temps,
                "q": q,
                "T_hot_surf": T_hot_surf,
                "T_metal_hot": T_metal_hot,
            }
            print(
                f"  {cfg_name:35s}  {T_hot_surf - 273.15:16.1f}  "
                f"{T_metal_hot - 273.15:13.1f}  {q / 1e3:10.2f}  {ok:>5}"
            )
        print()

    # -----------------------------------------------------------------------
    # Pressure budget
    # -----------------------------------------------------------------------
    dP_burner = 0.04 * P_cool
    print("--- Pressure Budget ---")
    print(f"  Burner dP/P = 4%  ->  dP_burner = {dP_burner / 1e3:.1f} kPa")
    print(f"\n  {'Channel':40s}  {'dP_cool [kPa]':>14}  {'dP_cool/dP_burner':>18}")
    print(f"  {'-' * 40}  {'-' * 14}  {'-' * 18}")
    for name, r in channels.items():
        frac = r.dP / dP_burner
        print(f"  {name:40s}  {r.dP / 1e3:14.2f}  {frac:18.3f}")

    # -----------------------------------------------------------------------
    # Plots
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(
        f"Combustor Liner Cooling  (COT={T_hot:.0f} K, PR={PR:.0f}, phi={phi_comb:.2f})",
        fontsize=13,
    )

    ch_names = list(channels.keys())
    cfg_names = list(wall_configs.keys())
    x_ch = np.arange(len(ch_names))
    width = 0.2

    # --- Plot 1: Metal hot-face temperature per channel × wall config ---
    ax = axes[0]
    for i, cfg in enumerate(cfg_names):
        T_vals = [results[ch][cfg]["T_metal_hot"] - 273.15 for ch in ch_names]
        ax.bar(x_ch + i * width, T_vals, width, label=cfg)
    ax.axhline(
        T_metal_limit - 273.15,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Limit {T_metal_limit - 273.15:.0f} °C",
    )
    ax.set_xticks(x_ch + width * 1.5)
    ax.set_xticklabels([n.split("(")[0].strip() for n in ch_names], fontsize=9)
    ax.set_ylabel("Metal hot-face temperature [°C]")
    ax.set_title("Metal Temperature")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3, axis="y")

    # --- Plot 2: Heat flux per channel × wall config ---
    ax = axes[1]
    for i, cfg in enumerate(cfg_names):
        q_vals = [results[ch][cfg]["q"] / 1e3 for ch in ch_names]
        ax.bar(x_ch + i * width, q_vals, width, label=cfg)
    ax.set_xticks(x_ch + width * 1.5)
    ax.set_xticklabels([n.split("(")[0].strip() for n in ch_names], fontsize=9)
    ax.set_ylabel("Wall heat flux [kW/m²]")
    ax.set_title("Heat Flux Through Wall")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3, axis="y")

    # --- Plot 3: Pressure drop comparison ---
    ax = axes[2]
    dP_vals = [channels[n].dP / 1e3 for n in ch_names]
    bars = ax.bar(x_ch, dP_vals, 0.5, color=["steelblue", "darkorange", "green"])
    ax.axhline(
        dP_burner / 1e3,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Burner dP = {dP_burner / 1e3:.1f} kPa",
    )
    ax.set_xticks(x_ch)
    ax.set_xticklabels([n.split("(")[0].strip() for n in ch_names], fontsize=9)
    ax.set_ylabel("Coolant channel dP [kPa]")
    ax.set_title("Pressure Drop")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")
    for bar, val in zip(bars, dP_vals, strict=True):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            val + 0.5,
            f"{val:.1f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    plt.tight_layout()
    if plt.get_backend().lower() == "agg":
        out_path = pathlib.Path(__file__).parent / "combustor_liner_cooling.png"
        plt.savefig(out_path, dpi=150)
        print(f"\nPlot saved to '{out_path}'")
    else:
        plt.show()


if __name__ == "__main__":
    main()
