#!/usr/bin/env python3
"""Combustor liner convective cooling design with bypass sweep.

Flow topology (counterflow convective cooling):

  Compressor exit (T2, P2)
       |
       +--[f_cool]--> liner cooling channels --> T_cool_exit
       |
       +--[1-f_cool]--> bypass
       |
       +--[fuel]
       |
       [mix: coolant_exit + bypass + fuel] --> combustor inlet
                                                    |
                                               complete_combustion
                                                    |
                                               COT, X_hot
                                                    |
                                               turbine inlet

Counterflow pressure budget:
  dP_total/P = dP_burner/P + dP_cool/P   (both in series)

Design sweep: vary bypass fraction f_cool to find minimum cooling flow
that keeps T_metal < 800 °C for three channel geometries:
  - Smooth
  - Rib-roughened  (Han et al. 1988, e/D=0.07, P/e=8, alpha=60°)
  - Dimpled        (Chyu et al. 1997, d/Dh=0.2, h/d=0.2, S/d=2)

Wall configurations compared at the design point:
  - Bare Inconel 718
  - Inconel 718 + YSZ APS TBC (fresh and aged 1000 h)
  - Inconel 718 + YSZ EB-PVD TBC

Key combaero tools used:
  - humid_air_composition          (realistic coolant composition)
  - mix + complete_combustion      (coupled combustor inlet state)
  - channel_smooth / channel_ribbed / channel_dimpled
  - wall_temperature_profile       (multi-layer T profile)
  - k_inconel718 / k_tbc_ysz      (temperature-dependent material k)
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
# Helpers
# ---------------------------------------------------------------------------


def hot_side_htc(T_hot: float, P: float, X: list, D_h: float, u: float) -> float:
    """Hot-gas-side HTC via Dittus-Boelter on combustor annulus."""
    rho = ca.density(T_hot, P, X)
    mu = ca.viscosity(T_hot, P, X)
    k = ca.thermal_conductivity(T_hot, P, X)
    Re = rho * u * D_h / mu
    Nu = ca.nusselt_dittus_boelter(Re, ca.prandtl(T_hot, P, X), heating=False)
    return ca.htc_from_nusselt(Nu, k, D_h)


def wall_temps(
    T_hot: float,
    T_cool: float,
    h_hot: float,
    h_cool: float,
    layers: list[tuple[float, float]],
) -> tuple[list[float], float]:
    """Interface temperatures and heat flux through a multi-layer wall."""
    t_over_k = np.array([t / k for t, k in layers])
    temps, q = ca.wall_temperature_profile(T_hot, T_cool, h_hot, h_cool, t_over_k)
    return list(temps), float(q)


def channel_result(
    name: str,
    T: float,
    P: float,
    X: list,
    u: float,
    D_h: float,
    L: float,
) -> ca.ChannelResult:
    """Dispatch to the right channel function by name."""
    if name == "Smooth":
        return ca.channel_smooth(T, P, X, u, D_h, L)
    if name == "Ribbed":
        return ca.channel_ribbed(T, P, X, u, D_h, L, e_D=0.07, P_e=8.0, alpha=60.0)
    if name == "Dimpled":
        return ca.channel_dimpled(T, P, X, u, D_h, L, d_Dh=0.20, h_d=0.20, S_d=2.0)
    raise ValueError(f"Unknown channel type: {name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    sp = SpeciesLocator.from_core()

    print("=" * 70)
    print("Combustor Liner Convective Cooling Design — Bypass Sweep")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Fixed operating conditions
    # -----------------------------------------------------------------------
    PR = 20.0
    P2 = PR * 101325.0  # Pa  compressor delivery pressure
    T2 = 667.0  # K   compressor exit (isentropic, PR=20)
    RH_in = 0.60  # relative humidity at compressor inlet
    dP_burner_frac = 0.04  # burner dP/P (4 %, typical modern GT)
    phi_overall = 0.50  # overall equivalence ratio -> COT ~1758 K
    D_ann = 0.30  # m   combustor annulus hydraulic diameter
    u_hot_ann = 25.0  # m/s mean hot-gas velocity in annulus

    X_cool = ca.humid_air_composition(288.15, 101325.0, RH_in)
    X_air = ca.standard_dry_air_composition()
    X_ch4 = sp.empty()
    X_ch4[sp.indices["CH4"]] = 1.0

    # Fuel flow fixed for phi_overall on total air (1 kg/s reference)
    mdot_total = 1.0
    air_ref = ca.Stream()
    air_ref.T, air_ref.P, air_ref.X, air_ref.mdot = T2, P2, X_air, mdot_total
    fuel_ref = ca.Stream()
    fuel_ref.T, fuel_ref.P, fuel_ref.X = 300.0, P2, X_ch4
    fuel_ref = ca.set_fuel_stream_for_phi(phi_overall, fuel_ref, air_ref)
    mdot_fuel = fuel_ref.mdot

    # -----------------------------------------------------------------------
    # Channel geometry
    # -----------------------------------------------------------------------
    D_ch = 0.008  # m   hydraulic diameter of each cooling channel
    L_ch = 0.35  # m   channel length (liner axial extent)
    N_ch = 8  # number of parallel channels per liner panel
    A_ch_total = N_ch * np.pi * (D_ch / 2) ** 2  # total flow area

    # -----------------------------------------------------------------------
    # Materials
    # -----------------------------------------------------------------------
    T_metal_limit = 1073.15  # K = 800 °C
    k_inconel = ca.k_inconel718(900.0)
    k_ysz_fresh = ca.k_tbc_ysz(1200.0, hours=0.0, is_ebpvd=False)
    k_ysz_aged = ca.k_tbc_ysz(1200.0, hours=1000.0, is_ebpvd=False)
    k_ysz_ebpvd = ca.k_tbc_ysz(1200.0, hours=0.0, is_ebpvd=True)
    t_wall = 0.003  # m  Inconel liner thickness
    t_tbc = 0.0003  # m  TBC thickness (300 µm)

    wall_configs = {
        "Inconel only": [(t_wall, k_inconel)],
        "Inconel + YSZ APS (fresh)": [(t_tbc, k_ysz_fresh), (t_wall, k_inconel)],
        "Inconel + YSZ APS (1000 h)": [(t_tbc, k_ysz_aged), (t_wall, k_inconel)],
        "Inconel + YSZ EB-PVD": [(t_tbc, k_ysz_ebpvd), (t_wall, k_inconel)],
    }

    print(f"\nAvailable materials : {ca.list_materials()}")
    print(f"Inconel 718  k      : {k_inconel:.2f} W/(m·K)  @ 900 K")
    print(f"YSZ APS fresh k     : {k_ysz_fresh:.3f} W/(m·K)  @ 1200 K")
    print(f"YSZ APS 1000h k     : {k_ysz_aged:.3f} W/(m·K)  @ 1200 K")
    print(f"YSZ EB-PVD k        : {k_ysz_ebpvd:.3f} W/(m·K)  @ 1200 K")
    print(f"\nLiner: {t_wall * 1e3:.1f} mm Inconel,  TBC: {t_tbc * 1e6:.0f} µm YSZ")
    print(f"Channels: {N_ch} parallel, D_h={D_ch * 1e3:.1f} mm, L={L_ch * 1e3:.0f} mm")
    print(f"Metal limit: {T_metal_limit - 273.15:.0f} °C")

    # -----------------------------------------------------------------------
    # Bypass sweep: f_cool = 0.05 … 0.35
    # -----------------------------------------------------------------------
    f_cool_vals = np.linspace(0.05, 0.35, 25)
    channel_names = ["Smooth", "Ribbed", "Dimpled"]

    # Storage: [channel][f_cool_idx] -> dict
    sweep: dict[str, list[dict]] = {n: [] for n in channel_names}

    rho_cool = ca.density(T2, P2, X_cool)

    for f_cool in f_cool_vals:
        mdot_cool = f_cool * mdot_total
        mdot_bypass = (1.0 - f_cool) * mdot_total
        u_cool = mdot_cool / (rho_cool * A_ch_total)

        # Combustor inlet: mix coolant exit (≈T2, small heating) + bypass + fuel
        cool_s = ca.Stream()
        cool_s.T, cool_s.P, cool_s.X, cool_s.mdot = T2, P2, X_cool, mdot_cool
        byp_s = ca.Stream()
        byp_s.T, byp_s.P, byp_s.X, byp_s.mdot = T2, P2, X_air, mdot_bypass
        fuel_s = ca.Stream()
        fuel_s.T, fuel_s.P, fuel_s.X, fuel_s.mdot = 300.0, P2, X_ch4, mdot_fuel
        mixed_in = ca.mix([cool_s, byp_s, fuel_s])
        burned = ca.complete_combustion(mixed_in.T, mixed_in.X, mixed_in.P)
        T_hot = burned.T
        X_hot = burned.X
        P_hot = P2 * (1.0 - dP_burner_frac)

        h_hot = hot_side_htc(T_hot, P_hot, X_hot, D_ann, u_hot_ann)

        for ch_name in channel_names:
            r = channel_result(ch_name, T2, P2, X_cool, u_cool, D_ch, L_ch)

            # Wall temperatures for each wall config
            wall_res = {}
            for cfg_name, layers in wall_configs.items():
                temps, q = wall_temps(T_hot, T2, h_hot, r.h, layers)
                T_hot_surf = temps[0]
                T_metal_hot = temps[1] if len(temps) > 2 else temps[0]
                wall_res[cfg_name] = {
                    "T_hot_surf": T_hot_surf,
                    "T_metal_hot": T_metal_hot,
                    "q": q,
                }

            dP_total_frac = (dP_burner_frac * P2 + r.dP) / P2 * 100.0

            sweep[ch_name].append(
                {
                    "f_cool": f_cool,
                    "u_cool": u_cool,
                    "h_cool": r.h,
                    "Re": r.Re,
                    "Nu": r.Nu,
                    "dP_cool": r.dP,
                    "dP_total_frac": dP_total_frac,
                    "T_hot": T_hot,
                    "h_hot": h_hot,
                    "wall": wall_res,
                }
            )

    # -----------------------------------------------------------------------
    # Print summary table at f_cool = 0.10, 0.15, 0.20
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print("Bypass Sweep Summary  (Inconel only, no TBC)")
    print(f"{'=' * 70}")
    hdr = f"  {'f_cool':>7}  {'u [m/s]':>8}  {'h_cool':>7}  {'COT [K]':>8}  "
    hdr += f"{'T_metal [°C]':>13}  {'dP_cool [kPa]':>14}  {'dP_tot/P [%]':>13}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    for ch_name in channel_names:
        print(f"\n  Channel: {ch_name}")
        for pt in sweep[ch_name]:
            T_m = pt["wall"]["Inconel only"]["T_metal_hot"]
            ok = "" if T_m < T_metal_limit else " !"
            print(
                f"  {pt['f_cool']:7.2f}  {pt['u_cool']:8.1f}  {pt['h_cool']:7.0f}"
                f"  {pt['T_hot']:8.1f}  {T_m - 273.15:13.1f}{ok}"
                f"  {pt['dP_cool'] / 1e3:14.2f}  {pt['dP_total_frac']:13.2f}"
            )

    # -----------------------------------------------------------------------
    # Find minimum f_cool to meet T_metal < 800 °C per channel × wall config
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print("Minimum bypass fraction to achieve T_metal < 800 °C")
    print(f"{'=' * 70}")
    print(f"  {'Channel':10s}  {'Wall config':35s}  {'f_cool_min':>10}  {'dP_tot/P [%]':>13}")
    print(f"  {'-' * 10}  {'-' * 35}  {'-' * 10}  {'-' * 13}")

    design_points: dict[str, dict[str, dict]] = {}
    for ch_name in channel_names:
        design_points[ch_name] = {}
        for cfg_name in wall_configs:
            f_min = None
            for pt in sweep[ch_name]:
                T_m = pt["wall"][cfg_name]["T_metal_hot"]
                if T_m < T_metal_limit:
                    f_min = pt["f_cool"]
                    dp = pt["dP_total_frac"]
                    design_points[ch_name][cfg_name] = pt
                    break
            if f_min is not None:
                print(f"  {ch_name:10s}  {cfg_name:35s}  {f_min:10.2f}  {dp:13.2f}")
            else:
                print(f"  {ch_name:10s}  {cfg_name:35s}  {'NOT MET':>10}  {'---':>13}")

    # -----------------------------------------------------------------------
    # Plots
    # -----------------------------------------------------------------------
    colors = {"Smooth": "steelblue", "Ribbed": "darkorange", "Dimpled": "green"}
    cfg_plot = "Inconel only"  # wall config shown in sweep plots

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Combustor Liner Cooling — Bypass Sweep\n"
        f"PR={PR:.0f}, phi={phi_overall:.2f}, {N_ch} channels D_h={D_ch * 1e3:.0f} mm, "
        f"dP_burner/P={dP_burner_frac * 100:.0f}%",
        fontsize=12,
    )

    # --- Plot 1: T_metal vs f_cool ---
    ax = axes[0, 0]
    for ch_name in channel_names:
        f_vals = [pt["f_cool"] for pt in sweep[ch_name]]
        T_vals = [pt["wall"][cfg_plot]["T_metal_hot"] - 273.15 for pt in sweep[ch_name]]
        ax.plot(f_vals, T_vals, color=colors[ch_name], label=ch_name)
    ax.axhline(
        T_metal_limit - 273.15,
        color="red",
        linestyle="--",
        label=f"Limit {T_metal_limit - 273.15:.0f} °C",
    )
    ax.set_xlabel("Bypass fraction f_cool [-]")
    ax.set_ylabel("Metal hot-face temperature [°C]")
    ax.set_title(f"Metal Temperature ({cfg_plot})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Plot 2: Total dP/P vs f_cool ---
    ax = axes[0, 1]
    for ch_name in channel_names:
        f_vals = [pt["f_cool"] for pt in sweep[ch_name]]
        dp_vals = [pt["dP_total_frac"] for pt in sweep[ch_name]]
        ax.plot(f_vals, dp_vals, color=colors[ch_name], label=ch_name)
    ax.axhline(
        dP_burner_frac * 100,
        color="gray",
        linestyle=":",
        label=f"Burner only {dP_burner_frac * 100:.0f}%",
    )
    ax.set_xlabel("Bypass fraction f_cool [-]")
    ax.set_ylabel("Total dP/P [%]")
    ax.set_title("Counterflow Pressure Budget\n(dP_burner + dP_cool) / P")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Plot 3: COT vs f_cool (coupling effect) ---
    ax = axes[1, 0]
    for ch_name in channel_names:
        f_vals = [pt["f_cool"] for pt in sweep[ch_name]]
        cot_vals = [pt["T_hot"] for pt in sweep[ch_name]]
        ax.plot(f_vals, cot_vals, color=colors[ch_name], label=ch_name)
    ax.set_xlabel("Bypass fraction f_cool [-]")
    ax.set_ylabel("Combustor outlet temperature [K]")
    ax.set_title("COT vs Bypass Fraction\n(more bypass → slightly lower COT)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Plot 4: T_metal at design f_cool_min per wall config ---
    ax = axes[1, 1]
    cfg_names = list(wall_configs.keys())
    x = np.arange(len(cfg_names))
    width = 0.25
    for i, ch_name in enumerate(channel_names):
        T_vals = []
        for cfg_name in cfg_names:
            if cfg_name in design_points[ch_name]:
                pt = design_points[ch_name][cfg_name]
                T_vals.append(pt["wall"][cfg_name]["T_metal_hot"] - 273.15)
            else:
                T_vals.append(float("nan"))
        ax.bar(x + i * width, T_vals, width, color=colors[ch_name], label=ch_name, alpha=0.85)
    ax.axhline(
        T_metal_limit - 273.15,
        color="red",
        linestyle="--",
        label=f"Limit {T_metal_limit - 273.15:.0f} °C",
    )
    ax.set_xticks(x + width)
    ax.set_xticklabels(
        [c.replace(" + ", "\n+\n").replace(" (", "\n(") for c in cfg_names],
        fontsize=7,
    )
    ax.set_ylabel("Metal temperature at f_cool_min [°C]")
    ax.set_title("Metal Temp at Minimum Bypass\n(per channel × wall config)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    if plt.get_backend().lower() == "agg":
        out_path = pathlib.Path(__file__).parent / "combustor_liner_cooling.png"
        plt.savefig(out_path, dpi=150)
        print(f"\nPlot saved to '{out_path}'")
    else:
        plt.show()


if __name__ == "__main__":
    main()
