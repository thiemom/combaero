#!/usr/bin/env python3
"""Toy industrial combustor cooling example.

Demonstrates a counterflow convective-cooled combustor can with an effusion
front plate and burner jets using coupled combustion, transport, and cooling
correlations.
"""

from __future__ import annotations

from dataclasses import dataclass

import combaero as ca
from combaero.species import SpeciesLocator


@dataclass(frozen=True)
class CombustorInputs:
    pressure_pa: float = 15.0e5
    air_inlet_temperature_k: float = 720.0
    fuel_inlet_temperature_k: float = 310.0
    mdot_air_total_kg_s: float = 9.5
    phi_primary: float = 0.78


@dataclass(frozen=True)
class LinerGeometry:
    can_diameter_m: float = 0.36
    can_length_m: float = 0.52
    wall_thickness_m: float = 0.0025
    wall_conductivity_w_mk: float = 18.0
    cooling_channel_dh_m: float = 0.011


@dataclass(frozen=True)
class WallLayer:
    name: str
    thickness_m: float
    conductivity_w_mk: float


def _make_air_stream(inputs: CombustorInputs) -> ca.Stream:
    stream = ca.Stream()
    stream.T = inputs.air_inlet_temperature_k
    stream.P = inputs.pressure_pa
    stream.X = ca.standard_dry_air_composition()
    stream.mdot = inputs.mdot_air_total_kg_s
    return stream


def _make_fuel_template(inputs: CombustorInputs, sp: SpeciesLocator) -> ca.Stream:
    fuel_x = sp.empty()
    fuel_x[sp.indices["CH4"]] = 0.92
    fuel_x[sp.indices["C2H6"]] = 0.05
    fuel_x[sp.indices["N2"]] = 0.03

    stream = ca.Stream()
    stream.T = inputs.fuel_inlet_temperature_k
    stream.P = inputs.pressure_pa
    stream.X = fuel_x
    stream.mdot = 0.0  # set by set_fuel_stream_for_phi
    return stream


def _hot_gas_side_htc(
    gas_temperature_k: float,
    pressure_pa: float,
    gas_x: list[float],
    velocity_m_s: float,
    hydraulic_diameter_m: float,
) -> tuple[float, float, float]:
    re = ca.reynolds(gas_temperature_k, pressure_pa, gas_x, velocity_m_s, hydraulic_diameter_m)
    pr = ca.prandtl(gas_temperature_k, pressure_pa, gas_x)
    k = ca.thermal_conductivity(gas_temperature_k, pressure_pa, gas_x)
    nu = ca.nusselt_dittus_boelter(re, pr, heating=False)
    h = ca.htc_from_nusselt(nu, k, hydraulic_diameter_m)
    return h, re, nu


def _coolant_side_htc(
    coolant_temperature_k: float,
    pressure_pa: float,
    coolant_x: list[float],
    velocity_m_s: float,
    channel_dh_m: float,
) -> tuple[float, float, float]:
    re = ca.reynolds(coolant_temperature_k, pressure_pa, coolant_x, velocity_m_s, channel_dh_m)
    pr = ca.prandtl(coolant_temperature_k, pressure_pa, coolant_x)
    k = ca.thermal_conductivity(coolant_temperature_k, pressure_pa, coolant_x)
    nu = ca.nusselt_dittus_boelter(re, pr, heating=True)
    h = ca.htc_from_nusselt(nu, k, channel_dh_m)
    return h, re, nu


def _slice_wall_temperature_profile(
    layers: list[WallLayer],
    interface_temperatures_k: list[float],
    slices_per_layer: int,
) -> list[tuple[float, float, str]]:
    samples: list[tuple[float, float, str]] = []
    x_cursor_m = 0.0

    for idx, layer in enumerate(layers):
        t_start = interface_temperatures_k[idx]
        t_end = interface_temperatures_k[idx + 1]

        for i in range(slices_per_layer + 1):
            if idx > 0 and i == 0:
                continue
            frac = i / slices_per_layer
            x_m = x_cursor_m + frac * layer.thickness_m
            t_k = t_start + frac * (t_end - t_start)
            samples.append((x_m, t_k, layer.name))

        x_cursor_m += layer.thickness_m

    return samples


def _combustor_pressure_loss_budget(
    inputs: CombustorInputs,
    geom: LinerGeometry,
    gas_temperature_k: float,
    gas_x: list[float],
    v_hot_m_s: float,
    v_cool_m_s: float,
    coolant_x: list[float],
) -> dict[str, float]:
    rho_hot = ca.density(gas_temperature_k, inputs.pressure_pa, gas_x)
    q_dyn_hot = 0.5 * rho_hot * v_hot_m_s**2

    # Toy assumption: burner swirler/mixer local loss coefficient
    k_mixer = 3.2
    dp_burner_mixer = k_mixer * q_dyn_hot

    dp_primary_friction, re_primary, f_primary = ca.pressure_drop_pipe(
        T=gas_temperature_k,
        P=inputs.pressure_pa,
        X=gas_x,
        v=v_hot_m_s,
        D=geom.can_diameter_m,
        L=0.70 * geom.can_length_m,
        roughness=2.0e-5,
        correlation="haaland",
    )

    k_dilution = 0.8
    dp_dilution_entry = k_dilution * q_dyn_hot

    dp_coolant_friction, re_cooling, f_cooling = ca.pressure_drop_pipe(
        T=inputs.air_inlet_temperature_k + 40.0,
        P=inputs.pressure_pa,
        X=coolant_x,
        v=v_cool_m_s,
        D=geom.cooling_channel_dh_m,
        L=geom.can_length_m,
        roughness=1.0e-5,
        correlation="haaland",
    )

    dp_hot_total = dp_burner_mixer + dp_primary_friction + dp_dilution_entry

    return {
        "dp_burner_mixer": dp_burner_mixer,
        "dp_primary_friction": dp_primary_friction,
        "dp_dilution_entry": dp_dilution_entry,
        "dp_hot_total": dp_hot_total,
        "dp_hot_total_pct": 100.0 * dp_hot_total / inputs.pressure_pa,
        "dp_coolant_friction": dp_coolant_friction,
        "dp_coolant_pct": 100.0 * dp_coolant_friction / inputs.pressure_pa,
        "re_primary": re_primary,
        "f_primary": f_primary,
        "re_cooling": re_cooling,
        "f_cooling": f_cooling,
    }


def main() -> None:
    inputs = CombustorInputs()
    geom = LinerGeometry()
    sp = SpeciesLocator.from_core()

    air = _make_air_stream(inputs)
    fuel_template = _make_fuel_template(inputs, sp)

    # Burner block: direct phi-targeted fuel stream helper.
    fuel = ca.set_fuel_stream_for_phi(inputs.phi_primary, fuel_template, air)

    mixed = ca.mix([fuel, air], P_out=inputs.pressure_pa)
    burned = ca.combustion_equilibrium(mixed.T, mixed.X, mixed.P)

    # Toy airflow split
    mdot_cool_front = 0.06 * air.mdot  # front plate effusion
    mdot_cool_liner = 0.24 * air.mdot  # counterflow liner channel
    mdot_primary = mixed.mdot - mdot_cool_front - mdot_cool_liner

    print("=" * 84)
    print("COUNTERFLOW + EFFUSION COOLED COMBUSTOR CAN (TOY INDUSTRIAL CASE)")
    print("=" * 84)
    print(f"Pressure: {inputs.pressure_pa / 1e5:.2f} bar")
    print(f"Primary-zone phi: {inputs.phi_primary:.2f}")
    print(f"Primary mixed flow: {mdot_primary:.3f} kg/s")
    print(f"Adiabatic flame temperature: {burned.T:.1f} K")
    print()

    # Geometry for flow estimates
    perimeter = 3.141592653589793 * geom.can_diameter_m
    hot_flow_area = 0.25 * 3.141592653589793 * geom.can_diameter_m**2

    # Assume counterflow channel area sized by a gap around liner
    cooling_gap_m = 0.004
    cooling_area = perimeter * cooling_gap_m

    # Hot-gas and coolant velocities
    rho_hot = ca.density(burned.T, inputs.pressure_pa, burned.X)
    v_hot = mdot_primary / (rho_hot * hot_flow_area)

    t_cool_in = inputs.air_inlet_temperature_k + 40.0
    rho_cool = ca.density(t_cool_in, inputs.pressure_pa, air.X)
    v_cool = mdot_cool_liner / (rho_cool * cooling_area)

    h_hot, re_hot, nu_hot = _hot_gas_side_htc(
        burned.T,
        inputs.pressure_pa,
        burned.X,
        v_hot,
        geom.can_diameter_m,
    )
    h_cool, re_cool, nu_cool = _coolant_side_htc(
        t_cool_in,
        inputs.pressure_pa,
        air.X,
        v_cool,
        geom.cooling_channel_dh_m,
    )

    # Front-plate effusion estimate from burner-side jets
    # Use correlation-friendly representative design point
    effusion_eta = ca.effusion_effectiveness(
        x_D=8.0,
        M=2.0,
        DR=1.75,
        porosity=0.055,
        s_D=6.0,
        alpha_deg=30.0,
    )

    # Baseline front plate: bare Haynes 230 metal wall
    k_haynes = ca.k_haynes230(min(max(900.0, t_cool_in), 1400.0))
    q_front_bare = ca.cooled_wall_heat_flux(
        T_hot=burned.T,
        T_coolant=t_cool_in,
        h_hot=h_hot,
        h_coolant=h_cool,
        eta=effusion_eta,
        t_wall=geom.wall_thickness_m,
        k_wall=k_haynes,
    )

    taw_front = ca.adiabatic_wall_temperature(burned.T, t_cool_in, effusion_eta)

    # If needed, upgrade wall with YSZ coating on top of Haynes substrate
    needs_upgrade = q_front_bare / 1e3 > 160.0
    front_layers: list[WallLayer]
    if needs_upgrade:
        k_tbc = ca.k_tbc_ysz(T=min(taw_front, 1690.0), hours=2000.0, is_ebpvd=False)
        front_layers = [
            WallLayer("YSZ coating (APS)", 0.00035, k_tbc),
            WallLayer("Haynes 230 substrate", geom.wall_thickness_m, k_haynes),
        ]
    else:
        front_layers = [
            WallLayer("Haynes 230 substrate", geom.wall_thickness_m, k_haynes),
        ]

    t_over_k_front = [layer.thickness_m / layer.conductivity_w_mk for layer in front_layers]
    wall_temps_front, q_front = ca.wall_temperature_profile(
        T_hot=taw_front,
        T_cold=t_cool_in,
        h_hot=h_hot,
        h_cold=h_cool,
        t_over_k=t_over_k_front,
    )
    wall_slices = _slice_wall_temperature_profile(
        front_layers, wall_temps_front, slices_per_layer=4
    )

    q_liner = ca.cooled_wall_heat_flux(
        T_hot=burned.T,
        T_coolant=t_cool_in,
        h_hot=h_hot,
        h_coolant=h_cool,
        eta=0.0,
        t_wall=geom.wall_thickness_m,
        k_wall=k_haynes,
    )

    dp_budget = _combustor_pressure_loss_budget(
        inputs=inputs,
        geom=geom,
        gas_temperature_k=burned.T,
        gas_x=burned.X,
        v_hot_m_s=v_hot,
        v_cool_m_s=v_cool,
        coolant_x=air.X,
    )

    area_front = 0.25 * 3.141592653589793 * geom.can_diameter_m**2
    area_liner = perimeter * geom.can_length_m

    qdot_front = q_front * area_front
    qdot_liner = q_liner * area_liner

    print("Hot side:")
    print(f"  Velocity: {v_hot:.2f} m/s")
    print(f"  Re: {re_hot:.3e}")
    print(f"  Nu: {nu_hot:.1f}")
    print(f"  h_hot: {h_hot:.1f} W/(m^2*K)")
    print()

    print("Counterflow coolant side:")
    print(f"  Velocity: {v_cool:.2f} m/s")
    print(f"  Re: {re_cool:.3e}")
    print(f"  Nu: {nu_cool:.1f}")
    print(f"  h_cool: {h_cool:.1f} W/(m^2*K)")
    print()

    print("Effusion front plate:")
    print(f"  Cooling effectiveness eta: {effusion_eta:.3f}")
    print(f"  Adiabatic wall temperature: {taw_front:.1f} K")
    print(f"  Bare Haynes wall heat flux: {q_front_bare / 1e3:.1f} kW/m^2")
    print(f"  Wall upgrade needed: {'yes (add YSZ coating)' if needs_upgrade else 'no'}")
    print(f"  Heat flux (cooled): {q_front / 1e3:.1f} kW/m^2")
    print(f"  Heat flux (uncooled baseline): {q_liner / 1e3:.1f} kW/m^2")
    print()

    print("Front-wall layer stack:")
    for layer in front_layers:
        print(
            f"  - {layer.name}: t={layer.thickness_m * 1e3:.2f} mm, "
            f"k={layer.conductivity_w_mk:.3f} W/(m*K)"
        )
    print("  Sliced wall temperature profile (from hot surface into wall):")
    for x_m, t_k, layer_name in wall_slices:
        print(f"    x={x_m * 1e3:6.3f} mm, T={t_k:7.2f} K ({layer_name})")
    print()

    print("Toy combustor pressure-loss budget:")
    print("  Assumption: burner mixer K = 3.2")
    print(f"  Burner mixer loss: {dp_budget['dp_burner_mixer'] / 1e3:.3f} kPa")
    print(f"  Primary friction loss: {dp_budget['dp_primary_friction'] / 1e3:.3f} kPa")
    print(f"  Dilution-entry local loss: {dp_budget['dp_dilution_entry'] / 1e3:.3f} kPa")
    print(
        f"  Hot-side total: {dp_budget['dp_hot_total'] / 1e3:.3f} kPa "
        f"({dp_budget['dp_hot_total_pct']:.3f}% of combustor pressure)"
    )
    print(f"  Cooling-channel friction loss: {dp_budget['dp_coolant_friction'] / 1e3:.3f} kPa")
    print(f"  Primary Re={dp_budget['re_primary']:.3e}, f={dp_budget['f_primary']:.4f}")
    print(f"  Cooling Re={dp_budget['re_cooling']:.3e}, f={dp_budget['f_cooling']:.4f}")
    print()

    print("Integrated heat pickup:")
    print(f"  Front plate heat load: {qdot_front / 1e6:.3f} MW")
    print(f"  Liner heat load: {qdot_liner / 1e6:.3f} MW")
    print(f"  Total wall heat load: {(qdot_front + qdot_liner) / 1e6:.3f} MW")


if __name__ == "__main__":
    main()
