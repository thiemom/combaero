#!/usr/bin/env python3
"""
Heat Transfer Measurement Example
==================================

Practical example: Inferring heat flux from thermocouple measurements
in a furnace wall with known fluid boundary conditions.

Scenario:
---------
A furnace wall separates hot combustion gases (inside) from ambient air (outside).
The wall consists of:
  - 10mm steel liner (k = 50 W/(m*K))
  - 100mm ceramic insulation (k = 0.5 W/(m*K))
  - 5mm steel casing (k = 50 W/(m*K))

Thermocouples are installed at:
  - TC1: Hot-side steel surface (edge 0)
  - TC2: Steel-ceramic interface (edge 1)
  - TC3: Embedded 20mm into ceramic (depth = 30mm from hot surface)
  - TC4: Cold-side casing surface (edge 3)

Given the measured temperatures and known fluid conditions, we calculate
the heat flux through the wall.
"""

import numpy as np

import combaero as ca

print("=" * 70)
print("FURNACE WALL HEAT TRANSFER ANALYSIS")
print("=" * 70)

# ============================================================================
# Wall geometry and materials
# ============================================================================
print("\n--- Wall Construction ---")

# Layer properties
layers = [
    {"name": "Steel liner", "thickness": 0.010, "k": 50.0},
    {"name": "Ceramic insulation", "thickness": 0.100, "k": 0.5},
    {"name": "Steel casing", "thickness": 0.005, "k": 50.0},
]

for i, layer in enumerate(layers):
    print(f"  Layer {i + 1}: {layer['name']}")
    print(f"           t = {layer['thickness'] * 1000:.0f} mm, k = {layer['k']:.1f} W/(m*K)")

thicknesses = np.array([L["thickness"] for L in layers])
conductivities = np.array([L["k"] for L in layers])
t_over_k = thicknesses / conductivities

total_thickness = sum(thicknesses)
print(f"\n  Total wall thickness: {total_thickness * 1000:.0f} mm")

# ============================================================================
# Fluid boundary conditions (known from process)
# ============================================================================
print("\n--- Fluid Boundary Conditions ---")

T_gas = 1200 + 273.15  # Hot combustion gas temperature [K]
T_amb = 25 + 273.15  # Ambient air temperature [K]

# Convective heat transfer coefficients (from correlations or experience)
h_gas = 150.0  # Hot gas side (forced convection) [W/(m^2*K)]
h_amb = 10.0  # Ambient side (natural convection) [W/(m^2*K)]

print(f"  Hot gas:     T = {T_gas - 273.15:.0f} degC, h = {h_gas:.0f} W/(m^2*K)")
print(f"  Ambient air: T = {T_amb - 273.15:.0f} degC, h = {h_amb:.0f} W/(m^2*K)")

# ============================================================================
# Theoretical analysis (what we expect)
# ============================================================================
print("\n--- Theoretical Analysis ---")

# Overall heat transfer coefficient
U = ca.overall_htc(np.array([h_gas, h_amb]), t_over_k)
print(f"  Overall HTC: U = {U:.3f} W/(m^2*K)")

# Expected heat flux
dT = T_gas - T_amb
q_theory = ca.heat_flux(U, dT)
print(f"  Expected heat flux: q = {q_theory:.1f} W/m^2")

# Temperature profile
temps_theory, q_check = ca.wall_temperature_profile(T_gas, T_amb, h_gas, h_amb, t_over_k)

print("\n  Expected temperature profile:")
edge_names = ["Steel liner (hot)", "Steel-ceramic", "Ceramic-casing", "Casing (cold)"]
for i, (name, T) in enumerate(zip(edge_names, temps_theory, strict=False)):
    print(f"    Edge {i} ({name}): {T - 273.15:.1f} degC")

# ============================================================================
# Simulated measurements (with some noise)
# ============================================================================
print("\n--- Thermocouple Measurements ---")

# Add realistic measurement noise (+/-2 degC)
np.random.seed(42)
noise = np.random.normal(0, 2, 4)

# Measured temperatures
T_meas = {
    "TC1": temps_theory[0] + noise[0],  # Edge 0: hot surface
    "TC2": temps_theory[1] + noise[1],  # Edge 1: steel-ceramic interface
    "TC3": None,  # Will calculate for embedded TC
    "TC4": temps_theory[3] + noise[3],  # Edge 3: cold surface
}

# TC3 is embedded 20mm into ceramic (total depth = 10mm steel + 20mm = 30mm)
depth_TC3 = 0.030  # 30mm from hot surface
R_to_TC3 = 1 / h_gas + 0.010 / 50.0 + 0.020 / 0.5  # convective + steel + 20mm ceramic
T_TC3_true = T_gas - q_theory * R_to_TC3
T_meas["TC3"] = T_TC3_true + noise[2]

print("  Measured temperatures:")
print(f"    TC1 (hot surface):      {T_meas['TC1'] - 273.15:.1f} degC")
print(f"    TC2 (steel-ceramic):    {T_meas['TC2'] - 273.15:.1f} degC")
print(f"    TC3 (20mm into ceramic): {T_meas['TC3'] - 273.15:.1f} degC")
print(f"    TC4 (cold surface):     {T_meas['TC4'] - 273.15:.1f} degC")

# ============================================================================
# Heat flux inference from measurements
# ============================================================================
print("\n--- Heat Flux Inference ---")

# Method 1: From edge temperatures
q_from_TC1 = ca.heat_flux_from_T_at_edge(T_meas["TC1"], 0, T_gas, T_amb, h_gas, h_amb, t_over_k)
q_from_TC2 = ca.heat_flux_from_T_at_edge(T_meas["TC2"], 1, T_gas, T_amb, h_gas, h_amb, t_over_k)
q_from_TC4 = ca.heat_flux_from_T_at_edge(T_meas["TC4"], 3, T_gas, T_amb, h_gas, h_amb, t_over_k)

print("  From edge thermocouples:")
print(
    f"    TC1 -> q = {q_from_TC1:.1f} W/m^2 (error: {(q_from_TC1 - q_theory) / q_theory * 100:+.1f}%)"
)
print(
    f"    TC2 -> q = {q_from_TC2:.1f} W/m^2 (error: {(q_from_TC2 - q_theory) / q_theory * 100:+.1f}%)"
)
print(
    f"    TC4 -> q = {q_from_TC4:.1f} W/m^2 (error: {(q_from_TC4 - q_theory) / q_theory * 100:+.1f}%)"
)

# Method 2: From embedded thermocouple (TC3)
q_from_TC3 = ca.heat_flux_from_T_at_depth(
    T_meas["TC3"], depth_TC3, T_gas, T_amb, h_gas, h_amb, thicknesses, conductivities
)
print("\n  From embedded thermocouple:")
print(
    f"    TC3 -> q = {q_from_TC3:.1f} W/m^2 (error: {(q_from_TC3 - q_theory) / q_theory * 100:+.1f}%)"
)

# Average of all measurements
q_avg = np.mean([q_from_TC1, q_from_TC2, q_from_TC3, q_from_TC4])
print(f"\n  Average from all TCs: q = {q_avg:.1f} W/m^2")
print(f"  Theoretical value:    q = {q_theory:.1f} W/m^2")

# ============================================================================
# Sensitivity analysis: Which TC location is most sensitive?
# ============================================================================
print("\n--- Measurement Sensitivity ---")
print("  (Temperature change per 1 W/m^2 change in heat flux)")

# Sensitivity = dT/dq = R (thermal resistance to that point)
R_to_edge = [1 / h_gas]  # Edge 0
for tk in t_over_k:
    R_to_edge.append(R_to_edge[-1] + tk)

for i, (name, R) in enumerate(zip(edge_names, R_to_edge, strict=False)):
    print(f"    Edge {i}: {R * 1000:.2f} mK per W/m^2 ({name})")

print("\n  -> Thermocouples deeper in the wall (higher R) are more sensitive")
print("     to heat flux changes but also more affected by material uncertainty.")

# ============================================================================
# Heat loss calculation
# ============================================================================
print("\n--- Heat Loss Estimate ---")

furnace_area = 50.0  # m^2 (example furnace surface area)
Q_loss = q_avg * furnace_area / 1000  # kW

print(f"  Furnace surface area: {furnace_area:.0f} m^2")
print(f"  Total heat loss: {Q_loss:.1f} kW")
print(f"  Annual energy cost @ $0.10/kWh: ${Q_loss * 8760 * 0.10:.0f}")

print("\n" + "=" * 70)
print("Analysis complete.")
print("=" * 70)
