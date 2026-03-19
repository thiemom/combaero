"""Simple two-pipe wall coupling example.

Demonstrates heat transfer between a hot pipe and a cold pipe through a shared wall.
This is a post-processing example showing how to use ConvectiveSurface and WallConnection.
"""

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, RibbedModel, SmoothModel, WallConnection

# ============================================================================
# Flow conditions
# ============================================================================

# Hot side (combustion products)
T_hot = 1200.0  # K
P_hot = 2.0e5  # Pa
u_hot = 50.0  # m/s
D_hot = 0.05  # m
L_hot = 0.5  # m

# Cold side (cooling air)
T_cold = 400.0  # K
P_cold = 2.5e5  # Pa
u_cold = 60.0  # m/s
D_cold = 0.01  # m
L_cold = 0.5  # m

# Wall properties (Hastelloy X)
t_wall = 0.002  # m
k_wall = 25.0  # W/(m*K)

# Composition (air for both sides)
X_air = cb.standard_dry_air_composition()

# ============================================================================
# Convective surfaces
# ============================================================================

# Hot side: smooth channel
hot_surface = ConvectiveSurface(
    area=3.14159 * D_hot * L_hot,  # pi * D * L
    model=SmoothModel(correlation="gnielinski"),
    heating=False,  # Hot gas is being cooled
)

# Cold side: ribbed channel for enhanced cooling
cold_surface = ConvectiveSurface(
    area=3.14159 * D_cold * L_cold,
    model=RibbedModel(e_D=0.05, pitch_to_height=10.0, alpha_deg=60.0),
    heating=True,  # Coolant is being heated
    Nu_multiplier=1.1,  # 10% enhancement from test data
)

# ============================================================================
# Compute HTCs and adiabatic wall temperatures
# ============================================================================

print("=" * 70)
print("Two-Pipe Wall Coupling Example")
print("=" * 70)
print()

# Hot side
h_hot, T_aw_hot, A_hot = hot_surface.htc_and_T(
    T=T_hot, P=P_hot, X=X_air, velocity=u_hot, diameter=D_hot, length=L_hot
)

print(f"Hot Side (Smooth Channel):")
print(f"  T_bulk = {T_hot:.1f} K")
print(f"  T_aw   = {T_aw_hot:.1f} K")
print(f"  h      = {h_hot:.1f} W/(m^2*K)")
print(f"  A_conv = {A_hot:.4f} m^2")
print()

# Cold side
h_cold, T_aw_cold, A_cold = cold_surface.htc_and_T(
    T=T_cold, P=P_cold, X=X_air, velocity=u_cold, diameter=D_cold, length=L_cold
)

print(f"Cold Side (Ribbed Channel):")
print(f"  T_bulk = {T_cold:.1f} K")
print(f"  T_aw   = {T_aw_cold:.1f} K")
print(f"  h      = {h_cold:.1f} W/(m^2*K)")
print(f"  A_conv = {A_cold:.4f} m^2")
print()

# ============================================================================
# Wall coupling
# ============================================================================

wall = WallConnection(
    id="liner_wall",
    element_a="hot_pipe",
    element_b="cold_pipe",
    wall_thickness=t_wall,
    wall_conductivity=k_wall,
)

Q, T_wall = wall.compute_coupling(h_hot, T_aw_hot, A_hot, h_cold, T_aw_cold, A_cold)

print(f"Wall Coupling:")
print(f"  Wall thickness = {t_wall*1000:.1f} mm")
print(f"  Wall conductivity = {k_wall:.1f} W/(m*K)")
print(f"  Effective area = {min(A_hot, A_cold):.4f} m^2")
print()

# Overall HTC
t_over_k = t_wall / k_wall
U = 1.0 / (1.0 / h_hot + t_over_k + 1.0 / h_cold)

print(f"Results:")
print(f"  Overall HTC (U) = {U:.1f} W/(m^2*K)")
print(f"  Wall temperature = {T_wall:.1f} K")
print(f"  Heat transfer rate = {Q:.1f} W")
print(f"  Heat flux = {Q/min(A_hot, A_cold):.1f} W/m^2")
print()

# Temperature drops
dT_hot_side = T_aw_hot - T_wall
dT_wall = T_wall - T_aw_cold  # Should be small for thin, conductive wall
dT_cold_side = T_wall - T_aw_cold

print(f"Temperature Drops:")
print(f"  Hot side (convection):  {dT_hot_side:.1f} K")
print(f"  Wall (conduction):      {dT_wall:.1f} K")
print(f"  Cold side (convection): {dT_cold_side:.1f} K")
print(f"  Total:                  {T_aw_hot - T_aw_cold:.1f} K")
print()

# Thermal resistances
R_hot = 1.0 / (h_hot * min(A_hot, A_cold))
R_wall = t_over_k / min(A_hot, A_cold)
R_cold = 1.0 / (h_cold * min(A_hot, A_cold))
R_total = R_hot + R_wall + R_cold

print(f"Thermal Resistances:")
print(f"  Hot side:  {R_hot:.4f} K/W ({R_hot/R_total*100:.1f}%)")
print(f"  Wall:      {R_wall:.4f} K/W ({R_wall/R_total*100:.1f}%)")
print(f"  Cold side: {R_cold:.4f} K/W ({R_cold/R_total*100:.1f}%)")
print(f"  Total:     {R_total:.4f} K/W")
print()

print("=" * 70)
print("Note: This is a post-processing example. For full network solver")
print("integration with Newton coupling, see the design document.")
print("=" * 70)
