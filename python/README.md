# CombAero Python bindings

This package provides Python bindings for the core CombAero C++ library:

- Thermodynamic and transport properties for gas mixtures
- Combustion helpers (oxygen requirement, complete combustion)
- Equivalence ratio and Bilger mixture-fraction utilities
- Humid air properties and standard dry-air composition
- Species common-name mapping
- Compressible flow (isentropic nozzle, quasi-1D, Fanno flow)
- Incompressible flow (Bernoulli, orifice, pipe pressure drop)
- Friction factor correlations (Colebrook, Haaland, Serghides)
- Heat transfer correlations (Dittus-Boelter, Gnielinski, Sieder-Tate)
- Integrated channel-flow solvers (`channel_smooth`, `channel_ribbed`, `channel_dimpled`, `channel_pin_fin`, `channel_impingement`)
- Advanced cooling correlations (rib enhancement, impingement, film cooling, effusion, pin fins, dimples)
- Materials database (Inconel 718, Haynes 230, SS316, Al 6061, YSZ TBC)
- Acoustics (duct modes, Helmholtz resonators, Q-factor screening)

The bindings are implemented with [pybind11](https://pybind11.readthedocs.io/) and built via [scikit-build-core](https://scikit-build-core.readthedocs.io/).

## Installation (from source)

The recommended workflow is from the repository root; see the top-level README for details. In short:

```bash
python -m build --wheel
python -m pip install dist/combaero-0.0.1-*.whl
```

## Usage

Example using NumPy arrays and the high-level `combaero` package:

```python
import numpy as np
import combaero as ca

# Get species index for N2 and O2
i_N2 = ca.species_index_from_name("N2")
i_O2 = ca.species_index_from_name("O2")

# Create air composition (79% N2, 21% O2)
X = np.zeros(ca.num_species(), dtype=float)
X[i_N2] = 0.79
X[i_O2] = 0.21

T = 300.0   # K
P = 101325.0  # Pa

print(f"cp = {ca.cp(T, X):.2f} J/(mol·K)")
print(f"density = {ca.density(T, P, X):.4f} kg/m³")
print(f"viscosity = {ca.viscosity(T, P, X):.2e} Pa·s")
```

### State-based API

The `State` class uses Pythonic properties for all attribute access:

```python
import numpy as np
import combaero as ca

# Create a State - direct property assignment
s = ca.State()
s.T = 300.0
s.P = 101325.0
s.X = X_air

# Or use fluent setters for chaining
s = ca.State()
s.set_T(300.0).set_P(101325.0).set_X(X_air)

# Access thermodynamic properties (all are properties, not methods)
print(f"cp = {s.cp:.2f} J/(mol·K)")
print(f"h = {s.h:.0f} J/mol")
print(f"rho = {s.rho:.4f} kg/m³")
print(f"gamma = {s.gamma:.3f}")
print(f"speed of sound = {s.a:.1f} m/s")

# Access transport properties
print(f"viscosity = {s.mu:.2e} Pa·s")
print(f"thermal conductivity = {s.k:.4f} W/(m·K)")
print(f"Prandtl = {s.Pr:.3f}")
```

### Combustion and Equilibrium

```python
import numpy as np
import combaero as ca

# Create a CH4 + air mixture
X = np.zeros(ca.num_species(), dtype=float)
X[ca.species_index_from_name("CH4")] = 0.095
X[ca.species_index_from_name("O2")] = 0.19
X[ca.species_index_from_name("N2")] = 0.715

# Adiabatic complete combustion
burned = ca.complete_combustion(T=300.0, X=X)
print(f"Adiabatic flame T: {burned.T:.0f} K")
print(f"Flame density: {burned.rho:.4f} kg/m³")

# WGS equilibrium (isothermal or adiabatic)
eq = ca.wgs_equilibrium_adiabatic(T=1500.0, X=burned.X)
print(f"Equilibrium T: {eq.T:.0f} K")
```

### Reforming Equilibrium (Rich Combustion)

For rich combustion products with unburned hydrocarbons, use reforming equilibrium
to compute the steam reforming + water-gas shift equilibrium:

```python
import numpy as np
import combaero as ca

# Rich combustion products (with unburned CH4, C2H6, etc.)
X = np.zeros(ca.num_species(), dtype=float)
X[ca.species_index_from_name("CH4")] = 0.02    # Unburned methane
X[ca.species_index_from_name("C2H6")] = 0.003  # Unburned ethane
X[ca.species_index_from_name("H2O")] = 0.20
X[ca.species_index_from_name("CO2")] = 0.08
X[ca.species_index_from_name("N2")] = 0.70

# General reforming equilibrium (handles all hydrocarbons)
# CnHm + n*H2O <-> n*CO + (n + m/2)*H2
# CO + H2O <-> CO2 + H2 (WGS)
eq = ca.reforming_equilibrium_adiabatic(T=2200.0, X=X)
print(f"Equilibrium T: {eq.T:.0f} K")
print(f"CO: {eq.X[ca.species_index_from_name('CO')]*100:.2f}%")
print(f"H2: {eq.X[ca.species_index_from_name('H2')]*100:.2f}%")

# Isothermal version (fixed temperature)
eq_iso = ca.reforming_equilibrium(T=2000.0, X=X)

# SMR+WGS (CH4 only, for backward compatibility)
eq_smr = ca.smr_wgs_equilibrium_adiabatic(T=2200.0, X=X)
```

**Available equilibrium functions:**
- `wgs_equilibrium(T, X, P)` - Isothermal WGS (CO + H2O ⇌ CO2 + H2)
- `wgs_equilibrium_adiabatic(T, X, P)` - Adiabatic WGS
- `smr_wgs_equilibrium(T, X, P)` - Isothermal SMR+WGS (CH4 only)
- `smr_wgs_equilibrium_adiabatic(T, X, P)` - Adiabatic SMR+WGS (CH4 only)
- `reforming_equilibrium(T, X, P)` - Isothermal reforming (all hydrocarbons)
- `reforming_equilibrium_adiabatic(T, X, P)` - Adiabatic reforming (all hydrocarbons)
- `combustion_equilibrium(T, X, P)` - **One-step**: combustion + reforming + WGS

**Tip:** Use `combustion_equilibrium` when starting from an unburned fuel+air mixture.
It combines complete combustion with reforming equilibrium in one call:

```python
# Unburned fuel + air -> equilibrium products in one step
X_unburned = ...  # CH4 + O2 + N2
result = ca.combustion_equilibrium(T=300.0, X=X_unburned)
print(f"Equilibrium T: {result.T:.0f} K")
```

### Stream Mixing

Mix multiple streams with mass and enthalpy balance:

```python
import numpy as np
import combaero as ca

# Create streams - use property assignment (Pythonic)
air = ca.Stream()
air.T = 400.0
air.P = 101325.0
air.X = X_air
air.mdot = 10.0

fuel = ca.Stream()
fuel.T = 300.0
fuel.P = 101325.0
fuel.X = X_fuel
fuel.mdot = 0.5

# Or use fluent setters for one-liners
fuel.set_T(300.0).set_P(101325.0).set_X(X_fuel).set_mdot(0.5)

# Mix streams (uses minimum inlet pressure by default)
mixed = ca.mix([air, fuel])
print(f"Mixed T: {mixed.T:.1f} K")
print(f"Mixed mdot: {mixed.mdot:.2f} kg/s")

# Or specify output pressure explicitly
mixed2 = ca.mix([air, fuel], P_out=150000.0)
```

### Compressible Flow

```python
import combaero as ca

X = ca.standard_dry_air_composition()

# Isentropic nozzle flow
T0, P0 = 500.0, 500000.0  # Stagnation conditions
P_back = 300000.0  # Back pressure
A_eff = 0.001  # Effective area [m²]

sol = ca.nozzle_flow(T0, P0, P_back, A_eff, X)
print(f"Mach: {sol.M:.3f}, mdot: {sol.mdot:.4f} kg/s, choked: {sol.choked}")

# Critical pressure ratio
P_ratio = ca.critical_pressure_ratio(T0, P0, X)
print(f"P*/P0 = {P_ratio:.4f}")

# Converging-diverging nozzle with axial profile
sol_cd = ca.nozzle_cd(
    T0=600.0, P0=1e6, P_exit=2e5,
    A_inlet=0.01, A_throat=0.005, A_exit=0.008,
    x_throat=0.1, x_exit=0.2, X=X
)
print(f"Mass flow: {sol_cd.mdot:.4f} kg/s")
for st in sol_cd.profile[::10]:  # Every 10th station
    print(f"x={st.x:.3f} m, M={st.M:.3f}, T={st.T:.1f} K")

# Fanno flow (adiabatic pipe with friction)
f = ca.friction_colebrook(Re=100000, e_D=0.0001)
sol_fanno = ca.fanno_pipe(T_in=400.0, P_in=5e5, u_in=100.0,
                          L=5.0, D=0.05, f=f, X=X)
print(f"Outlet P: {sol_fanno.outlet.P/1000:.1f} kPa")
```

### Transport Properties

```python
import combaero as ca

X = ca.standard_dry_air_composition()
T, P = 300.0, 101325.0

# Basic transport
print(f"Viscosity: {ca.viscosity(T, P, X):.2e} Pa·s")
print(f"Thermal conductivity: {ca.thermal_conductivity(T, P, X):.4f} W/(m·K)")
print(f"Prandtl: {ca.prandtl(T, P, X):.3f}")

# Dimensionless numbers
V, L = 10.0, 0.1  # Velocity [m/s], length scale [m]
Re = ca.reynolds(T, P, X, V, L)
Pe = ca.peclet(T, P, X, V, L)
print(f"Reynolds: {Re:.0f}, Peclet: {Pe:.0f}")
```

### Species Common Names

```python
import combaero as ca

# Formula to common name
print(ca.common_name("CH4"))  # "Methane"
print(ca.common_name("N2"))   # "Nitrogen"

# Common name to formula
print(ca.formula("Methane"))  # "CH4"
print(ca.formula("Water"))    # "H2O"

# Get full mappings
formula_map = ca.formula_to_name()  # dict: formula -> name
name_map = ca.name_to_formula()     # dict: name -> formula
```

### Incompressible Flow

For liquids and low-speed gas flows (Ma < 0.3):

```python
import combaero as ca

rho = 998.0  # Water density [kg/m3]

# Bernoulli equation
P2 = ca.bernoulli_P2(P1=200000, v1=2.0, v2=5.0, rho=rho)
v2 = ca.bernoulli_v2(P1=200000, P2=150000, v1=2.0, rho=rho)

# Orifice flow
mdot = ca.orifice_mdot(P1=200000, P2=100000, A=0.001, Cd=0.62, rho=rho)
A = ca.orifice_area(mdot=5.0, P1=200000, P2=100000, Cd=0.62, rho=rho)

# Pipe pressure drop (Darcy-Weisbach)
f = ca.friction_colebrook(Re=100000, e_D=0.0001)
dP = ca.pipe_dP(v=2.0, L=10.0, D=0.05, f=f, rho=rho)

# Hydraulic diameter
Dh = ca.hydraulic_diameter_rect(a=0.1, b=0.05)  # Rectangular duct
```

### Channel Flow (Heat Transfer + Pressure Drop)

Integrated solvers that return `ChannelResult` with h, Nu, Re, Pr, f, dP, Mach, T_aw, and q in one call:

```python
import combaero as ca

X = ca.standard_dry_air_composition()
T, P, u = 600.0, 20e5, 30.0  # K, Pa, m/s
D, L = 0.008, 0.35           # m

# Smooth pipe (Gnielinski baseline)
sol = ca.channel_smooth(T, P, X, u, D, L)
print(f"h = {sol.h:.0f} W/(m2*K), Nu = {sol.Nu:.1f}, dP = {sol.dP:.0f} Pa")

# Rib-roughened (Han et al. 1988)
sol_rib = ca.channel_ribbed(T, P, X, u, D, L,
                             e_D=0.07, pitch_to_height=8.0, alpha_deg=60.0)
print(f"Rib enhancement: Nu_rib/Nu_smooth = {sol_rib.Nu/sol.Nu:.2f}")

# Dimpled surface (Chyu et al. 1997)
sol_dim = ca.channel_dimpled(T, P, X, u, D, L, d_Dh=0.2, h_d=0.2, S_d=2.0)

# Pin-fin array (Metzger et al. 1982)
sol_pin = ca.channel_pin_fin(T, P, X, u,
                              channel_height=0.01, pin_diameter=0.003,
                              S_D=2.5, X_D=2.5, N_rows=4)

# With wall temperature (enables q calculation)
sol_q = ca.channel_smooth(T, P, X, u, D, L, T_wall=800.0)
print(f"q = {sol_q.q/1000:.1f} kW/m2")
```

**`ChannelResult` attributes:** `h`, `Nu`, `Re`, `Pr`, `f`, `dP`, `M`, `T_aw`, `q`

### Advanced Cooling Correlations

```python
import combaero as ca

# Rib enhancement factor - geometry-only (Han et al. 1988, used by channel_ribbed)
enh = ca.rib_enhancement_factor(e_D=0.05, pitch_to_height=10.0, alpha_deg=60.0)
fmul = ca.rib_friction_multiplier(e_D=0.05, pitch_to_height=10.0)

# High-Re direct empirical fit (Singh & Ekkad 2017, Re = 30k-400k)
# Nu/Nu0 = 54.9 * Re^(-0.1755) * (e/D)^0.38 * (P/e)^(-0.11) * F(alpha)
enh_hre = ca.rib_enhancement_factor_high_re(e_D=0.0625, pitch_to_height=10.0,
                                             alpha_deg=60.0, Re=200_000.0)
fmul_hre = ca.rib_friction_multiplier_high_re(e_D=0.0625, pitch_to_height=10.0)

# Thermal-hydraulic performance factor (Webb & Eckert 1972)
# eta = (Nu/Nu0) / (f/f0)^(1/3)  -- eta > 1: net improvement; eta < 1: pressure-drop generator
eta = ca.thermal_performance_factor(Nu_ratio=enh_hre, f_ratio=fmul_hre)

# Impingement Nusselt (Martin 1977 / Florschuetz 1981)
Nu = ca.impingement_nusselt(Re_jet=20000, Pr=0.7, z_D=6.0)
Nu_arr = ca.impingement_nusselt(Re_jet=20000, Pr=0.7, z_D=6.0, x_D=8.0, y_D=8.0)

# Film cooling effectiveness (Baldauf et al. 2002)
eta = ca.film_cooling_effectiveness(x_D=10.0, M=1.0, DR=1.5, alpha_deg=30.0)
eta_avg = ca.film_cooling_effectiveness_avg(x_D=10.0, M=1.0, DR=1.5,
                                             alpha_deg=30.0, s_D=3.0)

# Effusion effectiveness
eta_eff = ca.effusion_effectiveness(x_D=5.0, M=0.5, DR=1.5,
                                    porosity=0.2, s_D=3.0, alpha_deg=30.0)
```

### Materials Database

```python
import combaero as ca

# Superalloys
k = ca.k_inconel718(T=900)      # W/(m*K), valid 300-1200 K
k = ca.k_haynes230(T=1100)      # W/(m*K), valid 300-1400 K

# Structural alloys
k = ca.k_stainless_steel_316(T=600)  # W/(m*K)
k = ca.k_aluminum_6061(T=400)        # W/(m*K)

# Thermal barrier coating (YSZ) with sintering model
k_fresh = ca.k_tbc_ysz(T=1200, hours=0)       # As-sprayed
k_aged  = ca.k_tbc_ysz(T=1200, hours=1000)    # After 1000 h
k_ebpvd = ca.k_tbc_ysz(T=1200, is_ebpvd=True) # EB-PVD process

# List all available materials
print(ca.list_materials())
```
