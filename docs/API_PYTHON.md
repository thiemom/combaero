# CombAero Python API Reference

This document provides the high-level reference for the `combaero` Python package. For the C++ technical reference, see [API_CPP.md](API_CPP.md).

## Table of Contents
- [Core Thermodynamics](#core-thermodynamics)
- [Combustion & Equilibrium](#combustion--equilibrium)
- [Flow Regimes (Symmetric API)](#flow-regimes-symmetric-api)
- [Advanced Flow Functions](#advanced-flow-functions)
- [Comprehensive Orifice Functions](#comprehensive-orifice-functions)
- [Heat Transfer](#heat-transfer)
  - [Correlations](#correlations)
  - [Network Heat Transfer](#network-heat-transfer)
- [Acoustics](#acoustics)
- [Network Solver](#network-solver)
- [NetworkRunner (GUI-JSON Programmatic Driver)](#networkrunner-gui-json-programmatic-driver)
- [Geometry & Materials](#geometry--materials)
- [Advanced Thermodynamics](#advanced-thermodynamics)
- [Psychrometrics (Humid Air)](#psychrometrics-humid-air)

---

## Core Thermodynamics

### State Object
The `State` object is the central data structure, managing T [K], P [Pa], and composition.

```python
import combaero as cb

state = cb.State()
state.TPX = (300.0, 101325.0, "O2:0.21, N2:0.79")
# OR setting mass fractions directly
state.Y = y_vector

print(f"Mole fractions: {state.X}")
print(f"Mass fractions: {state.Y}")

print(f"Density: {state.rho} kg/m³")
print(f"Enthalpy: {state.h} J/mol")
print(f"Viscosity: {state.mu} Pa·s")
```

### MixtureState (Network Data)
The `MixtureState` struct is used for network node properties and carries both static and stagnation conditions.
**Units**: P [Pa], T [K], m_dot [kg/s].

```python
# Constructor: MixtureState(P, P_total, T, T_total, m_dot, Y)
node_state = cb.MixtureState(1e5, 1.1e5, 300, 300, 1.5, y_vec)

print(node_state.P)       # Static [Pa]
print(node_state.m_dot)   # Mass flow [kg/s]
print(node_state.density()) # Static density [kg/m³]
```

### Species Module
The `combaero.species` module provides high-level helpers for composition vectors.

```python
import combaero as cb

# Get standard compositions (returns numpy array)
X_air = cb.species.dry_air()        # Mole fractions
Y_air = cb.mole_to_mass(X_air)   # Mass fractions (standard for network)

# Humid air at T [K], P [Pa], RH [0-1]
Y_humid = cb.species.humid_air_mass(300, 101325, 0.5)

# Pure species or from mapping
Y_ch4 = cb.species.pure_species("CH4")
Y_mix = cb.species.from_mapping({"O2": 0.20, "N2": 78, "AR": 0.02})

# Explicit conversion
Y = cb.species.to_mass(X_air)
X = cb.species.to_mole(Y_air)
```

### Species Lookup (Legacy)
```python
# Low-level metadata from core
num = cb.num_species()
name = cb.species_name(0)
idx = cb.species_index_from_name("H2O")
```

---

## Combustion & Equilibrium

### Complete Combustion
```python
# ADIABATIC flame temperature and products
# smooth=True ensures a continuous Jacobian for root-finding near Phi=0
burned = cb.complete_combustion(state, smooth=True)
print(f"Adiabatic T: {burned.T} K")

# Isothermal products
products = cb.complete_combustion_isothermal(state, smooth=True)
```

### Chemical Equilibrium
```python
# T, P constant equilibrium
eq_state = cb.combustion_equilibrium(state, smooth=True)

# Water-gas shift
wgs = cb.wgs_equilibrium(state)
```

---

## Flow Regimes (Symmetric API)

CombAero provides symmetric submodules for incompressible and compressible flow.

```python
from combaero import incompressible as flow   # OR: from combaero import compressible as flow

sol = flow.channel_flow(T=400, P=2e5, X=air, u=10.0, L=2.0, D=0.05, f=0.02)
print(f"Pressure drop: {sol.dP} Pa")
print(f"Outlet Mach: {sol.M}") # nan for incompressible
```

### FlowSolution Fields
- `mdot`: Mass flow rate [kg/s]
- `v`: bulk velocity [m/s]
- `dP`: Pressure drop [Pa]
- `M`: Mach number
- `choked`: Boolean flag

---

## Advanced Flow Functions

### Inverse Solvers
```python
# Solve for unknown given mass flow rate
A_eff = cb.solve_A_eff_from_mdot(T0=300, P0=2e5, P_back=1e5, mdot_target=0.1, X=air)
P_back = cb.solve_P_back_from_mdot(T0=300, P0=2e5, A_eff=1e-4, mdot_target=0.1, X=air)
P0 = cb.solve_P0_from_mdot(T0=300, P_back=1e5, A_eff=1e-4, mdot_target=0.1, X=air)
```

### Utility Functions
```python
# Critical pressure ratio and Mach number
P_crit_ratio = cb.critical_pressure_ratio(T0=300, P0=2e5, X=air)
Mach = cb.mach_from_pressure_ratio(T0=300, P0=2e5, P=1.5e5, X=air)

# Mass flux and flow analysis
G = cb.mass_flux_isentropic(T0=300, P0=2e5, P=1.5e5, X=air)
L_max = cb.fanno_max_length(T_in=300, P_in=2e5, u_in=50, D=0.05, f=0.02, X=air)
```

### Fanno flow in Mach number
Integrated in M rather than x: the state is algebraic in M at fixed mass flux,
so the duct length is a smooth integral that vanishes at sonic -- exact L*,
no cutoff. Fed isentropically from a stagnation state:

```python
G_ch = cb.fanno_choked_mass_flux(Pt=2e5, Tt=300, X=air, L=1.0, D=0.02,
                                 roughness=1e-5, friction_model="haaland")  # kg/(m^2 s)
r = cb.fanno_duct(2e5, 300, 0.8 * G_ch, air, 1.0, 0.02, 1e-5, "haaland")
r.inlet.M, r.exit.M, r.exit.P, r.Pt_exit, r.L_star, r.choked
cb.fanno_length_between(G=300, Tt=300, M1=0.3, M2=1.0, X=air, D=0.02,
                        roughness=0.0, friction_model="fixed", f_multiplier=0.02)
```
Agrees with the x-march (`fanno_channel_rough`, 20000 steps) to 1e-9 at
constant f and to ~1e-8 with a Re-dependent f, the floor set by kinks in
mu(T).

### Rocket Nozzle Functions
```python
# Converging-diverging nozzle analysis
sol = cb.nozzle_cd(T0=300, P0=2e5, P_design=1e5, P_amb=101325,
                   A_inlet=1e-3, A_throat=5e-4, A_exit=2e-3,
                   x_throat=0.1, x_exit=0.2, X=air)

# Thrust calculations
thrust = cb.nozzle_thrust_cd(T0=300, P0=2e5, P_design=1e5, P_amb=101325,
                             A_inlet=1e-3, A_throat=5e-4, A_exit=2e-3,
                             x_throat=0.1, x_exit=0.2, X=air)
print(f"Thrust: {thrust.thrust} N")
print(f"Specific impulse: {thrust.specific_impulse} s")
```

---

## Comprehensive Orifice Functions

### Discharge Coefficient Correlations
```python
# Normed metering orifice (ISO 5167 family, a function of beta = d/D)
Cd = cb.Cd_sharp_thin_plate(geom, state)

# Discharge hole in a wall -- no pipe, no beta. Separate selector, separate
# geometry. Idelchik (1966) diagrams 4-17/4-18, valid Re 25 to 1e6.
hole = cb.DischargeHoleGeometry(d=1e-3, L=2e-3, r=0.0)
flow = cb.DischargeHoleState(Re=5e4)
Cd = cb.discharge_cd(cb.DischargeCdCorrelation.Idelchik1966Thick, hole, flow)

# Solver-facing form: (Cd, dCd/dRe, dCd/d(U1_over_Vi)), analytic throughout
Cd, dCd_dRe, dCd_dU = cb.discharge_cd_and_derivatives(
    cb.DischargeCdCorrelation.Idelchik1966Thick, hole, flow
)

# Inside the source's own Re and l/d range? A flag, never a selector.
cb.discharge_cd_in_range(cb.DischargeCdCorrelation.Lichtarowicz1965, hole, flow)
```

There is deliberately no auto-selection from geometry: it is how a
rounded-entry request used to come back as Stolz.

#### Plenum-to-plenum holes: McGreehan and Schotsch (1988)

The correlations above are ISO 5167 pipe metering -- a single hole IN a pipe
run. For a cooling transfer hole or a jet-plate hole discharging
plenum-to-plenum, use the composite chain of McGreehan and Schotsch (1988),
ASME J. Turbomachinery 110(2), 213-217, which covers Reynolds number, inlet
corner radius, orifice length and inlet crossflow:

```python
# A bare, sharp-edged jet plate, plate thickness = hole diameter
cb.mcgreehan_schotsch_1988_cd(Re=1e4, r_over_d=0.0, L_over_d=1.0)   # 0.785

# A radiused, long cooling hole
cb.mcgreehan_schotsch_1988_cd(Re=3e4, r_over_d=0.1, L_over_d=3.0)   # 0.889

# Fed from a duct rather than a plenum: supply-side crossflow
cb.mcgreehan_schotsch_1988_cd(Re=3e4, r_over_d=0.0, L_over_d=1.0,
                              U1_over_Vi=0.5)                       # 0.770
```

`U1_over_Vi` is the **inlet (approach, supply-side)** tangential velocity over
the ideal through-flow velocity, and defaults to 0 for a plenum-fed plate. It
is *not* a discharge-side crossflow ratio -- do not pass Florschuetz's
`Gc/Gj`, which is spent air in the impingement channel on the far face of the
plate. `V_i` must be built from static inlet conditions, not from a total
pressure that already contains the tangential velocity head.

Two behaviours that look like bugs and are not:

- **Not monotonic in `U1_over_Vi`.** `Cd` rises up to ~5.7% above its
  zero-crossflow value near `U1_over_Vi ~ 0.09` before falling away. This is
  drawn in the source's Fig. 4 and supported by Rohde (NASA TN D-5467), whose
  result 2 is that slanting an orifice into the flow increases `Cd`.
- **Constant below `Re = 1e4`.** That is the correlation's stated floor, and
  below it the underlying equation diverges (it reaches `Cd = 1.0` at
  `Re = 904`), so the value is held at the floor rather than extrapolated.

When the plate's zero-crossflow `Cd` is already known -- measured, or a
literature value such as Florschuetz's per-configuration Table 1 -- apply only
the crossflow correction:

```python
cb.mcgreehan_schotsch_1988_crossflow_cd(cd_base=0.79, U1_over_Vi=0.5)  # 0.794
cb.mcgreehan_schotsch_1988_crossflow_cd(cd_base=0.79, U1_over_Vi=2.0)  # 0.541
```

Note the first: at a 0.79 baseline, `U1/Vi = 0.5` still sits inside the rise,
because a higher baseline pushes the peak to larger `U1/Vi` (`R_v` carries a
`(Cd/0.6)^-3` factor). The decay is well established by `U1/Vi = 2`.

This is how the source uses Eq. (17) in its own validation. Its Figs. 5 and 6
anchor to "a set baseline point at `U1/Vi = 0`" taken from Rohde's
measurements (0.64, 0.73, 0.88), not to the chain's prediction for the same
geometry -- the chain runs 3-12% higher, which is the paper's own remark that
Rohde's "basic values are lower".

**Measured agreement.** Against Rohde's data as replotted in Fig. 6
(`t/d = 0.51`, baseline 0.64, 11 points), Eq. (17) scores bias `+5.95%`,
RMS `7.54%` -- reading high, and increasingly so with crossflow (`+3%` below
`U1/Vi = 0.4`, `+16%` at 1.41). Scored by
`validation/cooling/orifice_runner.py`.

#### Film cooling: Baldauf et al. (2002)

Laterally averaged adiabatic effectiveness downstream of one row of
cylindrical, streamwise-inclined holes -- valid from the ejection point, and
carrying the adjacent jet interaction rather than excluding it.

```python
eta = cb.film_effectiveness_baldauf_2002(
    x_over_D=20.0, M=2.0, P=1.2, alpha_deg=30.0, s_over_D=3.0, Tu=0.015
)

# Solver form: (eta, deta/dM, deta/dP), analytic
eta, deta_dM, deta_dP = cb.film_effectiveness_baldauf_2002_and_derivatives(
    20.0, 2.0, 1.2, 30.0, 3.0, 0.015
)
```

`eta` uses the same convention as `adiabatic_wall_temperature`, so they
compose. `deta/dM` changes sign through the lift-off peak, which is physical
and is what a solver needs to push blowing the right way.

Envelope: `M` 0.2-2.5, `P` 1.2-1.8, `s/D` 2-5, `alpha` 30-90 deg,
`Tu` 0.0035-0.075; the paper's RMS deviation is 5.5%. Eq. (31) is
implemented as printed and disagrees with the paper's own worked example --
see `validation/cooling/extractions/baldauf_2002_film_effectiveness.md`.

`b_0=<float>` replaces Eq. (31)'s `b_0` (default `None` = as printed), so the
paper's Table 4 reading can be tested. Its effect is governed by `x/D`, not
`M`: under 1% for `x/D <= 20` at any blowing, -28% by `x/D = 400`.

Outside the envelope the value is unchanged -- it still answers, finitely --
but each out-of-range input now raises a warning through the global handler,
one per parameter:

```python
with cb.suppress_warnings():          # or set_warning_handler to collect them
    eta = cb.film_effectiveness_baldauf_2002(20.0, 3.0, 1.0, 30.0, 7.37, 0.05)
# unsuppressed, that call reports M, P and s_over_D as extrapolated
```

#### Effusion plate internal heat transfer

```python
nu = cb.effusion_internal_nusselt(Re=3000.0, Pr=0.727,
                                  X_over_L=15.24 / 6.3, L_over_D=6.3 / 3.27)

# the two terms separately, if you want to see which dominates
nu_a = cb.effusion_approach_nusselt(3000.0, 0.727, 15.24 / 6.3)
nu_t = cb.effusion_throat_nusselt(3000.0, 0.727, 6.3 / 3.27)
r    = cb.mills_entry_length_factor(6.3 / 3.27)     # entry-length factor
```

`Re` is on the hole diameter and `Nu` on the **hole internal area**. For a
plate-area coefficient multiply by `k/D` and then by `A_h/A` with
`A_h = pi D L` and `A = X**2 - pi D**2/4` -- a factor of 3.46 for Andrews'
plate C, so not optional.

The approach term leads at low `Re` (exponent 0.476) and the throat
overtakes it (0.8); the crossover falls inside a typical effusion plate's
own operating range, which is why both are needed.

##### On the element

`EffusionPlateElement.internal_heat_transfer(state_in)` applies the above to
the panel's own geometry and appears in `diagnostics()` after a solve:

```python
diag = result["__element_diag__"]["panel"]
diag["h_plate_area"]   # W/m^2 K on the approach area, Andrews' definition
diag["h_hole_area"]    # W/m^2 K on the hole internal area
diag["area_ratio"]     # approach / hole internal, 3.46 for Andrews' plate C
diag["Nu_approach"], diag["Nu_throat"]
```

**Both coefficients are reported because they differ by `area_ratio`**, and
using one where the other is meant is a factor-of-three error. For a panel
energy balance: `Q = h_plate_area * area_approach_total * dT`.

The external film is **not** included here, and `overall_effectiveness()`
below explains why that is a measured decision rather than a gap.

##### The plate's own wall in a network (#471)

The plate IS the wall the coolant passes through, so it owns it: no
`ThermalWall` connects to it, and its diagnostics carry the wall solution
(`wall_heat_transfer`). The node it DISCHARGES into decides the gas side:

| Discharge into | Gas side |
|---|---|
| `MomentumChamberNode` (flow) | the chamber's own surface correlation is the unblown `h_gas_unblown` (so the chamber's surface model IS the baseline); its velocity gives `blowing_ratio`; `gas_augmentation` (the caller's, 1.0) scales h. `T_gas` is the approaching MAIN stream when the chamber has one -- its own state is the mixed outlet, diluted by the effused coolant |
| a plenum (state, no flow) | the imposed `gas_heat_flux` [W/m^2], 0 by default (adiabatic) |

```python
panel = EffusionPlateElement("eff", "coolant", "liner", hole_diameter=3.27e-3,
                             wall_thickness=6.3e-3, pitch=15.24e-3,
                             panel_length=0.152, panel_width=0.152,
                             wall_conductivity=20.0,      # one layer, t/k
                             gas_augmentation=1.0,        # the caller's
                             gas_film="none",             # or "baldauf_sellers"
                             gas_heat_flux=0.0)           # plenum discharge only
diag["T_wall_hot"], diag["T_wall_cold"], diag["q_wall"], diag["eta_overall"]
```

- **Energy-neutral for the network.** The wall heat leaves the gas and
  returns with the effusing coolant into the same node, so it is an output,
  not a source.
- **Film off by default.** `gas_film="baldauf_sellers"` adds C++'s
  `effusion_panel_film_effectiveness` (Baldauf per row, Sellers over the
  panel's rows; needs `panel_length`), flagged `film_extrapolated` outside
  Baldauf's envelope. Scored on Andrews 88-GT-290 the data refuse it offered
  alone; the missing physics is the gas-side augmentation.
- **Plenum-fed.** McGreehan-Schotsch's U1/Vi stays 0 and the coolant side is
  Andrews' still plenum; `coolant_crossflow_ignored` flags a supply node a
  channel runs through. A channel-fed liner segment is a separate element.
- With `wall_conductivity` large the wall reduces exactly to
  `overall_effectiveness` below.

##### Duct-fed: the effusion liner (#471)

A liner's coolant usually runs along a backside duct and bleeds through the
wall, so the holes see a SUPPLY-side crossflow. `EffusionLiner` builds that
configuration from N stations:

```python
from combaero.network import EffusionLiner

liner = EffusionLiner("ln", n_segments=4, length=0.152, width=0.152, duct_height=0.03,
                      hole_diameter=3.27e-3, wall_thickness=6.3e-3,
                      pitch_x=15.24e-3, pitch_y=15.24e-3)
liner.add_to(net, coolant_in="cin", coolant_out="cout", gas="liner_chamber")
liner.summarize(res["__element_diag__"])   # bleed, U1/Vi, Cd, P_backside, T_wall per station
```

- **Stations are Bassett bleeds** (`CrossflowSegmentElement`, kappa 0.75):
  the backside static pressure rises along the duct as it slows, so the
  downstream panels see more drive.
- **Each panel is duct-fed**: `EffusionPlateElement(crossflow_segments=(prev,
  next), crossflow_area=A)` takes McGreehan-Schotsch's supply-side
  `U1/Vi`, with `U1` the duct's mean velocity at the station and `Vi` the
  isentropic jet velocity from the station's STATIC state (C++
  `crossflow_velocity_ratio`). Only `McGreehanSchotsch` has a crossflow term,
  so a duct-fed panel refuses the others.
- **What sets U1/Vi** is mainly the duct's own pressure drop against the
  holes' drive, roughly `sqrt(dp_duct/dp_hole)`, not the duct size: a
  pressure-driven duct carries more flow when made bigger. Flags:
  `crossflow_cd_beyond_8pct` (> 0.2) and `crossflow_cd_degraded` (> 0.35),
  the Rohde-scored bounds.
- **Flagged, provisional:** the coolant-side heat transfer is still Andrews'
  plenum-fed correlation (`coolant_ht_plenum_assumed`), and the stations
  bleed through holes far smaller than the duct, beyond Bassett's measured
  psi = 1-3 (`station_psi_beyond_bassett`).

##### Overall effectiveness, as an output

```python
out = panel.overall_effectiveness(state_in, h_gas_unblown=250.0,
                                  T_gas=1600.0, U_gas=80.0,
                                  gas_augmentation=1.0, eta_film=0.0)
out["eta_overall"]       # (Tg - Tw)/(Tg - Tc)
out["T_wall"]            # K
out["T_adiabatic_wall"]  # K, what eta_film sets
out["resistance_ratio"]  # h_internal / h_gas
out["velocity_ratio"]    # jet / mainstream -- reported, never applied
out["blowing_ratio"], out["momentum_flux_ratio"]
```

```
h_gas = h_gas_unblown * gas_augmentation
eta   = (h_i + h_gas * eta_film) / (h_i + h_gas)
```

An overall-effectiveness **correlation** must never be the closure here
-- it already contains the internal convection this computes, and the two
would double-count. So `eta` is an output and the gas side is the
caller's input.

**A film is TWO numbers, and the signature says so.** The standard form
is `q = h_f (T_aw - T_w)` with `T_aw = Tg - eta_f (Tg - Tc)`, so a film
owes both a driving temperature (`eta_film`) and a conductance
(`gas_augmentation = h_f/h_0`). `film_effectiveness_baldauf_2002` gives
only the first -- an adiabatic wall passes no heat, so it measures no
coefficient. **Passing `eta_film` while leaving `gas_augmentation` at 1.0
over-predicts**; both default to their no-film values so that omitting
the pair is consistent rather than half-right.

**Mind what `h_gas_unblown` is measured against.** An augmentation ratio
is meaningless without its baseline. For Andrews 88-GT-290's rig the
required `gas_augmentation` is 3.1-4.3 against a fully developed
Dittus-Boelter and 1.6-2.5 against the same duct with a thermal-entry
correction -- a factor of 1.75 from that choice alone. Published
film-cooling ratios are usually referenced to a flat-plate turbulent
boundary layer at the same x, a third baseline again.

**`gas_augmentation` and `internal_Nu_multiplier` are yours, never
fitted here**, matching `Nu_multiplier` on `ConvectiveSurface`. At their
defaults the closure scores +3.1% on Andrews' effusion plate C and +23.7%
on plate B; setting `gas_augmentation = 2.06` closes plate B at G = 0.6
exactly. Note that `internal_Nu_multiplier = 0.486` also closes it -- by
halving the coolant side, which is the wrong direction, since that
correlation runs 10.4% LOW against Andrews' Fig. 8. Reaching the target
is not the same as being the right dial.

**How good `eta` is depends on a regime this cannot predict.** Andrews'
two plates need gas-side coefficients differing by 1.7x, and that ratio
survives any choice of film model. Pass `U_gas` and read the reported
ratios; if `velocity_ratio` exceeds about 1, treat `eta` as an upper
bound unless `gas_augmentation` already accounts for the jets. No
threshold is applied, because the data supports none -- see
`validation/cooling/extractions/andrews_effusion_overall_eta.md`.

#### Multi-row film superposition

```python
per_row = [cb.film_effectiveness_baldauf_2002(x, M, P, 30.0, 3.0, 0.015)
           for x in row_distances]

eta = cb.film_superposition_sellers(per_row)                  # Gao Eq. (1)
eta = cb.film_superposition_corrected(per_row, alphas)        # Gao Eq. (7)
```

`alphas` has one fewer entry than `per_row`: `alpha_j` is the fraction of the
film's temperature deficit surviving between row `j` and `j+1`. All ones is
plain Sellers.

Sellers **overestimates, and worsens as rows accumulate** -- which is why
`alpha` exists and why effusion cannot reuse a few-row film model unchanged.
`cb.mainstream_temperature_correction(r, a, b)` gives Gao's published form
for it. `a` and `b` are required rather than defaulted even though the paper *does*
print them (12 and 0.9465, section 4.3.1). Since `a r/(a r + 1) >= 0`, that
pair gives `alpha >= 0.9465` for **every** `r`, so it cannot reach the
0.69-0.85 per-row damping a tighter-pitched plate needs -- an argument that
holds however `r` is scaled. Eq. (5) carries no streamwise-spacing term,
and that is the variable Murray shows to dominate.

##### Compressibility: the expansion factor

The paper's mass flow is `W = Cd * Y * A * sqrt(2 rho_t1 (P_t1 - P_s2))`. `Y`
carries compressibility in an otherwise incompressible orifice equation:

```python
S = P_s2 / P_t1
cb.mcgreehan_schotsch_1988_expansion_factor(cd=0.80, S=0.7, gamma=1.4)  # 0.912
```

**Only for the incompressible formulation.** `regime='compressible'` already
solves the isentropic nozzle exactly (`nozzle_flow`, choked branch included),
and Eq. (5) reproduces that solve to 0.008% -- applying `Y` on top of it
corrects for compressibility twice, worth over 8% at `S = 0.7`.

The paper carries two forms because an orifice is not a nozzle, and blends
between them on `Cd`:

```python
cb.mcgreehan_schotsch_1988_expansion_orifice(S=0.7, gamma=1.4)  # 0.912, Eq. (4)
cb.mcgreehan_schotsch_1988_expansion_nozzle(S=0.7, gamma=1.4)   # 0.824, Eq. (5)
```

A sharp hole (`Cd < 0.82`) gets the orifice form, one rounded enough to behave
like a nozzle (`Cd > 0.94`) gets the nozzle form. Two departures from the
printed equations, both deliberate and both recorded as decisions D11/D12 in
the extraction:

- **The blend weight saturates smoothly.** The paper prints `X = 8.333(Cd -
  0.82)` with no upper bound, which reaches 1.45 for a well-rounded long hole
  and extrapolates past a nozzle. A hard clamp would give an exactly-zero
  `dY/dCd` outside the blend and a discontinuous jump of 0.734 at each knee.
  Pass `eps=0` for the paper's exact hard clamp; the default costs 0.47% in
  `Y` and keeps the derivative continuous and non-zero.
- **`Y` saturates at the critical pressure ratio.** Below `S*` the isentropic
  form predicts *decreasing* flow. The paper has no choked branch.

Provenance, the checks behind every constant, and the two errata found in
the paper are in
`validation/cooling/extractions/orifice_discharge_coefficient.md`.

### Geometry and Flow Analysis
```python
# Area calculations
area = cb.orifice_area(d=0.01)
area = cb.orifice_area_from_beta(beta=0.5, D=0.02)

# Cd-K conversions
K = cb.orifice_K_from_Cd(Cd=0.65)
Cd = cb.orifice_Cd_from_K(K=0.5)

# Flow state analysis
state = cb.orifice_flow_state(P1=2e5, P2=1e5, rho=1.2, mu=1.8e-5, d=0.01, D=0.02)
Re_d = cb.orifice_Re_d_from_mdot(mdot=0.1, d=0.01, mu=1.8e-5)
```

### Advanced Flow Functions
```python
# Thermodynamic orifice flow
sol = cb.orifice_flow_thermo(T=300, P=2e5, X=air, m_dot=0.1, area=1e-4, Cd=0.65)

# Impedance with flow effects
Z = cb.orifice_impedance_with_flow(mdot=0.1, area=1e-4, Cd=0.65, rho=1.2, c=340)

```

### Pressure and Flow Calculations
```python
# Pressure drop
dP = cb.orifice_dP(mdot=0.1, area=1e-4, Cd=0.65, rho=1.2)
dP = cb.orifice_dP_Cd(mdot=0.1, area=1e-4, Cd=0.65, rho=1.2)

# Mass flow
mdot = cb.orifice_mdot(P1=2e5, P2=1e5, area=1e-4, Cd=0.65, rho=1.2)
mdot = cb.orifice_mdot_Cd(P1=2e5, P2=1e5, area=1e-4, Cd=0.65, rho=1.2)

# Velocity and heat transfer
v = cb.orifice_velocity(mdot=0.1, area=1e-4, rho=1.2)
v = cb.orifice_velocity_from_mdot(mdot=0.1, area=1e-4, rho=1.2)
Q = cb.orifice_Q(mdot=0.1, h1=3e5, h2=2.8e5)
```

---

## Heat Transfer

### Correlations
```python
Nu = cb.nusselt_gnielinski(Re=1e5, Pr=0.7)
h = cb.htc_from_nusselt(Nu, k=0.026, L=0.05)

# Multi-layer wall temperature profile
temps, q = cb.wall_temperature_profile(T_hot=1200, T_cold=300, h_hot=200, h_cold=20, t_over_k=[0.01/50, 0.05/0.5])
```

### Network Heat Transfer

#### ConvectiveSurface
Defines convective heat transfer properties on network elements. Supports different channel models.

```python
import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, PinFinModel, RibbedModel, SmoothModel

# Smooth channel with Gnielinski correlation (default)
surface = ConvectiveSurface(
    area=np.pi * 0.04 * 1.0,  # Surface area [m²]
    model=SmoothModel(correlation="gnielinski")
)

# Ribbed walls (see "Ribbed Channels" below)
ribbed = ConvectiveSurface(
    area=2.5,
    model=RibbedModel(e_D=0.06, p_e=10.0, alpha_deg=90.0, W_H=1.0, n_ribbed_walls=2),
)

# Staggered pin-fin array (see "Pin-Fin Channels" below); area is the
# endwall's BASE (planform) area
pins = ConvectiveSurface(
    area=0.01,
    model=PinFinModel(pin_diameter=0.005, S_D=2.5, X_D=2.5, H_D=1.0, N_rows=10),
)
```

#### WallConnection
Thermal coupling between two network elements through a shared wall.

```python
from combaero.heat_transfer import WallConnection

# Couple hot and cold channels through a wall
wall = WallConnection(
    id="coupling_wall",
    element_a="hot_channel",      # First element ID
    element_b="cold_channel",     # Second element ID
    wall_thickness=0.002,      # Wall thickness [m]
    wall_conductivity=25.0,    # Wall thermal conductivity [W/(m·K)]
    contact_area=None          # Optional: override surface area
)

# Add to network
network.add_wall(wall)
```

Each side's heat goes into that element's `to_node` (or the node itself when
the wall names a node). A **boundary** cannot take it -- its state is fixed --
so heat a wall puts there leaves the network with the stream crossing that
boundary, and the result reports it:

```python
res["coupling_wall.Q"]              # W, side a -> side b
res["coupling_wall.Q_to_boundary"]  # W, the part that landed on boundaries (0 if none)
res["outlet.Q_wall_out"]            # W, wall heat leaving through boundary "outlet"
```

With it the network's energy balance closes:
`sum_out m h(T_upstream) + sum_boundaries Q_wall_out = sum_in m h(T_in)`
(`python/tests/test_energy_conservation_cooling.py` checks it to 1e-9 for
every cooling configuration).

- **Reversed flow (#481).** Streams follow the SIGN of an element's flow: a
  2-port element with `m_dot < 0` delivers into its `from_node`, at its
  `to_node`'s state, and a wall puts that side's heat into the node the flow
  goes into. Junction and tee elements keep their declared directions, which
  their own models enforce.
- **k(T) layers** (a `material` from `cb.list_materials()`) take their
  conductivity at the wall temperatures of the same evaluation, iterated to
  a fixed point -- not from the previous one.
- **h <= 0** from a correlation past its range is floored smoothly at
  `cb.WALL_HTC_KNEE` (1e-2 W/(m^2 K)); the wall never refuses it.
- With `thermal_coupling_enabled = False` the walls report nothing.

Every converged solve in the test suite is also checked for this balance
(`python/tests/conftest.py`, `_closure_check.py`).

#### Element Integration
Elements with convective surfaces support heat transfer calculations.

```python
from combaero.network import ChannelElement

# Channel with convective heat transfer
channel = ChannelElement(
    id="heated_channel",
    from_node="inlet",
    to_node="outlet",
    diameter=0.04,
    length=1.0,
    roughness=0.0,
    surface=surface  # ConvectiveSurface instance
)

# Get heat transfer coefficient and temperature
h, T = channel.htc_and_T(state)  # Returns (htc [W/(m²·K)], T [K])
```

#### Thermal Coupling Toggle
Enable/disable thermal coupling globally for debugging or fast solves.

```python
network = FlowNetwork()
network.thermal_coupling_enabled = True   # Enable wall heat transfer
network.thermal_coupling_enabled = False  # Disable (faster, no heat transfer)
```

---

## Acoustics

```python
# Acoustic properties bundle
props = cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=1.0)
print(f"SPL: {props.spl} dB")

# Cavity resonators
f_helm = cb.helmholtz_frequency(V=0.001, A_neck=1e-4, L_neck=0.01, c=340)

# Duct modes
tube = cb.Tube(L=1.0, D=0.1)
modes = cb.tube_axial_modes(tube, c=340, bc1=cb.BoundaryCondition.Closed, bc2=cb.BoundaryCondition.Open)
```

---

## Programmatic Network Design

For automated design and agents, `combaero.network` provides a pure-Python interface to build and solve networks without the GUI.

```python
from combaero.network import FlowNetwork, NetworkSolver, OrificeElement, PressureBoundary, PlenumNode
import combaero.species as species

# 1. Initialize Network
net = FlowNetwork()

# 2. Define Boundaries & Nodes
net.add_node(PressureBoundary("inlet", P_total=5e5, T_total=400, Y=species.dry_air_mass()))
net.add_node(PlenumNode("plenum", P=4e5, T=400))
net.add_node(PressureBoundary("outlet", P_total=1e5, T_total=300))

# 3. Add Elements
net.add_element(OrificeElement("feed", "inlet", "plenum", diameter=0.01, Cd=0.65))
net.add_element(OrificeElement("exhaust", "plenum", "outlet", diameter=0.015, Cd=0.62))

# 4. Solve
solver = NetworkSolver(net)
results = solver.solve(method="hybr")

# 5. Access Results
print(f"Plenum Pressure: {results.nodes['plenum'].P} Pa")
print(f"System Mass Flow: {results.elements['feed'].m_dot} kg/s")
```

> [!TIP]
> **Agent Hint**: Use `net.to_dict()` to get a JSON-serializable representation of the network, which can be saved or sent to the GUI.

---

## Network Solver

### What a solve reports

`NetworkSolver.solve` returns the solved unknowns plus a set of `__dunder__`
bookkeeping keys. Five of them say what happened:

| key | type | meaning |
|---|---|---|
| `__success__` | `bool` | converged **and** consistent. Unchanged; keep using it for "can I trust this result". |
| `__converged__` | `bool` | the root finder reached the residual tolerance, whatever the consistency checks then said |
| `__consistent__` | `bool \| None` | the elements' own physical-consistency verdict. **`None` means not checked** -- either the solve never converged, or no element in this network has a verifier. It is not a pass. |
| `__inconsistent_elements__` | `list[str]` | ids of the elements that rejected the solution |
| `__outcome__` | `SolveOutcome` | why it ended, in one machine-readable value |
| `__worst_residuals__` | `list[dict]` | the rows carrying the residual, largest first |

`__success__` is False in two quite different situations, and reading it alone
cannot tell them apart: Newton never got there, or Newton got there and a
junction rejected the root as unphysical. The second reports a *small*
`__final_norm__` next to `success=False`, which looks like a contradiction
until `__converged__` and `__consistent__` are read separately.

```python
from combaero.network import SolveOutcome

result = solver.solve()
if result["__outcome__"] == SolveOutcome.INCONSISTENT:
    print("converged, then rejected by", result["__inconsistent_elements__"])
elif not result["__success__"]:
    worst = result["__worst_residuals__"][0]
    print(f"{result['__outcome__']}: {worst['name']} carries {worst['residual']:.3e}")
```

`SolveOutcome` is a `StrEnum`, so it compares and serialises as a plain
string: `CONVERGED`, `INCONSISTENT`, `NO_PROGRESS`, `RESIDUAL_TOO_LARGE`,
`TIMEOUT`, `ERROR`, `NOT_CONVERGED`, `NO_UNKNOWNS`.

> [!TIP]
> Prefer `__outcome__` over matching on `__message__`. Part of that text comes
> from SciPy and can be reworded without notice; `__outcome__` is classified
> once, by the code that knows the answer.

Note there is no `__residual_norm__`; the norm is `__final_norm__`.


### Basic Usage
```python
from combaero.network import FlowNetwork, NetworkSolver, OrificeElement

graph = FlowNetwork()
graph.add_element(OrificeElement("ori1", "nodeA", "nodeB", Cd=0.6, diameter=0.011284))

solver = NetworkSolver(graph)
# init_strategy options: "default", "analytical_pt_prop", "homotopy",
# "continuation", "outletref_warmstart", "incompressible_warmstart"
# (deprecated).
# "default" auto-upgrades to "analytical_pt_prop" for networks that
# contain a MultiPortChamberBase; pass another strategy or an
# explicit x0 to opt out.
# "outletref_warmstart" solves the outlet-referenced incompressible
# proxy (densities at the downstream static) and warm-starts the
# compressible solve from it directly -- the auto-retry's seed as a
# primary strategy, for networks known to stall the cold path.
# auto_retry (default True): failed cold solves on compressible
# MPCE networks are retried once from an outlet-referenced
# incompressible warm start (densities at the downstream static).
# The primary attempt gets 40% of `timeout`, the retry the rest.
# A cold solve of a network with walls that still fails is then
# retried from the same network solved WITHOUT its walls -- the
# message says "Converged from the wall-free flow solution ...".
# An iterate the model refuses (an exception, NaN/inf) is answered with a
# residual ten times the last physical one, so the trust region rejects
# the step and shrinks; a failed solve lists the rejected probes in its
# message. Each evaluation repeats the state propagation until states read
# ahead of its order (wall back-edges, recirculation) settle, so the
# residual is a function of x alone.
# Stall detection ends doomed phases early: when hybr's best |F|
# plateaus far from tolerance it hands over to the LM fallback, and
# when the LM fallback plateaus too the attempt fails fast so the
# auto-retry gets the budget. Cold attempts only -- seeded solves
# (explicit x0 or a warm-start seed) are the last rung and keep
# their full budget.
results = solver.solve(method="hybr", init_strategy="homotopy")
```

### Network Nodes
```python
from combaero.network import (
    PressureBoundary, MassFlowBoundary, WallNode, PlenumNode,
    CombustorNode, MomentumChamberNode, MomentumBoundary, EnergyBoundary
)

# Boundary conditions
inlet = PressureBoundary("inlet", P_total=2e5, T_total=300)
outlet = PressureBoundary("outlet", P_total=1e5, T_total=300)
mass_in = MassFlowBoundary("mass_in", m_dot=0.1, T_total=400, Y=air)
mom_in = MomentumBoundary("mom_in", P_total=2e5, T_total=300, area=1e-4)

# Closed-end boundary (zero mass flow)
# Typical use: symmetry plane of a ring manifold, or any dead-end branch.
wall = WallNode("wall")

# Internal nodes
junction = PlenumNode("junction", P=1.5e5, T=350, Y=air)
combustor = CombustorNode("combustor", P=2e5, T=1800, Y=products)
momentum = MomentumChamberNode("momentum", P=2e5, T=1800, Y=products, regime="compressible")
energy = EnergyBoundary("energy", Q=50000)  # Heat addition [W]
loss = EnergyBoundary("loss", fraction=-0.05)  # removes 5% of the SENSIBLE enthalpy
```

- **`fraction`** scales the inflow's sensible enthalpy, `h(T) - h(298.15 K)`
  at its own composition (`cb.SENSIBLE_ENTHALPY_REF_T`); on a combustor, the
  products'. It is not a fraction of absolute enthalpy, which includes
  formation enthalpy (#481). The heat it applied is reported as
  `{node}.Q_fraction` [W].
- **`Q`** is spread over the node's real inflow. Only a stagnant node (below
  `cb.MIXER_HEAT_MDOT_FLOOR` = 1 mg/s) cannot take it; what it withholds is
  reported as `{node}.Q_withheld` [W], as is wall heat aimed at a `WallNode`.
- **A `MassFlowBoundary` between elements** is an injection: its `m_dot`
  joins the node at its own `Tt` and `Y`.

#### Merge chamber: main inlet plus side streams (#471)

A `MomentumChamberNode` with `main_inlet` declared accepts further inflows as
SIDE STREAMS, still with exactly one outlet. The chamber's axis is its outlet
direction; the main stream arrives along it; each side stream brings
`m u_jet cos(theta)` of axial momentum and discharges at the chamber's static
pressure. A constant-area control volume closes it:

    P_face A + m_main u_main + sum(m_s u_jet cos theta) = P A + m_out u_out

```python
liner = MomentumChamberNode("liner", main_inlet="duct")   # area inherited from "duct"
panel = EffusionPlateElement("eff", "coolant", "liner", hole_diameter=0.8e-3,
                             wall_thickness=2e-3, pitch=6e-3, panel_area=0.01,
                             angle_deg=30.0)               # jets at 30 deg to the gas
```

- **The node's `(P, Pt)` is the chamber (outlet) state.** The main-inlet
  element sees the MAIN-FACE state the impulse balance gives
  (`cb.chamber_merge_face_state`): exact for compressible flow, the face
  static pressure on the subsonic root and its Pt from the same
  entropy-based closure as the chamber. With no side streams the face is the
  chamber, so declaring a single-inlet chamber's main inlet changes nothing.
- **Jets are compressible too.** A side stream's impulse is `cb.jet_impulse`:
  the isentropic velocity from its supply (Pt, Tt) to the chamber pressure,
  or, choked, the sonic momentum plus its pressure thrust.
- **No loss coefficient.** The stagnation-pressure loss of mixing follows
  from momentum; transverse momentum is reacted by the walls and its kinetic
  energy dissipated. Normal injection costs the main stream its mixing loss;
  inclined jets push it, ejector-like.
- **Directions live only inside the junction.** Edges carry scalars and a
  channel is 1D along its own axis, so no angle is ever passed on.
- **One area.** The main inlet's port area is the chamber area, inherited
  from the main-inlet element; an explicit area that disagrees is refused.
- **The injector owns its jet.** Orifice-type elements carry
  `injection_angle_deg` (90 by default, no axial momentum);
  `EffusionPlateElement` uses its hole inclination `angle_deg`. The jet
  velocity is the element's own, from its supply's stagnation state.
- **Refused:** two inflows without `main_inlet`, a `main_inlet` that does not
  flow in, and any split (more than one outlet needs a loss closure: use a
  tee).
- In the GUI, inflows on the chamber's "s" handle are side streams; the one on
  its flow handle is the main inlet.

### Network Elements
```python
from combaero.network import (
    OrificeElement, ChannelElement, EffectiveAreaConnectionElement,
    LosslessConnectionElement, DiameterDischargeCoefficientConnectionElement,
    TeeJunctionElement, VortexElement, EffusionPlateElement,
    BorderCarnotLossElement,
)
from combaero.network.mpce_element import ConstantKTeeElement, MultiPortChamberElement

# Flow elements
orifice = OrificeElement("orifice", "node1", "node2", Cd=0.65, diameter=0.011284, regime="compressible")
channel = ChannelElement("channel", "node2", "node3", length=2.0, diameter=0.05, roughness=1e-4, regime="compressible")

# Effusion (multi-perforated) wall panel. Geometry is given the way a plate
# is designed -- pitch, hole diameter, wall thickness, inclination -- and the
# hole count follows from the panel area.
panel = EffusionPlateElement(
    "panel", "coolant", "gas",
    hole_diameter=0.6e-3, wall_thickness=1.0e-3,
    pitch=3.0e-3, panel_area=0.02 * 0.02, angle_deg=30.0,
)
# One panel is one coolant pressure and one gas pressure, so it cannot show
# coolant migration WITHIN itself. For a wall along a feed channel, hang a
# panel off each channel segment -- the node mass balance then gives
# m_channel_in = m_channel_out + m_effusion, i.e. coolant flow as f(x). The
# resolution is your choice of segment count.

# Connection elements
effective_area = EffectiveAreaConnectionElement("ea", "node3", "node4", diameter=0.015958)
lossless = LosslessConnectionElement("lossless", "node4", "node5")
area_cd = DiameterDischargeCoefficientConnectionElement("area_cd", "node5", "node6", diameter=0.011284, Cd=0.7)

# Three-port tee junction (Unified0D compressible model, Mynard & Valen-Sendstad 2015)
import math
# Merging tee: two inlets (straight_node, branch_node) -> one outlet (common_node)
tee_merge = TeeJunctionElement(
    "tee", common_node="outlet", straight_node="inlet_s", branch_node="inlet_b",
    theta=math.pi / 2.0, F_C=0.01, psi=1.0, tee_type="merging",
)
# Branching tee: one inlet (common_node) -> two outlets (straight_node, branch_node)
tee_branch = TeeJunctionElement(
    "tee", common_node="inlet", straight_node="outlet_s", branch_node="outlet_b",
    theta=math.pi / 2.0, F_C=0.01, psi=1.0, tee_type="branching",
)
# Solved unknowns: tee.m_dot_com (total), tee.m_dot_branch (branch arm)
# Straight flow is implicit: m_dot_straight = m_dot_com - m_dot_branch

# Momentum-CV junction (PDF spec, supersedes K-closure for n>3 manifolds and
# high-Mach / ejector behaviour; see docs/junction/momentum cv implementation guide.pdf).
# N port-MCNs feed a single junction; each lateral port carries a turning loss.
mc_com = MomentumChamberNode("mc_com", area=0.01)
mc_str = MomentumChamberNode("mc_str", area=0.01)
mc_bra = MomentumChamberNode("mc_bra", area=0.008)

# MultiPortChamberBase is the ABSTRACT BASE since 0.6.0 -- it owns the port
# machinery and no longer carries a residual. Use MultiPortChamberElement (the Mynard
# closure) or ConstantKTeeElement (fixed handbook K).
jct = MultiPortChamberElement(
    "jct",
    inlet_nodes=["mc_com"],
    outlet_nodes=["mc_str", "mc_bra"],
    inlet_angles_deg=[0.0],
    outlet_angles_deg=[0.0, 90.0],   # geometric branch angles per outlet port
    flow_direction="branch",         # 1 supplier + N-1 collectors
    # port_areas (length = N_in + N_out) inherited from connecting channels if not given
)
# Lateral port turning loss (sharp-edged 90 deg branch). Straight ports
# (delta_geom = 0) do not need a loss element.
loss_bra = BorderCarnotLossElement(
    "loss_bra", from_node="mc_bra", to_node="ch_bra_in",
    delta_geom_deg=90.0,
    # area inherited from neighbouring channel if not given
)
# Solved unknowns: jct.P_jct (junction static pressure). Per-port mass flows
# live on the connecting channels / loss elements; the junction reads them via
# the graph. The N+1 residuals = N port total-pressure relations
# (Pt_i - Pt_jct +/- K_i * q_dyn_com) + global sum-of-port-mdots = 0, computed
# in C++ (see mpce_residuals_and_jacobian in API_CPP.md).

# Constant-K junction (the "simplest model" tier): fixed handbook loss
# coefficients instead of the Mynard closure. K is referenced to the
# common-leg dynamic head (Idelchik convention) and does not vary with the
# flow split. K_ports keys are port indices (inlets first, then outlets);
# the common port (single inlet for 'branch', single outlet for 'merge')
# is ignored. Exposed in the GUI as the tee's "Constant K (handbook)"
# junction model.
from combaero.network.mpce_element import ConstantKTeeElement
jct_k = ConstantKTeeElement(
    "jct",
    inlet_nodes=["mc_com"],
    outlet_nodes=["mc_str", "mc_bra"],
    flow_direction="branch",
    K_ports={1: 0.4, 2: 1.0},   # K_straight, K_branch
)

# Rotating cavity / disc-pump pressure rise (Vatistas n-vortex model)
vortex = VortexElement(
    "vortex", "node_in", "node_out",
    r_c=0.02,       # vortex core radius [m]
    r_out=0.10,     # outer evaluation radius [m]
    r_in=0.0,       # inner evaluation radius [m] (0 = on-axis)
    omega_rpm=3000, # shaft speed [rpm]
    n=2.0,          # Vatistas shape parameter (>= 1, default n=2)
)
# Residual: Pt_out - Pt_in - dP_vortex = 0
# dP_vortex = vatistas_delta_p(r_out, ...) - vatistas_delta_p(r_in, ...)
# Gamma is derived from omega so that V_theta(r_c) = omega * r_c exactly (solid-body core).
# Diagnostics include: m_dot, omega_rpm, dP_vortex [Pa], Pt_in, Pt_out.
# If omega_rpm=0 the element is lossless (transparent); r_in=0 omits the on-axis term.

# Supersonic ejector across all operating regimes -- critical (double-choked),
# subcritical droop, and unchoked-primary subsonic jet pump -- in one C1 residual
# system (Huang 1999 entrainment ratio + Kracik & Dvorak 2016 mixing closure;
# see validation/ejector/README.md and OPERATING_REGIMES_DESIGN.md).
from combaero.network.ejector_element import EjectorElement
ejector = EjectorElement(
    "ej1", primary_node="mc_primary", secondary_node="mc_secondary", outlet_node="mc_outlet",
    throat_area=3.14e-5,      # primary nozzle throat area A_t [m^2]
    nozzle_exit_area=1.0e-4,  # primary nozzle exit area A_p1 [m^2], must exceed throat_area
    mixing_area=8.0e-4,       # constant-area mixing section A_3 [m^2], must exceed nozzle_exit_area
    recovery_efficiency=1.0,  # multiplies the lossless mixed stagnation pressure; not fitted to data
)
# Solved unknowns: the three port mass flows + ej1.P_jct, the single owned scalar,
# repurposed as the mixing-plane static pressure P_py. Four residual rows: primary
# C-D nozzle (choked or unchoked), blended entrainment, mass conservation, and a
# regime-dependent outlet/discharge closure. The regime is picked smoothly by two
# smootherstep weights (s_choke on mp/mdot_choked, s_sub on outlet.Pt/P_c*): in
# critical mode the outlet floats free below P_c* and P_jct is the diagnostic
# recovered stagnation; in the subcritical/jet-pump regimes the outlet is pinned
# to the ejector's discharge and P_py floats. Working fluid: combaero's real-fluid
# EOS (combustion/air species only, no halocarbon refrigerants); gamma/R evaluated
# live at the entrained-flow choking plane. The whole 4-row (f, J) is assembled in
# C++ (full analytic Jacobian, no FD fallback). diagnostics() reports the actual
# omega = ms/mp, a critical_mode flag, and (in critical mode) p_c_star_pa.
# verify_solution_consistent(sol) returns True unconditionally now that every
# regime is modeled.
```

### Flow-Area Inference

Components that need a reference flow area infer one from the topology when you
do not supply it. Inference searches the neighbourhood of the attached nodes:
channels first (so networks that resolved before resolve identically), then any
other area-bearing element (`AreaChangeElement`, `TeeJunctionElement`,
`BorderCarnotLossElement`, `PressureLossElement`), then adjacent nodes whose
area you set explicitly (`MomentumChamberNode`, `CombustorNode`). An area that
was itself auto-sized is never used as a source -- inheriting a placeholder
would launder one guess into another.

What happens when nothing can be inferred depends on what the area means:

| Component | Parameter | Unresolvable |
|---|---|---|
| `AreaChangeElement` | `F0`, `F1` | **raises** |
| `ChannelElement` | `diameter` | **raises** |
| `TeeJunctionElement` | `F_C` | **raises** |
| `BorderCarnotLossElement` | `area` | **raises** |
| `MomentumChamberNode` | `area` | nominal default |
| `CombustorNode` | `area` | nominal default |
| `PressureLossElement` | `area` | nominal default, warns if a head-loss correlation reads it |

The first group raises because the area *is* the geometry being modelled -- a
default there invents a contraction or expansion the network does not contain,
and the solver either stalls on it or converges to a confidently wrong answer.
The second group defaults because the area is only a dynamic-head reference and
a network can legitimately never read it.

```python
# Unresolvable -> a named error saying which parameter clears it
AreaChangeElement("expansion", "plenum_a", "plenum_b")
# ValueError: AreaChangeElement 'expansion': cannot determine upstream area F0
# -- no neighbouring channel, area-bearing element or node with a known area
# was found at 'plenum_a'. Set upstream area F0 explicitly ...

# Either fix clears it: set the value, or give a neighbour a known area
AreaChangeElement("expansion", "plenum_a", "plenum_b", F0=0.05, F1=0.08)
MomentumChamberNode("chamber", area=0.15)  # neighbour the area change can read
```

### Tuned Constants in the Junction Closure

`MultiPortChamberElement`'s Mynard closure carries three tuned constants, each named
and documented at its definition in `combaero.network._mynard2010`:

| constant | value | what it is | switch |
|---|---|---|---|
| `MYNARD_ETA_A0`, `MYNARD_ETA_A1` | 0.8, -0.2 | Mynard 2015 Eq 36 energy-transfer factor, CFD-fitted | `eta_scale` |
| `FLOW_RATIO_DAMPING` | 0.02 | regulariser from the Matlab reference, not physics | -- |
| `MultiPortChamberElement.DEFAULT_JOINING_ETRANSFER_ALPHA` | 0.2 | combaero's joining-side correction, fitted by `validation/junction/calibrate_etransfer.py` | pass `joining_etransfer_alpha=0.0` |

Under this repo's policy a tuned constant earns its place only by improving
agreement on the digitised validation data; the on/off tables live on issue
#271. `eta_scale` exists for that measurement:

```python
MultiPortChamberElement(..., eta_scale=1.0)   # default: the faithful Mynard port
MultiPortChamberElement(..., eta_scale=0.0)   # energy-transfer term off, for scoring
junction_loss_coefficient(U, A, theta, eta_scale=0.0)
```

Changing any of these values is a retune and needs a before/after table
before it lands, not a silent edit -- `python/tests/test_junction_tuned_constants.py`
pins the documented values.

### The Junction Soft-Barrier Weight

Distinct from the tuned constants above: this is a **numerical** parameter, not
physics, and it is derived rather than declared.

When `strict=False` and a port flows against its declared direction,
`MultiPortChamberElement` replaces the physics with a continuity residual plus a
one-sided quadratic penalty, `alpha * max(0, -e_i * mdot_i)^2`. The penalty
shares its row with the continuity relation, so it balances against a pressure
error rather than driving the offending flow to zero, and has a fixed point at
`slack* = sqrt(dP / alpha)`. A solve that reaches it parks there.

`alpha` therefore carries `Pa/(kg/s)^2` and cannot be a constant: the weight
needed scales as `1/m_ref^2`, two decades per decade of network size.
`NetworkSolver` derives it before each solve from the reference state it
already computes for seeding:

```python
alpha = P_ref / (BARRIER_SLACK_FRACTION * m_ref) ** 2   # f = 0.005
```

placing the fixed point at 0.5% of the reference mass flow whatever the
network's size. It is frozen for the solve, so the residual and Jacobian are
unchanged in form.

| name | where | meaning |
|---|---|---|
| `BARRIER_SLACK_FRACTION` | `combaero.network.mpce_element` | where the fixed point is placed, as a fraction of `m_ref` |
| `DEFAULT_SOFT_PENALTY_ALPHA` | same | fallback when no solver has supplied a weight, or the reference state is degenerate |
| `scaled_penalty_alpha(P_ref, m_ref)` | same | the derivation, exposed for testing |
| `MultiPortChamberElement.effective_penalty_alpha()` | element | the weight actually used, after precedence |

Precedence: an explicitly set `soft_penalty_alpha` always wins, then the
solver-supplied scale-aware weight, then the fallback. So the tuning knob keeps
working:

```python
element.soft_penalty_alpha = 5.0e7   # explicit: overrides the derived weight
```

Raising `alpha` shrinks the fixed point and never destabilises the solve --
the response is monotone and saturates -- so when in doubt, larger.

### Combustion Integration
```python
from combaero.network.combustion import (
    combustion_from_streams, combustion_from_phi, mix_streams, stoichiometric_products
)

# Create combustion from streams
fuel_stream = cb.Stream(m_dot=0.01, T=300, P_total=2e5, Y=fuel_y)
oxidizer_stream = cb.Stream(m_dot=0.2, T=300, P_total=2e5, Y=air_y)

result = combustion_from_streams([fuel_stream, oxidizer_stream], phi=1.0)
print(f"Products: {result.products.Y}")
print(f"Adiabatic T: {result.T_ad} K")

# Or from equivalence ratio
result = combustion_from_phi(phi=0.8, T=300, P=2e5, fuel="CH4", oxidizer="air")
```

### MixtureState for Network Data
```python
from combaero.network import MixtureState

# Constructor: MixtureState(P, P_total, T, T_total, m_dot, Y)
node_state = MixtureState(1e5, 1.1e5, 300, 300, 1.5, air_vec)

print(f"Static pressure: {node_state.P} Pa")
print(f"Total pressure: {node_state.P_total} Pa")
print(f"Mass flow: {node_state.m_dot} kg/s")
print(f"Static density: {node_state.density()} kg/m³")
```

### Advanced (f, J) Interface
High-accuracy analytical Jacobians are available for performance-critical solver loops.
```python
# Exact residuals and derivatives wrt (P_tot, P_static, T, Y)
res_ori = cb.orifice_residuals_and_jacobian(m_dot, P_tot, P_stat, T, Y, P_down, Cd, area)
res_chan = cb.channel_residuals_and_jacobian(m_dot, P_tot, P_stat, T, Y, P_down, L, D, roughness, correlation)
res_comb = cb.combustor_residuals_and_jacobians(m_dot, P_tot, P_stat, T, Y, Q_comb, method, smooth, pressure_loss_func)
res_plen = cb.plenum_residuals_and_jacobian(m_dot_vec, P_target, T_target, Y_target)
```

Combustor results return mapping specific to `(m_dot, P_total, T, Y)` state vectors.

### Tee Junction Interface

Two solver-facing APIs are available:

**Unified0D compressible model** (current, used by `TeeJunctionElement`): stagnation-pressure
residuals with mass-flow continuity, consistent with compressible duct elements. Supports
arbitrary per-branch gamma and R_gas.

**Legacy Bassett 2001 model**: empirical incompressible K tables retained for direct
low-level access; not used by the network solver.

#### Unified0D compressible interface

`BranchInput` holds per-branch thermodynamic state. `CompressibleTeeResult` holds the
two residuals and the full Jacobian.

```python
import combaero._core as _core
import math

# Build per-branch state (P_static [Pa], Pt [Pa], T [K], m_dot [kg/s],
#                         A [m^2], theta [rad], gamma_eff [-], R_gas [J/kg/K])
com = _core.BranchInput(P_static=1.95e5, Pt=2.0e5, T=400.0, m_dot=0.5,
                         A=0.01, theta=0.0, gamma_eff=1.4, R_gas=287.0)
str_b = _core.BranchInput(P_static=1.88e5, Pt=1.95e5, T=400.0, m_dot=0.3,
                           A=0.01, theta=0.0, gamma_eff=1.4, R_gas=287.0)
bra = _core.BranchInput(P_static=1.82e5, Pt=1.90e5, T=400.0, m_dot=0.2,
                         A=0.01, theta=math.pi/2, gamma_eff=1.4, R_gas=287.0)

# Branching: common=supplier, straight+branch=collectors
res = _core.compressible_branching_tee_rj(com=com, str=str_b, bra=bra)
# res.R_0, res.R_1           -- residuals [Pa]
# res.dR0_dPt_com, res.dR0_dPt_str, ...  -- full Jacobian fields

# Merging: straight+branch=suppliers, common=collector
res = _core.compressible_merging_tee_rj(com=com, str=str_b, bra=bra)
```

#### Legacy Bassett 2001 interface

```python
import combaero._core as _core
import math

Y = [0.0] * 15
Y[0] = 0.767  # N2
Y[1] = 0.233  # O2

# Check input validity (non-throwing)
status = _core.tee_check_inputs(q=0.5, psi=1.0, theta=math.pi/2)
# status.valid, status.q_in_range, status.psi_valid, status.theta_valid

# Raw K-coefficient functions (Bassett 2001 Table 2 + Eq 33/34 angle corrections).
# K2 and K5 take only q (psi/theta-independent per Eq 15); the rest take (q, psi, theta).
# Provided for diagnostics and validation against measured data.
k1  = _core.tee_K1(q, psi, theta)
k2  = _core.tee_K2(q)
k3  = _core.tee_K3(q, psi, theta)
k4  = _core.tee_K4(q, psi, theta)
k5  = _core.tee_K5(q)               # straight arm, separating tee
k6  = _core.tee_K6(q, psi, theta)   # branch arm, separating tee
k7  = _core.tee_K7(q, psi, theta)
k8  = _core.tee_K8(q, psi, theta)
k9  = _core.tee_K9(q, psi, theta)
k10 = _core.tee_K10(q, psi, theta)
k11 = _core.tee_K11(q, psi, theta)  # straight arm, joining tee
k12 = _core.tee_K12(q, psi, theta)  # branch arm, joining tee

# Blended K functions (always finite, smooth at q=0 and across topology reversal)
ks = _core.merging_tee_K_straight(q, psi, theta)
kb = _core.merging_tee_K_branch(q, psi, theta)

# Legacy solver residuals and full Jacobians
res = _core.merging_tee_residuals_and_jacobian(
    m_dot_com=0.5, m_dot_branch=0.2,
    dP0_straight=150.0, dP0_branch=200.0,
    P_static_com=2e5, T_com=600.0, Y_com=Y,
    theta=math.pi/2, psi=1.0, F_C=0.01,
)
# res.R_straight, res.R_branch
# res.dR_straight_d_mdot_com, res.dR_straight_d_mdot_branch
# res.dR_straight_dP_static_com, res.dR_straight_dT_com
# res.topology_valid, res.status  (CorrelationValidity.VALID or EXTRAPOLATED)

res = _core.branching_tee_residuals_and_jacobian(
    m_dot_com=0.5, m_dot_branch=0.15,
    dP0_straight=100.0, dP0_branch=180.0,
    P_static_com=1.5e5, T_com=500.0, Y_com=Y,
    theta=math.pi/3, psi=1.2, F_C=0.008,
)
```

---

## Vatistas n-Vortex Model

Reference: Vatistas, Kozel, Mih (1991), *Exp. Fluids* 11, 73–76.

Shape parameter `n >= 1` (default 2.0). `n=2` gives the best fit to most
experimental swirl data.  Closed-form antiderivatives for the pressure integral
are used for `n=1` and `n=2`; composite Simpson quadrature for general `n`.
All Jacobians are analytical.

### Free functions

```python
import combaero as cb

# Normalised shape functions (no units, r_bar = r / r_c)
v0   = cb.vatistas_v0_bar(r_bar, n)           # tangential velocity shape
dv0  = cb.vatistas_dv0_bar_drbar(r_bar, n)    # d(V0_bar)/d(r_bar)
vr   = cb.vatistas_vr_bar(r_bar, n)           # radial velocity shape (< 0)
dvr  = cb.vatistas_dvr_bar_drbar(r_bar, n)    # d(Vr_bar)/d(r_bar)
I    = cb.vatistas_pressure_integral(r_bar, n) # integral_0^r V0^2/r' dr'
dI   = cb.vatistas_d_pressure_integral_drbar(r_bar, n)  # V0_bar^2/r_bar (FTC)

# Dimensional functions (Gamma [m^2/s], r_c [m], rho [kg/m^3])
Vt   = cb.vatistas_v_theta(r, Gamma, r_c, n)  # [m/s]
Vt, dVt_dr, dVt_dG, dVt_drc = cb.vatistas_v_theta_and_jacobians(r, Gamma, r_c, n)

dP   = cb.vatistas_delta_p(r, rho, Gamma, r_c, n)  # P(r)-P(0) [Pa]
dP, ddP_dr, ddP_dG, ddP_drc = cb.vatistas_delta_p_and_jacobians(r, rho, Gamma, r_c, n)
```

### `VatistasVortex` class

```python
import numpy as np
import combaero as cb

vortex = cb.VatistasVortex(Gamma=0.5, r_c=0.05, n=2.0)

r = np.linspace(0, 0.2, 200)
V  = vortex.V_theta(r)          # [m/s], vectorised
dV = vortex.dV_theta_dr(r)      # [1/s], vectorised

Vmax = vortex.V_theta_max()     # peak tangential velocity at r = r_c

dP   = vortex.delta_P(r=0.1, rho=1.2)    # [Pa], scalar
dPbar = vortex.delta_P_bar(r=0.1)        # normalised [0, 1)
```

---

## NetworkRunner (GUI-JSON Programmatic Driver)

`NetworkRunner` loads a network saved from the GUI and runs it from Python
scripts or Jupyter notebooks — no web server required.

### Loading a network

```python
from gui.backend.runner import NetworkRunner

# From a file downloaded from the GUI
runner = NetworkRunner.from_file("my_network.json")

# From an already-loaded dict (normalises GUI camelCase keys automatically)
import json
with open("my_network.json") as f:
    runner = NetworkRunner.from_dict(json.load(f))
```

Node labels are read from the GUI's label field.  If a node has no label the
node ID is used as the fallback key for overrides and result lookup.

### Single solve with boundary-condition overrides

```python
# Override format: "<label>.<attribute>": value
result = runner.solve({
    "air_inlet.m_dot": 1.2,   # kg/s
    "air_inlet.Tt": 650.0,    # K
    "outlet.Pt": 200_000.0,   # Pa
})

print(result.success, result.final_norm)

# Scalar result by label.quantity (resolves label → node ID automatically)
T_combustor = result.get("combustor.T")         # K
phi          = result.get("combustor.phi")
m_dot_out    = result.get("hot_channel.m_dot")  # kg/s

# Full thermodynamic state dict for a node
state = result.node_state("combustor")
# {"T": 1738.0, "P": 195000.0, "Pt": 200000.0, "Tt": 1760.0, ...}

# Full-detail DataFrame (matches GUI CSV export, unit-annotated columns)
df = result.to_dataframe()
```

### Parametric sweep

```python
import pandas as pd

params = pd.DataFrame({
    "fuel_inlet.m_dot": [0.020, 0.025, 0.030, 0.035],
})

# Compact: one row per solve, only requested metrics
sweep_df = runner.sweep(params, metrics=["combustor.T", "combustor.phi"])
# columns: fuel_inlet.m_dot | combustor.T | combustor.phi | success | final_norm

# Full detail: complete to_dataframe() output for each solve, stacked
full_df = runner.sweep(params)
# columns: _sweep_index | fuel_inlet.m_dot | <all entity columns> | success | final_norm
```

### Injecting labels for unlabelled networks

GUI exports may omit node labels.  Inject them before constructing the runner:

```python
import json

with open("network.json") as f:
    schema = json.load(f)

label_map = {
    "node_1778698958490": "air_inlet",
    "node_1778698961265": "fuel_inlet",
    "node_1778699014636": "combustor",
}
for node in schema["nodes"]:
    if node["id"] in label_map:
        node["data"]["label"] = label_map[node["id"]]

runner = NetworkRunner.from_dict(schema)
```

### NetworkResult attributes

| Attribute | Type | Description |
|---|---|---|
| `success` | `bool` | True when solver converged |
| `message` | `str` | Human-readable solver status |
| `final_norm` | `float \| None` | Residual norm at convergence |

| Method | Returns | Description |
|---|---|---|
| `get(key)` | `float` | Scalar by `<id>.<qty>` or `<label>.<qty>` |
| `node_state(label)` | `dict` | Full state dict for a node |
| `to_dataframe()` | `DataFrame` | Unit-annotated full-detail export |
| `swap_boundary(label, new_type)` | `NetworkRunner` | Retype a boundary node (mass ↔ pressure) and return a new runner pre-seeded with the full solved state |

### Boundary condition swap

```python
result = runner.solve({"air_inlet.m_dot": 1.0})

# Re-run the same point with a pressure BC (Pt auto-populated from solved state)
p_runner = result.swap_boundary("air_inlet", "pressure_boundary")
result2 = p_runner.solve()  # starts from solved state, converges immediately

# Sweep inlet pressure around the design point
import pandas as pd, numpy as np
design_Pt = result.node_state("air_inlet")["Pt"]
params = pd.DataFrame({"air_inlet.Pt": design_Pt * np.linspace(0.90, 1.10, 11)})
sweep_df = p_runner.sweep(params, metrics=["combustor.T"])
```

Only `mass_boundary` ↔ `pressure_boundary` swaps are supported.  The returned
runner carries the full solved state (pressures, flows, junction states) as a
warm start so `solve()` requires no additional `init_strategy`.

---

## Geometry & Materials

### Geometry Classes
```python
# Basic geometry
tube = cb.Tube(L=1.0, D=0.1)
annulus = cb.Annulus(R_outer=0.1, R_inner=0.05)

# Can-annular geometry
can_geo = cb.CanAnnularFlowGeometry(
    N_cans=12, R_can=0.05, L_can=0.3,
    R_annulus_inner=0.06, R_annulus_outer=0.1
)
```

### Material Properties
```python
# List available materials
materials = cb.list_materials()

# Thermal conductivity [W/(m·K)]
k_al = cb.k_aluminum_6061
k_ss = cb.k_stainless_steel_316
k_inconel = cb.k_inconel718
k_haynes = cb.k_haynes230
k_tbc = cb.k_tbc_ysz
```

### Advanced Geometry Functions
```python
# Annular areas
area_ann = cb.annular_area(R_outer=0.1, R_inner=0.05)

# Residence times
t_res = cb.residence_time(V=0.01, Q=0.1)
t_res_ann = cb.residence_time_annulus(V=0.01, Q=0.1, R_outer=0.1, R_inner=0.05)
t_res_can = cb.residence_time_can_annular(V=0.01, Q=0.1, N_cans=12, R_can=0.05)
```

### Rib Correlations (parametrised)

Rib correlations are **data, not code**. A parameter set carries the
coefficients, their normalisers, the validity band, the stated accuracy, and
where all of it came from -- so a built-in set and a set you tuned on your own
rig are structurally distinguishable.

```python
import combaero as cb

s = cb.han_1988_orthogonal()          # Han (1988), 90 deg orthogonal ribs
g = cb.RibGeometry(e_D=0.047, p_e=10.0, W_H=1.0, alpha_deg=90.0)

r = cb.evaluate_rib(s, g, Re=10_000)
r.R, r.f, r.e_plus, r.G, r.St_r       # 3.2000, 0.04576, 71.1, 12.21, 0.00968
r.extrapolated                        # outside the set's advisory validity
```

**Two named sets, chosen explicitly.** They cover disjoint Reynolds-number
regimes reported by different papers, not one superseding the other -- so
there is no auto-switching between them:

```python
s90 = cb.han_1988_orthogonal()        # 90 deg, Re 10,000-60,000
s45 = cb.rallabandi_2009_high_re()    # 45 deg sharp ribs, Re 30,000-400,000

g = cb.RibGeometry(e_D=0.14, p_e=7.5, W_H=1.0, alpha_deg=45.0)
r = cb.evaluate_rib(s45, g, Re=100_000)
r.G                                    # 35.40 -- see han_ribbed_high_re.md
```

`rallabandi_2009_high_re` is 45 deg and sharp-edged ribs only: the source
reports round-edged ribs at the same conditions instead following Han's
correlation, but calls the agreement "coincidental" rather than physically
grounded, so it is not wired as an automatic edge-profile switch.

**A third set covers angled ribs generally**, `han_park_1988_angled` (Eq.
4.17/4.18, `alpha` 30-90 deg, `W/H` 1-4):

```python
angled = cb.han_park_1988_angled()
g = cb.RibGeometry(e_D=0.06, p_e=15.0, W_H=2.0, alpha_deg=60.0)
r = cb.evaluate_rib(angled, g, Re=30_000)
r.R, r.G                              # 3.2332, 18.4408
```

Its `R` and `G` carry genuine, unsmoothed discontinuities the source
states rather than a numerical artefact: `R`'s `(W/H)^m` term switches at
`alpha == 90 deg` (up to 62.5% at `W/H = 4`), and `G`'s `alpha`/`p_e`
exponents switch on whether the channel is square (`W/H == 1`, up to 27%
at low `alpha`). Neither the printed equations nor the extraction record a
smooth transition between the branches, so none is invented -- a caller
whose solver traverses either boundary exactly needs its own guard.

`RAlphaShape` and `GShapeModel` (both importable from `combaero`) are what
make this representable as data: setting `R_alpha_shape` to
`QuadraticAlpha` and populating `R_quad_c0/c1/c2` plus the two
`R_quad_WH_exponent_*` fields lets a user-supplied set describe its own
angle-dependent `R`, the same way `RibTerm`'s normaliser lets one describe
its own geometry dependence.

**A fourth set extends the same family to narrow channels**,
`han_1989_narrow_channel` (Eq. 4.19, `alpha` 30-90 deg, `W/H` 1/4-1):

```python
narrow = cb.han_1989_narrow_channel()
g = cb.RibGeometry(e_D=0.0625, p_e=15.0, W_H=0.3, alpha_deg=45.0)
r = cb.evaluate_rib(narrow, g, Re=30_000)
r.R, r.G
```

Its source has **no printed equation for `R`** below `W/H = 1` -- only a
drawn figure line -- so `R_alpha_shape = QuadraticAlphaTwoBand`'s
coefficients (`R_quad_c0/c1/c2` for `1/2 <= W/H < 1`,
`R_quad_narrow_c0/c1/c2` for `1/4 < W/H < 1/2`, switched by
`R_WH_band_boundary`) are this project's own fit to that line, and the
set's `provenance` is `Fitted` rather than `Extracted` -- see
`validation/cooling/extractions/han_ribbed.md`, item 40. `G`
(`GShapeModel.NarrowChannelAlphaSwitch`) *is* text-extracted: a genuine
~20% jump in its leading constant at `alpha == 90 deg`
(`G_narrow_C_alpha90` vs `G_narrow_C_off_axis`), plus a further `(W/H)`
correction to both the constant and the `e+` exponent itself below
`G_narrow_WH_band_boundary`.

**Rib shape, and nine configurations from one paper.** A set declares the
`RibShape` it was fitted to (`Transverse`, `Parallel`, `Crossed`, `V`,
`Lambda`; `Unspecified` binds nothing) because V, crossed and Lambda ribs at
the same angle are different geometries. `han_zhang_lee_1991` returns Han,
Zhang & Lee (1991)'s Table 2 row for a tested configuration -- one rig,
square channel, `e/D` 0.0625, `P/e` 10 -- and raises for anything else:

```python
vee = cb.han_zhang_lee_1991(cb.RibShape.V, 60.0)
g = cb.RibGeometry(e_D=0.0625, p_e=10.0, W_H=1.0, alpha_deg=60.0)
r = cb.evaluate_rib(vee, g, Re=30_000)
r.G, r.G_bar, r.has_G_bar           # its own printed G_bar, not 1.2 G
vee.symmetric                       # False: a V reversed is a Lambda
```

**Ratio-form sets** (`RibRatioSet`, #444) carry `Nu/Nu0` and `f/f0` as
multipliers on the smooth-duct baseline the source normalised by --
`Nu0_source` is required, because a ratio without its reference is not a
number. In the fitted range they reproduce the paper; below `Re_floor` they
hand over smoothly to Gnielinski by default, so Nu keeps a laminar floor as
Re -> 0 rather than collapsing with Dittus-Boelter:

```python
s = cb.RibRatioSet()
s.name, s.source = "my_rig", "measured 2026-10"
s.C_Nu, s.Nu_Re = 2.5, cb.RibTerm(-0.2, 1.0e4)       # Nu/Nu0 = 2.5 (Re/1e4)^-0.2
s.C_f, s.f_Re = 8.0, cb.RibTerm(0.1, 1.0e4)           # f/f0  = 8 (Re/1e4)^0.1
s.Nu0_source = cb.RatioBaseline(0.023, 0.8, 0.4)     # Dittus-Boelter
s.f0_source = cb.RatioBaseline(0.046, -0.2)
s.Re_floor = 1.0e4
cb.validate_rib_ratio_set(s)
r = cb.evaluate_rib_ratio(s, g, Re=30_000, Pr=0.7)
r.Nu, r.dNu_dRe, r.f, r.df_dRe, r.below_floor
```

Shipped: `cb.taslim_spring_1987(aspect_ratio_taslim, e_D)`, one set per
two-side configuration Taslim & Spring (1987) tested -- AR 0.5 (e/D 0.125,
0.250), 1.0 (0.083, 0.167), 3.5 (0.053, 0.107, 0.161), transverse ribs,
p/e 10; anything else raises. Their AR is height/width, so `W/H = 1/AR`.
Friction is a constant Fanning `f` over its own measured range
(`Re_floor_f`, which may differ from Nu's `Re_floor`):

```python
t = cb.taslim_spring_1987(1.0, 0.083)
g1 = cb.RibGeometry(e_D=0.083, p_e=10.0, W_H=1.0, alpha_deg=90.0)
cb.evaluate_rib_ratio(t, g1, Re=50_000, Pr=0.7).Nu
```

`RibRatioOptions.below_floor` picks the handover (`Gnielinski`,
`SourceBaseline`, or `User` with your own `user_Nu0`). It only acts below the
fitted range: inside it, Nu is the paper's whatever you choose -- matching a
specific rig across the range is what `Nu_multiplier` is for.

**Supply your own.** Real hardware needs it -- no published correlation is
precise enough for a specific rig:

```python
mine = cb.han_1988_orthogonal()
mine.name = "rig_3"
mine.source = "measured 2026-09"
mine.provenance = cb.RibProvenance.User   # not Extracted -- the claim differs
mine.C_G = 4.1
cb.validate_rib_set(mine)                 # rejects malformed sets, loudly
```

Three things worth knowing:

- **The normaliser is data.** `RibTerm(exponent, reference)` divides by
  `reference` before applying the exponent. `3.2 (p/e/10)^0.35` and
  `1.4294 (p/e)^0.35` are the same function; using one constant under the
  other convention is wrong by `10^0.35 = 2.24x`, uniformly, which never looks
  like a trend.
- **Validity is advisory.** `r.extrapolated` reports; nothing refuses. A band
  belongs to the source's rig, not to yours.
- **`St_r` is the ribbed side.** Combining it with the smooth walls is the
  caller's job, and depends on how many walls are ribbed.

Bad *parameters* raise from `validate_rib_set`. Bad *operating points* never
raise: reverse flow, zero flow and extreme values are guarded smoothly, because
a solver probes states that are not physical and a throw inside a residual
kills the solve.

### Ribbed Channels

```python
from combaero.network import ChannelElement, ConvectiveSurface, RibbedModel

surface = ConvectiveSurface(
    area=0.1,
    model=RibbedModel(
        e_D=0.06, p_e=10.0, alpha_deg=90.0,
        W_H=1.0,                  # aspect ratio: a Dh does not determine it
        n_ribbed_walls=2,         # 2 = two OPPOSITE walls, as Han measured
    ),
)
```

`n_ribbed_walls` accepts 1, 2 or 4. Three is rejected -- it has no unambiguous
geometry.

**Two asymmetries, both from the source rather than convenience:**

- **Friction needs no wall weighting.** The correlation's `f` is already the
  four-sided channel value, so the element uses it directly. Nothing multiplies
  pipe friction.
- **Heat transfer does.** The correlation gives the ribbed side; the smooth
  walls come from the base correlation and the channel average is the
  area-weighted combination. The result exposes `h_ribbed` and `h_smooth`
  separately, because the average hides a modelling choice worth seeing.

**A documented gap.** Plain smooth walls give `h_s/h_r` around 0.42 where Han's
own channel average implies 0.70 -- ribs enhance the adjacent smooth wall by
10-50% too, which no correlation here covers. The channel average is therefore
about **20% below** Han's measurement for the two-ribbed-wall square case.

```python
model = RibbedModel(..., smooth_wall_Nu_multiplier=1.67)   # reproduces Han
```

That knob exists because the ribbed side already has a better one: change `C_G`
on the parameter set, which records what you changed and why. The smooth walls
come from Gnielinski and have no set of their own.

### Pin-Fin Channels

`PinFinModel` puts a pin-fin array on a `ChannelElement`'s `ConvectiveSurface`
(#335). It uses the pin-fin sets below: heat transfer from `nu_set`
(default `metzger_1986_staggered_nu()`), the drop from `f_set` (default
`metzger_1982_staggered_friction()`), optionally transferred by a
`modifier`.

```python
from combaero.network import ChannelElement, ConvectiveSurface, PinFinModel

pins = PinFinModel(pin_diameter=0.005, S_D=2.5, X_D=2.5, H_D=1.0, N_rows=10,
                   k_pin=20.0)     # W/(m K), for the fin efficiency
ch = ChannelElement("te", "A", "B", length=0.125, diameter=0.009,
                    surface=ConvectiveSurface(area=0.01, model=pins))
```

- **The element's flow area is the unobstructed channel.** Vmax follows
  exactly from the array geometry (`A_min/A_frontal`), and
  `Re_D = rho Vmax D / mu`.
- **`ConvectiveSurface.area` is the endwall's base (planform) area.** The
  returned `h` is referenced to it and includes the fin efficiency of
  pins fed from each wall to mid-height:
  `h_eff = h (A_endwall_exposed + eta_fin A_pin) / A_base`. `h_array` and
  `eta_fin` on the result show both halves.
- **The drop is the correlation's own per-row form,**
  `dP = 2 rho Vmax^2 N f(Re_D)`, with an analytic Jacobian in mass flow,
  upstream T and P (through rho and mu). The code removed in 0.7.0 used
  `rho V^2/2` here, a factor of 4 low.
- **Inline arrays** default to `chyu_1990_friction(Inline)`, the only
  inline friction source in hand (one geometry: H/D 1, S/D = X/D = 2.5).
  Staggered friction is never substituted: a staggered `f_set` on an inline
  array raises. Heat transfer defaults to Chyu's inline set, or use
  `modifier=cb.chyu_1998_inline_over_staggered()` to transfer a staggered
  set.
- **Known gaps:**
  - `T_aw` is borrowed from the smooth correlation, as the impingement
    paths do.
  - `dh_dT` is the smooth result's relative property sensitivity scaled to
    `h`; the pin correlations expose no temperature derivative of their
    own.

### Impingement-Cooled Walls (jet-plate arrays)

A jet plate fed from a plenum, impinging on a target wall, the spent air
leaving through the gap as crossflow. **This is one configuration**, and
`ImpingementArray` builds it as one object (#465):

```python
from combaero.network import ImpingementArray, WallLayer

arr = ImpingementArray("ia", n_rows=10, d_jet=0.00254, xn_d=5.0, yn_d=4.0, z_d=2.0,
                       span=0.122, plate_thickness=0.00254)
arr.add_to(net, supply="supply", exit="exit")            # rows ia__p{i}, ia__x{i}, ia__c{i}
arr.add_wall(net, "w", hot_element="hot", layers=[WallLayer(0.001, 20.0)])
res = NetworkSolver(net).solve()
arr.summarize(res["__element_diag__"])   # m_dot, dP, Re_j/Gc_Gj ranges, rows_Nu, ...
```

`add_wall` puts one `ThermalWall` per row on that row's target footprint; the
hot side is evaluated once per row at its own state, so a hot side that changes
along the array needs its own segments. Heating the crossflow lowers its
density, which raises the jets' momentum cost and the supply non-uniformity
beyond Florschuetz's isothermal model (1.88 vs 1.61 on a 10-row (5,4,2) array
against 1200 K gas) -- the direction is physical; the isothermal source cannot
score it.

A bypass crossflow (#467) and spent air leaving through the target (#468) are
different configurations, with their own elements to come. The building
blocks below are public for hand-built networks, one pair per spanwise row:

```python
from combaero.network import ImpingementCrossflowElement, ImpingementPlateElement

d, span = 0.00254, 0.122
for i in range(1, 11):
    net.add_element(ImpingementPlateElement(
        f"p{i}", "supply", f"c{i}", d_jet=d, xn_d=5.0, yn_d=4.0, z_d=2.0,
        span=span, plate_thickness=d, row=i,   # row: diagnostics only
    ))
    net.add_element(ImpingementCrossflowElement(
        f"x{i}", f"c{i}", f"c{i+1}" if i < 10 else "exit",
        length=5.0 * d, height=2.0 * d, span=span,
    ))
```

```
supply plenum ---+--------------+--------------+
                 |              |              |
            [Plate 1]      [Plate 2]      [Plate 3]
                 |              |              |
crossflow:      c1 --[Cross]-- c2 --[Cross]-- c3 --[Cross]-- exit
```

**`ImpingementPlateElement`** is an orifice (`n_holes = round(span/(yn d))`
holes in parallel) whose flow is the row's jet flow, plus Florschuetz,
Truman and Metzger's (1981) heat transfer on the target footprint
`n_holes xn yn d^2`. Its **Gc/Gj comes from the network**, not the
uniform-supply closed form:

    Gc/Gj = (m_c / m_j) (pi/4) / ((yn/d)(z/d))

where `m_j` is the plate's own flow and `m_c` every other inflow to its
`to_node`, i.e. the crossflow approaching the row. A non-uniform supply, a row
of different geometry or a bleed therefore shows up in the heat transfer.
`h` depends on a neighbour's flow, so the element returns `dh_dsources` and
the solver relays it into the Jacobian (the `network_flow_inputs()` opt-in).
`T_aw` is the supply plenum temperature, Florschuetz's own reference.
The physics is C++'s `cb.jet_row_heat_transfer(set, cb.JetRowGeometry(...),
m_jet, m_crossflow, mu, k, Pr)`, which returns `h` and its analytic
derivatives with respect to both flows; the element only decides which flows
are crossflow.
- Hole Cd: default `'fixed'` at Florschuetz's 0.79 (measured 0.73-0.85);
  alternatives `'IdelchikThick'`, `'Lichtarowicz'` (long hole, l/d 2-10)
  and `'McGreehanSchotsch'` (supply-side `U1/Vi = 0`, never Gc/Gj).
  `Cd_in_range` reports whether the plate is inside the chosen source's range.
- Diagnostics: `n_holes`, `Re_j`, `Gc_Gj`, `Nu`, `htc`, `T_aw`,
  `surface_extrapolated`, and `Gc_Gj_closed_form` (Eq. 8) when `row` is given.
- Not modelled: the temperature sensitivity of the correlation's properties
  (`dh_dT = 0`) and the crossflow's own temperature.

**`ImpingementCrossflowElement`** is a `ChannelElement` (area `height *
span`) plus the momentum the jets cost: they enter with no streamwise
momentum, so each merge drops the static pressure by
`(m_b|m_b| - m_a|m_a|)/(rho A^2)`. The drop is split half-and-half between
the segments either side, so each crossflow node holds the static pressure at
its row centre. **Without this term every row sees the same pressure
difference and the supply stays uniform.** That puts Gc/Gj at twice Eq. 8 on
Florschuetz's strongest-crossflow geometry. The term is C++'s
`cb.station_half_drop` with `kappa = 0`, each station at the density of its
own node (see `CrossflowSegmentElement` below).

#### Crossflow segments: merge and bleed stations (#471)

`ImpingementCrossflowElement` is one case of `CrossflowSegmentElement`, a
channel segment whose end nodes are side-stream STATIONS. Over a station, with
the side stream's axial velocity `kappa * u_arriving`:

    dP_static = (m_b|m_b| - m_a|m_a| + kappa m_a (m_a - m_b)) / (rho A^2)

| `kappa` | station | source |
|---|---|---|
| `cb.STATION_KAPPA_MERGE_NORMAL` = 0 | jets merging normally (impingement) | Florschuetz's P + G^2/rho = const |
| `cb.STATION_KAPPA_BLEED_BASSETT` = 0.75 | bleed through the wall (effusion holes) | Bassett et al. (2001) K2/K5, Eq. 15, reproduced for every split |

```python
from combaero.network import CrossflowSegmentElement

seg = CrossflowSegmentElement("s1", "b1", "b2", length=0.02, area=4e-3, Dh=0.01,
                              from_kappa=cb.STATION_KAPPA_BLEED_BASSETT, prev_seg="s0",
                              to_kappa=cb.STATION_KAPPA_BLEED_BASSETT, next_seg="s2")
first = CrossflowSegmentElement("s0", "plenum", "b1", length=0.01, area=4e-3,
                                entry_K=0.5,          # reservoir entry, sharp
                                to_kappa=cb.STATION_KAPPA_BLEED_BASSETT, next_seg="s1")
```

- **Centred.** A segment carries half of the station at each end, each at its
  node's own density, so every node holds its station's mid static pressure.
  Station nodes must be `PlenumNode`s (Pt = P).
- **Ends are explicit.** `prev_seg`/`next_seg` name the neighbour segments;
  `entry_K` makes `from_node` a reservoir; a segment into a plenum with no
  `to_kappa` dumps its dynamic head (static continuity).
- **Bleed scored on Bassett's own measured K5** (Fig. 7c, 43 points): MAE
  0.039 in K overall, 0.021 at q >= 0.6 (the few-percent bleed of an effusion
  station). Bassett measured psi = 1-3; a station bleeding through holes far
  smaller than the duct is an extrapolation in psi, which K5 is independent
  of by derivation.

**One configuration only.** These elements cover spent air leaving down the
crossflow channel. A hand-wired network describing a different
configuration is refused at set-up, naming it:
- a bypass or initial crossflow feeding the chain (#467);
- spent air leaving through the target, impingement-effusion or film (#468);
- two plates merging at one node.

**Measured against the source** (`python/tests/test_impingement_plate_validation.py`):
- **Flow model:** over the 27 Fig. 6 geometries x 10 rows, the network
  reproduces Florschuetz's own 1D flow model (Eq. 8, and the cosh jet
  distribution behind it) to a max 3.7% in Gc/Gj (bias -0.5%) and 3.8% in
  Gj/Gj_mean.
- **Fig. 6:** the chain's Nu/Nu1 scores +4.0% bias and 6.9% MAE over 242
  points. The correlation at the paper's own abscissae scores +3.7% / 6.8%.
- **Row 1 vs Fig. 5:** absolute Nu1 matches as in #461.

### Impingement Channels

```python
from combaero.network import ConvectiveSurface, ImpingementModel, SingleJetImpingementModel

surface = ConvectiveSurface(
    area=0.01,   # this ROW's own target-plate footprint
    model=ImpingementModel(
        d_jet=0.002, xn_d=8.0, yn_d=6.0, z_d=2.0,
        row=3,       # 1-indexed, counting from upstream
    ),
)
```

For a jet plate fed from a plenum use the elements above: they carry the
row's own jet flow and read Gc/Gj from the network. `ImpingementModel` on a
`ChannelElement` remains for a row whose flow you impose yourself.

**One element models one spanwise row**, not a whole array -- there is no
single "channel Nu" for a jet array the way there is for a smooth or ribbed
duct, since downstream rows see progressively more crossflow than upstream
ones. Model a real array by chaining `n_rows` separate elements, each with
its own `row`; row 1 sees zero crossflow by definition (`Gc/Gj = 0`,
Florschuetz's own `Nu1`).

**The element's mass flow is this row's jet flow; `area` counts the
holes.** `area` is the row's target-plate footprint, so the row has
`n = area / (xn_d d_jet * yn_d d_jet)` holes, each carrying `m_dot / n`, and
`Re_j = 4 (m_dot / n) / (pi d_jet mu)` -- Florschuetz's jet mass velocity on
the hole area. Until #460 the jet flow was taken as `rho v * area`, which
made Re_j scale with the target area: driven with the paper's own flow, the
element then missed Florschuetz's Fig. 5 by +268% on average; it now matches
it to 5.2%, like the correlation. Whether every row gets the same total flow (a uniform-supply
approximation) or a row-dependent one (Florschuetz's own Eq. 7, deliberately
not implemented -- see `impingement_correlation.h`'s module comment) is a
choice for whoever assembles the chain, not something this element decides.

**Default convective area.** A `ChannelElement` whose surface was given no
`area` gets one that follows the model (#462): an impingement row gets its
target footprint matched to the channel's crossflow cross-section,
`A_flow * (xn/d) / (z/d)`, so the hole count follows from the channel; a single
jet gets Goldstein's averaging disc, `pi (R_D d_jet)^2`; every other surface
keeps the wetted wall `pi D L`. The channel diagnostics report the surface's
own `Re_surface` and `surface_extrapolated`.

**Pressure drop is not modelled here.** Model the jet plate's own orifice
loss with a proper `OrificeElement` upstream, using a real
discharge-coefficient correlation -- duplicating it inside this heat-transfer
model would risk double-counting it. `f`/`dP` reported by this element are
the crossflow's own plain smooth-duct values (Han's book gives no
impingement-specific friction correlation, unlike ribs' own `R`/`f`
relationship), borrowed to compute a `T_aw` Eq. 4.2/4.3's own recovery-factor
model would otherwise have to supply.

**A single free jet** (Goldstein, Behbahani and Heppelmann 1986) is a
separate model, since there is no array and no crossflow:

```python
surface = ConvectiveSurface(
    area=0.001,   # this jet's own target patch
    model=SingleJetImpingementModel(d_jet=0.003, L_D=7.75, R_D=5.0),
)
```

`R_D` is a real modelling choice, not a detail: the correlation's `Nu` is
the AVERAGE over the disc of radius `R = R_D d_jet` (Han's Eq. 4.1), so the
convective area should be that disc.

Both models expose the same wall-coupling derivatives (`dh_dmdot`, `dh_dT`,
`dT_aw_dmdot`, `dT_aw_dT`) and an `extrapolated` flag (array only -- a single
jet has no stated validity range to be outside of), matching `RibbedModel`'s
contract.

### Jet Impingement Correlations (parametrised)

Two independent regimes, matching
`validation/cooling/extractions/han_impingement.md`'s own split, wired into
the network elements above (`ImpingementModel`, `SingleJetImpingementModel`).

**Single jet** (Goldstein, Behbahani and Heppelmann, 1986): one free round
jet on a flat plate, no crossflow, no array.

```python
import combaero as cb

s = cb.goldstein_1986_single_jet()
nu_q = cb.single_jet_impingement_nu(
    s, cb.ImpingementThermalBC.ConstantHeatFlux, Re=25_000, L_D=7.75, R_D=5.0)
nu_t = cb.single_jet_impingement_nu(
    s, cb.ImpingementThermalBC.ConstantWallTemperature, Re=25_000, L_D=7.75, R_D=5.0)
nu_q, nu_t                            # 59.93, 55.71 -- Han's own check point (60, 56)
```

`L_D = 7.75` is the optimum spacing -- structural, not a separately fitted
fact: `A` itself is the numerator's value there, for any `Re` or `R_D`.

**Jet array with crossflow** (Florschuetz, Truman and Metzger, 1981): a
staggered or inline array fed from a common plenum, where every row's Nu is
degraded by the crossflow accumulated from every row upstream of it.

```python
inline = cb.florschuetz_1981_inline()

result = cb.jet_array_impingement_nu(
    inline, Re_j=10_000, Gc_Gj=0.3, Pr=0.7, xn_d=10.0, yn_d=6.0, z_d=2.0)
result.Nu, result.extrapolated        # 27.93, False
```

`Re_j` and `Gc_Gj` are that ROW's own local values, not an array mean --
Florschuetz correlates "the individual spanwise row jet Reynolds number".
`Gc_Gj` at a row has its own closed form (the paper's own Eq. 8), needing
only geometry and a jet-plate discharge coefficient:

```python
cb.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=cb.FLORSCHUETZ_1981_DEFAULT_CD, row=1)   # 0.0 -- Nu1's own definition
cb.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=cb.FLORSCHUETZ_1981_DEFAULT_CD, row=10)  # 0.404
```

It depends on `(yn/d)(z/d)` only, not `xn/d` -- the source states the flow
distribution is independent of streamwise hole spacing and hole pattern.
`FLORSCHUETZ_1981_DEFAULT_CD` (0.79) is the paper's own recommended default
absent a measured value, and remains the default. To predict `C_D` from the
plate's own geometry instead, pass `mcgreehan_schotsch_1988_cd` (see above);
for a plenum-fed plate leave its `U1_over_Vi` at 0, since Florschuetz's
`Gc/Gj` is discharge-side crossflow and must not be fed to it:

```python
C_D = cb.mcgreehan_schotsch_1988_cd(Re=1e4, r_over_d=0.0, L_over_d=t_over_d)
cb.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=C_D, row=10)
```

At `t/d = 1` the two agree to 0.6%, so this changes little; at `t/d = 0.5` the
correlation gives 0.688 and the predicted `Gc/Gj` rises by up to 12% at
downstream rows (about 1% in `Nu`). The ISO 5167 `Cd_sharp_thin_plate` family
in `orifice.h` remains inapplicable here -- it models a hole in a pipe run.

**Two named sets**, chosen explicitly by hole pattern -- their coefficients
differ, not just their validity:

```python
inl = cb.florschuetz_1981_inline()      # xn/d 5-15
stg = cb.florschuetz_1981_staggered()   # xn/d 5-10, genuinely tighter
```

Bad *parameters* raise from `validate_single_jet_set`/`validate_jet_array_set`.
Bad *operating points* never raise: reverse flow and out-of-range crossflow
ratios are guarded smoothly, matching the ribs' policy above.

Leading-edge/curved-surface impingement (Section 4.1.4) and the simpler,
looser forms (Eq. 4.6 Kercher-Tabakoff -- graphical, not closed-form anyway;
Eq. 4.7/4.8 -- Florschuetz's own less-tight alternate) are out of scope,
deferred per the extraction's I3.

### Pin-Fin Correlations (parametrised)

Pin-fin arrays (#335) as provenanced sets in each source's own form, with
exact geometry-only conversion to one canonical basis: `Re_D` on the pin
diameter and the velocity at the minimum flow area, `Nu_D = h D / k`, and
`f = dP / (2 rho Vmax^2 N)` (per row, so `dP = 2 rho Vmax^2 N f`). Heat
transfer and friction are separate sets, so any Nu set can pair with any
friction set.

```python
import combaero as cb

g = cb.PinFinGeometry(S_D=2.5, X_D=2.5, H_D=1.0, N_rows=10)   # staggered
nu = cb.evaluate_pin_fin_nu(cb.metzger_1986_staggered_nu(), g, Re_D=2.0e4, Pr=0.7)
f = cb.evaluate_pin_fin_friction(cb.metzger_1982_staggered_friction(), g, Re_D=2.0e4)
nu.Nu, nu.dNu_dRe, nu.extrapolated     # 91.78, 0.00317, False
f.f, f.extrapolated                    # 0.0755, False

# VanFossen works in its own D' basis; the result comes back canonical.
cb.evaluate_pin_fin_nu(cb.vanfossen_1982_staggered_nu(), g, 2.0e4, 0.7).Nu   # 90.91
```

The second call is flagged `extrapolated`: VanFossen's arrays had 4 rows, not
10. A set outside its validity box still evaluates; it never refuses.

| set | what | box |
|---|---|---|
| `metzger_1986_staggered_nu()` (default) | `0.135 Re_D^0.69 (X/D)^-0.34`, pin + endwall | H/D <= 3, 2 <= S/D <= 4, 1.5 <= X/D <= 5, Re 1e3-1e5, 10 rows |
| `metzger_1982_staggered_friction()` (default) | `0.317 Re^-0.132` / `1.76 Re^-0.318`, C1 blend at 1e4 | 0.5 <= H/D <= 6, 2 <= S/D <= 4 |
| `vanfossen_1982_staggered_nu()` | `0.153 Re_D'^0.685` on D' = 4V/A_t | H/D 0.5-2, 4 rows |
| `damerow_1972_staggered_friction()` | `2.06 (S/D)^-1.1 Re^-0.16`, per (N-1) rows | S/D 4.24-7.07 |
| `chyu_1998_nu(arrangement, surface)` | Han Table 4.7, staggered or inline, pin / endwall / total | S/D = X/D = 2.5, H/D = 1 |
| `chyu_1998_inline_over_staggered()` | inline/staggered Nu ratio, no friction ratio | Chyu's geometry |
| `chyu_1990_nu(arrangement, fillet=False)` | Chyu 1990 Table 2, pin surface, straight or fillet pins | S/D = X/D = 2.5, H/D = 1, Re 5e3-3e4 |
| `chyu_1990_friction(arrangement, fillet=False)` | Chyu 1990 Fig. 6, digitised and **Fitted**; inline straight = 0.1693/4, Re-independent | same, Re about 9e3-2.2e4 |
| `chyu_1990_fillet_over_straight(arrangement)` | fillet/straight Nu ratio (staggered 0.78-0.92), no friction ratio | same |

**Inline friction** comes from one rig at one geometry: Chyu (1990) Fig. 6,
digitised (points in `validation/cooling/data/chyu1990/`). Chyu's staggered
straight-pin friction sits 16-25% below Metzger's at the same geometry, a
disagreement between the labs that is reported, not reconciled.

**Fin efficiency** of conducting pins (each wall feeds the pin to mid-height):

```python
frac = cb.pin_fin_area_fractions(g)                 # per wall, over base area
e = cb.pin_fin_array_efficiency(h=2000.0, k_pin=20.0, D=1e-3, H=1e-3,
                                A_f_over_A_t=frac.pin_over_total)
e.eta_fin, e.eta_t, e.deta_t_dh                     # 0.968, 0.993, analytic
```

**Supply your own** by copying the default and saying so:

```python
mine = cb.metzger_1986_staggered_nu()
mine.name, mine.source = "rig_7", "measured 2026-10"
mine.provenance = cb.RibProvenance.User
mine.C = 0.150
cb.validate_pin_fin_nu_set(mine)
```

To use a set in a network, see **Pin-Fin Channels** (`PinFinModel`).

Not yet carried, and declared rather than assumed: row-count correction (a
geometry with a different `N_rows` is flagged), channel convergence, long pins
(H/D > 3), and fillets beyond Chyu's single geometry.

### Enhanced Cooling Surfaces

Removed in 0.7.0. The pin-fin, dimple, rib and impingement correlations could
not be traced to their cited sources -- the rib friction multiplier was 4-5x
below the only rib datum in the repository. Provenanced rib, jet
impingement and pin-fin correlations are re-added above (issues #334, #337,
#335); see issue #339 for the rebuild.

`channel_smooth` and the base convective correlations (Gnielinski,
Dittus-Boelter, Sieder-Tate, Petukhov) are unaffected, as are the user-set
`Nu_multiplier` and `f_multiplier` knobs on `ConvectiveSurface`.

---

## Advanced Thermodynamics

### Complete State Functions
```python
# Air properties bundle
air = cb.air_properties(T=300, P=101325, humidity=0.0)
print(f"Density: {air.density} kg/m³")
print(f"Viscosity: {air.viscosity} Pa·s")
print(f"Thermal conductivity: {air.thermal_conductivity} W/(m·K)")

# Complete thermodynamic state
thermo = cb.thermo_state(T=300, P=101325, X=air)
complete = cb.complete_state(T=300, P=101325, X=air)
```

### Derivatives
```python
# Temperature derivatives
dh_dT = cb.dh_dT(T=300, X=air)
ds_dT = cb.ds_dT(T=300, X=air)
dcp_dT = cb.dcp_dT(T=300, X=air)
dg_dT = cb.dg_over_RT_dT(T=300, X=air)
```

### Advanced Inverse Solvers
```python
# Molar-basis inverses (targets in J/mol or J/(mol*K))
T_from_h   = cb.calc_T_from_h(h_target=3e5, X=air)            # J/mol -> K
T_from_s   = cb.calc_T_from_s(s_target=7000, P=101325, X=air) # J/(mol*K), Pa -> K
T_from_cp  = cb.calc_T_from_cp(cp_target=29, X=air)           # J/(mol*K) -> K
T_from_u   = cb.calc_T_from_u(u_target=2.2e5, X=air)          # J/mol -> K

# Mass-specific inverses (targets in J/kg or J/(kg*K))
T_from_h_m = cb.calc_T_from_h_mass(h_mass_target=3e5, X=air)
T_from_s_m = cb.calc_T_from_s_mass(s_mass_target=7000, P=101325, X=air)
T_from_u_m = cb.calc_T_from_u_mass(u_mass_target=2.2e5, X=air)

# Flash calculations — solve T from two independent properties
T_from_sv  = cb.calc_T_from_sv_mass(s_mass_target=7000, v_mass_target=0.8, X=air)
T_from_sh  = cb.calc_T_from_sh_mass(s_mass_target=7000, h_mass_target=3e5, X=air)
```

### Flow Analysis
```python
# Bulk velocity and kinetic energy
v = cb.bulk_velocity(m_dot=0.5, rho=1.2, area=1e-3)   # m/s
KE = cb.kinetic_energy(v=v)                             # J/kg

# Dimensionless numbers
# From full mixture state (recomputes rho, mu internally):
Re = cb.reynolds(T=300, P=101325, X=air, V=v, L=0.05)
# From pre-computed state properties (use when CompleteState is available):
Re = cb.reynolds_from_state(rho=1.2, v=v, L=0.05, mu=1.8e-5)

# Stagnation conditions from static state and velocity
Tt, Pt = cb.stagnation_from_static(T=300, P=101325, v=v, X=air)  # (K, Pa)
```

---

## Psychrometrics (Humid Air)

```python
air = cb.HumidAir()
air.set_TP_RH(300.0, 101325.0, 0.5)

print(f"Humidity ratio: {air.humidity_ratio} kg/kg")
print(f"Dewpoint: {air.dewpoint} K")
print(f"Underlying state P: {air.state.P} Pa")
```
