# CombAero API Reference

This document provides detailed API reference for LLMs and developers.

## Units

All functions use a consistent SI-based unit system. See **[UNITS.md](UNITS.md)** for the complete reference.

**Key conventions:**
- Temperature: K (Kelvin)
- Pressure: Pa (Pascal)
- Mole fractions: mol/mol (dimensionless)
- Mass fractions: kg/kg (dimensionless)
- Thermodynamic properties: molar basis (J/mol, J/(mol·K))
- Compressible flow outputs: mass basis (J/kg, J/(kg·K))

**Programmatic unit queries (C++ and Python):**
```cpp
#include "units.h"
auto u = combaero::units::get_units("density");  // {input: "T: K, P: Pa, X: mol/mol", output: "kg/m^3"}
```
```python
import combaero
combaero.get_units("density")  # UnitInfo(input='T: K, P: Pa, X: mol/mol', output='kg/m^3')
```

## Species

The library uses a fixed set of 14 species. Use `species_index_from_name()` for robust lookup:

| Index | Formula | Common Name |
|-------|---------|-------------|
| 0 | N2 | Nitrogen |
| 1 | O2 | Oxygen |
| 2 | AR | Argon |
| 3 | CO2 | Carbon Dioxide |
| 4 | H2O | Water |
| 5 | CH4 | Methane |
| 6 | C2H6 | Ethane |
| 7 | C3H8 | Propane |
| 8 | IC4H10 | Isobutane |
| 9 | NC5H12 | n-Pentane |
| 10 | NC6H14 | n-Hexane |
| 11 | NC7H16 | n-Heptane |
| 12 | CO | Carbon Monoxide |
| 13 | H2 | Hydrogen |

### Thermodynamic Data

The library uses **NASA-9 polynomials** for thermodynamic properties:
- Valid temperature range: 200-6000 K (some species to 20000 K)
- 1-3 temperature intervals per species
- Higher accuracy than NASA-7, especially at high temperatures

### Species Data Files

- **`thermo_transport_data.h`**: Contains all species data including `species_names` vector, `species_index` unordered map, `molar_masses`, `nasa_coeffs`, `transport_props`, and `molecular_structures`. This file is auto-generated.
- **`common_names.h`**: Maps between formulas and human-readable names (in `combaero::` namespace):
  - `formula_to_name` map: `"CH4"` → `"Methane"`
  - `name_to_formula` map: `"Methane"` → `"CH4"`
  - `common_name(formula)` function: lookup with error handling
  - `formula(name)` function: inverse lookup with error handling
- **`thermo_data_generator/`**: Python utility that generates `thermo_transport_data.h` from NASA CEA data and Cantera YAML mechanisms.

**Note**: Always use `species_index_from_name("CH4")` rather than hardcoded indices for forward compatibility.

## Inverse Combustion Solvers

These functions find the fuel or oxidizer mass flow rate to achieve a target property in burned products. All use complete combustion to CO₂/H₂O (no WGS).

### Find Fuel Stream (oxidizer mdot fixed)

```cpp
// Target adiabatic flame temperature
Stream set_fuel_stream_for_Tad(
    double T_ad_target,           // Target Tad [K], must be > oxidizer T
    const Stream& fuel,           // Fuel stream (T, P, X set; mdot ignored)
    const Stream& oxidizer,       // Oxidizer stream (with mdot set)
    double tol = 1.0,             // Temperature tolerance [K]
    std::size_t max_iter = 100,   // Maximum bisection iterations
    bool lean = true,             // true: lean side (O2 excess), false: rich side
    double phi_max = 10.0         // Max equivalence ratio for search range
);

// Target O2 mole fraction (wet basis) - lean only
Stream set_fuel_stream_for_O2(double X_O2_target, const Stream& fuel, const Stream& oxidizer,
                               double tol = 1e-6, std::size_t max_iter = 100);

// Target O2 mole fraction (dry basis) - lean only
Stream set_fuel_stream_for_O2_dry(double X_O2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                   double tol = 1e-6, std::size_t max_iter = 100);

// Target CO2 mole fraction (wet basis) - lean only
Stream set_fuel_stream_for_CO2(double X_CO2_target, const Stream& fuel, const Stream& oxidizer,
                                double tol = 1e-6, std::size_t max_iter = 100);

// Target CO2 mole fraction (dry basis) - lean only
Stream set_fuel_stream_for_CO2_dry(double X_CO2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                    double tol = 1e-6, std::size_t max_iter = 100);
```

### Find Oxidizer Stream (fuel mdot fixed)

```cpp
// Target adiabatic flame temperature
Stream set_oxidizer_stream_for_Tad(
    double T_ad_target,           // Target Tad [K], must be > fuel T
    const Stream& fuel,           // Fuel stream (with mdot set)
    const Stream& oxidizer,       // Oxidizer stream (T, P, X set; mdot ignored)
    double tol = 1.0,             // Temperature tolerance [K]
    std::size_t max_iter = 100,   // Maximum bisection iterations
    bool lean = true,             // true: lean side (O2 excess), false: rich side
    double phi_max = 10.0         // Max equivalence ratio for search range
);

// Target O2 mole fraction (wet basis) - lean only
Stream set_oxidizer_stream_for_O2(double X_O2_target, const Stream& fuel, const Stream& oxidizer,
                                   double tol = 1e-6, std::size_t max_iter = 100);

// Target O2 mole fraction (dry basis) - lean only
Stream set_oxidizer_stream_for_O2_dry(double X_O2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                       double tol = 1e-6, std::size_t max_iter = 100);

// Target CO2 mole fraction (wet basis) - lean only
Stream set_oxidizer_stream_for_CO2(double X_CO2_target, const Stream& fuel, const Stream& oxidizer,
                                    double tol = 1e-6, std::size_t max_iter = 100);

// Target CO2 mole fraction (dry basis) - lean only
Stream set_oxidizer_stream_for_CO2_dry(double X_CO2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                        double tol = 1e-6, std::size_t max_iter = 100);
```

### Key Concepts

- **Lean combustion** (φ < 1): O₂ excess in products. Tad increases with fuel up to stoichiometric.
- **Rich combustion** (φ > 1): Fuel excess, no O₂ in products. Tad decreases with fuel beyond stoichiometric.
- **Tad ambiguity**: The same Tad can be achieved on both lean and rich sides. Use `lean` parameter to select.
- **O₂/CO₂ solvers**: Only support lean combustion (O₂ > 0 in products).
- **Dry basis**: Water vapor removed from products before computing mole fraction (for emission sampling).

### Python Usage

```python
import combaero as cb
import numpy as np

# Create streams
n = len(cb.standard_dry_air_composition())
X_CH4 = np.zeros(n)
X_CH4[6] = 1.0  # CH4 at index 6

fuel = cb.Stream()
fuel.T = 300.0
fuel.P = 101325.0
fuel.X = X_CH4

air = cb.Stream()
air.T = 300.0
air.P = 101325.0
air.X = cb.standard_dry_air_composition()
air.mdot = 10.0

# Find fuel for target Tad (lean side, default)
fuel_lean = cb.set_fuel_stream_for_Tad(1500.0, fuel, air)

# Find fuel for target Tad (rich side)
fuel_rich = cb.set_fuel_stream_for_Tad(1500.0, fuel, air, lean=False)

# Find fuel for target dry O2 (emission sampling)
fuel_o2 = cb.set_fuel_stream_for_O2_dry(0.12, fuel, air)

# Convert burned products to dry basis
mixed = cb.mix([fuel_lean, air])
burned = cb.complete_combustion(mixed.T, mixed.X, mixed.P)
X_dry = cb.convert_to_dry_fractions(burned.X)
```

## Combustion Functions

```cpp
// Complete combustion to CO2/H2O (adiabatic)
State complete_combustion(const State& in);

// Complete combustion to CO2/H2O (isothermal)
State complete_combustion_isothermal(const State& in);

// Get burned composition without solving temperature
std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double>& X);

// Set fuel stream for target equivalence ratio
Stream set_fuel_stream_for_phi(double phi, const Stream& fuel, const Stream& oxidizer);
```

## Equilibrium Functions

```cpp
// Water-gas shift equilibrium (isothermal)
State wgs_equilibrium(const State& in);

// Water-gas shift equilibrium (adiabatic)
State wgs_equilibrium_adiabatic(const State& in);

// Steam methane reforming + WGS equilibrium
State smr_wgs_equilibrium(const State& in);
State smr_wgs_equilibrium_adiabatic(const State& in);

// General reforming equilibrium (SMR + WGS + dry reforming)
State reforming_equilibrium(const State& in);
State reforming_equilibrium_adiabatic(const State& in);

// Combustion equilibrium (complete combustion + WGS)
State combustion_equilibrium(const State& in);
```

## Stream Mixing

```cpp
// Mix multiple streams with mass and enthalpy balance
Stream mix(const std::vector<Stream>& streams, double P_out = -1.0);
```

## Thermodynamic Properties

All functions take temperature [K] and mole fractions vector:

```cpp
double cp(double T, const std::vector<double>& X);      // J/(mol·K)
double h(double T, const std::vector<double>& X);       // J/mol
double s(double T, const std::vector<double>& X);       // J/(mol·K)
double cv(double T, const std::vector<double>& X);      // J/(mol·K)
double density(double T, double P, const std::vector<double>& X);  // kg/m³
double speed_of_sound(double T, const std::vector<double>& X);     // m/s
```

## Transport Properties

```cpp
double viscosity(double T, double P, const std::vector<double>& X);           // Pa·s
double thermal_conductivity(double T, double P, const std::vector<double>& X); // W/(m·K)
double prandtl(double T, double P, const std::vector<double>& X);             // dimensionless
double kinematic_viscosity(double T, double P, const std::vector<double>& X); // m²/s (ν = μ/ρ)
double thermal_diffusivity(double T, double P, const std::vector<double>& X); // m²/s (α = k/(ρ·cp))
double reynolds(double T, double P, const std::vector<double>& X, double V, double L);  // Re = ρVL/μ
double peclet(double T, double P, const std::vector<double>& X, double V, double L);    // Pe = VL/α
```

## Utility Functions

```cpp
// Mixture molecular weight
double mwmix(const std::vector<double>& X);  // [g/mol]

// Mole/mass fraction conversion
std::vector<double> mole_to_mass(const std::vector<double>& X);
std::vector<double> mass_to_mole(const std::vector<double>& Y);

// Remove water vapor and renormalize
std::vector<double> convert_to_dry_fractions(const std::vector<double>& X);

// Normalize fractions to sum to 1
std::vector<double> normalize_fractions(const std::vector<double>& X);

// Oxygen requirements
double oxygen_required_per_mol_fuel(const std::string& fuel);
double oxygen_required_per_kg_fuel(const std::string& fuel);
double oxygen_required_per_mol_mixture(const std::vector<double>& X);
double oxygen_required_per_kg_mixture(const std::vector<double>& X);

// Equivalence ratio
double equivalence_ratio_mole(const std::vector<double>& X_fuel, const std::vector<double>& X_ox,
                               const std::vector<double>& X_mix);
double equivalence_ratio_mass(const std::vector<double>& Y_fuel, const std::vector<double>& Y_ox,
                               const std::vector<double>& Y_mix);
```

## Humid Air

```cpp
// Standard dry air composition (N2, O2, Ar, CO2)
std::vector<double> standard_dry_air_composition();

// Humid air composition from relative humidity
std::vector<double> humid_air_composition(double T, double P, double RH);

// Dew point temperature
double dewpoint(double T, double P, const std::vector<double>& X);
```

## Compressible Flow

Compressible flow solvers for ideal gas with variable cp(T). Uses NASA polynomial fits for temperature-dependent properties.

### Isentropic Nozzle Flow

```cpp
// Forward problem: given geometry and pressures, find mass flow
CompressibleFlowSolution nozzle_flow(
    double T0, double P0,           // Stagnation conditions [K, Pa]
    double P_back,                  // Back pressure [Pa]
    double A_eff,                   // Effective flow area [m²] (= A * Cd)
    const std::vector<double>& X,   // Mole fractions
    double tol = 1e-8,
    std::size_t max_iter = 50
);

// Inverse problems: solve for unknown given mass flow
double solve_A_eff_from_mdot(double T0, double P0, double P_back, double mdot_target,
                              const std::vector<double>& X, ...);
double solve_P_back_from_mdot(double T0, double P0, double A_eff, double mdot_target,
                               const std::vector<double>& X, ...);
double solve_P0_from_mdot(double T0, double P_back, double A_eff, double mdot_target,
                           const std::vector<double>& X, ...);

// Utilities
double critical_pressure_ratio(double T0, double P0, const std::vector<double>& X, ...);
double mach_from_pressure_ratio(double T0, double P0, double P, const std::vector<double>& X, ...);
double mass_flux_isentropic(double T0, double P0, double P, const std::vector<double>& X, ...);
```

### Quasi-1D Nozzle Flow

Solves quasi-1D compressible flow through a nozzle with variable area A(x):
- Mass: ṁ = ρuA = const
- Energy: h + u²/2 = h₀ (adiabatic)
- Momentum: dp/dx + ρu·du/dx = 0

```cpp
// General form with area function A(x)
NozzleSolution nozzle_quasi1d(
    double T0, double P0, double P_exit,
    const AreaFunction& area_func,    // std::function<double(double)>
    double x_start, double x_end,
    const std::vector<double>& X,
    std::size_t n_stations = 100,
    double tol = 1e-8, std::size_t max_iter = 50
);

// With area as (x, A) pairs (linearly interpolated)
NozzleSolution nozzle_quasi1d(
    double T0, double P0, double P_exit,
    const std::vector<std::pair<double, double>>& area_profile,
    const std::vector<double>& X, ...);

// With polynomial area: A(x) = a₀ + a₁x + a₂x² + ...
NozzleSolution nozzle_quasi1d_poly(
    double T0, double P0, double P_exit,
    const std::vector<double>& area_coeffs,
    double x_start, double x_end,
    const std::vector<double>& X, ...);

// Converging-diverging nozzle with smooth cosine profile
NozzleSolution nozzle_cd(
    double T0, double P0, double P_exit,
    double A_inlet, double A_throat, double A_exit,
    double x_throat, double x_exit,
    const std::vector<double>& X, ...);
```

#### NozzleSolution Structure

```cpp
struct NozzleStation {
    double x, A, P, T, rho, u, M, h;  // Axial profile data
};

struct NozzleSolution {
    State inlet, outlet;              // Inlet/outlet states
    double mdot, h0, T0, P0;          // Flow parameters
    bool choked;                      // True if throat is sonic
    double x_throat, A_throat;        // Throat location and area
    std::vector<NozzleStation> profile; // Axial profile
};
```

### Fanno Flow (Adiabatic Pipe with Friction)

Solves compressible pipe flow with wall friction using RK4 integration:
- Mass: ρuA = const
- Energy: h + u²/2 = h₀ (adiabatic)
- Momentum: dp/dx = -f/(2D)·ρu²

```cpp
// Solve Fanno flow through a pipe segment
FannoSolution fanno_pipe(
    double T_in, double P_in,       // Inlet static conditions [K, Pa]
    double u_in,                    // Inlet velocity [m/s]
    double L, double D,             // Pipe length and diameter [m]
    double f,                       // Darcy friction factor [-]
    const std::vector<double>& X,   // Mole fractions
    std::size_t n_steps = 100,      // Integration steps
    bool store_profile = false      // Store axial profile
);

// Convenience overload with State
FannoSolution fanno_pipe(const State& inlet, double u_in, double L, double D, double f,
                          std::size_t n_steps = 100, bool store_profile = false);

// Find maximum pipe length before choking (L*)
double fanno_max_length(double T_in, double P_in, double u_in, double D, double f,
                         const std::vector<double>& X, double tol = 1e-6, ...);
```

### Result Structures

```cpp
struct CompressibleFlowSolution {
    State stagnation;    // Stagnation state (T0, P0, ...)
    State outlet;        // Outlet static state
    double v;            // Outlet velocity [m/s]
    double M;            // Outlet Mach number [-]
    double mdot;         // Mass flow rate [kg/s]
    bool choked;         // True if flow is choked (M = 1)
};

struct FannoStation {
    double x, P, T, rho, u, M, h, s;  // Axial profile data
};

struct FannoSolution {
    State inlet, outlet;              // Inlet/outlet states
    double mdot, h0, L, D, f;         // Flow parameters
    bool choked;                      // True if M reached 1
    double L_choke;                   // Length to choking [m]
    std::vector<FannoStation> profile; // Axial profile (if requested)
};
```

## Friction Factor Correlations

Darcy friction factor for turbulent pipe flow. Valid for Re > ~2300.

```cpp
// Haaland (1983) - explicit, ~2-3% accuracy vs Colebrook
// 1/√f = -1.8·log₁₀((ε/D/3.7)^1.11 + 6.9/Re)
double friction_haaland(double Re, double e_D);

// Serghides (1984) - explicit, <0.3% accuracy vs Colebrook
double friction_serghides(double Re, double e_D);

// Colebrook-White (1939) - implicit, reference standard
// 1/√f = -2·log₁₀(ε/D/3.7 + 2.51/(Re·√f))
// Solved iteratively with Haaland initial guess
double friction_colebrook(double Re, double e_D, double tol = 1e-10, int max_iter = 20);
```

### Parameters

- `Re`: Reynolds number [-] (must be > 0)
- `e_D`: Relative roughness ε/D [-] (must be ≥ 0)

### Usage Example

```cpp
double Re = 50000.0;
double e_D = 0.001;  // Relative roughness

double f_haaland = friction_haaland(Re, e_D);     // Fast explicit
double f_serghides = friction_serghides(Re, e_D); // More accurate explicit
double f_colebrook = friction_colebrook(Re, e_D); // Reference (iterative)

// Use in pressure drop: ΔP = f·(L/D)·(ρv²/2)
```

## Incompressible Flow (`incompressible.h`)

Analytical equations for incompressible (constant density) flow. These are valid for:
- Mach number < 0.3 (density change < 5%)
- Pressure ratio close to 1 (ΔP/P << 1)

For higher Mach numbers or large pressure ratios, use `compressible.h`.

### Bernoulli Equation

```cpp
// Solve for downstream pressure P2
// P1 + ½ρv1² + ρgh1 = P2 + ½ρv2² + ρgh2
double bernoulli_P2(double P1, double v1, double v2, double rho,
                    double dz = 0.0, double g = 9.80665);

// Solve for downstream velocity v2
double bernoulli_v2(double P1, double P2, double v1, double rho,
                    double dz = 0.0, double g = 9.80665);
```

### Orifice Flow

```cpp
// Mass flow rate: ṁ = Cd · A · √(2 · ρ · ΔP)
double orifice_mdot(double P1, double P2, double A, double Cd, double rho);

// Volumetric flow rate: Q = Cd · A · √(2 · ΔP / ρ)
double orifice_Q(double P1, double P2, double A, double Cd, double rho);

// Ideal velocity (no losses): v = √(2 · ΔP / ρ)
double orifice_velocity(double P1, double P2, double rho);

// Inverse: solve for area given mass flow
double orifice_area(double mdot, double P1, double P2, double Cd, double rho);

// Inverse: solve for pressure drop given mass flow
double orifice_dP(double mdot, double A, double Cd, double rho);
```

**Typical discharge coefficients (Cd):**
- Sharp-edge orifice: 0.60-0.65
- Rounded orifice: 0.95-0.98
- Nozzle: 0.95-0.99

### Pipe Flow (Darcy-Weisbach)

```cpp
// Pressure drop: ΔP = f · (L/D) · (ρ · v² / 2)
double pipe_dP(double v, double L, double D, double f, double rho);

// Pressure drop from mass flow rate
double pipe_dP_mdot(double mdot, double L, double D, double f, double rho);

// Velocity from mass flow: v = ṁ / (ρ · π · D² / 4)
double pipe_velocity(double mdot, double D, double rho);

// Mass flow from velocity: ṁ = ρ · v · π · D² / 4
double pipe_mdot(double v, double D, double rho);
```

### Hydraulic Utilities

```cpp
// Dynamic pressure (velocity head): q = ½ · ρ · v²
double dynamic_pressure(double v, double rho);

// Velocity from dynamic pressure: v = √(2 · q / ρ)
double velocity_from_q(double q, double rho);

// Hydraulic diameter: Dh = 4 · A / P_wetted
double hydraulic_diameter(double A, double P_wetted);

// Rectangular duct: Dh = 2·a·b / (a + b)
double hydraulic_diameter_rect(double a, double b);

// Annulus: Dh = D_outer - D_inner
double hydraulic_diameter_annulus(double D_outer, double D_inner);
```

### Usage Example

```cpp
#include "incompressible.h"
#include "friction.h"

// Water flow through pipe
double rho = 998.0;  // kg/m³
double D = 0.05;     // m
double L = 10.0;     // m
double v = 2.0;      // m/s

// Calculate Reynolds number and friction factor
double mu = 0.001;   // Pa·s (water viscosity)
double Re = rho * v * D / mu;
double e_D = 0.0001; // Relative roughness
double f = friction_colebrook(Re, e_D);

// Pressure drop
double dP = pipe_dP(v, L, D, f, rho);

// Orifice sizing
double mdot = pipe_mdot(v, D, rho);
double Cd = 0.62;
double dP_orifice = 50000.0;  // Pa
double A_orifice = orifice_area(mdot, 200000.0, 150000.0, Cd, rho);
```

## Orifice Cd Correlations (`orifice.h`)

Discharge coefficient (Cd) correlations for various orifice geometries. Cd relates actual to ideal flow:
```
mdot_actual = Cd * A * sqrt(2 * rho * dP)
```

**Compressible flow note**: For sharp-edged orifices, Cd is weakly dependent on Mach number. These correlations provide adequate Cd values into the choked-flow regime when combined with compressible mass-flow relations (see `compressible.h`).

### Geometry and State Structs

```cpp
struct OrificeGeometry {
    double d = 0.0;       // Orifice bore diameter [m]
    double D = 0.0;       // Pipe diameter [m]
    double t = 0.0;       // Plate thickness [m] (for thick plate)
    double r = 0.0;       // Inlet edge radius [m] (for rounded entry)

    double beta() const;          // Diameter ratio d/D [-]
    double area() const;          // Orifice area [m²]
    double t_over_d() const;      // Thickness ratio t/d [-]
    double r_over_d() const;      // Radius ratio r/d [-]
    bool is_valid() const;
};

struct OrificeState {
    double Re_D = 0.0;    // Pipe Reynolds number (based on D) [-]
    double dP = 0.0;      // Differential pressure [Pa]
    double rho = 0.0;     // Fluid density [kg/m³]
    double mu = 0.0;      // Dynamic viscosity [Pa·s]
};
```

### Cd Correlation Functions

```cpp
// Sharp thin-plate orifice (ISO 5167-2, Reader-Harris/Gallagher)
// Valid for: 0.1 <= beta <= 0.75, Re_D >= 5000, D >= 50mm
double Cd_sharp_thin_plate(const OrificeGeometry& geom, const OrificeState& state);

// Thick-plate orifice (Idelchik thickness correction)
// Valid for: 0 < t/d < ~3
double Cd_thick_plate(const OrificeGeometry& geom, const OrificeState& state);

// Rounded-entry orifice (Idelchik-based)
// Valid for: 0 < r/d <= 0.2
double Cd_rounded_entry(const OrificeGeometry& geom, const OrificeState& state);

// Auto-select correlation based on geometry
double Cd(const OrificeGeometry& geom, const OrificeState& state);
```

### Flow Calculations with Cd

```cpp
// Mass flow through orifice
double orifice_mdot(const OrificeGeometry& geom, double Cd, double dP, double rho);

// Pressure drop for given mass flow
double orifice_dP(const OrificeGeometry& geom, double Cd, double mdot, double rho);

// Solve for Cd from measurement
double orifice_Cd_from_measurement(const OrificeGeometry& geom,
                                    double mdot, double dP, double rho);
```

### Utility Functions

```cpp
namespace orifice {
    // Loss coefficient K from Cd: K = (1/Cd² - 1) * (1 - beta⁴)
    double K_from_Cd(double Cd, double beta);

    // Cd from loss coefficient K
    double Cd_from_K(double K, double beta);

    // Thickness correction factor (multiplies thin-plate Cd)
    double thickness_correction(double t_over_d, double beta);
}
```

### Usage Example

```cpp
#include "orifice.h"

// Define geometry: 50mm orifice in 100mm pipe
OrificeGeometry geom;
geom.d = 0.050;
geom.D = 0.100;

// Define flow state
OrificeState state;
state.Re_D = 100000.0;
state.dP = 10000.0;   // 10 kPa
state.rho = 1.2;      // kg/m³

// Get Cd (auto-selects thin-plate for this geometry)
double Cd = Cd(geom, state);  // ~0.61

// Calculate mass flow
double mdot = orifice_mdot(geom, Cd, state.dP, state.rho);

// For thick plate, set thickness
geom.t = 0.010;  // 10mm
double Cd_thick = Cd_thick_plate(geom, state);  // Higher than thin plate

// For rounded entry, set radius
geom.t = 0.0;
geom.r = 0.008;  // 8mm radius
double Cd_round = Cd_rounded_entry(geom, state);  // ~0.95-0.98
```

### Python Usage

```python
import combaero as cb

# Define geometry
geom = cb.OrificeGeometry()
geom.d = 0.050
geom.D = 0.100

# Define flow state
state = cb.OrificeState()
state.Re_D = 100000.0
state.dP = 10000.0
state.rho = 1.2

# Get Cd
Cd = cb.Cd_orifice(geom, state)

# Calculate mass flow
mdot = cb.orifice_mdot_Cd(geom, Cd, state.dP, state.rho)

# Convert between Cd and loss coefficient K
K = cb.orifice_K_from_Cd(Cd, geom.beta)
Cd_back = cb.orifice_Cd_from_K(K, geom.beta)
```

## Heat Transfer (`heat_transfer.h`)

Convective heat transfer correlations, overall heat transfer coefficients, and utilities for multi-layer wall analysis.

### Nusselt Number Correlations

All correlations are for turbulent internal pipe flow.

```cpp
// Dittus-Boelter (1930) - simple, widely used
// Nu = 0.023 * Re^0.8 * Pr^n  (n=0.4 heating, n=0.3 cooling)
// Valid: Re > 10000, 0.6 < Pr < 160, L/D > 10
double nusselt_dittus_boelter(double Re, double Pr, bool heating = true);

// Gnielinski (1976) - more accurate, wider range
// Valid: 2300 < Re < 5e6, 0.5 < Pr < 2000
double nusselt_gnielinski(double Re, double Pr);
double nusselt_gnielinski(double Re, double Pr, double f);  // with known friction factor

// Sieder-Tate (1936) - includes viscosity correction
// Valid: Re > 10000, 0.7 < Pr < 16700, L/D > 10
double nusselt_sieder_tate(double Re, double Pr, double mu_bulk, double mu_wall);

// Petukhov (1970) - high accuracy for smooth tubes
// Valid: 10^4 < Re < 5×10^6, 0.5 < Pr < 2000
double nusselt_petukhov(double Re, double Pr);
double nusselt_petukhov(double Re, double Pr, double f);
```

### Heat Transfer Coefficient

```cpp
// HTC from Nusselt number: h = Nu * k / D
double htc_from_nusselt(double Nu, double k, double D);

// State-based convenience (uses Gnielinski)
double nusselt_pipe(const State& state, double V, double D);
double htc_pipe(const State& state, double V, double D);
```

### Overall Heat Transfer Coefficient

For composite walls with convective boundaries and conductive layers:

```cpp
// General form: 1/U = Σ(1/h_i) + Σ(t_j/k_j)
// h_values: convective HTCs [W/(m²·K)]
// t_over_k: thickness/conductivity for each layer [m²·K/W]
double overall_htc(const std::vector<double>& h_values,
                   const std::vector<double>& t_over_k);

// Two-sided wall: hot fluid | wall layers | cold fluid
double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double>& t_over_k_layers);

// Multi-layer with fouling resistances
double overall_htc_wall_multilayer(double h_inner, double h_outer,
                                    const std::vector<double>& t_over_k_layers,
                                    double R_fouling_inner = 0.0,
                                    double R_fouling_outer = 0.0);
```

### Thermal Resistance

```cpp
// Convective: R = 1/(h·A)
double thermal_resistance(double h, double A);

// Conductive: R = t/(k·A)
double thermal_resistance_wall(double t, double k, double A);
```

### Log Mean Temperature Difference (LMTD)

```cpp
// General LMTD from temperature differences
double lmtd(double dT1, double dT2);

// Counter-flow heat exchanger
// dT1 = T_hot_in - T_cold_out, dT2 = T_hot_out - T_cold_in
double lmtd_counterflow(double T_hot_in, double T_hot_out,
                        double T_cold_in, double T_cold_out);

// Parallel-flow heat exchanger
// dT1 = T_hot_in - T_cold_in, dT2 = T_hot_out - T_cold_out
double lmtd_parallelflow(double T_hot_in, double T_hot_out,
                         double T_cold_in, double T_cold_out);
```

### Heat Rate and Flux

```cpp
// Heat transfer rate: Q = U * A * ΔT [W]
double heat_rate(double U, double A, double dT);

// Heat flux: q = U * ΔT [W/m²]
double heat_flux(double U, double dT);

// Inverse: solve for area or temperature difference
double heat_transfer_area(double Q, double U, double dT);
double heat_transfer_dT(double Q, double U, double A);
```

### Wall Temperature Profile

Compute temperatures at each layer interface through a multi-layer wall:

```cpp
// Returns vector of interface temperatures (N+1 values for N layers)
// Edge 0: hot-side wall surface
// Edge 1: interface between layer 0 and layer 1
// ...
// Edge N: cold-side wall surface
std::vector<double> wall_temperature_profile(
    double T_hot, double T_cold,      // Bulk fluid temperatures [K]
    double h_hot, double h_cold,      // Convective HTCs [W/(m²·K)]
    const std::vector<double>& t_over_k,  // Layer resistances [m²·K/W]
    double& q);                       // Output: heat flux [W/m²]

// Without heat flux output
std::vector<double> wall_temperature_profile(
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k);
```

### Heat Flux from Measured Temperature

Infer heat flux from thermocouple measurements:

```cpp
// From temperature at a wall edge (layer interface)
double heat_flux_from_T_at_edge(
    double T_measured, std::size_t edge_idx,
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k);

// From temperature at arbitrary depth (embedded thermocouple)
double heat_flux_from_T_at_depth(
    double T_measured, double depth_from_hot,
    double T_hot, double T_cold,
    double h_hot, double h_cold,
    const std::vector<double>& thicknesses,
    const std::vector<double>& conductivities);

// Infer bulk temperature from edge T and known heat flux
double bulk_T_from_edge_T_and_q(
    double T_measured, std::size_t edge_idx, double q,
    double h_hot, double h_cold,
    const std::vector<double>& t_over_k,
    const std::string& solve_for);  // "hot" or "cold"
```

### Temperature Sensitivity (Risk Assessment)

Partial derivatives for uncertainty propagation and risk assessment:

```cpp
// Sensitivity to bulk temperature changes
// Returns: ∂T_edge/∂T_hot [-] (0 to 1)
double dT_edge_dT_hot(std::size_t edge_idx, double h_hot, double h_cold,
                       const std::vector<double>& t_over_k);

// Returns: ∂T_edge/∂T_cold [-] (0 to 1)
double dT_edge_dT_cold(std::size_t edge_idx, double h_hot, double h_cold,
                        const std::vector<double>& t_over_k);

// Both at once: {∂T/∂T_hot, ∂T/∂T_cold}
std::pair<double, double> dT_edge_dT_bulk(
    std::size_t edge_idx, double h_hot, double h_cold,
    const std::vector<double>& t_over_k);

// Sensitivity to heat flux: ∂T_edge/∂q [K·m²/W]
double dT_edge_dq(std::size_t edge_idx, double h_hot,
                   const std::vector<double>& t_over_k);
```

**Key properties:**
- `∂T/∂T_hot + ∂T/∂T_cold = 1` (always, by linearity)
- Edges closer to hot side have higher `∂T/∂T_hot`
- `∂T/∂q` is negative (higher flux → lower edge T)

### NTU-Effectiveness Method

For heat exchanger analysis:

```cpp
// Number of Transfer Units: NTU = U * A / C_min
double ntu(double U, double A, double C_min);

// Capacity ratio: C_r = C_min / C_max
double capacity_ratio(double C_min, double C_max);

// Effectiveness for counter-flow HX
// ε = (1 - exp(-NTU*(1-C_r))) / (1 - C_r*exp(-NTU*(1-C_r)))
double effectiveness_counterflow(double NTU, double C_r);

// Effectiveness for parallel-flow HX
// ε = (1 - exp(-NTU*(1+C_r))) / (1 + C_r)
double effectiveness_parallelflow(double NTU, double C_r);

// Heat rate from effectiveness
// Q = ε * C_min * (T_hot_in - T_cold_in)
double heat_rate_from_effectiveness(double epsilon, double C_min,
                                     double T_hot_in, double T_cold_in);
```

### Usage Example

```cpp
#include "heat_transfer.h"

// Furnace wall: steel + ceramic + steel
std::vector<double> t_over_k = {
    0.010 / 50.0,   // 10mm steel (k=50)
    0.100 / 0.5,    // 100mm ceramic (k=0.5)
    0.005 / 50.0    // 5mm steel (k=50)
};

double h_hot = 150.0;   // Hot gas convection
double h_cold = 10.0;   // Natural convection

// Overall HTC
double U = overall_htc({h_hot, h_cold}, t_over_k);

// Temperature profile
double T_hot = 1473.15;  // 1200°C
double T_cold = 298.15;  // 25°C
double q;
auto temps = wall_temperature_profile(T_hot, T_cold, h_hot, h_cold, t_over_k, q);

// Risk assessment: if gas temperature spikes 100K
auto [dT_hot, dT_cold] = dT_edge_dT_bulk(1, h_hot, h_cold, t_over_k);
double delta_T_interface = dT_hot * 100.0;  // Interface temperature increase
```

### Python Usage

```python
import combaero as ca
import numpy as np

# Multi-layer wall
t_over_k = np.array([0.01/50, 0.1/0.5, 0.005/50])
h_hot, h_cold = 150.0, 10.0

# Temperature profile
temps, q = ca.wall_temperature_profile(500+273, 25+273, h_hot, h_cold, t_over_k)

# Heat flux from thermocouple at steel-ceramic interface
q_measured = ca.heat_flux_from_T_at_edge(temps[1], 1, 773, 298, h_hot, h_cold, t_over_k)

# Risk assessment
dT_hot, dT_cold = ca.dT_edge_dT_bulk(1, h_hot, h_cold, t_over_k)
print(f"If T_hot +100K → interface +{dT_hot*100:.1f}K")

# NTU-effectiveness
NTU = ca.ntu(U=100, A=10, C_min=500)
eps = ca.effectiveness_counterflow(NTU, C_r=0.5)
Q = ca.heat_rate_from_effectiveness(eps, 500, 400, 300)
```

---

## Acoustics

Functions for acoustic mode analysis, resonator design, and damping estimation.

### Geometry

```python
# Tube geometry
tube = ca.Tube(L=0.5, D=0.05)  # Length, Diameter [m]
tube.area()       # Cross-sectional area [m²]
tube.volume()     # Volume [m³]
tube.perimeter()  # Circumference [m]

# Annular geometry (gas turbine combustors)
annulus = ca.Annulus(L=0.3, D_inner=0.5, D_outer=0.7)
annulus.D_mean()        # Mean diameter [m]
annulus.gap()           # Annular gap [m]
annulus.area()          # Cross-sectional area [m²]
annulus.circumference() # Mean circumference [m]
```

### Acoustic Modes

```python
# Boundary conditions
ca.BoundaryCondition.Open    # Pressure release (p'=0)
ca.BoundaryCondition.Closed  # Rigid wall (u'=0)

# Tube axial modes
tube = ca.Tube(L=1.0, D=0.1)
c = 343.0  # Speed of sound [m/s]
modes = ca.tube_axial_modes(tube, c,
                             ca.BoundaryCondition.Closed,
                             ca.BoundaryCondition.Open,
                             n_max=5)

for m in modes:
    print(f"{m.label()}: {m.frequency:.1f} Hz")
# Output: 1L: 85.8 Hz, 2L: 257.3 Hz, ...

# Annulus modes (axial + azimuthal)
combustor = ca.Annulus(L=0.3, D_inner=0.5, D_outer=0.7)
modes = ca.annulus_modes(combustor, c=800,
                          ca.BoundaryCondition.Closed,
                          ca.BoundaryCondition.Open,
                          n_max=3, m_max=3)
# Returns: 1T, 1L, 1L1T, 2T, ...

# Utility functions
danger_modes = ca.modes_in_range(modes, f_min=200, f_max=600)
closest = ca.closest_mode(modes, f_target=500)
min_sep = ca.min_mode_separation(modes)
```

### Mean Flow Correction

```python
# Frequency shift due to mean flow (axial modes only)
f0 = 500.0  # Quiescent frequency [Hz]
M = 0.15    # Mach number

f_up = ca.axial_mode_upstream(f0, M)      # f+ = f0/(1-M)
f_down = ca.axial_mode_downstream(f0, M)  # f- = f0/(1+M)
f_up, f_down = ca.axial_mode_split(f0, M) # Both at once
```

### Resonator Frequencies

```python
# Helmholtz resonator
# f = (c/2π) * √(A/(V*L_eff))
f = ca.helmholtz_frequency(
    V=0.0005,           # Cavity volume [m³]
    A_neck=0.0003,      # Neck area [m²]
    L_neck=0.03,        # Neck length [m]
    c=343.0,            # Speed of sound [m/s]
    end_correction=0.85 # Flanged (default), use 0.6 for unflanged
)

# Quarter-wave and half-wave tubes
f_qw = ca.quarter_wave_frequency(L=0.25, c=343)  # c/(4L)
f_hw = ca.half_wave_frequency(L=0.5, c=343)      # c/(2L)
```

### Strouhal Number

```python
# Dimensionless frequency for vortex shedding, flame dynamics
St = ca.strouhal(f=100, L=0.01, u=5.0)  # St = fL/u
f = ca.frequency_from_strouhal(St=0.2, L=0.01, u=5.0)  # f = St*u/L
```

### Viscothermal Damping (Screening)

Educational screening tools for resonator Q estimation. Expect ±50% accuracy.
Real designs require experimental validation.

```python
# Boundary layer thicknesses
nu = 1.5e-5     # Kinematic viscosity [m²/s]
alpha = 2.1e-5  # Thermal diffusivity [m²/s]
f = 400.0       # Frequency [Hz]
gamma = 1.4     # Specific heat ratio

delta_nu = ca.stokes_layer(nu, f)      # δ_ν = √(2ν/ω)
delta_kappa = ca.thermal_layer(alpha, f)  # δ_κ = √(2α/ω)
delta_eff = ca.effective_viscothermal_layer(delta_nu, delta_kappa, gamma)

# Quality factor estimates
Q_tube = ca.tube_Q(L=0.2, D=0.025, nu=nu, alpha=alpha, gamma=gamma, f=f)
Q_helm = ca.helmholtz_Q(V=0.0005, A_neck=0.0003, L_neck=0.03,
                         nu=nu, alpha=alpha, gamma=gamma, f=200)

# Derived quantities
zeta = ca.damping_ratio(Q_tube)   # ζ = 1/(2Q)
df = ca.bandwidth(f, Q_tube)      # Δf = f/Q
```

### Residence Time

```python
# Time for fluid to pass through a volume
tau = ca.residence_time(V=0.001, Q=0.0001)  # τ = V/Q̇ [s]
tau = ca.residence_time_mdot(V=0.001, mdot=0.001, rho=1.2)  # τ = Vρ/ṁ

# With geometry structs
tau = ca.residence_time_tube(tube, Q=0.01)
tau = ca.residence_time_annulus(annulus, Q=0.1)

# Space velocity (inverse)
SV = ca.space_velocity(Q=0.01, V=0.001)  # SV = 1/τ [1/s]
```

### Complete Example: Combustor Acoustic Analysis

```python
import combaero as ca
import math

# Gas turbine annular combustor
combustor = ca.Annulus(L=0.3, D_inner=0.5, D_outer=0.7)
c = 800.0  # Hot products speed of sound [m/s]

# Get all modes
modes = ca.annulus_modes(combustor, c,
                          ca.BoundaryCondition.Closed,
                          ca.BoundaryCondition.Open,
                          n_max=3, m_max=3)

print("Acoustic modes:")
for m in modes[:5]:
    print(f"  {m.label():6s} {m.frequency:7.1f} Hz")

# Check for modes in instability-prone range
danger = ca.modes_in_range(modes, 200, 600)
print(f"\nModes in 200-600 Hz range: {[m.label() for m in danger]}")

# Design quarter-wave damper for 1T mode
f_1T = modes[0].frequency
L_damper = ca.quarter_wave_frequency.__doc__  # c/(4L) → L = c/(4f)
L_damper = c / (4 * f_1T)
print(f"\nQuarter-wave damper for {modes[0].label()}: L = {L_damper*100:.1f} cm")

# Estimate damper Q (screening)
nu, alpha, gamma = 3e-5, 4e-5, 1.3  # Hot gas properties
Q = ca.tube_Q(L_damper, D=0.03, nu=nu, alpha=alpha, gamma=gamma, f=f_1T)
print(f"Estimated Q ≈ {Q:.0f} (±50%, validate experimentally)")
```
