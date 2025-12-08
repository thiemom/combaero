# CombAero API Reference

This document provides detailed API reference for LLMs and developers.

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

### Species Data Files

- **`thermo_transport_data.h`**: Contains all species data including `species_names` vector, `species_index` unordered map, `molar_masses`, `nasa_coeffs`, `transport_props`, and `molecular_structures`. This file is auto-generated.
- **`common_names.h`**: Maps between formulas and human-readable names (in `combaero::` namespace):
  - `formula_to_name` map: `"CH4"` → `"Methane"`
  - `name_to_formula` map: `"Methane"` → `"CH4"`
  - `common_name(formula)` function: lookup with error handling
  - `formula(name)` function: inverse lookup with error handling
- **`thermo_data_generator/`**: Python utility that generates `thermo_transport_data.h` from source data (NASA polynomials, transport properties).

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
