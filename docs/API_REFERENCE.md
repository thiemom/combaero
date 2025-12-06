# CombAero API Reference

This document provides detailed API reference for LLMs and developers.

## Species

The library uses a fixed set of 12 species with 0-indexed positions:

| Index | Formula | Common Name |
|-------|---------|-------------|
| 0 | N2 | Nitrogen |
| 1 | O2 | Oxygen |
| 2 | H2 | Hydrogen |
| 3 | H2O | Water |
| 4 | CO | Carbon Monoxide |
| 5 | CO2 | Carbon Dioxide |
| 6 | CH4 | Methane |
| 7 | C2H6 | Ethane |
| 8 | C3H8 | Propane |
| 9 | C4H10 | Butane |
| 10 | Ar | Argon |
| 11 | He | Helium |

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
double viscosity(double T, const std::vector<double>& X);           // Pa·s
double thermal_conductivity(double T, const std::vector<double>& X); // W/(m·K)
double prandtl(double T, const std::vector<double>& X);             // dimensionless
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
