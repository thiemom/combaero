# CombAero

A C++ library for thermodynamic and combustion calculations for gas mixtures.

## Overview

This library provides tools for ideal-gas thermodynamic properties, combustion stoichiometry, transport properties, and humid air calculations. It includes functionality for:

- Thermodynamic property calculations for multi-species gas mixtures
- Combustion utilities (stoichiometric O₂ demand, complete combustion products)
- Saturation vapor pressure calculations for water and ice
- Humid air composition and property calculations
- Mole fraction utilities
- Python bindings for these utilities via the `combaero` package

## Features

- **NASA-9 thermodynamic polynomials** (200-6000 K, some species to 20000 K)
- Ideal-gas mixture properties: Cp, H, S, density, speed of sound
- Gibbs free energy helpers `g_over_RT` and `dg_over_RT_dT`
- Combustion tools: O₂ demand per fuel/mixture, complete combustion to CO₂/H₂O
- Equivalence ratio and Bilger mixture fraction utilities
- **Inverse combustion solvers**: find fuel or oxidizer flow for target Tad, O₂, or CO₂
- **Compressible flow**: isentropic nozzle flow and Fanno flow (adiabatic pipe with friction)
- **Friction factors**: Colebrook-White, Haaland, and Serghides correlations
- **Orifice Cd correlations**: ISO 5167, thick-plate, rounded-entry (Idelchik/Bohl)
- **Heat transfer correlations**: Dittus-Boelter, Gnielinski, Sieder-Tate, Petukhov
- **Acoustics**: Tube/annulus modes, Helmholtz resonators, Q factor screening
- Accurate saturation vapor pressure using Hyland-Wexler equations
  - For ice (-80°C to 0°C): maximum relative error ≤ 0.023%
  - For water vapor (0°C to 80°C): maximum relative error ≤ 0.0057%
- Dry air and humid air utilities (humidity ratio, compositions, dew point)
- Mole fraction normalization and wet/dry conversion utilities
- Transport properties: viscosity, thermal conductivity, Prandtl number

## Library Structure

The library is organized into focused modules:

| Header | Description |
|--------|-------------|
| `state.h` | `State` and `Stream` structs for thermodynamic state representation |
| `thermo.h` | Thermodynamic properties (cp, h, s, density, etc.) and physical constants |
| `transport.h` | Transport properties (viscosity, thermal conductivity, Prandtl) |
| `combustion.h` | Combustion calculations (O₂ demand, equivalence ratio, complete combustion) |
| `equilibrium.h` | Chemical equilibrium solver (WGS partial equilibrium) |
| `compressible.h` | Compressible flow (isentropic nozzle, Fanno pipe flow) |
| `friction.h` | Darcy friction factor correlations (Colebrook, Haaland, Serghides) |
| `heat_transfer.h` | Nusselt correlations, LMTD for heat exchangers |
| `geometry.h` | Hydraulic diameter, Tube/Annulus geometry, residence time |
| `incompressible.h` | Incompressible flow (Bernoulli, pipe flow, orifice) |
| `orifice.h` | Orifice discharge coefficient (Cd) correlations |
| `acoustics.h` | Acoustic modes, Helmholtz/quarter-wave resonators, Q factors |
| `humidair.h` | Humid air properties and saturation vapor pressure |
| `utils.h` | Utility functions (mixture property printing) |

### State-based API

The library provides `State` and `Stream` structs for representing thermodynamic states:

```cpp
struct State {
    double T = 298.15;           // Temperature [K]
    double P = 101325.0;         // Pressure [Pa]
    std::vector<double> X;       // Mole fractions [-]
    
    // Thermodynamic properties
    double mw() const;           // Molecular weight [g/mol]
    double cp() const;           // Specific heat [J/(mol·K)]
    double h() const;            // Enthalpy [J/mol]
    double rho() const;          // Density [kg/m³]
    // ... and more (s, cv, u, R, gamma, a)
    
    // Transport properties
    double mu() const;           // Dynamic viscosity [Pa·s]
    double k() const;            // Thermal conductivity [W/(m·K)]
    double nu() const;           // Kinematic viscosity [m²/s]
    double Pr() const;           // Prandtl number [-]
    double alpha() const;        // Thermal diffusivity [m²/s]
    
    // Setters with chaining
    State& set_T(double T);
    State& set_P(double P);
    State& set_X(const std::vector<double>& X);
};

struct Stream {
    State state;
    double mdot = 0.0;           // Mass flow rate [kg/s]
};
```

Combustion and equilibrium functions accept and return `State` objects:

```cpp
State in;
in.set_T(300.0).set_P(101325.0).set_X(X_unburned);

// Adiabatic complete combustion
State burned = complete_combustion(in);
std::cout << "Adiabatic flame T: " << burned.T << " K\n";
std::cout << "Flame density: " << burned.rho() << " kg/m³\n";
std::cout << "Flame viscosity: " << burned.mu() << " Pa·s\n";

// WGS equilibrium (isothermal or adiabatic)
State eq_iso = wgs_equilibrium(in);
State eq_ad = wgs_equilibrium_adiabatic(in);
```

### Stream Mixing

Mix multiple streams with mass and enthalpy balance:

```cpp
Stream air, fuel;
air.state.set_T(400.0).set_P(101325.0).set_X(X_air);
air.mdot = 10.0;  // kg/s

fuel.state.set_T(300.0).set_P(101325.0).set_X(X_fuel);
fuel.mdot = 0.5;  // kg/s

// Mix streams (uses minimum inlet pressure by default)
Stream mixed = mix({air, fuel});

// Or specify output pressure explicitly
Stream mixed2 = mix({air, fuel}, 150000.0);
```

## Building the Project

### Prerequisites

- C++17 compatible compiler
- CMake 3.10 or higher

### Build Instructions

```bash
# Clone the repository
git clone https://github.com/thiemom/combaero.git
cd combaero

# Create a build directory
mkdir -p build
cd build

# Configure and build
cmake ..
make

# Run tests
make test
```

## Examples

The `examples` directory contains sample applications demonstrating the library's functionality:

- `thermo_example`: Basic thermodynamic calculations
- `humidair_example`: Humid air property calculations

## Development workflow

### Environment setup (direnv + Python)

These steps assume you have `direnv` and `uv` installed.

```bash
git clone https://github.com/thiemom/combaero.git
cd combaero

# Allow direnv to manage the root .venv
direnv allow

# direnv will create .venv/ and uv will bootstrap Python tooling
which python   # should point to .venv/bin/python inside the repo
```

The root `.envrc` manages a single virtual environment at `./.venv` for:

- Building the Python wheel via `scikit-build-core`
- Running Python tests with `pytest`
- Running data generator scripts in `thermo_data_generator/`
- Running `scripts/generate_units_md.py` to regenerate `docs/UNITS.md`

### Build C++ and run C++ tests

```bash
mkdir -p build
cd build

cmake ..
cmake --build .

# Run all C++ tests
ctest --output-on-failure

# Or invoke the main test executable directly
./tests/combaero_tests

# Accuracy tests for saturation vapor pressure
./tests/test_ice_equation_accuracy
./tests/test_water_equation_accuracy
```

### Build and test the Python wheel

From the repository root (with direnv active):

```bash
# Optionally clean old wheels
rm -f dist/combaero-*.whl

# Build the wheel
python -m build --wheel

# Install (or reinstall) the wheel into the root .venv
python -m pip install --force-reinstall dist/combaero-0.0.1-*.whl

# Run Python tests
python -m pytest python/tests
```

## Testing

The project uses GoogleTest and CTest for C++ tests, and pytest for Python tests:

- Core thermodynamic and transport property tests
- Combustion, equivalence-ratio, and Bilger mixture-fraction tests
- Humid air and saturation pressure tests
- Accuracy validation against Hyland-Wexler reference values
- Python tests for the `combaero` API (thermo/transport, humid air, combustion helpers, species common names)
- **Cantera validation tests**: Comprehensive validation against Cantera for combustion, mixing, and transport properties

See the *Development workflow* section above for concrete commands.

### Cantera Validation Tests

A comprehensive test suite validates CombAero's combustion, mixing, and transport calculations against Cantera as a reference implementation. These tests are integrated into the pre-commit workflow.

```bash
cd cantera_validation_tests

# Install dependencies
poetry install

# Run all validation tests
poetry run pytest -v

# Run specific test categories
poetry run pytest test_combustion_validation.py -v
poetry run pytest test_mixing_validation.py -v
poetry run pytest test_transport_validation.py -v

# Run in parallel
poetry run pytest -n auto
```

The validation tests cover:
- **Combustion**: Adiabatic flame temperature, product composition, equivalence ratio
- **Mixing**: Stream mixing with enthalpy balance, density calculations
- **Transport**: Viscosity, thermal conductivity, Prandtl number at various conditions

See [`cantera_validation_tests/README.md`](cantera_validation_tests/README.md) and [`cantera_validation_tests/INTEGRATION_GUIDE.md`](cantera_validation_tests/INTEGRATION_GUIDE.md) for detailed documentation.

## Documentation

- **[API Reference](docs/API_REFERENCE.md)**: Detailed function signatures and usage examples (useful for LLMs)
- **[Units Reference](docs/UNITS.md)**: SI unit system used throughout the library (auto-generated from `units_data.h`)

## Acknowledgments

This library was developed with significant assistance from AI coding agents, including Claude (Anthropic) for the majority of implementation and testing, and ChatGPT (OpenAI) for early water-gas shift equilibrium guidance. These tools made it possible to build a comprehensive thermodynamic library in a fraction of the time it would otherwise require.

## License

[MIT License](LICENSE)
