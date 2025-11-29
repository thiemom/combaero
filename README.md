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

- Ideal-gas mixture properties: Cp, H, S, density, speed of sound
- Gibbs free energy helpers `g_over_RT` and `dg_over_RT_dT`
- Combustion tools: O₂ demand per fuel/mixture, complete combustion to CO₂/H₂O
- Equivalence ratio and Bilger mixture fraction utilities
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
| `thermo.h` | Thermodynamic properties (cp, h, s, density, etc.) and physical constants |
| `transport.h` | Transport properties (viscosity, thermal conductivity, Prandtl) |
| `combustion.h` | Combustion calculations (O₂ demand, equivalence ratio, Bilger) |
| `equilibrium.h` | Chemical equilibrium solver (WGS partial equilibrium) |
| `humidair.h` | Humid air properties and saturation vapor pressure |
| `utils.h` | Utility functions (mixture property printing) |

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

See the *Development workflow* section above for concrete commands.

## License

[MIT License](LICENSE)
