# CombAero

A C++ library for thermodynamic and combustion calculations for gas mixtures.

## Overview

This library provides tools for ideal-gas thermodynamic properties, combustion stoichiometry, transport properties, and humid air calculations. It includes functionality for:

- Thermodynamic property calculations for multi-species gas mixtures
- Combustion utilities (stoichiometric O₂ demand, complete combustion products)
- Saturation vapor pressure calculations for water and ice
- Humid air composition and property calculations
- Mole fraction utilities

## Features

- Ideal-gas mixture properties: Cp, H, S, density, speed of sound
- Gibbs free energy helpers `g_over_RT` and `dg_over_RT_dT`
- Combustion tools: O₂ demand per fuel/mixture, complete combustion to CO₂/H₂O
- Accurate saturation vapor pressure using Hyland-Wexler equations
  - For ice (-80°C to 0°C): maximum relative error ≤ 0.023%
  - For water vapor (0°C to 80°C): maximum relative error ≤ 0.0057%
- Dry air and humid air utilities (humidity ratio, compositions, dew point)
- Mole fraction normalization and wet/dry conversion utilities

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

## Testing

The project uses GoogleTest and CTest to validate the accuracy of calculations:

- Core thermodynamic and transport property tests
- Combustion and mixture-fraction tests
- Humid air and saturation pressure tests
- Accuracy validation against Hyland-Wexler reference values

From the build directory you can run:

```bash
# Run all tests via CTest
ctest --output-on-failure

# Run only the main combaero test suite
ctest -R combaero_tests --output-on-failure

# Or invoke the test executable directly
./tests/combaero_tests

# Accuracy tests for saturation vapor pressure
./tests/test_ice_equation_accuracy
./tests/test_water_equation_accuracy
```

## License

[MIT License](LICENSE)
