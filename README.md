# Thermo

A C++ library for thermodynamic calculations with a focus on humid air properties.

## Overview

This library provides tools for calculating thermodynamic properties of gases, with a particular focus on humid air calculations using the Hyland-Wexler equations. It includes functionality for:

- Thermodynamic property calculations
- Saturation vapor pressure calculations for water and ice
- Humid air composition and property calculations
- Mole fraction utilities

## Features

- Accurate saturation vapor pressure calculations using Hyland-Wexler equations
  - For ice (-80°C to 0°C): Maximum relative error ≤ 0.023%
  - For water vapor (0°C to 80°C): Maximum relative error ≤ 0.0057%
- Dry air composition handling (N₂, O₂, Ar, CO₂)
- Humidity ratio calculations
- Mole fraction normalization and conversion utilities

## Building the Project

### Prerequisites

- C++17 compatible compiler
- CMake 3.10 or higher

### Build Instructions

```bash
# Clone the repository
git clone https://github.com/yourusername/thermo.git
cd thermo

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

The project includes comprehensive tests to validate the accuracy of calculations:

- Core thermodynamic property tests
- Humid air calculation tests
- Mole fraction utility tests
- Accuracy validation against Hyland-Wexler reference values

Run the tests with:

```bash
cd build
./tests/thermo_tests
```

For specific accuracy tests:

```bash
./tests/test_ice_equation_accuracy
./tests/test_water_equation_accuracy
```

## License

[MIT License](LICENSE)
