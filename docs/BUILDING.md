# Building CombAero

This document provides detailed instructions for building the CombAero C++ library and its Python bindings.

## Prerequisites

- **C++ Compiler**: A C++17 compatible compiler (GCC 8+, Clang 7+, or MSVC 2019+).
- **CMake**: Version 3.10 or higher.
- **Python**: Version 3.12 or higher (for Python bindings).
- **uv**: Recommended for Python dependency management.

## C++ Library Build

To build the C++ library and its native tests:

```bash
# Create a build directory
mkdir -p build
cd build

# Configure and build
cmake ..
cmake --build .

# Run all C++ tests
ctest --output-on-failure
```

### Build Options

- `-DCOMBAERO_BUILD_TESTS=ON/OFF`: Build C++ unit tests (default: ON).
- `-DCOMBAERO_BUILD_EXAMPLES=ON/OFF`: Build C++ examples (default: ON).
- `-DCMAKE_BUILD_TYPE=Release/Debug`: Specify build configuration.

## Python Bindings Build

The Python package `combaero` uses `scikit-build-core` to compile the C++ extensions.

### Environment Setup

We recommend using the repository-local virtual environment:

```bash
# Create/update .venv and install dev tooling
./scripts/bootstrap.sh

# Activate shell
source .venv/bin/activate
```

### Building the Wheel

```bash
# Build the wheel using build
python -m build --wheel

# Install the wheel locally
python -m pip install dist/combaero-*.whl
```

## Troubleshooting

### Windows (MSVC)
Ensure you use the "Developer Command Prompt for VS" or "x64 Native Tools Command Prompt".

### macOS
Ensure Xcode Command Line Tools are installed: `xcode-select --install`.
