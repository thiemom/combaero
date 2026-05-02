# CombAero

**High-performance C++17 engine for aerospace thermodynamics, combustion, and fluid network simulation.**

![CombAero Banner](.github/assets/combaero_banner.png)

CombAero provides aerospace engineers and researchers with a fast, accurate, and extensible toolset for modeling complex gas mixtures, chemical equilibrium, and integrated cooling systems.

## Key Capabilities

- **🚀 Thermodynamics**: NASA-9 polynomials for multi-species gas mixtures (200-20000 K).
- **🔥 Combustion**: Complete combustion, chemical equilibrium (WGS, SMR), and inverse solvers.
- **💨 Fluid Dynamics**: Compressible flow (nozzle, Fanno), friction correlations, and orifice models.
- **❄️ Cooling**: Advanced correlations for rib enhancement, impingement, film cooling, and pin fins.
- **🕸️ Network Solver**: Fast-path native C++ solver for large-scale fluid-thermal networks.
- **🐍 Python Native**: High-level Python bindings for rapid prototyping and data analysis.

## Quick Start

### Python (Recommended)
Install the pre-built wheel:
```bash
pip install combaero
```
Calculate adiabatic flame temperature:
```python
import combaero as cb
state = cb.State().set_TPX(300, 101325, "CH4:1, O2:2, N2:7.52")
burned = cb.complete_combustion(state.T, state.X)
print(f"Adiabatic Flame Temperature: {burned.T:.2f} K")
```

### C++
Include the headers in your project:
```cpp
#include <combaero/combustion.h>
#include <combaero/state.h>

combaero::State in;
in.set_T(300.0).set_P(101325.0).set_X(X_mixture);
auto burned = combaero::complete_combustion(in);
```

## Documentation

- **[API Reference (C++)](docs/API_CPP.md)**: Technical reference for C++ developers and agents.
- **[API Reference (Python)](docs/API_PYTHON.md)**: High-level Python usage and examples.
- **[Units Guide](docs/UNITS.md)**: SI unit system and physical constants.
- **[CombAero GUI](gui/README.md)**: Interactive desktop application for network design.

## Technical Details

- **[Building CombAero](docs/BUILDING.md)**: Prerequisites and build instructions.
- **[Development Workflow](docs/DEVELOPMENT.md)**: Testing, styling, and contribution guide.
- **[Packaging Strategy](docs/PACKAGING.md)**: How we decouple core physics from the GUI.

## License
[MIT License](LICENSE)
