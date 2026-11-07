# CombAero

**C++17 thermodynamics, combustion and flow-network library with Python bindings, for gas turbine combustor and cooling design.**

## Key Capabilities

- **🚀 Thermodynamics**: NASA-9 polynomials for multi-species gas mixtures (200-6000 K; to 20 000 K for N2, O2, Ar, CO2, CO, H2).
- **🔥 Combustion**: Complete combustion, partial equilibrium (reforming and water-gas shift, no dissociation), and inverse solvers.
- **💨 Fluid Dynamics**: Compressible flow (nozzle, Fanno, Rayleigh), friction correlations, and orifice models.
- **❄️ Cooling**: Correlations for rib enhancement, impingement, film and effusion cooling, and pin fins.
- **🕸️ Network Solver**: Fluid-thermal networks solved with Newton-type methods (MINPACK) on C++ element kernels with analytic Jacobians.
- **🐍 Python Native**: High-level Python bindings for rapid prototyping and data analysis.

Correlations are checked against their source's own data and, separately, against other sources; how is in the [Validation Policy](https://github.com/thirox-lab/combaero/blob/main/docs/VALIDATION_POLICY.md).

## Quick Start

### Python (Recommended)
Install the pre-built wheel (Python 3.12+; we recommend [uv](https://docs.astral.sh/uv/)):
```bash
uv pip install combaero
```

Calculate adiabatic flame temperature:
```python
import combaero as cb
state = cb.State().set_TPX(300, 101325, "CH4:1, O2:2, N2:7.52")
burned = cb.complete_combustion(state.T, state.X)
print(f"Adiabatic Flame Temperature: {burned.T:.0f} K")
```

Example scripts, including worked networks from papers, are in [`python/examples`](https://github.com/thirox-lab/combaero/tree/main/python/examples).

### C++
Add `include/` to your include path and link the library ([Building CombAero](https://github.com/thirox-lab/combaero/blob/main/docs/BUILDING.md)):
```cpp
#include <combustion.h>
#include <state.h>
#include <thermo.h>

// Mole fractions in the library's species order: CH4 + 2 (O2 + 3.76 N2).
std::vector<double> X(combaero::num_species(), 0.0);
X[combaero::species_index_from_name("CH4")] = 1.0;
X[combaero::species_index_from_name("O2")] = 2.0;
X[combaero::species_index_from_name("N2")] = 7.52;

combaero::State in;
in.set_T(300.0).set_P(101325.0).set_X(X);
auto burned = combaero::complete_combustion(in);
```

## GUI

For an interactive drag-and-drop network designer, install **[combaero-gui](https://pypi.org/project/combaero-gui/)**; its Examples menu loads the worked examples:
```bash
uv pip install combaero-gui  # then run: combaero-gui
```

## Documentation

- **[API Reference (Python)](https://github.com/thirox-lab/combaero/blob/main/docs/API_PYTHON.md)**: Python usage and examples.
- **[API Reference (C++)](https://github.com/thirox-lab/combaero/blob/main/docs/API_CPP.md)**: Reference for C++ developers.
- **[Units Guide](https://github.com/thirox-lab/combaero/blob/main/docs/UNITS.md)**: SI unit system and physical constants.
- **[Building CombAero](https://github.com/thirox-lab/combaero/blob/main/docs/BUILDING.md)**: Prerequisites and build instructions.
- **[Development Workflow](https://github.com/thirox-lab/combaero/blob/main/docs/DEVELOPMENT.md)**: Testing, styling, and contribution guide.
- **[Changelog](https://github.com/thirox-lab/combaero/blob/main/CHANGELOG.md)**: What changed in each release.

## License
[MIT License](https://github.com/thirox-lab/combaero/blob/main/LICENSE)
