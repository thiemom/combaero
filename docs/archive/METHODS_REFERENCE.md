# Methods & Computational Reference

This document tracks the numerical methods, profiling benchmarks, and computational costs of the CombAero Python bindings. It is intended to guide the architectural design of the Python-level Network Solver, ensuring that inner loops remain performant and Jacobians remain stable.

## 1. Physics Kernel Latency (Profiling)

Evaluating network solver architectures requires understanding the latency of the underlying basic physics evaluations. A network with $N$ elements and $E$ edges, solved via a finite-difference dense Jacobian, requires $O(N \times E)$ element residual evaluations per Newton step. Consequently, identifying bottlenecks in the CombAero C++ bindings is critical for building a fast `scipy.optimize.root` inner loop.

The following table details the execution time per call (in microseconds) for fundamental physics kernels when called from Python. (Tested using `timeit` on Python 3.12, macOS Apple Silicon context).

| Operation | Time (μs) | Time (ns) | Implementation Note |
| :--- | ---: | ---: | :--- |
| **Object Allocation** | | | |
| `s = cb.State()` | 0.497 μs | 497 ns | Python object initialization to C++ pointer |
| `s.TPX = (T, P, X)` | 0.787 μs | 787 ns | Vector translation across pybind11 |
| **Polynomial Evaluations** | | | |
| `cb.cp(T, X)` | 0.691 μs | 691 ns | Horner method on NASA-9 bases |
| `cb.h(T, X)` | 1.016 μs | 1016 ns | NASA-9 bases integration |
| **Equation of State / Acoustics** | | | |
| `cb.density(T, P, X)` | 0.660 μs | 660 ns | Ideal gas law + molar mass summation |
| `cb.speed_of_sound(T, X)` | 0.728 μs | 728 ns | Dependent on $c_p$ evaluation |
| **Transport Properties** | | | |
| `cb.thermal_conductivity(T, P, X)`| 0.912 μs | 912 ns | Collision integrals & mixing rules |
| `cb.viscosity(T, P, X)` | 0.928 μs | 928 ns | |
| `cb.transport_state(T, P, X)` | 1.318 μs | 1318 ns | Combined bundle evaluated at once |
| **Incompressible Flow Elements** | | | |
| `cb.orifice_flow_thermo(Cd)` | 1.382 μs | 1382 ns | Ideal equation + internal Density |
| `cb.channel_flow_rough(Haaland)` | 1.458 μs | 1458 ns | Explicit friction factor + Density/Viscosity |
| `cb.channel_flow_rough(Colebrook)` | 1.678 μs | 1678 ns | Iterative Implicit friction + Density/Viscosity |
| **Compressible ODE Integration** | | | |
| `cb.fanno_channel_rough(10 steps)`| 30.480 μs | 30480 ns | Heavy due to repeated `solve_T_from_energy` Newton sub-loop |
| `cb.fanno_channel_rough(100 steps)`| 290.194 μs | 290194 ns| Far too slow for dense finite-differences |

## 2. API Architectural Implications

### 2.1 Pass Primitives, Not Objects
A single call to instantiate a `cb.State()` object takes approximately `~0.5 μs`, which is nearly as expensive as evaluating the ideal gas law (`0.66 μs`).
**To maximize inner-loop speed**, network nodes should NOT pass `combaero.State()` objects to edge elements during solver iterations. Edges should be called directly with pure numerical primitives: `(T, P, X, mdot)`.
A pure python `dataclass` or NumPy array containing these primitives will be substantially faster than wrapping and unwrapping C++ bindings.

### 2.2 PyBind11 Return Type Latency

When designing the custom solver C++ API, the structure of the data returned across the pybind11 boundary is as critical as the passing of input arguments. A second profiling pass was conducted specifically comparing overheads:

| Return Type | Example Scenario | Time (μs) | Observation |
| :--- | :--- | ---: | :--- |
| **Primitive (`float`)** | `state.T` | 0.113 μs | Instantaneous return. |
| **Tuple (`float`, `float`, `list`)** | `state.TPX` | 0.214 μs | Packing a C++ `std::tuple` into a Python `tuple` is blazingly fast. |
| **Pass-by-ref (numpy array)** | `test_pass_by_ref` | 0.257 μs | Mutating a pre-allocated numpy `py::array_t` requires buffer locking/unlocking, making it slightly slower than pure tuple creation. |
| **List (`std::vector`)** | `mole_to_mass` | 0.766 μs | Array copying incurs slight latency. |
| **Object (Custom Struct)**| `channel_flow_rough` | 1.449 μs | Binding and returning custom C++ classes/structs (`IncompressibleFlowSolution`) adds over 1.2 μs of pure overhead. |

**API Implication:** The dedicated solver C++ methods must return `std::tuple<double, double>` representing `(f, J)` (value and exact jacobian). Attempting to wrap these results into a custom 'SolverElementResult' class object or dictionary drops the API performance by ~600% due to pybind11 conversion overheads. Interestingly, pure tuple returns actually beat pass-by-reference `py::array_t` mutations due to numpy's buffer lifecycle overheads.

### 2.3 Incompressible vs Compressible Cost Scaling
Incompressible calculations (`channel_flow_rough`, `orifice_flow_thermo`) are highly optimized. Because they utilize `v_ideal` rather than iterative Newton loops to find friction factors, they execute in **~1.5 μs** (including the heavy object return overhead noted above). A network containing 100 incompressible channels requires less than `0.2` milliseconds to evaluate its total residuals once.
