# Network & GUI Integration Review

**Date**: 2026-04-11
**Scope**: `python/combaero/network/`, `gui/backend/`, `gui/frontend/`

## Design Principles

1. **C++ backend for all calculations and derivatives** — Python code delegates all
   thermodynamic calculations, derivatives, and Jacobian entries to the C++ backend.
2. **Smooth handling of non-physical states** — the solver gracefully handles transient
   non-physical states during iteration by smoothly bounding them or returning penalty
   residuals, rather than throwing exceptions that terminate the solver.
3. **Efficient post-solve state extraction** — after convergence, all node and element
   states should be computed once via C++ `CompleteState` bundles, not re-derived
   multiple times.

---

## Findings

### HIGH impact

#### H1 — Redundant evaluations post-solve
- **File**: `solver.py:1356, 1385, 1422`
- **Issue**: After `root()` converges, `_residuals(final_x)` is called to check the
  norm (line 1356). This triggers a full `_propagate_states` + residual + Jacobian
  evaluation, then discards the Jacobian. Then `_propagate_states(final_x)` is called
  *again* (line 1422) for post-solve dict assembly. For the fallback path (line 1385),
  a third full evaluation occurs.
- **Cost**: 2-3 unnecessary full evaluations after convergence.
- **Fix**: Track the final residual norm from the last `residuals_wrapper` call (the
  `best_res_norm` tracker already exists for best-iterate). Remove the standalone
  `_propagate_states` at line 1422 since the final `_residuals` call already populated
  `_derived_states`.

#### H2 — GUI re-computes all states the solver already computed
- **File**: `gui/backend/main.py:93-123`
- **Issue**: The solver's `sol_dict` already contains all diagnostic keys
  (`{nid}.h`, `{nid}.rho`, `{nid}.mach`, `{eid}.mach_in`, etc.) computed via
  `cb.complete_state`. The GUI backend ignores these and re-calls
  `node.diagnostics(state_obj)` and `elem.diagnostics(mix_in, mix_out)` for every
  node and element, triggering a second (or third) `complete_state` call each.
- **Cost**: ~3x redundant `complete_state` C++ calls per node, ~2x per element.
- **Fix**: Extract diagnostics directly from `sol_dict` instead of re-calling
  diagnostics methods. The solver result already has everything the GUI needs.

#### H3 — GUI uses C++ MixtureState for Python diagnostics (type mismatch)
- **File**: `gui/backend/main.py:93`
- **Issue**: `cb.MixtureState` is the **C++ pybind11** struct, but
  `node.diagnostics()` expects the **Python dataclass** `MixtureState` from
  `components.py`. This works by duck-typing today but the `.X` property caching
  and any future method additions will not work on the C++ object.
- **Fix**: Eliminated when H2 is fixed (GUI stops constructing MixtureState objects).

### MEDIUM impact

#### M1 — `_residuals()` always computes and discards Jacobian
- **File**: `solver.py:1063-1068`
- **Issue**: `_residuals()` calls `_residuals_and_jacobian()` and discards the sparse
  Jacobian. Used at lines 1273, 1356, 1385 for dimension/norm checks. The Jacobian is
  ~50% of the computational work.
- **Fix**: Add a `_residuals_only()` method or a flag to skip COO/CSR matrix
  construction when only residuals are needed.

#### M2 — `extract_complete_states()` exists but is unused
- **File**: `solver.py:1514-1593`
- **Issue**: This method computes `CompleteState` for all nodes — exactly the pattern
  needed for efficient post-solve extraction. Neither `solve()` nor the GUI uses it.
- **Fix**: Integrate into the post-solve path: after convergence, call
  `complete_state(T, P, X)` once per node/element, attach to results.

#### M3 — Missing Jacobian coupling for auto_Cd
- **File**: `components.py:1080-1108` (`OrificeElement._effective_Cd`)
- **Issue**: Reynolds number is computed in Python to pass to `cb.Cd_orifice`. The
  returned Cd is treated as constant in the Jacobian (no `dCd/dRe` sensitivity).
  Missing `dCd/dm_dot` and `dCd/dT` terms.
- **Fix**: Either move the Cd-and-Jacobian computation into a single C++ wrapper, or
  document the Jacobian truncation as a known limitation. Impact on convergence is
  typically small (Cd varies slowly with Re).

#### M4 — `PipeElement.htc_and_T` missing `rho > 0` guard
- **File**: `components.py:1629-1630`
- **Issue**: `rho = state.density(); u = state.m_dot / (rho * self.area)` will throw
  `ZeroDivisionError` if rho <= 0 during non-physical intermediate states. The penalty
  wrapper catches this, but `MomentumChamberNode.htc_and_T` (line 655-656) already has
  the correct guard pattern: `u = m_dot / (rho * self.area) if rho > 0 and self.area > 0 else 0.0`.
- **Fix**: Add same guard to `PipeElement.htc_and_T`.

#### M5 — `mass_to_mole` called repeatedly for the same composition
- **File**: `gui/backend/main.py:99`, `components.py:358-361`
- **Issue**: For each node, `mass_to_mole(Y)` is called independently in: (a) the
  `MixtureState.X` property, (b) every `complete_state` call, (c) the GUI's explicit
  `cb.mass_to_mole(y_vals)`. Total: 2-3 redundant conversions per node per evaluation.
- **Fix**: Eliminated when H2 is fixed. For the solver path, the `MixtureState.X`
  caching already helps — ensure it is used consistently.

### LOW impact

#### L1 — Legacy `combustion.py` `mix_streams()` duplicates C++ mixing
- **File**: `combustion/combustion.py:22-106`
- **Issue**: Performs mass-weighted mixing, enthalpy balancing, and unit conversion in
  Python. The solver path correctly uses `cb.mixer_from_streams_and_jacobians`.
  `mix_streams` is a legacy convenience function with no Jacobian.
- **Fix**: Deprecate or redirect to `cb.mixer_from_streams_and_jacobians`.

#### L2 — Legacy `combustion_from_streams()` Python-side post-processing
- **File**: `combustion/combustion.py:109-227`
- **Issue**: Efficiency scaling, pressure drop, and derived property computation done
  in Python without Jacobians. Solver's `CombustorNode.compute_derived_state` uses
  the proper C++ path.
- **Fix**: Deprecate or annotate as "standalone convenience, not for solver use".

#### L3 — `MixtureState.X` property can throw on invalid Y
- **File**: `components.py:354-361`
- **Issue**: `cb.mass_to_mole()` is called without a try/except. If Y values are
  somehow all-zero or cause a C++ error, this throws. The penalty wrapper in
  `residuals_wrapper` catches this, but it adds noise to solver diagnostics.
- **Fix**: Add try/except with fallback to equal-fraction default, or pre-validate Y.

#### L4 — Double `_propagate_states` at solver startup
- **File**: `solver.py:1261-1273`
- **Issue**: `_propagate_states(x0_use)` is called explicitly (line 1262), then
  immediately `_residuals(x0_use)` (line 1273) calls it again internally.
- **Fix**: Remove the standalone call; the `_residuals` call handles it.

#### L5 — `stoichiometry.py` and `combustion.py` raise on non-physical states
- **File**: `combustion/combustion.py:49-56,160-161`, `combustion/stoichiometry.py:81-91`
- **Issue**: `ValueError` raised on conditions like `m_total <= 0` or mismatched
  pressures that can occur during solver iteration. Currently not on the solver
  hot-path, but would break if wired into a solver node.
- **Fix**: Low priority — these are legacy functions not used by the solver.

### TRIVIAL

#### T1 — Hardcoded pi in `PipeElement.__init__`
- **File**: `components.py:1451`
- **Issue**: `self.area = 3.1415926535 * (diameter / 2) ** 2` — should use `math.pi`.
- **Fix**: Replace with `math.pi`.

---

## Part C — Thermal Wall Refactor: Stackable Multi-Layer Walls

### Current State

The current `WallConnection` is a **single-layer** coupler:
- Scalar `wall_thickness` [m] and `wall_conductivity` [W/(m·K)]
- Computes `t_over_k = thickness / conductivity` and passes to C++
- C++ `wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b, t_over_k, A)`
  uses `U = 1 / (1/h_a + t/k + 1/h_b)` — a single resistance layer

The convective HTC comes from `ConvectiveSurface.htc_and_T()` on elements that have
a `surface` attribute (`PipeElement`, `MomentumChamberNode`, `CombustorNode`). The
base `NetworkElement.htc_and_T()` returns `None`.

### Desired State

1. **Stackable wall layers** — e.g., TBC → bond coat → metal substrate, each with its
   own thickness and conductivity. The combined thermal resistance is
   `R_wall = Σ(t_i / k_i)`.
2. **Walls as network building blocks** — analogous to flow elements that connect
   between nodes. A wall connects a hot-side element to a cold-side element.
3. **Convective HTC only from physically meaningful elements** — only elements that
   have transport states and velocity (pipe, orifice, momentum chamber, combustor)
   contribute convective boundary conditions. Plenums and boundaries do not.

### C++ Already Supports Multi-Layer

The C++ backend already has the multi-layer overloads:

```cpp
// heat_transfer.h
double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double>& t_over_k_layers);

double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double>& t_over_k_layers,
                        double R_fouling);

std::vector<double> wall_temperature_profile(
    double T_hot, double T_cold, double h_hot, double h_cold,
    const std::vector<double>& t_over_k, double& q);
```

However, `wall_coupling_and_jacobian` currently takes a **scalar** `t_over_k`, not a
vector. This needs a multi-layer overload.

### Refactor Plan

#### C.1 — Replace `WallConnection` with `WallLayer` + `ThermalWall`

```python
@dataclass
class WallLayer:
    """Single material layer in a composite wall."""
    thickness: float          # [m]
    conductivity: float       # [W/(m·K)]
    label: str = ""           # e.g. "TBC", "bond_coat", "IN718"

    @property
    def t_over_k(self) -> float:
        return self.thickness / self.conductivity


@dataclass
class ThermalWall:
    """Multi-layer wall connecting two convective elements.

    Layers are ordered hot-side (element_a) to cold-side (element_b).
    """
    id: str
    element_a: str            # hot-side element/node ID
    element_b: str            # cold-side element/node ID
    layers: list[WallLayer]   # ordered hot → cold
    contact_area: float | None = None   # override [m^2]
    R_fouling: float = 0.0    # additional fouling [m^2·K/W]

    @property
    def t_over_k_total(self) -> float:
        return sum(layer.t_over_k for layer in self.layers) + self.R_fouling
```

**Migration**: `WallConnection(id, a, b, thickness=t, conductivity=k)` becomes
`ThermalWall(id, a, b, layers=[WallLayer(t, k)])`. Backward-compatible factory or
`from_single_layer()` classmethod.

#### C.2 — Extend C++ `wall_coupling_and_jacobian` for multi-layer

Add an overload that accepts `std::vector<double> t_over_k_layers`:

```cpp
WallCouplingResult wall_coupling_and_jacobian(
    double h_a, double T_aw_a,
    double h_b, double T_aw_b,
    const std::vector<double>& t_over_k_layers,
    double A,
    double R_fouling = 0.0);
```

Internally: `R_wall = Σ(t_over_k_layers) + R_fouling`. The Jacobian structure is
unchanged — `dQ/dh_a`, `dQ/dh_b`, `dQ/dT_aw_a`, `dQ/dT_aw_b` depend only on
`R_total = 1/h_a + R_wall + 1/h_b`, which is a scalar sum regardless of number of
layers.

**Additionally**: Add `T_wall_profile` output (vector of interface temperatures) for
post-solve diagnostics, using the existing `wall_temperature_profile` C++ function.

#### C.3 — Restrict convective HTC to transport-capable components

Currently `NetworkElement.htc_and_T()` returns `None` by default and only `PipeElement`,
`MomentumChamberNode`, and `CombustorNode` override it. This already provides the
correct restriction — elements without transport states return `None` and walls skip
them.

The current design is correct but the check is implicit (duck-typing on `htc_and_T`
returning `None`). To make this explicit:

```python
class NetworkElement(ABC):
    @property
    def has_convective_surface(self) -> bool:
        """Whether this element can provide convective HTC."""
        return hasattr(self, 'surface') and self.surface.area > 0.0
```

This clarifies intent without breaking the existing wall solver logic.

#### C.4 — Solver integration changes

The solver's `_evaluate_walls_for_node` currently computes `t_over_k` as:

```python
t_over_k = wall.wall_thickness / wall.wall_conductivity
```

With multi-layer walls, this becomes:

```python
t_over_k_layers = [layer.t_over_k for layer in wall.layers]
# Call new C++ overload
wall_result = cb.wall_coupling_and_jacobian(
    h_a, T_aw_a, h_b, T_aw_b, t_over_k_layers, A_eff, wall.R_fouling
)
```

The Jacobian relay structure (`dQ/dh_a * dh_a/dmdot`, etc.) is **unchanged** because
the wall's internal layer count does not affect the outer Jacobian chain — only the
scalar `R_total` matters.

#### C.5 — Post-solve: wall temperature profile

After convergence, call `cb.wall_temperature_profile(T_hot, T_cold, h_hot, h_cold,
t_over_k_layers)` once per wall to get interface temperatures. Add to `sol_dict`:

```
{wall_id}.T_interface[0]   # hot-surface temp
{wall_id}.T_interface[1]   # TBC/bond_coat interface
{wall_id}.T_interface[2]   # bond_coat/metal interface
{wall_id}.T_interface[3]   # cold-surface temp
```

#### C.6 — GUI schema additions

```python
class WallLayerData(BaseModel):
    thickness: float         # [m]
    conductivity: float      # [W/(m·K)]
    label: str = ""

class ThermalWallData(BaseModel):
    layers: list[WallLayerData]
    contact_area: float | None = None
    R_fouling: float = 0.0
```

### Summary of Changes

| Component | Change | Effort |
|-----------|--------|--------|
| `components.py` | Add `WallLayer`, rename `WallConnection` → `ThermalWall`, add `layers` list | Small |
| `graph.py` | Update `add_wall` / `from_dict` / `to_dict` for new schema | Small |
| C++ `heat_transfer.h/cpp` | Add multi-layer `wall_coupling_and_jacobian` overload | Small (sum reduction) |
| C++ `_core.cpp` (pybind) | Expose new overload | Small |
| `solver.py` | Update `_evaluate_walls_for_node` to pass layer vector | Small |
| `solver.py` | Add wall temperature profile to post-solve | Small |
| `schemas.py` | Add `WallLayerData`, `ThermalWallData` | Small |
| `graph_builder.py` | Parse wall layers from GUI schema | Small |
| Tests | Update wall coupling tests, add multi-layer test | Medium |
| Backward compat | Factory `ThermalWall.from_single_layer(id, a, b, t, k)` | Small |

**Total effort**: Medium. The C++ multi-layer infrastructure already exists. The main
work is the Python dataclass refactor and wiring through the solver and GUI. The
Jacobian structure is unchanged because layer count only affects the scalar `R_total`.

---

## Recommended Implementation Order

### Phase 1 — Solver efficiency (no API change)
1. [x] **H1 + L4** — Eliminate redundant evaluations in solver (Completed in 4e1b2d0, fix in 1340974)
2. [x] **M1** — Add `_residuals_only` method (Completed in 837cda2)

### Phase 2 — GUI state extraction
3. [x] **H2 + H3 + M5** — Refactor GUI to use solver results directly (Completed in 13d2842)
4. [x] **M2** — Integrate `extract_complete_states` into post-solve path (Completed in 13d2842, fix in 1c2c2e3)

### Phase 3 — Robustness fixes
5. [x] **M4** — Add `rho > 0` guard to `PipeElement.htc_and_T` (Completed in dafa6d3, polish in 045763c)
6. [x] **T1** — Standardize Pi constants (Completed in dafa6d3)
7. [x] **L3** — Add defensive guard to `MixtureState.X` (Completed in dafa6d3)

### Phase 4 — Thermal wall refactor (Part C)
8. [x] **C.1** — `WallLayer` + `ThermalWall` Python dataclasses (Completed in eb3dd21, fix in cc249c2)
9. [x] **C.2** — C++ multi-layer `wall_coupling_and_jacobian` overload + pybind (Completed in eb3dd21, fix in cc249c2)
10. [x] **C.4** — Solver wiring for multi-layer walls (Completed in eb3dd21, fix in cc249c2)
11. [x] **C.5** — Post-solve wall temperature profile (Completed in eb3dd21, fix in cc249c2)
12. [x] **C.3** — Explicit `has_convective_surface` property (Completed in 2b0289e, fix in aa94c32)
13. [x] **C.6** — GUI schema + graph builder for wall layers (Completed in eb3dd21, fix in cc249c2)

### Phase 5 — Tuning factor completeness (Part E)
14. [x] **E.2.1** — Add `Nu_multiplier`/`f_multiplier` to `heat_transfer.py` wrappers (Completed in 1878894)
15. [x] **E.2.2** — Wire multipliers through GUI schemas + graph builder (Completed in 1878894)

### Phase 6 — Thermal wall refactor advanced (Parts C + F)
16. **F.2.1** — `FilmCoolingModifier` on `ThermalWall` hot side
17. **F.2.2** — `EffusionCoolingModifier` (T_aw + flow coupling)
18. **F.3** — Expose non-smooth channel models in GUI (schema + builder)

### Phase 7 — Missing pybind wrappers (Part G)
19. **G.2b** — Bind `calc_T_from_s_mass` (isentropic process solver — **unblocks H**)
20. **G.2c-d** — Bind `calc_T_from_u_mass`, `calc_T_from_sv_mass`, `calc_T_from_sh_mass`
21. **G.2a** — Bind `dh_dT`, `ds_dT`, `dcp_dT` (thermo derivatives)
22. **G.8a** — Wrap `get_material` for generic k(T) lookup by name

### Phase 8 — Turbomachinery stages (Part H)
23. **H.1a** — C++ `compressor_stage_pr_and_jacobian` + `turbine_stage_pr_and_jacobian`
24. **H.1b** — C++ `compressor_stage_capacity_and_jacobian` + `turbine_stage_capacity_and_jacobian`
25. **H.2** — C++ central-FD accuracy tests for all 4 functions (8 test cases, tol < 1e-6)
26. **H.3** — Pybind bindings for `StageResult` + `StageCapacityResult` structs
27. **H.4** — Python `CompressorStageElement` + `TurbineStageElement` (both modes)
28. **H.5** — C++ `shaft_balance_residual_and_jacobian` + central-FD test
29. **H.6** — Python `ShaftConnection`
30. **H.7** — GUI schemas + graph builder (mode dropdown)

### Phase 9 — Cleanup
31. **M3** — Document or fix auto_Cd Jacobian truncation
32. **L1 + L2** — Deprecate legacy combustion functions
33. **L5** — Guard legacy combustion raises (only if they become solver-facing)

---

## Architectural Note: Ideal Post-Solve Flow

After all fixes, the post-solve data flow should be:

```
solve()
  -> root() converges -> final_x
  -> _propagate_states(final_x)          # once, populates _derived_states
  -> for each node:
       cb.complete_state(T, P, X)        # once per node, returns CompleteState
  -> for each element:
       cb.complete_state(T_in, P_in, X)  # once per element
  -> sol_dict contains all unknowns + derived T/Y + CompleteState fields
  -> GUI reads sol_dict directly — zero re-computation
```

---

## Part D — PyPI Publishing Readiness

Current state assessed against [PyPA packaging guide](https://packaging.python.org/)
and [scikit-build-core best practices](https://scikit-build-core.readthedocs.io/).

### D.1 — pyproject.toml metadata (HIGH — PyPI rejects without this)

The current `pyproject.toml` is **minimal**. PyPI requires or strongly recommends:

```toml
[project]
name = "combaero"
version = "0.2.0"
description = "Combustion and aerothermal tools (C++ core with Python bindings)"
readme = "README.md"
requires-python = ">=3.9"
license = "MIT"
authors = [
    { name = "...", email = "..." },
]
keywords = ["combustion", "thermodynamics", "heat-transfer", "gas-turbine", "CFD"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Programming Language :: C++",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Topic :: Scientific/Engineering",
]
dependencies = ["numpy>=1.24", "scipy>=1.10"]

[project.urls]
Homepage = "https://github.com/.../combaero"
Documentation = "https://..."
Repository = "https://github.com/.../combaero"
Changelog = "https://github.com/.../combaero/blob/main/CHANGELOG.md"
```

**Missing today**: `license`, `authors`, `keywords`, `classifiers`, `[project.urls]`.
PyPI will accept without them, but the listing will look empty/unprofessional.

### D.2 — Cross-platform wheel builds with cibuildwheel (HIGH)

The current CI builds C++ and runs tests, but **does not build wheels**. For PyPI you
need pre-built wheels for at least:

- **Linux**: x86_64 (manylinux), aarch64
- **macOS**: x86_64, arm64 (Apple Silicon)
- **Windows**: x86_64

[cibuildwheel](https://cibuildwheel.readthedocs.io/) is the standard tool. Add to
`pyproject.toml`:

```toml
[tool.cibuildwheel]
build = "cp39-* cp310-* cp311-* cp312-* cp313-*"
skip = "*-musllinux_*"
test-requires = "pytest"
test-command = "pytest {project}/python/tests/ -x"

[tool.cibuildwheel.linux]
archs = ["x86_64", "aarch64"]
before-build = ""  # cmake is included in manylinux images

[tool.cibuildwheel.macos]
archs = ["x86_64", "arm64"]

[tool.cibuildwheel.windows]
archs = ["AMD64"]
```

Add a CI workflow (e.g. `.github/workflows/publish.yml`) triggered on tags:

```yaml
name: Publish
on:
  push:
    tags: ["v*"]
jobs:
  build-wheels:
    uses: pypa/cibuildwheel@v2
  publish:
    needs: build-wheels
    uses: pypa/gh-action-pypi-publish@release/v1
```

### D.3 — Type stubs / `py.typed` marker (MEDIUM)

The package has no `py.typed` marker file and no `.pyi` type stubs. For a
scientific library with a pybind11 C++ extension, this means:

- Users get no autocomplete or type checking for `combaero._core` functions
- `mypy` / `pyright` cannot validate calls

**Fix**:
- Add `python/combaero/py.typed` (empty marker file)
- Generate stubs with [pybind11-stubgen](https://github.com/sizmailov/pybind11-stubgen):
  `pybind11-stubgen combaero._core -o python/combaero/`
- Include the `.pyi` files in the wheel via `scikit-build` config

### D.4 — Version management (MEDIUM)

Version is hardcoded in `pyproject.toml` (`version = "0.2.0"`). For PyPI releases,
**every upload must have a unique version**. Options:

- **Manual**: bump in `pyproject.toml` before each release (current approach — fragile)
- **`setuptools-scm`** or **`scikit-build-core`'s built-in git versioning**:
  derive version from git tags automatically

Recommended: use `scikit-build-core`'s `[tool.scikit-build.metadata.version]` with
git tags:

```toml
[project]
dynamic = ["version"]

[tool.scikit-build.metadata]
version.provider = "scikit_build_core.metadata.setuptools_scm"

[tool.setuptools_scm]
```

Or simpler: keep manual but add a CI check that the tag matches `pyproject.toml`.

### D.5 — CHANGELOG.md (MEDIUM)

No changelog exists. PyPI best practice is a `CHANGELOG.md` tracking releases.
This is also linked from `[project.urls]`. Format:
[Keep a Changelog](https://keepachangelog.com/).

### D.6 — Source distribution (sdist) includes build files (LOW)

Without a `MANIFEST.in` or explicit sdist config, `scikit-build-core` includes
everything. This means the sdist may contain `build/`, `.venv/`, `gui/`,
`diagnostic_plots/`, etc.

Add to `pyproject.toml`:

```toml
[tool.scikit-build]
sdist.include = ["src/**", "include/**", "python/**", "CMakeLists.txt",
                 "README.md", "LICENSE", "pyproject.toml"]
sdist.exclude = ["build", ".venv", "gui", "diagnostic_*", "benchmarks",
                 "cantera_validation_tests", "thermo_data_generator", "dist"]
```

### D.7 — `__all__` cleanup in `__init__.py` (LOW)

The `__init__.py` exports **~300 symbols** in `__all__`. This is fine but:
- Two different `MixtureState` classes exist: `combaero.MixtureState` (C++ pybind11)
  and `combaero.network.MixtureState` (Python dataclass). This will confuse users.
- Many internal solver helpers are exported (`orifice_mdot_and_jacobian`,
  `density_and_jacobians`, etc.) — consider whether these belong in the public API
  or under a `combaero._solver_tools` namespace.

**Recommendation**: Audit `__all__` and split into public API vs internal helpers.
Resolve the `MixtureState` name collision.

### D.8 — Network/GUI should not be in the PyPI package (LOW)

The `network/` subpackage and `gui/` are application-layer code, not library code.
For a PyPI release:
- `combaero` (the library) goes on PyPI
- `combaero.network` could be a separate `combaero-network` package or stay bundled
  as an optional extra
- `gui/` should **never** be in the PyPI wheel — it's a standalone application

Current `wheel.packages = ["python/combaero"]` includes `network/` since it's a
subpackage. Options:
1. Keep it bundled (simpler, what numpy/scipy do with sub-modules)
2. Split `combaero-network` as a separate distribution with `depends = ["combaero"]`

**Recommendation**: Keep bundled for now (option 1), but ensure `gui/` is excluded
from wheels.

### D.9 — Python version support floor (LOW)

`requires-python = ">=3.9"` — Python 3.9 is EOL (Oct 2025). For a new PyPI package,
`>=3.10` is reasonable. `>=3.11` if you want to drop old manylinux CI complexity.
The `ruff.toml` already targets `py312`.

### D.10 — License file in wheel (TRIVIAL)

The `LICENSE` file exists but isn't explicitly included in the wheel metadata. Add:

```toml
[project]
license = {file = "LICENSE"}
```

Or with the newer spec:
```toml
license = "MIT"
license-files = ["LICENSE"]
```

### Summary Table

| # | Item | Severity | Effort | Status |
|---|------|----------|--------|--------|
| D.1 | pyproject.toml metadata | High | Small | Missing authors, classifiers, urls |
| D.2 | cibuildwheel CI for multi-platform wheels | High | Medium | No wheel CI exists |
| D.3 | Type stubs + py.typed | Medium | Small | No stubs, no marker |
| D.4 | Version management (git tags) | Medium | Small | Hardcoded version |
| D.5 | CHANGELOG.md | Medium | Small | Does not exist |
| D.6 | sdist include/exclude | Low | Small | May ship build artifacts |
| D.7 | `__all__` audit + MixtureState collision | Low | Medium | 300+ exports, name clash |
| D.8 | network/gui packaging scope | Low | Decision | Keep bundled vs split |
| D.9 | Python version floor | Low | Trivial | 3.9 is EOL |
| D.10 | License in wheel metadata | Trivial | Trivial | Missing `license-files` |

---

## Part E — Tuning Factor Wiring Audit

**Goal**: Ensure `Nu_multiplier` and `f_multiplier` are wired end-to-end from GUI
through network components into C++ channel functions, and that they participate in
Jacobians correctly.

### E.1 — Wiring Status by Layer

| Layer | Nu_multiplier | f_multiplier | Notes |
|-------|:---:|:---:|-------|
| C++ `channel_smooth` | **YES** | **YES** | Default 1.0, applied inside correlation |
| C++ `channel_ribbed` | **YES** | **YES** | Applied on top of rib enhancement |
| C++ `channel_dimpled` | **YES** | **YES** | Applied on top of dimple enhancement |
| C++ `channel_pin_fin` | **YES** | **YES** | Applied on top of pin-fin correlation |
| C++ `channel_impingement` | **YES** | **YES** | Applied on top of impingement Nu |
| pybind `_core.cpp` | **YES** | **YES** | All 5 channel functions bound with both params |
| Python `ConvectiveSurface.htc_and_T()` | **YES** | **YES** | Passed from `self.Nu_multiplier`/`self.f_multiplier` |
| Python `heat_transfer.py` wrappers | **NO** | **NO** | `smooth()`, `ribbed()`, `dimpled()`, `pin_fin()`, `impingement()` do not accept multipliers |
| `PipeElement.htc_and_T()` | **YES** | **YES** | Delegates to `self.surface.htc_and_T()` |
| `MomentumChamberNode.htc_and_T()` | **YES** | **YES** | Delegates to `self.surface.htc_and_T()` |
| `CombustorNode.htc_and_T()` | Implicit | Implicit | Uses `self.surface` which carries multipliers |
| Solver `_evaluate_walls_for_node` | **YES** | **YES** | Calls `obj.htc_and_T(state)` → `ConvectiveSurface` |
| GUI `schemas.py` | **NO** | **NO** | No multiplier fields on any element schema |
| GUI `graph_builder.py` | **NO** | **NO** | Hardcodes `ConvectiveSurface(area=conv_area)` with defaults |
| GUI frontend | **NO** | **NO** | No UI controls for multipliers |

### E.2 — Findings

#### E.2.1 — `heat_transfer.py` wrappers drop multipliers (MEDIUM)

The high-level Python wrappers (`smooth()`, `ribbed()`, `dimpled()`, `pin_fin()`,
`impingement()`) in `heat_transfer.py` do **not** accept or forward `Nu_multiplier`
or `f_multiplier`. These are standalone convenience functions (not on the solver path),
but any user calling e.g. `combaero.heat_transfer.smooth(...)` cannot apply tuning.

```python
# heat_transfer.py line 96-108
def smooth(...) -> ChannelResult:
    return _channel_smooth(
        T, P, list(X), u, D, L,
        T_wall=T_wall, correlation=correlation,
        heating=heating, mu_ratio=mu_ratio, roughness=roughness,
        # Nu_multiplier and f_multiplier NOT passed
    )
```

**Fix**: Add `Nu_multiplier: float = 1.0` and `f_multiplier: float = 1.0` keyword
args to all 5 wrapper functions and forward them.

#### E.2.2 — GUI does not expose multipliers (LOW — future work)

The GUI `PipeData`/`OrificeData` schemas and `graph_builder.py` do not wire
`Nu_multiplier`/`f_multiplier`. The `ConvectiveSurface` is constructed with defaults:

```python
# graph_builder.py line 161
surface=ConvectiveSurface(area=conv_area)
# → Nu_multiplier=1.0, f_multiplier=1.0 (defaults)
```

This is acceptable for the initial GUI (no heat transfer tuning exposed yet). When
thermal coupling is fully exposed in the GUI, add:

```python
class PipeData(BaseModel):
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
```

And in `graph_builder.py`:
```python
surface=ConvectiveSurface(
    area=conv_area,
    Nu_multiplier=data.Nu_multiplier,
    f_multiplier=data.f_multiplier,
)
```

#### E.2.3 — Jacobian correctness of multipliers

In C++, the multipliers are applied as `h = Nu_multiplier * h_base` and
`f = f_multiplier * f_base`. The `ChannelResult` Jacobians (`dh_dmdot`, `dh_dT`)
are computed **after** the multiplier is applied, so the chain rule is correct:

```
dh/dmdot = Nu_multiplier * dh_base/dmdot
```

This is the correct design (multiplier inside C++, not Python post-correction) per
the design doc. **No Jacobian issue.**

#### E.2.4 — Channel model selection not exposed in GUI (LOW)

The GUI hardcodes `SmoothModel()` for pipes. Users cannot select ribbed, dimpled,
pin-fin, or impingement models from the GUI. This is future GUI work (Part F overlap).

### E.3 — Summary

The multiplier wiring is **correct through the solver path** (Python network →
`ConvectiveSurface` → C++ `channel_*`). The two gaps are:

1. **`heat_transfer.py` convenience wrappers** don't forward multipliers (MEDIUM fix)
2. **GUI** doesn't expose multipliers yet (LOW, future GUI feature)

---

## Part F — Missing Correlations: C++ Available but Not Wired to Network

### F.1 — Inventory: C++ Cooling Functions

The C++ library in `cooling_correlations.h` provides these standalone functions:

| C++ Function | Pybind | Network Model | GUI |
|---|:---:|:---:|:---:|
| `rib_enhancement_factor` | YES | via `channel_ribbed` | NO |
| `rib_enhancement_factor_high_re` | YES | via `channel_ribbed` | NO |
| `rib_friction_multiplier` | YES | via `channel_ribbed` | NO |
| `rib_friction_multiplier_high_re` | YES | via `channel_ribbed` | NO |
| `thermal_performance_factor` | YES | — (post-processing) | NO |
| `impingement_nusselt` | YES | via `channel_impingement` | NO |
| `film_cooling_effectiveness` | YES | **NO network model** | NO |
| `film_cooling_effectiveness_avg` | YES | **NO network model** | NO |
| `film_cooling_multirow_sellers` | YES | **NO network model** | NO |
| `effusion_effectiveness` | YES | **NO network model** | NO |
| `effusion_discharge_coefficient` | YES | **NO network model** | NO |
| `pin_fin_nusselt` | YES | via `channel_pin_fin` | NO |
| `pin_fin_friction` | YES | via `channel_pin_fin` | NO |
| `dimple_nusselt_enhancement` | YES | via `channel_dimpled` | NO |
| `dimple_friction_multiplier` | YES | via `channel_dimpled` | NO |
| `adiabatic_wall_temperature` | YES | **NO network model** | NO |
| `cooled_wall_heat_flux` | YES | **NO network model** | NO |

And `solver_interface.h` Jacobian-equipped versions:

| C++ Function | Pybind | Used in Network |
|---|:---:|:---:|
| `pin_fin_nusselt_and_jacobian` | YES | Inside `channel_pin_fin` |
| `pin_fin_friction_and_jacobian` | YES | Inside `channel_pin_fin` |
| `dimple_nusselt_enhancement_and_jacobian` | YES | Inside `channel_dimpled` |
| `dimple_friction_multiplier_and_jacobian` | YES | Inside `channel_dimpled` |
| `rib_enhancement_factor_high_re_and_jacobian` | YES | Inside `channel_ribbed` |
| `impingement_nusselt_and_jacobian` | YES | Inside `channel_impingement` |
| `film_cooling_effectiveness_and_jacobian` | YES | **NOT used in network** |
| `effusion_effectiveness_and_jacobian` | YES | **NOT used in network** |

### F.2 — Missing Network Models (Not Yet Wired)

Five C++ correlation families are exposed via pybind but have **no** corresponding
`ChannelModel` subclass in the network:

#### F.2.1 — Film Cooling

- `film_cooling_effectiveness(x_D, M, DR, alpha_deg)` → local eta
- `film_cooling_effectiveness_avg(x_D, M, DR, alpha_deg, s_D)` → laterally averaged
- `film_cooling_multirow_sellers(row_positions, eval_xD, M, DR, alpha_deg)` → multi-row

These modify the *adiabatic wall temperature* (`T_aw`) seen by the hot-side wall, not
the convective HTC. They are typically used as **boundary condition modifiers** rather
than channel-flow models.

**Network model**: A `FilmCoolingModel` could modify `T_aw` on the hot side:
```
T_aw_film = T_hot - eta * (T_hot - T_coolant)
```
This would sit on `WallConnection` as a hot-side modifier, not as a `ChannelModel`.

#### F.2.2 — Effusion Cooling

- `effusion_effectiveness(x_D, M, DR, porosity, s_D, alpha_deg)` → eta for dense arrays
- `effusion_discharge_coefficient(Re_d, P_ratio, alpha_deg, L_D)` → Cd for effusion holes

Similar to film cooling — modifies `T_aw` on the hot side. Additionally, the mass
flow through effusion holes couples the coolant supply pressure to the hot-side
pressure (like an orifice array).

**Network model**: More complex — requires coupling coolant mass flow (through
effusion holes) to both thermal and flow networks. Could be modeled as:
1. Orifice array for flow coupling (using `effusion_discharge_coefficient`)
2. `T_aw` modifier for thermal coupling (using `effusion_effectiveness`)

#### F.2.3 — Standalone Helpers

- `adiabatic_wall_temperature(T_hot, T_coolant, eta)` — convenience for T_aw from eta
- `cooled_wall_heat_flux(T_hot, T_coolant, h_hot, h_coolant, eta, t_wall, k_wall)` —
  single-pass combined heat flux
- `thermal_performance_factor(Nu_ratio, f_ratio)` — post-processing diagnostic

These are building blocks, not channel models. They don't need network model classes.

### F.3 — Channel Models Wired but Not Exposed in GUI

All 5 `ChannelModel` subclasses exist in Python but the GUI only uses `SmoothModel`:

| Model | Python Class | GUI Schema | GUI Frontend |
|-------|-------------|:---:|:---:|
| Smooth | `SmoothModel` | Implicit only | NO (no controls) |
| Ribbed | `RibbedModel` | **NO** | NO |
| Dimpled | `DimpledModel` | **NO** | NO |
| Pin-fin | `PinFinModel` | **NO** | NO |
| Impingement | `ImpingementModel` | **NO** | NO |

### F.4 — Recommended Implementation for Missing Models

#### Phase F-I: Film/Effusion as T_aw modifiers (new concept)

Film and effusion cooling are not channel-flow models — they modify the adiabatic wall
temperature seen by the `WallConnection`. This requires:

1. Add a `CoolingModifier` concept to `ThermalWall` (Part C refactor):
```python
@dataclass
class FilmCoolingModifier:
    """Modifies hot-side T_aw via film cooling effectiveness."""
    x_D: float           # downstream distance / hole diameter
    M: float             # blowing ratio
    DR: float            # density ratio
    alpha_deg: float     # injection angle [deg]
    s_D: float = 3.0     # hole spacing / diameter
    multirow: bool = False
    row_positions_xD: list[float] | None = None

@dataclass
class EffusionCoolingModifier:
    """Modifies hot-side T_aw via effusion cooling effectiveness."""
    x_D: float
    M: float
    DR: float
    porosity: float
    s_D: float
    alpha_deg: float
```

2. In `ThermalWall.compute_coupling()`:
```python
if self.hot_side_modifier:
    eta = cb.film_cooling_effectiveness(...)
    T_aw_a = cb.adiabatic_wall_temperature(T_aw_a, T_coolant, eta)
```

3. The Jacobian chain extends: `dQ/dT_aw_a * dT_aw_a/deta * deta/dM`.
   The C++ `film_cooling_effectiveness_and_jacobian` already provides `deta/dM`.

#### Phase F-II: Expose channel model selection in GUI

Add `htc_model` selection to pipe/element schemas:
```python
class PipeData(BaseModel):
    htc_model: str = "smooth"  # "smooth"|"ribbed"|"dimpled"|"pin_fin"|"impingement"
    htc_params: dict = {}      # model-specific parameters
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
```

Map to the correct `ChannelModel` subclass in `graph_builder.py`.

### F.5 — Summary Table

| # | Item | Severity | Effort | Status |
|---|------|----------|--------|--------|
| E.2.1 | `heat_transfer.py` wrappers missing multipliers | Medium | Small | 5 functions to fix |
| E.2.2 | GUI multiplier exposure | Low | Small | Future GUI feature |
| E.2.4 | GUI channel model selection | Low | Medium | Future GUI feature |
| F.2.1 | Film cooling network model | Medium | Medium | New `CoolingModifier` concept needed |
| F.2.2 | Effusion cooling network model | Medium | Medium | Flow+thermal coupling |
| F.3 | Non-smooth models in GUI | Low | Small | Schema + builder wiring |

---

## Part G — Full C++ → Pybind → Python → GUI Wiring Audit

Systematic screen of every public C++ header for functions/types that are:
- **Not bound** in pybind (`_core.cpp`)
- **Bound** but not re-exported in `__init__.py`
- **Exported** in Python but not wired to network/GUI

Assessment criterion: *Would an engineering student or practising engineer benefit
from having this exposed in the Python API or GUI?*

### G.1 — `solver_interface.h`

All solver-interface functions are internal building blocks consumed by the Python
network solver. They are **correctly not exposed** as public API — users should not
call `orifice_residuals_and_jacobian` directly. The only user-facing items are the
Jacobian-equipped correlation functions, all of which are already bound.

| Function | Pybind | `__init__` | Network | GUI | Assessment |
|----------|:---:|:---:|:---:|:---:|------------|
| `orifice_mdot_and_jacobian` | YES | YES | Solver-internal | — | Correct: internal |
| `orifice_residuals_and_jacobian` | YES | YES | Solver-internal | — | Correct: internal |
| `pipe_residuals_and_jacobian` | YES | YES | Solver-internal | — | Correct: internal |
| `lossless_pressure_and_jacobian` | YES | YES | Solver-internal | — | Correct: internal |
| `orifice_compressible_*` | YES | YES | Solver-internal | — | Correct: internal |
| `pipe_compressible_*` | YES | YES | Solver-internal | — | Correct: internal |
| `nusselt_and_jacobian_*` (4 variants) | YES | YES | — | — | OK: useful for advanced users |
| `friction_and_jacobian_*` (5 variants) | YES | YES | — | — | OK: useful for advanced users |
| `pin_fin_nusselt_and_jacobian` | YES | YES | — | — | OK |
| `pin_fin_friction_and_jacobian` | YES | YES | — | — | OK |
| `dimple_*_and_jacobian` | YES | YES | — | — | OK |
| `rib_enhancement_factor_high_re_and_jacobian` | YES | YES | — | — | OK |
| `impingement_nusselt_and_jacobian` | YES | YES | — | — | OK |
| `film_cooling_effectiveness_and_jacobian` | YES | YES | — | — | OK |
| `effusion_effectiveness_and_jacobian` | YES | YES | — | — | OK |
| `density_and_jacobians` | YES | YES | — | — | OK |
| `enthalpy_and_jacobian` | YES | YES | — | — | OK |
| `viscosity_and_jacobians` | YES | YES | — | — | OK |
| `mach_number_and_jacobian_v` | YES | YES | — | — | OK |
| `T_adiabatic_wall_and_jacobian_v` | YES | YES | — | — | OK |
| `T0_from_static_and_jacobian_M` | YES | YES | — | — | OK |
| `P0_from_static_and_jacobian_M` | YES | YES | — | — | OK |
| `adiabatic_T_complete_and_jacobian_T` | YES | YES | — | — | OK |
| `adiabatic_T_equilibrium_and_jacobians` | YES | YES | — | — | OK |
| `mixer_from_streams_and_jacobians` | YES | YES | Solver-internal | — | Correct |
| `adiabatic_T_*_from_streams` | YES | YES | Solver-internal | — | Correct |
| `combustor_residuals_and_jacobians` | YES | YES | Solver-internal | — | Correct |
| `momentum_chamber_residual_and_jacobian` | YES | YES | Solver-internal | — | Correct |

**Verdict: solver_interface.h — fully wired, no gaps.**

### G.2 — `thermo.h`

| Function | Pybind | `__init__` | Assessment |
|----------|:---:|:---:|------------|
| `cp`, `h`, `s`, `cv`, `u` (molar) | YES | YES | OK |
| `cp_mass`, `h_mass`, `s_mass`, `cv_mass`, `u_mass` | YES | YES | OK |
| `density`, `molar_volume`, `specific_gas_constant` | YES | YES | OK |
| `isentropic_expansion_coefficient`, `speed_of_sound` | YES | YES | OK |
| `air_properties` | YES | YES | OK |
| `thermo_state`, `complete_state` | YES | YES | OK |
| `calc_T_from_h`, `calc_T_from_s`, `calc_T_from_cp` | YES | YES | OK |
| `calc_T_from_h_mass` | YES | YES | OK |
| `dh_dT`, `ds_dT`, `dcp_dT`, `dg_over_RT_dT` | **NO** | **NO** | **MISSING** — useful for sensitivity analysis |
| `calc_T_from_u` (molar internal energy) | **NO** | **NO** | Low value (rare use case) |
| `calc_T_from_s_mass` | **NO** | **NO** | **MISSING** — isentropic process solver |
| `calc_T_from_u_mass` | **NO** | **NO** | **MISSING** — isochoric process solver |
| `calc_T_from_sv_mass` | **NO** | **NO** | **MISSING** — (s, v) flash |
| `calc_T_from_sh_mass` | **NO** | **NO** | **MISSING** — (s, h) flash |
| `nusselt_pipe` (State-based overload) | **NO** | **NO** | Low value (htc_pipe covers this) |
| State-based overloads (`cp(State&)` etc.) | **NO** | **NO** | Low value (internal convenience) |

**Key gaps:**
- **`dh_dT`, `ds_dT`, `dcp_dT`** — Useful for students studying thermodynamic
  sensitivity and for engineers doing uncertainty propagation. Medium value.
- **`calc_T_from_s_mass`** — Essential for isentropic expansion/compression
  calculations (turbine/compressor modelling). **High value for students.**
- **`calc_T_from_u_mass`** — Isochoric heating (e.g., constant-volume combustion).
  Medium value.
- **`calc_T_from_sv_mass`, `calc_T_from_sh_mass`** — State-point flash calculations.
  Medium value for advanced thermodynamic analysis.

### G.3 — `transport.h`

| Function | Pybind | `__init__` | Assessment |
|----------|:---:|:---:|------------|
| `viscosity` | YES | YES | OK |
| `thermal_conductivity` | YES | YES | OK |
| `prandtl` | YES | YES | OK |
| `kinematic_viscosity` | YES | YES | OK |
| `thermal_diffusivity` | YES | YES | OK |
| `reynolds` | YES | YES | OK |
| `peclet` | YES | YES | OK |
| `transport_state` | YES | YES | OK |
| `omega22`, `linear_interp` | **NO** | **NO** | Correct: internal implementation details |

**Verdict: transport.h — fully wired, no gaps.**

### G.4 — `friction.h`

| Function | Pybind | `__init__` | Assessment |
|----------|:---:|:---:|------------|
| `friction_haaland` | YES | YES | OK |
| `friction_serghides` | YES | YES | OK |
| `friction_colebrook` | YES | YES | OK |
| `friction_petukhov` | YES | YES | OK |
| `pipe_roughness` | YES | YES | OK |
| `standard_pipe_roughness` | YES | YES | OK |

**Verdict: friction.h — fully wired, no gaps.**

### G.5 — `orifice.h`

| Function | Pybind | `__init__` | Assessment |
|----------|:---:|:---:|------------|
| `OrificeGeometry` struct | YES | YES | OK |
| `OrificeState` struct | YES | YES | OK |
| `CdCorrelation` enum | YES | YES | OK |
| `Cd_sharp_thin_plate` | YES | YES | OK |
| `Cd_thick_plate` | YES | YES | OK |
| `Cd_rounded_entry` | YES | YES | OK |
| `Cd` (auto-select) | YES | YES | OK |
| `orifice::Cd_ReaderHarrisGallagher` | YES | YES | OK |
| `orifice::Cd_Stolz` | YES | YES | OK |
| `orifice::Cd_Miller` | YES | YES | OK |
| `orifice::thickness_correction` | YES | YES | OK |
| `orifice::Cd_rounded` | YES | YES | OK |
| `orifice::K_from_Cd` | YES | YES | OK |
| `orifice::Cd_from_K` | YES | YES | OK |
| `orifice_mdot` | YES | YES | OK |
| `orifice_dP` | YES | YES | OK |
| `orifice_Cd_from_measurement` | YES | YES | OK |
| `solve_orifice_mdot` | YES | YES | OK |
| `OrificeFlowResult` struct | YES | YES | OK |
| `orifice_flow` | YES | YES | OK |
| `orifice_flow_state` (convenience) | YES | YES | OK |
| `orifice_velocity_from_mdot` | YES | YES | OK |
| `orifice_area_from_beta` | YES | YES | OK |
| `beta_from_diameters` | YES | YES | OK |
| `orifice_Re_d_from_mdot` | YES | YES | OK |
| `expansibility_factor` | YES | YES | OK |
| `make_correlation` (factory) | **NO** | **NO** | Low value: polymorphic OO pattern, not Pythonic |
| `make_user_correlation` | **NO** | **NO** | Low value: users can pass Cd directly |
| `make_tabulated_correlation` | **NO** | **NO** | **MISSING** — useful for test-data driven Cd. Medium. |

**Verdict: orifice.h — nearly complete. `make_tabulated_correlation` is the only
notable gap for users with empirical test data.**

### G.6 — `heat_transfer.h`

| Function | Pybind | `__init__` | Assessment |
|----------|:---:|:---:|------------|
| `nusselt_dittus_boelter` | YES | YES | OK |
| `nusselt_gnielinski` (both overloads) | YES | YES | OK |
| `nusselt_sieder_tate` | YES | YES | OK |
| `nusselt_petukhov` (both overloads) | YES | YES | OK |
| `htc_from_nusselt` | YES | YES | OK |
| `overall_htc` | YES | YES | OK |
| `overall_htc_wall` (all overloads) | YES | YES | OK |
| `thermal_resistance` | YES | YES | OK |
| `thermal_resistance_wall` | YES | YES | OK |
| `lmtd`, `lmtd_counterflow`, `lmtd_parallelflow` | YES | YES | OK |
| `heat_rate`, `heat_flux`, `heat_transfer_area`, `heat_transfer_dT` | YES | YES | OK |
| `wall_temperature_profile` | YES | YES | OK |
| `heat_flux_from_T_at_edge` | YES | YES | OK |
| `heat_flux_from_T_at_depth` | YES | YES | OK |
| `bulk_T_from_edge_T_and_q` | YES | YES | OK |
| `dT_edge_dT_hot`, `dT_edge_dT_cold`, `dT_edge_dT_bulk`, `dT_edge_dq` | YES | YES | OK |
| `ntu`, `capacity_ratio` | YES | YES | OK |
| `effectiveness_counterflow`, `effectiveness_parallelflow` | YES | YES | OK |
| `heat_rate_from_effectiveness` | YES | YES | OK |
| `ChannelResult` struct | YES | YES | OK |
| `WallCouplingResult` struct | YES | YES | OK |
| `wall_coupling_and_jacobian` | YES | YES | OK |
| `channel_smooth/ribbed/dimpled/pin_fin/impingement` | YES | YES | OK |
| `htc_pipe` (composite tuple overload) | YES | YES | OK |
| `nusselt_pipe` (State-based) | **NO** | **NO** | Low value — `htc_pipe` covers the use case |
| `htc_pipe` (State-based) | YES | YES | OK |
| Laminar constants `NU_LAMINAR_CONST_T/Q` | **NO** | **NO** | Low: trivial constants (3.66, 4.36) |

**Verdict: heat_transfer.h — fully wired. No meaningful gaps.**

### G.7 — `combustion.h`

| Function | Pybind | `__init__` | Assessment |
|----------|:---:|:---:|------------|
| `oxygen_required_per_mol_fuel/kg_fuel` | YES | YES | OK |
| `oxygen_required_per_mol_mixture/kg_mixture` | YES | YES | OK |
| `equivalence_ratio` (elemental) | YES | YES | OK |
| `equivalence_ratio_mass` (elemental) | YES | YES | OK |
| `fuel_lhv_molar`, `fuel_lhv_mass` | YES | YES | OK |
| `dryair_required_per_mol_fuel/kg_fuel` | YES | YES | OK |
| `dryair_required_per_mol_mixture/kg_mixture` | YES | YES | OK |
| `equivalence_ratio_mole` (fuel/ox streams) | YES | YES | OK |
| `set_equivalence_ratio_mole/mass` | YES | YES | OK |
| `bilger_*` functions (all) | YES | YES | OK |
| `complete_combustion_to_CO2_H2O` (both) | YES | YES | OK |
| `combustion_state` (both overloads) | YES | YES | OK |
| `combustion_state_from_streams` (both) | YES | YES | OK |
| `set_fuel_stream_for_*` (5 variants) | YES | YES | OK |
| `set_oxidizer_stream_for_*` (5 variants) | YES | YES | OK |
| `complete_combustion` (State-based) | **NO** | **NO** | **MISSING** — convenience for quick State→State combustion |
| `complete_combustion_isothermal` (State) | **NO** | **NO** | Low value: niche use case |

**Verdict: combustion.h — nearly complete. `complete_combustion(State)` is a useful
convenience but not critical since `combustion_state` covers the same physics.**

### G.8 — `materials.h`

| Function | Pybind | `__init__` | Assessment |
|----------|:---:|:---:|------------|
| `k_inconel718` | YES | YES | OK |
| `k_haynes230` | YES | YES | OK |
| `k_stainless_steel_316` | YES | YES | OK |
| `k_aluminum_6061` | YES | YES | OK |
| `k_tbc_ysz` | YES | YES | OK |
| `list_materials` | YES | YES | OK |
| `get_material` (returns `ThermalMaterial`) | **NO** | **NO** | **MISSING** — generic k(T) lookup by name. Medium value. |
| `ThermalMaterial` struct | **NO** | **NO** | Would need wrapping of `std::function<double(double)>` |

**Verdict: materials.h — mostly wired. `get_material` would be valuable for generic
wall-layer specification by material name, but requires wrapping a C++ functor.**

### G.9 — Summary of All Missing Items

| # | Header | Function(s) | Usefulness | Effort | Notes |
|---|--------|------------|-----------|--------|-------|
| G.2a | thermo.h | `dh_dT`, `ds_dT`, `dcp_dT` | Medium | Small | Sensitivity analysis, uncertainty propagation |
| G.2b | thermo.h | `calc_T_from_s_mass` | **High** | Small | Isentropic T after expansion/compression |
| G.2c | thermo.h | `calc_T_from_u_mass` | Medium | Small | Isochoric process solver |
| G.2d | thermo.h | `calc_T_from_sv_mass`, `calc_T_from_sh_mass` | Medium | Small | State-point flash calculations |
| G.5a | orifice.h | `make_tabulated_correlation` | Medium | Medium | Test-data driven Cd(beta, Re) tables |
| G.7a | combustion.h | `complete_combustion(State)` | Low | Small | Quick State→State convenience |
| G.8a | materials.h | `get_material` (generic k(T) lookup) | Medium | Medium | Needs functor wrapping |

### G.10 — Recommended Priority

1. **`calc_T_from_s_mass`** — High value for turbomachinery students. Trivial to bind
   (same pattern as `calc_T_from_h_mass`).
2. **`calc_T_from_u_mass`**, **`calc_T_from_sv_mass`**, **`calc_T_from_sh_mass`** —
   Complete the inverse-solver family. Small effort each.
3. **`dh_dT`**, **`ds_dT`**, **`dcp_dT`** — Expose thermodynamic derivatives for
   educational use. Trivial pybind wrappers.
4. **`get_material`** — Useful when Part C (multi-layer walls) lands; allows
   specifying wall materials by name instead of hardcoded k(T).
5. **`make_tabulated_correlation`** — Useful for experimental validation workflows.

---

## Part H — Turbomachinery Stage Elements (Compressor / Turbine)

The missing `calc_T_from_s_mass` (G.2b) and related inverse thermo solvers are the
key C++ building blocks for adding compressor and turbine stage elements to the
network solver. This section sketches the design.

### H.1 — Thermodynamic Model

A single compressor or turbine stage with isentropic efficiency `eta_is`:

```
Given:  T1, P1, X   (inlet stagnation state from upstream node)
        PR           (pressure ratio P2/P1 — design parameter)
        eta_is       (isentropic efficiency — design parameter)

Step 1: P2 = P1 * PR
Step 2: T2s = calc_T_from_s_mass(s_mass(T1, X, P1), P2, X)   ← isentropic outlet
Step 3: h2  = h1 + (h2s - h1) / eta_is                        (compressor)
         h2  = h1 - eta_is * (h1 - h2s)                        (turbine)
Step 4: T2  = calc_T_from_h_mass(h2, X)                        ← actual outlet
Step 5: W   = mdot * (h2 - h1)                                 [W] shaft work
```

All five C++ functions (`s_mass`, `h_mass`, `calc_T_from_s_mass`, `calc_T_from_h_mass`,
`cp_mass`) already exist with proper Newton solvers. Only `calc_T_from_s_mass` is
missing from pybind (the single blocker for this feature).

### H.2 — Operating Modes

A turbomachinery stage needs three operating modes to cover different analysis types:

| Mode | Inputs (fixed) | Solved for | Use case |
|------|----------------|-----------|----------|
| **Prescribed PR** | `pressure_ratio`, `eta_is` | — | Single design-point "what-if" |
| **Capacity** | `throat_area` (or `capacity`), `eta_is` | `PR` emerges from network | System matching, gas turbine cycles |
| **Performance map** | `CompressorMap` / `TurbineMap` | `PR`, `eta_is` both from map | Off-design analysis |

### H.2a — Capacity Model (Key for System Matching)

**Problem with fixed PR**: if PR is prescribed, the mass flow is fully determined by
the rest of the network. The stage is passive — it can't "push back" or limit flow.
This prevents modelling compressor–turbine matching, surge, or choking.

**Capacity** is the standard turbomachinery coupling parameter:

```
Capacity ≡ ṁ·√(T₀_in) / P₀_in    [kg·√K / (s·Pa)]
```

For a **choked turbine nozzle**, capacity is fixed by throat geometry:
```
capacity_choked = A_throat · √(γ/R) · (2/(γ+1))^((γ+1)/(2(γ-1)))
```

For a **compressor**, capacity at design point is a characteristic of the machine.

In capacity mode, the element's residual constrains mass flow, not pressure ratio:
```
Residual:  ṁ·√(T₀_in) / P₀_in  -  capacity  =  0
```
The pressure ratio is NOT prescribed — it emerges from the downstream node pressures
set by the rest of the network. The element observes the actual PR = P₀_out / P₀_in
and uses it to compute T_out via the isentropic model.

This is architecturally identical to how `OrificeElement` constrains ṁ from ΔP via
`ṁ = Cd·A·√(2ρΔP)` — the orifice doesn't prescribe ΔP, it relates ṁ to the
pressure difference the network provides.

### H.2b — Network Element Design

Following the existing `NetworkElement` pattern (cf. `PipeElement`, `OrificeElement`):

```python
class CompressorStageElement(NetworkElement):
    """Single compressor stage with isentropic efficiency."""

    def __init__(self, id, from_node, to_node,
                 eta_is: float = 0.85,
                 # --- Mode selection (exactly one must be set) ---
                 pressure_ratio: float | None = None,   # Mode 1: fixed PR
                 capacity: float | None = None,          # Mode 2: ṁ√T/P capacity
                 throat_area: float | None = None,       # Mode 2 alt: compute capacity from geometry
                 # performance_map: CompressorMap | None = None,  # Mode 3: future
                 shaft_id: str | None = None):
        ...

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def residuals(self, state_in, state_out):
        if self.pressure_ratio is not None:
            # Mode 1: Prescribed PR — residual enforces pressure ratio
            P2 = state_in.P_total * self.pressure_ratio
            res = [state_out.P_total - P2]
            jac = {0: {
                f"{self.from_node}.P_total": -self.pressure_ratio,
                f"{self.to_node}.P_total": 1.0,
            }}
        else:
            # Mode 2: Capacity — residual enforces mass flow capacity
            # C++ provides: capacity_residual_and_jacobian(m_dot, T1, P1, X, capacity)
            res_cpp = cb.compressor_capacity_residuals_and_jacobian(
                state_in.m_dot, state_in.T, state_in.P_total, state_in.X,
                self._capacity,
            )
            res = [res_cpp.residual]
            jac = {0: {
                f"{self.id}.m_dot": res_cpp.d_res_d_mdot,
                f"{self.from_node}.P_total": res_cpp.d_res_dP_total,
                f"{self.from_node}.T": res_cpp.d_res_dT,
            }}

        # Both modes: compute outlet T from actual PR seen by the element
        # (forward-propagated via compute_derived_state on downstream node)
        return res, jac
```

The `TurbineStageElement` is identical except the efficiency formula
flips: `h2 = h1 - eta_is * (h1 - h2s)`. In capacity mode a choked turbine nozzle
is the most common case — `throat_area` directly gives the choking capacity.

### H.2c — Corrected Quantities (Diagnostics)

For post-solve diagnostics the element reports standard corrected quantities:
```python
def diagnostics(self, state_in, state_out):
    T_ref, P_ref = 288.15, 101325.0   # ISA sea-level
    theta = state_in.T_total / T_ref
    delta = state_in.P_total / P_ref
    PR_actual = state_out.P_total / state_in.P_total

    return {
        "PR": PR_actual,
        "eta_is": self.eta_is,
        "W_specific": h2 - h1,
        "corrected_flow": state_in.m_dot * sqrt(theta) / delta,
        "capacity": state_in.m_dot * sqrt(state_in.T_total) / state_in.P_total,
    }
```

### H.3 — Temperature Propagation

**Key insight**: outlet temperature `T2` is *derived*, not a solver unknown. It follows
the same pattern as `CombustorNode.compute_derived_state()`:

The downstream `MomentumChamberNode` (or `PlenumNode`) calls
`compute_derived_state()` which receives upstream states. For a turbomachinery
element, the upstream state arriving at the downstream node already has the correct
T2 computed via the isentropic model. This means:

- The element computes `T2_calc` during residual evaluation
- The downstream node's `compute_derived_state` receives `T2_calc` as
  `state.T_total` from the element's exit
- **No additional solver unknown** for temperature — it's forward-propagated

This is architecturally identical to how `CombustorNode` propagates `T_ad` without
making it a solver unknown.

### H.4 — C++ Residuals + Jacobians (Required)

**Architecture rule**: every element wired into the Newton solver must provide its
residuals and analytical Jacobians from C++, with a matching C++ finite-difference
accuracy test (cf. `PipeCompressibleDerivatives`, `OrificeDerivatives` patterns in
`tests/test_solver_interface.cpp`). Python-side FD Jacobians are not acceptable for
production solver elements.

Two C++ functions are needed — one per operating mode:

#### H.4a — Prescribed-PR mode

```cpp
struct StageResult {
    double T2;           // actual outlet stagnation temperature [K]
    double T2s;          // isentropic outlet temperature [K]
    double W_specific;   // specific work h2 - h1 [J/kg]

    // Residual: P_total_out - P1 * PR = 0
    double residual;
    double d_res_dP_total_in;    // = -PR
    double d_res_dP_total_out;   // = 1

    // Derived-state Jacobians (for T propagation through mixer)
    double dT2_dT1;      // ∂T2/∂T1
    double dT2_dP1;      // ∂T2/∂P1
    double dW_dT1;       // ∂W/∂T1
    double dW_dP1;       // ∂W/∂P1
};

StageResult compressor_stage_pr_and_jacobian(
    double T1, double P1, const std::vector<double>& X,
    double P_total_out, double PR, double eta_is);

StageResult turbine_stage_pr_and_jacobian(
    double T1, double P1, const std::vector<double>& X,
    double P_total_out, double PR, double eta_is);
```

#### H.4b — Capacity mode

The capacity residual constrains mass flow:

```cpp
struct StageCapacityResult {
    double T2;           // actual outlet stagnation temperature [K]
    double T2s;          // isentropic outlet temperature [K]
    double W_specific;   // specific work h2 - h1 [J/kg]
    double PR_actual;    // observed P_total_out / P_total_in

    // Residual: m_dot * sqrt(T1) / P1 - capacity = 0
    double residual;
    double d_res_d_mdot;         // = sqrt(T1) / P1
    double d_res_dT1;            // = m_dot / (2 * sqrt(T1) * P1)
    double d_res_dP_total_in;    // = -m_dot * sqrt(T1) / P1^2

    // Derived-state Jacobians (for T propagation — same chain rule as PR mode)
    double dT2_dT1;
    double dT2_dP1;
    double dT2_dP_total_out;     // via PR_actual = P_out/P_in
    double dW_dT1;
    double dW_dP1;
};

StageCapacityResult compressor_stage_capacity_and_jacobian(
    double m_dot, double T1, double P1, const std::vector<double>& X,
    double P_total_out, double capacity, double eta_is);

StageCapacityResult turbine_stage_capacity_and_jacobian(
    double m_dot, double T1, double P1, const std::vector<double>& X,
    double P_total_out, double capacity, double eta_is);
```

**Key difference**: in capacity mode the residual depends on `m_dot`, `T1`, and
`P_total_in` (3 coupled unknowns), whereas PR mode only depends on pressures.
The T2 Jacobians now also include `∂T2/∂P_total_out` because the observed PR
changes when the downstream pressure floats.

#### H.4c — Required C++ tests

```cpp
TEST(CompressorStagePR_Derivatives, FDAccuracy) {
    // Central FD test for every Jacobian entry, tolerance < 1e-6 relative
    // Points: PR=2..30, eta=0.7..0.95, T=250..500K, air + combustion products
}

TEST(CompressorStageCapacity_Derivatives, FDAccuracy) {
    // Same FD methodology, perturb m_dot, T1, P1, P_total_out independently
}

TEST(TurbineStagePR_Derivatives, FDAccuracy) { ... }
TEST(TurbineStageCapacity_Derivatives, FDAccuracy) { ... }
```

#### H.4d — Analytical Jacobian derivation

The chain rule through the isentropic solver:
```
dT2/dT1 = (dh2/dh1) · (dT2/dh2)
         where dh2/dh1 = 1 - 1/eta_is + (1/eta_is) · (dh2s/dh1)  [compressor]
         and   dh2s/dh1 = cp(T2s) · (dT2s/dT1)
         and   dT2s/dT1 = cp(T1) / cp(T2s)   (from ds = cp·dT/T - R·dP/P = 0 at const P2)
         and   dT2/dh2  = 1 / cp(T2)
```

For the capacity residual, Jacobians are simple closed-form:
```
d_res/d_mdot = √T1 / P1
d_res/dT1    = m_dot / (2·√T1·P1)
d_res/dP1    = -m_dot·√T1 / P1²
```

All analytically tractable — no FD needed for the Jacobians themselves.

### H.5 — Network Topology Patterns

**Simple turbojet** (student exercise):
```
PressureBoundary(inlet) → Pipe → CompressorStage → CombustorNode → TurbineStage → Nozzle → PressureBoundary(exit)
                                      ↑                                  ↓
                                      └──────── ShaftConnection ─────────┘
                                              (W_compressor = W_turbine)
```

**Shaft coupling** adds a global constraint: `W_compressor + W_turbine = 0`.
This is a single scalar residual linking two elements — implementable as a
`ShaftConnection` similar to `WallConnection` for thermal coupling.

**Multi-spool** extends to multiple shafts, each with its own balance equation.

### H.6 — GUI Representation

| Component | GUI Schema | Parameters |
|-----------|-----------|------------|
| `CompressorStageElement` | `CompressorData` | `mode`, `pressure_ratio`\|`capacity`\|`throat_area`, `eta_is`, `shaft_id` |
| `TurbineStageElement` | `TurbineData` | `mode`, `pressure_ratio`\|`capacity`\|`throat_area`, `eta_is`, `shaft_id` |
| `ShaftConnection` | `ShaftData` | `connected_elements[]`, `mechanical_efficiency` |

In the GUI canvas, compressor and turbine would be rendered as distinct node shapes
(trapezoids per turbomachinery convention). The shaft connection would be a dashed
line between coupled stages.

### H.7 — Off-Design Performance Maps (Future)

For off-design analysis, `pressure_ratio` and `eta_is` become functions of corrected
mass flow and corrected speed (compressor map) or expansion ratio (turbine map).
This would follow the `make_tabulated_correlation` pattern from orifice.h:

```python
CompressorStageElement(
    ...,
    performance_map=cb.CompressorMap(
        corrected_flow=[...],    # corrected mass flow points
        corrected_speed=[...],   # corrected speed lines
        PR_table=[[...]],        # PR(corrected_flow, corrected_speed)
        eta_table=[[...]],       # eta_is(corrected_flow, corrected_speed)
    ),
)
```

This is a significant feature extension and should follow the basic constant-PR
element.

### H.8 — Implementation Dependencies

| Step | Description | Depends On | Effort |
|------|------------|-----------|--------|
| **H.0** | Bind `calc_T_from_s_mass` in pybind | G.2b | Small |
| **H.1a** | C++ `*_stage_pr_and_jacobian` (compressor + turbine) | H.0 | Medium |
| **H.1b** | C++ `*_stage_capacity_and_jacobian` (compressor + turbine) | H.0 | Medium |
| **H.2** | C++ FD accuracy tests for all 4 functions (8 test cases) | H.1a, H.1b | Small |
| **H.3** | Pybind bindings for `StageResult` + `StageCapacityResult` | H.2 | Small |
| **H.4** | Python `CompressorStageElement` + `TurbineStageElement` (both modes) | H.3 | Medium |
| **H.5** | C++ `shaft_balance_residual_and_jacobian` + FD test | H.4 | Medium |
| **H.6** | Python `ShaftConnection` | H.5 | Small |
| **H.7** | GUI schemas + graph builder (mode dropdown) | H.4 | Small |
| **H.8** | Performance maps (off-design, Mode 3) | H.4, G.5a | Large |

### H.9 — Summary

The existing C++ thermo backend (`s_mass`, `h_mass`, `calc_T_from_h_mass`) plus the
missing `calc_T_from_s_mass` provides everything needed for turbomachinery stages.
The network architecture already supports the pattern (forward T propagation via
`compute_derived_state`, pressure residuals from elements). The implementation
sequence is:

1. **Bind `calc_T_from_s_mass`** — trivial, unblocks everything
2. **C++ `*_stage_pr_and_jacobian`** — prescribed-PR mode (simplest, validates thermo chain)
3. **C++ `*_stage_capacity_and_jacobian`** — capacity mode (enables system matching)
4. **C++ FD accuracy tests** for all 4 functions (central FD, < 1e-6 relative)
5. **Pybind + Python `NetworkElement` subclasses** with mode dispatch
6. **C++ `shaft_balance_residual_and_jacobian`** + FD test, then Python `ShaftConnection`

Every solver-facing function must follow the established discipline:
**C++ implementation → C++ FD test → pybind binding → Python element → integration test.**

The **capacity mode** is what makes this feature genuinely useful: it enables
compressor–turbine matching, choked turbine nozzle modelling, and coupled gas turbine
cycle analysis. The prescribed-PR mode alone is only a teaching tool.

This would make the network solver capable of modelling complete gas turbine cycles,
which is high-value for both students (thermodynamic cycle analysis) and engineers
(secondary air system + engine cycle coupling).
