# CombAero Network Solver Architecture

This document tracks the topological strategies and mathematical numerical methods for solving interconnected flow element networks in Python.

## 1. Open Source Solver Ecosystem

Developing a monolithic Newton solver from scratch is ill-advised for large systems. The following permissive, open-source Python numerical optimizers serve as architectural templates or potential backends:

| Framework | License | Paradigm | Pros / Cons for CombAero |
| :--- | :--- | :--- | :--- |
| **SciPy (`scipy.optimize.root`)** | BSD | Black-box array function | **Pros**: Native to NumPy, `hybrd` (MINPACK) is incredibly robust, easy to deploy without heavy installations, handles finite-difference automatically. <br>**Cons**: No intrinsic sparsity intelligence; gradients scale horribly for $N > 1000$ networks without supplying an explicit sparse Jacobian. |
| **OpenMDAO** | Apache 2.0 | Multidisciplinary systems | **Pros**: Designed explicitly for connecting "components" together in graphs, handles analytic gradients (MAUD algorithm) flawlessly. <br>**Cons**: Steep learning curve, component API requires writing dedicated derivatives which we lack for empirical friction factors. |
| **CasADi** | LGPL | Symbolic AD & NLP | **Pros**: Incredibly fast, handles sparsity intrinsically via automatic differentiation, integrates seamlessly with IPOPT solver. <br>**Cons**: Requires formulating equations as CasADi symbolic expressions, not opaque C++ compiled functions. Cannot digest our pre-compiled C++ CombAero wheel directly. |
| **GEKKO** | MIT | Differential Algebraic | **Pros**: Solves coupled DAEs gracefully. <br>**Cons**: Same as CasADi; it needs symbolic/internal representations of the equations, making opaque C++ bindings difficult to embed. |

**Verdict:** Our physics kernels are hidden inside a compiled C++ PyBind11 wheel. Because of this opacity to symbolic Automatic Differentiation (AD), frameworks like CasADi and GEKKO cannot see "inside" the functions to compute gradients.
We must either use **SciPy's `root` solver with finite-differences** (acceptable for Phase 1/2) or migrate to **OpenMDAO** treating our elements as "Explicit Components" relying on numerical, sparse Jacobian approximations.

## 2. Graph State Architecture

The user identified a crucial trajectory: the network must evolve from pure cold flow to reacting multicomponent flow with heat transfer.
The internal data structure mapped by the solver must minimize the conversion cost between the "1D Array" the solver wants, and the physical parameters the C++ functions need.

### Phase 1: Incompressible Flow Only
* **Variables**: Fluids are isothermal, non-reacting, with predefined properties.
* **Nodes (Junctions)**: Carry unknown **P** (1 DOF per node).
* **Edges (Elements)**: Carry unknown **$\dot{m}$** (1 DOF per edge).
* **Data Structure**: A flat 1D NumPy array `[P1, P2... P_n, mdot1... mdot_m]`.
* **State Object**: Since composition (X) and temperature (T) are constant, creating C++ `State` objects is unnecessary. The network explicitly passes `(T_const, P_i, X_const, mdot_i)` down into the C++ kernels.

### Phase 2: Reacting Flow (Combustion)
* **Variables**: Fluids undergo adiabatic composition shifts. The system must conserve mass, momentum, and **energy** (enthalpy), and species must be perfectly tracked.
* **Nodes**: Must mix arbitrary enthalpy streams to find a unified downstream temperature and composition. Unknowns per node: **P, h, X** ($2+N_{species}$ DOFs per node).
* **Edges**: Transport state. **$\dot{m}$** remains the primary edge unknown.
* **State Object**: Iterations evaluate `T = f(h, P, X)` heavily. Passing properties as decoupled primitives `(h, P, X)` into C++ is strictly required to avoid excessive object initialization latency. The `MixtureState` tracking class should remain a pure Python `dataclass` wrapping views into the monolithic NumPy state array.

### Phase 3: Flow, Combustion, and Heat Transfer
* **Variables**: Enthalpy is no longer conserved; walls leak energy $\dot{Q}$.
* **Edges**: $\dot{Q}$ becomes a required residual at each element boundary. Elements might need discrete spatial integration (1D Fanno / Rayleigh) to properly evaluate temperature profiles.
* **State Object**: Nodes treat `h` exactly as Phase 2. Edges must dynamically resolve local temperature gradients to evaluate wall fluxes.

## 3. Jacobian Stability

The speed and reliability of `scipy.optimize.root(method='hybrd')` hinges entirely on the quality of the Jacobian matrix. While SciPy can approximate this via finite differences ($\mathbf{f(x + \Delta x)}$), this approach is fundamentally fragile to numerical noise and dimensionally scales poorly ($O(N)$ evaluations per matrix construction).

**3.1 Analytical Jacobians (The Gold Standard)**
As noted in the roadmap, Phase 1 focuses on closed-form algebraic fluid mechanics (e.g., Orifice Flow: $\dot{m} \propto \sqrt{\Delta P}$). These equations possess exact, continuous analytical derivatives:
* $\frac{\partial \dot{m}}{\partial P_{in}} = \frac{1}{2} C_d A \sqrt{\frac{2 \rho}{P_{in} - P_{back}}}$
* $\frac{\partial h}{\partial T} = c_p$

Furthermore, virtually all of CombAero's empirical correlations are formulated as smooth, explicit power-laws or logarithms perfectly suited for exact chain-rule differentiation:
* **Friction (Petukhov):** $f = (0.79 \ln Re - 1.64)^{-2} \Rightarrow \frac{\partial f}{\partial Re} = -1.58 \cdot Re^{-1} \cdot (0.79 \ln Re - 1.64)^{-3}$
* **Heat Transfer (Dittus-Boelter):** $Nu = 0.023 Re^{0.8} Pr^{0.4} \Rightarrow \frac{\partial Nu}{\partial Re} = 0.0184 Re^{-0.2} Pr^{0.4}$
* **Pressure Loss ($K$-factor):** $\Delta P = \zeta \frac{1}{2} \rho v^2 \Rightarrow \frac{\partial \Delta P}{\partial v} = \zeta \rho v$

The solver architecture must prioritize computing these exact analytical gradients in Python and feeding them directly to `scipy.optimize.root(jac=...)`. This provides maximum mathematical stability, renders finite-difference noise irrelevant, and dramatically accelerates the convergence sequence.

**3.2 C++ Solver API: (f, J) Optimizations**
To avoid duplicating derivative logic across Python and C++ (violating DRY), the C++ physics kernels themselves should expose a dedicated solver API. For instance, rather than Python manually calculating $\partial \dot{m} / \partial P$ for an orifice, a C++ function like `orifice_mdot_derivative(dP, ...)` returns a struct or `std::tuple`:
`(value, derivative_wrt_target)`

This allows `scipy.optimize.root` to call the Python binding exactly once per block, unpack the value and analytical gradient instantaneously, and inject them into the sparse Newton solver. Where strictly analytical derivatives are impossible (e.g., highly implicit implicit functions), the C++ layer internally evaluates a tightly bounded central-difference numerical derivative at compiled C++ speeds, returning the identical `(f, J)` signature. This permanently encapsulates finite difference noise-handling inside the core rather than exposing the global Python solver to it.

**3.3 C++ Safety Protections**
When analytical gradients are too complex to derive (e.g., implicit friction factor networks) and the solver falls back to finite differences, sudden jumps or exceptions during gradient probing will immediately destroy the linesearch. This is precisely why the CombAero C++ engine returns robust approximations or clamps values at extreme boundaries instead of aggressively throwing `std::out_of_range` exceptions. The solver *must* be allowed to blindly probe unrealistic pressures (e.g. $P < P_{back}$) to calculate a local gradient pointing it back towards physical reality without the program crashing.

**Sparsity pattern**: Nodes in a pipe network represent a sparse graph map. Pressure $P_i$ at node $i$ only affects the flow $\dot{m}_{ij}$ in elements physically connected to it. The solver architecture must eventually pre-evaluate this sparsity pattern (e.g., via `scipy.sparse`); otherwise, finite-differencing a 500-node network requires $O(N^2)$ independent full-system evaluations to construct the dense $J$ matrix, resulting in catastrophic scaling delays.

---

## 4. Energy Conservation Implementation Plan

This section details the concrete implementation steps for adding energy conservation to the network solver, covering **mixing**, **combustion**, and **heat transfer**.

> **Design Rule: Solver → pybind → `solver_interface.h`**
>
> The `NetworkSolver` and element `residuals()` methods must call Python functions that interface with C++ via `include/solver_interface.h`. This ensures:
> 1. All physics evaluations return both value and analytical Jacobian `(f, J)` in a single call
> 2. Derivative logic is centralized in C++, not duplicated in Python
> 3. Finite-difference fallbacks are encapsulated at C++ speed
>
> **Exception:** Convenience wrappers in `components.py` (e.g., `AreaDischargeCoefficientConnectionElement`) may use higher-level C++ APIs like `cb.mix()` or `cb.htc_pipe()` for non-solver purposes (e.g., post-processing, initialization). However, residual and Jacobian evaluation in the solver loop must use `solver_interface.h` functions.
>
> **Performance Note:** Even for convenience wrappers, if they are called within the solver loop, the additional latency from higher-level APIs must be proven negligible via profiling before adoption. The solver may evaluate residuals thousands of times per solve; any per-call overhead compounds significantly.

### 4.1 Current State Summary

| Capability | Status |
|:---|:---|
| Mass conservation at nodes | ✅ Implemented |
| Momentum (pressure-flow) residuals | ✅ Implemented |
| Sparse analytical Jacobian assembly | ✅ Implemented |
| Energy conservation at nodes | ❌ Not implemented |
| Species tracking through network | ❌ Partial (X stored but not solved) |
| Heat transfer in elements | ❌ Stub only (C++ `htc_pipe`, `nusselt_*` exposed) |
| Combustion nodes | ❌ Stub only (C++ `mix`, `complete_combustion` exposed) |

### 4.2 Governing Equations

**Node Energy Balance (Mixing)**

At each interior node, the steady-state energy balance is:

$$\sum_{\text{in}} \dot{m}_i h_i = \sum_{\text{out}} \dot{m}_j h_j$$

where $h$ is the mass-specific stagnation enthalpy $h_0 = h(T) + \frac{v^2}{2}$.

For a plenum node ($v \approx 0$), this simplifies to:

$$\sum_{\text{in}} \dot{m}_i h(T_i, X_i) = \sum_{\text{out}} \dot{m}_j h(T_{\text{node}}, X_{\text{node}})$$

The mixed outlet temperature $T_{\text{node}}$ and composition $X_{\text{node}}$ are unknowns.

**Species Conservation (Mixing)**

$$X_{\text{node},k} = \frac{\sum_{\text{in}} \dot{m}_i X_{i,k}}{\sum_{\text{in}} \dot{m}_i}$$

This is an explicit algebraic equation (not a residual) that can be evaluated directly once mass flows are known.

**Element Heat Transfer**

For a pipe element with wall heat flux:

$$\dot{Q} = \dot{m} (h_{\text{out}} - h_{\text{in}}) = h_{\text{conv}} A_{\text{wall}} \Delta T_{\text{lm}}$$

where $\Delta T_{\text{lm}}$ is the log-mean temperature difference.

**Combustion**

For a combustor node, the outlet state satisfies:

$$h_{\text{out}}(T_{\text{ad}}, X_{\text{products}}) = h_{\text{in}}(T_{\text{in}}, X_{\text{reactants}})$$

with $X_{\text{products}}$ determined by the combustion model (complete or equilibrium).

### 4.3 Implementation Steps

#### Step 1: Extend `MixtureState` with Enthalpy Helpers

Add convenience methods to `MixtureState`:

```python
def stagnation_enthalpy(self) -> float:
    """h0 = h(T) + 0.5 * v^2, where v = m_dot / (rho * A)."""
    return cb.h_mass(self.T_total, self.X)

def total_enthalpy_flux(self) -> float:
    """m_dot * h0 [W]."""
    return self.m_dot * self.stagnation_enthalpy()
```

#### Step 2: Add Temperature Unknown to Interior Nodes

Modify `PlenumNode.unknowns()` to include temperature:

```python
def unknowns(self) -> list[str]:
    return [f"{self.id}.P", f"{self.id}.P_total", f"{self.id}.T"]
```

Update `_get_node_state()` in solver to unpack `T` from the solution vector.

#### Step 3: Add Energy Residual to Node Residuals

In `NetworkSolver._residuals_and_jacobian()`, after the mass conservation residual, add:

```python
# Energy Conservation: Sum(m_dot_in * h_in) - Sum(m_dot_out * h_out) = 0
energy_res_idx = len(res)
h_in_total = 0.0
h_out_total = 0.0

for elem in upstream_elems:
    state_up = self._get_node_state(self.network.nodes[elem.from_node], x)
    m_dot_elem = ...  # from element unknown
    h_in_total += m_dot_elem * state_up.stagnation_enthalpy()

for elem in downstream_elems:
    m_dot_elem = ...
    h_out_total += m_dot_elem * state.stagnation_enthalpy()

res.append(h_in_total - h_out_total)
```

**Jacobian entries:**
- $\partial R_E / \partial T_{\text{node}} = -\sum_{\text{out}} \dot{m}_j \cdot c_{p,\text{mass}}(T_{\text{node}}, X_{\text{node}})$
- $\partial R_E / \partial \dot{m}_j = h_j$ (for each connected element)

Use `cb.enthalpy_and_jacobian(T, X)` to obtain $(h, \partial h / \partial T)$.

#### Step 4: Species Mixing Strategy

**For solver residuals (implicit):** Use `cb.enthalpy_and_jacobian(T, X)` from `solver_interface.h` to compute enthalpy terms and their derivatives. The energy residual (Step 3) handles enthalpy balance implicitly.

**For post-processing / initialization (explicit):** The convenience function `cb.mix(streams, P_out)` from `state.h` can be used outside the solver loop to compute mixed states:

```python
def compute_mixed_state(self, node_id: str, x: np.ndarray) -> cb.Stream:
    """Compute mixed state at a node using C++ enthalpy-balanced mixing.

    NOTE: This is a convenience wrapper for post-processing/initialization.
    Solver residuals use cb.enthalpy_and_jacobian() per the design rule.
    """
    upstream_elems = self.network.get_upstream_elements(node_id)
    streams = []

    for elem in upstream_elems:
        m_dot = ...  # from x
        state_up = self._get_node_state(self.network.nodes[elem.from_node], x)
        stream = cb.Stream()
        stream.set_T(state_up.T).set_P(state_up.P).set_X(state_up.X).set_mdot(m_dot)
        streams.append(stream)

    if streams:
        return cb.mix(streams)  # Returns Stream with mixed T, X, mdot
    return None
```

#### Step 5: Use `solver_interface.h` Heat Transfer Functions in `PipeElement`

Per the design rule, residual evaluation must use `solver_interface.h` functions that return `(value, Jacobian)`. Use `cb.nusselt_and_jacobian_gnielinski(Re, Pr, f)` (or the appropriate correlation variant) rather than the convenience `cb.htc_pipe()`.

Add heat transfer residual when `htc_model != "none"`:

```python
def residuals(self, state_in, state_out):
    # ... existing pressure drop residual ...

    if self.htc_model != "none" and self.t_wall is not None:
        v = abs_mdot / (rho * A)
        k = cb.thermal_conductivity(state_in.T, state_in.P, state_in.X)
        Pr = cb.prandtl(state_in.T, state_in.P, state_in.X)

        # Use solver_interface.h function for (Nu, dNu/dRe)
        corr = cb.nusselt_and_jacobian_gnielinski(Re, Pr, f)
        Nu, dNu_dRe = corr.result
        h_conv = Nu * k / self.diameter

        # Log-mean temperature difference
        dT_in = self.t_wall - state_in.T
        dT_out = self.t_wall - state_out.T
        if abs(dT_in - dT_out) > 1e-6:
            dT_lm = (dT_in - dT_out) / np.log(dT_in / dT_out)
        else:
            dT_lm = dT_in  # Avoid log(1) singularity

        # Heat transfer rate
        A_wall = np.pi * self.diameter * self.length
        Q_dot = h_conv * A_wall * dT_lm

        # Energy residual: m_dot * (h_out - h_in) - Q_dot = 0
        h_in, dh_in_dT = cb.enthalpy_and_jacobian(state_in.T, state_in.X).result
        h_out, dh_out_dT = cb.enthalpy_and_jacobian(state_out.T, state_out.X).result
        res.append(m_dot * (h_out - h_in) - Q_dot)

        # Jacobian entries use dNu_dRe, dh_dT from solver_interface.h
```

This ensures all derivatives are computed analytically in C++ via `solver_interface.h`.

#### Step 6: Implement `CombustorNode` Energy Balance

Replace the placeholder in `CombustorNode.residuals()`:

```python
def residuals(self, state: MixtureState) -> tuple[list[float], dict]:
    # Get inlet enthalpy from upstream element
    h_in = ...  # from upstream state
    X_in = ...  # from upstream state

    # Compute adiabatic flame temperature and products
    if self.method == "complete":
        T_ad, X_products = cb.adiabatic_T_complete(T_in, P, X_in)
    else:
        T_ad, X_products = cb.adiabatic_T_equilibrium(T_in, P, X_in)

    # Residual: node temperature matches adiabatic flame temperature
    res = [state.T - T_ad]

    # Pressure loss residual
    res.append(state.P_total - (1 - self.pressure_loss_frac) * P_total_in)

    return res, jac
```

### 4.4 Jacobian Strategy for Energy Terms

| Term | Derivative | Source |
|:---|:---|:---|
| $h(T, X)$ | $\partial h / \partial T = c_{p,\text{mass}}$ | `cb.enthalpy_and_jacobian(T, X)` |
| $\rho(T, P, X)$ | $\partial \rho / \partial T$, $\partial \rho / \partial P$ | `cb.density_and_jacobians(T, P, X)` |
| $\mu(T, P, X)$ | $\partial \mu / \partial T$ | `cb.viscosity_and_jacobians(T, P, X)` |
| $Nu(Re, Pr, f)$ | $\partial Nu / \partial Re$ | `cb.nusselt_and_jacobian_*(Re, Pr, ...)` |

All required Jacobian APIs are already implemented and tested in `solver_interface.cpp`.

### 4.5 Testing Strategy

1. **Unit tests for energy residuals**: Verify $\sum \dot{m} h = 0$ at mixing nodes
2. **Adiabatic pipe**: Confirm $T_{\text{out}} = T_{\text{in}}$ when `htc_model="none"`
3. **Isothermal wall**: Verify $T_{\text{out}} \to T_{\text{wall}}$ for long pipes
4. **Combustor validation**: Compare $T_{\text{ad}}$ against Cantera for known fuel/air mixtures
5. **Jacobian accuracy**: Extend `test_network_jacobian.py` to include energy terms

### 4.6 Phase Summary

| Phase | Scope | Key Unknowns per Node |
|:---|:---|:---|
| **Current** | Isothermal, fixed composition | $P$, $P_{\text{total}}$ |
| **Phase 2** | Energy conservation, mixing | $P$, $P_{\text{total}}$, $T$ |
| **Phase 3** | Heat transfer in elements | $P$, $P_{\text{total}}$, $T$ (+ element $T_{\text{out}}$) |
| **Phase 4** | Reacting flow, species tracking | $P$, $P_{\text{total}}$, $T$, $X$ (14 species) |
