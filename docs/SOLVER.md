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
