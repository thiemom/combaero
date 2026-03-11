# Combustion Elements — Design and API

This document covers the stoichiometric functions, combustion result dataclass,
and element design for the network solver's combustion and mixing elements.

---

## Solver Integration Design Rules

> **Design Rule: Solver → pybind → `solver_interface.h`**
>
> All combustion element `residuals()` methods must call Python functions that interface
> with C++ via `include/solver_interface.h`. This ensures:
> 1. All physics evaluations return both value and analytical Jacobian `(f, J)` in a single call
> 2. Derivative logic is centralized in C++, not duplicated in Python
> 3. Finite-difference fallbacks are encapsulated at C++ speed
>
> See `docs/SOLVER.md` Section 4 for the complete design rule.

### Derivative Requirements for Solver Stability

All functions called within the solver loop must satisfy:

1. **Smoothness**: Functions must be continuous and differentiable (C¹) over the
   expected operating range. Discontinuities or kinks cause Newton solver divergence.

2. **Non-throwing**: Functions must never throw exceptions during gradient probing.
   The solver may evaluate at physically unrealistic states (e.g., T < 0, P < 0)
   while computing finite-difference gradients. Return clamped/extrapolated values
   instead of throwing.

3. **Bounded gradients**: Avoid infinite or near-infinite derivatives. Clamp or
   regularize near singularities (e.g., division by zero, log(0)).

4. **Analytical Jacobians preferred**: Where possible, use `solver_interface.h`
   functions that return `(value, derivative)` tuples. This eliminates finite-difference
   noise and improves convergence.

---

## Supporting Types

> **Canonical definitions** for `MixtureState` and `NetworkElement` live in
> `NETWORK_ROADMAP.md`. The combustion functions use the Phase 1/2 subset:
> `P`, `T`, `m_dot`, `X` — `P_total` and `T_total` are not required here.

---

## Core Dataclass: `CombustionResult`

The single return type from all combustion and mixing functions. Contains the
complete thermodynamic state of the burned/mixed stream — everything the network
solver needs to continue propagating state downstream.

> **Note:** The existing C++ `CombustionState` struct (reactants + products +
> phi + mixture_fraction) lives in `include/state.h` alongside `State`,
> `Stream`, `TransportState`, `ThermoState`, `AirProperties`, and
> `CompleteState`. All state types are consolidated there as the single
> source of truth for the state design pattern.

```python
from dataclasses import dataclass

@dataclass
class CombustionResult:
    # Composition
    X: list[float]          # mole fractions [14 species], sums to 1
    Y: list[float]          # mass fractions [14 species], sums to 1
    mw: float               # mixture molecular weight [g/mol]

    # Thermodynamic state — STATIC conditions (low-velocity plenum assumption)
    # For high-velocity combustors (M > 0.3), convert to static via:
    #   T_static = combaero.T_from_stagnation(T_total, M, X)
    #   P_static = combaero.P_from_stagnation(P_total, T_total, M, X)
    T: float                # static temperature [K]
    P: float                # static pressure [Pa]
    m_dot: float            # total mass flow [kg/s]

    # Derived properties (computed via CombAero)
    h: float                # specific enthalpy [J/mol]
    cp: float               # heat capacity [J/(mol·K)]
    rho: float              # density [kg/m³]
    gamma: float            # isentropic expansion coefficient [-]
    a: float                # speed of sound [m/s]

    # Combustion diagnostics
    phi: float              # equivalence ratio [-]  (0 if no fuel)
    T_adiabatic: float      # adiabatic flame temperature [K]
    eta: float              # combustion efficiency applied [-]
    Q_released: float       # heat released [J/s] = [W]

    # Note: T_out = T_in + η·(T_adiabatic − T_in) is a linear blend
    # approximation. It is not enthalpy-consistent for large η deviations
    # from 1. For accurate off-design predictions, invert h_mix via
    # calc_T_from_h() with a scaled Q_released instead.
```

---

## Function 1: `stoichiometric_products`

Computes burned gas composition from reactant compositions using atom balance
(C, H, O, N conservation). Assumes complete combustion. No iteration required.

```python
def stoichiometric_products(
    X_fuel: list[float],
    X_oxidiser: list[float],
    phi: float,
) -> list[float]:
    """
    Compute burned gas mole fractions from atom balance.

    Assumes complete combustion:
      - All C → CO2
      - All H → H2O
      - Excess O2 remains
      - N2 passes through inert

    Parameters
    ----------
    X_fuel : list[float]
        Mole fractions of fuel stream [14 species].
    X_oxidiser : list[float]
        Mole fractions of oxidiser stream [14 species].
    phi : float
        Equivalence ratio. phi < 1: lean, phi = 1: stoichiometric, phi > 1: rich.

    Returns
    -------
    list[float]
        Burned gas mole fractions [14 species], normalised to sum to 1.

    Notes
    -----
    Uses molecular_structures from thermo_transport_data.h via CombAero
    for atomic composition (C, H, O, N atoms per species).
    Rich combustion (phi > 1) produces CO and H2 as incomplete products —
    for accurate rich mixture composition use equilibrium chemistry (Cantera).
    """
```

**Implementation approach:**

1. Compute stoichiometric fuel/air ratio `FAR_stoich` from atom balance
2. Scale fuel stream by `phi * FAR_stoich`
3. Combine fuel + oxidiser → reactant atom counts [C, H, O, N]
4. Assign products: C→CO2, H→H2O (pairs), remaining O→O2, N→N2
5. For rich (φ > 1): distribute excess fuel to CO and H₂ using the
   water-gas shift assumption (all excess C → CO, all excess H → H₂).
   This is an approximation — for accurate rich or high-T dissociation
   use Cantera equilibrium chemistry.
6. Normalise to mole fractions

---

## Function 2: `mix_streams`

Mixes two streams by mass flow weighting. No reaction — pure mixing of
composition and enthalpy. Used by `MixingElement` and as the first step
in `CombustionElement`.

> **⚠️ Chamber Automatic Mixing Requirement**
>
> All chamber-type nodes (`PlenumNode`, `MomentumChamberNode`, and any future
> chamber variants) with multiple upstream connections must **automatically**
> act as mixing nodes. This is not optional — the solver must enforce:
>
> 1. **Species mass conservation**: $X_{\text{node},k} = \frac{\sum_{\text{in}} \dot{m}_i X_{i,k}}{\sum_{\text{in}} \dot{m}_i}$
> 2. **Energy conservation**: $\sum_{\text{in}} \dot{m}_i h_i = \sum_{\text{out}} \dot{m}_j h_j$
>
> The C++ `cb.mix(streams)` function implements both balances correctly and should
> be used for computing the mixed state. For solver residuals, decompose into
> `cb.enthalpy_and_jacobian()` calls per the design rule (see SOLVER.md Section 4).
>
> **Implication**: Every interior chamber node with >1 upstream element implicitly
> adds energy and species residuals to the system. The solver must detect this
> topology and inject the appropriate conservation equations automatically.
>
> **Applies to**: `PlenumNode`, `MomentumChamberNode`, `CombustorNode`, and any
> `NetworkNode` subclass representing a finite-volume chamber where streams mix.

```python
def mix_streams(
    state_1: MixtureState,
    state_2: MixtureState,
) -> CombustionResult:
    """
    Mix two streams conserving mass, species, and enthalpy.

    The mixed temperature is found by solving:
        h_mix = (m1·h1 + m2·h2) / (m1 + m2)
    using combaero.calc_T_from_h().

    Parameters
    ----------
    state_1 : MixtureState
        First stream (P, T, m_dot, X).
    state_2 : MixtureState
        Second stream (P, T, m_dot, X).
        Pressure must be compatible with state_1 (solver enforces this).

    Returns
    -------
    CombustionResult
        Mixed stream state. phi=0, eta=1, Q_released=0.
    """
```

---

## Pressure-Loss Correlation Hook (C++ / Python)

The C++ functions `combustion_state()` and `combustion_state_from_streams()` accept an
optional user-supplied pressure-loss callable via `PressureLossCorrelation`.

### Why a hook instead of a fixed `delta_P_frac`?

Real combustors have pressure losses that depend on operating conditions (φ, T_ad,
fuel splits). The hook lets the user supply any correlation without modifying CombAero.

### Loop-free design

All fields in `PressureLossContext` are **loop-free** — they depend only on inlet
conditions and combustion outputs, never on the unknown outlet pressure. This makes
the hook safe to call inside a Newton solver without creating an implicit loop.

### Derivative Requirements for User-Supplied Pressure Loss

> **⚠️ Stability Warning**: User-supplied pressure-loss correlations are a primary
> candidate for solver instability. The hook is called on every residual evaluation
> (potentially thousands of times per solve).

User-supplied `PressureLossCorrelation` callables **must** satisfy:

1. **Smoothness (C¹)**: The correlation must be continuous and differentiable with
   respect to all `PressureLossContext` fields. Avoid `if/else` branches that create
   discontinuities — use smooth blending functions instead.

2. **Non-throwing**: Never raise exceptions. The solver may probe extreme values
   (e.g., θ < 0, φ > 10) during gradient estimation. Return clamped values:
   ```python
   def my_loss(ctx):
       theta_safe = max(0.0, min(ctx.theta, 5.0))  # Clamp to valid range
       return 0.02 + 0.005 * theta_safe
   ```

3. **Bounded output**: Return values must stay in [0, 1). Values ≥ 1 would imply
   P_out ≤ 0, causing solver failure. Clamp explicitly:
   ```python
   def my_loss(ctx):
       raw = 0.02 + 0.01 * ctx.theta
       return min(raw, 0.5)  # Never exceed 50% loss
   ```

4. **Analytical derivative (recommended)**: If the solver requires Jacobians for
   the pressure-loss term, the correlation should also return `∂(ΔP/P)/∂θ`,
   `∂(ΔP/P)/∂φ`, etc. Future API may support a `(value, jacobian_dict)` return
   signature for this purpose.

5. **Performance**: The hook is called in the inner solver loop. Avoid expensive
   operations (file I/O, network calls, large allocations). Pure arithmetic only.

| Context field | Loop-free? | Reason |
|---|---|---|
| `state_in` (T, P, X) | ✅ | Upstream, known before solve |
| `phi` | ✅ | From ṁ_fuel/ṁ_air, both known |
| `T_ad` | ✅ | From h-balance, no P dependence (complete combustion) |
| `X_products` | ✅ | From atom balance, no P dependence |
| `theta = T_ad/T_in − 1` | ✅ | Dimensionless temperature rise, no P dependence |
| `mdot_fuel`, `mdot_air` | ✅ | Fixed stream inputs |
| Outlet P, ρ | ❌ | **Not provided** — would create a loop |

For equilibrium combustion `T_ad` and `X_products` shift slightly with P (~1% for
typical 4% ΔP/P). The hook uses P_in for the equilibrium call — this is acceptable
for engineering accuracy.

### C++ usage

```cpp
auto result = combustion_state(
    X_fuel, X_ox, phi, T_reactants, P, "",
    CombustionMethod::Complete,
    [](const PressureLossContext& c) {
        return 0.02 + 0.005 * c.theta;   // 2% base + 0.5% per unit θ
    }
);
// result.products.thermo.P = P * (1 - hook_return_value)
```

### Python usage

```python
def my_loss(ctx):
    return 0.02 + 0.005 * ctx.theta

result = cb.combustion_state(X_fuel, X_ox, phi=0.8, T_reactants=700, P=500000,
                              pressure_loss=my_loss)

# Available context fields:
# ctx.phi, ctx.T_ad, ctx.theta, ctx.T_in, ctx.P_in, ctx.X_in, ctx.X_products
# ctx.mdot_fuel, ctx.mdot_air  (populated by combustion_state_from_streams)
```

---

## Function 3: `combustion_from_streams`

The primary combustion function for the network solver. Takes a mixed inlet
state (air + fuel already mixed via `mix_streams`) and returns the burned
state. This is the function `CombustionElement` calls internally.

```python
def combustion_from_streams(
    state_air: MixtureState,
    state_fuel: MixtureState,
    eta: float = 1.0,
    delta_P_frac: float = 0.04,
) -> CombustionResult:
    """
    Compute burned gas state from separate air and fuel streams.

    Steps:
      1. Mix air + fuel streams (mass flow weighted)
      2. Compute equivalence ratio phi from atom balance
      3. Compute burned composition via stoichiometric_products()
      4. Compute adiabatic flame temperature via calc_T_from_h()
      5. Apply combustion efficiency: T_out = T_in + eta*(T_adiabatic - T_in)
      6. Apply pressure drop: P_out = P_in * (1 - delta_P_frac)

    Parameters
    ----------
    state_air : MixtureState
        Oxidiser stream (P, T, m_dot, X).
    state_fuel : MixtureState
        Fuel stream (P, T, m_dot, X). m_dot is the fuel mass flow rate.
    eta : float
        Combustion efficiency [0, 1]. Default 1.0 (complete combustion).
    delta_P_frac : float
        Fractional total pressure drop across burner. Default 0.04 (4%).

    Returns
    -------
    CombustionResult
        Complete burned gas state including diagnostics.
    """
```

---

## Function 4: `combustion_from_phi`

Convenience function for cases where equivalence ratio is the natural input
(preliminary design, parametric studies) rather than separate mass flows.

```python
def combustion_from_phi(
    state_air: MixtureState,
    X_fuel: list[float],
    phi: float,
    eta: float = 1.0,
    delta_P_frac: float = 0.04,
) -> CombustionResult:
    """
    Compute burned gas state from equivalence ratio.

    Derives fuel mass flow from phi and FAR_stoich, then delegates
    to combustion_from_streams().

    Parameters
    ----------
    state_air : MixtureState
        Oxidiser stream.
    X_fuel : list[float]
        Fuel composition (mole fractions).
    phi : float
        Equivalence ratio. phi=1: stoichiometric.
    eta : float
        Combustion efficiency.
    delta_P_frac : float
        Fractional pressure drop.

    Returns
    -------
    CombustionResult
        Complete burned gas state.
    """
```

---

## Network Element: `CombustionElement`

Wraps the combustion functions behind the `NetworkElement` interface. The
solver calls `residuals()` — the element handles all chemistry internally.

> **⚠️ Design Rule Compliance**: The example below shows the conceptual residual
> structure. For actual solver integration, the implementation must:
> 1. Use `solver_interface.h` functions (e.g., `cb.enthalpy_and_jacobian()`) for
>    all thermodynamic evaluations
> 2. Return both residuals and analytical Jacobians: `(list[float], dict[int, dict[str, float]])`
> 3. Avoid calling high-level convenience functions like `combustion_from_streams()`
>    directly in the residual evaluation — instead, decompose into primitive
>    `solver_interface.h` calls that provide derivatives

```python
@dataclass
class CombustionElement(NetworkElement):
    eta: float = 1.0
    delta_P_frac: float = 0.04
    fuel_bc: MassFlowBoundary = field(default=None)  # fixed fuel stream; must be set before solve

    def residuals(self,
                  state_in: MixtureState,
                  state_out: MixtureState) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Conceptual structure — actual implementation uses solver_interface.h
        # functions with Jacobians (see SOLVER.md Section 4)

        # Energy balance: h_in + Q_comb = h_out
        h_in, dh_in_dT = cb.enthalpy_and_jacobian(state_in.T, state_in.X).result
        h_out, dh_out_dT = cb.enthalpy_and_jacobian(state_out.T, state_out.X).result

        # ... compute T_ad, X_products via C++ combustion functions ...

        res = [
            state_out.T     - T_ad,           # energy balance
            state_out.P     - P_out,          # pressure drop
            state_out.m_dot - m_dot_total,    # mass balance
            # species balance residuals...
        ]

        jac = {
            0: {f"{self.to_node}.T": 1.0, ...},  # ∂R_energy/∂T_out, etc.
            1: {f"{self.to_node}.P": 1.0, ...},
            # ... Jacobian entries for all residuals ...
        }

        return res, jac

    def n_equations(self) -> int:
        return 3 + 14  # T, P, m_dot, 14 species
```

---

## Network Element: `MixingElement`

Pure mixing, no reaction. Used for dilution air, cooling flow injection,
stream merging.

> **⚠️ Design Rule Compliance**: Same requirements as `CombustionElement` —
> use `solver_interface.h` functions with Jacobians, not `mix_streams()`.
> The C++ `cb.mix()` function can be used for post-processing but not in
> the solver residual loop (see SOLVER.md Section 4).

```python
@dataclass
class MixingElement(NetworkElement):
    """Mixes two inlet streams into one outlet. No reaction."""

    def residuals(self,
                  state_in_1: MixtureState,
                  state_in_2: MixtureState,
                  state_out: MixtureState) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Energy balance: m1*h1 + m2*h2 = m_out*h_out
        h1, dh1_dT = cb.enthalpy_and_jacobian(state_in_1.T, state_in_1.X).result
        h2, dh2_dT = cb.enthalpy_and_jacobian(state_in_2.T, state_in_2.X).result
        h_out, dh_out_dT = cb.enthalpy_and_jacobian(state_out.T, state_out.X).result

        m1, m2 = state_in_1.m_dot, state_in_2.m_dot
        m_out = m1 + m2

        # Energy residual
        energy_res = m1 * h1 + m2 * h2 - m_out * h_out

        res = [
            energy_res,
            state_out.P - min(state_in_1.P, state_in_2.P),  # pressure (simplified)
            state_out.m_dot - m_out,
            # species balance: m1*X1 + m2*X2 = m_out*X_out
        ]

        jac = {
            0: {
                f"{self.to_node}.T": -m_out * dh_out_dT,
                # ... other Jacobian entries ...
            },
            # ...
        }

        return res, jac

    def n_equations(self) -> int:
        return 3 + 14
```

---

## What CombAero Already Provides

All thermodynamic calls in the combustion functions map directly to existing
CombAero functions.

### For Solver Residuals (with Jacobians via `solver_interface.h`)

| Need | Solver function | Returns |
|---|---|---|
| Enthalpy + derivative | `cb.enthalpy_and_jacobian(T, X)` | `(h, ∂h/∂T)` |
| Density + derivatives | `cb.density_and_jacobians(T, P, X)` | `(ρ, ∂ρ/∂T, ∂ρ/∂P)` |
| Viscosity + derivative | `cb.viscosity_and_jacobians(T, P, X)` | `(μ, ∂μ/∂T)` |

### For Convenience / Post-Processing

| Need | CombAero function |
|---|---|
| Mixture enthalpy | `combaero.h(T, X)` |
| Invert enthalpy → T | `combaero.calc_T_from_h(h_target, X)` |
| Density | `combaero.density(T, P, X)` |
| Speed of sound | `combaero.speed_of_sound(T, X)` |
| Heat capacity | `combaero.cp(T, X)` |
| Gamma | `combaero.isentropic_expansion_coefficient(T, X)` |
| Mole → mass fractions | `combaero.mole_to_mass(X)` |
| Molecular weight | `combaero.mwmix(X)` |
| Equivalence ratio from composition | `combaero.equivalence_ratio_mole(X, X_fuel, X_ox)` |
| Stream mixing | `combaero.mix(streams, P_out)` |
| Atomic composition | `molecular_structures` in `thermo_transport_data.h` |

---

## Missing: `stoichiometric_products` Implementation

This is the only function not yet in CombAero. It requires:

1. Reading atomic composition (C, H, O, N) per species from
   `molecular_structures` in `thermo_transport_data.h`
2. Pure atom-balance arithmetic — no iteration, no equilibrium chemistry
3. ~50–80 lines of Python

**Module layout** (all network solver code lives outside the C++ core):

```
combaero_network/
    combustion/
        stoichiometry.py      # stoichiometric_products()
        combustion.py         # combustion_from_streams(), combustion_from_phi()
        mixing.py             # mix_streams()
    elements/
        combustion_element.py # CombustionElement
        mixing_element.py     # MixingElement
    core/
        network_element.py    # NetworkElement base class, MixtureState
```

For rich mixtures (φ > 1) or high-temperature dissociation, the atom-balance
approximation breaks down. In those cases, delegate to Cantera's equilibrium
solver via the existing `cantera_validation_tests` infrastructure.

---

## Validation Targets

| Test case | Expected result | Reference |
|---|---|---|
| CH4/air, phi=1.0, T_in=700K, P=5bar | T_ad ≈ 2230 K | Cantera GRI-Mech |
| CH4/air, phi=0.6, T_in=700K, P=5bar | T_ad ≈ 1800 K | Cantera GRI-Mech |
| H2/air, phi=1.0, T_in=300K, P=1bar | T_ad ≈ 2380 K | Cantera GRI-Mech |
| Dilution mixing: 50% hot + 50% cold | T_mix = enthalpy-weighted mean | Analytical |
| Two-stream mix, same composition | X_out = X_in, T_out = weighted mean | Analytical |
