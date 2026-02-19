# Combustion Elements — Design and API

This document covers the stoichiometric functions, combustion result dataclass,
and element design for the network solver's combustion and mixing elements.

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

    # Thermodynamic state
    T: float                # temperature [K]
    P: float                # pressure [Pa]
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

```python
@dataclass
class CombustionElement(NetworkElement):
    eta: float = 1.0
    delta_P_frac: float = 0.04
    fuel_bc: MassFlowBoundary = field(default=None)  # fixed fuel stream; must be set before solve

    def residuals(self,
                  state_in: MixtureState,
                  state_out: MixtureState) -> list[float]:
        result = combustion_from_streams(
            state_air=state_in,
            state_fuel=self.fuel_bc.state,
            eta=self.eta,
            delta_P_frac=self.delta_P_frac,
        )
        return [
            state_out.T     - result.T,       # energy balance
            state_out.P     - result.P,       # pressure drop
            state_out.m_dot - result.m_dot,   # mass balance
            *(state_out.X[i] - result.X[i]   # species balance
              for i in range(14)),
        ]

    def n_equations(self) -> int:
        return 3 + 14  # T, P, m_dot, 14 species
```

---

## Network Element: `MixingElement`

Pure mixing, no reaction. Used for dilution air, cooling flow injection,
stream merging.

```python
@dataclass
class MixingElement(NetworkElement):
    """Mixes two inlet streams into one outlet. No reaction."""

    def residuals(self,
                  state_in_1: MixtureState,
                  state_in_2: MixtureState,
                  state_out: MixtureState) -> list[float]:
        result = mix_streams(state_in_1, state_in_2)
        return [
            state_out.T     - result.T,
            state_out.P     - result.P,
            state_out.m_dot - result.m_dot,
            *(state_out.X[i] - result.X[i] for i in range(14)),
        ]

    def n_equations(self) -> int:
        return 3 + 14
```

---

## What CombAero Already Provides

All thermodynamic calls in the combustion functions map directly to existing
CombAero functions:

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
