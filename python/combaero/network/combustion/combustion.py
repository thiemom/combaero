"""Combustion functions for network solver.

UPDATED: Now uses proper C++ Jacobian functions via pybind11:
- adiabatic_T_complete_and_jacobian_T()
- adiabatic_T_equilibrium_and_jacobians()

These functions return (T_ad, dT_ad_dT_in, X_products) with analytical Jacobians,
complying with SOLVER.md Section 4 design rules.
"""

from typing import TYPE_CHECKING

import combaero._core as core

from .combustion_result import CombustionResult
from .stoichiometry import stoichiometric_products

if TYPE_CHECKING:
    from combaero.network.components import MixtureState


def mix_streams(
    state_1: "MixtureState",
    state_2: "MixtureState",
) -> CombustionResult:
    """
    Mix two streams conserving mass, species, and enthalpy.

    The mixed temperature is found by solving:
        h_mix = (m1*h1 + m2*h2) / (m1 + m2)
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
    import combaero as cb

    # Validate inputs
    if state_1.P != state_2.P:
        raise ValueError("Streams must have compatible pressures for mixing")

    m1, m2 = state_1.m_dot, state_2.m_dot
    m_total = m1 + m2

    if m_total <= 0:
        raise ValueError("Total mass flow must be positive")

    # Mass-flow weighted composition
    Y_mix = [(m1 * y1 + m2 * y2) / m_total for y1, y2 in zip(state_1.Y, state_2.Y, strict=False)]

    # Mass-weighted enthalpy
    h1 = state_1.enthalpy()
    h2 = state_2.enthalpy()
    h_mix = (m1 * h1 + m2 * h2) / m_total

    # Find temperature that gives this enthalpy
    # Convert mass enthalpy [J/kg] to molar enthalpy [J/mol] because cb.calc_T_from_h expects molar enthalpy
    X_mix = cb.mass_to_mole(Y_mix)
    mw = cb.mwmix(X_mix)
    h_mix_molar = h_mix * mw / 1000.0

    try:
        T_mix = cb.calc_T_from_h(h_mix_molar, X_mix)
    except Exception as e:
        # Fallback to mass-weighted average if enthalpy inversion fails
        # This should be removed once the enthalpy inversion is working properly
        m1, m2 = state_1.m_dot, state_2.m_dot
        m_total = m1 + m2
        T_mix = (m1 * state_1.T + m2 * state_2.T) / m_total
        # Log warning for debugging
        print(f"Warning: calc_T_from_h failed, using fallback: {e}")

    # Compute derived properties
    mw = cb.mwmix(X_mix)
    rho = cb.density(T_mix, state_1.P, X_mix)
    cp = cb.cp(T_mix, X_mix)
    gamma = cb.isentropic_expansion_coefficient(T_mix, X_mix)
    a = cb.speed_of_sound(T_mix, X_mix)

    return CombustionResult(
        X=X_mix,
        Y=Y_mix,
        mw=mw,
        T=T_mix,
        P=state_1.P,
        m_dot=m_total,
        h=h_mix,
        cp=cp,
        rho=rho,
        gamma=gamma,
        a=a,
        phi=0.0,  # No combustion
        T_adiabatic=T_mix,  # No temperature rise from combustion
        eta=1.0,  # Perfect mixing efficiency
        Q_released=0.0,  # No heat release
    )


def combustion_from_streams(
    state_air: "MixtureState | dict",
    state_fuel: "MixtureState | dict",
    method: str = "complete",
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
    state_air : MixtureState or dict
        Oxidiser stream (P, T, m_dot, Y).
    state_fuel : MixtureState or dict
        Fuel stream (P, T, m_dot, Y). m_dot is the fuel mass flow rate.
    method : str
        Combustion method: 'complete' or 'equilibrium'. Default 'complete'.
    eta : float
        Combustion efficiency [0, 1]. Default 1.0 (complete combustion).
    delta_P_frac : float
        Fractional total pressure drop across burner. Default 0.04 (4%).

    Returns
    -------
    CombustionResult
        Complete burned gas state including diagnostics.
    """
    import combaero as cb

    # Convert dictionaries to MixtureState objects if needed
    if isinstance(state_air, dict):
        # Import here to avoid circular import
        from combaero.network.components import MixtureState as LocalMixtureState

        state_air = LocalMixtureState(**state_air)
    if isinstance(state_fuel, dict):
        # Import here to avoid circular import
        from combaero.network.components import MixtureState as LocalMixtureState

        state_fuel = LocalMixtureState(**state_fuel)

    # Validate inputs
    if state_air.P != state_fuel.P:
        raise ValueError("Air and fuel streams must have the same pressure")

    # Step 1: Mix air + fuel streams
    mixed = mix_streams(state_air, state_fuel)

    # Step 2: Compute equivalence ratio (using mole fractions for the API)
    phi = cb.equivalence_ratio_mole(
        cb.mass_to_mole(mixed.Y), cb.mass_to_mole(state_fuel.Y), cb.mass_to_mole(state_air.Y)
    )

    # Step 3: Compute burned composition
    X_products = stoichiometric_products(
        cb.mass_to_mole(state_fuel.Y), cb.mass_to_mole(state_air.Y), phi
    )

    # Step 4: Compute adiabatic flame temperature using proper Jacobian function
    # Use the new C++ function that returns (T_ad, dT_ad_dT_in, X_products)
    # Note: We need to pass the method parameter explicitly since mixed is a CombustionResult
    combustion_method = method  # Use the method parameter passed to this function

    if combustion_method == "complete":
        T_adiabatic, dT_ad_dT_in, X_products = core.adiabatic_T_complete_and_jacobian_T(
            mixed.T, mixed.P, mixed.X
        )
    else:  # equilibrium
        T_adiabatic, dT_ad_dT_in, X_products = core.adiabatic_T_equilibrium_and_jacobians(
            mixed.T, mixed.P, mixed.X
        )

    # Store the analytical Jacobian for potential use in residuals
    # TODO: Pass this to the network element for proper Jacobian assembly

    # Step 5: Apply combustion efficiency
    T_out = mixed.T + eta * (T_adiabatic - mixed.T)

    # Step 6: Apply pressure drop
    P_out = mixed.P * (1 - delta_P_frac)

    # Compute derived properties at output conditions
    rho = cb.density(T_out, P_out, X_products)
    h = cb.h_mass(T_out, X_products)
    cp = cb.cp(T_out, X_products)
    gamma = cb.isentropic_expansion_coefficient(T_out, X_products)
    a = cb.speed_of_sound(T_out, X_products)
    mw = cb.mwmix(X_products)
    Y_products = cb.mole_to_mass(X_products)

    # Compute heat released: enthalpy difference at same T_in between reactants and products
    Q_released = mixed.m_dot * (mixed.h - cb.h_mass(mixed.T, X_products))

    return CombustionResult(
        X=X_products,
        Y=Y_products,
        mw=mw,
        T=T_out,
        P=P_out,
        m_dot=mixed.m_dot,
        h=h,
        cp=cp,
        rho=rho,
        gamma=gamma,
        a=a,
        phi=phi,
        T_adiabatic=T_adiabatic,
        eta=eta,
        Q_released=Q_released,
    )


def combustion_from_phi(
    state_air: "MixtureState",
    X_fuel: list[float],
    phi: float,
    method: str = "complete",
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
    method : str
        Combustion method: 'complete' or 'equilibrium'. Default 'complete'.
    eta : float
        Combustion efficiency.
    delta_P_frac : float
        Fractional pressure drop.

    Returns
    -------
    CombustionResult
        Complete burned gas state.
    """
    import combaero as cb

    # Create fuel stream with appropriate mass flow for target phi
    # Use CombAero to set the fuel mass flow

    # Create streams for CombAero's set_fuel_stream_for_phi function
    fuel_stream = cb.Stream()
    fuel_stream.set_T(state_air.T).set_X(X_fuel)

    air_stream = cb.Stream()
    air_stream.set_T(state_air.T).set_P(state_air.P).set_X(cb.mass_to_mole(state_air.Y)).set_mdot(
        state_air.m_dot
    )

    # Set fuel stream mass flow to achieve target phi
    fuel_stream_phi = cb.set_fuel_stream_for_phi(phi, fuel_stream, air_stream)

    # Create a simple state dictionary for fuel instead of using cb.MixtureState
    # to avoid circular import
    fuel_state_dict = {
        "P": state_air.P,
        "P_total": state_air.P_total,
        "T": state_air.T,
        "T_total": state_air.T_total,
        "m_dot": fuel_stream_phi.mdot,
        "Y": cb.mole_to_mass(X_fuel),
    }

    # Delegate to combustion_from_streams using the state dictionary
    # We'll create proper MixtureState objects inside combustion_from_streams
    return combustion_from_streams(
        state_air,
        fuel_state_dict,
        method=state_air.method if hasattr(state_air, "method") else "complete",
        eta=eta,
        delta_P_frac=delta_P_frac,
    )
