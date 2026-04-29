"""Combustion functions for network solver.

All thermodynamic calculations delegate to C++ via pybind11:
- mixer_from_streams_and_jacobians()  -- adiabatic stream mixing
- adiabatic_T_complete_and_jacobian_T_from_streams()  -- complete combustion
- adiabatic_T_equilibrium_and_jacobians_from_streams() -- equilibrium combustion
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import combaero._core as core

from .combustion_result import CombustionResult

if TYPE_CHECKING:
    from combaero.network.components import NetworkMixtureState


def mix_streams(
    state_1: NetworkMixtureState,
    state_2: NetworkMixtureState,
) -> CombustionResult:
    """Mix two streams conserving mass, species, and enthalpy.

    Delegates to mixer_from_streams_and_jacobians (C++ via pybind11).
    """
    import combaero as cb

    streams = [
        core.MassStream(state_1.m_dot, state_1.Tt, state_1.Pt, list(state_1.Y)),
        core.MassStream(state_2.m_dot, state_2.Tt, state_2.Pt, list(state_2.Y)),
    ]
    mix = core.mixer_from_streams_and_jacobians(streams)

    T = mix.T_mix
    P = mix.P_total_mix
    Y = list(mix.Y_mix)
    X = list(cb.mass_to_mole(Y))
    mw = cb.mwmix(X)

    return CombustionResult(
        X=X,
        Y=Y,
        mw=mw,
        T=T,
        P=P,
        m_dot=state_1.m_dot + state_2.m_dot,
        h=cb.h_mass(T, X),
        cp=cb.cp(T, X),
        rho=cb.density(T, P, X),
        gamma=cb.isentropic_expansion_coefficient(T, X),
        a=cb.speed_of_sound(T, X),
        phi=0.0,
        T_adiabatic=T,
        eta=1.0,
        Q_released=0.0,
    )


def combustion_from_streams(
    state_air: NetworkMixtureState,
    state_fuel: NetworkMixtureState,
    method: str = "complete",
    eta: float = 1.0,
    delta_P_frac: float = 0.04,
) -> CombustionResult:
    """Compute burned gas state from separate air and fuel streams.

    Steps:
      1. Mix streams in C++ to get inlet state (T_in, Y_mix)
      2. Compute adiabatic flame temperature and product composition in C++
      3. Apply combustion efficiency: T_out = T_in + eta*(T_ad - T_in)
      4. Apply pressure drop: P_out = P_in * (1 - delta_P_frac)

    Parameters
    ----------
    state_air : NetworkMixtureState
        Oxidiser stream (Pt, Tt, m_dot, Y).
    state_fuel : NetworkMixtureState
        Fuel stream (Pt, Tt, m_dot, Y).
    method : str
        'complete' or 'equilibrium'. Default 'complete'.
    eta : float
        Combustion efficiency [0, 1]. Default 1.0.
    delta_P_frac : float
        Fractional total pressure drop. Default 0.04.
    """
    import combaero as cb

    streams = [
        core.MassStream(state_air.m_dot, state_air.Tt, state_air.Pt, list(state_air.Y)),
        core.MassStream(state_fuel.m_dot, state_fuel.Tt, state_fuel.Pt, list(state_fuel.Y)),
    ]
    P = state_air.P

    mix_in = core.mixer_from_streams_and_jacobians(streams)
    T_in = mix_in.T_mix
    X_mix_in = list(cb.mass_to_mole(list(mix_in.Y_mix)))

    if method == "complete":
        res = core.adiabatic_T_complete_and_jacobian_T_from_streams(streams, P)
    else:
        res = core.adiabatic_T_equilibrium_and_jacobians_from_streams(streams, P)

    T_ad = res.T_mix
    Y_products = list(res.Y_mix)
    X_products = list(cb.mass_to_mole(Y_products))

    phi = cb.equivalence_ratio_mole(
        X_mix_in,
        list(cb.mass_to_mole(list(state_fuel.Y))),
        list(cb.mass_to_mole(list(state_air.Y))),
    )

    T_out = T_in + eta * (T_ad - T_in)
    P_out = P * (1.0 - delta_P_frac)
    m_dot = state_air.m_dot + state_fuel.m_dot
    mw = cb.mwmix(X_products)

    Q_released = m_dot * (cb.h_mass(T_in, X_mix_in) - cb.h_mass(T_in, X_products))

    return CombustionResult(
        X=X_products,
        Y=Y_products,
        mw=mw,
        T=T_out,
        P=P_out,
        m_dot=m_dot,
        h=cb.h_mass(T_out, X_products),
        cp=cb.cp(T_out, X_products),
        rho=cb.density(T_out, P_out, X_products),
        gamma=cb.isentropic_expansion_coefficient(T_out, X_products),
        a=cb.speed_of_sound(T_out, X_products),
        phi=phi,
        T_adiabatic=T_ad,
        eta=eta,
        Q_released=Q_released,
    )


def combustion_from_phi(
    state_air: NetworkMixtureState,
    X_fuel: list[float],
    phi: float,
    method: str = "complete",
    eta: float = 1.0,
    delta_P_frac: float = 0.04,
) -> CombustionResult:
    """Compute burned gas state from equivalence ratio.

    Derives fuel mass flow from phi and FAR_stoich, then delegates
    to combustion_from_streams().

    Parameters
    ----------
    state_air : NetworkMixtureState
        Oxidiser stream.
    X_fuel : list[float]
        Fuel mole fractions.
    phi : float
        Equivalence ratio. phi=1: stoichiometric.
    method : str
        'complete' or 'equilibrium'. Default 'complete'.
    eta : float
        Combustion efficiency.
    delta_P_frac : float
        Fractional pressure drop.
    """
    import combaero as cb
    from combaero.network.components import NetworkMixtureState

    air_cb = cb.Stream()
    air_cb.set_T(state_air.T).set_P(state_air.Pt).set_X(
        list(cb.mass_to_mole(list(state_air.Y)))
    ).set_mdot(state_air.m_dot)

    fuel_cb = cb.Stream()
    fuel_cb.set_T(state_air.T).set_X(list(X_fuel))
    fuel_cb = cb.set_fuel_stream_for_phi(phi, fuel_cb, air_cb)

    state_fuel = NetworkMixtureState(
        P=state_air.P,
        Pt=state_air.Pt,
        T=state_air.T,
        Tt=state_air.Tt,
        m_dot=fuel_cb.mdot,
        Y=list(cb.mole_to_mass(list(X_fuel))),
    )

    return combustion_from_streams(
        state_air,
        state_fuel,
        method=method,
        eta=eta,
        delta_P_frac=delta_P_frac,
    )
