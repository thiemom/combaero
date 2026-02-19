"""Compressible flow submodule.

Provides a symmetric high-level API that mirrors ``combaero.incompressible``.
All functions accept ``(T, P, X, ...)`` and return a :class:`FlowSolution`.

Swapping the import is sufficient to switch flow regimes::

    # compressible
    from combaero import compressible as flow
    sol = flow.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05)

    # incompressible (same call signature)
    from combaero import incompressible as flow
    sol = flow.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05)

Functions
---------
pipe_flow
    Fanno (adiabatic compressible) pipe flow with a constant friction factor.
pipe_flow_rough
    Fanno pipe flow with friction factor from roughness + local Re.
nozzle_flow
    Isentropic nozzle flow.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from ._core import (
    CompressibleFlowSolution,
    FannoSolution,
)
from ._core import (
    fanno_pipe as _fanno_pipe,
)
from ._core import (
    fanno_pipe_rough as _fanno_pipe_rough,
)
from ._core import (
    nozzle_flow as _nozzle_flow,
)
from ._flow_solution import FlowSolution


def _fanno_to_flow_solution(sol: FannoSolution, store_profile: bool = False) -> FlowSolution:
    """Convert a ``FannoSolution`` to the unified ``FlowSolution``.

    ``sol`` must have been computed with ``store_profile=True`` so that the
    last ``FannoStation`` is available for outlet M and velocity.
    """
    dP = sol.inlet.P - sol.outlet.P
    last = sol.profile[-1]
    return FlowSolution(
        regime="compressible",
        mdot=sol.mdot,
        v=last.u,
        dP=dP,
        Re=sol.Re_in,
        rho=sol.inlet.rho,
        f=sol.f_avg,
        M=last.M,
        T_out=sol.outlet.T,
        P_out=sol.outlet.P,
        h0=sol.h0,
        choked=sol.choked,
        L_choke=sol.L_choke if sol.choked else math.nan,
        profile=list(sol.profile) if store_profile else [],
    )


def _nozzle_to_flow_solution(sol: CompressibleFlowSolution) -> FlowSolution:
    """Convert a ``CompressibleFlowSolution`` (nozzle) to ``FlowSolution``."""
    return FlowSolution(
        regime="compressible",
        mdot=sol.mdot,
        v=sol.v,
        M=sol.M,
        T_out=sol.outlet.T,
        P_out=sol.outlet.P,
        h0=sol.stagnation.h,
        choked=sol.choked,
    )


def pipe_flow(
    T: float,
    P: float,
    X: Sequence[float],
    *,
    u: float,
    L: float,
    D: float,
    f: float,
    n_steps: int = 100,
    store_profile: bool = False,
) -> FlowSolution:
    """Compressible (Fanno) pipe flow with a constant Darcy friction factor.

    Integrates the Fanno flow equations (mass + energy + momentum) along the
    pipe using RK4.  The friction factor is held constant throughout.

    Parameters
    ----------
    T:
        Inlet static temperature [K].
    P:
        Inlet static pressure [Pa].
    X:
        Mole fractions [-].
    u:
        Inlet bulk velocity [m/s].
    L:
        Pipe length [m].
    D:
        Pipe inner diameter [m].
    f:
        Darcy friction factor [-] (constant along pipe).
    n_steps:
        Number of RK4 integration steps.
    store_profile:
        If ``True``, the returned ``FlowSolution.profile`` contains a list of
        ``FannoStation`` objects along the pipe axis.

    Returns
    -------
    FlowSolution
        ``regime="compressible"``.  ``Cd`` is not applicable and carries
        ``nan``.
    """
    sol = _fanno_pipe(T, P, u, L, D, f, list(X), n_steps=n_steps, store_profile=True)
    return _fanno_to_flow_solution(sol, store_profile=store_profile)


def pipe_flow_rough(
    T: float,
    P: float,
    X: Sequence[float],
    *,
    u: float,
    L: float,
    D: float,
    roughness: float = 0.0,
    correlation: str = "haaland",
    n_steps: int = 100,
    store_profile: bool = False,
) -> FlowSolution:
    """Compressible (Fanno) pipe flow with roughness-based variable friction.

    At each RK4 stage the local friction factor is recomputed from the local
    Reynolds number and wall roughness using the specified correlation.

    Parameters
    ----------
    T:
        Inlet static temperature [K].
    P:
        Inlet static pressure [Pa].
    X:
        Mole fractions [-].
    u:
        Inlet bulk velocity [m/s].
    L:
        Pipe length [m].
    D:
        Pipe inner diameter [m].
    roughness:
        Absolute wall roughness [m].  ``0.0`` for a hydraulically smooth pipe.
    correlation:
        Friction factor correlation.  One of ``"haaland"`` (default),
        ``"serghides"``, ``"colebrook"``.
    n_steps:
        Number of RK4 integration steps.
    store_profile:
        If ``True``, the returned ``FlowSolution.profile`` contains a list of
        ``FannoStation`` objects along the pipe axis.

    Returns
    -------
    FlowSolution
        ``regime="compressible"``.  ``f`` carries the length-averaged friction
        factor.  ``Re`` carries the inlet Reynolds number.
    """
    sol = _fanno_pipe_rough(
        T,
        P,
        u,
        L,
        D,
        roughness,
        list(X),
        correlation=correlation,
        n_steps=n_steps,
        store_profile=True,
    )
    return _fanno_to_flow_solution(sol, store_profile=store_profile)


def nozzle_flow(
    T: float,
    P: float,
    X: Sequence[float],
    *,
    P_back: float,
    A_eff: float,
    tol: float = 1e-8,
    max_iter: int = 50,
) -> FlowSolution:
    """Isentropic nozzle flow.

    Parameters
    ----------
    T:
        Stagnation temperature [K].
    P:
        Stagnation pressure [Pa].
    X:
        Mole fractions [-].
    P_back:
        Back (downstream) pressure [Pa].
    A_eff:
        Effective throat area [m^2].
    tol:
        Convergence tolerance for the Mach number iteration.
    max_iter:
        Maximum number of Newton iterations.

    Returns
    -------
    FlowSolution
        ``regime="compressible"``.  ``dP``, ``Re``, ``rho``, ``f``, ``Cd``,
        and ``L_choke`` are not applicable and carry ``nan``.
    """
    sol = _nozzle_flow(T, P, P_back, A_eff, list(X), tol=tol, max_iter=max_iter)
    return _nozzle_to_flow_solution(sol)
