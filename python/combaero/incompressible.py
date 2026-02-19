"""Incompressible flow submodule.

Provides a symmetric high-level API that mirrors ``combaero.compressible``.
All functions accept ``(T, P, X, ...)`` and return a :class:`FlowSolution`.

Swapping the import is sufficient to switch flow regimes::

    # incompressible
    from combaero import incompressible as flow
    sol = flow.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05)

    # compressible (same call signature)
    from combaero import compressible as flow
    sol = flow.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05)

Functions
---------
pipe_flow
    Darcy-Weisbach pipe flow with a user-supplied friction factor.
pipe_flow_rough
    Darcy-Weisbach pipe flow with friction factor from roughness + Re.
orifice_flow
    Incompressible orifice flow (Bernoulli + Cd).
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from ._core import (
    IncompressibleFlowSolution,
    orifice_flow_thermo,
)
from ._core import (
    pipe_flow as _pipe_flow,
)
from ._core import (
    pipe_flow_rough as _pipe_flow_rough,
)
from ._flow_solution import FlowSolution


def _to_flow_solution(sol: IncompressibleFlowSolution, Cd: float = math.nan) -> FlowSolution:
    """Convert an ``IncompressibleFlowSolution`` to the unified ``FlowSolution``."""
    return FlowSolution(
        regime="incompressible",
        mdot=sol.mdot,
        v=sol.v,
        dP=sol.dP,
        Re=sol.Re,
        rho=sol.rho,
        f=sol.f,
        Cd=Cd,
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
) -> FlowSolution:
    """Incompressible pipe flow with a constant Darcy friction factor.

    Evaluates density from ``(T, P, X)`` internally and applies the
    Darcy-Weisbach equation.

    Parameters
    ----------
    T:
        Static temperature [K].
    P:
        Static pressure [Pa].
    X:
        Mole fractions [-].
    u:
        Bulk flow velocity [m/s].
    L:
        Pipe length [m].
    D:
        Pipe inner diameter [m].
    f:
        Darcy friction factor [-].

    Returns
    -------
    FlowSolution
        ``regime="incompressible"``.  ``M``, ``T_out``, ``P_out``, ``h0``,
        ``choked``, ``L_choke``, and ``Cd`` are not applicable and carry
        ``nan`` / ``False``.
    """
    sol = _pipe_flow(T, P, list(X), u=u, L=L, D=D, f=f)
    return _to_flow_solution(sol)


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
) -> FlowSolution:
    """Incompressible pipe flow with roughness-based friction factor.

    Evaluates density and dynamic viscosity from ``(T, P, X)`` internally,
    computes the Reynolds number, derives the Darcy friction factor from the
    specified correlation, and applies the Darcy-Weisbach equation.

    Parameters
    ----------
    T:
        Static temperature [K].
    P:
        Static pressure [Pa].
    X:
        Mole fractions [-].
    u:
        Bulk flow velocity [m/s].
    L:
        Pipe length [m].
    D:
        Pipe inner diameter [m].
    roughness:
        Absolute wall roughness [m].  ``0.0`` for a hydraulically smooth pipe.
    correlation:
        Friction factor correlation.  One of ``"haaland"`` (default),
        ``"serghides"``, ``"colebrook"``, ``"petukhov"``.

    Returns
    -------
    FlowSolution
        ``regime="incompressible"``.
    """
    sol = _pipe_flow_rough(
        T, P, list(X), u=u, L=L, D=D, roughness=roughness, correlation=correlation
    )
    return _to_flow_solution(sol)


def orifice_flow(
    T: float,
    P: float,
    X: Sequence[float],
    *,
    P_back: float,
    A: float,
    Cd: float = 1.0,
) -> FlowSolution:
    """Incompressible orifice flow.

    Evaluates density from ``(T, P, X)`` internally and applies the
    incompressible orifice equation::

        mdot = Cd * A * sqrt(2 * rho * (P - P_back))

    Parameters
    ----------
    T:
        Upstream static temperature [K].
    P:
        Upstream static pressure [Pa].
    X:
        Mole fractions [-].
    P_back:
        Downstream static pressure [Pa].
    A:
        Orifice (throat) area [m^2].
    Cd:
        Discharge coefficient [-].  Default ``1.0``.

    Returns
    -------
    FlowSolution
        ``regime="incompressible"``.  ``f`` carries ``Cd`` for convenience
        (matching the raw ``IncompressibleFlowSolution`` convention).
        ``Cd`` is also populated explicitly.
    """
    sol = orifice_flow_thermo(T, P, list(X), P_back=P_back, A=A, Cd=Cd)
    return _to_flow_solution(sol, Cd=Cd)
